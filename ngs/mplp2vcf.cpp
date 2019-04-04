#include <limits>
#include <algorithm>
#include <ctype.h>
#include <string>
#include <iostream>
#include <sstream>
#include <map>
#include <set>
#include <vector>
#include <iterator>

extern "C" {
#include "sam.h"
#include "bgzf.h"
#include "faidx.h"
}


// -----------------------------------------------------------
// -- Structure to store and merge statistics for each sample
// -----------------------------------------------------------
struct fr_pair:std::pair<int,int>{
  int min() const{return std::min(first,second);}
  int sum() const{return first+second;}
  fr_pair& operator+=(const fr_pair& b) {
    first += b.first;
    second += b.second;
    return *this;
  }
};



struct sample_gt {
  typedef std::map<std::string,fr_pair> gt_map_t;
  fr_pair dp;
  gt_map_t h;
  std::vector<std::string> alt;
  
  // default constructor
  sample_gt() {
  }
  
  // Initialize structure with a set of pileup records
  sample_gt(const bam_pileup1_t* first,const bam_pileup1_t* last) {
    // Compute DP, DPF, DPR
    dp.first = std::count_if(first,last,[](const bam_pileup1_t& o){return !bam_is_rev(o.b);});
    dp.second = std::count_if(first,last,[](const bam_pileup1_t& o){return bam_is_rev(o.b);});

    // Update counts for each alternative found
    for(;first != last;first++) {
      std::stringstream ss;
      if (first->indel<0) {
        ss << "<DEL" << -first->indel << '>';
        if (bam_is_rev(first->b)) h[ss.str()].second++; else h[ss.str()].first++;
      } else if (first->indel>0) {
        for (int j = 0;j <= first->indel; j++) {
          char c = seq_nt16_str[bam_seqi(bam_get_seq(first->b), first->qpos + j)];
          ss << c;
        }
        if (bam_is_rev(first->b)) h[ss.str()].second++; else h[ss.str()].first++;
      } else if (!first->is_del) {
        char c = seq_nt16_str[bam_seqi(bam_get_seq(first->b), first->qpos)];
        ss << c;
        if (bam_is_rev(first->b)) h[ss.str()].second++; else h[ss.str()].first++;
      }
    }
  }
  
  // Merge data of two libraries
  sample_gt& operator+=(const sample_gt& b) {
    this->dp += b.dp;
    for(auto kv:b.h) {
      auto& e = h[kv.first];
      e.first += kv.second.first;
      e.second += kv.second.second;
    }
    return *this;
  }
  
  // 
  void update_alt() {
    alt.clear();
    for(const auto& a:h) alt.push_back(a.first);
    std::sort(alt.begin(),alt.end(),[this](const std::string& a,const std::string& b){return h[a].min()>h[b].min();});
  }

  int ac(const std::string& a) const {
    auto i = h.find(a);
    return (i==h.end())?0:i->second.min();
  }

  double af(const std::string& a,double missing=0.5) const {
    return (dp.min()==0)?missing:(((double) ac(a))/dp.min());
  }

  int ac2() const { 
    return alt.size()>1?ac(alt[1]):0; 
  }
  
  int acf(const std::string& a) const {
    auto i = h.find(a);
    return (i==h.end())?0:i->second.first;
  }

  int acr(const std::string& a) const {
    auto i = h.find(a);
    return (i==h.end())?0:i->second.second;
  }
  
  int acfr(const std::string& a) const {
    auto i = h.find(a);
    return (i==h.end())?0:i->second.sum();
  }

  std::string acf(const std::vector<std::string>& alt) const {
    std::stringstream ss;
    auto a = alt.begin();
    if (a!=alt.end()) {
      ss << acf(*a);
      for(a++;a!=alt.end();a++) {
        ss << ',' << acf(*a);
      }
    }
    return ss.str();
  }
  
  std::string acr(const std::vector<std::string>& alt) const {
    std::stringstream ss;
    auto a = alt.begin();
    if (a!=alt.end()) {
      ss << acr(*a);
      for(a++;a!=alt.end();a++) {
        ss << ',' << acr(*a);
      }
    }
    return ss.str();
  }

  std::string acfr(const std::vector<std::string>& alt) const {
    std::stringstream ss;
    auto a = alt.begin();
    if (a!=alt.end()) {
      ss << acfr(*a);
      for(a++;a!=alt.end();a++) {
        ss << ',' << acfr(*a);
      }
    }
    return ss.str();
  }

};




// -----------------------------------------------------------------------------------------------
// -- The structure used to pass data to mplp_func.
// -- It stores filename, file pointer, and header structure of each input BAM file
// -----------------------------------------------------------------------------------------------
struct mplp_aux_t {
  const char* fn;
  BGZF* fp;
  bam_hdr_t *hdr;

  mplp_aux_t(const char* fn):fn(fn) {
    fp = bgzf_open(fn,"r");
    hdr = bam_hdr_read(fp);
  }
  
  ~mplp_aux_t() {
    bam_hdr_destroy(hdr);
    bgzf_close(fp);
  }
};


// -- This function is called at each iteration and can be used to process read data or skip reads
static int mplp_func(void *data, bam1_t *b) {
  mplp_aux_t *ma = (mplp_aux_t*) data;
  return bam_read1(ma->fp, b);
}

void usage(const char* argv0) {
	fprintf(stderr,"usage:%s ref.fasta in1.bam ...\n",argv0);
	exit(1);
}



// ----------------------------------------------
// -- The Main function
// ----------------------------------------------
int main(int argc, const char **argv) {

  if (argc<3) usage(argv[0]);
  faidx_t* fai = fai_load(argv[1]);
  if (!fai) usage(argv[0]);

  // -- initialize the structure storing the parameters
  std::vector<mplp_aux_t*> data;
  for (int i=2; i<argc; i++) data.push_back(new mplp_aux_t(argv[i]));
  bam_mplp_t iter = (bam_mplp_t) bam_mplp_init(data.size(), mplp_func,(void**) &data[0]);
  
  // -- print VCF header
  std::cout << "##fileformat=VCFv4.1" << std::endl;
  std::cout << "##ALT=<ID=DELx,Description=\"Deletion of length x\">" << std::endl;
  std::cout << "##ALT=<ID=REF,Description=\"Alternate is equal to REF\">" << std::endl;
  std::cout << "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">" << std::endl;
  std::cout << "##INFO=<ID=AC2,Number=1,Type=Integer,Description=\"Number of read supporting the second most frequent allele on both strand. A value of 1 means the second most frequent allele is observed at least 1 time on each strand.\">" << std::endl;
  std::cout << "##INFO=<ID=AFD,Number=1,Type=Float,Description=\"Maximum difference of allele frequency with consensus.\">" << std::endl;
  std::cout << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << std::endl;
  std::cout << "##FORMAT=<ID=DPF,Number=1,Type=Integer,Description=\"Forward Read Depth\">" << std::endl;
  std::cout << "##FORMAT=<ID=DPR,Number=1,Type=Integer,Description=\"Reverse Read Depth\">" << std::endl;
  std::cout << "##FORMAT=<ID=ACF,Number=A,Type=Integer,Description=\"Forward Allele Count\">" << std::endl;
  std::cout << "##FORMAT=<ID=ACR,Number=A,Type=Integer,Description=\"Reverse Allele Count\">" << std::endl;
  std::cout << "##FORMAT=<ID=ACFR,Number=A,Type=Integer,Description=\"ACF+ACR\">" << std::endl;
  std::cout << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
  for(auto i=data.begin();i!=data.end();i++) std::cout << '\t' << (*i)->fn;
  std::cout << std::endl;
  
  // -- print information for each position
  int tid,pos;
  int ref_tid=-1,ref_len=0;
  const char* ref = 0;
  std::vector<const bam_pileup1_t *> plp(data.size());
  std::vector<int> n_plp(plp.size());
  bam_mplp_set_maxcnt(iter,100000);
  while(bam_mplp_auto(iter, &tid, &pos, &n_plp[0], &plp[0]) > 0) {
    // -- load reference sequence
    if (tid != ref_tid) {
      if (ref) delete[] ref;
      ref_tid=tid;
      ref = faidx_fetch_seq(fai, data[0]->hdr->target_name[ref_tid], 0, 0x7fffffff, &ref_len);
    }
    std::string ref_bp(1,(ref && pos<ref_len)?ref[pos]:'N');
    
    // Fill genotypes
    std::vector<sample_gt> genotypes;
    for(size_t i=0;i<plp.size();i++) genotypes.push_back(sample_gt(plp[i],plp[i] + n_plp[i]));
    
    // Merge all genotypes
    sample_gt ALT;
    for(const auto& gt:genotypes) ALT += gt;
    ALT.update_alt();
        
    // Compute INFO.AC2 field
    int ac2 = ALT.ac2();
    //if (ac2<3) continue;

    // Compute INFO.AFD field
    std::vector<double> af_delta;
    for(auto& a:ALT.alt) {
      auto i = std::max_element(genotypes.begin(),genotypes.end(),[&a](const sample_gt& x,const sample_gt& y){return x.af(a,-1)<y.af(a,-1);});
      auto j = std::min_element(genotypes.begin(),genotypes.end(),[&a](const sample_gt& x,const sample_gt& y){return x.af(a,2)<y.af(a,2);});
      af_delta.push_back(i->af(a,1) - j->af(a,0));
    }
    auto max_af_delta_ptr = std::max_element(af_delta.begin(),af_delta.end());
    double max_af_delta = (max_af_delta_ptr==af_delta.end())?-1.0:*max_af_delta_ptr;


    // Update genotypes alternatives separately
    for(auto& gt:genotypes) gt.update_alt();
    
    // Print VCF mandatory fields
    std::cout << data[0]->hdr->target_name[tid] << '\t' << pos+1 << "\t.\t" << ref_bp << '\t';
    auto a=ALT.alt.begin();
    if (a!=ALT.alt.end()) {
      std::cout << (*a==ref_bp?"<REF>":*a);
      for(a++;a!=ALT.alt.end();a++) std::cout << ',' << (*a==ref_bp?"<REF>":*a);
    }
    std::cout << "\t.\t.\tAC2=" << ac2 << ";AFD=" << max_af_delta;
    
    // -- print VCF FORMAT and genotypes
    std::cout << "\tGT:DPF:DPR:ACF:ACR:ACFR";
    for(auto& gt:genotypes) {
      int GT_idx = 0;
      if (gt.alt.size()>0) {
        auto i = find(ALT.alt.begin(),ALT.alt.end(),gt.alt[0]);
        GT_idx = std::distance(ALT.alt.begin(),i) + 1;
      }
      std::cout << '\t';
      if (GT_idx<1) std::cout << '.'; else std::cout << GT_idx;
      std::cout << ':' << gt.dp.first << ':' << gt.dp.second << ':' << gt.acf(ALT.alt) << ':' << gt.acr(ALT.alt) << ':' << gt.acfr(ALT.alt);
    }
    std::cout << std::endl;
  }
  
  // -- clean memory
  bam_mplp_destroy(iter);
  fai_destroy(fai);
  for (auto i=data.begin(); i!=data.end();i++) delete *i;
  
  return 0;
}


