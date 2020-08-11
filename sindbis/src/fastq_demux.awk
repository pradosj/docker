#!/usr/bin/awk -f
#usage: zcat file.fastq.gz | awk -v indexmap=index_map.tsv -f fastq_demux.awk -

BEGIN{
  FS=OFS="\t";
  while(getline < indexmap) {
    map[$1]=$4 "_" $3 OFS $5
  }
}
NR%4==2 {
  idx = substr($0,89,6);
  umi = substr($0,95,12);
  barcode = substr($0,1,32);
  if (idx in map) {
    print map[idx],idx,umi,barcode;
  } else {
    print "Unmapped_NNNNNN",0,idx,umi,barcode;
  }
}
