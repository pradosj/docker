

docker run -p 80:80 -it -v $(PWD)/jbrowse_data/:/usr/local/apache2/htdocs/data/ -v $(PWD)/mimol_data/:/mnt/mimol/ jbrowse
open http://`docker-machine ip default`

./bin/prepare-refseqs.pl --fasta /mnt/mimol/ncbi_genomes/GCF_001281145.1_ASM128114v1/GCF_001281145.1_ASM128114v1_genomic.fna
./bin/flatfile-to-json.pl --gff /mnt/mimol/ncbi_genomes/GCF_001281145.1_ASM128114v1/GCF_001281145.1_ASM128114v1_genomic.gff --trackType CanvasFeatures --trackLabel "SA564 GFF"


