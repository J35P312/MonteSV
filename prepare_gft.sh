zgrep "	gene	" Homo_sapiens.GRCh37.87.gtf.gz | grep "_coding" | bedtools sort -i /dev/stdin > Homo_sapiens.GRCh37.87.clean.gtf
bgzip Homo_sapiens.GRCh37.87.clean.gtf
tabix -p gff Homo_sapiens.GRCh37.87.clean.gtf.gz
