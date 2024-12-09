gtf="$HOME/fyp/reference/gtf/gencode/gencode.vM25.annotation.gtf"

mkdir ~/fyp/test/name_map
### gene_id to gene_name
grep 'gene_id' $gtf | awk -F 'gene_id \"' '{print $2}' |awk -F '\"' '{print $1}' >gene_id_tmp
grep 'gene_id' $gtf | awk -F 'gene_name \"' '{print $2}' |awk -F '\"' '{print $1}' >gene_name_tmp
paste gene_id_tmp gene_name_tmp >last_tmp
uniq last_tmp >~/fyp/test/name_map/g2s_vm25_gencode.txt
rm *_tmp

### transcript_id to gene_name
grep 'transcript_id' $gtf | awk -F 'transcript_id \"' '{print $2}' |awk -F '\"' '{print $1}' >gene_id_tmp
grep 'transcript_id' $gtf | awk -F 'gene_name \"' '{print $2}' |awk -F '\"' '{print $1}' >gene_name_tmp
paste gene_id_tmp gene_name_tmp >last_tmp
uniq last_tmp >~/fyp/test/name_map/t2s_vm25_gencode.txt
rm *_tmp
