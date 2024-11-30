############################
echo -e  " \n \n \n 333#  Align !!! hisat2 !!!\n \n \n "
date
########set#### ###set#### ###set####   
index='~/fyp/reference/index/hisat/mm10/genome'

mkdir  -p   ~/fyp/test/align/flag
cd ~/fyp/test/align/
pwd

cat ~/fyp/test/idname | while read id;
do
         echo "333#  ${id} ${id} ${id}  is on the hisat2 Working !!!"
################ paired ###############################         
            hisat2 -t -p 12 -x  $index \
       -1 ~/fyp/test/clean/${id}_1_val_1.fq.gz  \
       -2  ~/fyp/test/clean/${id}_2_val_2.fq.gz  -S  ${id}.sam
######################################################
################Single################################
#              hisat2 -t -p 12 -x  $index  \
#                     -U ~/fyp/test/clean/${id}_trimmed.fq.gz \
#                     -S  ~/fyp/test/clean/${id}.sam  
########################################################           

### sam2bam and remove sam ###
    echo -e  " ${id} sam2bam and remove sam   "
    samtools sort -@ 12 -o  ~/fyp/test/align/${id}_sorted.bam  ${id}.sam
    rm ${id}.sam
done

#### samtools index and flagstat ####
echo -e  " \n \n \n samtools index and flagstat \n  " 
cd -p  ~/fyp/test/align/flag
pwd

#ls ~/fyp/test/align/*.bam | xargs  -i  samtools index  -@  12  {}    
ls ~/fyp/test/align/*.bam  | while read id ;\
do
        samtools flagstat -@ 12 $id > $(basename $id '.bam').flagstat
done
multiqc ./

echo -e " \n \n \n  333# ALL  Work Done !!!\n \n \n "
date
