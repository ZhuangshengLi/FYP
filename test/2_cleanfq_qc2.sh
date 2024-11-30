##############################################
echo -e " \n \n \n 222# Clean ! trim_galore !!! \n \n \n"
date
mkdir ~/fyp/test/clean/
cd ~/fyp/test/clean/
pwd

##########single###########################################################################
#cd ..
#ls ~/fyp/test/raw/fq/*.f* | while read id 
#do 
#       trim_galore -q 25 -j 4 --phred33 --length 35 --stringency 3 \
#               --gzip  -o ~/fyp/test/clean/  $id 
#done
#
##########paired###########################################################################
# First store the file paths and filenames of file_1 and file_2 respectively, then merge them into a two-column format and save as config.
  ls ~/fyp/test/raw/fq/*_1*  >1
  ls ~/fyp/test/raw/fq/*_2*  >2
  paste 1 2 >config
  cat config | while read id
  do
    arr=($id)
    fq1=${arr[0]}
    fq2=${arr[1]}
   trim_galore  -j 4  -q 25  --phred33 --length 35 --stringency 3 \
               --paired --gzip -o ~/fyp/test/clean/  $fq1 $fq2
  done
###########################################################################################                  


echo -e "\n \n \n 222# qc2 check the cleaned results!!! \n \n \n"
mkdir  ~/fyp/test/clean/qc2
cd ~/fyp/test/clean/qc2
pwd

ls ~/fyp/test/clean/*f*.gz | xargs fastqc -t 12 -o   ./
multiqc   ./

echo -e " \n 222# ALL  Work Done !!! \n "
date
