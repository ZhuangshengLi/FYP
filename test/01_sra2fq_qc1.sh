###########################################  
date 
echo  -e "\n \n \n  111#  move files !!! \n \n \n  " 
cd ~/fyp/test/raw/sra/
cat ~/fyp/test/idname | while read id
do 
mv $id/*  ./ 
rm -rf $id/ 
done 
date 


echo  -e "\n \n \n  111#  sra>>>fq !!! \n \n \n  "
mkdir -p ~/fyp/test/raw/fq/
cd ~/fyp/test/raw/fq/
pwd
ls  ~/fyp/test/raw/sra/*.sra |while read id 
do
echo " PROCESS $(basename $id) "
fasterq-dump -3 -e 12 -O ./ $id
pigz -p 12   ~/fyp/test/raw/fq/*q
done
date


echo -e " \n \n \n  111# qc 1 !!! \n \n \n "       
mkdir ~/fyp/test/raw/qc1/
cd  ~/fyp/test/raw/qc1/
pwd
ls ~/fyp/test/raw/fq/* | xargs fastqc -t 12 -o  ./
multiqc ./

echo -e  " \n 111#  ALL  Work Done!!! \n "
date
