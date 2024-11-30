# fetch the raw datasets
echo -e "\n \n \n prefetch sra !!! \n \n \n " 
date 
mkdir -p ~/fyp/test/raw/sra/ 
cd ~/fyp/test/raw/sra/ 
pwd 

cat  ~/fyp/test/idname | while read id ;
do       
  ( prefetch -O ./ $id & ) 
done   
