###################################################
echo -e  '\n \n \n ### salmon quant is Working !!! \n \n'
###set##set###set#############
index="$HOME/fyp/reference/index/salmon/grcm38/"

mkdir ~/fyp/test/salmon
cd ~/fyp/test/salmon
pwd
cat ~/fyp/test/idname | while read id 
do
            echo " '\n  !!!!!! Processing sample $id !!!!! '\n" 
########single#############################################
# salmon quant -i $index -l A \
#             -r ~/fyp/test/clean/$(basename $id)_trimmed.fq.gz \
#             -p 12  -o  $(basename $id)_quant
#######paired#############################################
salmon quant -i $index -l A \
 -1 ~/fyp/test/clean/${id}_1_val_1.fq.gz \
 -2 ~/fyp/test/clean/${id}_2_val_2.fq.gz  \
 -p 12  --output ${id}_quant
###############################################################     
done

multiqc ./

echo -e " \n \n \n !!!!ALL WORK DONE !!!!!!!!!!!!!!!!!!!!! \n"
date
