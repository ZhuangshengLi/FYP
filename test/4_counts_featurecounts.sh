###########################################
echo -e "\n \n \n444# Count #featureCounts !!! \n \n \n"
date
#####set####set###set
gtf="$HOME/fyp/reference/gtf/gencode/gencode.vM10.annotation.gtf"
#gtf="$HOME/fyp/reference/gtf/UCSC/mm10.refGene.gtf.gz"

mkdir  -p  ~/fyp/test/counts
cd  ~/fyp/test/counts/
pwd

######## single##########################################################################
#featureCounts -T 12 -a $gtf  -o counts.txt  ~/fyp/test/align/*.bam        
####### paired ###########################################################################
featureCounts -T  12  -p  -a  $gtf  -o  counts.txt  ~/fyp/test/align/*.bam
#######################################################################################            
# Copy script
cp counts.txt ~/fyp/test/R

### Generate html files to display the results
multiqc ./

echo -e " \n \n \n ALL WORK DONE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  \n "
date
