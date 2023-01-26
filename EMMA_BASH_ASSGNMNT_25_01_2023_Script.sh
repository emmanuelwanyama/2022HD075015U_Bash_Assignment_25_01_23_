#!/usr/bin/bash 
clear
echo -e "\n \e[1;34;5m This is a script for Manuplation of SAM and VCF files with: \e[0m \n 
\e[35m 1. BCFtools #VCF file Manipulation . \n
 2. samtools #SAM file manipulation. \n
 \t\t\t\t \e[5mBy Wanyama Emmanuel Okuku 2022/HD07/5015U\n\e[0m"

read -p "Please Provide the VCF file for Analysis:    " VCF_FILE
read -p "Please Provide the SAM file for Analysis:    " SAM_FILE

echo -e "\e[1;32m \n\nThe Provided Files are: \n ${VCF_FILE} \n ${SAM_FILE}"

echo -e "\e[1;35m\n The following is the response the question with answers Provided:\n\n"

echo -e "\e[1;5;35m\n\n\n ##Manipulating VCF files##\e[0m"                               


echo -e "\e[1;33m2.What does the header section of the file contain"
cat ${VCF_FILE} | grep "^##" > Qtn1#Bcftools_header_all_output.out
cat ${VCF_FILE} | grep "^##"  | less   > Qtn1#Bcftools_header_output_less.out
echo -e "\e[1;34m"
echo -e 'Code Used: \n  cat ${VCF_FILE}  | grep "^##"  #To get the header like lines in the vcf file. \n cat ${VCF_FILE}  | grep "^##"  | less   #To View few lines at a time'


echo -e "\e[1;33m3. How many variants are in the file"
bcftools query -l ${VCF_FILE}  | head -10  
bcftools query -l ${VCF_FILE}   | wc -l
echo -e "\e[1;34m"
echo -e 'Code Used: \n bcftools query -l ${VCF_FILE}  | head -10 #To view the files \n bcftools query -l ${VCF_FILE}   | wc -l #Counting the Files'


echo -e "\e[1;33m4. How many variants are in the file"
bcftools view ${VCF_FILE}  | grep -v "^#" | wc -l
echo -e '\e[1;34m'
echo -e 'Code Used: \n bcftools view ${VCF_FILE}  | grep -v "^#" | wc -l'

echo -e "\e[1;33m5.How would you extract the chromosome, position, QualByDepth and RMSMappingQuality fields? Save the output to a tab-delimited file."
bcftools query -f '%CHROM\t%POS\t%QUAL\t%MQ\n' ${VCF_FILE}  > Qtn5#Bcftools_output.txt
echo -e '\e[1;34m'
echo -e 'Code Used: \n bcftools query -f '%CHROM\t%POS\t%QUAL\t%MQ\n' ${VCF_FILE}  > Qtn5#Bcftools_output.txt'


echo -e "\e[1;33m6.Extract data that belongs to chromosomes 2,4 and MT."
#To Compress the vcf file
bcftools view ${VCF_FILE}  -Oz -o ${VCF_FILE}.gz   #To Compress the vcf file
#Then indexed thus to be able to filter out the required data of specific chromosomes. 
bcftools index ${VCF_FILE}.gz #To index the compressed vcf file 
bcftools view -r 2,4,MT ${VCF_FILE}.gz > Qtn6#Chrom_2_4_MT_output.vcf  #To filter  data
echo -e '\e[1;34m'
echo -e 'Code Used: \n bcftools view ${VCF_FILE}  -Oz -o ${VCF_FILE}.gz \n bcftools index ${VCF_FILE} .gz \n bcftools view -r 2,4,MT ${VCF_FILE}.gz > Qtn6#Chrom_2_4_MT_output.vcf'


echo -e "\e[1;33m7.Print out variants that do not belong to chr20:1-30000000."
bcftools view -R ^20:1-30000000 ${VCF_FILE}  > Qtn7#output.vcf
echo -e '\e[1;34m'
echo -e 'Code Used: \n bcftools view -R ^20:1-30000000 ${VCF_FILE}  > Qtn7#output.vcf'


echo -e "\e[1;33m8.Extract variants that belong to SRR13107019."
bcftools view -s SRR13107019 ${VCF_FILE}  > Qtn8#SRR13107019_output.vcf 
bcftools query -l Qtn8#SRR13107019_output.vcf  #To verify 
echo -e '\e[1;34m'
echo -e 'Code Used: \n bcftools view -s SRR13107019 ${VCF_FILE}  > Qtn8#SRR13107019_output.vcf  \n bcftools query -l SRR13107019_output.vcf #To verify'

echo -e "\e[1;33m9.Filter out variants with a QualByDepth above 7."
bcftools filter -i 'QUAL/FORMAT/DP > 7' ${VCF_FILE}  > Qtn9#Qua_depth_FORMAT_output.vcf 
bcftools filter -i 'QUAL/INFO/DP > 7' ${VCF_FILE}  > Qtn9#Qual_depth_INFO_output.vcf
echo -e '\e[1;34m'
echo -e 'Code Used: \n bcftools filter -i 'QUAL/FORMAT/DP > 7' ${VCF_FILE}  > Qtn9#Qua_depth_FORMAT_output.vcf /n bcftools filter -i 'QUAL/INFO/DP > 7' ${VCF_FILE}  > Qtn9#Qual_depth_INFO_output.vcf'

echo -e "\e[1;33m10.How many contigs are referred to in the file. Check the header section"
bcftools view -h ${VCF_FILE}  | grep -c "^##contig" > Qtn10#contig.txt
echo -e '\e[1;34m'
echo -e 'Code Used: \n bcftools view -h ${VCF_FILE}  | grep -c "^##contig"' > Qtn10#contig.txt


echo -e "\e[1;33m12.Extract data on the read depth of called variants for sample SRR13107018"
bcftools view -s SRR13107018 ${VCF_FILE}  | bcftools query -f '%DP\n' > Qnt12#SRR13107018_output.txt
echo -e '\e[1;34m'
echo -e 'Code Used: \n bcftools view -s SRR13107018 ${VCF_FILE}  | bcftools query -f '%DP\n' > Qnt12#SRR13107018_output.txt'

echo -e "\e[1;33m13.Extract data on the allele frequency of alternate alleles. Combine this data with the chromosome and position of the alternate allele"
bcftools query -f '%CHROM\t%POS\t%AF\n' ${VCF_FILE}  > Qnt13#Chrom_POS_AF_output.txt
echo -e '\e[1;34m'
echo -e 'Code Used: \n bcftools query -f '%CHROM\t%POS\t%AF\n' ${VCF_FILE}  > Qnt13#Chrom_POS_AF_output.txt'



echo -e "\e[1;5;35m\n\n\n ##Manipulating SAM files##\e[0m\n\n\n"                               


echo -e "\e[1;33m3.How many samples are in the file."
grep '@RG' ${SAM_FILE} | awk '{print $2}' | awk -F':' '{print $2}'| sort | uniq | wc -l
echo -e '\e[1;34m'
echo -e 'Code Used: \n grep '@RG' ${SAM_FILE} | awk '{print $2}' | awk -F':' '{print $2}'| sort | uniq | wc -l'

echo -e "\e[1;33m4. How many alignments are in the file."
samtools view -S -b ${SAM_FILE} > sample.bam
samtools flagstat sample.bam
echo -e '\e[1;34m'
echo -e 'Code Used: \n samtools view -S -b ${SAM_FILE} > sample.bam \n samtools flagstat sampe.bam > Qn4.bam_ouput.txt'

echo -e "\e[1;33m5.Get summary statistics for the alignments in the file."
samtools stats sample.bam > Qtn#5Samtools_sample.stats
samtools depth sample.bam > Qtn#5Samtools_sample.depth
samtools view -c -F4 sample.bam
echo -e '\e[1;34m'
echo -e 'Code Used: \n samtools stats sample.bam > sample.stats \n samtools depth sample.bam > sample.depth \n samtools view -c -F4 sample.bam'

echo -e "\e[1;33m6.Count the number of fields in the file."
cut -f1- ${SAM_FILE} | awk '{print NF;}' | sort -nu | tail -n1 
echo -e '\e[1;34m'
echo -e 'Code Used: \n cut -f1- ${SAM_FILE} | awk '{print NF;}' | sort -nu | tail -n1'

echo -e "\e[1;33m7.Print all lines in the file that have @SQ and sequence name tag beginning with NT_"
grep '^@SQ.*NT_' ${SAM_FILE} > Qtn7#Samtools_lines_@SQsequence.txt
echo -e '\e[1;34m'
echo -e 'Code Used: \n grep '^@SQ.*NT_' ${SAM_FILE} > Qtn7#Samtools_lines_@SQsequence.txt'

echo -e "\e[1;33m8.Print all lines in the file that have @RG and LB tag beginning with Solex"
awk '$1=="@RG" && $3=="LB:Solexa"' ${SAM_FILE}
echo -e '\e[1;34m'
echo -e 'Code Used: \n awk '$1=="@RG" && $3=="LB:Solexa"' ${SAM_FILE}'

echo -e "\e[1;33m9.Extract primarily aligned sequences and save them in another file"
grep -v 'XA:Z:' ${SAM_FILE} > Qtn9#Samtools_primary_alignments.sam
echo -e '\e[1;34m'
echo -e 'Code Used: \n grep -v 'XA:Z:' ${SAM_FILE} > Qtn9#primary_alignments.sam'


echo -e "\e[1;33m10.Extract alignments that map to chromosomes 1 and 3. Save the output in BAM format"
samtools view -S -b ${SAM_FILE} > sample.bam
samtools view -b sample.bam chr1 chr3 > Qtn10#Samtoolschromosomes_1_3.bam
echo -e '\e[1;34m'
echo -e 'Code Used: \n samtools view -S -b ${SAM_FILE} > sample.bam \n samtools view -b sample.bam chr1 chr3 > Qtn10#Samtoolschromosomes_1_3.bam'


echo -e "\e[1;33m11.How would you obtain unmapped reads from the file"
samtools view -b -f 4 sample.bam > Qtn11Samtools#unmapped_reads.bam
echo -e '\e[1;34m'
echo -e 'Code Used: \n samtools view -b -f 4 sample.bam > Qtn11#Samtoolsunmapped_reads.bam'


echo -e "\e[1;33m12.How many reads are aligned to chromosome 4"
samtools index sample.bam
samtools view -c -F 4 sample.bam 
echo -e '\e[1;34m'
echo -e 'Code Used: \n samtools index sample.bam \n samtools view -c -F 4 sample.bam'


echo -e "\e[1;33m13.Comment of the second and sixth column of the file"
awk '{print $2 "\t" $6}' ${SAM_FILE} > Qnt14Samtools#Sec_Six.txt
awk '{ for (i = 0; i < length($6); i++) { if (substr($6, i, 1) == "M") { matches++ } if (substr($6, i, 1) == "I") { insertions++ } if (substr($6, i, 1) == "D") { deletions++ } } } END { print "matches: " matches "\ninsertions: " insertions "\ndeletions: " deletions }' ${SAM_FILE}
echo -e '\e[1;34m'
echo -e 'Code Used: \n awk '{print $2 "\t" $6}' ${SAM_FILE} \n awk '{ for (i = 0; i < length($6); i++) { if (substr($6, i, 1) == "M") { matches++ } if (substr($6, i, 1) == "I") { insertions++ } if (substr($6, i, 1) == "D") { deletions++ } } } END { print "matches: " matches "\ninsertions: " insertions "\ndeletions: " deletions }' ${SAM_FILE}'


echo -e "\e[1;33m14.Extract all optional fields of the file and save them in “optional_fields.txt”"
awk '{for(i=12;i<=NF;i++) printf("%s\t",$i); printf("\n")}' ${SAM_FILE} > Qnt14Samtools_optional_fields.txt
echo -e '\e[1;34m'
echo -e 'Code Used: \n awk '{for(i=12;i<=NF;i++) printf("%s\t",$i); printf("\n")}' ${SAM_FILE} > Qnt14Samtools_optional_fields.txt'