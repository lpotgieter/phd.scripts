for f1 in /media/lizel/82AEF34EAEF3396D/chapter2/2_vcf/45_pop/*GT.FORMAT
do
#vcftools --vcf $VCFNAME.vcf --extract-FORMAT-info GT
cat ${f1} | cut -f 3- > ${f1}.2
head -n 1 ${f1}.2 | sed "s/\t/\n/g" > ${f1}.ind
sed "s/\t//g" ${f1}.2 | tail -n +2 > ${f1}.geno
sone
