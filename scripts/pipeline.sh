##make file subset.sh
#!/bin/bash
echo 'Start'
tabix -h zombies.vcf.gz

for chr in {1..22}
do
    echo "Processing chromosome $i"
    bgzip -c "N1511_fil_ch$chr.phased.vcf" > "N1511_fil_ch$chr.phased.vcf.gz"
    tabix -h "N1511_fil_ch$chr.phased.vcf.gz"

    bcftools isec -Ob zombies.vcf.gz "N1511_fil_ch$chr.phased.vcf.gz" -p "intersected_ch$chr"
    tabix -h "intersected_ch$chr"/0002.bcf.gz
    tabix -h "intersected_ch$chr"/0003.bcf.gz
    bcftools view -s "GA002786,GA003530,GA002674,GA000130,GA002973" -Ou "intersected_ch$chr"/0003.bcf.gz | bcftools merge "intersected_ch$chr"/0002.bcf.gz -Ob > "intersected_ch$chr"/merged.bcf.gz
    tabix -h "intersected_ch$chr"/merged.bcf.gz

done

bcftools concat -n intersected_ch1/merged.bcf.gz intersected_ch2/merged.bcf.gz intersected_ch3/merged.bcf.gz intersected_ch4/merged.bcf.gz intersected_ch5/merged.bcf.gz intersected_ch6/merged.bcf.gz intersected_ch7/merged.bcf.gz intersected_ch8/merged.bcf.gz intersected_ch9/merged.bcf.gz intersected_ch10/merged.bcf.gz intersected_ch11/merged.bcf.gz intersected_ch12/merged.bcf.gz intersected_ch13/merged.bcf.gz intersected_ch14/merged.bcf.gz intersected_ch15/merged.bcf.gz intersected_ch16/merged.bcf.gz intersected_ch17/merged.bcf.gz intersected_ch18/merged.bcf.gz intersected_ch19/merged.bcf.gz intersected_ch20/merged.bcf.gz intersected_ch21/merged.bcf.gz intersected_ch22/merged.bcf.gz -Ob > full.bcf.gz