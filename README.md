# pylae
Python for local ancestry estimation

### Data preparation stage:
(will be performed by script itself in future)

1. In case we have .bed .bim .fam files, we need to convert to vcf using plink:
```bash
plink2 --bfile <bfile_prefix> --recode vcf --out <vcf_file>
```

2. Calculate snp frequencies for population groups using bcftools. 
User groups defined in file `configs/vcf_groups.txt`:

```bash
cat <vcf_file> | bcftools view -c 1 -Ou | bcftools +fill-tags -Ou -- -S configs/vcf_groups.txt -t AF | bcftools query -H -f "%CHROM %POS %ID %AF_<group> %AF_Mediterranean %AF_NativeAmerican %AF_NorthEastAsian %AF_NorthernEuropean %AF_Oceanian %AF_SouthAfrican %AF_SouthEastAsian %AF_SouthWestAsian %AF_SubsaharanAfrican\n" > <group>.<sample>.txt
```

In case vcf file is (b)gzipped use samtools tabix.

### Script usage:
Currently supported modes: bayes, fb, softmax.
Note: fb is around 20 times slower.

```bash
python3 src/process_individuals.py --mode fb --window-len 200  <group>.<sample>.txt
```

Example pipeline:
```bash
plink2 --bfile America.QuechuaCandelaria_3.txt_GENO --recode vcf --out America.QuechuaCandelaria_3_GENO

cat America.QuechuaCandelaria_3_GENO.vcf | bcftools view -c 1 -Ou | bcftools +fill-tags -Ou -- -S vcf_groups.txt -t AF | bcftools query -H -f "%CHROM %POS %ID %AF_QuechuaCandelaria_3 %AF_Mediterranean %AF_NativeAmerican %AF_NorthEastAsian %AF_NorthernEuropean %AF_Oceanian %AF_SouthAfrican %AF_SouthEastAsian %AF_SouthWestAsian %AF_SubsaharanAfrican\n" > "QuechuaCandelaria_3.GA002786.txt"

python3 src/process_individuals.py --mode fb --window-len 200  "QuechuaCandelaria_3.GA002786.txt"
```

## Estimated performance:
for vcf file with around 120k SNPs.
|mode|exec time, min| ? |  
|--|--|--|  
|fb   | 20 | ? |  
|bayes| 1  | ? |
|softmax| 0.5 | ? |



### Files explanation
As a result of the pipeline we get 3 files:
1. `<group>_<mode>_<window-len>_predictions.csv`   
Csv file with a list of most probable population in each window.

2. `<group>_<mode>_<window-len>_snp_prob.tsv`   
Tsv (tab-separated) file with a list of all SNPs and probabilities that it came from each population.

3. `<group>_<mode>_<window-len>_stats.csv`  
Csv file with statistics that shows the fraction of windows assigned to each population.
