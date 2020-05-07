# pylae
Python for local ancestry estimation

## Requirements and installation:

* Python 3.5+ is required
* bcftools
* (optionally) plink / plink2

Installing python requirements:
```bash
pip3 install -r requirements.txt
```
## Usage:
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
Currently supported mode: bayes.

```bash
python3 src/bayesian_pipeline.py --sample <sample_name> --admixtures <admixture_vectors_file> --window-len 50 <group>.<sample>.txt
```

### Example pipeline:
```bash
plink2 --bfile America.QuechuaCandelaria_3.txt_GENO --recode vcf --out America.QuechuaCandelaria_3_GENO

cat America.QuechuaCandelaria_3_GENO.vcf | bcftools view -c 1 -Ou | bcftools +fill-tags -Ou -- -S vcf_groups.txt -t AF | bcftools query -H -f "%CHROM %POS %ID %AF_QuechuaCandelaria_3 %AF_Mediterranean %AF_NativeAmerican %AF_NorthEastAsian %AF_NorthernEuropean %AF_Oceanian %AF_SouthAfrican %AF_SouthEastAsian %AF_SouthWestAsian %AF_SubsaharanAfrican\n" > "QuechuaCandelaria_3.GA002786.txt"

python3 src/bayesian_pipeline.py --window-len 50  "QuechuaCandelaria_3.GA002786.txt"
```


### Files explanation
As a result of the pipeline we get 3 files:
1. `<group>_<mode>_<window-len>_predictions.csv`   
Csv file with a list of most probable population in each window.

2. `<group>_<mode>_<window-len>_snp_prob.tsv`   
Tsv (tab-separated) file with a list of all SNPs and probabilities that it came from each population.

3. `<group>_<mode>_<window-len>_stats.csv`  
Csv file with statistics that shows the fraction of windows assigned to each population.

Depending on your needs you might need only one file or all of them.

## Algorithm explanation
Algorithm can be split into 4 stages:  
* Data preparation 
* Calculating probabilities of assigning each SNP to populations using naive bayes algorithm.  
* Choosing best population for each window with selected length (in SNPs).  
In this slog (p). Then this information (I) is summed in each window and the window 
is assigned to population with max I. Pop = argmax(I)
* Calculating fraction of windows assigned to each population.




## Modes explanation
### 1. Bayes
Probability of assigning snp to population is calculated according to the Bayes formula:  
<img src="https://render.githubusercontent.com/render/math?math=P(Population | SNP) = \frac{P(SNP | Population) \cdot P(Population)}{P(SNP)}">  
Here,   
<img src="https://render.githubusercontent.com/render/math?math=P(SNP | Population)"> can be estimated as frequency of SNP in selected Population.  
<img src="https://render.githubusercontent.com/render/math?math=P(Population) = \frac{1}{|Populations|}"> - we take prior population probabilities from admixture vectors.  
<img src="https://render.githubusercontent.com/render/math?math=P(SNP)"> can be estimated as average frequency of SNP among all populations or samples.  

