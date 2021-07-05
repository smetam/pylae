# PyLAE
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
Currently supported modes: 
#### Bayes:

```bash
python3 src/bayesian_pipeline.py --sample <sample_name> --admixtures <admixture_vectors_file> --window-len 50 <group>.<sample>.txt
```
#### Bayes viterbi (used in the paper):

```
python3 src/bayes_viterbi.py --sample <sample_name> --admixtures <admixture_vectors_file> --window-len 50 <group>.<sample>.txt -m 

```

`-m` option is used to switch "merged" window mode (windows will overlap by 1 SNP)


### Example pipeline:
```bash
plink2 --bfile sample.txt_GENO --recode vcf --out sample

cat sample.vcf | bcftools view -c 1 -Ou | bcftools +fill-tags -Ou -- -S vcf_groups.txt -t AF | bcftools query -H -f "%CHROM %POS %ID %AF_QuechuaCandelaria_3 %AF_Mediterranean %AF_NativeAmerican %AF_NorthEastAsian %AF_NorthernEuropean %AF_Oceanian %AF_SouthAfrican %AF_SouthEastAsian %AF_SouthWestAsian %AF_SubsaharanAfrican\n" > "population.sample.txt"

python3 src/bayesian_pipeline.py --window-len 50  "population.sample.txt"
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

Using PyLAE with different genomes and/or sets of markers
A different set of putative ancestral populations or a different set of markers require 
additional processing. First, we need to collect a database of putatively un-admixed individuals. 
If there is an existing validated set of ancestry informative features, these markers should run the 
admixture in supervised mode. For each self-reported ancestry, samples should be clustered 
based on their admixture profiles to identify subgroups within each self-reported ancestry. These 
subgroups are then examined using information about the studied population's history, and the 
most representative subset is retained. Then, putative ancestral populations (from 15 to 20 
individuals per group) are generated for every ancestry. The validity and stability of the ancestral 
populations are evaluated using 1) PCA, 2) leave-one-out supervised admixture, and 3) by 
application of supervised admixture to the original datase

Algorithm can be split into 4 stages:  
* Data preparation 
* Calculating probabilities of assigning each SNP to populations using naive bayes algorithm.  
* Choosing best population for each window with selected length (in SNPs).  
In this slog (p). Then this information (I) is summed in each window and the window 
is assigned to population with max I. Pop = argmax(I)
* Calculating fraction of windows assigned to each population.


## Preprint about the method
https://www.biorxiv.org/content/10.1101/2020.11.13.380105v1
