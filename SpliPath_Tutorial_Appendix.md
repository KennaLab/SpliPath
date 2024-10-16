
##  SpliPath Tutorial Appendix

### 1. Use RegTools to generate splice junctions from RNAseq alignment for input to SpliPath.

#### Install [RegTools](https://regtools.readthedocs.io/en/latest/):
```{sh}
git clone https://github.com/griffithlab/regtools
cd regtools/
mkdir build
cd build/
cmake ..
make
```

#### Command to extract junction data from split-reads with minimum 6bp anchors in aligned RNAseq BAM files:
```{sh}
regtools junctions extract -a 6 -s 1 -o sample.junc sample.bam
```
Option ```-s``` in ```regtools junctions extract``` means strand specificity of RNA library preparation, please change it according to the input RNAseq samples.

The output ```sample.junc``` (e.g. in folder ```SpliPath/Example_data/regtools_junction/``` ) files will be the input of both LeafCutterMD splicing outlier analysis and SpliPath.


### 2. Run LeafCutterMD to identify sample with outlier splicing signal for each intron.

#### Install [LeafCutter](https://davidaknowles.github.io/leafcutter/index.html):
```{sh}
git clone https://github.com/davidaknowles/leafcutter
```
```{r}
install.packages("rstan")
devtools::install_github("stan-dev/rstantools")
```
If encountered problem in installing, see [LeafCutter](https://davidaknowles.github.io/leafcutter/index.html) documentation for extended information.

#### Run LeafCutterMD splicing outlier analysis for samples of each tissue.
```{sh}
tissues=($(cat example_RNAseq_meta.txt | tail -n+2 | cut -f4 | sort | uniq ))
for t in ${tissues[@]}; do
  # Make a list of BED file paths for each tissue
  awk -F"\t" -v t=${t} '$4==t {print $6}' example_RNAseq_meta.txt > leafcutter/sample_${t}_path.list
  # Clustering the junctions
  python2 leafcutter_cluster_regtools.py -r leafcutter/ -j leafcutter/sample_${t}_path.list -o example_LeafCutter_${t} -s True --maxintronlen 500000 -m 10 -M 1 -p 0.000001
  # Splicing junction outlier analysis
  Rscript leafcutterMD.R --num_threads 4 -c 10 -s 500 -o leafcutter/example_LeafCutter_${t}_outlier leafcutter/example_LeafCutter_${t}_perind_numers.counts.gz
done
```

### 3. Predict splice altering effect of DNA variants using SpliceAI and dbscSNV.

There are pre-computed SpliceAI scores (https://github.com/Illumina/SpliceAI), and pre-computed dbscSNV scores (http://www.liulab.science/dbscsnv.html) in VCF format. Or, users can run SpliceAI with customized parameters.
The VCF files contain splice-altering predictions or any other annotations should be first transformed to a table (e.g. by function [```rvat::vcfInfo2Table```](https://kennalab.github.io/rvat/reference/vcfInfo2Table.html)), where the rows are the variants and columns are CHROM, POS, ID, REF, ALT, and any other annotations (For SpliceAI, the columns are "spliceaiDS_AG", "spliceaiDS_AL", "spliceaiDS_DG", "spliceaiDS_DL", "spliceaiDP_AG", "spliceaiDP_AL", "spliceaiDP_DG", "spliceaiDP_DL", and etc. For dbscSNV, the columns are "ada_score", "rf_score", and etc.). 
For example,
```{r}
### R session
library(rvat)
vcfInfo2Table(vcf = "SpliceAI_precomputed_score.vcf", output = "SpliceAI_precomputed_score.txt", splitMultiallelic = TRUE)
```

Or, run SpliceAI using customized gene annotation and parameters, for example:
```{sh}
pip install spliceai
# or
conda install -c bioconda spliceai

spliceai \
-I example_rvatData.vcf \
-O example_rvatData_SpliceAI_prediction.vcf \
-D 500 \
-R Homo_sapiens.GRCh38.dna.primary_assembly.fa \
-A Ensembl_GRCh38.98.basic.transcript.txt
```

### 4. Extract annotated introns from genome annotation file.
The GTFtools can be downloaded in https://www.genemine.org/gtftools.php
```{sh}
python gtftools.py -i Intron_GRCh38.98.bed Homo_sapiens.GRCh38.98.gtf  
```