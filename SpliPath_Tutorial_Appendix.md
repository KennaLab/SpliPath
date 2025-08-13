
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

### 3. Predict splice altering effect of DNA variants using SpliceAI and Pangolin.

#### SpliceAI

There are pre-computed SpliceAI scores (https://github.com/Illumina/SpliceAI) in VCF format. Or, users can run SpliceAI with customized parameters.
The VCF files contain splice-altering predictions or any other annotations should be first converted to a table (e.g. by function [```rvat::vcfInfo2Table```](https://kennalab.github.io/rvat/reference/vcfInfo2Table.html)), where the rows are the variants and columns are CHROM, POS, ID, REF, ALT, and any other annotations (For SpliceAI, the columns are "spliceaiDS_AG", "spliceaiDS_AL", "spliceaiDS_DG", "spliceaiDS_DL", "spliceaiDP_AG", "spliceaiDP_AL", "spliceaiDP_DG", "spliceaiDP_DL", and etc.). 
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

In the output VCF, there might be multiple predictions for one variant in a line. The following code can be applied to convert the SpliceAI output VCF to table
NB: This code assumes the ID column in the VCF is VAR_id from genome GDB
```{sh}
(echo -e "VAR_id\tCHROM\tPOS\tREF\tALT\tspliceaiSYMBOL\tspliceaiDS_AG\tspliceaiDS_AL\tspliceaiDS_DG\tspliceaiDS_DL\tspliceaiDP_AG\tspliceaiDP_AL\tspliceaiDP_DG\tspliceaiDP_DL";
cat example_rvatData_SpliceAI_prediction.vcf | grep -v '^#' | grep -v ',' | awk '$8!="."' | awk -F'[|\t]' '{print $3, $1, $2, $4, $5, $9, $10, $11, $12, $13, $14, $15, $16, $17}' OFS='\t' | awk '$10!="."' ;
cat example_rvatData_SpliceAI_prediction.vcf | grep -v '^#' | grep ',' | awk -F'[,\t]' '{print $3, $1, $2, $4, $5, $8} {for (i = 9; i <= 100; i++) print $3, $1, $2, $4, $5, "SpliceAI="$i }' OFS="\t" | awk '$6!="SpliceAI="' | awk -F'[|\t]' '{print $1, $2, $3, $4, $5, $7, $8, $9, $10, $11, $12, $13, $14, $15}' OFS='\t' | awk '$10!="."' 
) > example_rvatData_SpliceAI_prediction.txt
```

#### Pangolin

Run Pangolin:
```{sh}
pangolin -m False -d 5000 example_rvatData.vcf hs38DH.fa gencode.v32.chr_patch_hapl_scaff.basic.annotation.db example_rvatData_pangolin_prediction
```
See https://github.com/tkzeng/Pangolin for input genome annotation files.

Convert output VCF to table
NB: This code assumes the ID column in the VCF is VAR_id from genome GDB
```{sh}
(echo -e "VAR_id\tCHROM\tPOS\tREF\tALT\tgene_id\tpos_increase\tscore_increase\tpos_decrease\tscore_decrease"
grep -v '^#' example_rvatData_pangolin_prediction.vcf | awk '$8!="."' | awk -F'[|:=\t]' '{print $3, $1, $2, $4, $5, $9, $10, $11, $12, $13}' OFS='\t' ) > example_rvatData_pangolin_prediction.txt
```

### 4. Extract annotated introns from genome annotation file.
The GTFtools can be downloaded in https://www.genemine.org/gtftools.php
```{sh}
python gtftools.py -i Intron_GRCh38.98.bed Homo_sapiens.GRCh38.98.gtf  
```