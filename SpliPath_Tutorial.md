#  Splicing Pathology (SpliPath)

## Introduction

SpliPath is designed to identify and functionally cluster splice-altering variants in WGS that have similar effects on RNA splicing. We refer to these variants as collapsed splicing quantitative trait loci (csQTLs). SpliPath aims to address two main difficulties in explaining missing heritability in rare disorders: effective method to functionally interpret genetic variants and increase statistical power to establish associations between rare variants and phenotypes. First, SpliPath links the prediction of SpliceAI with reference transcriptomics data to identify genetic variants that induce splice changes actually occuring in disease-relevant transcriptomics profiles. Second, SpliPath aggregates variants with similar functional consequences into csQTLs for more powerful genetic association analyses.

## Installation

### Clone repo (recommended)

Clone repo
```{sh}
git clone https://github.com/KennaLab/SpliPath.git
```
Open an R session and install using devtools 
```{r}
devtools::install("SpliPath") 
# if the SpliPath directory is not in your current directory,  provide the full path to the cloned repo
```

### Install dependencies

* [RVAT](https://github.com/KennaLab/rvat)
```{r}
remotes::install_github("kennalab/rvat")
```

## Tutorial
This tutorial shows how to nominate csQTL using either paired and unpaired WGS and RNAseq datasets. First, we annotate splicing junctions outliers derived from LeafCutter outlier analysis, or, a simple list of query splicing junction. Then, we putatively link them with genetic variants based on SpliceAI prediction of their effect on altering splice sites. Finally, we show how to use SpliPath data browser to visualize and manually inspect candidate csQTL. Considering the size of human genome and transcriptome, we recommend to perform the data analyses (Step 1-4.1) using a high performance computing cluster, while data visualization (Step 4.2) can be performed on any standard laptop. Considering the size of human genome and transcriptome, we recommend to perform the data analyses (Step 1-4.1) using a high performance computing cluster, while data visualization (Step 4.2) can be performed on any standard laptop.

* The example input and output files can be found in ```SpliPath/Example_data/```.

### Step 1. Data preparation

To nominate csQTL, SpliPath takes 1. reference splicing junctions and 2. genomic variants as inputs:

#### 1. Splicing junctions

If the reference splicing junctions to nominate csQTL are from LeafCutterMD outlier analysis, the following data should be available:
1. LeafCutter splicing outlier analysis output (e.g. ```SpliPath/Example_data/example_LeafCutter_outlier_pVals.txt``` )
(Two separated meta files are need for DNA & RNA sequencing samples is because there can be multiple RNAseq samples for a subject. SpliPath match the paired samples by 'SubjectID'.)
2. Splicing junction data in BED format (e.g. ```SpliPath/Example_data/regtools_junction/*``` )
3. A tab-delimited metadata table of the RNA sequencing samples. It should at least contain 'SubjectID', 'SampleID', 'Group', 'Tissue' columns and a 'Path' column which specify the paths to the splicing junction BED files of each sample (e.g. ```SpliPath/Example_data/example_RNAseq_meta.txt``` )

Alternatively, a simple list of splicing junctions can also be used to nominate csQTL:
1. A BED file containing the coordinates of the splicing junctions (e.g. ```SpliPath/Example_data/example_query_junction.bed``` ). 

#### 2. Genomic variants
1. Genomic variants in VCF format (e.g.  ```SpliPath/Example_data/example_rvatData.vcf``` ). It should includes all samples in the metadata table.
2. SpliceAI prediction (e.g. ```SpliPath/Example_data/example_chr1_varAnno_SpliceAI.txt``` & ```SpliPath/Example_data/example_chr1_varAnno_dbscSNV.txt```. There are pre-computed SpliceAI scores (https://github.com/Illumina/SpliceAI), and pre-computed dbscSNV scores (http://www.liulab.science/dbscsnv.html) in VCF format.)
3. A tab-delimited metadata table of the DNA sequencing samples. It should at least contain 'SubjectID', 'SampleID', and 'Group' (e.g. ```SpliPath/Example_data/example_DNA_meta.txt``` ). 

See [Appendix](https://github.com/KennaLab/SpliPath/blob/main/SpliPath_Tutorial_Appendix.md) for the commands of running RegTools, LeafCutter and SpliceAI for generating the input files above.

We use [RVAT](https://github.com/kkenna/rvat) to store variant annotations, individual genotypes and other required data in SQLite based “.gdb” files. This format has been optimized to support efficient execution of SpliPath functions as well as a range of other rare variant analyses.
```{r}
library(SpliPath)
setwd("SpliPath/Example_data/")

vcfpath <- "example_rvatData.vcf"
vcfpheno <- "example_DNA_meta.txt"
gdbpath <- "example_rvatData.gdb"
buildGdb(vcf = vcfpath,  output = gdbpath, overWrite = TRUE)

# Upload sample information
mygdb = gdb(gdbpath)
cohort_pheno <- read.table(vcfpheno, header = TRUE, sep = "\t")
cohort_pheno$IID = cohort_pheno$SampleID
uploadCohort(object = mygdb, name = "pheno", value = cohort_pheno) # The 'name' parameter is the name of the cohort table in the GDB
```

The variants annotation tables (SpliceAI predictions) can be uploaded into GDB by 
```{r}
spliceai_pred_path <- "example_chr1_varAnno_SpliceAI.txt"
spliceai_pred <- read.table(spliceai_pred_path, header = TRUE, sep = "\t")
uploadAnno(object = mygdb, name = "SpliceAI",value = spliceai_pred) # The 'name' parameter is the name of the variants annotation table in the GDB
```
Users can also upload other variants annotation to GDB and include them into analysis in Step 3.

#### (Optional) Customize reference files
The default reference transcripts annotations are from Ensembl GRCh38.98 GTF file and Snaptron annotation, which have been installed in SpliPath.

Users can create customized reference file using ```prepareGenoRef``` function:
```{r}
prepareGenoRef(gtf_file = "Homo_sapiens.GRCh38.98.gtf.gz", intron_bed = "Intron_GRCh38.98.sort.bed.gz", output_dir = "new_reference")
```
The argument ```intron_bed``` of ```prepareGenoRef``` function is a BED file contain coordinates (chromosome, start position, end position, strand), transcript id, and gene id information of introns. This file can be generated by [GTFtools] (https://www.genemine.org/gtftools.php) from GTF file, see [Appendix 4](https://github.com/KennaLab/SpliPath/blob/main/SpliPath_Tutorial_Appendix.md)

The output reference files will be written into the provided ```output_dir``` directory. The four output files are 
1) "Gene.bed": coordinates of genes; 
2) "Exon_proteincoding.bed": exons in protein coding transcripts; 
3) "Intron_proteincoding.bed": introns in protein coding transcripts; 
4) "Intron_noncoding.bed" introns in noncoding transcripts.

Or, users can generate reference files following the file format of the default reference files in ```SpliPath/inst/extdata/Reference``` directory. And name them with "Gene.bed", "Exon_proteincoding.bed", "Intron_proteincoding.bed" and "Intron_noncoding.bed".

After generating reference files, specify the directory to the generated reference files in ```reference``` argument in functions of Step 2-4. 

### Step 2. Annotate splicing junctions

Considering computational capacity, we highly recommend performing Step 2-3 per chromosome.

After data preparation, SpliPath starts with annotating the splicing junctions. The annotation include mapping the junctions to genes, finding unannotated splice sites and classifying splice events. Excluding possible gene fusion events and artefacts where the two splice sites of a junction are mapped two genes or neither of them are annotated. 

If the splicing junctions are from LeafCutterMD outlier analysis, function ```annotateLeafCutterJunc``` is applied:
* ```tissues_leafcutter_pvals_file``` lists paths to the LeafCutter analysis P values of each tissue type. Each element is a LeafCutter file path named under a tissue. The tissues should be the same with those in the RNAseq sample metadata file. e.g.
```{r}
tissues_leafcutter_pvals_file <- list(motor_cortex = "leafcutter/example_LeafCutter_motor_cortex_outlier_pVals.txt", 
                                     frontal_cortex = "leafcutter/example_LeafCutter_frontal_cortex_outlier_pVals.txt", 
                                     cervical = "leafcutter/example_LeafCutter_cervical_outlier_pVals.txt", 
                                     lumbar = "leafcutter/example_LeafCutter_lumbar_outlier_pVals.txt")

annotateLeafCutterJunc(rna_meta = "example_RNAseq_meta.txt", 
                       chrom = 1, 
                       tissues_leafcutter_pvals_file,
                       max_LeafCutter_Pval = 0.05,
                       min_nr_tissue = 1,
                       output_prefix = "example_chr1", 
                       reference = "Default")
```

The function outputs three files: 
1. example_chr1.txt.gz contains junction read count table merged from all the samples. The unique junctions are in rows and samples are in columns; 
2. example_chr1_anno.txt.gz contains the annotation of the junctions in the LeafCutter files. 
3. example_chr1_leafcutter_outliers.txt.gz contains the annotation of the unannotated junction outliers, filtered by the provided max_LeafCutter_Pval and min_nr_tissue thresholds.

Alternatively, if a simple list of splicing junctions will be used for further csQTL annotation, ```annotateQueryJunc``` is applied.
```{r}
annotateQueryJunc(query_junc = "example_query_junction.bed", 
                  output_prefix = "example_chr1", 
                  reference = "Default")
```

The function outputs one file: example_chr1_query_junc_anno.txt.gz, which contains the annotation of the query splicing junctions 

### Step 3. Nominate csQTL

The csQTLs are annotated by function ```collapsedsQTL```. Different modes will be used for paired and unpaired WGS - RNAseq data analyses.

The arguments are 
1. Junction and variants inputs:
* ```junc_anno_file``` is the file containing the annotations of LeafCutter outlier junctions or query junction, generated in Step 2.
* ```tissues_leafcutter_pvals_file``` lists paths to the LeafCutter analysis P values of each tissue type. Each element is a LeafCutter file path named under a tissue. The tissues should be the same with those in the RNAseq sample metadata file. e.g.
* ```gdb_path``` is the GDB file of DNA variants generated in Step 1. The table of phenotype data should be specified in ```cohort_name```.
* ```rna_meta``` and ```wgs_meta``` is the file of phenotype / metadata of the RNA and DNA sequencing samples. Both of them should have 'SubjectID', 'SampleID', and 'Group' to specify which individual the sample is collected and which phenotype group the individuals belong to. In addition, 'rna_meta' should have 'Tissue' column, which specify the tissue source of the RNAseq sample. The two files are not required for unpaired data analyses 

2. Thresholds to nominates csQTLs
* ```min_nr_tissue``` requires LeafCutter P value of a junction should be under the threshold in at least this number of tissues
* ```max_allele_freq``` is the maximum allele frequency to nominate ultra rare sQTL candidate 
* ```splice_prediction``` is the splice-altering prediction to annotate the DNA variants. The names of elements in the list are prediction tools, they should be in GDB database tables. The items are the scoring fields in the GDB database prediction tables A splice prediction tool may provide more then one scores. The maximum of the given scores will be used as the prediction of the corresponding tool for each variants. Default to:
* ```min_splice_score``` is the min splice score of each prediction tools to nominate ur-sQTL candidate.
* ```spliceai_default_reference``` suggests whether SpliceAI predictions are pre-computed or SpliceAI default annotation were used when running SpliceAI

#### Paired WGS and RNAseq analyses

For paired WGS and RNAseq analysis, the argument ```junc_anno_file``` takes the output of ```annotateLeafCutterJunc``` function.

```{r}
collapsedsQTL(paired_data = TRUE,
              junc_anno_file = "example_chr1_leafcutter_outliers.txt.gz",
              gdb_path = "example_rvatData.gdb",
              cohort_name = "pheno",
              rna_meta = "example_RNAseq_meta.txt",
              wgs_meta = "example_DNA_meta.txt",
              max_allele_freq = 0.05,
              reference = "Default",
              spliceai_prediction=list(SpliceAI = c("spliceaiDS_AG", "spliceaiDS_AL", "spliceaiDS_DG", "spliceaiDS_DL")),
              min_spliceai_score = 0.2,
              spliceai_default_reference = T,
              output_prefix = "example_chr1")
```
The nominated csQTLs are written in example_chr1_csQTL_paired_csQTL.txt.gz file.
One row in the the data.frame is an genetic variant that mapped to a unannotated junction in a subject's tissues. The LeafCutter P values and read count in tissues are in columns "Outlier_P_tissue" and "Reads_tissue". "match_by_threshold" column show whether the SpliceAI predicted splicing consequences of a variant match a unannotated junction; "pred_cryptic_exon" column show whether SpliceAI predict a cryptic exon (gain of an unannotated donor site and an unannotated acceptor site at the same time). "csQTL_candidate" column suggest whether the variant and unannotated junction is a sQTL candidate, based on the provided thresholds.


#### Unpaired WGS and RNAseq analyses

For unpaired WGS and RNAseq analysis, the argument ```junc_anno_file``` can take the output of both ```annotateLeafCutterJunc``` and ```annotateQueryJunc``` function.
```{r}
collapsedsQTL(junc_anno_file = "example_chr1_leafcutter_outliers.txt.gz",
              paired_data = FALSE,
              gdb_path = "example_rvatData.gdb",
              cohort_name = "pheno",
              max_allele_freq = 0.05,
              reference = "Default",
              spliceai_prediction=list(SpliceAI = c("spliceaiDS_AG", "spliceaiDS_AL", "spliceaiDS_DG", "spliceaiDS_DL")),
              min_spliceai_score = 0.2,
              spliceai_default_reference = T,
              output_prefix = "example_chr1")
```
The nominated csQTLs are written in example_chr1_csQTL_unpaired_csQTL.txt.gz file.

```{r}
collapsedsQTL(paired_data = FALSE,
              junc_anno_file = "example_chr1_query_junc_anno.txt.gz",
              gdb_path = "example_rvatData.gdb",
              cohort_name = "pheno",
              max_allele_freq = 0.05,
              reference = "Default",
              spliceai_prediction=list(SpliceAI = c("spliceaiDS_AG", "spliceaiDS_AL", "spliceaiDS_DG", "spliceaiDS_DL")),
              min_spliceai_score = 0.2,
              spliceai_default_reference = T,
              output_prefix = "example_chr1_query_junc")
```
The output is written in example_chr1_query_junc_csQTL_unpaired_csQTL.txt.gz file.

### Step 4. Run SpliPath data browser

#### Step 4.1 Prepare data for data browser

To prepare data for visualization, first, we separate the junction files in per gene and put them into a directory (e.g. browser_junction_paired) by:

#### Paired WGS and RNAseq analyses
For paired data analyses, the argument ```junc_count_files``` and ```junc_anno_files``` take the output files of function ```annotateLeafCutterJunc```; the argument ```cs_qtl_files``` take the output of ```collapsedsQTL``` function:
```{r}
browserData(paired_data = TRUE,
            junc_count_files = c("example_chr1.txt.gz"), 
            junc_anno_files = c("example_chr1_anno.txt.gz"), 
            cs_qtl_files = c("example_chr1_paired_csQTL.txt.gz"), 
            tissues_leafcutter_pvals_file = tissues_leafcutter_pvals_file,
            output_dir = "browser_junction_paired")
```

#### Unpaired WGS and RNAseq analyses
For unpaired data analyses, the argument ```junc_anno_file``` can take the output of both ```annotateLeafCutterJunc``` and ```annotateQueryJunc``` function; the argument ```cs_qtl_files``` take the output of ```collapsedsQTL``` function:
```{r}
browserData(paired_data = FALSE,
            junc_anno_files = c("example_chr1_leafcutter_anno.txt.gz"), 
            cs_qtl_files = c("example_chr1_unpaired_csQTL.txt.gz"), 
            output_dir = "browser_junction_unpaired")
```

```{r}
browserData(paired_data = FALSE,
            junc_anno_files = c("example_chr1_query_junc_anno.txt.gz"), 
            cs_qtl_files = c("example_chr1_query_junc_unpaired_csQTL.txt.gz"), 
            output_dir = "browser_query_junction_unpaired")
```

* If Step 2 and 3 were performed per chromosome, arguments ```junc_count_files```, ```junc_anno_files```, and ```cs_qtl_files``` should be vectors of matched file names for all the chromosomes. 

In the output directory, there will be four types of files:
1) "geneID_geneName_intron_count.txt.gz" are junction-sample read count table per gene,
2) "geneID_geneName_leafcutter_pval.txt.gz" are junction-sample LeafCutter P values table per gene,
3) "geneID_geneName_csQTL.txt.gz" are the mapped ur-sQTL candidates per gene.  

The other files necessary for SpliPath browser are (the input or output files in Step 1-3):
1) GDB file (e.g. example_rvatData.gdb)
Two sample mata data only necessary for visualizing paired data:
2) DNA sequencing sample meta data file (e.g. example_DNA_meta.txt)
3) RNAseq sample meta data file (e.g. example_RNAseq_meta.txt)

#### Step 4.2 Launch data browser

If user preformed the data analysis in a remote computer cluster and prefer to launch data browser on a local computer, user can download the files mentioned in Step 4.1 to local and then run ```runDataBroswer``` function:

#### Paired WGS and RNAseq analyses
```{r}
# 'data_dir' should be the same with 'output_dir' argument in browserData function.

library(SpliPath)
runDataBroswer(paired_data = TRUE,
               rna_meta = "example_RNAseq_meta.txt", 
               dna_meta = "example_DNA_meta.txt", 
               gdb_path = "example_rvatData.gdb",
               data_dir = "browser_junction_paired",
               reference = "Default")
```

#### Unpaired WGS and RNAseq analyses
```{r}
runDataBroswer(paired_data = FALSE,
               gdb_path = "example_rvatData.gdb",
               data_dir = "browser_query_junction_unpaired",
               reference = "Default")

```

The following graphic tutorial show how to use SpliPath to anwser splicing pathology questions:



