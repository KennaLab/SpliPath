#  Splicing Pathology (SpliPath)

## Introduction

SpliPath is developed to nominate ultra-rare splicing quantitative trait loci (ur-sQTL) in genomewide analyses of paried DNAseq and RNA splicing datasets. It serves to mitigate key challenges of distinguishing rare bona fide aberrant splicing events from abundant artefects, and to identify causal genetic factors responsible for these events. SpliPath provides genomewide splice junction annotation, ur-sQTL nomination and an interactive Shiny-based data visualization and analysis tool. It enables users to look beyond consensus splice sites for causal genetic variants in deeper intronic regions, and can be applied to analyse any dataset that includes matched DNA and RNA sequencing data for donors believed to harbour rare pathogenic splice-altering variants. 

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

This tutorial shows how to nominate ur-sQTL by integrating outlier splicing analysis (LeafCutterMD) with splice-altering variants predictions (SpliceAI and dbscSNV) of paired DNA-RNAseq samples. First, we annotate novel splice junctions. Then, we putatively link them with predicted splice-altering variants and nominate ultra rare ur-sQTL by finding splice-altering variants that were linked to novel splicing junctions with moderate or strong evidence. Finally, we show how to use SpliPath data browser to view splicing evidence and facilitate linking rare splice-altering variants with consequent rare aberrant splicing events. Considering the size of human genome and transcriptome, we recommend to perform the data analyses (Step 1-4.1) using a high performance computing cluster, while data visualization (Step 4.2) can be performed on any standard laptop.

* The example input and output files can be found in ```SpliPath/Example_data/``` (download by ```git clone https://github.com/KennaLab/SpliPath.git``` if the package was installed by ```remotes::install_github``` ).

### Step 1. Data preperation

SpliPath performs paired DNA-RNAseq analysis using:
1 A tab-delimited metadata table of the RNA sequencing samples. It should at least contain 'SubjectID', 'SampleID', 'Group', 'Tissue' columns and a 'Path' column which specify the paths to the splicing junction BED files of each sample (e.g. ```SpliPath/Example_data/example_RNAseq_meta.txt``` )
2. Splicing junction data in BED format (e.g. ```SpliPath/Example_data/regtools_junction/*``` )
3. LeafCutter splicing outlier analysis output (e.g. ```SpliPath/Example_data/example_LeafCutter_outlier_pVals.txt``` )

4. A tab-delimited metadata table of the DNA sequencing samples. It should at least contain 'SubjectID', 'SampleID', and 'Group' (e.g. ```SpliPath/Example_data/example_DNA_meta.txt``` )
5. Genomic variants in VCF format (e.g.  ```SpliPath/Example_data/example_rvatData.vcf``` ). It should includes all samples in the metadata table.
6. SpliceAI and dbscSNV prediction (e.g. ```SpliPath/Example_data/example_chr1_varAnno_SpliceAI.txt``` & ```SpliPath/Example_data/example_chr1_varAnno_dbscSNV.txt```. There are pre-computed SpliceAI scores (https://github.com/Illumina/SpliceAI), and pre-computed dbscSNV scores (http://www.liulab.science/dbscsnv.html) in VCF format.)

See [Appendix](https://github.com/KennaLab/SpliPath/blob/main/SpliPath_Tutorial_Appendix.md) for the commands of running RegTools, LeafCutter and SpliceAI for generating the input files above.


##### Convert genomic variants VCF to GDB
We use [RVAT](https://github.com/kkenna/rvat) to build GDB from VCF file.
```{r}
library(SpliPath)
library(rvat)

vcfpath <- "example_rvatData.vcf"
vcfpheno <- "example_DNA_meta.txt"
gdbpath <- "example_rvatData.gdb"
buildGdb(vcf = vcfpath,  output = gdbpath, overWrite = TRUE)

# Upload sample information
mygdb = gdb(gdbpath)
cohort_pheno <- read.table(vcfpheno, header = TRUE, sep = "\t")
cohort_pheno$IID = cohort_pheno$SampleID
uploadCohort(object = mygdb, name = "pheno",value = cohort_pheno) # The 'name' parameter is the name of the cohort table in the GDB
```

The variants annotation tables (SpliceAI & dbscSNV predictions) can be uploaded into GDB by 
```{r}
spliceai_pred_path <- "example_chr1_varAnno_SpliceAI.txt"
spliceai_pred <- read.table(spliceai_pred_path, header = TRUE, sep = "\t")
uploadAnno(object = mygdb, name = "SpliceAI",value = spliceai_pred) # The 'name' parameter is the name of the variants annotation table in the GDB

dbscsnv_pred_path <- "example_chr1_varAnno_dbscSNV.txt"
dbscsnv_pred <- read.table(dbscsnv_pred_path, header = TRUE, sep = "\t")
uploadAnno(object = mygdb, name = "dbscSNV", value = dbscsnv_pred) 
```
Users can also upload other variants annotation to GDB and include them into analysis in Step 3.

##### (Optional) Customize reference files
The default reference genomic annotations are from Ensembl GRCh38.98 GTF file and Snaptron annotation. Users can create customized reference file using ```prepareGenoRef``` function.

The argument ```intron_bed``` of ```prepareGenoRef``` function is a BED file contain coordinates (chromosome, start position, end position, strand), transcript id, and gene id information of introns. This file can be generated by [GTFtools] (https://www.genemine.org/gtftools.php) from GTF file, see [Appendix 4](https://github.com/KennaLab/SpliPath/blob/main/SpliPath_Tutorial_Appendix.md)

```{r}
prepareGenoRef(gtf_file = "Homo_sapiens.GRCh38.98.gtf.gz", intron_bed = "Intron_GRCh38.98.sort.bed.gz", output_dir = "new_reference")
```
The output reference files will be written into the provided ```output_dir``` directory. The four output files are 
1) "Gene.bed": coordinates of genes; 
2) "Exon_proteincoding.bed": exons in protein coding transcripts; 
3) "Intron_proteincoding.bed": introns in protein coding transcripts; 
4) "Intron_noncoding.bed" introns in noncoding transcripts.

Or, users can generate reference files following the file format of the default reference files in ```SpliPath/inst/extdata/Reference``` directory. And name them with "Gene.bed", "Exon_proteincoding.bed", "Intron_proteincoding.bed" and "Intron_noncoding.bed".

After generating reference files, specify the directory to the generated reference files in ```reference``` argument in functions of Step 2-4. 


### Step 2. Annotate splicing junctions

Considering computational capacity, we highly recommend performing Step 2-4 per chromosome.

After data preparation, SpliPath starts with annotating the splicing junctions. The annotation include mapping the junctions to genes, finding novel splice sites, finding novel junctions and annotate splice events. Excluding possible gene fusion events and artefacts where the two splice sites of a junction are mapped two genes or neither of them are annotated. 

```{r}
annotateJunc(rna_meta = "example_RNAseq_meta.txt", chrom = 1, output_prefix = "example_chr1", reference = "Default")
```

The function outputs two files: 
1. example_chr1.txt.gz, which is junction read count table merged from all the samples. The unique junctions are in rows and samples are in columns; 
2. example_chr1_anno.txt.gz, which is the annotation of the unique junctions. 

### Step 3. Map potential splice-altering variants to novel junctions and nominate ur-sQTL

The ur-sQTL candidates is nominated by function ```mapVariantJunc```.

The arguments are 
1. Junction and variants inputs:
* ```junc_anno_file``` and ```read_count_file``` are the files generated in Step 2.
* ```tissues_leafcutter_pvals_file``` lists paths to the LeafCutter analysis P values of each tissue type. Each element is a LeafCutter file path named under a tissue. The tissues should be the same with those in the RNAseq sample metadata file.
* ```gdb_path``` is the GDB file of DNA variants generated in Step 1. The table of phenotype data should be specified in ```cohort_name```.
* ```rna_meta``` and ```wgs_meta``` is the file of phenotype / metadata of the RNA and DNA sequencing samples. Both of them should have 'SubjectID', 'SampleID', and 'Group' to specify which individual the sample is collected and which phenotype group the individuals belong to. In addition, 'rna_meta' should have 'Tissue' column, which specify the tissue source of the RNAseq sample.
* ```as_annotated``` defines annotation sources (should be in colnames in junc_anno_file) that were combined to determine whether a junction is novel. 

2. Thresholds to nominate ur-sQTL
* ```max_LeafCutter_Pval``` is the maximum LeafCutter P value of a novel junction to nominate ultra rare sQTL candidate
* ```min_nr_tissue``` requires LeafCutter P value of a junction should be under the threshold in at least this number of tissues
* ```source_allele_freq``` is the source of the variants allele frequency. If set to "cohort", the allele frequency is calcuated from the cohort provided by the argument "cohort_name"
* ```max_allele_freq``` is the aximum allele frequency to nominate ultra rare sQTL candidate 
* ```splice_prediction``` is the splice-altering prediction to annotate the DNA variants. A splice prediction tool may provide more then one score. The   
* ```min_splice_score``` is the min splice score of each prediction tools to nominate ur-sQTL candidate.
* ```SpliceAI_default_reference``` suggests whether SpliceAI predictions are pre-computed or SpliceAI default annotation were used when running SpliceAI

```{r}
var2junc = mapVariantJunc( junc_anno_file = "example_chr1_anno.txt.gz", 
                           read_count_file = "example_chr1.txt.gz", 
                           tissues_leafcutter_pvals_file = list(motor_cortex = "leafcutter/example_LeafCutter_motor_cortex_outlier_pVals.txt", 
                                                                 frontal_cortex = "leafcutter/example_LeafCutter_frontal_cortex_outlier_pVals.txt", 
                                                                 cervical = "leafcutter/example_LeafCutter_cervical_outlier_pVals.txt", 
                                                                 lumbar = "leafcutter/example_LeafCutter_lumbar_outlier_pVals.txt"), 
                           gdb_path = "example_rvatData.gdb", 
                           cohort_name = "pheno", 
                           rna_meta = "example_RNAseq_meta.txt", 
                           wgs_meta = "example_DNA_meta.txt",
                           reference = "Default",
                           as_annotated = c("in.ensembl", "in.snaptron"),
                           max_LeafCutter_Pval = 0.05,
                           min_nr_tissue = 1,
                           source_allele_freq = "cohort", 
                           max_allele_freq = 0.05, 
                           splice_prediction=list(SpliceAI = c("spliceaiDS_AG", "spliceaiDS_AL", "spliceaiDS_DG", "spliceaiDS_DL"), 
                                                  dbscSNV = c("rf_score", "ada_score")), 
                           min_splice_score = list(SpliceAI = 0.2, dbscSNV = 0.7),
                           SpliceAI_default_reference = T,
                           output_prefix = "example_chr1")
                                  
```
The output data.frame is also written in example_chr1_crossref.txt.gz. One row in the the data.frame is an annotated DNA variant that mapped to a novel junction in a subject's tissues. The LeafCutter P values and read count in tissues are in columns "LeafCutter_pVal_tissue" and "Read_tissue". "SpliceAI_pred_match_junction" column show whether the SpliceAI predicted splicing consequences of a variant match a novel junction; "SpliceAI_pred_cryptic_exon" column show whether SpliceAI predict a cryptic exon. "sQTL_candidate" column suggest whether the variant and novel junction is a sQTL candidate, based on the provided thresholds.

### Step 4. Run SpliPath data browser

#### Step 4.1 Input data of data browser

To prepare data for visualization, first, separate the junction files in per gene and put them into a directory (e.g. browser_junction) by:
```{r}
browserData(junc_count = c("example_chr1.txt.gz"), 
            junc_anno = c("example_chr1_anno.txt.gz"), 
            var2junc = c("example_chr1_crossref.txt.gz"), 
            tissues_leafcutter_pvals_file = list(motor_cortex = "leafcutter/example_LeafCutter_motor_cortex_outlier_pVals.txt", 
                                                 frontal_cortex = "leafcutter/example_LeafCutter_frontal_cortex_outlier_pVals.txt", 
                                                 cervical = "leafcutter/example_LeafCutter_cervical_outlier_pVals.txt", 
                                                 lumbar = "leafcutter/example_LeafCutter_lumbar_outlier_pVals.txt"),
            output_dir = "browser_junction")
```
* If Step 2 and 3 were performed per chromosome, input files of ```junc_count```, ```junc_anno```, ```var2junc```, should be vectors of matched file names for all the chromosomes. 

In the output directory, there will be four types of files:
1) "number_of_novel_junction_per_subject.txt" records the total number of novel junctions per RNAseq sample,
2) "geneID_geneName_intron_count.txt.gz" are junction-sample read count table per gene,
3) "geneID_geneName_leafcutter_pval.txt.gz" are junction-sample LeafCutter P values table per gene,
4) "geneID_geneName_crossref.txt.gz" are the mapped ur-sQTL candidates per gene.  

The other files necessary for SpliPath browser are (the input or output files in Step 1-3):
1) RNAseq sample meta data file (e.g. example_RNAseq_meta.txt)
2) DNA sequencing sample meta data file (e.g. example_DNA_meta.txt)
3) GDB file (e.g. example_rvatData.gdb)

#### Step 4.2 Launch data browser

If user preformed the data analysis in a remote computer cluster and prefer to launch data browser on a local computer, user can download the files mentioned in Step 4.1 to local and then run ```runDataBroswer``` function:
```{r}
# 'data_dir' should be the same with 'output_dir' argument in browserData function.

library(SpliPath)
library(ggVennDiagram)
runDataBroswer(rna_meta = "example_RNAseq_meta.txt", 
               dna_meta = "example_DNA_meta.txt", 
               gdb_path = "example_rvatData.gdb",
               data_dir = "browser_junction",
               reference = "Default")
```

The following graphic tutorial show how to use SpliPath to anwser splicing pathology questions:
![](https://github.com/KennaLab/SpliPath/blob/main/inst/SpliPath_Browser/www/Guide.png)


