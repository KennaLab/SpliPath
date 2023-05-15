#  Splicing Pathology (SpliPath)

## Introduction

SpliPath is developed to nominate novel pathogenic mutation hotspots in genomewide analyses of paried DNAseq and RNA splicing datasets. It serves to mitigate key challenges of distinguishing rare bona fide aberrant splicing events from abundant technical artefects, and to identify causal genetic factors responsible for these events. SpliPath provides genomewide splice junction annotation, mutation hotspots nomination and an interactive Shiny-based data visualization and analysis tool. It enables users to look beyond consensus splice sites for causal genetic variants in deeper intronic regions, and can be applied to analyse any dataset that includes matched DNA and RNA sequencing data for donors believed to harbour rare pathogenic splice-altering variants (sRV). 

## Installation

### Install from github

```{sh}
remotes::install_github("KennaLab/SpliPath")
```

### Clone repo

Clone repo

```{sh}
git clone https://github.com/KennaLab/SpliPath.git
```
Open an R session and install using devtools 

```{r}
devtools::install("SpliPath") # if the SpliPath directory is not in your current directory,  provide the full path to the cloned repo
```

## Tutorial

This tutorial shows how to nominate splice variants hotspots by paired DNA-RNAseq analyses. We start from preparing splicing junction data from aligned RNAseq BAM file and DNA variants data from VCF file. Next, we annotate novel splice junctions and putatively link them with predicted splice-altering variants. Then, we nominate splice variants hotspots by finding splice-altering variants that were linked to novel splicing junctions with moderate or strong evidence. Finally, we show how to use SpliPath data browser to view splicing evidence and facilitate linking rare splice-altering variants with consequent rare aberrant splicing events. Considering the size of human genome and transcriptome, we recommend to perform the data analyses (Step 1-7.1) using a high performance computing cluster, while data visualization (Step 7.2) can be performed on any standard laptop.

* The example input and output files can be found in ```SpliPath/Example_data/``` .

### Step 1. Data preperation

SpliPath performs paired DNA-RNAseq analysis using splicing junction data in BED format and genomic variants in GDB format. 

#### Step 1.1 Splicing junction data from RNAseq alignment

We recommend [RegTools](https://regtools.readthedocs.io/en/latest/) to extract junction data from split-reads with minimum 6bp anchors in aligned RNAseq BAM files, using the command:
```{sh}
regtools junctions extract -a 6 -s 1 -o ${bamfile}.junc ${bamfile} 
awk -F'[\t|,]' '$6 != "?" {print $1,$2+$13,$3-$14,$4,$5,$6}' OFS='\t' ${bamfile}.junc > ${bamfile}.bed
```
Option ```-s``` in ```regtools junctions extract``` means strand specificity of RNA library preparation, please change it according to the input RNAseq samples.

A tab-delimited metadata table of the RNAseq samples should be prepared for further use. It should have at least "SubjectID", "SampleID", "Group", and "Tissue" columns which specify the individuals the samples were collected, the phenotype group the individuals belong to, and the tissue source the sample were collected. 

#### Step 1.2 Genomic variants in GDB format

We use [RVAT](https://github.com/kkenna/rvat) to create GDB from VCF file. Please follow the [Setting up and populating a gdb](https://kkenna.github.io/rvat/articles/basics.html) steps in RVAT tutorial. 

Here, we use the GDB in rvatData pacakge as an example. 
```{r}
library(rvat)
library(rvatData)
library(SummarizedExperiment)

# Load the example GDB object.
gdb_path = rvat_example("rvatData.gdb")
gdb = gdb(gdb_path)
```

A tab-delimited metadata table of the DNA sequencing samples should be prepared for further use. It should have at least "SubjectID" and "SampleID" columns which specify the individuals the samples were collected. The SampleIDs in the metadata table should be in the "IID" field in the GDB cohort table. 

#### Step 1.3 SplicAI and dbscSNV for splice-altering prediction

There are pre-computed SpliceAI scores (https://github.com/Illumina/SpliceAI), and pre-computed dbscSNV scores (http://www.liulab.science/dbscsnv.html) in VCF format. Or, users can run SpliceAI with customized parameters.
The output VCF files contain splice-altering predictions or any other annotations should be first transformed to a table, where the rows are the variants and columns are CHROM, POS, ID, REF, ALT, and any other annotations (For SpliceAI, the columns are "spliceaiDS_AG", "spliceaiDS_AL", "spliceaiDS_DG", "spliceaiDS_DL", "spliceaiDP_AG", "spliceaiDP_AL", "spliceaiDP_DG", "spliceaiDP_DL", and etc. For dbscSNV, the columns are "ada_score", "rf_score", and etc.). 

The table can be uploaded into GDB by 
```{r}
spliceai_pred_path <- "example_chr1_varAnno_SpliceAI.txt" # example variant annotation file in ```SpliPath/Example_data/```
spliceai_pred <- read.table(spliceai_pred_path, header = TRUE, sep = "\t")
uploadAnno(object = gdb, name = "SpliceAI",value = spliceai_pred) # The 'name' parameter is the name of the variants annotation table in the GDB

dbscsnv_pred_path <- "example_chr1_varAnno_dbscSNV.txt"
dbscsnv_pred <- read.table(dbscsnv_pred_path, header = TRUE, sep = "\t")
uploadAnno(object = gdb, name = "dbscSNV", value = dbscsnv_pred) 
```
Users can also upload other variants annotation to GDB and include them into analysis in Step 5.

### Step 2. Merge splicing junction data to junction-sample data table

After data preparation, start SpliPath analysis by merging splicing junctions per sample to a data frame. 
Considering computational capacity, we highly recommend to perform Step 2-7 per chromosome. E.g. chromosome 1
```{r}
library(SpliPath)
# 'sample_path_file' is a tab-delimited file which contain two columns: 1, unique RNAseq sample ids; 2, paths to the junction BED files (Step 1.1) of the corresponding sample.

sample_path_file = "sample_path.txt"
mergeJuncData(sample_path_file, chrom=1, output_prefix = "example")
```
The output junction-sample data table will be written in example_chr1.txt.gz file. The unique junctions are in rows and samples are in columns. The first four columns are the coordinates of junctions (chr, start, end, and strand). 

### Step 3. Annotate splicing junctions

The annotation include mapping the junctions to genes, finding novel splice sites, finding novel junctions and annotate splice events. Excluding possible gene fusion events and artefacts where the two splice sites of a junction are mapped two genes or neither of them are annotated. 
```{r}
annotateJunc(merged_junc = "example_chr1.txt.gz", output_prefix = "example_chr1", reference = "Default")
```
The output junction annotation will be written in example_chr1_anno.txt.gz file. 

* The default reference genomic annotations are from Ensembl GRCh38.98 GTF file and Snaptron annotation. Users can create customized reference file using ```prepareGenoRef``` function.
* The argument ```intron_bed``` of ```prepareGenoRef``` function is a BED file contain coordinates (chromosome, start position, end position, strand), transcript id, and gene id information of introns. This file can be generated by [GTFtools] (https://www.genemine.org/gtftools.php) from GTF file: python gtftools.py -i introns.bed ensembl.gtf  
* After generating reference files, specify the directory to the generated reference files in ```reference``` argument in ```annotateJunc``` function and other functions in Step 4-7. 

```{r}
# gtf_file is an Ensembl genomic annotation GTF file

prepareGenoRef(gtf_file = "ensembl.gtf", intron_bed = "introns.bed", output_dir = "new_reference")
```
The output reference files will be written into the provided ```output_dir``` directory. The four output files are 
1) "Gene.bed": coordinates of genes; 
2) "Exon_proteincoding.bed": exons in protein coding transcripts; 
3) "Intron_proteincoding.bed": introns in protein coding transcripts; 
4) "Intron_noncoding.bed" introns in noncoding transcripts.

Or, users can generate reference files following the file format of the default reference files in ```SpliPath/inst/extdata/Reference``` directory. And name them with "Gene.bed", "Exon_proteincoding.bed", "Intron_proteincoding.bed" and "Intron_noncoding.bed".

### Step 4. Calculate PSI

Next, the percentage spliced in (PSI) of the novel junctions are calculated.
```{r}
# 'junc_file' is the name of junction-sample file, the output of Step 2.
# 'annotation' is the annotation source that included in defining whether the junction is "novel" or "annotated". Should be in the column names in "example_chr1_anno.txt.gz" file.

novel_junc <- calJuncPSI(junc_file = "example_chr1.txt.gz", 
                         junc_anno="example_chr1_anno.txt.gz", 
                         output_prefix = "example_chr1", 
                         annotation = c("in.ensembl", "in.snaptron"))
```
The PSI values of the novel junctions in each sample will be written in example_chr1_norm.txt.gz. The retured ```novel_junc``` object is a 7-column data.frame containing the information of the novel junctions: chr, start, end, strand, gene.id, gene.name, and event.

### Step 5. Map potential splice-altering variants to novel junctions

First, find the annotated junctions and flanking exons that were disrupted by the novel junctions, where the cis-acting splice-altering variants most likely to lie in. 
```{r}
var_region <- juncVariantRegion(novel_junc, output_prefix = "example_chr1", reference = "Default")
```
The output will be written in ./example_chr1_variant_region.txt.gz. This file is the input of the following mapping potential splice-altering variants to novel junctions step:
```{r}
# 'count_file' is the file of junction-sample read count generated by Step 2.
# 'psi_file' is the file of novel junction PSIs generated by Step 3.
# 'gdb_path' is the GDB file of DNA variants generated by Step 1. The table of phenotype data should be specified in 'chort_name'.
# 'rna_meta' and 'wgs_meta' is the file of phenotype / metadata of the RNA and DNA sequencing samples. Both of them should have 'SubjectID', 'SampleID', and 'Group' to specify which individual the sample is collected and which phenotype group the individuals belong to. In addition, 'rna_meta' should have 'Tissue' column, which specify the tissue source of the RNAseq sample.
# 'tissue' specifies RNAseq of which tissue source will be included into analysis.
# 'splice_prediction' is the splice-altering prediction to annotate the DNA variants. Default: list(SpliceAI = c("spliceaiDS_AG", "spliceaiDS_AL", "spliceaiDS_DG", "spliceaiDS_DL"), dbscSNV = c("rf_score", "ada_score")). The names of elements in the list are prediction tools, they should be in GDB database tables. The items are the scoring fields in the GDB database prediction tables. The maximum of the given scores will be used as the prediction of the corresponding tool for each variants. E.g. The SpliceAI score of a variant will be the maximum of spliceaiDS_AG, "pliceaiDS_AL, spliceaiDS_DG, and spliceaiDS_DL. 

tissue_var2junc <- mapVariantJunc(junc_var_region = var_region, 
                                  psi_file = "example_chr1_norm.txt.gz", 
                                  count_file = "example_chr1.txt.gz", 
                                  gdb_path = rvat_example("rvatData.gdb"), 
                                  cohort_name = "pheno", 
                                  rna_meta = "example_RNAseq_meta.txt", 
                                  wgs_meta = "example_DNA_meta.txt", 
                                  tissues = c("tissue1", "tissue2", "tissue3", "tissue4"), 
                                  allele_freq = "gdb",
                                  splice_prediction=list(SpliceAI = c("spliceaiDS_AG", "spliceaiDS_AL", "spliceaiDS_DG", "spliceaiDS_DL"),
                                                         dbscSNV = c("rf_score", "ada_score")), 
                                  output_prefix = "example_chr1")
```
The returned ```tissue_var2junc``` is a data.frame. One row in the the data.frame is an annotated DNA variant that mapped to a novel junction in a subject's tissues. The junction PSI and read count in tissues are in columns "PSI_tissue" and "Read_tissue". The data.frame is also written into example_chr1_crossref.txt.gz

* (Optional) The following function compares whether the SpliceAI predicted splice site gain or loss match the observed novel junctions.
```{r}
spliceai_match <- matchPredNObs(gdb_path = rvat_example("rvatData.gdb"), 
                                tissue_var2junc, 
                                min_spliceai_score = 0.2, 
                                reference = "Default")
```
In the returned ```spliceai_match```, one row in the the data.frame is an annotated DNA variant that mapped to a novel junction in a subject's tissues. "match_by_threshold" column show whether the SpliceAI predicted splicing consequences of a variant match a novel junction; "pred_cryptic_exon" column show whether SpliceAI predict a cryptic exon.  

### Step 6. Nominate mutation hotspots

Users can nominate mutation hotspots by setting thresholds for junction splicing signal, variants splice-altering scores, allele frequency.
```{r}
# 'min_psi' is the minimum PSI required for observed novel junction.
# 'min_read' is the minimum number of reads required to support the observed novel junction. 
# 'min_read_include' means also include junctions supported by more than min_read_include reads even of their PSI are under min_psi.
# 'tissues' specifies the tissues included in the analysis.
# 'min_nr_tissue' means user require novel junction observation concordant across min_nr_tissue tissues.
# 'min_splice_pred' is a list of the minimum prediction score of prediciton tools. Default: list(SpliceAI = 0.2, dbscSNV = 0.7).
# 'splicai_pred_match' Optional. If the matching of SpliceAI prediction and observed novel junctions are provided (see Step 5), the "match_by_threshold" and "pred_cryptic_exon" column will be add into the output. 
# 'blacklist' Optional. A genomic region blacklist. If provided, the regions will be removed.

nominateHotspot(tissue_var2junc, 
                min_psi = 0.2, 
                min_read = 5, 
                min_read_include = 20, 
                tissues = c("tissue1", "tissue2", "tissue3", "tissue4"), 
                min_nr_tissue = 2, 
                min_splice_pred = list(SpliceAI = 0.2, dbscSNV = 0.7), 
                max_allele_freq = 0.05,
                splicai_pred_match = spliceai_match, 
                output_prefix = "example_chr1")
```
The paired DNA variants and novel junctions will be written in example_chr1_hotspot.txt.gz. In this file, column "Hotspot" shows whether they are nominated as mutation hotspots based on provided thresholds.

### Step 7. Run SpliPath data browser

#### Step 7.1 Input data of data browser

To prepare data for visualization, first, separate the junction files in per gene and put them into a directory (e.g. browser_junction) by:
```{r}
browserData(junc_count = c("example_chr1.txt.gz"), 
            junc_psi = c("example_chr1_norm.txt.gz"), 
            junc_anno = c("example_chr1_anno.txt.gz"), 
            hotspot = c("example_chr1_hotspot.txt.gz"), 
            output_dir = "browser_junction")
```
* If Step 2-6 were performed per chromosome, input files of ```junc_count```, ```junc_psi```, ```junc_anno```, ```hotspot```, should be vectors of matched file names for all the chromosomes. 

In the output directory, there will be four types of files:
1) "number_of_novel_junction_per_subject.txt" records the total number of novel junctions per RNAseq sample,
2) "geneID_geneName_intron_count.txt.gz" are junction-sample read count table per gene,
3) "geneID_geneName_intron_psi.txt.gz" are junction-sample PSI table per gene,
4) "geneID_geneName_intron_sRV_candidate.txt.gz" are the candidate mutation hotspots per gene.  

The other files necessary for SpliPath browser are (the input or output files in Step 2-6):
1) RNAseq sample meta data file (e.g. example_RNAseq_meta.txt)
2) DNA sequencing sample meta data file (e.g. example_DNA_meta.txt)
3) GDB file (e.g. rvatData.gdb)
4) Gene reference file (e.g. Gene_GRCh38.98.bed)
5) Exon reference file (e.g. Exon_GRCh38.98.proteincoding.bed)

#### Step 7.2 Launch data browser

After the files mentioned above are prepared, user can run SpliPath data browser by:
```{r}
# 'data_dir' should be the same with 'output_dir' argument in browserData function.

runDataBroswer(rna_meta = "example_RNAseq_meta.txt", 
              dna_meta = "example_DNA_meta.txt", 
              gdb_path = rvat_example("rvatData.gdb"),
              data_dir = "browser_junction",
              reference = "Default")
```
If user preformed the data analysis in a remote computer cluster and prefer to launch data browser on a local computer, user can download the files mentioned in Step 7.1 to local and then run ```runDataBroswer``` function.

The following graphic tutorial show how to use SpliPath to anwser splicing pathology questions:
![](https://github.com/yanwang271/SpliPath/blob/main/SpliPath_Browser/www/Guide.png)


