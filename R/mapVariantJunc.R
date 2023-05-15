#' Map potential splice-altering DNA variants to novel junctions in tissues
#' 
#' For each individual, map the DNA variants that of splice-altering potential to a novel junction in this individual's tissues.
#'
#' @param junc_var_region string or data.frame The output file or returned data.frame of juncVariantRegion function. The file contains the variants mapping region of each novel junction. 
#' @param psi_file string The output file of calJuncPSI function, which contain junction-sample PSI table.
#' @param count_file string The output file of mergeJuncData funciton, which contain junction-sample read count table.
#' @param gdb_path string The path to GDB database that store the genotype data from whole genome sequencing.
#' @param cohort_name string The table name in GDB database that contain sample metadata
#' @param rna_meta string RNAseq sample metadata file. The file should contain at least "SampleID", "SubjectID", "Tissue" columns.
#' @param wgs_meta string WGS sample metadata file. The file should contain at least "SampleID", "SubjectID" columns.
#' @param tissues vector potential splice-altering DNA variants to novel junctions in given tissues.
#' @param allele_freq list or string Allele frequency of the variants. Default: "gdb", the allele frequency among the whole samples in the provided GDB. Or use other allele frequency annotations of the variants recorded in the GDB tables, e.g. list(gdb_table = "gnomAD", af_field = "AF"). The gdb_table and af_field are the name of variant annotation table and field where the allele frequency is in the GDB.  
#' @param splice_prediction list The splice-altering prediction to annotate the DNA variants. Default: list(SpliceAI = c("spliceaiDS_AG", "spliceaiDS_AL", "spliceaiDS_DG", "spliceaiDS_DL"), dbscSNV = c("rf_score", "ada_score")). 
#' The names of elements in the list are prediction tools, they should be in GDB database tables. The items are the scoring fields in the GDB database prediction tables. The maximum of the given scores will be used as the prediction of the corresponding tool for each variants. E.g. The SpliceAI score of a variant will be the maximum of spliceaiDS_AG, "pliceaiDS_AL, spliceaiDS_DG, and spliceaiDS_DL. 
#' @param output_prefix string The mapping results will be written in output_prefix_map.txt.gz
#'
#' @import rvat
#' @import RSQLite
#' @import SummarizedExperiment
#' @import dplyr
#' @import tidyr
#' @import bedr
#' @returns data.frame One row in the the dataframe is an annotated DNA variant that mapped to a novel junction in a subject's tissues. The junction PSI and read count in tissues are in columns "PSI_tissue" and "Read_tissue".
#' @export 
mapVariantJunc <-
function(junc_var_region, psi_file, count_file, gdb_path, cohort_name, rna_meta, wgs_meta, tissues, 
         allele_freq = "gdb",
         splice_prediction=list(SpliceAI = c("spliceaiDS_AG", "spliceaiDS_AL", "spliceaiDS_DG", "spliceaiDS_DL"), 
                                dbscSNV = c("rf_score", "ada_score")), 
         output_prefix = ""){

  if (!is.data.frame(junc_var_region)){
    junc_var_region = read.table(junc_var_region, header=T, sep="\t", stringsAsFactors = F, check.names = F)
  }
  junc_var_region = junc_var_region[, c("chr", "region.start", "region.end", "strand", "junc", "event", "gene.name", "gene.id")]
  colnames(junc_var_region)[2:3] = c("start", "end")
  junc_var_region[, c("start", "end")] = sapply(junc_var_region[, c("start", "end")], FUN = as.numeric)
  junc_var_region = junc_var_region[order(junc_var_region$chr, junc_var_region$start), ]
  
  # Get variants and annotation from GDB
  gdb = rvat::gdb(gdb_path)
  
  pred_field = paste(paste0(names(splice_prediction), ".*"), collapse = ", ")
  query = sprintf("select VAR_id, var.CHROM, var.POS, var.REF, var.ALT, %s from var ", pred_field)
  for (item in names(splice_prediction)){
    query = paste0(query, sprintf(" left join %s using ('VAR_id')", item))
  }
  query = paste0(query, sprintf(" where var.CHROM in ('%s') ", paste(c(paste0('chr', unique(junc_var_region$chr)), sub("^chr", "", unique(junc_var_region$chr))), collapse = "', '")))
  
  vars = RSQLite::dbGetQuery(gdb, query)
  vars[is.na(vars)] = 0
  for (item in names(splice_prediction)){
    vars[, splice_prediction[[item]]] = sapply(vars[, splice_prediction[[item]]], as.numeric)
    vars[[item]] = apply(vars[, splice_prediction[[item]]], 1, max)
    vars = vars[, !colnames(vars) %in% splice_prediction[[item]]]
  }
  vars = vars[apply(vars[, names(splice_prediction)], 1, max) > 0, ]
  
  if (!grepl(junc_var_region$chr[1], 'chr', fixed = T)){
    vars$CHROM = sub("^chr", "", vars$CHROM)
  }else{
    if (!grepl(vars$CHROM [1], 'chr', fixed = T)){
      vars$CHROM = paste0("chr", vars$CHROM)
    }
  }
  vars$POS = as.numeric(vars$POS)
  vars$start = vars$POS - 1
  vars$DNA_variant = paste0(vars$CHROM, ":g.", vars$POS, vars$REF, ">", vars$ALT)
  vars = vars[, c("CHROM", "start", "POS", names(splice_prediction), "VAR_id", "DNA_variant")]
  colnames(vars) = c("chr", "start", "end", names(splice_prediction), "VAR_id", "DNA_variant")
  vars = vars[order(vars$chr, vars$start), ]
  
  ### Map variant to aberrant splicing region
  var2junc = bedr::bedr(engine = "bedtools", 
                         params = "-d -t all", 
                         input = list(a = vars, b = junc_var_region ), 
                         method = "closest", check.chr = FALSE , check.merge=F, check.valid = F)
  colnames(var2junc) = c("chr", "pos_1", "pos", names(splice_prediction), "VAR_id", "DNA_variant", 
                         "region.chr", "region.start", "region.end", "region.strand", "junc", "event", "gene.name", "gene.id", "dist")
  rm(junc_var_region, vars)
  var2junc = subset(var2junc, dist == 0, select = -pos_1)
  var2junc$var.region = paste(var2junc$region.chr, var2junc$region.start, var2junc$region.end, sep=":")

  rna_meta = read.table(rna_meta, sep="\t", header = T, stringsAsFactors = F)
  rnaid2subject = lapply(split(rna_meta$SubjectID, rna_meta$SampleID), FUN=unique)
  wgs_meta = read.table(wgs_meta, sep="\t", header = T, stringsAsFactors = F)
  wgsid2subject = lapply(split(wgs_meta$SubjectID, wgs_meta$SampleID), FUN=unique)
  subject2wgsid = lapply(split(wgs_meta$SampleID, wgs_meta$SubjectID), FUN=unique)
  
  geno = rvat::getGT(gdb, VAR_id=unique(var2junc$VAR_id), cohort=cohort_name)
  if (allele_freq == "gdb"){
    AF = data.frame(rvat::getAF(geno))
  }else if (class(allele_freq) == "list"){
    AF = rvat::getAnno(gdb, VAR_id = unique(var2junc$VAR_id), table = allele_freq[["gdb_table"]], fields = allele_freq[["af_field"]])
  }
  colnames(AF) = c("AF")
  
  geno = as.data.frame(SummarizedExperiment::assays(geno)$GT, check.names=F)
  geno[is.na(geno)] = 0
  geno = geno[, colnames(geno) %in% wgs_meta$SampleID]
  colnames(geno) = unname(unlist(wgsid2subject[colnames(geno)]))
  rm(wgs_meta, wgsid2subject, subject2wgsid)
  
  psi = read.table(psi_file, header=T, sep="\t", stringsAsFactors = F, check.names = F)
  psi = psi[!duplicated(psi[, c("chr", "start", "end", "strand", "gene.id")]), ]
  rownames(psi) = paste(psi$chr, psi$start, psi$end, psi$strand, psi$gene.id, sep=":")
  
  read_count = read.table(count_file, sep="\t", header=T, stringsAsFactors = F, check.names = F)
  row.names(read_count) = paste(read_count$chr, read_count$start, read_count$end, read_count$strand, sep=":")
  
  tissue_map = data.frame()
  for (tissue in tissues){
    if (nrow(tissue_map) == 0){
      tissue_map = tissueVar2Junc(rna_meta, tissue, var2junc, psi, read_count, geno,  names(splice_prediction))
    }else{
      tissue_map = dplyr::full_join(tissue_map, 
                             tissueVar2Junc(rna_meta, tissue, var2junc, psi, read_count, geno,  names(splice_prediction)))
    }
  }
  tissue_map[is.na(tissue_map)] = 0
  if (nrow(tissue_map) == 0){warning("No variant was mapped to novel junctions")}
  
  # Add AF column
  tissue_map$AF = AF[match(as.character(tissue_map$VAR_id), rownames(AF)), "AF"]
  
  tissue_map = data.frame(lapply(tissue_map, as.character),  stringsAsFactors=F)
  
  output_file = gzfile(sprintf("%s_crossref.txt.gz", output_prefix), "w")
  write.table(tissue_map, output_file, sep="\t", quote=F, row.names = F)
  close(output_file)
  
  tissue_map
}

