#' Map potential splice-altering DNA variants to novel junctions in one tissue 
#' 
#' For each individual, map the DNA variants that of splice-altering potential to a novel junction in this individual's tissue.  
#' 
#' @returns data.frame One row in the the dataframe is an annotated DNA variant that mapped to a novel junction in a subject's tissues. The junction PSI and read count in the are in columns "PSI_tissue" and "Read_tissue".
#' @export 
getVarCarrier <-
  function(rna_meta, wgs_meta, varid_ls, mygdb, cohort_name){
    ### Get geno table
    rna_meta = read.table(rna_meta, sep="\t", header = T, stringsAsFactors = F)
    wgs_meta = read.table(wgs_meta, sep="\t", header = T, stringsAsFactors = F)
    wgsid2subject = lapply(split(wgs_meta$SubjectID, wgs_meta$SampleID), FUN=unique)
    
    geno = rvat::getGT(mygdb, VAR_id=unique(varid_ls), cohort=cohort_name)
    geno = rvat::flipToMinor(geno)
    # if (source_allele_freq == "cohort"){
    #   AF = data.frame(rvat::getAF(geno))
    # }else if (class(source_allele_freq) == "list"){
    #   AF = rvat::getAnno(mygdb, VAR_id = unique(varid_ls), table = source_allele_freq[["gdb_table"]], fields = source_allele_freq[["af_field"]])
    # }
    # colnames(AF) = c("AF")
    
    geno_carrier = rvat::getCarriers(geno)
    geno_carrier = geno_carrier[geno_carrier$IID %in% wgs_meta$SampleID, ]
    geno_carrier$SubjectID = unname(unlist(wgsid2subject[geno_carrier$IID]))
    # geno_carrier$AF = AF[as.character(geno_carrier$VAR_id), "AF"]
    
    geno_carrier[, c("VAR_id", "SubjectID", "genotype")]
  }
