#' Map potential splice-altering DNA variants to novel junctions in one tissue 
#' 
#' For each individual, map the DNA variants that of splice-altering potential to a novel junction in this individual's tissue.  
#' 
#' @returns data.frame One row in the the dataframe is an annotated DNA variant that mapped to a novel junction in a subject's tissues. The junction PSI and read count in the are in columns "PSI_tissue" and "Read_tissue".
#' @export 
tissueVar2Junc <-
function(rna_meta, tissue, var2junc, psi, read_count, geno, var_pred){
  rna_meta = rna_meta[rna_meta$Tissue == tissue, ]
  rna_meta = rna_meta[!duplicated(rna_meta$SubjectID), ]
  rna_meta = rna_meta[rna_meta$SubjectID %in% colnames(geno), ]

  raw_junc_coord = do.call(rbind, strsplit(var2junc$junc, split=":EN"))[,1]
  tissue_read_count = read_count[raw_junc_coord, rna_meta$SampleID]
  colnames(tissue_read_count) = rna_meta$SubjectID
  
  tissue_psi = psi[var2junc$junc, rna_meta$SampleID]
  colnames(tissue_psi) = rna_meta$SubjectID
  tissue_var = geno[as.character(var2junc$VAR_id), colnames(tissue_psi)]
  tissue_var[tissue_var > 0] = 1
  tissue_var2psi = tissue_var * tissue_psi 
  tissue_map_idx = which(tissue_var2psi > 0, arr.ind = T)
  rm(tissue_var2psi)
  
  var_ls = do.call(rbind, strsplit(row.names(tissue_var)[unname(tissue_map_idx[,"row"])], split = "[.]"))[,1]
  junc_ls = do.call(rbind, strsplit(row.names(tissue_psi)[unname(tissue_map_idx[,"row"])], split = "[.]"))[,1]
  tissue_map = data.frame(VAR_id=var_ls, 
                          DNA_variant = var2junc[match(var_ls, var2junc$VAR_id), "DNA_variant"], 
                          SubjectID=colnames(tissue_psi)[unname(tissue_map_idx[,"col"])], 
                          Coordinates_of_novel_junc=junc_ls, 
                          Event = var2junc[match(junc_ls, var2junc$junc), "event"],
                          Gene=var2junc[match(junc_ls, var2junc$junc), "gene.name"], 
                          Gene_id=var2junc[match(junc_ls, var2junc$junc), "gene.id"], 
                          Var_region = var2junc[match(junc_ls, var2junc$junc), "var.region"], 
                          Strand = var2junc[match(junc_ls, var2junc$junc), "region.strand"]
                          )
  for (pred in var_pred){
    tissue_map[[pred]] = var2junc[match(var_ls, var2junc$VAR_id), pred]
  }
  tissue_map[[paste0("PSI_", tissue)]]=tissue_psi[tissue_map_idx]
  tissue_map[[paste0("Read_", tissue)]]=tissue_read_count[tissue_map_idx]
  
  tissue_map
}
