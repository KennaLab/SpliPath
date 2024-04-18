#' Map potential splice-altering DNA variants to novel junctions in one tissue 
#' 
#' For each individual, map the DNA variants that of splice-altering potential to a novel junction in this individual's tissue.  
#' 
#' @returns data.frame One row in the the dataframe is an annotated DNA variant that mapped to a novel junction in a subject's tissues. The junction PSI and read count in the are in columns "PSI_tissue" and "Read_tissue".
#' @export 
tissueVar2Junc <-
function(rna_meta, tissue, var2junc, psi, read_count, leafcutter_pvals, geno, var_pred){
  rna_meta = rna_meta[rna_meta$Tissue == tissue, ]
  rna_meta = rna_meta[!duplicated(rna_meta$SubjectID), ]
  rna_meta = rna_meta[rna_meta$SubjectID %in% colnames(geno), ]
  
  var2junc$junc_coord = do.call(rbind, strsplit(var2junc$junc, split=":EN"))[,1]
  # Remove junctions that are not in the leafcutter test results
  var2junc = var2junc[var2junc$junc_coord %in% rownames(leafcutter_pvals),  ]  
  
  tissue_read_count = read_count[var2junc$junc_coord, rna_meta$SampleID]
  colnames(tissue_read_count) = rna_meta$SubjectID
  
  tissue_leafcutter_pvals = leafcutter_pvals[var2junc$junc_coord, rna_meta$SampleID]
  colnames(tissue_leafcutter_pvals) = rna_meta$SubjectID
  
  # tissue_psi = psi[var2junc$junc, rna_meta$SampleID]
  # colnames(tissue_psi) = rna_meta$SubjectID
  
  tissue_var = geno[as.character(var2junc$VAR_id), colnames(tissue_read_count)]
  tissue_var[tissue_var > 0] = 1
  tissue_var2rc = tissue_var * tissue_read_count
  tissue_map_idx = which(tissue_var2rc > 0, arr.ind = T)
  rm(tissue_var2rc)
  
  var_ls = do.call(rbind, strsplit(row.names(tissue_var)[unname(tissue_map_idx[,"row"])], split = "[.]"))[,1]
  junc_ls = do.call(rbind, strsplit(row.names(tissue_read_count)[unname(tissue_map_idx[,"row"])], split = "[.]"))[,1]
  tissue_map = data.frame(VAR_id=var_ls, 
                          DNA_variant = var2junc[match(var_ls, var2junc$VAR_id), "DNA_variant"], 
                          SubjectID=colnames(tissue_read_count)[unname(tissue_map_idx[,"col"])], 
                          Coordinates_of_novel_junc=junc_ls, 
                          Event = var2junc[as.numeric(unname(tissue_map_idx[,"row"])), "event"],
                          Gene=var2junc[as.numeric(unname(tissue_map_idx[,"row"])), "gene.name"],
                          Gene_id=var2junc[as.numeric(unname(tissue_map_idx[,"row"])), "gene.id"],
                          Var_region = var2junc[as.numeric(unname(tissue_map_idx[,"row"])), "var.region"],
                          Strand = var2junc[as.numeric(unname(tissue_map_idx[,"row"])), "region.strand"]
                          
                          # Event = var2junc[match(junc_ls, var2junc$junc), "event"],
                          # Gene=var2junc[match(junc_ls, var2junc$junc), "gene.name"], 
                          # Gene_id=var2junc[match(junc_ls, var2junc$junc), "gene.id"], 
                          # Var_region = var2junc[match(junc_ls, var2junc$junc), "var.region"], 
                          # Strand = var2junc[match(junc_ls, var2junc$junc), "region.strand"]
                          )
  for (pred in var_pred){
    tissue_map[[pred]] = var2junc[match(var_ls, var2junc$VAR_id), pred]
  }
  #tissue_map[[paste0("PSI_", tissue)]]=tissue_psi[tissue_map_idx]
  tissue_map[[paste0("Read_", tissue)]]=tissue_read_count[tissue_map_idx]
  tissue_map[[paste0("LeafCutter_pVal_", tissue)]]=tissue_leafcutter_pvals[tissue_map_idx]
  tissue_map[[paste0("PSI_", tissue)]] = 0
  
  tissue_map
}
