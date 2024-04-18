### Simulate regtools BED from junc.txt

template = read.table("regtools_junc/CGND-HRA-00660.Aligned.sortedByCoord.out.bam.junc", sep="\t", header=F, stringsAsFactors = F)
template = tidyr::separate(template, col="V11", into=c("V13", "V14"), sep=",", remove=F)
template$start = template$V2 + as.integer(template$V13)
template$end = template$V3 - as.integer(template$V14)
rownames(template) = paste(template$V1, template$start, template$end, template$V6, sep=":")

junc_mtx = read.table("SpliPath-test/example_chr1.txt.gz", sep="\t", header=T, stringsAsFactors = F)
sample_list = colnames(junc_mtx)[5:ncol(junc_mtx)]
rownames(junc_mtx) = paste(junc_mtx$chr, junc_mtx$start, junc_mtx$end, junc_mtx$strand, sep=":")

junc_mtx$V13 = sample(6:124, nrow(junc_mtx), replace=T)
junc_mtx$V14 = sample(6:124, nrow(junc_mtx), replace=T)
junc_mtx$V11 = paste(junc_mtx$V13, junc_mtx$V14, sep=",")
junc_mtx$V2 = junc_mtx$start - junc_mtx$V13
junc_mtx$V3 = junc_mtx$end + junc_mtx$V14

junc_mtx$V4 = template$V4[1:nrow(junc_mtx)]
junc_mtx$V9 = template$V9[1:nrow(junc_mtx)]
junc_mtx$V10 = template$V10[1:nrow(junc_mtx)]
junc_mtx$V12 = template$V12[1:nrow(junc_mtx)]

for (s in sample_list){
  s_bed = junc_mtx[, c("chr", paste0("V", 2:4), s, "strand", paste0("V", c(2:3, 9:12)))]
  s_bed = s_bed[s_bed[s] != 0, ]
  write.table(s_bed, sprintf("SpliPath-test/regtools_junction/%s.Aligned.sorted.junc", s), row.names = F, col.names=F, quote=F, sep="\t")
}