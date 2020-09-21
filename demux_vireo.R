#Run cellSNP on cluster, for each lane
#cellSNP -s outs/possorted_genome_bam.bam -b outs/filtered_gene_bc_matrices/GRCh38/barcodes.tsv -o s3.cellsnp_out.vcf -R common_SNP_files/genome1K.commonSNPsLiftedOver.hg38.vcf.gz -p 10 --minMAF 0.1 --minCOUNT 20

library(cardelino)
library(Seurat, lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library/v2")

cell_data <- load_cellSNP_vcf("orgs1-6.vcf.gz", 
                              max_other_allele = NULL, 
                              min_count = 0, min_MAF = 0)
summary(cell_data)
ids <- vireo(cell_data = cell_data, n_donor = 2)
table(ids$assigned$donor_id)
saveRDS(ids, "orgs1-6.vireo.ids.rds")

seur = readRDS("clusteredSeur.rds")
seur = SetAllIdent(seur, "orig.ident")
#org2 = SubsetData(seur, ident.use =2)

ids_new = ids$assigned
#rownames(ids_new) = paste0("2_",substr(ids$assigned$cell,1,nchar(ids$assigned$cell)-2))

#get cell names - merge vcf has a different way to solve cell name conflicts than Seurat does
cellnames = ids_new$cell
for (index in 1:length(cellnames)) {
  cell = cellnames[index]
  f = substr(cell,1,1)
  if (f %in% c("A","G","C","T")) {
    for (i in 1:6) {
      new = paste0(i, "_", substr(cell,1,nchar(cell)-2))
      if (paste0(gsub("_",":",new),"-1") %in% cellnames) { break }
      if (new %in% seur@cell.names) {
        cellnames[index] = new
        break
      }
    }
  } else {
    new = paste0(f, "_", substr(cell,3,nchar(cell)-2))
    if (new %in% seur@cell.names) {
      cellnames[index] = new
    }
  }
}
rownames(ids_new) = cellnames

ids_new_sub = ids_new[seur@cell.names,] #there are 80k in vireo, only 63k in Seurat due to filtering

seur@meta.data$donor_id = NA
seur@meta.data[rownames(ids_new_sub),"donor_id"] = ids_new_sub$donor_id
save.image("seur.org1-6assignments.RData")

seur = SetAllIdent(seur, "donor_id")
pdf("vireo_org1-6_assignments.pdf")
TSNEPlot(seur, do.label=T)
dev.off()


line1 = SubsetData(seur, cells.use = seur@cell.names[seur@meta.data$donor_id=="donor1"])
saveRDS(line1, "line1.rds")
line2 = SubsetData(seur, cells.use = seur@cell.names[seur@meta.data$donor_id=="donor2"])
saveRDS(line2, "line2.rds")
