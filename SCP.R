#Code for uploading data to single cell portal

library(Seurat)
library(Matrix)

seur = readRDS("DK_GM.harmonized.rds")
sampleName = "GM_wt_harmonized"
savedir = "portal"
dir.create(savedir)

keep=c("nFeature_RNA","nCount_RNA","percent.mito","percent.ribo","dataset","broadCellType","organoid","treat")
line1="NAME\tnGene\tnUMI\tpercentMito\tpercentRibo\tdataset\tCellType\tOrganoid\tCondition"
line2="TYPE\tnumeric\tnumeric\tnumeric\tnumeric\tgroup\tgroup\tgroup\tgroup"	

meta=seur@meta.data[,keep]
saveMeta=paste(savedir,"/meta_all.txt",sep="")
write(c(line1,line2),saveMeta,sep="\n")
write.table(meta,saveMeta,sep="\t",col.names=F,quote=F,append=T)

tsn=seur@reductions$tsne@cell.embeddings
saveTSNE=paste(savedir,"/tsne_",sampleName,".txt",sep="")
line1="NAME	X	Y"
line2="TYPE	numeric	numeric"
write(c(line1,line2),saveTSNE,sep="\n")
write.table(tsn,saveTSNE,sep="\t",col.names=F,quote=F,append=T)


saveData=paste(savedir,"/expression_",sampleName,".txt",sep="")
writeMM(GetAssayData(seur, slot = "data"),saveData)
genes = rownames(GetAssayData(seur, slot = "data"))
barcodes = colnames(GetAssayData(seur, slot = "data"))
write.table(genes,paste0(savedir,"/expression_",sampleName,"_genes.txt"),col.names=F,row.names=F,quote=F)
write.table(barcodes, paste0(savedir,"/expression_",sampleName,"_barcodes.txt"),col.names=F,row.names=F,quote=F)

#python ~/Desktop/UsefulCode/SortSparseMatrix.py expression_Name.txt
