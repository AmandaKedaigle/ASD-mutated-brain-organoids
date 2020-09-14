#Amanda Kedaigle 7/12/18
#grown from Sean's dir10X function in load_Seurat.R

library(Seurat)
library(stringr)
library(Matrix)
library(writexl)

dir = '../round2_data/'
outdir = '../round2_results/'
dir.create(outdir)

#loads all 10X lanes from a given directory
dir10X<-function(dir="",outdir="",minGenes=500,regress=c("nCount_RNA"))
{
  print(paste("ls ",dir,"/[1-3]/*/outs/filt* | grep : | sed 's/://g'",sep=""))
  lst=system(paste("ls ",dir,"/[1-3]/*/outs/filt* | grep : | sed 's/://g'",sep=""),intern=T)
  print(lst)
  
  print("Read in!")
  dat=Read10X(lst)
  print(head(dat))
  print("Fix colnames!")
  cols=colnames(dat)
  cols_new=c()
  for(col in cols){
    start=str_sub(col,1,1)
    cur=col
    if(start %in% c("A","T","G","C")){cur=paste("1_",cur,sep="")}
    cols_new<-c(cols_new,cur)
  }
  colnames(dat)=cols_new

  print(paste("Dims: ",toString(dim(dat))))
  
  print("Make object!")
  seur<-CreateSeuratObject(dat,"Seurat",min.features=minGenes)
  
  # The % of UMI mapping to MT-genes is a common QC metric.
  mito.genes <- grep(pattern = "^MT-", x = rownames(x = seur), value = TRUE)
  percent.mito <- Matrix::colSums(GetAssayData(seur, slot='counts')[mito.genes, ]) / Matrix::colSums(GetAssayData(seur, slot="counts"))
  seur[['percent.mito']] <- percent.mito
  ribo.genes <- grep("^RP[S,L]",rownames(seur), value = TRUE)
  percent.ribo <- Matrix::colSums(GetAssayData(seur, slot='counts')[ribo.genes, ]) / Matrix::colSums(GetAssayData(seur, slot="counts"))
  seur[['percent.ribo']] <- percent.ribo
  pdf(paste0(outdir,'QC_Rplots.pdf'))
  par(mfrow = c(1, 3))
  print(FeatureScatter(object = seur, feature1 = "nCount_RNA", feature2 = "percent.mito"))
  print(FeatureScatter(object = seur, feature1 = "nCount_RNA", feature2 = "percent.ribo"))
  print(FeatureScatter(object = seur, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
  dev.off()
  
  print("Normalize data!")
  seur = NormalizeData(seur, normalization.method = "LogNormalize", scale.factor=1000000)
  
  print("Get variable genes!")
  seur<-FindVariableFeatures(seur,selection.method='mean.var.plot')
  
  print("Regress out!")
  seur<-ScaleData(seur,features=VariableFeatures(seur),vars.to.regress=regress)
  
  print("Run PCA!")
  seur<-RunPCA(seur,features = VariableFeatures(seur),verbose=F)
  
  pdf(paste0(outdir,'pca_Rplots.pdf'))
  print(VizDimLoadings(object = seur, dims = 1:2))
  print(DimPlot(object = seur, dim.1 = 1, dim.2 = 2))
  DimHeatmap(object = seur, dims = 1:9, cells = 500, balanced = TRUE)
  print(ElbowPlot(object = seur, 30))
  dev.off()
}

dir10X(dir=dir, outdir=outdir)

#pause here - pick an appropriate numPCs from pca_Rplots.pdf ElbowPlot
numPCs=25

#Cluster & run dim reduction using numPCs
seur = FindNeighbors(seur, dims=1:numPCs)
seur = FindClusters(seur, resolution=1.0)
seur[[paste0("clusts_",numPCs,"PCs")]] = Idents(seur)
seur <- RunTSNE(object = seur, dims = 1:numPCs)

saveRDS(seur, "clustered_seur.rds")

#Get marker genes for each cluster, to be used in assigning cell types
Idents(seur) <- paste0("clusts_",numPCs,"PCs")
dir.create("markerListsExcels")

markers = Signac::VeniceAllMarkers(seur, only.pos=T)
write.table(markers, file="markerListsVenice.txt")
for (i in levels(markers$cluster)){
  write_xlsx(markers[markers$cluster==i,], path=paste0("markerListsExcels/Cluster",i,"_markerLists.xlsx"))
}
