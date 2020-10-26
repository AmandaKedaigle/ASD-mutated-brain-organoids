#Takes a Seurat Object and performs differential expression analysis 
#By splitting the object into seperate samples (i.e organoids) and using DESeq2,
#which accounts for noise between samples and then looks for DEGs that overcome that noise
#based on method by Sean Simmons

library(Seurat)
library(DESeq2)
library(writexl)
library(readxl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(cowplot)
library(viridis)


#Function that performs the DESeq calculation by summing counts in each sample
combineDE<-function(seur,id,condition="treat",base="wt",combineOn="orig.ident",
                    minCells=20, #Minimum # of cells that must be in this cluster per sample to keep that sample
                    minBatches=2, #Minimum # of samples you can have per condition. Cannot be lower than 2.
                    minReads=10, #Mininum # of reads to have total per gene to calculate DE for that gene
                    genes=c(),  #Genes to consider for DE analysis, if you don't want to use all expressed genes.
                    form="" #design formula for DESeq2 if you want to include more variables in addition to "condition"
                    )
{
  print("Subsample")
  seur<-subset(seur,idents=c(id))
  
  Idents(seur)=combineOn
  genes.use=rownames(GetAssayData(seur,slot="counts"))
  if(length(genes)>0){genes.use=genes}
  
  print("Combine data per sample")
  data.all=data.frame(row.names = genes.use)
  for(i in levels(Idents(seur))) {
    temp.cells=WhichCells(seur,ident=i)
    if (length(temp.cells)==1) data.temp=(GetAssayData(seur,slot="counts")[genes.use,temp.cells])
    if (length(temp.cells)>1) data.temp=apply(GetAssayData(seur,slot="counts")[genes.use,temp.cells],1,sum)
    data.all=cbind(data.all,data.temp)
    colnames(data.all)[ncol(data.all)]=i
  }
  
  print("Filter samples for minimum cells")
  keepOrgs=names(summary(Idents(seur)))[summary(Idents(seur))>minCells]
  numOrg=length(keepOrgs)
  print(paste("Keeping", numOrg, "samples"))
  data.all=data.all[,keepOrgs]

  extraColumns<-strsplit(form,"+",fixed=T)[[1]]
  val=seur@meta.data[,c(condition,combineOn,extraColumns)]
  val=val[!duplicated(val[,2]),]
  rownames(val)=val[,2]
  keepBatch=as.character(val[keepOrgs,1])
  levels = levels(factor(keepBatch))
  if(length(levels)<2) {
    print("Not enough batches per treatment group with minimum # of cells!")
    return(NULL)
  }
  for (level in levels) {
    if(sum(keepBatch==level)<2) {
      print("Not enough batches per treatment group with minimum # of cells!")
      return(NULL)
    }
  }
  
  print("Save meta data")
  colDat=factor(keepBatch)
  if (base != "") { colDat = relevel(colDat, ref=base)}
  colDat=data.frame(colDat)
  colnames(colDat)="condition"
  rownames(colDat)=colnames(data.all)
  colDat[keepOrgs,extraColumns]=val[keepOrgs,extraColumns]
  
  print("Run DESeq2")
  design= ~ condition
  if(nchar(form)>0){
    design=as.formula(paste("~",form,"+ condition",sep=""))
  }
  print(design)
    
  dds <- DESeqDataSetFromMatrix(countData = data.all,colData = colDat,design = design)
  dds <- dds[ rowSums(counts(dds) > minReads)>=2, ]
  dds <- DESeq(dds)
  out=data.frame(results(dds))
  out=out[order(out$pvalue),]
  return(out)
}

#DEseq for pseudotime subsets of cells
#I made a metadata column in each object called "psuedoSub"(sic) with T/F based on whether that cell is included in the trajectory of interest for Fig. 3
for (dataset in c("SUV_mito210_1m_d35/",
                  "PTEN_mito210_3m/",
                  "HUES66_3mon/")) {
  
  #Load Seurat Object
  seur = readRDS(paste0(dataset,"celltypesSeur.rds"))
  
  #Set Ident to cluster or celltypes, whatever groups you want to seperate before detecting DEGs in that group
  Idents(seur) = "psuedoSub"
  
  #Set condition to the metadata column you want DEGs between, and base to the "wildtype" or base level of that column
  condition = "treat"
  base = "wt"
  
  #set combineOn to the metadata column that contains the samples (i.e. different organoids)
  combineOn = "orig.ident"
  
  dir=paste0(dataset,"wtvMut_DEGs_pseudoSub")
  
  #This loop will run DE analysis for each cluster and save a .xlsx file for each!
  id=TRUE
  degs <- combineDE(seur, id=id, condition=condition, base=base, combineOn=combineOn)
  degs$gene = rownames(degs)
  if (length(degs)>0) {
    write_xlsx(degs, path=paste0(dir,".DEGs.xlsx"), format_headers = T)
  }
}

#Read in DEGs for each and find overlaps and GO terms
allGenesUp = list()
allGenesDown = list()
datasets =  c("~/Documents/DorsalKadoshima/SUV_mito210_1m_d35/wtAndMut/",
              "~/Documents/DorsalKadoshima/HUES66_3mon/wtAndMut/",
              "~/Documents/DorsalKadoshima/PTEN_mito210_3m/")
names = c("SUV_d35","CHD8_3m","PTEN_3m")
dat = data.frame()
for (d in 1:length(datasets)) {
  dir=paste0(datasets[[d]],"wtvMut_DEGs_pseudoSub")
  xls = paste0(dir,".DEGs.xlsx")
  if (file.exists(xls)) {
    res = read_xlsx(xls)
    res = res[!is.na(res$pvalue),]
    res$dataset=names[[d]]
    genesUp = res[res$padj<0.05 & res$log2FoldChange>0,"gene"]
    genesDown = res[res$padj<0.05 & res$log2FoldChange<0,"gene"]
    allGenesUp[[names[[d]]]] = genesUp$gene
    allGenesDown[[names[[d]]]] = genesDown$gene
    dat = rbind(dat,res)
  }
}

ckU = compareCluster(allGenesUp, fun="enrichGO", OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="BP")
ckD = compareCluster(allGenesDown, fun="enrichGO", OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="BP")

ckU2 = clusterProfiler::simplify(ckU,cutoff=0.7, by="p.adjust", select_fun=min)
ckD2 = clusterProfiler::simplify(ckD,cutoff=0.7, by="p.adjust", select_fun=min)

dotplot(ckU2) + theme(axis.text.x=element_text(angle=60, hjust=1)) + 
  ggtitle("Upregulated in mutant") + scale_color_viridis(option="D", end=0.8)
ggsave("pseudo-GO.compareDatasets-Up.pdf")

dotplot(ckD2) + theme(axis.text.x=element_text(angle=60, hjust=1)) + 
  ggtitle("Downregulated in mutant") + scale_color_viridis(end=0.8)
ggsave("pseudo-GO.compareDatasets-Down.pdf")

#Get genes that overlap
overlap2 = intersect(allGenesUp$SUV_d35, allGenesUp$CHD8_3m)
overlap2 = union(overlap2, intersect(allGenesDown$SUV_d35, allGenesDown$CHD8_3m))
datSlim = data.frame()
for (gene in overlap2) {
  dats = dat[dat$gene==gene,]
  datSlim = rbind(datSlim, dats)
}
datSlim = datSlim[datSlim$dataset!="PTEN_3m",]
datSlim$dataset = factor(datSlim$dataset, levels=c("SUV_d35","CHD8_3m"))
dss = datSlim[datSlim$dataset=="SUV_d35",]
geneOrder = unique(dss[order(dss$log2FoldChange),]$gene)
datSlim$gene = factor(datSlim$gene, levels=geneOrder)
ggplot(datSlim, aes(x=gene, y=log2FoldChange, fill=log2FoldChange)) +
  geom_bar(stat="identity", position=position_dodge()) +
  facet_grid(dataset ~ ., scales="free") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=0)
ggsave("DEGs.SUVandCHD8overlap-rb.pdf")
