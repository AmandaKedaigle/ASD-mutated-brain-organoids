library(Seurat)
library(ggplot2)
library(harmony)

#Load individual Seurat object, merge, and run Harmony batch correction
object.names = c("GM_1m_demuxed","GM_3mon","GM_6mon")
object.list = list()
for (n in object.names) {
  seur = readRDS(paste0("indiv_datasets/v3wt.",n,".seur.rds"))
  seur$dataset = n
  object.list[[n]] = seur
}
combined = merge(object.list[[1]], object.list[2:length(object.list)], add.cell.ids=1:length(object.list), merge.data = T)
rm(object.list)
combined = FindVariableFeatures(combined, selection.method='mean.var.plot')
combined <- ScaleData(combined, features=VariableFeatures(combined),vars.to.regress=c("nCount_RNA", "CC.Difference"))
combined <- RunPCA(combined)

combined = RunHarmony(combined, group.by.vars = "dataset", theta=2) #theta=Diversity penalty parameter. Default=2. theta=0 does not encourage any diversity.

combined = RunTSNE(combined, reduction="harmony", dims=1:30)
#combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:30) 
#combined = FindClusters(combined)
#combined$harmonizedClusts_PC30 = Idents(combined)

saveRDS(combined, "DK_GM.harmonized.rds")

#updating & harmonizing Cell Types for this object
combined$CellType = NA
types1m = readRDS("GM_1m_demuxed/wtv3/celltypes6-11-20.rds")
rownames(types1m) = paste0("1_", rownames(types1m))
combined$CellType[combined$dataset=="GM_1m_demuxed"] = 
  types1m[rownames(combined@meta.data[combined$dataset=="GM_1m_demuxed",]),"CellType"]
types3m = readRDS("GM_3mon/wtAndMut/celltypes6-17-20.rds")
rownames(types3m) = paste0("2_", substr(rownames(types3m),4,nchar(rownames(types3m))))
combined$CellType[combined$dataset=="GM_3mon"] = 
  types3m[rownames(combined@meta.data[combined$dataset=="GM_3mon",]),"CellType"]
types6m = readRDS("GM_6mon/wtAndMut/celltypes6-17-20.rds")
rownames(types6m) = paste0("3_", substr(rownames(types6m),4,nchar(rownames(types6m))))
combined$CellType[combined$dataset=="GM_6mon"] = 
  types6m[rownames(combined@meta.data[combined$dataset=="GM_6mon",]),"CellType"]
combined$CellType[combined$CellType=="IP"] = "IPCs"
combined$CellType[combined$CellType=="Unknown"] = "PNs"

combined$broadCellType = combined$CellType
combined$broadCellType[combined$CellType %in% c("CFuPNs","CPNs","Cajal Retzius","Newborn DL neurons","Newborn PN")] = "PNs"
combined$broadCellType[combined$CellType %in% c("Immature Interneurons")] = "INs"
combined$broadCellType[combined$CellType %in% c("aRG","oRG","oRG/Astroglia")] = "RG"
combined$broadCellType[combined$CellType %in% c("Cycling Ventral Prog.")] = "Cycling Progenitors"


saveRDS(combined, "DK_GM.harmonized.rds")


#Enrichment of ASD risk genes from Satterstrom et. al.
TADAgenes = c("CHD8","SCN2A","SYNGAP1","ADNP","FOXP1","POGZ","ARID1B","KMT5B","DYRK1A","SLC6A1",
              "GRIN2B","PTEN","SHANK3","MED13L","GIGYF1","CHD2","ANKRD11","ANK2","ASH1L","TLK2","DNMT3A",
              "DEAF1","CTNNB1","KDM6B","DSCAM","SETD5","KCNQ3","SRPRA","KDM5B","WAC","SHANK2","NRXN1",
              "TBL1XR1","MYT1L","BCL11A","RORB","RAI1","DYNC1H1","DPYSL2","AP2S1","KMT2C","PAX5","MKX",
              "GABRB3","SIN3A","MBD5","MAP1A","STXBP1","CELF4","PHF12","TBR1","PPP2R5D","TM9SF4",
              "PHF21A","PRR12","SKI","ASXL3","SPAST","SMARCC2","TRIP12","CREBBP","TCF4","CACNA1E",
              "GNAI1","TCF20","FOXP2","NSD1","TCF7L2","LDB1","EIF3G","PHF2","KIAA0232","VEZF1","GFAP",
              "IRF2BPL","ZMYND8","SATB1","RFX3","SCN1A","PPP5C","TRIM23","TRAF7","ELAVL3","GRIA2",
              "LRRC4C","CACNA2D3","NUP155","KMT2E","NR3C2","NACC1","PTK7","PPP1R9B","GABRB2","HDLBP",
              "TAOK1","UBR1","TEK","KCNMA1","CORO1A","HECTD4","NCOA1","DIP2A")

combined = AddModuleScore(combined, features=list(TADAgenes))
combined$ASDmoduleScore = combined$Cluster1
combined$Cluster1 = NULL
FeaturePlot(combined, "ASDmoduleScore") + ggtitle("ASD Gene Module Score")

#Enrichment of DEPs
prots = readRDS("../Proteomics/allSigProteins.rds")
combined = AddModuleScore(combined, features=list(prots))
combined$DEPmoduleScore = combined$Cluster1
combined$Cluster1 = NULL
FeaturePlot(combined, "DEPmoduleScore") + ggtitle("DEP Gene Module Score")

#Figures
levels =  c("RG","Cortical Hem","Cycling Progenitors","IPCs","Subcortical","PNs","INs")
cols = c("#41ae76","#bebada","#bbcc33","#fdb462","#77aadd","#aa4499","#332288")
cols = cols[match(levels(factor(combined$broadCellType)),levels)]
DimPlot(combined, group.by = "broadCellType",cols=cols,pt.size = 0.2) + NoAxes() + NoLegend()
ggsave("Fig1b1.tiff")

DimPlot(combined, group.by = "dataset", pt.size=0.2,cols = c("#ddaa33","#bb5566","#004488")) + NoAxes() + NoLegend()
ggsave("Fig1b2.tiff")

cols = cols[match(c("INs","PNs","IPCs","Subcortical","RG","Cycling Progenitors","Cortical Hem"),levels(factor(combined$broadCellType)))]
VlnPlot(combined, "ASDmoduleScore", group.by = "broadCellType", pt.size = 0,sort=T,cols=cols) + 
  ggtitle("ASD Gene Module Score") + NoLegend() + labs(x="")
ggsave("Fig1c.pdf")

cols = cols[match(c("PNs","IPCs","Cycling Progenitors","INs","Subcortical","RG","Cortical Hem"),
                  c("INs","PNs","IPCs","Subcortical","RG","Cycling Progenitors","Cortical Hem"))]
VlnPlot(combined, "DEPmoduleScore", group.by = "broadCellType", pt.size = 0,sort=T,cols=cols) + 
  ggtitle("DEP Gene Module Score") + NoLegend() + labs(x="") + stat_summary(fun=median, geom="crossbar")
ggsave("DEPs.pdf")
