#Load Monocle
library(monocle3)
library(Seurat)
library(cowplot)
library(ggridges)

#Load data
setwd("PTEN_mito210_3m")
seur = readRDS("celltypesSeur.rds")
fd = data.frame("gene_short_name" = rownames(GetAssayData(seur)))
rownames(fd) = rownames(GetAssayData(seur))
cds <- new_cell_data_set(GetAssayData(seur,slot = "counts"),
                         cell_metadata = seur@meta.data,
                         gene_metadata = fd)


#Run Monocle
cds <- preprocess_cds(cds, num_dim = 27) #use the same # of top PCs as used for clustering
cds <- reduce_dimension(cds,umap.fast_sgd=TRUE)
cds <- cluster_cells(cds, partition_qval = 0.5) #force few partitions
cds <- learn_graph(cds)


#figures
levels =  c("aRG","Cajal Retzius","Cortical Hem","Cycling Progenitors","IPCs",
            "Newborn Deep Layer Neurons","Newborn PNs","Subcortical","Unknown",
            "CFuPNs","CPNs","oRG","PNs","Immature INs",
            "oRG/Astroglia","Ventral Progenitors", "Cycling Ventral Prog.")
cols = c("#41ae76","#ee8866","#bebada","#bbcc33","#fdb462",
         "#f768a1","#fa9fb5","#77aadd","darkgray",
         "#cc6677","#882255","#225522","#aa4499","#332288",
         "#009988","#5C61FF","#B87ACF")
cols = cols[match(levels(factor(colData(cds)$CellType)),levels)]
p1 = plot_cells(cds, color_cells_by = "CellType", 
                cell_size = 0.5,
                label_cell_groups=FALSE,
                label_leaves=FALSE,
                label_branch_points=FALSE,
                labels_per_group = 0) + scale_color_manual(values=cols) + theme_void()
colData(cds)$clusts_27PCs = factor(colData(cds)$clusts_27PCs)
p2=plot_cells(cds, color_cells_by = "clusts_27PCs",show_trajectory_graph = F,
              label_groups_by_cluster = F,group_label_size = 4,cell_size = 0.5) + theme_void()
cds = order_cells(cds) #manually chosen to be the bottom in cycling progenitors
p3 = plot_cells(cds,show_trajectory_graph = T,
              color_cells_by = "pseudotime",
              label_cell_groups=FALSE,
              label_leaves=FALSE,
              label_branch_points=FALSE,
              graph_label_size=1.5,
              cell_size = 0.5) + theme_void()
p4=plot_cells(cds, color_cells_by = "treat", show_trajectory_graph = F, label_cell_groups = F,cell_size = 0.5) + theme_void()


png("Monocle3.png",res=125,width=2000,height=400)
plot_grid(p1,p2,p3,p4,ncol=4,rel_widths = c(0.3,0.3,0.2,0.2))
dev.off()

saveRDS(cds, "Monocle3-cds.rds")

df = as.data.frame(colData(cds))
df$pseudotime = pseudotime(cds)
df = df[is.finite(df$pseudotime),]
png("Monocle3.ridges.png",res=125,height=160,width=400)
ggplot(df[df$CellType=="CPNs",], aes(x=pseudotime, y=treat)) + 
  geom_density_ridges(aes(fill=treat)) +
  theme_classic()
dev.off()
