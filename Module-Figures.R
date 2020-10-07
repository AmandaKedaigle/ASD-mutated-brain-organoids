library(Seurat)
library(monocle3)
library(reshape2)
library(ggplot2)
library(dplyr)
library(ggbeeswarm)
library(ggpubr)
library(rstatix)
library(lme4)
library(patchwork)

setwd("HUES66_3mon/")

#Load modules
modules = readRDS("modules.CHD8_HUES66_3m_rep2_genes_in_module_pow12.rds")
seur = readRDS("celltypesSeur.rds")

#subset to cell used to calculate modules (one partition of the pseudotime)
cds = readRDS("Monocle3-INsSub-cds.rds")
sub = subset(seur, cells = rownames(colData(cds)))

#subset to equal number of cells per organoid
Idents(sub) = "orig.ident"
sub = subset(sub, cells=WhichCells(sub,downsample = min(table(sub$orig.ident))))
sub = AddModuleScore(sub, modules)
colnames(sub@meta.data)[(ncol(sub@meta.data)-length(modules)+1):ncol(sub@meta.data)] = names(modules)


#Plot distributions (inspired by Super Plots Lord et. al. 2020)
df = sub@meta.data[,c("orig.ident", "treat", names(modules))]
df$orig.ident = factor(df$orig.ident)
df$treat = factor(df$treat, levels=c("wt","mut"))
dfm = melt(df)
ReplicateAverages <- dfm %>% group_by(variable,treat,orig.ident) %>% summarize(mean(value))
colnames(ReplicateAverages)[[4]] = 'value'
ggplot(dfm, aes(x=treat,y=value)) + 
  geom_violin(aes(fill=treat)) + 
  stat_summary(fun=median, geom="crossbar")+
  scale_fill_manual(values=c("#00BFC4","#F8766D"))+
  geom_beeswarm(data=ReplicateAverages, size=3, cex = 4)  + 
  theme_classic() +
  facet_grid(. ~ variable) +
  labs(y="Module Expression Score", x=NULL, color="Organoid Mean Value")


#Calculate significance with mixed model
ret=c()
for (i in 3:(length(modules)+2)) {
  colnames(df)[i] = "col"
  res=lmer(col ~ treat + (1|orig.ident),data=df)
  res2=lmer(col ~ (1|orig.ident),data=df)
  colnames(df)[i] = names(modules)[i-2]
  coef=fixef(res)["treatmut"]
  OR=exp(coef)
  CI=confint(res,parm="treatmut",method="Wald")
  CI_OR=exp(CI)
  anov=anova(res,res2)
  pv=anov$"Pr(>Chisq)"[2]
  resf=as.numeric(c(coef,OR,CI,CI_OR,pv))
  ret=rbind(ret,resf)
}
ret=data.frame(ret)
rownames(ret)=names(modules)
colnames(ret)=c("coef","OR","CI_coef_low","CI_coef_high","CI_OR_low","CI_OR_high","pval")
ret = ret[rownames(ret)!="grey",]
ret$adj.pval = p.adjust(ret$pval, method="BH")

write.table(ret,file="modules-lmm-pvals.txt",sep="\t",quote=F,row.names=T)

#Figure with p values added
pval = data.frame(
  variable = rownames(ret),
  group1 = "wt",
  group2 = "mut",
  p = ret$pval,
  p.adj = signif(ret$adj.pval, digits=2)
)

dfm = dfm[dfm$variable!="grey",]
ReplicateAverages = ReplicateAverages[ReplicateAverages$variable!="grey",]
pval = add_y_position(pval, data=dfm, formula = value~treat)
pval = add_significance(pval)

ggplot(dfm, aes(x=treat,y=value)) + 
  geom_violin(aes(fill=treat)) + 
  stat_summary(fun=median, geom="crossbar")+
  scale_fill_manual(values=c("#00BFC4","#F8766D"))+
  geom_beeswarm(data=ReplicateAverages, size=3, cex = 4)  +
  theme_classic() +
  facet_grid(. ~ variable, labeller = label_wrap_gen(width=15)) +
  labs(y="Module Expression Score", x=NULL, color="Organoid Mean Value") +
  stat_pvalue_manual(pval, label="{p.adj}{p.adj.signif}", hide.ns = T)

#Feature plots - 
#We want the UMAP from the psuedotimes so the shapes match in the figure
#first, need to recalculate cells since analyses were downsampled in different ways
cds = readRDS("monocle3-cds-sub.rds")
sub = subset(seur, cells = rownames(colData(cds))) #downsample to the ones on the UMAP - which were downsampled per treatment instead of per org
sub = AddModuleScore(sub, modules)
colnames(sub@meta.data)[(ncol(sub@meta.data)-length(modules)+1):ncol(sub@meta.data)] = names(modules)
sub@reductions$tsne@cell.embeddings = cds@int_colData@listData$reducedDims@listData$UMAP[rownames(sub@meta.data),]
colnames(sub@reductions$tsne@cell.embeddings) = c("tSNE_1","tSNE_2")

modules = within(modules, rm("grey"))
p = list()
for (m in 1:length(modules)) {
  p[[m]] = FeaturePlot(sub, names(modules)[[m]],pt.size=0.3) + 
    scale_color_gradient2(low="goldenrod",mid="red",high="red4",midpoint=2) + NoAxes()
}
wrap_plots(p)
ggsave("modules-Featureplots.tiff")
