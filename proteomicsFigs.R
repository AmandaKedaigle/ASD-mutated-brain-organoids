#Amanda Kedaigle 6/22/19

library(readxl)
library(reshape2)
library(ggplot2)
library(patchwork)
library(ggpubr)

setwd("Proteomics/")

res = as.data.frame(read_xlsx("PTEN_2timepts/Results.organoids.update.0524.xlsx",sheet = 3))
colnames(res)[1:2] = c("Gene", "MT.vs.WT")
res$`Significantly changed`= NULL
res$exp = "PTEN"
resTemp = read.table("CHD8.Results.txt", header = T)
colnames(resTemp)[2] = "MT.vs.WT"
resTemp$exp = "CHD8"
res = rbind(res, resTemp)
resTemp = read_xlsx("SUV.strict.Isoform.Results.MT.vs.WT.100820.xlsx", sheet=1)
colnames(resTemp)[1] = "MT.vs.WT"
colnames(resTemp)[6] = "Gene"
resTemp$Significance = NULL
resTemp$exp="SUV"
res = rbind(res,resTemp)

rm(resTemp)


res$log.adj.pval = -1*log10(res$adj.P.Val)
protSig = unique(res$Gene[res$adj.P.Val<0.1])
saveRDS(res[res$adj.P.Val<0.1,], "allSigResults.rds")
saveRDS(res, "allResults.rds")


#PTEN 1m correlation to 2m ------

res = as.data.frame(read_xlsx("PTEN_2timepts/Results.organoids.update.0524.xlsx",sheet = 1))
res$exp = "WT"
colnames(res)[2] = "logFC"
resTemp = as.data.frame(read_xlsx("PTEN_2timepts/Results.organoids.update.0524.xlsx",sheet = 2))
resTemp$exp = "MT"
colnames(resTemp)[2] = "logFC"
res = rbind(res, resTemp)
rm(resTemp)
res$exp = factor(res$exp, levels=c("WT","MT"))
res$abslogFC = abs(res$logFC)

ggplot(res, aes(x=exp, y=abslogFC, fill=exp)) +
  geom_violin() + theme_minimal() +
  stat_compare_means(paired=T,label.y = 3.5,label.x = 1.3) +
  stat_summary(fun=median, geom="crossbar") +
  ylab("Absolute Value of logFC Day 35 vs. 70") + xlab(NULL) + theme(legend.position = "none")
ggsave("PTENproteomics-1to2m-allProts-paired.pdf")

resSub = res[res$Protein %in% protSig,]
ggpaired(resSub, x="exp", y = "abslogFC", id="Protein",
         color="exp", line.color="gray", line.size=0.1, palette=c("#00BFC4","#F8766D")) +
  stat_compare_means(paired=T,label.y = 2.3,label.x = 1.3) +
  ylab("Absolute Value of logFC Day 35 vs. 70") + xlab(NULL) + theme(legend.position = "none")
ggsave("PTENproteomics-1to2m-allSigProts-paired.pdf")

#PTEN logFCs WT to MT at 35days versus WT 35 to 70days -------
setwd("PTEN_2timepts/")
data = list()
data[["WT70v35"]] = as.data.frame(read_xlsx("Results.organoids.update.0524.xlsx",sheet = 1))
rownames(data[["WT70v35"]]) = data[["WT70v35"]]$Protein
data[["MT70v35"]] = as.data.frame(read_xlsx("Results.organoids.update.0524.xlsx",sheet = 2))
rownames(data[["MT70v35"]]) = data[["MT70v35"]]$Protein
data[["MTvWT35"]] = as.data.frame(read_xlsx("Results.organoids.update.0524.xlsx",sheet = 3))
rownames(data[["MTvWT35"]]) = data[["MTvWT35"]]$Protein
data[["MTvWT70"]] = as.data.frame(read_xlsx("Results.organoids.update.0524.xlsx",sheet = 4))
rownames(data[["MTvWT70"]]) = data[["MTvWT70"]]$Protein

proteins = data[["MTvWT35"]]$Protein[data[["MTvWT35"]]$`Significantly changed`]
dat = data.frame("Protein"= proteins)
dat$MTvWT35 = data[["MTvWT35"]][proteins, "MT_35d.vs.WT_35d"]
dat$MTvWT70 = data[["MTvWT70"]][proteins, "MT_70d.vs.WT_70d"]
dat$WT70v35 = data[["WT70v35"]][proteins, "WT_70d.vs.WT_35d"]
dat$MT70v35 = data[["MT70v35"]][proteins, "MT_70d.vs.MT_35d"]

protOrder = dat[order(dat$MTvWT35),"Protein"]

dat = melt(dat)
pval = apply(dat,1,function(x) data[[x["variable"]]][x["Protein"],"adj.P.Val"])
dat$pval = -1*log10(pval)
dat$Protein = factor(dat$Protein, levels=protOrder)

#diverging bars
p1=ggplot(dat[dat$variable %in% c("MTvWT35"),], aes(x=Protein, y=value, fill=value)) + 
  geom_bar(stat="identity") +
  #scale_fill_gradientn(values=c(0,0.125,1), colors=c("gray","pink","red"),limits=c(0,8),name="-log(adj. p value)") +
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=0,name="Fold Change", limits=c(-2,2.5)) +
  coord_flip() +
  theme_classic() +
  theme(legend.position = "bottom")+
  xlab(NULL) + ylab("FC MT v WT at 35d")
p2= ggplot(dat[dat$variable %in% c("WT70v35"),], aes(x=Protein, y=value,fill=value)) +
  geom_bar(stat="identity") +
  #scale_fill_gradientn(values=c(0,0.125,1), colors=c("gray","pink","red"),limits=c(0,8),name="-log(adj. p value)") +
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=0,name="Fold Change",limits=c(-2,2.5)) +
  coord_flip() +
  theme_classic() +
  theme(legend.position = "none",axis.text.y = element_blank(), 
        axis.line.y=element_blank(), axis.ticks.y=element_blank()) +
  xlab(NULL) + ylab("FC 35d v 70d in WT")

p1 + p2 + plot_layout(guides="collect") & theme(legend.position = "bottom")
ggsave("PTENproteomicsBarPlot.pdf")

#PTEN correlations mutant effect vs time effect --------
#From Kalliopi Tsafou

results <- read.csv("Organoid_proteomics/figures/organoids.normalized.imputed.0624.csv")
results$LogFC.MT.vs.WT.d35 <- results$MT_35-results$WT_35
results$LogFC.d70.vs.d35 <- results$WT_70-results$WT_35

#Here we make the plot that shows the correlation with the non-significant PTEN+/- /WT
subresults <- results[which(results$Significant.MT.vs.WT.35d=="FALSE"),]
correlation = stats::cor(subresults$LogFC.MT.vs.WT.d35, subresults$LogFC.d70.vs.d35)
c1 <- cor.test(subresults$LogFC.MT.vs.WT.d35, subresults$LogFC.d70.vs.d35)
corr.index <- paste("r =", signif(correlation, digits = 4))
p1 <- ggplot(subresults, aes(x=LogFC.MT.vs.WT.d35, y=LogFC.d70.vs.d35)) + 
  geom_abline(color = "grey50", linetype = 1)+
  geom_point(alpha = 0.3, size=3, shape=16, fill="grey") + 
  stat_cor(method = "pearson", label.x = -3, label.y = 4) +
  #  geom_smooth(method='lm') +
  scale_x_continuous(limits = c(-3,4), name = expression(MT~vs.~WT~log[2](FC)))+
  scale_y_continuous(limits = c(-3,4), name = expression(d70~vs.~d35~log[2](FC)))+
  ggtitle(corr.index)+
  theme_bw()+
  theme(aspect.ratio=4/4)+
  theme(axis.text.x = element_text(color="black", size=10)) +
  labs(x = "LogFC.MT.vs.WT.d35", y = "LogFC.d70.vs.d35")+
  labs(title = "PTEN+/- /WT\nNon significant")
p1

#Here we make the plot that shows the correlation with the significant PTEN+/- /WT
subresults <- results[which(results$Significant.MT.vs.WT.35d=="TRUE"),]
write.table(subresults, file="organoid.results.pten.maturation.tsv", sep="\t'", quote=F, row.names = F)
correlation = stats::cor(subresults$LogFC.MT.vs.WT.d35, subresults$LogFC.d70.vs.d35)
c2 <- cor.test(subresults$LogFC.MT.vs.WT.d35, subresults$LogFC.d70.vs.d35)
corr.index <- paste("r =", signif(correlation, digits = 4))
p2 <- ggplot(subresults, aes(x=LogFC.MT.vs.WT.d35, y=LogFC.d70.vs.d35)) + 
  geom_abline(color = "grey50", linetype = 1)+
  geom_point(alpha = 0.3, size=3, shape=16, fill="grey") + 
  #  geom_smooth(method='lm') +
  stat_cor(method = "pearson", label.x = -3, label.y = 4) +
  scale_x_continuous(limits = c(-3,4), name = expression(MT~vs.~WT~log[2](FC)))+
  scale_y_continuous(limits = c(-3,4), name = expression(d70~vs.~d35~log[2](FC)))+
  ggtitle(corr.index)+
  theme_bw()+
  theme(aspect.ratio=4/4)+
  theme(axis.text.x = element_text(color="black", size=10)) +
  labs(x = "LogFC.MT.vs.WT.d35", y = "LogFC.d70.vs.d35")+
  labs(title = "PTEN+/- /WT\nSignificant")
p2

#we make the dataframe with the cor. estimates for the PTEN+/- /WT significant or non significant sets
#and make the respective plot
df <- data.frame(Set=c("PTEN-/- /WT\nNon significant", "PTEN+/- /WT\nSignificant"),
                 cor=c(c1$estimate, c2$estimate),
                 ymin=c(c1$conf.int[1], c2$conf.int[1]),
                 ymax=c(c1$conf.int[2], c2$conf.int[2]))
p3 <- ggplot(data=df) +
  geom_pointrange(aes(x=Set, y=cor, ymin=ymin, ymax=ymax)) +
  theme_bw()+
  theme(aspect.ratio=4/2)+
  labs(x = "", y = "Correlation estimates")+
  labs(title = "Correlation") +
  theme(axis.text.x = element_text(color="black", size=10, angle = 45, hjust = 1))
p3

p.all <- p1|p2|p3

p.all <- p.all +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 14), 
    axis.text = element_text(size = 12),
    axis.title = element_text(color="black", size=12),
    axis.text.y = element_text(color="black", size=10)
  )
p.all

ggsave("correlationPlot.pdf", plot=p.all, width=10)
