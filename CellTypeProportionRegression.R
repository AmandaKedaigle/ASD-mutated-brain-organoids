library(lme4)
library(Seurat)
library(reshape2)


#Load data
seur = readRDS("celltypesSeur.rds")
meta=seur@meta.data[,c("CellType","treat","orig.ident")] #treat contains wtVMut info, and orig.ident contains organoid identity
meta$treat = factor(meta$treat, levels=c("wt","mut")) #set wt as base
meta$orig.ident = factor(paste0(meta$treat,"_",meta$orig.ident))
meta = meta[meta$CellType != "Unknown",] #remove unknown cells

#Change cell type names to be compatiable with regression
meta$CellType = factor(gsub(" ",".",meta$CellType))
meta$CellType = factor(gsub("-",".",meta$CellType))
meta$CellType = factor(gsub("/",".",meta$CellType))

lmm_celltype<-function(meta,form="(1|orig.ident)",celltype="CellType",test="treat") {
  dat = meta
  levs=levels(dat[,celltype])
  ret=c()
  for(i in levs) {
    typ=as.character(i)
    dat[paste("Celltype",typ,sep="_")]=as.numeric(dat[,celltype]==i) #column of 1s for celltype of interest, 0s otherwise
    
    form_cur=as.formula(paste(paste("Celltype",typ,sep="_"),paste(test,form,sep=" + "),sep="~ ")) #model with treatment as fixed effect
    form_cur2=as.formula(paste(paste("Celltype",typ,sep="_"),form,sep="~ ")) #model without treatment included
    print(form_cur)
    res=glmer(form_cur,data=dat,family="binomial")
    print(form_cur2)
    res2=glmer(form_cur2,data=dat,family="binomial")
    vari=paste(test,as.character(levels(dat[,test]))[2],sep="")
    
    #In the first model, get coefficient related to treatment (wt-vs-mut) effect
    coef=fixef(res)[vari]
    OR=exp(coef)
    print("start CI")
    CI=confint(res,parm=vari,method="Wald")
    CI_OR=exp(CI)
    
    #compare two models to get p value for each cell type
    anov=anova(res,res2)
    pv=anov$"Pr(>Chisq)"[2]
    res=as.numeric(c(coef,OR,CI,CI_OR,pv))
    ret=rbind(ret,res)
  }
  ret=data.frame(ret)
  rownames(ret)=as.character(levs)
  colnames(ret)=c("coef","OR","CI_coef_low","CI_coef_high","CI_OR_low","CI_OR_high","pval")
  return(ret)
}

ret = lmm_celltype(meta)
ret$adj.pval = p.adjust(ret$pval, method="BH") #get adjusted p value amount of regressions we ran

save(ret,file="CellComposition.LMM.NoUnk.Robj")
write.table(ret,file="CellComposition.LMM.NoUnk.txt",sep="\t",quote=F,row.names=T)



#For figure, combine results from several datasets
load("PTEN_mito210_1m/CellComposition.LMM.Robj")
ret1m1 = ret
ret1m1$dataset = "PTEN_35div_batch1"
ret1m1$CellType = rownames(ret1m1)
load("PTEN_mito210_1m_r2/CellComposition.LMM.Robj")
ret1m2 = ret
ret1m2$dataset = "PTEN_35div_batch2"
ret1m2$CellType = rownames(ret1m2)
load("PTEN_mito210_3m/CellComposition.LMM.Robj")
ret3m = ret
ret3m$dataset = "PTEN_3mon"
ret3m$CellType = rownames(ret3m)
ret = rbind(ret1m1, ret1m2, ret3m)
ret$CellType = factor(ret$CellType, levels = rev(c("Cycling.Progenitors","aRG","IPCs","Newborn.PNs",
                                               "Newborn.DL.PNs","Subcortical","Cortical.Hem","Cajal.Retzius",
                                               "Newborn.Neurons","oRG","PNs","CFuPNs","CPNs",
                                               "Cycling.IN.Progenitors","IN.Progenitors","Immature.INs"))) #order of y axis
ret$sig = ret$adj.pval<0.05
ggplot(ret,aes(x=dataset,y=CellType,fill=coef,size=-log(adj.pval,10),color=sig)) +
  geom_point(shape=21) + scale_fill_gradient2(low="blue",mid="white",high="red",name="Regression Coeff.") +
  scale_color_manual(values=c("gray","black")) +
  xlab("") + ylab("") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("CellComposition.LMM.allPTEN.sig.pdf")
