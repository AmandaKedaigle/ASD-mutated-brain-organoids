#Getting pMM PPI-weighted distance score between gene sets 
#see Yoon et. al. BMC Genomics 2019, code is based on theirs

library(GScluster)

#download human STRING PPI
#download.file(url='https://github.com/unistbig/GScluster-Data/raw/master/human/string.RData', dest='humanPPI.RData')
load("humanPPI.RData")

#Get indices of proteinSet in PPI
pidx = function(gs) {
  unlist(sapply(gs, function(j) {
           a <- which(j == rownames(string))
           if (length(a)) { a }
         }), use.names = FALSE)
}

#pMM score function
DOP <- function(gs1, gs2, gsi1, gsi2, PPI, Alpha = 1) {
  interact <- function(idx1, idx2, PPI) {
    if (length(idx1) == 0 | length(idx2) == 0) {
      return(0)
    }
    return(sum(PPI[idx1, idx2]))
  }
  
  score <- function(A, B, IA, IB, PPI, Alpha) {
    w <- min(length(A), length(B)) / (length(A) + length(B))
    a <- length(A)
    i <- length(intersect(A, B))
    uniqA <- setdiff(IA, IB)
    uniqB <- setdiff(IB, IA)
    common <- intersect(IA, IB)
    nom <- w * interact(uniqA, common, PPI) + interact(uniqA, uniqB, PPI)
    denom <- w * i + length(setdiff(B, A))
    1 - (i + Alpha * nom / denom) / a
  }
  s1 <- score(gs1, gs2, gsi1, gsi2, PPI, Alpha)
  s2 <- score(gs2, gs1, gsi2, gsi1, PPI, Alpha)
  
  return(min(s1, s2))
}

#load DEP sets
res = readRDS("allSigResults.rds")
gs1 = res[res$exp=="SUV","Gene"]
gs2 = res[res$exp=="CHD8","Gene"]
gs3 = res[res$exp=="PTEN","Gene"]

suv.chd8=DOP(gs1, gs2, pidx(gs1), pidx(gs2), string) #0.93
suv.pten=DOP(gs1, gs3, pidx(gs1), pidx(gs3), string) #0.89
chd8.pten=DOP(gs2, gs3, pidx(gs2), pidx(gs3), string) #0.98

#Get background distribution of scores for proteins detected in organoids
allRes = readRDS("allResults.rds")
scores = data.frame(suv.chd8=NA,suv.pten=NA,chd8.pten=NA)
pairs = list(c("SUV","CHD8"),c("SUV","PTEN"),c("CHD8","PTEN"))
for (i in 1:1000) {
  for (pair in 1:3) {
    gs1.rand = sample(allRes[allRes$exp==pairs[[pair]][[1]],"Gene"], length(gs1))
    gs2.rand = sample(allRes[allRes$exp==pairs[[pair]][[2]],"Gene"], length(gs2))
    score = DOP(gs1.rand, gs2.rand, pidx(gs1.rand), pidx(gs2.rand), string)
    scores[i,pair] = score
  }
}
p1 = sum(scores$suv.chd8<suv.chd8)/1000
p2 = sum(scores$suv.pten<suv.pten)/1000
p3 = sum(scores$chd8.pten<chd8.pten)/1000

scores = melt(scores)
actual = data.frame(variable=c("suv.chd8","suv.pten","chd8.pten"), val=c(suv.chd8, suv.pten, chd8.pten), p = c(p1,p2,p3))
ggplot(scores, aes(x=variable, y=value)) +
  geom_boxplot(aes(fill="Random Protein Sets")) + scale_fill_manual(values="white") +
  geom_point(data=actual, aes(x=variable, y=val,col="Differential Protein Sets"), shape=23, size=4, fill="#ee6677") +
  annotate("text",x=actual$variable, y=rep(1,3), label=paste0("p = ",actual$p)) +
  theme_minimal() +
  labs(y="Distance between Protein Sets", x="") +
  theme(legend.title=element_blank())
ggsave("pMM-distance-boxplots.pdf")
