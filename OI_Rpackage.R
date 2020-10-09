library(PCSF)

setwd("~/Documents/DorsalKadoshima/Proteomics/")

data("STRING")
ppi <- construct_interactome(STRING)

prot = readRDS("allSigResults.rds")
terminals = abs(prot$MT.vs.WT)
names(terminals) = prot$Gene
subnet <- PCSF(ppi, terminals, w = 0.8, b = 15, mu = 0.01) #one network
plot(subnet)

subnet <- PCSF_rand(ppi,terminals, n = 10, r = 0.1, w = 0.8, b = 15,mu = 0.01) #multiple networks in parameter sweep
res <- enrichment_analysis(subnet)
plot(res$subnet)
write_graph(res$subnet, "OI_network.gml", format = "gml")
