load("../data/mouse_msigdb/mouse_c2_v5p2.rdata")
load("../data/mouse_msigdb/mouse_c3_v5p2.rdata")
load("../data/mouse_msigdb/mouse_c4_v5p2.rdata")
load("../data/mouse_msigdb/mouse_c5_v5p2.rdata")
load("../data/mouse_msigdb/mouse_c6_v5p2.rdata")
load("../data/mouse_msigdb/mouse_c7_v5p2.rdata")
load("../data/mouse_msigdb/mouse_H_v5p2.rdata")

library(piano)
library(annotables)

# Convert ImmGen symbols
tV <- grcm38[["symbol"]]
names(tV) <- as.character(as.numeric(gsub("ENSMUSG", "", grcm38[["ensgene"]])))

immGen <- read.csv("../data/mm10.refGenes.2016.1018.csv")
igenes <- immGen[,3]

fixx <- function(xx){
  return(unname(xx[!is.na(xx) & xx %in% igenes]))
}

mm.c2 <- lapply(Mm.c2, function(pathway) fixx(tV[pathway]))
mm.c3 <- lapply(Mm.c3, function(pathway) fixx(tV[pathway]))
mm.c4 <- lapply(Mm.c4, function(pathway) fixx(tV[pathway]))
mm.c5 <- lapply(Mm.c5, function(pathway) fixx(tV[pathway]))
mm.c6 <- lapply(Mm.c6, function(pathway) fixx(tV[pathway]))
mm.c7 <- lapply(Mm.c7, function(pathway) fixx(tV[pathway]))
mm.H  <- lapply(Mm.H,  function(pathway) fixx(tV[pathway]))

save(mm.c2, mm.c3, mm.c4, mm.c5, mm.c6, mm.c7, mm.H, file = "../output/immgenMSIGDBlists.rda")
