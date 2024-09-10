rm(list=ls())

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/TCGA")

library(OmicCircos)

data("TCGA.BC.cnv.2k.60")
data("TCGA.BC.gene.exp.2k.60")
data("TCGA.BC.sample60")
data("TCGA.BC_Her2_cnv_exp")

pvalue <- -1*log10(TCGA.BC_Her2_cnv_exp[,5])
pvalue <- cbind(TCGA.BC_Her2_cnv_exp[,c(1:3)], pvalue)

Her2.i <- which(TCGA.BC.sample60[,2] == "Her2")
Her2.n <- TCGA.BC.sample60[Her2.i, 1] 
Her2.j <- which(colnames(TCGA.BC.cnv.2k.60) %in% Her2.n) 
cnv <- TCGA.BC.cnv.2k.60[,c(1:3, Her2.j)]
cnv.m <- cnv[,c(4:ncol(cnv))]
cnv.m[cnv.m > 2] <- 2                 
cnv.m [cnv.m < -2] <- -2
cnv <- cbind(cnv[,1:3], cnv.m ) 

Her2.j <- which(colnames(TCGA.BC.gene.exp.2k.60) %in% Her2.n)
gene.exp <- TCGA.BC.gene.exp.2k.60[,c( 1:3, Her2.j) ] 

colors <- rainbow(10, alpha = 0.5)
par(mar = c(2, 2, 2, 2)) 

plot(c(1, 800), c(1, 800), type ="n", axes = FALSE, xlab = "", ylab = "", main = "")

gene.exp <- top500
cnv <- luad_cnv_circos

circos(R=400, cir = "hg18", W=4, type = "chr", print.chr.lab = TRUE, scale = TRUE)
circos(R=300, cir = "hg18", W=100, mapping = gene.exp, col.v = 4, type = "heatmap2", cluster = T, col.bar = TRUE , lwd = 0.1, col = "blue")
circos(R=220, cir = "hg18", W=80, mapping = cnv, col.v = 4, type = "ml3", B=FALSE, lwd =1, cutoff = 0) 
circos(R=140, cir = "hg18", W=80, mapping = pvalue, col.v = 4, type = "l", B=TRUE, lwd=1, col = colors[1])

# Prepare input data #
sel_sample <- luad_cnv_seg[luad_cnv_seg$sample %in% colnames(luad_cnv),]
sel_sample$sample %>% unique()
x <- luad_cnv_seg$chr
luad_cnv_seg$chr <- paste("chr",x, sep = "")
s1 <- sel_sample[ (sel_sample$sample %in% c("TCGA-44-6147")),]
luad_seg_list <- list(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12)
save(luad_seg_list, file = "relationship_plot/luad_seg_list.Rdata")

# Select most variable genes #
rna_log_normalized <- luad_rna %>% as.matrix() %>% DESeq2::rlog()
topVarGenes <- head(order(matrixStats::rowVars(rna_log_normalized), decreasing = TRUE), 300)
top500 <- rna_log_normalized[topVarGenes, ]
rna_zscore <- t(scale(t(top500)))
top500 <- rna_zscore %>% as.data.frame() %>% dplyr::select(13:24)
luad_cnv_circos <- luad_cnv_circos[luad_cnv_circos$NAME %in% rownames(top500),]
s2_cnv_anno <- s2_cnv_anno[s2_cnv_anno$NAME %in% s1_cnv_anno$NAME,]
top500 <- top500[top500$NAME %in% luad_cnv_circos$NAME, ]

# Formatting cnv table #
s2_cnv_anno <- setNames(s2_cnv_anno, c("chr","start","po","TCGA-50-6595","sample","NAME")) %>% dplyr::select(1,3,4,6)
s2_cnv_anno <-s2_cnv_anno %>% dplyr::select(3,4)

luad_cnv_circos$chr <- stringr::str_sub(luad_cnv_circos$chr,4,5)
colnames(top500) <- stringr::str_sub(colnames(top500),1,12)

rownames_idx <- match(luad_cnv_circos$NAME, rownames(top500))
top500 <- top500[rownames_idx,]
save(luad_cnv_circos, file = "circos_plot/luad_cnv_circos.Rdata")


