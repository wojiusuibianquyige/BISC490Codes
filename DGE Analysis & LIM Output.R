
BiocManager::install("statmod")
library(TCGAbiolinks)
library(limma)
library(edgeR)
library(statmod)
#基因表达矩阵
targets <- read.csv("/Users/yiyang/Desktop/GBM/Files/Rawcounts.csv", check.names = F)[,-1]
rownames(targets) <- targets[,1]
targets <- targets[,-1]
targets <- as.numeric(targets)

#分组
group <- c(rep("NT",5),rep("TP",169))
#构建DGEList对象
dgelist <- DGEList(counts = targets, group = group)
#过滤低表达率基因，进行标准化
keep <- rowSums(cpm(dgelist) > 1 ) >= 2 #过滤
dgelist <- dgelist[keep, ,keep.lib.sizes = FALSE]
#TMM标准化
dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')  #TMM标准化
#绘制MDS
plotMDS(dgelist_norm, col = rep(c('red', 'blue'), each = 5), dim = c(1, 2))
#估算离散值
design <- model.matrix(~group)    #构建分组矩阵
dge <- estimateDisp(dgelist_norm, design, robust = TRUE) #估算离散值
#绘制BCV估算离散值
plotBCV(dge) #作图查看

#基因差异分析
fit <- glmFit(dge, design, robust = TRUE)     #拟合模型
lrt <- glmLRT(fit)   #统计检验
topTags(lrt)

dge_de <- decideTestsDGE(lrt, adjust.method = 'fdr', p.value = 0.05)  #查看默认方法获得的差异基因
summary(dge_de)

plotMD(lrt, status = dge_de, values = c(1, -1), col = c('blue', 'red'))#作图观测
abline(h = c(-1, 1), col = 'gray', lty = 2)

# export the standard Limma output 
v <- voom(dge, design, plot=TRUE)
fit <- lmFit(v, design)
fit <- eBayes(fit)
topTable(fit, coef=ncol(design), n=500)
LIM <- as.data.frame(topTable(fit, coef=ncol(design), n=500))
LIM <- cbind(rownames(LIM),LIM)
colnames(LIM) <- c("ENSEMBL", "logFC", "AveExpre", "T", "Pvalue", "Adj.P.Value","B")

#convert gene ID
ensemblID <- LIM$RefID
gene.df <- bitr(ensemblID, fromType = "ENSEMBL",
                toType = "REFSEQ",
                OrgDb = org.Hs.eg.db)
LIM <- merge(gene.df,LIM,by="ENSEMBL")
LIM <- LIM[,-1]
write.csv(LIM,"/Users/yiyang/Desktop/GBM/Files/LIM.csv")
save.image("/Users/yiyang/Desktop/GBM/GBM/DGE Analysis & LIM Output.Rdata")





















