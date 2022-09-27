#!/apps/R/gnu/9.1/3.6.3/bin/Rscript

options(warn = -1)

suppressMessages( library(ggplot2) )
suppressMessages( library(dplyr) )
suppressMessages( library(DESeq2) )

# read in the data 
argv = commandArgs(TRUE)

# counts
counts = read.table(argv[1], header = TRUE)
outdir = dirname(argv[1])

# coldata
coldata = read.table(argv[2], header = TRUE)
row.names(coldata) = coldata$sample
coldata = subset(coldata, select = c("condition"))

# comparisons
comparisons = read.table(argv[3], header = FALSE)
names(comparisons) = c("x", "y")

# outname
outname = argv[4]

# select cols that are non-numeric/numeric
tmp = data.frame(counts, row.names = 1)
names_df = select_if(tmp, Negate(is.numeric))
cts = select_if(tmp, is.numeric)

# remove rows that have no counts
x = rowSums(cts) >= 1
cts_filt = cts[x,]

# set condition equal to coldata condition
condition = coldata$condition

# DESeq
deobj <- DESeqDataSetFromMatrix(countData = cts_filt, colData = coldata, design = ~condition)
dds <- DESeq(deobj)
dds_counts = counts(dds, normalized = TRUE)

df <- cbind(x = rownames(dds_counts), dds_counts)
names(df)[names(df) == 'x'] = paste0(names_df[1])
rownames(df) <- 1:nrow(df)

write.table( merge(names_df, df, by = paste0(names_df[1]), all.y = TRUE), paste0(outname,"-normalized-counts.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# PCA plot to assess variance within sample groups and between sample groups 
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

PCA = ggplot(pcaData, aes(PC1, PC2, color=name)) +
  geom_point(size=5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme(aspect.ratio = 1) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

pdf(paste0(outname,"-PCA.pdf"), height = 8, width = 8)
PCA
dev.off() 

#results and plot for each comparison
for(i in 1:nrow(comparisons)){

  sample_x = as.character(comparisons[i,1])
  sample_y = as.character(comparisons[i, 2])

  curr = paste0(sample_x,"-vs-",sample_y)

  xdf = subset(coldata, condition == paste0(sample_x))
  xnames = rownames(xdf)

  ydf = subset(coldata, condition == paste0(sample_y))
  ynames = rownames(ydf)

  res = results(dds, contrast = c("condition",paste0(sample_y),paste0(sample_x)))

  c = as.data.frame(counts(dds,normalized = TRUE))
  keep = as.vector(rbind(paste0(xnames), paste0(ynames)))
  c_sub = subset(c, select = keep)

  resdata = merge(as.data.frame(res), c_sub, by = 'row.names', sort = FALSE)
  names(resdata)[1] <- paste0(names_df[1])

  resdata[paste0(sample_x)] = rowMeans( resdata[ , xnames] )
  resdata[paste0(sample_y)] = rowMeans( resdata[ , ynames] )
  
  resdata = merge(resdata, geneids, by = paste0(names_df[1]), all.x = TRUE)

  deseq_cols = c("baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")
  cols = c(names_df, deseq_cols)

  ord = c(cols, deseq_cols, list(xnames)[[1]], list(ynames)[[1]], sample_x, sample_y)
  write.table(resdata[,ord], paste0(outname,"-",curr,"-deseq-results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

}