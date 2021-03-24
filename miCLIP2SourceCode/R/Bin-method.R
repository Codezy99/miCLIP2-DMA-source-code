## Differential Methylation Analysis (DMA) function

## import function and packages
#' @import dplyr
#' @import rtracklayer
#' @import DESeq2
#' @import GenomicFeatures
#' @import GenomicRanges
#' @import S4Vectors
#' @import methods
#'
#' @title Bin-based method
#'
#' @description A function to run the Bin-based differential methylation analysis.
#'
#' @author You Zhou, Kathi Zarnack
#'
#' @param peaks A GRanges object which contains the single nucleotide peaks from the peak calling
#'   result of pureCLIP.
#' @param pathToHTseqCount A path way to all the result of HTSeq-count folder.
#' @param sampleCondition A vector which contains the condition of experiment. For example,
#'   sampleCondition <- c(rep("KO",3), rep("WT",3)).
#' @param txDB If the user has assigned the peaks to the gene, please keep txDB=NA and save the gene
#'   ID of each peaks under the column "gene_id". If the user has not assigned the peaks yet, this
#'   parameter should point to a txDB object of the annotation file.
#' @param condition1 A vector which indicates the colname of the control or WT group in the input
#'   peaks objects. For example, condition1 <- c("read_WT_1", "read_WT_2", "read_WT_3").
#' @param condition2 A vector which indicates the colname of the KO or other treatment in the input
#'   peaks objects. For example, condition2 <- c("read_KO_1", "read_KO_2", "read_KO_3").
#' @param binwidth The bin width of gene expression changes for grouping the peaks.


binBased_DMA <- function(peaks, pathToHTseqCount, sampleCondition, txDB=NA,
                         condition1, condition2, binwidth=0.3){

  ## evaluate the log2 fold change for genes
  sampleFiles <- grep("count",list.files(pathToHTseqCount),value=TRUE)
  sampleCondition <- sampleCondition

  sampleTable <- data.frame(sampleName = sampleFiles,
                            fileName = sampleFiles,
                            condition = sampleCondition)
  ddsHTSeq <- DESeq2::DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                                 directory = pathToHTseqCount,
                                                 design= ~ condition)

  dds <- DESeq2::DESeq(ddsHTSeq)
  res.dds <- lfcShrink(dds, contrast=c("condition","KO","WT"), type='normal')
  res.dds <- res.dds[res.dds$baseMean != 0, ]

  ## Bin-based method
  if (!is.na(txDB)) {
    ## assign peaks
    gene <- genes(txDB)
    gene <- gene[order(width(gene), decreasing = TRUE)]
    o <- findOverlaps(peaks, gene, select = "first")
    peaks$gene_id <- gene$gene_id[o]
  }
  peaks$log2FoldChange <- res.dds$log2FoldChange[match(peaks$gene_id,
    rownames(res.dds))] %>% round(.,4)

  peaks <- peaks[!is.na(peaks$log2FoldChange)]
  res.bin <- data.frame(baseMean = 0, log2FoldChange = 0, lfcSE = 0,
                               stat = 0, pvalue= 0, padj = 0, ID_peaks = "x")
  peaks$ID_peaks <- paste0(seqnames(peaks), ",", start(peaks), ",", end(peaks), ",",
                           strand(peaks))
  i = min(res.dds$log2FoldChange)-0.01
  ## the bin-based DMA
  while(i <= (max(peaks$log2FoldChange) + binwidth)){

    tryCatch({
      temp1 <- peaks[peaks$log2FoldChange >= i & peaks$log2FoldChange < i + binwidth]

      df <- data.frame(mcols(temp1)[colnames(mcols(temp1)) %in% condition1])
      df2 <- data.frame(mcols(temp1)[colnames(mcols(temp1)) %in% condition2])
      df <- cbind(df,df2)
      rownames(df) <- temp1$ID_peaks

      coldata3 <- data.frame(condition = as.factor(c(rep("condition1", length(condition1)),
                                                     rep("condition2", length(condition2)))))
      dds.x <- DESeqDataSetFromMatrix(countData = df,
                                      colData = coldata3,
                                      design = ~condition)

      dds.x <- DESeq(dds.x)
      res.x <- results(dds.x, contrast=c("condition","condition2","condition1"))
      res.x$ID_peaks <- rownames(res.x)
      res.bin <- rbind(res.bin, as.data.frame(res.x))
    },
    warning = function(msg){
      message(paste0(msg))
      return(NA)
    },
    error = function(msg){
      message(paste0(msg,"\n"))
      return(NULL)
    })

    i = i+binwidth
  }
  res.bin <- res.bin[-1,]
  return(res.bin)
}
