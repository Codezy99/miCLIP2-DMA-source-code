# miCLIP2-DMA-source-code
Author: You Zhou, Kathi Zarnack    
    
A function for running the bin-based differential methylation analysis (DMA)     

## Example of usage:

path2htseq <- "/path/to/the/htseqresult/"    
conditionS <- c(rep("KO",3), rep("WT",3))    
condition1 <- c("read_WT_1", "read_WT_2", "read_WT_3")    
condition2 <- c("read_KO_1", "read_KO_2", "read_KO_3")    

result <- binBased_DMA(peaks, path2htseq, conditionS, txDB=NA, condition1=condition1,
                  condition2=condition2)
