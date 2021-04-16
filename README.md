# miCLIP2-DMA-source-code
Author: You Zhou, Kathi Zarnack    
    
A function for running the bin-based differential methylation analysis (DMA)     

## Installation
install.packages("path_to/miCLIP2SourceCode_1.0.1.tar.gz", repos = NULL, 
type="source")    

## Description
In order to learn about the features of genuine m6A sites in the miCLIP2 data, 
we sought to extract all miCLIP2 peaks that significantly changed in the Mettl3 
KO mESCs. However, changes at individual peaks were overshadowed by massive 
shifts in gene expression in Mettl3 KO cells, with more than 2,809 genes 
altered at least 2-fold in comparison to WT mESCs (false discovery 
rate [FDR] ≤ 0.01). These massive shifts in the underlying transcript 
abundances meant that miCLIP2 read counts at individual peaks could not be 
compared directly. In order to overcome this shortcoming, we tested several 
strategies for differential methylation analysis to account for the substantial 
gene expression changes in the Mettl3 KO cells. Best performance was achieved 
with the bin-based approach, here we provide a function to do the 
`bin-based differential methylation analysis` for any miCLIP2 data.

## Workflow for running the binBased_DMA function
The `binBased_DMA` function requires the gene counting result by `htseq-count`  
and the single nucleotide peaks that output by the `pureCLIP` with the 
truncation signal miCLIP2 data as input.   

After user collecting the output from the `htseq-count` and `pureCLIP`, user 
can follow the following workflow to complete the 
`bin-based differential methylation analysis`.

1) Import the bed file of the peaks as a `GRanges` object.
2) Assign the truncation signal to the peaks. For doing that, one option 
could be the function `truncationAssignment` in the package 
[m6Aboost](https://github.com/ZarnackGroup/m6Aboost).
3) Generate the TxDb annotation object for the experiment.
4) Fill the parameters and run the `binBased_DMA` function. 
5) The `binBased_DMA` exports an GRanges object that contains the result of 
the peak changes.

## Example of usage:
path2htseq <- "/path/to/the/htseqresult/"    
conditionS <- c(rep("KO",3), rep("WT",3))    
condition1 <- c("read_WT_1", "read_WT_2", "read_WT_3")    
condition2 <- c("read_KO_1", "read_KO_2", "read_KO_3")    

result <- binBased_DMA(peaks, path2htseq, conditionS, txDB=NA, 
                       condition1=condition1, condition2=condition2)

