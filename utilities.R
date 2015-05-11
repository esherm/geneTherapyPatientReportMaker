#' Create a word bubble for genes weighted by abundance.
#' @param gene_abundance, a data frame of two cloumns
#'        gene_name: name of the genes, chr
#'        gene_abundance: weight, integer 
#' @param max_num_words
#' @note it plots a word bubble
#' @note this function plots the top max_num_words words sorted by gene_abundance
#' @examples
#' abstracts <- getAbstracts("22238265")
#' counts <- cleanAbstracts(abstracts)
#' generate_word_bubble(counts) ## error to be expected
#' genecount=data.frame(gene_name=counts$word, gene_abundance=counts$freq)
#' generate_word_bubble(genecount)
#' generate_word_bubble(genecount, max_num_words=10)
generate_word_bubble <- function(gene_abundance, max_num_words=500) {
    
    if(! require('PubMedWordcloud')) stop("Need PubMedWordcloud package")
    
    ## check if gene_abundance has the correct columns
    if ( !all( c("gene_name", "gene_abundance")
              %in% colnames(gene_abundance) ) ) stop(
"Expecting two columns:
     gene_name as character
     gene_abundance as integer")
    
    counts <- data.frame(word=gene_abundance$gene_name,
                         freq=gene_abundance$gene_abundance )
    
    plotWordCloud(counts,
                  scale=c(3,0.5),
                  min.freq=1, max.words=max_num_words, 
                  rot.per = 0, 
                  colors=c(colSets("Set1")[-6],colSets("Paired")))
    
}

