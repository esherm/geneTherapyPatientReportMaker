library(stringr)

sanitize <- function(string) {
  result <- gsub("&", "\\&", string, fixed = TRUE)
  result <- gsub("_", "\\_", result, fixed = TRUE)
  result
}

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

#' convert date format to number of days.
#'
#' the format is "^[mdy]\d+\.\d*$" for example: m10, y7, d3.5, m3., m3.0.
#' @param vector of dates to convert
#' @return vector of days
mdy_to_day <- function(dates) {
    stopifnot(check_date_format(dates))
    mdy_letter <- get_mdy_letter(dates)
    mdy_value <- get_mdy_value(dates)
    letter_value <- data.frame(letter=mdy_letter, value=mdy_value)
    convert_to_day(letter_value)
}

#' the format is "^[mdy]\d+\.?\d*$" 
#' @return TRUE if format is correct else FALSE
check_date_format <- function(dates) {
    pattern <- "^[mdy]\\d+\\.?\\d*$"
    all(str_detect(dates, pattern))
}

# only get m, d or y
get_mdy_letter <- function(dates) {
    date_letter <- "^[mdy]"
    str_extract(dates, date_letter)
}

# only get number(ignore m, d, y)
get_mdy_value <- function(dates) {
    date_number <- "\\d+\\.?\\d*$" 
    mdy_value <- str_extract(dates, date_number) 
    as.numeric(as.character(mdy_value))
}

#' lookup table used in convert_to_day f
mdy_to_day_lookup <- data.frame(
    letter=c('d', 'm', 'y'),
    days=c(1, 30.5, 366)
)

#' convert to days
#' @param datafram with 2 cols: letter(d, m, y) and value(numeric).
#' @return days
convert_to_day <- function(letter_value) {
    let_val_days <- merge(letter_value, mdy_to_day_lookup) # days here per unit
    let_val_days$value * let_val_days$days  
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

