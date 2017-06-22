source("/Users/yoona96/Dropbox/YoonA_Simon_ChipSeq/rcode/helper_functions.R")

get_chrbychr_read_count_ratio <- function(samfile, confile){
  info <- get_info(samfile, confile)
  
  sam_pileup <- info[1]
  con_pileup <- info[2]
  chr <- info[3]
  chr_size <- info[4]
  
  r_matrix <- matrix(ncol = 2, nrow = length(chr_size))
  colnames(r_matrix) <- c("chromosome", "r")
  
  for (i in 1:length(chr)){
    cur_chr <- chr[i]

    cur_chr_sam <- sam_pileup[sam_pileup$seqnames == cur_chr, ]$count
    cur_chr_con <- con_pileup[con_pileup$seqnames == cur_chr, ]$count
    
    sam_cur_count <- sum(cur_chr_sam)
    con_cur_count <- sum(cur_chr_con)
    
    cur_r <- sam_cur_count/con_cur_count
    r_matrix[i, ] <- c(cur_chr, cur_r)
  }
  return(r_matrix)
}
