source("/Users/yoona96/Dropbox/YoonA_Simon_ChipSeq/rcode/helper_functions.R")

get_bin_coverage<- function(samfile, confile, bin = 100, outname){
  
  info <- get_info(samfile, confile)
  
  sam_pileup <- info[1]
  con_pileup <- info[2]
  chr <- info[3]
  chr_size <- info[4]
  
  total_bin_size <- 0
  for (i in 1: length(chr_size)){
    total_bin_size <- total_bin_size + chr_size[i] %/% 100 + 1
  }
  sam_con_count <- data.table(matrix(0, ncol = 2, nrow = total_bin_size))
  colnames(sam_con_count) <- c("sam_count", "con_count")
  cur_index <- 1
  
  for (i in 1:length(chr)){
    cur_chr <- chr[i]
    cur_chr_size <- chr_size[i]
    
    sam_con_count_with_index <- create_chr_bin_ratio_table(sam_con_count, sam_pileup,
                      con_pileup, cur_chr, cur_chr_size, bin, overlap = NULL,
                      rate = 1, cur_index = cur_index)
    sam_con_count <- sam_con_count_with_index[1]
    cur_index <- sam_con_count_with_index[2]
    
  }
  plot_bin_coverage(sam_con_count, outname)
}
