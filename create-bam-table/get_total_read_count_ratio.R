source("/Users/yoona96/Dropbox/YoonA_Simon_ChipSeq/rcode/helper_functions.R")

get_total_read_count_ratio <- function(samfile, confile){
  info <- get_info(samfile, confile)
  sam_pileup <- info[1]
  con_pileup <- info[2]
  
  r <- sam_count/con_count
  return(r)
}
