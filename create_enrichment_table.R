source("/Users/yoona96/Dropbox/YoonA_Simon_ChipSeq/rcode/helper_functions.R")

create_enrichment_table<- function(samfile, confile, overlap = 50, bin = 100, outname, r){
  
  info <- get_info(samfile, confile)
  
  sam_pileup <- info[1]
  con_pileup <- info[2]
  chr <- info[3]
  chr_size <- info[4]
  
  for (i in 1:length(chr)){
    if (is.matrix(r)){
      rate <- r[i,2]
    } else{
      rate <- r
    }
    cur_chr <- chr[i]
    cur_chr_size <- chr_size[i]
    
    chr_bin_size <- cur_chr_size %/% overlap + 1
    sam_con_count <- data.table(matrix(0, ncol = 4, nrow = chr_bin_size))
      
    sam_con_count_with_index <- create_chr_bin_ratio_table(sam_con_count, sam_pileup,
                      con_pileup, cur_chr, cur_chr_size, bin, overlap, rate,
                      cur_index = NULL)
    sam_con_count <- sam_con_count_with_index[1]
    
    colnames(sam_con_count) <- c("range","con_count", "sam_count", "enrichment")
    write.table(sam_con_count, paste(outname, "_",cur_chr), sep = '\t', row.names = FALSE, col.names = FALSE)
  }
}
