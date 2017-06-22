require(data.table)
require(Rsamtools)
require(rtracklayer)
require(bamsignals)

create_average_o_gram_file <- function(samfile, confile, gff, annotation = 'ARS', binsize = 100, bin = 50, outname){
  gff_data <- read.table(gff, fill = T)
  gff_ann <- subset(gff_data, V3 == annotation)
  ann_num <- nrow(gff_ann)
  
  info <- get_info(samfile, confile)
  sam_pileup <- info[1]
  con_pileup <- info[2]
  chr <- info[3]
  chr_size <- info[4]
  
  aog_table <- data.table(matrix(0, ncol = bin * 2 + 1, nrow = ann_num))
  colnames(aog_table) <- c(paste("UP", seq(bin, 1, -1)), "CENTRE", paste("DOWN", seq(1, bin, 1)))

  for (i in 1 : ann_num){
      print(paste(annotation, " : Line ", i))
      cur_chr <- gff_ann[i, 1]
      cur_chr_size <- chr_size[as.character(cur_chr)]
      
      up_starts <- seq(gff_ann[i,5], min(cur_chr_size, gff_ann[i,5] + binsize * (bin-1)), binsize)
      down_starts <- seq(gff_ann[i,4], max(1, gff_ann[i,4] - binsize * (bin-1)), -binsize)
      
      cur_chr_sam <- sam_pileup[sam_pileup$seqnames == cur_chr, ]
      cur_chr_con <- con_pileup[con_pileup$seqnames == cur_chr, ]
      sam_position <- cur_chr_sam$pos
      con_position <- cur_chr_con$pos
      
      centre_start <- gff_ann[i,4]
      centre_end <- gff_ann[i,5]
      centre_sam_sub <- sum(cur_chr_sam[sam_position >= centre_start & sam_position < centre_end,]$count)
      centre_con_sub <- sum(cur_chr_con[con_position >= centre_start & con_position < centre_end,]$count)
      
      aog_table[i, bin + 1] <- centre_sam_sub / centre_con_sub 
      for (j in 1:bin){
        print(paste("Examining base: ",up_starts[j], "and", down_starts[j]))
  
        up_end <- up_starts[j] + binsize
        up_sam_sub <- sum(cur_chr_sam[sam_position >= up_starts[j] & sam_position < up_end,]$count)
        up_con_sub <- sum(cur_chr_con[con_position >= up_starts[j] & con_position < up_end,]$count)
        if (up_con_sub != 0){
          aog_table[i, bin + 1 - j] <- up_sam_sub / up_con_sub
        } else{
          aog_table[i, bin + 1 - j] <- NA
        }
        
        down_end <- down_starts[j] - binsize
        down_sam_sub <- sum(cur_chr_sam[sam_position < down_starts[j] & sam_position >= down_end,]$count)
        down_con_sub <- sum(cur_chr_con[con_position < down_starts[j] & con_position >= down_end,]$count)
        if (down_con_sub != 0){
          aog_table[i, bin + 1 +j] <- down_sam_sub / down_con_sub
        } else{
          aog_table[i, bin + 1 +j] <- NA
        }
    }
  }
  write.table(aog_table, paste(outname,"_",annotation), sep = '\t', row.names = FALSE, col.names = TRUE)
  
}

