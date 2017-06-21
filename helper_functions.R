require(data.table)
require(Rsamtools)
require(rtracklayer)
require(bamsignals)
library(bedr)

#sam <- "/Users/yoona96/Desktop/bam_file_for_Yoona/Brown_MBS28_TK.bam"
#con <- "/Users/yoona96/Desktop/bam_file_for_Yoona/Brown_MBS24_TK.bam"
#gffFile <- "/Users/yoona96/Dropbox/YoonA_Simon_ChipSeq/Genomic_Features/saccharomyces_cerevisiae_R64-2-1_20150113.gff"

#sambam <- open(BamFile(sam))
#conbam <- open(BamFile(con))

get_info <- function(samfile, confile){
  pup <- PileupParam(max_depth=8000, min_mapq=0, min_base_quality=0)
  
  sbp <- ScanBamParam()
  sam_pileup <- as.data.table(pileup(samfile, scanBamParam=sbp, pileUpParam=pup))
  con_pileup <- as.data.table(pileup(confile, scanBamParam=sbp, pileUpParam=pup))
  
  chr <- seqlevels(samfile)
  chr_size <- seqlengths(samfile)
  
  info <- c(sam_pileup, con_pileup, chr, chr_size)
  
  return (info)
}

plot_bin_coverage <- function(sam_con_count, outname){
  lm(as.numeric(unlist(sam_con_count[, 1])) ~ 
            as.numeric(unlist(sam_con_count[, 2])) +
            I(as.numeric(unlist(sam_con_count[, 2]))^2))

  min_count <- max(log10(min(sam_con_count)),1)
  max_count <- log10(max(sam_con_count))
  sam_con_count[sam_con_count == 0] <- 1

  pdf(outname)
  plot(log10(as.numeric(unlist(sam_con_count[, 2]))), 
     log10(as.numeric(unlist(sam_con_count[, 1]))),
     xlab = "Control read count", ylab = "Sample read count",
     xlim = c(min_count, max_count), ylim = c(min_count, max_count))
  lines(sort(log10(as.numeric(unlist(sam_con_count[, 2])))),
      log10(predict(fit))[order(log10(as.numeric(unlist(sam_con_count[, 2]))))], 
      col = "red", lwd=2)
  abline(lm(log10(as.numeric(unlist(sam_con_count[, 1])))
          ~ log10(as.numeric(unlist(sam_con_count[, 2])))), col = "blue")
  dev.off()
}




create_chr_bin_ratio_table <- function(sam_con_count,sam_pileup, con_pileup, 
                      cur_chr, cur_chr_size, bin, overlap = NULL, 
                      rate = NULL, cur_index = NULL){
  
  print(paste(cur_chr,cur_chr_size))
  cur_chr_sam <- sam_pileup[sam_pileup$seqnames == cur_chr, ]
  cur_chr_con <- con_pileup[con_pileup$seqnames == cur_chr, ]
  if (is.null(overlap)){
    starts <- seq(1, cur_chr_size, bin)
  } else{
    starts <- seq(1, cur_chr_size, overlap)
  }
  sam_position <- cur_chr_sam$pos
  con_position <- cur_chr_con$pos
  for (j in 1:length(starts)){
    if (j %% 100 == 0){
      print(paste("Examining base: ",starts[j]))
    }
    end <- starts[j] + bin
    sam_sub <- cur_chr_sam[sam_position >= starts[j] & sam_position < end,]$count
    con_sub <- cur_chr_con[con_position >= starts[j] & con_position < end,]$count/rate
    
    if (ncol(sam_con_count) == 4){
      sam_con_count[i, 1] <- starts[j] -1 
      sam_con_count[i, 3] <- sum(sam_sub)
      sam_con_count[i, 2] <- sum(con_sub)
      sam_con_count[i, 4] <- sam_con_count[cur_index, 3] / sam_con_count[i, 2]
    } else if (ncol(sam_con_count) == 2){
      sam_con_count[cur_index, 1] <- sum(sam_sub)
      sam_con_count[cur_index, 2] <- sum(con_sub)
      cur_index <- cur_index + 1
    }
  }
  return (sam_con_count, cur_index)
}
