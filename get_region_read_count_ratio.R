source("/Users/yoona96/Dropbox/YoonA_Simon_ChipSeq/rcode/helper_functions.R")


get_region_read_count_ratio <- function(samfile, confile, gff, annotation = 'ARS', distance = 100000){
  
  gff_data <- read.table(gff, fill = T)
  gff_ann <- subset(gff_data, V3 == annotation)
  ann_num <- nrow(gff_ann)
  pup <- PileupParam(max_depth=8000, min_mapq=0, min_base_quality=0)
  
  ranges <- c()
  
  chr <- seqlevels(samfile)
  chr_size <- seqlengths(samfile)
  
  for (i in 1 : ann_num){
    chr <- gff_ann[i, 1]
    start <- max(1,gff_ann[i, 4] - distance)
    end <- min(gff_ann[i,5] + distance, chr_size[as.character(chr)])
    range <- paste(chr,":",start, "-", end,sep="")
    ranges <- c(ranges, range)
  }
  
  if (check.binary("bedtools")) {
    ranges_sort   <- bedr.sort.region(ranges, check.chr = FALSE)
    ranges_merged <- bedr.merge.region(ranges_sort, check.chr = FALSE)
  }
  
  total_sam_count <- 0
  total_con_count <- 0
  
  for (i in 1: length(ranges_merged)){
    range_start <- ranges_merged[i]
    cur_chr <- substr(range_start,0, regexpr(":", range_start)-1)
    
    if (i==1){
      start_pt <- 1
      end_pt <- max(as.integer(substr(range_start, regexpr(":", range_start)+1, regexpr("-", range_start)-1))-1, 1)
      
      which <- GRanges(cur_chr, IRanges(start = start_pt, end = end_pt))
      sbp <- ScanBamParam(which = which)
      total_sam_count <- total_sam_count + sum(pileup(samfile, scanBamParam=sbp, pileUpParam=pup)$count)
      total_con_count <- total_con_count + sum(pileup(confile, scanBamParam=sbp, pileUpParam=pup)$count)
    }
    
    start_pt <- as.integer(substr(range_start, regexpr("-", range_start)+1, nchar(range_start)))+1
    
    if (i == length(ranges_merged)){
      end_pt <- chr_size[as.character(cur_chr)]
      
      which <- GRanges(cur_chr, IRanges(start = start_pt, end = end_pt))
      sbp <- ScanBamParam(which = which)
      total_sam_count <- total_sam_count + sum(pileup(samfile, scanBamParam=sbp, pileUpParam=pup)$count)
      total_con_count <- total_con_count + sum(pileup(confile, scanBamParam=sbp, pileUpParam=pup)$count)
    } else{
      range_end <- ranges_merged[i+1]
      cur_chr_next <- substr(range_end,0, regexpr(":", range_end)-1)
      
      if (cur_chr == cur_chr_next){
        end_pt <- as.integer(substr(range_end, regexpr(":", range_end)+1, regexpr("-", range_end)-1))-1
        
        which <- GRanges(cur_chr, IRanges(start = start_pt, end = end_pt))
        sbp <- ScanBamParam(which = which)
        total_sam_count <- total_sam_count + sum(pileup(samfile, scanBamParam=sbp, pileUpParam=pup)$count)
        total_con_count <- total_con_count + sum(pileup(confile, scanBamParam=sbp, pileUpParam=pup)$count)
      } else{
        end_pt <- chr_size[as.character(cur_chr)]
        
        which <- GRanges(cur_chr, IRanges(start = start_pt, end = end_pt))
        sbp <- ScanBamParam(which = which)
        total_sam_count <- total_sam_count + sum(pileup(samfile, scanBamParam=sbp, pileUpParam=pup)$count)
        total_con_count <- total_con_count + sum(pileup(confile, scanBamParam=sbp, pileUpParam=pup)$count)
        
        start_pt2 <- 1
        end_pt2 <- max(as.integer(substr(range_end, regexpr(":", range_end)+1, regexpr("-", range_end)-1))-1, 1)
        
        which <- GRanges(cur_chr, IRanges(start = start_pt2, end = end_pt2))
        sbp <- ScanBamParam(which = which)
        total_sam_count <- total_sam_count + sum(pileup(samfile, scanBamParam=sbp, pileUpParam=pup)$count)
        total_con_count <- total_con_count + sum(pileup(confile, scanBamParam=sbp, pileUpParam=pup)$count)
      }
    }
  }
  r <- total_sam_count/total_con_count
  return(r)
}
