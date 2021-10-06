library(GenomicInteractions)
library(GenomicRanges)
library(tidyverse)
library(fitdistrplus)
options(scipen = 999)

peak2gr <- function(dat, col_names=T, delim="\t", comment="#",keep.extra.columns = F){
  bed <- readr::read_delim(dat, 
                           col_names = col_names, 
                           delim=delim,
                           comment = comment)
  colnames(bed)[1:3] <- c("chr","start","end")
  gr <- GenomicRanges::makeGRangesFromDataFrame(bed, keep.extra.columns = keep.extra.columns)
  return(gr)
}

get_loop_bg <- function(loop){
  #prepare distal interaction
  con <- isInteractionType(loop, "distal", "enhancer") |
        isInteractionType(loop, "distal", "promoter") |
    isInteractionType(loop, "distal", "ctcf")
  
  #background
  bg <- loop[con]
  bg$dis <- calculateDistances(bg)
  return(bg)
}

get_element_loop <- function(loop){
  #only get loops between ctcf, promoter, enhancer
  con <- (isInteractionType(loop, "enhancer","enhancer")) | (isInteractionType(loop, "enhancer","promoter")) | 
        (isInteractionType(loop, "enhancer","ctcf")) | (isInteractionType(loop, "promoter","promoter")) |
        (isInteractionType(loop, "promoter","ctcf")) | (isInteractionType(loop, "ctcf","ctcf"))
  
  element_loop <- loop[con]
  element_loop$dis <- calculateDistances(element_loop)
  return(element_loop)
}

get_sig_loop <- function(bg, element_loop, resolution=10000){
  require(fitdistrplus)
  #calculate pvalue for each distance
  dis <- (2:(10^6/resolution))*resolution
  
  for(d in dis){
    print(d)
    #get background for distance
    bg_d <- bg[bg$dis == d,]
    #calculate parameters for the distance
    param <- mledist(bg_d$counts, distr = "nbinom")
    
    #ks.test
    #ks.test(bg_d$counts, "pnbinom", size=param$estimate[1], 
    #        mu=param$estimate[2])
    
    #calculate significance
    if(d == 10^6){
      element_loop_d <- element_loop[element_loop$dis >= d,]
    }else{
      element_loop_d <- element_loop[element_loop$dis == d,]
    }
    
    element_loop_d$p.value <- pnbinom(element_loop_d$counts, 
                             size=param$estimate[1], 
                             mu=param$estimate[2],
                             lower.tail = F)
    
    #merge loop data for all distance
    if(!exists("all_element_loop")){
      all_element_loop <- element_loop_d
    }else{
      all_element_loop <- c(all_element_loop, element_loop_d)
    }
  }
  
  all_element_loop$fdr <- p.adjust(all_element_loop$p.value, method = "BH")
  return(all_element_loop)
}


if(T){
  resolution=10000
  #1. prepare genomic partion
  bin <- peak2gr("mm10_10000.bed")
  
  #2. read ctcf region
  ctcf <- peak2gr("ctcf.bed")
  
  ctcf.bin <- subsetByOverlaps(bin, ctcf)
  
  #jsut shrink region to the center 1 bp
  ctcf.bin.c <- GRanges(seqnames(ctcf.bin), 
                        IRanges(start(ctcf.bin) + resolution/2, 
                                start(ctcf.bin) + resolution/2 + 1))
  
  ctcf.ny.bin <- unique(c(shift(ctcf.bin, shift = -resolution), 
                          shift(ctcf.bin, shift = resolution)))
  
  ctcf.ny.bin.c <- GRanges(seqnames(ctcf.ny.bin), 
                           IRanges(start(ctcf.ny.bin) + resolution/2, 
                                   start(ctcf.ny.bin) + resolution/2 +1))
  
  rm(ctcf, ctcf.bin, ctcf.ny.bin)
  
  
  #3. read promoter region
  pro <- read.table("mm10_rna_info.tsv", 
                    header=T,sep="\t", stringsAsFactors = F)
  
  pro1 <- GRanges(paste0("chr",pro$chromosome_name),
                  IRanges(pro$transcription_start_site, 
                          pro$transcription_start_site+1))
  
  pro.bin <- subsetByOverlaps(bin, pro1)
  
  pro.bin.c <- GRanges(seqnames(pro.bin), 
                       IRanges(start(pro.bin) + resolution/2, 
                               start(pro.bin) + resolution/2 + 1))
  
  pro.ny.bin <- unique(c(shift(pro.bin, shift = -resolution), 
                         shift(pro.bin, shift = resolution)))
  
  pro.ny.bin.c <- GRanges(seqnames(pro.ny.bin), 
                          IRanges(start(pro.ny.bin) + resolution/2, 
                                  start(pro.ny.bin) + resolution/2 + 1))
  
  rm(pro, pro1, pro.ny.bin)
  
  
  #4. read enhancer region
  enh <- peak2gr("enhancer.bed")
  
  enh.bin <- subsetByOverlaps(bin, enh)
  
  enh.bin.c <- GRanges(seqnames(enh.bin), 
                       IRanges(start(enh.bin) + resolution/2, 
                               start(enh.bin) + resolution/2 + 1))
  
  #enhancer nearby regions
  enh.ny.bin <- unique(c(shift(enh.bin, shift = -resolution), 
                         shift(enh.bin, shift = resolution)))
  
  enh.ny.bin.c <- GRanges(seqnames(enh.ny.bin), 
                          IRanges(start(enh.ny.bin) + resolution/2, 
                                  start(enh.ny.bin) + resolution/2 + 1))
}

#5. load loop files treated by HiCcompare 
total_loop <- readr::read_delim("ctrl_shCTCF.10k.filter.tsv", 
                                col_names = F, delim = "\t")

#calculate significant loops for ctrl
if(T){
  ctrl <- total_loop %>% as_tibble %>% 
    dplyr::select(-'X8', -'X9') %>% 
    dplyr::mutate(across(7, round)) %>%
    dplyr::filter(X7>1)
  
  anchor1 <- ctrl %>% dplyr::select(1,2,3) %>% 
    rename("chr"=X1, "start"=X2, "end"=X3) %>% 
    GRanges()
  
  anchor2 <- ctrl %>% dplyr::select(4,5,6) %>% 
    rename("chr"=X4, "start"=X5, "end"=X6) %>% 
    GRanges()

  counts <- ctrl %>% pull(7)
  
  ctrl_loop <- GenomicInteractions(anchor1 = anchor1, 
                                   anchor2 = anchor2,
                                   counts = counts)
  rm(ctrl, anchor1, anchor2, counts)
  
  #6. loop annoation
  annotation.features = list(ctcf=ctcf.bin.c,
                             promoter=pro.bin.c,
                             enhancer=enh.bin.c,
                             ctcf_ny=ctcf.ny.bin.c,
                             promoter_ny=pro.ny.bin.c,
                             enhancer_ny=enh.ny.bin.c)
  
  suppressWarnings(annotateInteractions(ctrl_loop, 
                                        annotation.features))
  
  #annotationFeatures(loop)
  #categoriseInteractions(loop_dat)
  #plotInteractionAnnotations(loop)
  
  ctrl_bg <- get_loop_bg(ctrl_loop)
  ctrl_element_loop <- get_element_loop(ctrl_loop)
  ctrl_sig_loop <- get_sig_loop(ctrl_bg, ctrl_element_loop)
  
  rm(ctrl_loop, ctrl_bg, ctrl_element_loop)
}

#calculate significant loops for treatment
if(T){
  treat <- total_loop %>% 
    dplyr::select(-X7, -X9) %>% 
    dplyr::mutate(across(7, round)) %>% 
    dplyr::filter(X8>1)
  
  anchor1 <- treat %>% dplyr::select(1,2,3) %>% 
    rename("chr"=X1, "start"=X2, "end"=X3) %>% 
    GRanges()
  
  anchor2 <- treat %>% dplyr::select(4,5,6) %>% 
    rename("chr"=X4, "start"=X5, "end"=X6) %>% 
    GRanges()
  
  counts <- treat %>% pull(7)
  
  treat_loop <- GenomicInteractions(anchor1 = anchor1, 
                                   anchor2 = anchor2,
                                   counts = counts)
  
  rm(treat, anchor1, anchor2, counts)

  #6. loop annoation
  annotation.features = list(ctcf=ctcf.bin.c,
                             promoter=pro.bin.c,
                             enhancer=enh.bin.c,
                             ctcf_ny=ctcf.ny.bin.c,
                             promoter_ny=pro.ny.bin.c,
                             enhancer_ny=enh.ny.bin.c)
  
  suppressWarnings(annotateInteractions(treat_loop, 
                                        annotation.features))
  
  #annotationFeatures(loop)
  #categoriseInteractions(loop_dat)
  #plotInteractionAnnotations(loop)
  
  treat_bg <- get_loop_bg(treat_loop)
  treat_element_loop <- get_element_loop(treat_loop)
  treat_sig_loop <- get_sig_loop(treat_bg, treat_element_loop)
  
  rm(treat_loop, treat_bg, treat_element_loop)
}

#keep ctrl and treat loop information
for(q in c(0.05, 0.01)){
  print(q)
  #for ctrl
  print("For ctrl")
  name1 = "ctrl"
  ctrl_loop_kp <- ctrl_sig_loop %>% 
    as.data.frame %>% 
    dplyr::filter(fdr < {q}) %>%
    dplyr::select(c("seqnames1","start1","end1",
                    "seqnames2","start2","end2",
                    "counts","p.value","fdr"))
  write.table(ctrl_loop_kp, file=paste0(name1,".fdr_",q,".tsv"),
              sep="\t", quote=F, row.names = F)
  
  ctrl_loop_kp1 <- ctrl_loop_kp %>% 
    mutate(info = paste0(seqnames2,":",
                         start2,"-",
                         end2, ",",
                         counts)) %>%
    dplyr::select(1:3, info)
  
  write.table(ctrl_loop_kp1, file=paste0(name1,".washu.fdr_",q,".tsv"),
              col.names = F, row.names = F, sep="\t", quote=F)
  
  #for treatment
  print("For treatment")
  name2 = "shCTCF"
  treat_loop_kp <- treat_sig_loop %>% 
    as.data.frame %>% 
    dplyr::filter(fdr < {q}) %>%
    dplyr::select(c("seqnames1","start1","end1",
                    "seqnames2","start2","end2",
                    "counts","p.value","fdr"))
  write.table(treat_loop_kp, file=paste0(name2,".fdr_",q,".tsv"),
              sep="\t", quote=F, row.names = F)
  
  #for ctrl and treatment
  print("For ctrl and treatment")
  name3 = "ctrl_shCTCF"
  if(!exists("idx")){
    idx <- with(total_loop, 
                paste0(X1, X2, X3,
                       X4, X5, X6))
  }
  
  ctrl_idx <- with(ctrl_loop_kp, 
                   paste0(seqnames1, start1, end1,
                          seqnames2, start2, end2))
  
  treat_idx <- with(treat_loop_kp, 
                     paste0(seqnames1, start1, end1,
                            seqnames2, start2, end2))
  
  ctrl_treat_idx <- unique(c(ctrl_idx, treat_idx))
  
  ctrl_treat_loop <- total_loop[idx %in% ctrl_treat_idx, ]
  
  write.table(ctrl_treat_loop, file=paste0(name3, ".fdr_",q,".tsv"),
              sep="\t", quote=F, row.names = F)
}

  



