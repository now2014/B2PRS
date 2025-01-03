#!/usr/bin/env Rscript

library(BEDMatrix)
library(data.table)
library(parallel)

check.plink.files <- function(cur.chrom, bed.file.pattern='./chr@.bed'){
  if(cur.chrom==23 || cur.chrom=='23' || cur.chrom=='X'){
    bed.file.23 <- gsub('@', '23', bed.file.pattern)
    bed.file.X <- gsub('@', 'X', bed.file.pattern)
    cur.chrom <- 23
    if(file.exists(bed.file.23)){
      bed.file <- bed.file.23
    }else if(file.exists(bed.file.X)){
      bed.file <- bed.file.X
    }else{
      bed.file <- NA
    }
  }else{
    cur.chrom <- as.numeric(cur.chrom)
    bed.file <- gsub('@', cur.chrom, bed.file.pattern)
    if(!file.exists(bed.file)){
      bed.file <- NA
    }
  }

  ## check .bim file
  bim.file <- gsub('.bed$', '.bim', bed.file)
  fam.file <- gsub('.bed$', '.fam', bed.file)
  if(!file.exists(bim.file)) bim.file <- NA
  if(!file.exists(fam.file)) fam.file <- NA
  if(any(is.na(c(bed.file, bim.file, fam.file)))){
    return(NULL)
  }
  return(bed.file)
}

calc.PRS <- function(kept.snps, bed.file, fill.missing=c('mean', 'EAF')){
  ## bed.file.pattern: path to the bed file
  ## kept.snps: a data.frame with columns of 'rsid', 'EA', 'OA', 'EAF', and 'b1', 'b2', 'b3', ...
  ## fill.missing: a string of 'EAF' and/or 'mean' to fill missing genotype with 'EAF' column or mean frequency from genotype data
  fill.missing <- match.arg(fill.missing)
  
  ## remove duplicated rsids
  idx.dup <- duplicated(kept.snps$rsid)
  if(any(idx.dup)){
    message('Duplicated rsids in kept.snps are removed.')
    kept.snps <- kept.snps[!idx.dup, ]
  }
  row.names(kept.snps) <- kept.snps$rsid

  direc <- rep(c(1, -1), each=nrow(kept.snps))
  names(direc) <- c(
    paste0(kept.snps$rsid, ':', kept.snps$EA, ':', kept.snps$OA),
    paste0(kept.snps$rsid, ':', kept.snps$OA, ':', kept.snps$EA)
  )

  bim.file <- gsub('.bed$', '.bim', bed.file)
  bim <- fread(bim.file, drop=c(1, 3, 4), col.names=c('rsid', 'A1', 'A2'))
  bim$i <- seq_len(nrow(bim))
  bim <- bim[bim$rsid %in% kept.snps$rsid, ]
  ## remove duplicated rsids
  idx.dup <- duplicated(bim$rsid)
  if(any(idx.dup)){
    message('Duplicated rsids in .bim file are removed.')
    bim <- bim[!idx.dup, ]
  }
  bim$snp <- paste0(bim$rsid, ':', bim$A1, ':', bim$A2)
  bim <- bim[bim$snp %in% names(direc), ]
  bim$direc <- direc[bim$snp]

  ## remove duplicated rsids
  idx.dup <- duplicated(bim$rsid)
  if(any(idx.dup)){
    message('Duplicated rsids in .bim file are removed.')
    bim <- bim[!idx.dup, ]
  }

  i.snps <- bim$i
  if(nrow(bim)==0){
    warning('No SNPs for PRS are found in the bed file: ', bed.file)
    return(NULL)
  }
  ## read bed file
  bed <- suppressMessages(BEDMatrix(bed.file, simple_names = TRUE))
  X <- as.matrix(bed[, i.snps])  # N samples x M SNPs
  colnames(X) <- bim$rsid
  rownames(X) <- rownames(bed)
  idx.change.col <- bim$direc==-1
  X[, idx.change.col] <- 2 - X[, idx.change.col]

  ## fill missing genotypes
  if(fill.missing=='EAF'){
    fill.cnt <- kept.snps$EAF * 2
    names(fill.cnt) <- kept.snps$rsid
    fill.cnt <- fill.cnt[bim$rsid]
  }else if(fill.missing=='mean'){
    fill.cnt <- colMeans(X, na.rm=TRUE)
    fill.cnt[is.na(fill.cnt)] <- 0
  }else{
    stop('fill.missing should be either "EAF" or "mean".')
  }
  fill.mat <- t(t(is.na(X)) * fill.cnt)
  X[is.na(X)] <- 0
  X <- X + fill.mat

  b <- kept.snps[colnames(X), ]
  b <- b[, setdiff(colnames(b), c('rsid', 'EA', 'OA', 'EAF', 'POS38')), drop=FALSE] # M SNPs x K traits
  b <- as.matrix(b)
  b[is.na(b)] <- 0

  PRS <- X %*% b # N samples x K traits
  return(PRS)
}

match.PRS.samples <- function(PRS, samples){
  if(is.null(PRS)){
    PRS <- matrix(0, nrow=length(samples), ncol=1)
    rownames(PRS) <- samples
  }
  PRS <- PRS[intersect(samples, rownames(PRS)), 1, drop=FALSE]
  PRS.match <- matrix(NA, nrow=length(samples), ncol=1)
  rownames(PRS.match) <- samples
  PRS.match[rownames(PRS), ] <- PRS
  PRS.match[is.na(PRS.match)] <- mean(PRS.match, na.rm=TRUE)
  return(PRS.match)
}

PRS.chr <- function(cur.chrom, kept.snps, kept.samples=NULL, bed.file.pattern='./chr@.bed', fill.missing=c('mean', 'EAF'), chunk.size=1000, mc.cores=1){
  kept.snps <- kept.snps[kept.snps$CHROM==cur.chrom, -c(1, 2)]

  bed.file <- check.plink.files(cur.chrom, bed.file.pattern)
  if(is.null(bed.file)){
    warning('PLINK files (.bed + .bim + .fam) are missing for chr', cur.chrom, '. Skip this chromosome.')
    return(NULL)
  }
  if(is.null(kept.samples)){
    kept.samples <- fread(gsub('bed$', 'fam', bed.file), select=1, header=FALSE)
  }
  kept.samples <- as.character(unlist(kept.samples))
  bim.rsids <- fread(gsub('bed$', 'bim', bed.file), select=2, header=FALSE)
  bim.rsids <- as.character(unlist(bim.rsids))

  kept.snps <- kept.snps[kept.snps$rsid %in% bim.rsids, ]
  M <- nrow(kept.snps)
  if(M==0){
    warning('No SNPs for PRS are found on chromosome ', cur.chrom)
    return(NULL)
  }
  chunks <- seq_len(ceiling(M/chunk.size))

  PRS <- parallel::mclapply(
    chunks,
    function(i){
      idx.start <- (i-1) * chunk.size + 1
      idx.end <- min(i * chunk.size, M)
      kept.snps.cur <- kept.snps[idx.start:idx.end, ]
      message('Calculating PRS for ', idx.start, '-', idx.end, ' SNPs on chromosome ', cur.chrom)
      PRS.cur <- calc.PRS(kept.snps.cur, bed.file, fill.missing)
      match.PRS.samples(PRS.cur, kept.samples)
    },
    mc.cores=mc.cores
  )
  if(all(sapply(PRS, is.null))){
    return(NULL)
  }
  PRS <- Reduce('+', PRS)
  PRS <- matrix(PRS, ncol=1)
  rownames(PRS) <- kept.samples
  return(PRS)
}

PRS.wg <- function(
  EAF.b.rds='./PRS-b/EAF.b.rds',
  snp.info.rds='./PRS-b/snp.info.rds',
  bed.file.pattern='./chr@.bed',
  fill.missing=c('mean', 'EAF'),
  chunk.size=1000, mc.cores=1, PRS.name='PRS'){
  #### PRS for whole genome ####
  ## EAF.b.rds: path to the RDS file with EAF and b for PRS calculation
  ## snp.info.rds: path to the RDS file with SNP information
  ## bed.file.pattern: a pattern of bed file names, e.g., './chr@.bed', where '@' is replaced by chromosome number
  ## fill.missing: a string of 'EAF' and/or 'mean' to fill missing genotype with 'EAF' column or mean frequency from genotype data
  ## chunk.size: max number of SNPs to be read from bed file in each chunk
  ## mc.cores: number of cores for parallel computing
  ## PRS.name: name of the PRS column in the output

  kept.snps <- cbind(readRDS(snp.info.rds), readRDS(EAF.b.rds))
  chroms <- unique(kept.snps$CHROM)

  bed.files <- lapply(chroms, check.plink.files, bed.file.pattern)
  idx.null <- sapply(bed.files, is.null)
  if(any(idx.null)){
    missing.chroms <- paste0(chroms[idx.null], collapse=', ')
    warning('Skip chromosomes [', missing.chroms, '] due to missing PLINK files.')
    chroms <- chroms[!idx.null]
    bed.files <- bed.files[!idx.null]
  }

  kept.samples <- fread(gsub('bed$', 'fam', bed.files[1]), select=1, header=FALSE)
  kept.samples <- as.character(unlist(kept.samples))

  PRS <- NULL
  for(cur.chrom in chroms){
    PRS.cur <- PRS.chr(cur.chrom, kept.snps, kept.samples, bed.file.pattern, fill.missing, chunk.size, mc.cores)
    if(is.null(PRS.cur)){
      next
    }
    PRS.cur <- match.PRS.samples(PRS.cur, kept.samples)
    if(is.null(PRS)){
      PRS <- PRS.cur
    }else{
      PRS <- PRS + PRS.cur
    }
  }
  if(is.null(PRS)){
    warning('No SNP is found for PRS calculation.')
    return(NULL)
  }
  colnames(PRS) <- PRS.name
  return(PRS)
}
