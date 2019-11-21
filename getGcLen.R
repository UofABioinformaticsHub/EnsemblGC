#' @title Get the GC content and Length for a transcriptome
#'
#' @description Get the GC content and Length for an Ensembl transcriptome
#' 
#' @param sp The requested species. Must be specified in standard binomial 
#' nomenclature, with the genus capitalised and species in lower case. 
#' An underscore should be used instead of a space between names. For example,
#' sp = "Homo_sapiens".
#' @param rls The Ensembl release as an integer or character, e.g. rls = "98" 
#' @param bld The genome build associated with the corresponding Ensembl 
#' release, e.g. bld = "GRCh38"
#' @param db The Ensembl database for the requested species.
#' @param dir The location to download the file to. Defaults to tempdir()
#' 
#' @details 
#' This will download the requested transcriptome (i.e. cDNA.fa) from Ensembl
#' to a temporary location and calculate both the GC content and length of
#' every transcript in the downloaded file.
#' 
#' @return  
#' A tibble with GC content and length information per transcript.
#' 
#' @author Lachlan Baer & Steve Pederson
#' 
#' @examples 
#' library(tidyverse)
#' library(RCurl)
#' library(Biostrings)
#' gc <- getGCLen("Notechis_scutatus", 98, "TSX10Xv2-PRI")
#' 
#' @import stringr
#' @import dplyr
#' @import tidyr
#' @import magrittr
#' @import tibble
#' @import RCurl
#' @import Biostrings
getGCLen <- function(
  sp, 
  rls, 
  bld, 
  db = c("ensembl", "plants"), 
  seq_type = c("cdna", "cds"), 
  dir = tempdir(), 
  write_rds = TRUE
) {
  
  ## As this function exists outside of a package, the above imports are not 
  # going to be performed. Put them here and remove later if moving to a pkg
  reqPkg <- c(
    "stringr", "dplyr", "tidyr", "magrittr", "tibble", "RCurl", "Biostrings"
  )
  notLoaded <- setdiff(reqPkg, (.packages()))
  if (length(notLoaded)) {
    instPkg <- rownames(installed.packages())
    notInstalled <- setdiff(notLoaded, instPkg)
    if (length(notInstalled)) {
      if (!"BiocManager" %in% instPkg) install.packages("BiocManager")
      BiocManager::install(notInstalled)
    }
    sapply(notLoaded, library, character.only = TRUE)
  }
  
  ## Check all arguments
  if (!all(.checkSpecies(sp))) 
    stop("Species must be specified as 'Genus_species'")
  if (is.na(as.integer(rls)))
    stop("Ensembl release not specified correctly")
  seq_type <- match.arg(seq_type)
  bld <- match.arg(bld)
  db <- match.arg(db)
  stopifnot(dir.exists(dir))
  
  ## Define the file
  fl <- paste(sp, bld, seq_type, "all.fa.gz", sep = ".")
  ## Define the url. This is structured differently depending on database.
  if (db == "ensembl") {
    rootUrl <- paste0(
      "ftp://ftp.ensembl.org/pub/release-", 
      rls, 
      "/", 
      "fasta"
    )
  }
  if (db == "plants") {
    rootUrl <- paste0(
      "ftp://ftp.ensemblgenomes.org/pub/release-", 
      rls, 
      "/",
      "plants/",
      "fasta"
    )
  }
  ## Now form the complete url & download
  url <- paste(rootUrl, str_to_lower(sp), seq_type, fl, sep = "/")
  stopifnot(url.exists(url))
  localFa <- file.path(dir, fl)
  tryCatch(download.file(url, localFa))
  
  tbl <- .faToGC(localFa, bld)
  if (write_rds) saveRDS(tbl, paste0(sp, ".", bld, ".", rls, ".rds"))
  ## Silently return the object
  tbl
  
}

#' @description For use on local fa files
#' @param faFile Local Fasta File for processing
#' @param bld The genome build associated with the correspondig Ensembl 
#' release, e.g. bld = "GRCh38"
.faToGC <- function(faFile, bld){
  
  ## Import as a DNAStringSet
  tr <- readDNAStringSet(faFile)
  ## Split the sequence headers into an 7 column matrix using the 
  ## space delimiters. Beyond the first 6 columns, the sequences headers are 
  ## too variable both within and between species.
  ids <- str_split_fixed(names(tr), " ", 7)
  colnames(ids) <- c(
    "seq_id", "seq_type", "location", "gene_id", "gene_biotype", 
    "transcript_biotype", "misc"
  )
  df <- as_tibble(ids)
  ## Now create a data.frame (tibble) with all of the relevant information
  df <- dplyr::mutate(
    df,
    length = width(tr),
    gc = rowSums(alphabetFrequency(tr, baseOnly = TRUE)[,c("G", "C")]) / length,
    location = str_remove(location, paste0(".+", bld, ":")),
    gene_id = str_remove(gene_id, "gene:"),
    gene_biotype = str_remove(gene_biotype, "gene_biotype:"),
    transcript_biotype = str_remove(transcript_biotype, "transcript_biotype:"),
    ## Add the gene symbol, making sure there are NAs where it's not included
    gene_name =  str_extract(misc, "gene_symbol:[^ ]+"),
    gene_name = str_replace(gene_name, "gene_symbol:(.+)", "\\1")
  )
  df <- dplyr::select(df, -misc)
  df <- dplyr::arrange(df, seq_id) 
  ## Separate the genomic location into separate columns
  df <- tidyr::separate(
    data = df, 
    col = location, 
    into = c("seqnames", "start", "end", "strand"),
    sep = ":"
  ) 
  df <- dplyr::mutate(df, strand = ifelse(strand == 1, "+", "-"))
  df <- tidyr::unite(
    data = df,
    col = location,
    c("seqnames", "start", "end", "strand"),
    sep = ":"
  )
  
  ## Most species have ENS identifiers, but some don't
  ## If they do, they will have version numbers after transcript & gene ids
  ## which can be removed. Others (such as C elegans) won't
  hasENS <- all(grepl("ENS", df$seq_id))
  if (hasENS) {
    df$tx_id <- str_remove(df$seq_id, "\\.[0-9]+$")
    df$gene_id <- str_remove(df$gene_id, "\\.[0-9]+$")
  }
  
  ## Setup as a tibble then export as an RDS
  tbl <- as_tibble(df)
  tbl
}

.checkSpecies <- function(x){
  
  ## Check the species is specified correctly
  c(
    chkLength = length(x) == 1,
    chkCaps = str_split(x, "")[[1]][[1]] %in% LETTERS,
    chkUnderscore = str_detect(x, "_") == 1,
    chkNoSpace = str_detect(x, " ") == 0
  )
}
