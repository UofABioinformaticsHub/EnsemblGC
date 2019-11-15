#' Generate the GC Content for:
#' - Bos taurus
#' - Danio rerio
#' - Gallus gallus
#' - Homo sapiens
#' - Mus musculus
#' - Rattus norvegicus
#' 
#' All transcriptomes use the primary genome build for Ensembl release 95
#' 

library(dplyr)
library(tidyverse)
library(Biostrings)
library(AnnotationHub)
library(here)

setwd(here::here("Release95"))

spList <- list(
  Bs = c(sp = "Bos_taurus", bld = "ARS-UCD1.2", rls = "95"),
  Dr = c(sp = "Danio_rerio", bld = "GRCz11", rls = "95"),
  Gg = c(sp = "Gallus_gallus", bld = "GRCg6a", rls = "95"),
  Hs = c(sp = "Homo_sapiens", bld = "GRCh38", rls = "95"),
  Mm = c(sp = "Mus_musculus", bld = "GRCm38", rls = "95"),
  Rn = c(sp = "Rattus_norvegicus", bld = "Rnor_6.0", rls = "95")
)

getGCLen <- function(x) {
  sp <- x[1]
  bld <- x[2]
  rls <- x[3]
  url <- paste0(
    "ftp://ftp.ensembl.org/pub/release-", 
    rls, 
    "/fasta/", 
    str_to_lower(sp), 
    "/cdna/", 
    sp,
    ".",
    bld, 
    ".cdna.all.fa.gz"
  )
  download.file(url, file.path(tempdir(), "temp.fa"))
  tr <- file.path(tempdir(), "temp.fa") %>% readDNAStringSet()
  df <- tibble(
    id = names(tr), 
    length = width(tr), 
    gc = rowSums(alphabetFrequency(tr, baseOnly = TRUE)[,c("G", "C")]) / length
  ) %>% 
    mutate(
      tx_id = str_extract(id, "ENS[A-Z]+[0-9]+"),
      gene_id = str_extract(id, "ENS[A-Z]*G[0-9]+"),
      location = str_extract(id, paste0("[a-z]+:", bld, ":.+:[0-9]+:[0-9]+:-?[0-9]+")),
      location = str_remove(location, paste0("[a-z]+:", bld, ":")),
      gene_symbol =  str_extract(id, "gene_symbol:.+ desc"),
      gene_symbol = str_remove(gene_symbol, "gene_symbol:"),
      gene_symbol = str_remove(gene_symbol, " desc")
    ) %>%
    dplyr::select(-id) %>%
    dplyr::arrange(tx_id, desc(length)) %>%
    separate("location", c("seqnames", "start", "end", "strand"), sep = ":") %>%
    mutate(strand = ifelse(strand == 1, "+", "-")) %>%
    as.data.frame() %>%
    column_to_rownames("tx_id")
  ah <- AnnotationHub()
  ahRls <- ah %>%
    subset(species == str_replace(sp, "_", " ")) %>%
    subset(dataprovider == "Ensembl") %>%
    subset(rdataclass == "EnsDb")
  ens <- ah[[
    ahRls$ah_id[
      ahRls$title == grep(
        paste0("Ensembl ", rls), ahRls$title, value = TRUE
      )
      ]
    ]]
  gr <- makeGRangesFromDataFrame(
    df, 
    seqinfo = seqinfo(ens),
    keep.extra.columns = TRUE
  )
  saveRDS(gr, paste0(sp, ".", bld, ".", rls, ".rds"))
}

lapply(spList, getGCLen)
