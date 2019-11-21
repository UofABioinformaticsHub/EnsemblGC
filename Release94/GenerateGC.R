library(tidyverse)
library(Biostrings)
library(RCurl)
library(here)

## Load the main function
source(file.path(here::here(), "getGcLen.R"))

## Set the Ensembl release
rls <- 94

## Set the working directory in a portable manner
setwd(here::here(paste0("Release", rls)))

## Define the species
sp <- c(
  Bs = "Bos_taurus",
  Ce = "Caenorhabditis_elegans", 
  Dr = "Danio_rerio", 
  Dm = "Drosophila_melanogaster",
  Gg = "Gallus_gallus", 
  Hs = "Homo_sapiens", 
  Mm = "Mus_musculus", 
  Oa = "Ovis_aries", 
  Rn = "Rattus_norvegicus", 
  Sc = "Saccharomyces_cerevisiae"
)

bld <- c(
  Bs = "UMD3.1",
  Ce = "WBcel235", 
  Dr = "GRCz11", 
  Dm = "BDGP6",
  Gg = "Gallus_gallus-5.0", 
  Hs = "GRCh38", 
  Mm = "GRCm38", 
  Oa = "Oar_v3.1", 
  Rn = "Rnor_6.0", 
  Sc = "R64-1-1"
)

db <- c(
  Bs = "ensembl",
  Ce = "ensembl", 
  Dr = "ensembl", 
  Dm = "ensembl",
  Gg = "ensembl", 
  Hs = "ensembl", 
  Mm = "ensembl", 
  Oa = "ensembl", 
  Rn = "ensembl", 
  Sc = "ensembl"
)

## Form a list of arguments
spList <- names(sp) %>%
  lapply(
    function(x){
      list(sp = sp[[x]], rls = rls, bld = bld[[x]], db = db[[x]])
    }
  )

## Run everything
lapply(spList, function(x){do.call(getGCLen, x)})
