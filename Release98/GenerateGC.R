library(tidyverse)
library(Biostrings)
library(AnnotationHub)
library(GenomicRanges)
library(RCurl)
library(here)

## Load the main function
source(file.path(here::here(), "getGcLen.R"))

## Set the Ensembl release
rls <- 98

## Set the working directory in a portable manner
setwd(here::here(paste0("Release", rls)))

## Define the Annotation Hub object to extract the EnsDb objects from
ah <- subset(AnnotationHub(), rdataclass == "EnsDb")

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

## Form a list of arguments
spList <- lapply(sp, function(x){list(sp = x, rls = rls, ah = ah)})

## Run everything
lapply(spList, function(x){do.call(getGCLen, x)})
