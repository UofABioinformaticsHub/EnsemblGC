library(tidyverse)
library(Biostrings)
library(RCurl)
library(here)

## Load the main function
source(file.path(here::here(), "getGcLen.R"))

## Set the Ensembl release
rls <- 40

## Set the working directory in a portable manner
setwd(here::here(paste0("PlantsRelease", rls)))

## Define the species
sp <- c(
  At = "Arabidopsis_thaliana",
  Ta = "Triticum_aestivum", 
  Hv = "Hordeum_vulgare"
)

bld <- c(
  At = "TAIR10",
  Ta = "IWGSC", 
  Hv = "Hv_IBSC_PGSB_v2"
)

db <- c(
  At = "plants",
  Ta = "plants", 
  Hv = "plants"
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
