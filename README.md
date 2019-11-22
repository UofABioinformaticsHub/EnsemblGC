# Ensembl_GC

**This site is best viewed and used using the github pages site:**
https://uofabioinformaticshub.github.io/Ensembl_GC/

This repository holds precalculated GC content for select species from multiple ensembl releases.
All RDS files contain the GC content and length for the transcripts as defined at time of release, as tibbles. Links are provided below.

To load these directly into your R session you can execute something similar to the following lines of code.
This example will import gc content information for the *Homo sapiens* transcriptome based on Ensembl release 98.
Simply copy the url for your desired transcriptome in place of the human one, and this should work. 

```
library(GenomicRanges)
con <- url("https://uofabioinformaticshub.github.io/Ensembl_GC/Release98/Homo_sapiens.GRCh38.98.rds")
gc <- readRDS(con)
```

Alternatively if you prefer a more 'tidy' approach.

```
library(GenomicRanges)
library(magrittr)
gc <- url("https://uofabioinformaticshub.github.io/Ensembl_GC/Release98/Homo_sapiens.GRCh38.98.rds") %>% 
    readRDS()
```

The script used to generate all objects can be found in the respective folders of each release on the main repository https://github.com/UofABioinformaticsHub/Ensembl_GC

## Ensembl 98
Last generated 21 November, 2019.
- [Bos taurus](Release98/Bos_taurus.ARS-UCD1.2.98.rds)
- [Caenorhabditis elegans](Release98/Caenorhabditis_elegans.WBcel235.98.rds)
- [Danio rerio](Release98/Danio_rerio.GRCz11.98.rds)
- [Drosophila melanogaster](Release98/Drosophila_melanogaster.BDGP6.22.98.rds)
- [Gallus gallus](Release98/Gallus_gallus.GRCg6a.98.rds)
- [Homo sapiens](Release98/Homo_sapiens.GRCh38.98.rds)
- [Mus musculus](Release98/Mus_musculus.GRCm38.98.rds)
- [Ovis aries](Release98/Ovis_aries.Oar_v3.1.98.rds)
- [Rattus norvegicus](Release98/Rattus_norvegicus.Rnor_6.0.98.rds)
- [Saccharomyces cerevisiae](Release98/Saccharomyces_cerevisiae.R64-1-1.98.rds)

## EnsemblPlants 45
Last generated 22 November, 2019.
- [Arabidopsis thaliana](PlantsRealease45/Arabidopsis_thaliana.TAIR10.45.rds)
- [Hordeum vulgare](PlantsRealease45/Hordeum_vulgare.IBSC_v2.45.rds)
- [Triticum aestivum](PlantsRealease45/Triticum_aestivum.IWGSC.45.rds)

## Ensembl 97
Last generated 21 November, 2019.
- [Bos taurus](Release97/Bos_taurus.ARS-UCD1.2.97.rds)
- [Caenorhabditis elegans](Release97/Caenorhabditis_elegans.WBcel235.97.rds)
- [Danio rerio](Release97/Danio_rerio.GRCz11.97.rds)
- [Drosophila melanogaster](Release97/Drosophila_melanogaster.BDGP6.22.97.rds)
- [Gallus gallus](Release97/Gallus_gallus.GRCg6a.97.rds)
- [Homo sapiens](Release97/Homo_sapiens.GRCh38.97.rds)
- [Mus musculus](Release97/Mus_musculus.GRCm38.97.rds)
- [Ovis aries](Release97/Ovis_aries.Oar_v3.1.97.rds)
- [Rattus norvegicus](Release97/Rattus_norvegicus.Rnor_6.0.97.rds)
- [Saccharomyces cerevisiae](Release97/Saccharomyces_cerevisiae.R64-1-1.97.rds)

## EnsemblPlants 44
Last generated 22 November, 2019.
- [Arabidopsis thaliana](PlantsRealease44/Arabidopsis_thaliana.TAIR10.44.rds)
- [Hordeum vulgare](PlantsRealease44/Hordeum_vulgare.IBSC_v2.44.rds)
- [Triticum aestivum](PlantsRealease44/Triticum_aestivum.IWGSC.44.rds)

## Ensembl 96
Last generated 21 November, 2019.
- [Bos taurus](Release96/Bos_taurus.ARS-UCD1.2.96.rds)
- [Caenorhabditis elegans](Release96/Caenorhabditis_elegans.WBcel235.96.rds)
- [Danio rerio](Release96/Danio_rerio.GRCz11.96.rds)
- [Drosophila melanogaster](Release96/Drosophila_melanogaster.BDGP6.22.96.rds)
- [Gallus gallus](Release96/Gallus_gallus.GRCg6a.96.rds)
- [Homo sapiens](Release96/Homo_sapiens.GRCh38.96.rds)
- [Mus musculus](Release96/Mus_musculus.GRCm38.96.rds)
- [Ovis aries](Release96/Ovis_aries.Oar_v3.1.96.rds)
- [Rattus norvegicus](Release96/Rattus_norvegicus.Rnor_6.0.96.rds)
- [Saccharomyces cerevisiae](Release96/Saccharomyces_cerevisiae.R64-1-1.96.rds)

## EnsemblPlants 43
Last generated 22 November, 2019.
- [Arabidopsis thaliana](PlantsRealease43/Arabidopsis_thaliana.TAIR10.43.rds)
- [Hordeum vulgare](PlantsRealease43/Hordeum_vulgare.IBSC_v2.43.rds)
- [Triticum aestivum](PlantsRealease43/Triticum_aestivum.IWGSC.43.rds)

## Ensembl 95
Last generated 21 November, 2019.
- [Bos taurus](Release95/Bos_taurus.ARS-UCD1.2.95.rds)
- [Caenorhabditis elegans](Release95/Caenorhabditis_elegans.WBcel235.95.rds)
- [Danio rerio](Release95/Danio_rerio.GRCz11.95.rds)
- [Drosophila melanogaster](Release95/Drosophila_melanogaster.BDGP6.22.95.rds)
- [Gallus gallus](Release95/Gallus_gallus.GRCg6a.95.rds)
- [Homo sapiens](Release95/Homo_sapiens.GRCh38.95.rds)
- [Mus musculus](Release95/Mus_musculus.GRCm38.95.rds)
- [Ovis aries](Release95/Ovis_aries.Oar_v3.1.95.rds)
- [Rattus norvegicus](Release95/Rattus_norvegicus.Rnor_6.0.95.rds)
- [Saccharomyces cerevisiae](Release95/Saccharomyces_cerevisiae.R64-1-1.95.rds)

## EnsemblPlants 42
Last generated 22 November, 2019.
- [Arabidopsis thaliana](PlantsRealease42/Arabidopsis_thaliana.TAIR10.42.rds)
- [Hordeum vulgare](PlantsRealease42/Hordeum_vulgare.IBSC_v2.42.rds)
- [Triticum aestivum](PlantsRealease42/Triticum_aestivum.IWGSC.42.rds)

## Ensembl 94
Last generated 21 November, 2019.
- [Bos taurus](Release94/Bos_taurus.UMD3.1.94.rds)
- [Caenorhabditis elegans](Release94/Caenorhabditis_elegans.WBcel235.94.rds)
- [Danio rerio](Release94/Danio_rerio.GRCz11.94.rds)
- [Drosophila melanogaster](Release94/Drosophila_melanogaster.BDGP6.22.94.rds)
- [Gallus gallus](Release94/Gallus_gallus.Gallus_gallus-5.0.94.rds)
- [Homo sapiens](Release94/Homo_sapiens.GRCh38.94.rds)
- [Mus musculus](Release94/Mus_musculus.GRCm38.94.rds)
- [Ovis aries](Release94/Ovis_aries.Oar_v3.1.94.rds)
- [Rattus norvegicus](Release94/Rattus_norvegicus.Rnor_6.0.94.rds)
- [Saccharomyces cerevisiae](Release94/Saccharomyces_cerevisiae.R64-1-1.94.rds)

## EnsemblPlants 41
Last generated 22 November, 2019.
- [Arabidopsis thaliana](PlantsRealease41/Arabidopsis_thaliana.TAIR10.41.rds)
- [Hordeum vulgare](PlantsRealease41/Hordeum_vulgare.IBSC_v2.41.rds)
- [Triticum aestivum](PlantsRealease41/Triticum_aestivum.IWGSC.41.rds)

## Ensembl 93
Last generated 21 November, 2019.
- [Bos taurus](Release93/Bos_taurus.UMD3.1.93.rds)
- [Caenorhabditis elegans](Release93/Caenorhabditis_elegans.WBcel235.93.rds)
- [Danio rerio](Release93/Danio_rerio.GRCz11.93.rds)
- [Drosophila melanogaster](Release93/Drosophila_melanogaster.BDGP6.22.93.rds)
- [Gallus gallus](Release93/Gallus_gallus.Gallus_gallus-5.0.93.rds)
- [Homo sapiens](Release93/Homo_sapiens.GRCh38.93.rds)
- [Mus musculus](Release93/Mus_musculus.GRCm38.93.rds)
- [Ovis aries](Release93/Ovis_aries.Oar_v3.1.93.rds)
- [Rattus norvegicus](Release93/Rattus_norvegicus.Rnor_6.0.93.rds)
- [Saccharomyces cerevisiae](Release93/Saccharomyces_cerevisiae.R64-1-1.93.rds)

## EnsemblPlants 40
Last generated 22 November, 2019.
- [Arabidopsis thaliana](PlantsRealease40/Arabidopsis_thaliana.TAIR10.40.rds)
- [Hordeum vulgare](PlantsRealease40/Hordeum_vulgare.Hv_IBSC_PGSB_v2.40.rds)
- [Triticum aestivum](PlantsRealease40/Triticum_aestivum.IWGSC.40.rds)

## Ensembl 92
Last generated 21 November, 2019.
- [Bos taurus](Release92/Bos_taurus.UMD3.1.92.rds)
- [Caenorhabditis elegans](Release92/Caenorhabditis_elegans.WBcel235.92.rds)
- [Danio rerio](Release92/Danio_rerio.GRCz11.92.rds)
- [Drosophila melanogaster](Release92/Drosophila_melanogaster.BDGP6.22.92.rds)
- [Gallus gallus](Release92/Gallus_gallus.Gallus_gallus-5.0.92.rds)
- [Homo sapiens](Release92/Homo_sapiens.GRCh38.92.rds)
- [Mus musculus](Release92/Mus_musculus.GRCm38.92.rds)
- [Ovis aries](Release92/Ovis_aries.Oar_v3.1.92.rds)
- [Rattus norvegicus](Release92/Rattus_norvegicus.Rnor_6.0.92.rds)
- [Saccharomyces cerevisiae](Release92/Saccharomyces_cerevisiae.R64-1-1.92.rds)

## Ensembl 91
Last generated 21 November, 2019.
- [Bos taurus](Release91/Bos_taurus.UMD3.1.91.rds)
- [Caenorhabditis elegans](Release91/Caenorhabditis_elegans.WBcel235.91.rds)
- [Danio rerio](Release91/Danio_rerio.GRCz10.91.rds)
- [Drosophila melanogaster](Release91/Drosophila_melanogaster.BDGP6.22.91.rds)
- [Gallus gallus](Release91/Gallus_gallus.Gallus_gallus-5.0.91.rds)
- [Homo sapiens](Release91/Homo_sapiens.GRCh38.91.rds)
- [Mus musculus](Release91/Mus_musculus.GRCm38.91.rds)
- [Ovis aries](Release91/Ovis_aries.Oar_v3.1.91.rds)
- [Rattus norvegicus](Release91/Rattus_norvegicus.Rnor_6.0.91.rds)
- [Saccharomyces cerevisiae](Release91/Saccharomyces_cerevisiae.R64-1-1.91.rds)

## Ensembl 90
Last generated 21 November, 2019.
- [Bos taurus](Release90/Bos_taurus.UMD3.1.90.rds)
- [Caenorhabditis elegans](Release90/Caenorhabditis_elegans.WBcel235.90.rds)
- [Danio rerio](Release90/Danio_rerio.GRCz10.90.rds)
- [Drosophila melanogaster](Release90/Drosophila_melanogaster.BDGP6.22.90.rds)
- [Gallus gallus](Release90/Gallus_gallus.Gallus_gallus-5.0.90.rds)
- [Homo sapiens](Release90/Homo_sapiens.GRCh38.90.rds)
- [Mus musculus](Release90/Mus_musculus.GRCm38.90.rds)
- [Ovis aries](Release90/Ovis_aries.Oar_v3.1.90.rds)
- [Rattus norvegicus](Release90/Rattus_norvegicus.Rnor_6.0.90.rds)
- [Saccharomyces cerevisiae](Release90/Saccharomyces_cerevisiae.R64-1-1.90.rds)