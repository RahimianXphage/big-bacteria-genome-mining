# big-bacteria-genome-mining
Genomic survey of 82 giant bacteria (Beggiatoa, Thioploca, Thiomargarita, Epulonipiscium) revealing ARGs, BGCs, prophages, plasmids, and defensive-metabolic islands
## Overview

This repository contains data from the manuscript:

**"Genomic diversification of giant bacteria: Antimicrobial resistance, secondary metabolites, and phage-host coevolution in free-living and symbiotic lineages"**

Mohammadreza Rahimian

This study presents a comprehensive genomic survey of 82 assemblies from six genera of giant bacteria: *Beggiatoa*, *Thioploca*, *Thiomargarita*, *Candidatus* Epulonipiscium, *Candidatus* Epulonipiscioides, and *Candidatus* Parepulonipiscium.

## Repository Contents

### /plasmids/
- `plasmids.fasta` — Full nucleotide sequences of 10 identified plasmids
- `plasmids_annotation.xlsx` — Functional annotation of plasmid genes

### /prophages/
- `prophages.fasta` — Full nucleotide sequences of 11 quality-filtered prophages
- `prophages_annotation.xlsx` — Functional annotation of prophage genes

### /msa_phage_proteins/
Multiple sequence alignments (FASTA format) of core phage proteins used for HyPhy co-evolutionary analysis

## Methods Summary

- **ARG Prediction:** DeepARG v1.0.4, ABRicate v1.0.1 (AMRFinderPlus, CARD, Resfinder, ARG-ANNOT, MEGARES)
- **BGC Prediction:** antiSMASH v8.0.4 (strict mode)
- **Defense Systems:** DefenseFinder v2.0.0, CRISPR-Cas++ (CRISPRCasFinder, CRISPRCasMeta)
- **Prophage Prediction:** PHASTEST v1.0.1, VirSorter2 v2.2.4, PhaMer, geNomad v1.11.1
- **Plasmid Prediction:** geNomad v1.11.1
- **Quality Control:** CheckV v1.5 (>90% completeness, <5% contamination)
- **Phylogenetics:** ViPTree v4.0, IQ-TREE v2.4.0
- **Coevolutionary Analysis:** HyPhy v2.5.96 (GARD, BGM)

## Citation

If you use this data, please cite:

Rahimian M. (2026). Genomic diversification of giant bacteria: Antimicrobial resistance, secondary metabolites, and phage-host coevolution in free-living and symbiotic lineages. [Journal Name], Volume(Issue), Pages. DOI: [DOI]

## Contact

Mohammadreza Rahimian
Email: rahimianmohammadreza66@gmail.com

## Acknowledgments

NCBI Genome database provides genome assemblies. All bioinformatics tools used in this study are acknowledged in the manuscript.
