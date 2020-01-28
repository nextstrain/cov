# Code and Data for Human Coronavirus Builds

## About this repository

This respository currently includes the data (as of 23 Jan 20) and code to run Nextstrain `augur` builds for Coronaviruses 229E, NL63, SARS, and Betacoronavirus1. These builds focus on the human transmission of these viruses, so do not necessarily include all available samples, focusing instead on the part of the tree which includes the majority of human infections.

This respository uses a basic form of 'Databasing' to make adding new sequences from ViPR easier.
This means you can always download the full latest metadata file from ViPR, yet only sequences that have not already been seen will be downloaded and aligned again, and then simply added to the sequences that already exist.

I recommend running the `database` rule(s) first to ensure this works ok, then proceed to running the `augur`/`build` rules.

# How to Run

## Initial Database Setup
To help make initial runs faster, all the data on ViPR as of 23-Jan-20 is included in this respository. 
To avoid having the whole thing start from scratch (including downloading and aligning ~600 sequences), after cloning this repo, rename the `genbank-for-database` folder to `genbank`.
(Keeping this separate means that others pushing won't affect your own database...)

Then run:
```
snakemake -t make_database_all
```
This will (hopefully) set all the files so that Snakemake won't try to re-build the databases from scratch. 
Do this **before** trying to add any new sequences.

## Adding New Sequences
Adding new sequences is easy.
Navigate to ViPR and use the search function, ensuring that 'Full Genome' is ticked at the bottom of the search.
Navigate first to 'Coronaviridae' and 'Coronavirinae', then "Select all" in the following categories:
* In 'Betacoronavirus':
    * For **Betacoronavirus 1** select 'Betacoronavirus 1' (about 250)
    * For **SARS** select 'Severe acute respiratory syndrome-related coronavirus' (about 330) and 'SARS-related coronavirus' (about 9)
* In 'Alphacoronavirus':
    * For **NL63** select 'Human coronavirus NL63' (about 70)
    * For **229E** select 'Human coronavirus 229E' (about 40)

Once each search results load, select all sequences, and download the samples in *tab delimited format*.
Once downloaded, move them into the `data` folder, and name them according to the specified format, setting the correct date (this makes debug easier).
Go into the `Snakefile` and at the top of the `files` rule, change the input file to match the one you just created. 

_(To be extra clear - you need to do each of these downloads separately - one for SARS, one for NL63... etc!)_

This should trigger Snakemake to build a new database for whatever files you've added.
You can build all the databases with `snakemake make_database_all` or build them individually with the `snakemake make_database_229E` (and similar) rules (see Snakefile for all rules).

## Nextstrain/Augur Build
Once the database is ready, you can run the Nexstrain or `augur` build.
You can run this for all the viruses using just `snakemake`, or for individual ones with `snakemake build_229E` (and similar - see Snakefile for all).
Different viruses have different rules for what is filtered and what paramters are used.
These are usually displayed while running, but see the Snakefile and/or below for more detail.

For the 229E and Betacoronavirus1 builds, you can also run analyses that split the genome into the replicase gene and all remaining genes, and run a separate analysis for each, using the commands `snakemake build_229E_tangle` and `snakemake build_beta_tangle`.
This allows you to create a tanglegram using the two sections of the genome.

All output is stored in the `cov/auspice` folder, to make it easier to copy/view them all together.

### Adding extra data
Sometimes, data is not available from Genbank or ViPR automatically, but is available by manually reading the Genbank file or from reading associated publications.
This -- in particular `host` and `date` information -- can be very useful to add manually.
This pipeline allows this - add extra metadata to the (virus)/`data/extra_meta_`(virus)`.tsv` file, following the format of the files that already have data in them.
Ensure the values are separated by tabs, not spaces!

# Notes on Each Build

For each build, see a list of the manually excluded/included sequences in the respective (virus)`/config/exclude_`(virus)`.txt` and `include_`(virus)`.txt`.

## 229E
### Filtering
Filters only to Human samples
### Mutation rate
Estimates mutation rate - currently coming out to ~2.6 x 10^-4 subs/site/yr.

## SARS
### Filtering
Excludes all Bat and Mouse samples. Some Mouse were with humans but all are 'mouse-adapted strains' and thus not appropriate to include in a human-focused run.
### Mutation rate
Estimates mutation rate - currently coming out to ~4.0 x 10^-4 subs/site/yr.
### Other notable exclusions
* MK062179 through MK062184 are modified viruses for vaccine trials
* GU553365 was a human isolate from 2003, passaged through money lungs and then Vero E6 cells
* 'VeroE6' in strain name are passaged and divergent
* JX163923 to JX163928 were isolated from vaccine-challenged ferrets
* Those with with "wtic" and "ExonN1 mutant" in the name are cultured from engineered viruses

## Betacoronavirus1
### Filtering
Filtered to only human and chimp samples.
### Mutation rate
Estimates mutation rate - currently coming out to ~2.8 x 10^-4 subs/site/yr.
### Other notable exclusions
* Strains ATCC VR-759 (accessions AY391777 and AY585228) are the 'prototype strain' (ATCC) sequences. These were originally isolated in 1967 England, but were propogated on cell culture.
* AY585229 is an early (2004) sequence of Betacoronavirus 1 isolated in 2001, but was sequenced same time as one of the prototypes above. It has little divergence from the prototype strains, despite the fact that it should be ~40 years diverged!

## NL63
### Filtering
Not filtered by any host.
However, sequences are masked to exclude bases 20,500 to 22,300, where there seems to be a lot of recombination, to help get a better tree.
### Mutation rate
Uses set mutation rate of 4.59 x 10^-4 as estimated for MERS by Dudas et al. 2018. eLife.
### Other notable exclusions
AY518894 and KF530106 seem to be recombinants in the Spike protein region.
### Other notes
This tree is particularly difficult. There seems to be a lot of recombination.
Unsure how to improve the tree further.

KT381875, KU521535, and KX179500 seem to be the same sample? They are not currently excluded but perhaps should be? ([Reference](https://www.ncbi.nlm.nih.gov/pubmed/27799635))

AY518894 is an isolate from 1988 but was passaged in cells before sequencing.

## Notes on NL63 and Betacoronavirus1 -- Seattle sequences
Both runs have sequences from [this study](https://academic.oup.com/jid/article/216/2/203/3858443) and were sequenced in Seattle. 
There are multiple samples from the same patients, and these cluster in both runs.
NL63 has N07 and N06 samples, each patient's samples cluster.
Betacoronavirus1 has N07, N08, and N09 samples (though N08-33B is misnamed and should be N09-33B), and patient's samples cluster.
However, in **BOTH RUNS** these samples are *under-diverged* - they have fewer mutations than expected. 
This may indicate something is wrong with the sequencing or something else.
Currently they are not removed, but possibly should be?
*Update: These sequences are now removed from the Betacoronavirus1 run*