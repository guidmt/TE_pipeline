# TE_pipeline
 
# How to run:
1) Install bpipe
2) Configure the parameters of *.bpipe files 
3) Run:

```bash
bpipe  step1_te_pipeline.bpipe
bpipe  step2_te_pipeline.bpipe
```
## Run background:
```bash
tmux new -s step1_te_pipeline
bpipe  step1_te_pipeline.bpipe
ctrl+b+d

recover session:
tmux attach-session -t step1_te_pipeline
```

# Dependencies:
- bpipe (last version, http://docs.bpipe.org)
- blastn (last version)
- seqkit (last version, https://github.com/shenwei356/seqkit)
- R: readr, data.table
- Perl
- repeatmasker
- abyssis [to check]

# Updates:

14/05/2020
- debug step2
- debug step1
- upload generate_longer_reads_from_blasted_db_archaicGMTOK.R

07/05/2020
- update step2 bpipe
- add filtering_5prime_from_blast_archaic_v2.R
- add generate_longer_reads_from_blasted_db_archaic.R
- improve speed some stages step1 

28/04/2020
- added and removed stages step2_te_pipeline.bpipe

10/04/2020:
- created draft step2 pipeline: step2_te_pipeline.bpipe
- filtering_5prime_and_filled_sites_from modern_specific_3prime.R > converted in bash 

09/04/2020:
- seqkit now split the reference genome by chromosome
- changed script fasta_reference_portion2.pl, add argument (path reference), to identiy the folder with the fasta for each chromosome and extract the sequences. 

08/04/2020: 
- add stages pipeline
- deduplicated_modern_flankings_and_filled_sites.R > converted bash 
- filtering_archaic_vs_modern_3prime.R >  converted bash 

# To do list:
General: create a config file with parameters (to check). To implement usage of arguments for each function in each stage. Optimize some chunks of code of step1_te_pipeline.bpipe

07/05/2020:
- check script R stages 2 (debug)
- add stages step2 - seqkit 

10/04/2020:
- step2: optimize R scripts, add stages pipeline
- add the possibility to upload the parameters from a tab both from step1 and step2
- improvers some steps step1

07/04/2020-Now: 
- filtering_5prime_and_filled_sites _from modern_specific_3prime  DONE
- to complete step1_te_pipeline.bpipe DONE 
- deduplicated_modern_flankings_and_filled_sites.R -> to improve DONE
- filtering_archaic_vs_modern_3prime.R  -> transform bash script DONE
- filtering_5prime_and_filled_sites _from modern_specific_3prime  -> to improve
