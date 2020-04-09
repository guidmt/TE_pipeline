# TE_pipeline
 
# How to run:
1) Install bpipe
2) Configure the parameters of *.bpipe files 
3) Run:

```bash
bpipe  step1_te_pipeline.bpipe
```

# Dependencies:
- bpipe (last version, http://docs.bpipe.org)
- blastn (last version)
- seqkit (last version, https://github.com/shenwei356/seqkit)
- R: readr, data.table
- Perl
- repeatmasker
- abyssis

# Update:
- add stages pipeline
- deduplicated_modern_flankings_and_filled_sites.R > converted bash command

# To do list:
General: create a config file with parameters (to check). To implement usage of arguments for each function in each stage. Optimize some chunks of code of step1_te_pipeline.bpipe

07/04/2020: 
- Complete step1_te_pipeline.bpipe
- deduplicated_modern_flankings_and_filled_sites.R -> to improve
- filtering_archaic_vs_modern_3prime.R  -> transform bash script
- filtering_5prime_and_filled_sites _from modern_specific_3prime  -> to improve
