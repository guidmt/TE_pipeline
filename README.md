# TE_pipeline
 
# How to run:
1) Install bpipe
2) Configure the parameters of *.bpipe files 
3) Run:

```bash
bpipe  step1_te_pipeline.bpipe
```

![Denisova](https://upload.wikimedia.org/wikipedia/commons/thumb/d/da/Известная_на_весь_Мир_Денисова_пещера._01.jpg/1920px-Известная_на_весь_Мир_Денисова_пещера._01.jpg=200x200)

# Dependencies:
- bpipe (last version, http://docs.bpipe.org)
- blastn (last version)
- seqkit (last version, https://github.com/shenwei356/seqkit)
- R: readr, data.table
- Perl
- repeatmasker
- abyssis

# Update:
07/04/2020: 
- add stages pipeline
- deduplicated_modern_flankings_and_filled_sites.R > converted bash 
- filtering_archaic_vs_modern_3prime.R >  converted bash 

# To do list:
General: create a config file with parameters (to check). To implement usage of arguments for each function in each stage. Optimize some chunks of code of step1_te_pipeline.bpipe

08/04/2020:
- to complete step1_te_pipeline.bpipe
- filtering_5prime_and_filled_sites _from modern_specific_3prime  -> to improve

07/04/2020: 
- to complete step1_te_pipeline.bpipe
- deduplicated_modern_flankings_and_filled_sites.R -> to improve DONE
- filtering_archaic_vs_modern_3prime.R  -> transform bash script DONE
- filtering_5prime_and_filled_sites _from modern_specific_3prime  -> to improve
