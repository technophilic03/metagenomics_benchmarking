# This file is for generating execution time plots

## Bracken MISSING

## Centrifuge
## More files than there are results, also why are there multiple execution times
# In no error folders, the correct files are the _noerr
ctf_time <- list.files("profilers/Centrifuge/outputs", 
                       pattern = "timing.txt", recursive = TRUE,
                       full.names = TRUE)


krk2_time <- list.files("profilers/Kraken2/outputs", 
                       pattern = "timing.txt", recursive = TRUE,
                       full.names = TRUE)

ps <- list.files("profilers/PathoScope2/outputs", 
                 pattern = "timing.txt", recursive = TRUE,
                 full.names = TRUE)
