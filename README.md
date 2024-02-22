#### These scripts and files are used to analyze fitness and gas exchange measurements from the manuscript "Elevated [CO2] and temperature augment gas exchange while depressing the fitness of a montane forb"
* fit.csv contains fitness data from the growth chamber experiment
* phys.csv contains physiological measurements obtained with the LICOR
* Sco_ecophys.txt contains LICOR data for the field data in Colorado

Gas_Exchange_updated_glmmTMB.R uses phys.csv and Sco_ecophys.txt

fitness_updated.R uses fit.csv

- - - - 
Updated: 2/22/24 removed old files and old R script. Previous Rscripts used lme4 to analyze data. New R files use glmmTMB, which allows us to pull confidence intervals for ggplot more effectively. 
