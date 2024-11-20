## These scripts and files are used to analyze fitness and gas exchange measurements in Boechera stricta from the manuscript "Elevated [CO2] and temperature augment gas exchange while depressing the fitness of a montane forb"

Manuscript is available at https://nph.onlinelibrary.wiley.com/doi/full/10.1111/nph.19765

## Table of Contents 
* [Datasets](#Datasets)  
   * [Fitness](#Fitness_gc)
   * [Physiology](#Phys_gc)
   * [Garden Physiology](#garden_phys)
* [Scripts](#Scripts)

### Datasets <a name="Datasets"></a>
The following are summary files from the growth chamber or field datasets. Metadata is listed within each expansion.

#### Growth chamber fitness <a name="Fitness_gc"></a>
 [Fit.csv](fit.csv) contains summarized fitness data from the growth chamber experiment  
  <details> 
    <summary> Growth chamber fitness metadata</summary>
    
  Column | Description
  -------| ----------
  Treatment | T1-T4 are codes to represent one of the 4 combinations of CO2 and temperature treatments
  Temperature | High - future projected temperatures; Low - contemporary based on average data from Crested Butte, CO.
  CO2 | High - future 650 ppm; Low - contemporary (~450 ppm due to scrubbing constraints; see manuscript for details)
  Tray | Block number ranging from 1-8
  Number | Unique ID number for each specimen
  Full_ID | ID number that incorporates treatment, tray, and number.
  Genotype | Accession of B. stricta samples; See supplementary table 1 for details on locality information around Gothic, CO
  Population | Population of origin for each genotype
  Elevation | Source elevation of each population, in meters
  Transplant_Date | Date in which seedling was transplanted into conetainers
  OD_transplant | Transplant date converted to ordinal date
  Lifetime_Fruited | Binary with 1 indicating fruits produced during the lifetime of the plant
  Overwinter_survival | Binary with 1 indicating plant survived 8 week vernalization period
  Season_survival | Binary with 1 indicating plant survived until the end of the second season when they were harvested
  Mature_silique_number | Total number of mature siliques (fruits) produced during the lifetime of the plant
  Mature_length_siliques | Total length of all mature siliques across the lifetime
  Cohort | 1 or 2 indicate germination cohort
  </details>
        
#### Growth chamber physiology data <a name="Phys_gc"></a>
[Phys.csv](phys.csv) contains leaf-level physiological measurements obtained with the LICOR 6800.   
See the [LICOR 6800 summary of symbols](https://www.licor.com/env/support/LI-6800/topics/symbols.html#gasex) for more details about the measured traits

<details>
<summary> Growth Chamber physiology metadata </summary> 
  
  Column | Description
  -------| ----------
  Treatment | T1-T4 are codes to represent one of the 4 combinations of CO2 and temperature treatments
  Temperature | High - future projected temperatures; Low - contemporary based on average data from Crested Butte, CO.
  CO2 | High - future 650 ppm; Low - contemporary (~450 ppm due to scrubbing constraints; see manuscript for details)
  Tray | Block number ranging from 1-8
  Number | Unique ID number for each specimen
  Full_ID | ID number that incorporates treatment, tray, and number.
  Genotype | Accession of B. stricta samples; See supplementary table 1 for details on locality information around Gothic, CO
  Population | Population of origin for each genotype
  Elevation | Source elevation of each population, m
  Cohort | 1 or 2 indicate germination cohort
  Licor_E | Evapotranspiration rate mol m-2 s-1
  Licor_A | Assimilation rate µmol m-2 s-1
  Licor_Ca | Ambient CO2 µmol mol-1
  Licor_Ci | Intercellular CO2 µmol mol-1
  Licor_gsw | Stomatal conductance to water vapor mol m-2 s-1
  Licor_LeafArea | Leaf area for each sample, in cm2
  Licor_dryLeafMass | Leaf sample dry mass, in g

</details>

#### Schofield garden physiology data <a name="garden_phys"></a>
[Sco_ecophys.txt](Sco_ecophys.txt) contains LICOR 6400 data for the Schofield garden at the Rocky Mountain Biological Laboratory.  
    See [LICOR 6400 Table 3.2](https://www.licor.com/env/support/LI-6400/topics/guided-tour-2-new-measurements.html) for more details on measurements and units
   <details>
     <summary> Schofield metadata </summary>
     
   Column | Description
   -------| ----------
   plant number | Unique ID for each individual
   Quad | Block in the garden
   Row | Row number within the block
   Column | Column number within the block
   Genotype | Population of the individual - equivalent to Population in other files
   Elevation | Source elevation of the population in meters
   Dist | Distance from original source population in meters
   Photosynthesis | Assimilation rate after leaf area adjustments
   Stomatal conductance | gsw 
   Ci | Intercellular [CO2]
   transpiration | Evapotranspiration rate, E
   Ci/Ca | Calculated ratio of intercellular to ambient [CO2]
   VpdL | Vapor pressure deficit based on leaf temperature
   WUE | A/E
  
 </details>   

### Scripts <a name="Scripts"></a>
Each of these scripts can be run independently from one another to produce ANOVA tables and figures from the manuscript.

[Gas_Exchange_updated_glmmTMB.R](Gas_Exchange_updated_glmmTMB.R) uses phys.csv and Sco_ecophys.txt to analyze photosynthesis, transipiration, stomatal conductance and WUE.

[fitness_updated.R](fitness_updated.R) uses fit.csv to analyze probability of reproduction, probability of survival, and fecundity. 

- - - - 
Updated: 2/22/24 removed old files and old R script. Previous Rscripts used lme4 to analyze data. New R files use glmmTMB, which allows us to pull confidence intervals for ggplot more effectively. 
