![](https://github.com/MFEh2o/BluegillMorphologyDataDOC/blob/main/resources/headerFish.jpg)

# Bluegill Shape Data
Relating lake DOC to Geometric Morphometric Analyses of bluegill captured from study lakes.

### Abstract
This repository contains code and data related to Chelsea Bishop's project on Bluegill morphometrics. The data can be found in the data/ folder, apart from the MFE database, which must be downloaded separately (see **Project setup**, below). Code is in the code/ folder. The [renv](https://rstudio.github.io/renv/articles/renv.html)/ folder preserves package versions used in this analysis, including older versions of packages that may not be current on CRAN. The resources/ folder has various relevant literature documents diagrams/notes. The figures/ folder contains the raw material for figures, which were then modified in Photoshop by CB before inclusion in the manuscript.

### Project setup
1. Clone this repository to your computer using the green Download Code button.
2. **Download MFE database** The MFE database, which these analyses draw from, is too large to be stored on GitHub. It is stored in a [Figshare repository](https://caryinstitute.figshare.com/articles/dataset/MFE_database_Data_from_ecosystem_ecology_research_by_Jones_Solomon_and_collaborators_on_the_ecology_and_biogeochemistry_of_lakes_and_lake_organisms_in_the_Upper_Midwest_USA/7438598/5). This analysis draws from MFEdb_20210423, which is version 4.5.4 of the MFE database and version 5 on Figshare (the previous link will take you to the correct Figshare version). There are lots of other supporting documents on the Figshare repository, but you only need to download the one called MFEdb_20210423.db. Once you've downloaded the database file, *Put this file into the "data/" folder of the BluegillMorphologyDataDOC cloned folder.*
3. **Download the NHD watershed shapefiles** These are used to assign major watershed basins to each study lake, by means of a spatial intersection with lake lat/long coordinates. You can download the shapefiles from [this Figshare repository](https://caryinstitute.figshare.com/articles/dataset/Morphometry_of_Bluegill_sunfish_Lepomis_macrochirus_varies_with_lake_dissolved_organic_carbon_concentration/14529303). In that repo, the README.txt file describes where the data comes from. Then, there are two zip files, "region4.zip" and "region7.zip". Download both of those and unzip them; you will now have two folders by the same names. *Move the "region4/" and "region7/" folders into the "data/inputs/" folder of the BluegillMorphologyDataDOC cloned folder.*
4. Make sure you have both R and RStudio installed. If you don't, install them here: [R](https://www.r-project.org/) [RStudio](https://www.rstudio.com/products/rstudio/download/)
5. Open the RStudio project, "BluegillMorphologyDataDOC.Rproj". This will open a new session of RStudio. All the scripts in this project are written with file paths that will work within this "project" setup. For more information on how RStudio projects work, see [this article](https://support.rstudio.com/hc/en-us/articles/200526207-Using-Projects).
6. This project is built with the `renv` package to manage package versions, to make sure that the code will still run even when the pacakges get upgraded in the future. Therefore, when you first open the RStudio project, you will probably see a message like this in the console (the example has a file path specific to Chris's computer--you would see a different file path):

```
* Project 'C:/Users/solomonc/OneDrive/Documents/Research projects/Microbes to Micropterus/Bishop bluegill morphometry/BluegillMorphologyDataDOC' loaded. [renv 0.12.5]
* The project library is out of sync with the lockfile.
* Use `renv::restore()` to install packages recorded in the lockfile.
```
As instructed, type `renv::restore()` in the console to install and load the R package versions that are needed for this project.

If you have specific R package versions installed elsewhere on your computer, don't worry--`renv` is set up to store local copies of package versions in the project folder itself. So all of your other R code in other projects on your computer will not be affected. For more information about how `renv` works, see [this article](https://rstudio.github.io/renv/articles/renv.html).

*Heads up*: if you're using a Windows computer, especially if you haven't used `renv` before, you might run into an error like this when you run `renv::restore()`:

```
Installing Matrix [1.3-2] ...
    FAILED
Error installing package 'Matrix':
==================================
* installing to library 'C:/Users/solomonc/AppData/Local/Temp/RtmpI9oIou/renv-staging-2eb03eb0751d'
* installing *source* package 'Matrix' ...
** package 'Matrix' successfully unpacked and MD5 sums checked
** using staged installation
** libs
Warning in system(paste(cmd, "shlib-clean")) : 'make' not found
Warning in system(cmd) : 'make' not found
ERROR: compilation failed for package 'Matrix'
* removing 'C:/Users/solomonc/AppData/Local/Temp/RtmpI9oIou/renv-staging-2eb03eb0751d/Matrix'
Error: install of package 'Matrix' failed
```
This error usually means that your computer isn't configured with the necessary tools to install packages from source, instead of from pre-compiled binaries, which `renv` needs to do if you're trying to install an older package version. To fix this, install the `Rtools` package from [this site](https://cran.r-project.org/bin/windows/Rtools/) (make sure to choose the appropriate link, according to whether your computer has a 32-bit or 64-bit processor).

Once you've installed Rtools, try `renv::restore()` again. You may have to restart R and RStudio first.

7. Run the analyses in order, from the beginning. 

Here is a diagram illustrating the reproducible data and code pipeline (click to make it larger). The pinkish ovals are R scripts that need to be run. So, start by running compileLakeDOC.R; then run hucIDforLakes.R, etc. (You don't have to run chelseaDataEnter_FISH_INFO_cleanup_gh92_gh29_gh93.R, shown in a lighter pink).

![](https://docs.google.com/drawings/d/e/2PACX-1vTEP9I5EuLxzX9hGX9sKzH35NbnyMQLg7ndA6maboKz3uW1_UmA13QyY7cssFAbMCt5Q2UDPTgbF9kv/pub?w=4116&h=947)

### Directory structure
- code/ contains final R code files, represented in the schema above.
- data/ contains final data files, including minimal necessary inputs and outputs. 
     - This analysis draws from MFEdb_20210423.db, which is version 4.5.4 of the database. See above for download instructions.
- figures/ contains saved images that are either final or near-final versions of the figures needed for the manuscript. These are generated in **analysis.R**. When further editing was needed, these were loaded into Photoshop and modified. Because some figures contain multiple parts and I wanted to leave room for arranging those parts as needed, I've created a separate folder for each figure, containing one or multiple image files. The 'allometry/' subfolder contains some allometry plots created in the various **create....R** scripts; these were not included in the final manuscript.
- renv/ contains data about package versions used in preparing this project, to facilitate re-running the analyses at a later date even if packages have been updated.
- resources/ contains various diagrams and literature that were helpful in performing this analysis.

```
.
├── BluegillMorphologyDataDOC.Rproj
├── README.md
├── code
│   ├── analysis.R
│   ├── compileLakeDOC.R
│   ├── createEyeWidths.R
│   ├── createGillRakers.R
│   ├── createIdentifiersUpdated.R
│   ├── createPecFinAngles.R
│   ├── createPecFins.R
│   ├── dbUtil.R
│   ├── defs.R
│   ├── hucIDforLakes.R
│   └── plotFishLateral.R
├── data
│   ├── MFEdb_20210423.db (** must download separately; see above)
│   ├── inputs
│   │   ├── FULL_2018_TPS_FILE_UPDATED_09-25-19.TPS
│   │   ├── Full_body_links.txt
│   │   ├── TotalReplicatesCoords.TPS
│   │   ├── fishBodyPhotos_fileNames.txt
│   │   ├── region4/ (** must download separately; see above)
│   │   ├── region7/ (** must download separately; see above)
│   │   └── ntl41_v1_0.csv
│   └── outputs
│       ├── Gill_Rakers_2018_Final.csv
│       ├── Identifiers_Update_2020.txt
│       ├── Lake_Info_2020wBasins.csv
│       ├── PecFinAnglesFINAL.csv
│       ├── PecFinDataNovemberFINAL.csv
│       ├── ReplicatesIdentifiers.txt
│       ├── eyewidthsFINAL.csv
│       ├── lakeInfo.csv
│       └── univariateModelSummary.csv
├── figures
│   ├── allometryPlots/
│   ├── eyeWidths/
│   ├── fishShapes_pc1_pc2/
│   ├── gillRakers/
│   ├── gradientLegend.pdf
│   ├── legend.pdf
│   └── pecFins/
├── renv/
├── renv.lock
└── resources/
```
### Raw data sources (files in data/ and data/inputs/)
- **MFEdb_20210423.db** See download instructions above.
- **FULL_2018_TPS_FILE_UPDATED_09-25-19.TPS** measurements taken in tpsDIG and exported
- **Full_body_links.txt** created manually by Chelsea Bishop, based on knowledge of how the body landmarks interact. This file serves as a precursor to geometric morphometric analyses.
- **TotalReplicatesCoords.TPS** measurements taken in tpsDIG and exported
- **fishBodyPhotos_fileNames.txt** created by running `ls()` on the folder "fishBodyPhotos", stored in the MFE Box drive under "MFE/notoriousBLGPhotoArchive".
- **region4 and region7 shapeiles** See download instructions above.
- **ntl41_v1_0.csv** Downloaded from [here](https://lter.limnology.wisc.edu/dataset/biocomplexity-north-temperate-lakes-lter-coordinated-field-studies-chemical-limnology-2001-2) (the "Download csv" option)

### Contributors
Repository created by Chelsea Bishop and Kaija Gahm in 2020-21. 

Data were collected by Chelsea Bishop with the help of Alex Ross, Colin Dassow, Matthew Farragher, Henry Chung and undergrads of the University of Notre Dame in 2018. Data and scripts were reorganized and reviewed by Kaija Gahm in consultation with Chelsea Bishop, and Kaija prepared the data reproducibility pipeline.

Geometric Morphometrics Analyses were conducted with the help of Madlen Stange. 
Fish morphometric data are contributed to the [MFE database](https://figshare.com/articles/MFE_database_Data_from_ecosystem_ecology_research_by_Jones_Solomon_and_collaborators_on_the_ecology_and_biogeochemistry_of_lakes_and_lake_organisms_in_the_Upper_Midwest_USA/7438598). 
