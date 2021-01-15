# BluegillShapeData
DOC-Geometric Morphometrics Analyses

### Abstract
This repo contains code and data related to Chelsea Bishop's project on Bluegill morphometrics.

### Directory structure

```
.
├── README.md
├── watershedCalculations/  # code and files used by Stuart to assign lakes to watershed basins
    ├── hucIDforLakes.R
    ├── Lake_Info_2020.csv
    ├── Lake_Info_2020.xlsx
    └── Lake_Info_2020wBasins.csv
├── Bishop_NotoriousBLG.xlsx  #excel file containing all capture and merisitic data that is included in MFE database
├── Code For SVG Images for Inkscape.R   # summarized code for all analyses to create figures. May just be repeat of what is already in larger analyses R files so might delete. 
├── Code_and_Results_BLG_DOC_v2020Mar15.Rmd  # code for geometric morphometrics analysis of full body images, includes getting eye width values from interlandmark distance of landmarks on eye. Also includes how to create backtransform morphospace figure from Olsen, A.M. 2017 paper in Functional Ecology. Some model trials as well. 
├── FULL_2018_TPS_FILE_UPDATED_09-25-19.TPS  # tps file of landmark x,y coordinates of 19 landmarks placed on each photo of sampled fish. To be uploaded in geometric morphometrics analyses
├── FULL_body_links.txt.  # tells geomorph in R how to connect landmarks. To be uploaded in geometric morphometrics analyses
├── Gill_Rakers_2018_Final.csv 
├── Identifiers_Update_2020.txt  # identifiers for geometric morphometrics analyses in code_and_results..
├── Identifiers_Update_2020.xlsx
```

| Data file      | Description | Output from | Input to |
| -------------- | ----------- | ----------- | -------- |
| angle_var.csv  | | |R_UnivariateAnalyses.Rmd|
| Bishop_NotoriousBLG.xlsx ||| |
| eye_variance.csv||||
| eyedata2020_modelSE_and_intercepts.csv||||
| eyewidthsFINAL.csv ||||
| findata2020_modelSE_and_intercepts.csv ||||
| FULL_2018_TPS_FILE_UPDATED_09-25-19.TPS ||||
| Full_body_links.txt ||||
| Gill_Rakers_2018_Final.csv ||||
| Identifiers_Update_2020.txt ||||
| Identifiers_Update_2020.xlsx ||||
| MFEdb_20201125.db | MFE database from which some of the data is pulled. ||||
| Pec Fin Angles.csv ||||
| PecFinDataNovember.csv ||||
| PecFins2018.TPS ||||
| pecL-w_variance.csv ||||
| proc_shape_variance.csv ||||
| raker_variances.csv ||||
| rakerdata2020_modelSE_and_intercepts2.csv ||||
| rakerLength_variances.xlsx ||||
| Ratio Identifiers Update.txt ||||
| Ratio Identifiers Update.xlsx ||||
| Replicates Identifiers.txt ||||
| Total Replicates Coords.TPS ||||
| TPS_Identifiers_Final_October.txt ||||


### Contributions
Repository created by Chelsea Bishop and Kaija Gahm in 2020. Data were collected by Chelsea Bishop with the help of Alex Ross, Colin Dassow, Matthew Farragher, Henry Chung and undergrads of the University of Notre Dame in 2018. Geometric Morphometrics Analyses were conducted with the help of Madlen Stange. Data will eventually be contributed to the [MFE database](https://figshare.com/articles/MFE_database_Data_from_ecosystem_ecology_research_by_Jones_Solomon_and_collaborators_on_the_ecology_and_biogeochemistry_of_lakes_and_lake_organisms_in_the_Upper_Midwest_USA/7438598). 
