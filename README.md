
# Code and Data for: "Age-Based Risk Estimates for *C9orf72*<sup>RE</sup>-related Diseases: Theoretical Developments and Added Value for Genetic Counselling"

This repository contains the R code and data accompanying the paper **"Age-Based Risk Estimates for *C9orf72*<sup>RE</sup>-related Diseases: Theoretical Developments and Added Value for Genetic Counselling"** by D.M. de Vienne and D. de Vienne. 

The aim of this study is to provide a quantitative framework for genetic counselling in ALS/FTD families, by estimating carrier probability and disease risks in the next years for relatives of carriers of the *C9orf72*<sup>RE</sup> mutation, based on age-dependent penetrance; 

An online simulator implementing these calculations is available to genetic consellors online at https://lbbe-shiny.univ-lyon1.fr/ftd-als/.

## Contents

- `article/` contains data and code related to the article: 
  - **`Age-dep-proba.R`**: R script for data processing, probability computations, and visualization of age-dependent carrier probabilities.
  - **`C9_Data.csv`**: Copy of the dataset from Murphy et al. (2017).
- **`app.R`**, **`report-template.Rmd`**, **`rinstall.txt`**, **`df_onset.csv`** and the `www/` folder concern the online *shiny* application. 

## Citation

If you use this code, please cite the accompanying paper:

> D.M. de Vienne & D. de Vienne. 2025. Improving Genetic Counselling for C9orf72-mediated ALS/FTD with Age-Based Risk Estimates. *MedrXiv*
>  

## References

Murphy NA, Arthur KC, Tienari PJ, Houlden H, Chiò A, Traynor BJ. Age-related penetrance of the C9orf72 repeat expansion. Sci Rep. 2017;7(1):2116. Published 2017 May 18. doi:10.1038/s41598-017-02364-1

## License

This repository is distributed under the [MIT License](LICENSE).