# EQAPOL Case Study

This subdirectory contains all the functions required to replicate the data analysis on the EQAPOL datasets.
The datasets consist of paired comma-separate-value (CSV) files (CMV_pp65 and Costim) for 11 individuals. 
Running the `R` functions above will generate Figures 4 and 5, and values in Section 6 in the main text, as well as, 
Table S1 and Figures S1-S9 in the Supplementary Material of the accompanying paper: 
Pura, Chan, Xie (2019) Multiple Testing Embedded in an Aggregation Tree to Identify
where Two Distributions Differ (https://arxiv.org/pdf/1906.07757.pdf)

To reproduce the above items, 1) install the `TEAM` package from `CRAN` : `install.packages("TEAM")`, 2) source `case_study_functions.R`, and 3) run the code in `run_EQAPOL_case_study.R`.
Additional details are contained in the header of each file. 

