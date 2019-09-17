# TEAM
The `R` package `TEAM` and source file TEAM2D.R identify local 1D (one marker) and 2D regions (two markers), respectively, that differ between case and control distributions.

`TEAM` accompanies the paper:
   
    Pura J., Chan C., Xie J. 2019. Multiple Testing Embedded in an Aggregation Tree to Identify 
    where Two Distributions Differ. https://arxiv.org/abs/1906.07757

## Installation
You can install `TEAM` from CRAN with:

    install.packages("TEAM")

Please see the `TEAM` reference guide in CRAN for details.

TEAM2D.R takes in two dataframes from case (df1) and control (df0). Each dataframe contains two columns, one for each marker. The user inputs the number of breaks, m, to be created in each direction, as well as the FDR level, alpha, and number of layers, L. To use `TEAM2D.R` simply source the file in the working environment.

The Simulation folder contains the code to reproduce the simulations in the paper above.

The EQAPOL folder contains the code to reproduce the data analysis on the 11 individuals, as mentioned in the paper above. The case (CMV pp65) and control (Costim) flow cytometry samples for the 11 individuals can be accessed via https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/D141FU.
