# TEAM
The functions TEAM.R and TEAM2D.R identify local 1D (one marker) and 2D regions (two markers) that differ between case and control distributions.

TEAM.R takes in two vectors of observations from case (x1) and control (x0). Each vector corresponds to the intensities for a single marker. Additionally, the procedure takes in 2^K bins (default K=3) and L layers (default 3). The FDR level is set by alpha (default 0.05)

TEAM2D.R takes in two dataframes from case (df1) and control (df0). Each dataframe contains two columns, one for each marker. The user inputs the number of breaks, m, to be created in each direction, as well as the FDR level, alpha, and number of layers, L.

The EQAPOL folder contains the code for the data analysis on 11 individuals. The case and control flow cytometry samples for the 11 individuals can be accessed via https://duke.box.com/s/vex1m5ogwyt966adldb27tk5snslqns8
