# READ-ME / OVERVIEW

The data and code provided here accompay the paper "Ensuring identifiability in hierarchical mixed effects Bayesian models" 
by Kiona Ogle and Jarrett Barber, published in Ecological Appilcations. The following R files are provided:

1) "Data S1.RData" -- this contains the simulated data that are used with the models contained within "Code A2.R" and "Code A3.R"
2) "Code A2.R" -- R script code for implementing the random effect regression summarized in Equation 5 in the manuscript, which loads and prepares data and implements the model(s) in JAGS. The code does not implement any reparameterization or computational solutions to solving the non-identifiabilty of the overall intercept and random effects terms, but it does include options for specifying a relatively non-informative Gamma prior or a weakly informative folded-Cauchy prior for the random effects variance-related parameter. The model is applied to the three datasets contained in the “Data S1.RData” object.
3) "Code A3.R" -- R script code for implementing Solutions 1 (hierarchical centering, HC), 2 (sum-to-zero, SZ), and 4 (post-sweeping, PS) with the random effect regression summarized in Equation 5 in the manuscript.
4) "Code A4.R" -- R script code for implementing Solutions 1 (HC), 2 (SZ), and 4 (PS) with the nested random effects regression summarized in Equation 11 in the manuscript.
5) "Code A5.R" -- R script code for implementing Solutions 2 (SZ) and 4 (PS) with the Bayesian regression involving crossed random effects, as summarized in Equation 14 in the manuscript.
