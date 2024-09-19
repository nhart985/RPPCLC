# RPPCLC

R implementation of the RPPCLC model proposed by Hartman, N. and He, K.

Scientific Registry of Transplant Recipients
https://www.srtr.org/about-the-data/technical-methods-for-the-program-specific-reports/

Link to Datasets:

https://www.srtr.org/reports/program-specific-reports/

Nu_Estimation.R: R script for running simulations to estimate the cluster-level confounding coefficient using our robust privacy-preserving cluster-level confounding model. 
    
CRE.R: R script for running simulations to compare our robust privacy-preserving cluster-level confounding model with the patient-level robust correlated random effects model. 

Flagging.R: R script for running simulations to compare our proposed inference methods (Frequentist and Pseudo-Bayesian) with naive versions that do not account for cluster-level confounding. 
    
Real.R: R script for evaluating United States kidney transplant programs using summary statistics from the Scientific Registry of Transplant Recipients and the proposed cluster-level confounding adjustments.
