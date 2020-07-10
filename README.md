##################################################################
################### READ ME for code for #########################
################### Gupta, Vats:   ###############################
############## Estimating Monte Carlo variance from multiple Markov chains######################
##################################################################

There is a running script "minimum.R" that produces all the 
Figures (1-6) and produces the contents of Tables 1 and 2. However, 
the script "minimum.R" loads objects produced by other codes
that take a while to run. Thus, the output files are provided as well in this folder. 
These output files are copies of the output files inside each example folder.

####
"multinomial_running_sqroot_2_m_2_n_10000_plot_step_100_.Rdata" is produced by "Multinomial/multinomial.R". This may take about
5 hours to run.

####
"MVG_CDE_rho_0.5_m_5_n_10000_.Rdata", "MVG_CDE_rho_0.5_m_10_n_10000_.Rdata", "MVG_CDE_rho_0.999_m_5_n_10000_.Rdata" and 
"MVG_CDE_rho_0.999_m_10_n_10000_.Rdata" are produced by "multivariate_Gibbs_normal/cove_det_ess.R" This may take about 40 minutes to run.

####
"MVG_running_rho_0.5_m_5_n_10000_.Rdata", "MVG_running_rho_0.5_m_10_n_10000_.Rdata", "MMVG_running_rho_0.999_m_5_n_10000_.Rdata" and 
"MVG_running_rho_0.999_m_10_n_10000_.Rdata" are produced by "multivariate_Gibbs_normal/multivariate_Gibbs_normal.R" This may take about 50 minutes to run.

####
"rosenbrock_CDE_8_8_m_5_n_1e+05_.Rdata" and "rosenbrock_CDE_8_8_m_10_n_1e+05_.Rdata" are produced by "rosenbrock/cove_det_ess.R". This may take about 3 hours to run.

####
"rosenbrock_running_ESS_8_8_m_5_n_1e+05_.Rdata" and "rosenbrock_running_ESS_8_8_m_10_n_1e+05_.Rdata" are produced by "rosenbrock/rosenbrock.R". This may take about 7 hours to run.


