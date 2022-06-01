Supplementary information / reproducible research files for the manuscript 
Title: "Incorporating infectious duration-dependent transmission into Bayesian epidemic models"

Authors: Caitlin Ward, Grant D. Brown, and Jacob J. Oleson

Authors of the code: Caitlin ward (caitlin.ward@ucalgary.ca)

The code was written/evaluated in R with the following software versions:

Simulations and data analysis were run in parallel on a Linux server with software versions:



This folder contains the following data and files that can be used to reproduce all analysis and 
figures of the manuscript.
It contains three subfolders containing the following files:

./ebola_analysis:
	run_models.R
	An R script that performs the analysis of the EVD data reported in the manuscript (section 4).
	The script loads the EVD data from the ABSEIR package, runs three MCMC chains for 
	each of the six models fitted, and compiles and saves summary statistics of the posterior 
	distribution, WAIC, and Gelman-Rubin convergence diagnostics. Each model was run in parallel
	on 3 cores of a linux server, but should also work on Windows. Results for each model are saved
	in ./output/ since computation of each model takes several hours. For faster computation, the
	number of iterations could be reduced by changing the niter variable.

	./output/
	A folder containing the results of run_models.R for each of the six models

	get_priors_inits.R
	An R script sourced by run_models.R which provides a function to create the initial values and
	priors used for each of the six models.

	post_processing.R
	An R script sourced by run_models.R which provides a function to take the model output from the
	three chains (posterior iterations, calculated WAIC), and combine and summarize
	aspects necessary for the final results: Gelman-Rubin, posterior summaries of parameters, 
	IDD curves, and R0(t), and WAIC.

	ebola_results.R
	An R script that takes the RDS files from the ./output/ subfolder and creates the results
	yields exactly the results presented in the manuscript (figures, tables, and in-text information)

./simulation_study:

./scalability_analysis:
	run_scalability.R
	An R script that performs the simulations to assess the scalability of models to larger epidemics
	(section 5). For increasing population sizes, epidemics are simulated, and 100 iterations of the 
	exponential, path-specific, and IDD models are timed using microbenchmark. This code yields 
	exactly the results reported in Figure 8 of the manuscript after running on the author's PC. 
	Timings may vary based on the machine used to run the script, but the trends seen in Figure 8 
	should not change. After timing each model on the increasingly large epidemics, the results
	are saved in the current folder as timePerIter.rds.

	scalability_results.R
	An R script that imports the RDS produced by run_scalability.R and creates Figure 8 from the
	manuscript.
	
	

