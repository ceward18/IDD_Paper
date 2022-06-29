Supplementary information / reproducible research files for the manuscript 
Title: "Incorporating infectious duration-dependent transmission into Bayesian epidemic models"

Authors: Caitlin Ward, Grant D. Brown, and Jacob J. Oleson

Authors of the code: Caitlin ward (caitlin.ward@ucalgary.ca)

The code was written/evaluated in R with the following software versions:

Simulations and data analysis were run in parallel on a Linux server with software versions:



This folder contains the following data and files that can be used to reproduce all analysis and 
figures of the manuscript.

It contains one file with helper functions:

helper_functions.R

It contains three subfolders containing the following files:

./ebola_analysis:
	run_models.R
	An R script that performs the analysis of the EVD data reported in the manuscript (section 4).
	The script loads the EVD data from the ABSEIR package, runs three MCMC chains for 
	each of the six models fitted, and compiles and saves summary statistics of the posterior 
	distribution, WAIC, and Gelman-Rubin convergence diagnostics. Each model was run in parallel
	on 3 cores of a linux server, but should also work on Windows. Results for each model are saved
	in ./output/ since computation of each model takes several hours. For faster computation, the
	number of iterations could be reduced by changing the niter variable, however, this may
	results in unconverged models and results which are not identical to those in the manuscript.

	./output/
	A folder containing the results of run_models.R for each of the six models

	get_priors_inits.R
	An R script sourced by run_models.R which provides a function to create the initial values and
	priors used for each of the six models. Priors are as specified in Supplemental Table 4.

	post_processing.R
	An R script sourced by run_models.R which provides a function to take the model output from the
	three chains (posterior iterations, calculated WAIC), and combine and summarize
	aspects necessary for the final results: Gelman-Rubin, posterior summaries of parameters, 
	IDD curves, and R0(t), and WAIC.

	ebola_results.R
	An R script that takes the RDS files from the ./output/ subfolder and yields exactly the
	results presented in the manuscript/supplemental material (figures, tables, and in-text 4
	information).

./simulation_study:
	data_generation.R
	An R scripts that is used to generated the simulated epidemics (section 3). Epidemics are 
	simulated in four scenarios: PS, IDD Peak, IDD Exp, IDD Logit. Stored data are saved
	in four RDS files in the ./data/ folder. Parameters of each model are as described in 
	Supplemental Table 1.

	./data/
	A folder containing the output of data_generation.R for each of the four data generating
	scenarios.

	get_priors_inits.R
	An R script sourced by model_fitting_knownE.R and model_fitting_estimatedE.R which provides
	a function to create the initial values and priors used for each of the six models. Priors
	are as specified in Supplemental Table 2.

	model_fitting_knownE.R
	An R script that performs the IDD analyses with known exposure times for the results 
	presented in Figures 2 and 3. Simulations are submitted and saved in batches of 100. 
	Each model was run in parallel on 3 cores of a linux server, but should also work on Windows.
	Results from each batch are saved in ./batch_output/ since computation takes several hours 
	per batch.  For faster computation, the number of models fit per bach could be reduced by
	changing the batchSize variable.
 
	model_fitting_estimatedE.R
	An R script that performs the main simulation analyses using all six models. As each model
	can take several hours to run, models are run individually. Each model was run in parallel 
	on 3 cores of a linux server, but should also work on Windows. Results from each batch are 
	saved in ./batch_output/. 

	post_processing.R
	An R script sourced by model_fitting_knownE.R and model_fitting_estimatedE.R which provides 
	a function to take the model output from the three chains (posterior iterations, calculated 
	WAIC, time to run), and combine and summarize aspects necessary for the final results: 
	Gelman-Rubin, posterior summaries of parameters, IDD curves, and R0(t), MCMC Efficiency and WAIC.

	./batch_output/
	A folder containing the outputs of model_fitting_knownE.R and model_fitting_estimatedE.R.

	simulation_results.R
	An R script that takes the RDS files from the ./batch_output/ subfolder and yields exactly 
	the results presented in the manuscript/supplemental material (figures, tables, and in-text
	information)
	

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
	
