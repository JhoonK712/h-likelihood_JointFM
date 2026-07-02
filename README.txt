README

The code is an R-based program to implement model fitting and computation.
The model is an h-likelihood-based model used for analyzing various types of clustered bivariate survival data (i.e., bivariate censored data, semi-competing risks data, and competing risks data).
The baseline hazards for each type of bivariate survival data can be visualized in figures, as the M-spline technique is employed.
The estimated values and 95% CIs for each frailty can be identified.
Real-world datasets (i.e., kidney data, bladder cancer data, and ovarian cancer data) are used as examples.

Detailed information on the R code and its execution can be found below.

-------------------------------------------------------------------------------------------------------------------------------

* Information on the R version and all packages used to produce these results is provided in sessionInfo().txt.
* The analysis relies on a master code and source files containing functions. Please save and run all R scripts in the same folder.

1. Execution Guide
	- Enter the folder path of the R scripts into the setwd() function in Master_code.R, then execute the source() commands for each table and figure.
	- After configuring the path via setwd(), you can also run each table or figure R script directly.
	* Defining the path in Master_code.R once eliminates the need for any additional path adjustments.

2. Guide to Key R Codes
	- Master_code.R: This is the master execution script designed to replicate every table and figure in the paper.
	- Source_Simulation.R: A source file dedicated to running simulations.
		- The corresponding outputs from this script include: 
			- Table 2
			- Table B.1
	- Source_Realdata.R: A source file for real data analysis.
		- It incorporates varying settings for the association parameter alpha (alpha = 1, 0, -1$).
		- It features predefined setups for kappa1 and kappa2 (kappa_1, kappa_2 = 1, 1e+8, 1e+15) for sensitivity analysis purposes.
		- The corresponding outputs from this script include: 
			- Table 4, Table 5, Table 6
			- Table C.1, Table C.2, Table C.3, Table C.4, Table C.5, Table C.6
			- Figure 1, Figure 2, Figure 3
			- Figure C.1, Figure C.2, Figure C.3

3. Description of Other R Codes
	- Regarding the R codes below, all tables are generated as Excel files in a format similar to the paper, and all figures are saved as PNG files.
		- Table 2, Table B.1.R - The simulation results on the proposed estimation method over 200 replications on three types of clustered bivariate survival data (Table 2), as well as the number of errors occurred during 200 simulation replications for clustered bivariate survival data (Table B.1), are estimated.
		- Table 4.R - The results for each package were estimated using the kidney infection data.
		- Table 5.R - The results for each package were estimated using the ovarian cancer data.
		- Table 6.R - The results for each package were estimated using the bladder cancer multicenter data.
		- Table C.1.R - This is to conduct a sensitivity analysis on how the results of the proposed method change as alpha varies in the kidney infection data.
		- Table C.2.R - This is to conduct a sensitivity analysis on how the results of the proposed method change as alpha varies in the ovarian cancer data.
		- Table C.3.R - This is to conduct a sensitivity analysis on how the results of the proposed method change as alpha varies in the bladder cancer multicenter data.
		- Table C.4.R - This is to compare how the results of the proposed method change as kappa1 and kappa2 vary in the kidney infection data.
		- Table C.5.R - This is to compare how the results of the proposed method change as kappa1 and kappa2 vary in the ovarian cancer data.
		- Table C.6.R - This is to compare how the results of the proposed method change as kappa1 and kappa2 vary in the bladder cancer multicenter data.
		- Figure 1, Figure C.1.R - Visual representation of frailties for each patient (Figure 1) and the baseline hazard function (Figure C.1) for each event under the proposed method using the kidney infection data.
		- Figure 2, Figure C.2.R - Visual representation of frailties for each group (Figure 2) and the baseline hazard function (Figure C.2) for each event under the proposed method using the ovarian cancer data.
		- Figure 3, Figure C.3.R - Visual representation of frailties for each centre (Figure 3) and the baseline hazard function (Figure C.3) for each event under the proposed method using the bladder cancer multicenter data.
	













