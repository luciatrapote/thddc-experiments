README: Reproducibility of tHDDC Experiments

This repository contains the source code, datasets, and scripts necessary to reproduce the experimental results presented in the paper. The experiments are organized into two main blocks: Synthetic Data Simulations (Section 5) and Real Data Experiments (Section 6).

1. Synthetic Data Simulations (Section 5)
This section evaluates the performance of the proposed tHDDC algorithm under controlled environments across three distinct scenarios. 
To assess its effectiveness, tHDDC is compared against two key benchmarks:

	-tclust: A well-established method for robust clustering that handles outliers but may struggle with the "curse of dimensionality" in very high-dimensional spaces.

	-HDDC: The state-of-the-art for parsimonious high-dimensional clustering. While highly efficient in high dimensions, it lacks inherent robustness against contamination and outliers.

The goal of these simulations is to demonstrate that tHDDC bridges the gap between these two approaches, offering both high-dimensional efficiency and robust performance.

Simulation Scripts (.R)
	sim_scen_1.R: Script for Scenario 1. It includes the data generation process and the code to reproduce Table 1 (first row) and Figures 3(a), 4(a), 5(a), 6(a), and 7(a).

	sim_scen_2.R: Script for Scenario 2. It includes the data generation process and the code to reproduce Table 1 (second row) and Figures 3(b), 4(b), 5(b), 6(b), and 7(b).

	sim_scen_3.R: Script for Scenario 3. It includes the data generation process and the code to reproduce Table 1 (third row) and Figures 3(a), 4(a), 5(c), 6(c), and 7(c).

Data Objects (.RData)
To ensure exact reproduction of the figures in the manuscript, the original R data objects (generated with fixed seeds) are provided:

	sim_scen_1.RData: Pre-generated dataset for Scenario 1.

	sim_scen_2.RData: Pre-generated dataset for Scenario 2.

	sim_scen_3.RData: Pre-generated dataset for Scenario 3.

2. Real Data Experiments: Digits Data (Section 6)
This section illustrates the robustness of the tHDDC method using a modified version of the United States Postal Service (USPS) digit recognition dataset. The experiment compares tHDDC against the parsimonious state-of-the-art HDDC and the robust benchmark tclust.

Dataset: data_ZZZ_df.csv: A specialized high-dimensional dataset derived from the USPS digits, consisting of:

		- 1,756 original observations: Hand-written digits 3, 5, and 8 (collectively known as the USPS358 subset, available in the MBCbook R package).

		- 195 sampled outliers: Observations randomly selected from the USPS dataset (available in the Rdimtools R package) consisting of digits other than 3, 5, and 8; specifically, observations from digits 0, 2, 4, and 9 were introduced to act as natural outliers.

		- 45 artificial outliers: Generated following specific geometric patterns as detailed in Section 6.

Scripts and Results
	sim_digits_data.R: This script generates the confusion matrices and performance metrics presented in Table 2.

		Note on Reproducibility: The script utilizes the HDclassif package to implement the standard HDDC algorithm. Due to the internal stochastic nature of certain functions within HDclassif, the HDDC results may show minor numerical variations upon re-run. However, results for the proposed tHDDC and other methods are fully reproducible.

	fig_digits_data.R: Dedicated script to reproduce all figures and visualizations presented in Section 6.

	tHDDC_results_uX.rds: A saved R object containing the final classification and parameter estimation results for the tHDDC method, provided for direct inspection and figure generation without requiring a full re-run of the algorithm.


3. Core Algorithms and Computational Efficiency
In addition to the simulation and figure scripts, this repository provides the source code for the algorithms in their parallelized versions.

Algorithm Scripts (.R)
	tHDDC_paralelizado.R: Implementation of the proposed tHDDC algorithm.

	tclust_paralelizado.R: A parallelized version of the tclust algorithm.

Two-Stage Strategy for Concentration Steps
To ensure a fair comparison of execution times during the synthetic simulations and to optimize computational resources, both scripts implement a two-stage concentration strategy:

	- Stage 1 (Initialization): The algorithm begins with a large number of random starts (nstart). Each start is subjected to only a few preliminary concentration steps (cstep1, e.g., 2 or 4).

	- Stage 2 (Refinement): Only the nkeep most "promising" initializations—those yielding the highest values of the objective function or BIC—are retained. These candidates are then iterated until convergence or until a maximum number of final concentration steps (cstep2) is reached.

This approach significantly improves efficiency by avoiding the full iteration of poor initializations. It allows for a more thorough exploration of the parameter space within the same computational budget, ensuring that the performance gains of tHDDC are assessed alongside its practical scalability.

Software Requirements:
To execute these scripts, ensure the following R packages are installed:

	-tclust

	-HDclassif

	-foreach / doParallel (for parallel processing)

	-ggplot2 / patchwork (for visualization)

	-clue (for cluster agreement metrics)

	-doRNG (for reproducible parallel seeds)
