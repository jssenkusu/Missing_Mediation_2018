R files reproduce results in the paper "Missing data and mediation analysis in logitudinal studies"

***** FILES

datagen_mix_mar_20180213.R :	Generates data from a linear mixed model
diffthensim_20180213.R     :	Imputes missing data by first taking the difference
				M_ij - M_i1 and then imputes missing values.
simthendiff_20180213.R     :	Imputes missing data by first imputes missing values
				in M_ij and Y_ij and then takes the difference
				M_ij - M_i1.
est_tables_20180213.R      :	Summarizes results and generates syntaxt that can be
				compiled in latex 

***** TO RUN THE SIMULATION STUDY

- Make sure the R scripts are all in one directory.
- Run the imputation R scripts. Point to the directory where the data generation 
  R script in located. Make sure the output is directed to the 
  same directory where other scripts are located. 
- This may take several hours.
- The results will be saved in your directory.
- Run the R script that summarizes the results. It will generate latex
  syntax that can be compiled to view the results.  

***** R PACKAGES NEEDED

foreach, doParallel, utils, MASS, pan, mice, mitools, lme4, mlmmm, 
reshape, simglm, xtable