/*************************************************************************************************************************************************************************************
VERSION CHANGES:
	CHANGE FROM VERSION 8 TO VERSION 9: VERSION 8 SAMPLES WITH REPLACEMENT SEPARATELY IN MEN AND WOMEN (bsample WITH strata OPTION) FOR SCENARIOS 0 AND 1. 
	                                    VERSION 9 DOES NOT SAMPLE SEPARATELY WITHIN STRATA FOR ANY OF THE SCENARIOS.
										
										VERSION 8 FITS THE SUBGROUP ANALYSIS MODEL ASSUMING THE SAME VARIANCE PARAMETER FOR MEN AND WOMEN.
										VERSION 9 FITS THE MODEL ONLY TO MEN ALLOWING FOR THE VARIANCE PARAMETER TO DIFFER BETWEEN MEN AND WOMEN
										
										VERSION 9 SIMPLIED THE BOOTSTRAP VARIANCE EQUATION IN LINE WITH THE MANUSCRIPT. EQUIVALENT TO THE BOOTSTRAP VARIANCE
 									              EQUATION IN VERSION 8

ARGS
	- scenario		: SITUATION WE ARE TESTING
	- sampleSize	: SIZE OF EACH DATASET
	- numberOfDatasets : NUMBER OF DATASETS IN THE SIMULATION STUDY
	- percentObserved : % OF OBSERVED WEIGHT VALUES
	- ratio	: omegaM/omegaF = ratio (ONLY NEEDED FOR SCENARIO 2)
	- distribution : CHOICE OF NON-NORMAL DISTRIBUTIONS FOR SCENARIO 4. IF SCENARIO 4 IS NOT SELECTED THEN THE DISTRUBTIONS ARE SET TO BE NORMAL.
    - B : NUMBER OF BOOTSTRAPPED DATASETS 
	- M : NUMBER OF IMPUTATIONS PER BOOTSTRAPPED DATASET

SCENARIOS
	0 - PROBABILITY DISTRIBUTIONS FOR ALL VARIABLES ARE THE SAME FOR MEN AND WOMEN
		ALL CONTINUOUS VARIABLES ARE NORMALLY DISTRIBUTED
		ONLY MEN HAVE INCOMPLETE DATA
		IMPUTATION MODEL IS CORRECT
		ANALYSIS MODEL IS CORRECT
	
	1 - PROBABILITY DISTRIBUTION FOR ALL VARIABLES ARE THE SAME FOR MEN AND WOMEN
		ALL CONTINUOUS VARIABLES ARE NORMALLY DISTRIBUTED
		MEN AND WOMEN HAVE MISSING DATA
		IMPUTATION MODEL IS CORRECT
		ANALYSIS MODEL IS INCORRECT (MORE GENERAL)
	
	2 - PROBABILITY DISTRIBUTION OF WEIGHT AND LOGINSINDEX HAVE DIFFERENT VARIANCES FOR MEN AND WOMEN
		ALL CONTINUOUS VARIABLES ARE NORMALLY DISTRIBUTED
		MEN AND WOMEN HAVE MISSING DATA
		IMPUTATION MODEL IS INCORRECT
		ANALYSIS MODEL IS INCORRECT
		
	3 - PROBABILITY DISTRIBUTION FOR MEN AND WOMEN HAVE DIFFERENT MEANS BUT THE SAME VARIANCES
		ALL CONTINUOUS VARIABLES ARE NORMALLY DISTRIBUTED
		MEN AND WOMEN HAVE MISSING DATA
		THE TRUE ASSOCIATION OF WEIGHT AND LOGINSINDEX ARE THE SAME IN MEN AND WOMEN
		IMPUTATION MODEL IS CORRECT
		ANALYSIS MODEL DOES NOT ASSUME THE ASSOCIATION OF WEIGHT AND LOGINSINDEX IS THE SAME IN MEN AND WOMEN

	4 - PROBABILITY DISTRIBUTIONS FOR ALL VARIABLES ARE THE SAME FOR MEN AND WOMEN
		THE VARIABLES weight AND loginsindex HAVE MODERATE NON-NORMAL DISTRIBUTIONS
		MEN AND WOMEN HAVE MISSING DATA
		IMPUTATION MODEL IS INCORRECT
		ANALYSIS MODEL IS INCORRECT
		
	5 - PROBABILITY DISTRIBUTIONS FOR ALL VARIABLES ARE THE SAME FOR MEN AND WOMEN
		THE VARIABLES weight AND loginsindex HAVE SEVERE NON-NORMAL DISTRIBUTIONS
		MEN AND WOMEN HAVE MISSING DATA
		IMPUTATION MODEL IS INCORRECT
		ANALYSIS MODEL IS INCORRECT
		
	6 - PROBABILITY DISTRIBUTION OF WEIGHT AND LOGINSINDEX HAVE DIFFERENT VARIANCES FOR MEN AND WOMEN
		ALL CONTINUOUS VARIABLES ARE NORMALLY DISTRIBUTED
		MEN AND WOMEN HAVE MISSING DATA
		IMPUTATION MODEL IS INCORRECT
		ANALYSIS MODEL IS INCORRECT
******************************************************************************************************************************************************************************************/
version 12.1

args scenario numberOfDatasets percentObserved B M v

* OUTCOMES OF THE IMPUTATION AND ANALYSIS MODELS, RESPECTIVELY
local Ioutcome = `"weight"'
local Aoutcome = `"loginsindex"'

if `scenario'== 0 {
	local pathway "Data\subgroup\"

	* COVARIATES OF THE IMPUTATION MODEL
	local Icovars `"age height loginsindex"'

	* COVARIATES OF THE ANALYSIS MODEL
	local Acovars `"constant age weight"'
	local A_b `"(b_rep_imp[1,1]) (b_rep_imp[1,2]) (b_rep_imp[1,3])"'
	local b_matrixIndex `"[1,1] [1,2] [1,3]"'
}
else if `scenario'== 1 {
	local pathway "Data\subgroup\"
	
	* COVARIATES OF THE IMPUTATION MODEL
	local Icovars `"age height loginsindex"'

	* COVARIATES OF THE ANALYSIS MODEL
	local Acovars `"constant age weight"'
	local A_b `"(b_rep_imp[1,1]) (b_rep_imp[1,2]) (b_rep_imp[1,3])"'
	local b_matrixIndex `"[1,1] [1,2] [1,3]"'
}
else if `scenario'== 2 {
	local pathway "Data\hetero_quadruple\"

	* COVARIATES OF THE IMPUTATION MODEL
	local Icovars `"sex age height loginsindex"'

	* COVARIATES OF THE ANALYSIS MODEL
	local Aoutcome = `"loginsindex"'
	local Acovars `"constant sex age weight"'
	local A_b `"(b_rep_imp[1,1]) (b_rep_imp[1,2]) (b_rep_imp[1,3]) (b_rep_imp[1,4])"'
	local b_matrixIndex `"[1,1] [1,2] [1,3] [1,4]"'
}
else if `scenario'== 3 {
	local pathway "Data\omitted\"
	
	* COVARIATES OF THE IMPUTATION MODEL
	local Icovars `"sex age height loginsindex"'

	* COVARIATES OF THE ANALYSIS MODEL
	local Aoutcome = `"loginsindex"'
	local Acovars `"constant sex age weight weightXsex"'
	local A_b `"(b_rep_imp[1,1]) (b_rep_imp[1,2]) (b_rep_imp[1,3]) (b_rep_imp[1,4]) (b_rep_imp[1,5])"'
	local b_matrixIndex `"[1,1] [1,2] [1,3] [1,4] [1,5]"'
}
else if `scenario'==4 {
	local pathway "Data\moderate\" 
	
	* COVARIATES OF THE IMPUTATION MODEL
	local Icovars `"age height loginsindex"'

	* COVARIATES OF THE ANALYSIS MODEL
	local Aoutcome = `"loginsindex"'
	local Acovars `"constant age weight"'
	local A_b `"(b_rep_imp[1,1]) (b_rep_imp[1,2]) (b_rep_imp[1,3])"'
	local b_matrixIndex `"[1,1] [1,2] [1,3]"'
}
else if `scenario'==5 {
	local pathway "Data\severe\" 
	
	* COVARIATES OF THE IMPUTATION MODEL
	local Icovars `"age height loginsindex"'

	* COVARIATES OF THE ANALYSIS MODEL
	local Aoutcome = `"loginsindex"'
	local Acovars `"constant age weight"'
	local A_b `"(b_rep_imp[1,1]) (b_rep_imp[1,2]) (b_rep_imp[1,3])"'
	local b_matrixIndex `"[1,1] [1,2] [1,3]"'
}
else if `scenario'== 6 {
	local pathway "Data\hetero_quarter\"

	* COVARIATES OF THE IMPUTATION MODEL
	local Icovars `"sex age height loginsindex"'

	* COVARIATES OF THE ANALYSIS MODEL
	local Aoutcome = `"loginsindex"'
	local Acovars `"constant sex age weight"'
	local A_b `"(b_rep_imp[1,1]) (b_rep_imp[1,2]) (b_rep_imp[1,3]) (b_rep_imp[1,4])"'
	local b_matrixIndex `"[1,1] [1,2] [1,3] [1,4]"'
}
else {
	display "error"
}
local numberOfIcovars : list sizeof Icovars
local numberOfAcovars : list sizeof Acovars

display "Simulation study for scenario " `scenario'
display "The " `numberOfIcovars' " covariates of the imputation model are " "`Icovars'"
display "The " `numberOfAcovars' " covariates of the analysis model are " "`Acovars'"

* SET TEMPORARY NAMES, VARIABLES AND FILES	
tempname postfile_bootstrap postfile_dataset ci_normal ci_percentile ci_bc bstar mlparameters
tempvar xb
tempfile file_incomplete file_bootstrap file_dataset 

* SET THE LISTS
local b_list `"b_full b_cca b_boot"'
local var_list `"var_full var_cca var_boot msb_boot msw_boot dfb_boot dfw_boot"'
local df_list `"df_boot"'
local lower_list `"lowerT_full lowerT_cca lowerT_boot lowerN_boot"'
local upper_list `"upperT_full upperT_cca upperT_boot upperN_boot"'

* GENERATE THE variable lists AND post lists
* PART 1
foreach item of local b_list {
	foreach var of local Acovars {
		local tempList "`item'_`var'"
		local part1VarList : list part1VarList | tempList
	}
	
	foreach index of local b_matrixIndex  {
		local tempList `"(`item'`index')"'
		local part1PostList : list part1PostList | tempList
	}
}
foreach item of local var_list {
	foreach var of local Acovars {
		local tempList "`item'_`var'"
		local part1VarList : list part1VarList | tempList
	}
	
	foreach index of local b_matrixIndex  {
		local tempList `"(`item'`index')"'
		local part1PostList : list part1PostList | tempList
	}
}
di "`part1PostList'"
* PART 2
foreach item of local df_list {
	foreach var of local Acovars {
		local tempList "`item'_`var'"
		local part2VarList : list part2VarList | tempList
	}
	
	foreach index of local b_matrixIndex  {
		local tempList `"(`item'`index')"'
		local part2PostList : list part2PostList | tempList
	}
}
foreach item of local lower_list {
	foreach var of local Acovars {
		local tempList "`item'_`var'"
		local part2VarList : list part2VarList | tempList
	}
	
	foreach index of local b_matrixIndex  {
		local tempList `"(`item'`index')"'
		local part2PostList : list part2PostList | tempList
	}
}
foreach item of local upper_list {
	foreach var of local Acovars {
		local tempList "`item'_`var'"
		local part2VarList : list part2VarList | tempList
	}
	
	foreach index of local b_matrixIndex  {
		local tempList `"(`item'`index')"'
		local part2PostList : list part2PostList | tempList
	}
}
di "`part2PostList'"

* DETERMINE SAMPLE SIZE OF A SIMULATED DATASET
quietly use "`pathway'dataset_1", clear
quietly count
local sampleSize = r(N)

matrix b_boot = J(1,`numberOfAcovars',.)
matrix var_boot = J(1,`numberOfAcovars',.)
matrix df_boot = J(1,`numberOfAcovars',.)
matrix lowerT_boot = J(1,`numberOfAcovars',.)
matrix upperT_boot = J(1,`numberOfAcovars',.)
matrix lowerT_full = J(1,`numberOfAcovars',.)
matrix upperT_full = J(1,`numberOfAcovars',.)
matrix lowerT_cca = J(1,`numberOfAcovars',.)
matrix upperT_cca = J(1,`numberOfAcovars',.)
matrix msb_boot = J(1,`numberOfAcovars',.)
matrix msw_boot = J(1,`numberOfAcovars',.)
matrix dfb_boot = J(1,`numberOfAcovars',.)
matrix dfw_boot = J(1,`numberOfAcovars',.)
matrix lowerN_boot = J(1,`numberOfAcovars',.)
matrix upperN_boot = J(1,`numberOfAcovars',.)

* STORE THE SEED NUMBERS FOR GENERATING AND PROCESSING EACH DATASET - USEFUL WHEN AN ERROR OCCURS
* AT THE END OF THE LOG FILE THE RESULTS ARE PRINTED
capture log close
log using "Log files\vonHippel_BMLMI_scenario`scenario'_`sampleSize'obs_observed`percentObserved'_M`M'_B`B'_V`v'", replace

postfile `postfile_dataset' dataset df_full df_cca `part1VarList' `part2VarList' numberofreplicates runtime using `file_dataset', replace
forvalues dataset=1(1)`numberOfDatasets' {

	display "vonHippel_BMLMI: scenario `scenario'; `percentObserved'% observed; B=`B'; M=`M': Processing dataset `dataset'"  
	display "current seed is " c(seed)
	
quietly {

	* OPEN DATASET
	use "`pathway'dataset_`dataset'", clear
	capture drop _*
	
	* COLLECT THE FULL DATA ESTIMATES	
		* COLLECT THE DEGREES OF FREEDOM FOR THE COMPLETED DATASET
	if `scenario'<2 regress `Aoutcome' `Acovars' if sex==0, noconstant
	else regress `Aoutcome' `Acovars', noconstant
	matrix b_full = e(b)
	matrix var_full = vecdiag(e(V))
	local nu_com = e(df_r)
	local studentT = invttail(`nu_com', 0.025)
	local column 0
	foreach covariate of local Acovars {
		local column = `column'+1
		matrix lowerT_full[1,`column'] = b_full[1,`column'] - sqrt(var_full[1,`column'])*`studentT' 
		matrix upperT_full[1,`column'] = b_full[1,`column'] + sqrt(var_full[1,`column'])*`studentT' 
	}

	* SET THE WEIGHT MEASUREMENTS TO BE MISSING
	replace weight =. if R==0
	replace weightXsex =. if R==0
	
	* COLLECT THE CCA DATA ESTIMATES	
	if `scenario'<2 regress `Aoutcome' `Acovars' if sex==0, noconstant
	else regress `Aoutcome' `Acovars', noconstant
	matrix b_cca = e(b)
	matrix var_cca = vecdiag(e(V))
	local nu_cca = e(df_r)
	local studentT = invttail(`nu_cca', 0.025)
	local column 0
	foreach covariate of local Acovars {
		local column = `column'+1
		matrix lowerT_cca[1,`column'] = b_cca[1,`column'] - sqrt(var_cca[1,`column'])*`studentT' 
		matrix upperT_cca[1,`column'] = b_cca[1,`column'] + sqrt(var_cca[1,`column'])*`studentT' 
	}
	
	save `file_incomplete', replace
	
	timer clear 1
	timer on 1
	
	* BOOTSTRAP THEN IMPUTE  
	postfile `postfile_bootstrap' replicate imputation `Acovars' using `file_bootstrap', replace
	forvalues replicate=1(1)`B' {
		
		use `file_incomplete', clear
		
		* STEP 1: SAMPLE WITH REPLACEMENT THE INCOMPLETE DATASET
		quietly bsample
				
		* STEP 2: OBTAIN THE ML ESTIMATES OF THE IMPUTATION MODEL
		regressML weight `Icovars'
		matrix `mlparameters' = r(parameters)
		local num_betas = colsof(`mlparameters') - 1
		matrix `bstar' = `mlparameters'[1,1..`num_betas']
		local rmsestar = `mlparameters'[1,`num_betas'+1]
	
		* STEP 3: IMPUTE M TIMES
		forvalues imputation=1(1)`M' {
		
			* IMPUTE BY SAMPLING CONDITIONAL DISTRIBUTION
			capture drop `xb'
			matrix score `xb' = `bstar'
			quietly replace weight = `xb' + `rmsestar'*rnormal() if R==0
			
			* PASSIVELY IMPUTE THE VARIABLES DERIVED FROM WEIGHT
			quietly replace weightXsex = weight*sex
			
			* APPLY THE ANALYSIS MODEL  
			if `scenario'<2 regress `Aoutcome' `Acovars' if sex==0, noconstant
			else regress `Aoutcome' `Acovars', noconstant
			matrix b_rep_imp = e(b)
	
			post `postfile_bootstrap' (`replicate') (`imputation') `A_b'
		} // END OF imputation FOR-LOOP
	} // END OF replicate FOR-LOOP
	postclose `postfile_bootstrap'

	* STEP 3: CALCULATE BOOTSTRAP POINT ESTIMATES AND VARIANCE ESTIMATES 
	use `file_bootstrap', clear
	
	quietly count
	local numberofreplicates = r(N)
	
	local column 0
	foreach var of local Acovars {
		local column = `column' + 1
				
		* POINT ESTIMATE
		summarize `var'
		matrix b_boot[1,`column'] = r(mean)
		
		if `M'>1 {
			* ONE-WAY ANALYSIS OF VARIANCE
			anova `var' replicate
			local MSB = e(mss)/e(df_m)
			local MSW = e(rss)/e(df_r)
			matrix msb_boot[1,`column'] = `MSB'
			matrix dfb_boot[1,`column'] = e(df_m)
			matrix msw_boot[1,`column'] = `MSW'
			matrix dfw_boot[1,`column'] = e(df_r)
			
			* VARIANCE ESTIMATE
			local E = (`B'+1)/(`B'*`M')
			matrix var_boot[1,`column'] = `E'*`MSB' - `MSW'/`M'	
			
			* DEGREES OF FREEDOM	
			local F = (`E'*`MSB' - `MSW'/`M')^2
			local G = 1/(`B'*(`M'^2)*(`M'-1))
			local H = (1/(`B-1'))*(`E'^2)*(`MSB'^2) + `G'*`MSW'^2
			local nu_boot = `F'/`H'
		}	
		else {
			* VARIANCE ESTIMATE
			summarize `var'
			
			* DEGREES OF FREEDOM	
			matrix var_boot[1,`column'] = r(Var)
			local nu_boot = `B'-1
		}	
		matrix df_boot[1,`column'] = `nu_boot'
			
		* LIMITS OF THE CONFIDENCE INTERVAL
		local studentT = invttail(`nu_boot', 0.025)
		matrix lowerT_boot[1,`column'] = b_boot[1,`column'] - sqrt(var_boot[1,`column'])*`studentT' 
		matrix upperT_boot[1,`column'] = b_boot[1,`column'] + sqrt(var_boot[1,`column'])*`studentT'	
	
		* NORMAL-BASED CONFIDENCE INTERVAL
		matrix lowerN_boot[1,`column'] = b_boot[1,`column'] - sqrt(var_boot[1,`column'])*1.96 
		matrix upperN_boot[1,`column'] = b_boot[1,`column'] + sqrt(var_boot[1,`column'])*1.96
	}
	
	timer off 1
	quietly timer list
    	
	post `postfile_dataset' (`dataset') (`nu_com') (`nu_cca') `part1PostList' `part2PostList' (`numberofreplicates') (r(t1))
	} // END OF QUIETLY	
} // END OF dataset FOR-LOOP	
postclose `postfile_dataset'
log close
	
use `file_dataset', clear
if c(version) >=13 saveold "Results\vonHippel_BMLMI_scenario`scenario'_`sampleSize'obs_observed`percentObserved'_M`M'_B`B'_V`v'", version(12) replace
else save "Results\vonHippel_BMLMI_scenario`scenario'_`sampleSize'obs_observed`percentObserved'_M`M'_B`B'_V`v'", replace

* END OF DO FILE
