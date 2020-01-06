/***********************************************************************************************************************************************************************
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
		THE VARIABLES weight AND loginsindex HAVE NON-NORMAL DISTRIBUTIONS
		MEN AND WOMEN HAVE MISSING DATA
		IMPUTATION MODEL IS INCORRECT
		ANALYSIS MODEL IS INCORRECT
****************************************************************************************************************************************************/
* REQUIRES ROW MATRICES AS INPUT
capture mata: mata drop mim4()
mata:
	void mim4(string matrix estimates, string matrix variances, real scalar nu_com)
{
	real matrix all, all_var, A, A_var, withinVar, betweenVar, coefficients, totalVar, lowerCI, upperCI, df
	real scalar m, cols, nu, gamma_m, nu_obs, nu_star, studentT, column

	A = st_matrix(estimates)
	A_var = st_matrix(variances)
	m = rows(A)
	cols = cols(A)
		
	// INITIALIZE THE MATRICES THAT WILL HOLD THE RESULTS
	coefficients = J(1,cols,0)
	totalVar = J(1,cols,0)
	withinVar = J(1,cols,.)
	betweenVar = J(1,cols,.)
	df = J(1,cols,0)
	lowerCI = J(1,cols,0)
	upperCI = J(1,cols,0)
	
	// MEAN OF THE m ESTIMATES
	coefficients[1,.] = mean(A)

	// CALCULATION OF THE BETWEEN, WITHIN VARIANCE AND TOTAL VARIANCE, DF AND CONFIDENCE LIMITS
	withinVar[1,.] = mean(A_var)
	
	if (m==1) {
		totalVar = withinVar
		studentT = invttail(nu_com, 0.025)
		lowerCI = coefficients - sqrt(totalVar)*studentT
		upperCI = coefficients + sqrt(totalVar)*studentT
	}
	else {
		betweenVar[1,.] = diagonal(variance(A))'

		for(column=1; column<=cols; column++) {
			totalVar[1,column] = withinVar[1,column] + (1+(1/m))*betweenVar[1,column]
				
			// CALCULATION OF DEGREES OF FREEDOM WITH REFINEMENT FOR SMALL SAMPLES (BARNARD AND RUBIN 1999)
			// SEE EQUATIONS 10.15, AND 10.16, PAGE 211, LITTLE AND RUBIN 2002 
			nu = (m-1)*(1 + (m/(m+1))*(withinVar[1,column]/betweenVar[1,column]))^2 
			gamma_m = ((1 + (1/m))*betweenVar[1,column])/(withinVar[1,column]+(1+(1/m))*betweenVar[1,column])
			nu_obs = (1-gamma_m)*((nu_com+1)/(nu_com+3))*nu_com 
			nu_star = 1/((1/nu) + (1/nu_obs))
			df[1,column] = nu_star
			
			// DETERMINE THE VALUE OF THE 2.5 PERCENTILE POINT OF THE T-DISTRIBUTION WITH nu_star DEGREES OF FREEDOM
			studentT = invttail(nu_star, 0.025)
			lowerCI[1,column] = coefficients[1,column] - sqrt(totalVar[1,column])*studentT
			upperCI[1,column] = coefficients[1,column] + sqrt(totalVar[1,column])*studentT
		}
	}
	
	// RETURNS THE CONTENTS OF MATRIX results SO THAT IT IS ACCESSIBLE OUTSIDE MATA AS r(imputes)
	st_rclear()
	st_matrix("r(coefficients)", coefficients)
	st_matrix("r(betweenVar)", betweenVar)
	st_matrix("r(withinVar)", withinVar)
	st_matrix("r(totalVar)", totalVar)
	st_matrix("r(df)", df)
	st_matrix("r(lowerCI)", lowerCI)
	st_matrix("r(upperCI)", upperCI)
}
end
********************************************************************************************************************************************************************************************************
version 12.1

args scenario numberOfDatasets percentObserved B M v

if `scenario'== 0 {
	local pathway "Data\subgroup\"

	* NAME VARIABLES FOR THE IMPUTATION MODEL
	local Ioutcome = `"weight"'
	local Icovars `"age height loginsindex"'

	* NAME THE VARIABLE FOR THE ANALYSIS MODEL
	local Aoutcome = `"loginsindex"'
	local Acovars `"constant age weight"'
	local A_b `"(b_imp_rep[1,1]) (b_imp_rep[1,2]) (b_imp_rep[1,3])"'
	local A_wv `"(wv[1,1]) (wv[1,2]) (wv[1,3])"'
	local b_matrixIndex `"[1,1] [1,2] [1,3]"'
}
else if `scenario'== 1 {
	local pathway "Data\subgroup\"
	
	* NAME VARIABLES FOR THE IMPUTATION MODEL
	local Ioutcome = `"weight"'
	local Icovars `"age height loginsindex"'

	* NAME THE VARIABLE FOR THE ANALYSIS MODEL
	local Aoutcome = `"loginsindex"'
	local Acovars `"constant age weight"'
	local A_b `"(b_imp_rep[1,1]) (b_imp_rep[1,2]) (b_imp_rep[1,3])"'
	local A_wv `"(wv[1,1]) (wv[1,2]) (wv[1,3])"'
	local b_matrixIndex `"[1,1] [1,2] [1,3]"'
}
else if `scenario'== 2 {
	local pathway "Data\hetero_quadruple\"

	* NAME VARIABLES FOR THE IMPUTATION MODEL
	local Ioutcome = `"weight"'
	local Icovars `"sex age height loginsindex"'

	* NAME THE VARIABLE FOR THE ANALYSIS MODEL
	local Aoutcome = `"loginsindex"'
	local Acovars `"constant sex age weight"'
	local A_b `"(b_imp_rep[1,1]) (b_imp_rep[1,2]) (b_imp_rep[1,3]) (b_imp_rep[1,4])"'
	local A_wv `"(wv[1,1]) (wv[1,2]) (wv[1,3]) (wv[1,4])"'
	local b_matrixIndex `"[1,1] [1,2] [1,3] [1,4]"'
}
else if `scenario'== 3 {
	local pathway "Data\omitted\"
	
	* NAME VARIABLES FOR THE IMPUTATION MODEL
	local Ioutcome = `"weight"'
	local Icovars `"sex age height loginsindex"'

	* NAME THE VARIABLE FOR THE ANALYSIS MODEL
	local Aoutcome = `"loginsindex"'
	local Acovars `"constant sex age weight weightXsex"'
	local A_b `"(b_imp_rep[1,1]) (b_imp_rep[1,2]) (b_imp_rep[1,3]) (b_imp_rep[1,4]) (b_imp_rep[1,5])"'
	local A_wv `"(wv[1,1]) (wv[1,2]) (wv[1,3]) (wv[1,4]) (wv[1,5])"'
	local b_matrixIndex `"[1,1] [1,2] [1,3] [1,4] [1,5]"'
}
else if `scenario'==4 {
	local pathway "Data\moderate\" 
	
	* NAME VARIABLES FOR THE IMPUTATION MODEL
	local Ioutcome = `"weight"'
	local Icovars `"age height loginsindex"'

	* NAME THE VARIABLE FOR THE ANALYSIS MODEL
	local Aoutcome = `"loginsindex"'
	local Acovars `"constant age weight"'
	local A_b `"(b_imp_rep[1,1]) (b_imp_rep[1,2]) (b_imp_rep[1,3])"'
	local A_wv `"(wv[1,1]) (wv[1,2]) (wv[1,3])"'
	local b_matrixIndex `"[1,1] [1,2] [1,3]"'
}
else if `scenario'==5 {
	local pathway "Data\severe\" 
	
	* NAME VARIABLES FOR THE IMPUTATION MODEL
	local Ioutcome = `"weight"'
	local Icovars `"age height loginsindex"'

	* NAME THE VARIABLE FOR THE ANALYSIS MODEL
	local Aoutcome = `"loginsindex"'
	local Acovars `"constant age weight"'
	local A_b `"(b_imp_rep[1,1]) (b_imp_rep[1,2]) (b_imp_rep[1,3])"'
	local A_wv `"(wv[1,1]) (wv[1,2]) (wv[1,3])"'
	local b_matrixIndex `"[1,1] [1,2] [1,3]"'
}
else if `scenario'== 6 {
	local pathway "Data\hetero_quarter\"

	* NAME VARIABLES FOR THE IMPUTATION MODEL
	local Ioutcome = `"weight"'
	local Icovars `"sex age height loginsindex"'

	* NAME THE VARIABLE FOR THE ANALYSIS MODEL
	local Aoutcome = `"loginsindex"'
	local Acovars `"constant sex age weight"'
	local A_b `"(b_imp_rep[1,1]) (b_imp_rep[1,2]) (b_imp_rep[1,3]) (b_imp_rep[1,4])"'
	local A_wv `"(wv[1,1]) (wv[1,2]) (wv[1,3]) (wv[1,4])"'
	local b_matrixIndex `"[1,1] [1,2] [1,3] [1,4]"'
}
if `scenario'== 7 {
	local pathway "Data\subgroup\"

	* NAME VARIABLES FOR THE IMPUTATION MODEL
	local Ioutcome = `"weight"'
	local Icovars `"age height loginsindex"'

	* NAME THE VARIABLE FOR THE ANALYSIS MODEL
	local Aoutcome = `"loginsindex"'
	local Acovars `"constant age weight"'
	local A_b `"(b_imp_rep[1,1]) (b_imp_rep[1,2]) (b_imp_rep[1,3])"'
	local A_wv `"(wv[1,1]) (wv[1,2]) (wv[1,3])"'
	local b_matrixIndex `"[1,1] [1,2] [1,3]"'
}
if `scenario'== 8 {
	local pathway "Data\perfect\"

	* NAME VARIABLES FOR THE IMPUTATION MODEL
	local Ioutcome = `"weight"'
	local Icovars `"age height loginsindex"'

	* NAME THE VARIABLE FOR THE ANALYSIS MODEL
	local Aoutcome = `"loginsindex"'
	local Acovars `"constant age weight"'
	local A_b `"(b_imp_rep[1,1]) (b_imp_rep[1,2]) (b_imp_rep[1,3])"'
	local A_wv `"(wv[1,1]) (wv[1,2]) (wv[1,3])"'
	local b_matrixIndex `"[1,1] [1,2] [1,3]"'
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
tempname postfile_bootstrap postfile_dataset postfile_withinvarboot
tempfile file_MI file_bootstrap file_dataset file_withinvarboot

* SET THE LISTS
local b_list `"b_full b_cca b_MI"' 
local var_list `"var_full var_cca var_MI varB_MI varW_MI var_boot varW_boot"'
local df_list `"df_MI dfRS_boot dfBR_boot"'
local lower_list `"lowerT_full lowerT_cca lowerT_BR_MI lowerT_RS_boot lowerT_BR_boot"'
local upper_list `"upperT_full upperT_cca upperT_BR_MI upperT_RS_boot upperT_BR_boot"'

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

matrix wv = J(1,`numberOfAcovars',.)
matrix lowerT_BR_MI = J(1,`numberOfAcovars',.)
matrix upperT_BR_MI = J(1,`numberOfAcovars',.)
matrix varW_boot = J(1,`numberOfAcovars',.)
matrix var_boot = J(1,`numberOfAcovars',.)
matrix dfRS_boot = J(1,`numberOfAcovars',.)
matrix dfBR_boot = J(1,`numberOfAcovars',.)
matrix lowerT_RS_boot = J(1,`numberOfAcovars',.)
matrix upperT_RS_boot = J(1,`numberOfAcovars',.)
matrix lowerT_BR_boot = J(1,`numberOfAcovars',.)
matrix upperT_BR_boot = J(1,`numberOfAcovars',.)

* MATRICES TO STORE RUBIN'S IMPUTATION RESULTS
matrix b_m = J(`M',`numberOfAcovars',.)
matrix V_m = J(`M',`numberOfAcovars',.)

matrix lowerT_full = J(1,`numberOfAcovars',.)
matrix upperT_full = J(1,`numberOfAcovars',.)
matrix lowerT_cca = J(1,`numberOfAcovars',.)
matrix upperT_cca = J(1,`numberOfAcovars',.)

* STORE THE SEED NUMBERS FOR GENERATING AND PROCESSING EACH DATASET - USEFUL WHEN AN ERROR OCCURS
* AT THE END OF THE LOG FILE THE RESULTS ARE PRINTED
capture log close
log using "Log files\SandH_method2_scenario`scenario'_`sampleSize'obs_observed`percentObserved'_V`v'", replace

set seed 070818

postfile `postfile_dataset' dataset df_full df_cca `part1VarList' `part2VarList' runtime using `file_dataset', replace
forvalues dataset=1(1)`numberOfDatasets' {

	display "SandH_method2: scenario `scenario'; `percentObserved'% observed; B=`B'; M=`M': Processing dataset `dataset'"  
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

	replace weight =. if R==0
	replace weightMen =. if R==0
	replace weightWomen =. if R==0
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
	
	timer clear 1
	timer on 1
	* STEP 1: IMPUTE M TIMES
		* CHECKED USING dryrun
	ice weight age height loginsindex sex weightXsex weightMen weightWomen, m(`M') saving(`file_MI', replace) ///
		eq(weight:`Icovars') ///
		passive(weightXsex:weight*sex\weightMen:weight-weightXsex\weightWomen:weightXsex)
	
		* OBTAIN THE MI POINT ESTIMATE; ALSO RECORD OTHER FEATURES OF RUBIN'S RULES
	use `file_MI', clear
	forvalues imputation=1(1)`M' {
		if `scenario'<2 regress `Aoutcome' `Acovars' if sex==0 & _mj==`imputation', noconstant
		else regress `Aoutcome' `Acovars' if _mj==`imputation', noconstant
		matrix b_m[`imputation',1] = e(b)
		matrix V_m[`imputation',1] = vecdiag(e(V))
	}
	mata: mim4("b_m", "V_m", `nu_com')
	matrix b_MI = r(coefficients)
	matrix var_MI = r(totalVar)
	matrix varB_MI = r(betweenVar)
	matrix varW_MI = r(withinVar)	
	matrix df_MI = r(df)
	matrix lowerT_BR_MI = r(lowerCI)
	matrix upperT_BR_MI = r(upperCI)
	
	* STEP 2: BOOTSTRAP THE IMPUTED DATASETS INDEPENDENTLY
		* FOR EACH IMPUTED DATASET POST THE WITHIN-IMPUTATION VARIANCE
	postfile `postfile_withinvarboot' imputation `Acovars' using `file_withinvarboot', replace
	forvalues imputation=1(1)`M' {	
		postfile `postfile_bootstrap' replicate `Acovars' using `file_bootstrap', replace
		forvalues replicate=1(1)`B' {	
			use `file_MI', clear
			keep if _mj==`imputation'
		
			* BOOTSTRAP THE IMPUTED DATASET  
			quietly bsample
					
			* APPLY THE ANALYSIS MODEL AND POST THE ESTIMATES 
			if `scenario'<2 regress `Aoutcome' `Acovars' if sex==0, noconstant
			else regress `Aoutcome' `Acovars', noconstant
			matrix b_imp_rep = e(b)			
			post `postfile_bootstrap' (`replicate') `A_b'
			
		} // END OF replicate FOR-LOOP
		postclose `postfile_bootstrap'
		
		* CALCULATE THE M WITHIN-IMPUTATION VARIANCES
		use `file_bootstrap', clear 
		local column 0
		foreach covariate of local Acovars {
			local column =`column'+1
			
			summarize `covariate'
			matrix wv[1,`column'] = r(Var)
		}		
		post `postfile_withinvarboot' (`imputation') `A_wv'
		
	} // END OF imputation FOR-LOOP
	postclose `postfile_withinvarboot'

	* STEP 3: CALCULATE RUBIN'S TOTAL VARIANCE AND T-INTERVAL USING THE BOOTSTRAP ESTIMATES
		* WITHIN VARIANCE AND TOTAL VARIANCE
	use `file_withinvarboot', clear
	local column 0
	foreach covariate of local Acovars {
		local column = `column'+1
		
		* WITHIN-VARIANCE
		summarize `covariate'
		matrix varW_boot[1,`column'] = r(mean)
		local Wvar_boot = varW_boot[1,`column']
		
		if `M'==1 {
			* TOTAL VARIANCE
			matrix var_boot[1,`column'] = `Wvar_boot'
			matrix dfRS_boot[1,`column'] = `nu_com'
			matrix dfBR_boot[1,`column'] = `nu_com'
			local studentT_RS = invttail(`nu_com', 0.025)
			local studentT_BR = invttail(`nu_com', 0.025)
		}
		else {
			* BETWEEN-VARIANCE
			local Bvar_MI = varB_MI[1,`column']
			
			* TOTAL VARIANCE
			matrix var_boot[1,`column'] = `Wvar_boot' + (1+(1/`M'))*`Bvar_MI'
					
			* RUBIN AND SCHENKER 1986 DF FOR SMALL SAMPLE
			local nu_RS = (`M'-1)*(1 + (`M'/(`M'+1))*(`Wvar_boot'/`Bvar_MI'))^2
			matrix dfRS_boot[1,`column'] = `nu_RS'
			
			* STUDENT-T USING RS DEGREES OF FREEDOM
			local studentT_RS = invttail(`nu_RS', 0.025)
		
			* BARNARD AND RUBIN'S 1999 DEGREES OF FREEDOM 
			local gamma_m = ((1 + (1/`M'))*`Bvar_MI')/(`Wvar_boot'+(1+(1/`M'))*`Bvar_MI')
			local nu_obs = (1-`gamma_m')*((`nu_com'+1)/(`nu_com'+3))*`nu_com' 
			local nu_star = 1/((1/`nu_RS') + (1/`nu_obs'))	
			matrix dfBR_boot[1,`column'] = `nu_star'
			
			* STUDENT-T USING BR DEGREES OF FREEDOM
			local studentT_BR = invttail(`nu_star', 0.025)
		}
		
		* CONFIDENCE LIMITS
			* USING RUBIN AND SCHENKER DF
		matrix lowerT_RS_boot[1,`column'] = b_MI[1,`column'] - sqrt(var_boot[1,`column'])*`studentT_RS' 
		matrix upperT_RS_boot[1,`column'] = b_MI[1,`column'] + sqrt(var_boot[1,`column'])*`studentT_RS' 
			* USING BARNARD AND RUBIN DF
		matrix lowerT_BR_boot[1,`column'] = b_MI[1,`column'] - sqrt(var_boot[1,`column'])*`studentT_BR' 
		matrix upperT_BR_boot[1,`column'] = b_MI[1,`column'] + sqrt(var_boot[1,`column'])*`studentT_BR' 
	}		
	timer off 1
	quietly timer list
    	
	post `postfile_dataset' (`dataset') (`nu_com') (`nu_cca') `part1PostList' `part2PostList' (r(t1))
	} // END OF QUIETLY	
} // END OF dataset FOR-LOOP	
postclose `postfile_dataset'
log close
	
use `file_dataset', clear
if c(version) >=13 saveold "Results\SandH_method2_scenario`scenario'_`sampleSize'obs_observed`percentObserved'_M`M'_B`B'_V`v'", version(12) replace
else save "Results\SandH_method2_scenario`scenario'_`sampleSize'obs_observed`percentObserved'_M`M'_B`B'_V`v'", replace

* END OF DO FILE
