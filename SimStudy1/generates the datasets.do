/***********************************************************************************************************************************************************************************************
DO FILE GENERATES THE DATASETS FOR ALL THREE SCENARIOS  

USES THE BARRY CAERPHILLY FOLLOW-UP STUDY ON THE 951 YOUNG ADULTS

ARGS
	- size		: NUMBER OF OBSERVATIONS IN THE DATASET
	- mu1		: SEX VECTOR OF MULTIVARIATE NORMAL DISTRIBUTION age,height | sex	
	- kappa1	: SEX COEFFICIENT OF weight|constant,age,sex,height 
	- omegaW	: VARIANCE TERM FOR WOMEN IN REGRESSION MODEL loginsindex|constant,age,sex,weight
	- theta1 	: SEX COEFFICIENT OF loginsindex|cons,age,sex,weight
	- lambdaW	: VARIANCE TERM FOR WOMEN IN REGRESSION MODEL weight|cons,age,sex,height 
	- notNormal	: SPECIFIES WHICH NON-NORMAL ERROR DISTRIBUTIONS ARE REQUIRED FOR weight AND loginsindex 
					- IF NOT AN INTEGER NUMBER 1 TO 17 THEN ASSUMES THE DISTRIBUTIONS ARE NORMAL
	
DATE CREATED: 1st JULY 2018
UPDATED: 4TH AUGUST 2018
********************************************************************************************************************************************************************************************/
args scenario sampleSize probabilityObserved numberOfDatasets

tempname mu1 kappa1 theta1 lambdaW omegaW

quietly do "Do files\deriving data model parameters.do"

* error_weight ~ N(0,1) AND error_loginsinsex ~ N(0,1) UNLESS OTHERWISE SPECIFIED BELOW
local gen_ew "rnormal()"
local gen_el "rnormal()"

* CALCULATE THE PARAMETERS FOR THE SPECIFIC SCENARIO
if `scenario'== 1 {
	local pathway "Data\subgroup\"
	matrix `mu1' = (0\0)
	scalar `kappa1' = 0
	scalar `omegaW' = omegaM
	scalar `theta1' = 0
	scalar `lambdaW' = lambdaM		
}
else if `scenario'== 2 {
	local pathway "Data\hetero_quadruple\"
	matrix `mu1' = mu1
	scalar `kappa1' = kappa1
	scalar `omegaW' = omegaM*4
	scalar `theta1' = theta1
	scalar `lambdaW' = lambdaM*4
}
else if `scenario'== 3 {
	local pathway "Data\omitted\"
	matrix `mu1' = mu1
	scalar `kappa1' = kappa1
	scalar `omegaW' = omegaM
	scalar `theta1' = theta1
	scalar `lambdaW' = lambdaM
}
else if `scenario'==4 {
	* error_weight ~ exp(N(0,1/4^2)) AND error_loginsinsex ~ exp(N(0,1/4^2))
	local pathway "Data\moderate\" 
	local gen_ew "exp(0.25*rnormal())"
	local gen_el "exp(0.25*rnormal())"
	matrix `mu1' = (0\0)
	scalar `kappa1' = 0
	scalar `omegaW' = omegaM
	scalar `theta1' = 0
	scalar `lambdaW' = lambdaM
}
else if `scenario'==5 {
	* error_weight ~ student(df=3) AND error_loginsinsex ~ exp(N(0,1))
	local pathway "Data\severe\" 
	local gen_ew "rt(3)"
	local gen_el "exp(rnormal())"
	matrix `mu1' = (0\0)
	scalar `kappa1' = 0
	scalar `omegaW' = omegaM
	scalar `theta1' = 0
	scalar `lambdaW' = lambdaM
} 
else if `scenario'==6 {
	local pathway "Data\hetero_quarter\" 
	matrix `mu1' = mu1
	scalar `kappa1' = kappa1
	scalar `omegaW' = omegaM*0.25
	scalar `theta1' = theta1
	scalar `lambdaW' = lambdaM*0.25
}
else if `scenario'==8 {		// SAME AS SCENARIO 1 EXCEPT MISSINGNESS WILL BE IN MEN AND WOMEN (SEE BELOW)
	local pathway "Data\perfect\"
	matrix `mu1' = (0\0)
	scalar `kappa1' = 0
	scalar `omegaW' = omegaM
	scalar `theta1' = 0
	scalar `lambdaW' = lambdaM		
}

forvalues dataset=1(1)`numberOfDatasets' {
	clear
	set obs `sampleSize'

	* GENERATE ERROR VARIABLES FOR weight AND loginsindex
	gen error_weight = `gen_ew'
	gen error_loginsindex = `gen_el'

	* CENTRE THE ERROR VARIABLES, WITH UNIT VARIANCE
	summarize error_weight
	replace error_weight =  (error_weight - r(mean))/r(sd)
	summarize error_loginsindex
	replace error_loginsindex =  (error_loginsindex - r(mean))/r(sd)
	
	* DRAW SEX FROM A BINOMIAL DISTRIBUTION
	gen sex = rbinomial(1, pi)

	* DRAW FROM THE MULTIVARIATE NORMAL DISTRIBUTION 
	* age   | sex  ~ N(mu0 + mu1*sex, Sigma)
	* height
	drawnorm age height, mean(mu0) cov(Sigma) 
	tempname mean
	matrix `mean' = mu0 + `mu1'
	drawnorm age1 height1, mean(`mean') cov(Sigma) 
	replace age = age1 if sex==1
	replace height = height1 if sex==1
	drop age1 height1

	gen weight = kappa0 + kappa2*age + kappa3*height + sqrt(omegaM)*error_weight if sex==0
	replace weight = kappa0 + `kappa1' + kappa2*age + kappa3*height + sqrt(`omegaW')*error_weight if sex==1

	gen loginsindex = theta0 + theta2*age + theta3*weight + sqrt(lambdaM)*error_loginsindex if sex==0
	replace loginsindex = theta0 + `theta1' + theta2*age + theta3*weight + sqrt(`lambdaW')*error_loginsindex if sex==1

	gen constant = 1

	* GENERATES SUBGROUP VARIABLES
	gen consMen = constant
	replace consMen = 0 if sex==1
	gen consWomen = constant
	replace consWomen = 0 if sex==0

	gen weightMen = weight if sex==0
	replace weightMen = 0 if sex==1
	gen weightWomen = weight if sex==1
	replace weightWomen = 0 if sex==0

	gen ageMen = age if sex==0
	replace ageMen = 0 if sex==1
	gen ageWomen = age if sex==1
	replace ageWomen = 0 if sex==0

	* GENERATE INTERACTION TERM
	gen weightXsex = weight*sex

	* SIMULATE MISSING DATA - MCAR
		* - FOR SCENARIOS 0 AND 1 ONLY MEN HAVE MISSING WEIGHT MEASUREMENTS
		* - OTHERWISE MEN AND WOMEN HAVE MISSING WEIGHT MEASUREMENTS
	if `scenario'==0 | `scenario'==1 {
		gen R = runiform()< `probabilityObserved' if sex==0
	}
	else {
		gen R = runiform()< `probabilityObserved'
	}

	save "`pathway'dataset_`dataset'", replace
} // END OF dataset FOR-LOOP

* END OF DO FILE
