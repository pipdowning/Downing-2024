# Downing-2024
Code and data for: Downing PA. 2024. Michener’s group-size paradox in cooperatively breeding birds. Journal of Evolutionary Biology, 37: 353-359.

Number of supplementary items: seven
1. MP_R_Code.R
2. MP_Table_S1.csv
3. MP_Table_S2.csv
4. MP_Table_S3.csv
5. MP_Data_Extraction.txt
6. Cooperative_bird_literature.pdf
7. MP_Supp_Methods_Figures.pdf


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

File name: MP_R_Code.R

This R script contains all the code needed to replicate the analyses (including packages and functions).
Models cited in the manuscript are located in this code.
Lettering in the code matches the lettering used in the methods section of the manuscript.

- Data manipulation (lines 25 to 100)
- Part A. per capita and total reproductive success across group sizes (lines 106 to 256)
- Part B. per capita and total reproductive success and the frequency of group sizes (lines 262 to 434)
- Part C. study effort (lines 440 to 627)

Note that the code includes two-stage models fit using restricted maximum likelihood in the metafor R package to investigate the sensitivity of the results to Bayesian methods (results reported in Table S3).
See Viechtbauer, W. (2010) https://www.metafor-project.org/doku.php/tips:two_stage_analysis


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

File name: MP_Table_S1.csv

This csv document contains:\
	+ reproductive success estimates (annual number of fledglings) in different sizes groups for 26 species\
	+ the data are organised in long format (137 rows)\
	+ read into R (some of the column headings will need changing to match the R code)\
	+ column descriptions:\
		A. common name = English name of each species\
		B. latin binomial = scientific name of each species (matches the Jetz et al. nomenclature)\
		C. measure = the metric used to measure annual reproductive success in each study\
		D. mean = estimate of mean reproductive success\
		E. standard_dev = standard deviation of the mean reproductive success estimate\
		F. variance = variance of the mean reproductive success estimate\
		G. N = sample size used to estimate mean reproductive success\
		H. group_size = the group size pertaining to each mean reproductive success and % of groups estimate\
		I. source = where the mean reproductive success estimates were extracted from in the study\
		J. reference = study from which the mean reproductive success estimates were extracted\
		K. % of groups = the frequency of each group size in nature\
		L. N = sample size used to estimate % of groups \
		M. source = where the % of groups estimates were extracted from in the study\
		N. reference = study from which % of groups estimates were extracted\


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

File name: MP_Table_S2.csv

This csv document contains:\
	+ the frequency of different group sizes in natural populations for 23 species\
	+ data are absent for three species in Table S1: chestnut-crowned babbler, toucan barbet, Pygmy Nuthatch\
	+ the data in columns D and E are entered in columns L and K in Table S1, but only for group size with reproductive success estimates\
	+ the data are organised in long format (183 rows)\
	+ read into R (some of the column headings will need changing to match the R code)\
	+ column descriptions:\
		A. common name = English name of each species\
		B. latin binomial = scientific name of each species (matches the Jetz et al. nomenclature)\
		C. group_size = the group size pertaining to each % of groups estimate\
		D. N = sample size used to estimate % of groups\
		E. % of groups = the frequency of each group size in nature\
		F. source = where the % of groups estimates were extracted from in the study\
		G. reference = study from which % of groups estimates were extracted


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

File name: MP_Table_S3.csv

This csv document contains:\
	+ parameter estimates from statistical models (see MP_R_Code.R for further information)\
	+ 34 rows = output from 8 statistical models\
	+ includes results from the two-stage models (see MP_Supp_Methods_Figures.pdf)\
	+ variable names:\
		beta = parameter estimate (slope) from the model\
		lwr CI = lower 95% credible interval from the posterior distribution of the model\
		upr CI = upper 95% credible interval from the posterior distribution of the model


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

File name: MP_Data_Extraction.txt

This plain text document contains details of the figures, tables, and text fragments from which data were obtained, and the calculations used to pool data, organised alphabetically by latin name.


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

File name: Cooperative_bird_literature.pdf

A database of 879 studies on 182 species of cooperatively breeding birds, organised alphabetically by latin name.
Comprises all known studies (published research, MSc and PhD theses, monographs and edited volumes) with data on the breeding biology and fitness parameters of cooperative birds.


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

File name: MP_Supp_Methods_Figures.pdf

This PDF document contains further details on the methods used in the study (prior specification and model convergence, quadratic effects and study effort) and supplementary figures S1 to S5.


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
