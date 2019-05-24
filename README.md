# Rcode-JASA2019
R code to reproduce all simulation results in "Estimating Optimal Dynamic Treatment Regimes with Survival Outcomes" in JASA (2019)

## Content
- plot_for_paper.R: to reproduce Figure 1 in the manuscript
- Supplemental_simulations.Rnw: reproduce tables and figures in Supplementary Material C
- Supplemental_asymptotic_variance.Rnw: reproduce tables and figures in Supplementary Material B
- Asymptotic_variance:
	- Asymptotic_variance.R: run this script to reproduce the results presented in
	  Section 4 of Supplementary Material B
	- 'Results': outputs from "Asymptotic_variance.R"
	- child_asymp_XXX.Rnw: separate .Rnw files with Tables 1-3 of 
	  Supplementary Material B
- Bias_precision:
	- bias_precision.R: run this script to reproduce the results presented in
	  Section 2 of Supplementary Material C
	- 'Independent30', 'Independent60', 'X1dependent30', 'X1dependent60', 
	  'X1X2dependent30' and 'X1X2dependent60': outputs from "bias_precision.R"
	- bias_precision_XXX.png: plots in Section 2 of Supplementary Material C
	- reproduce_plot_bias_precision.R: run this script to reproduce the plots
	  bias_precision_XXX.png
	- child_XXX.Rnw: separate .Rnw files to reproduce the Tables in Section 2 of
	  Supplementary Material C
- Comparison_survival_time:
	- Comparison_survival_time.R: run this script to reproduce the results presented in
	  Section 4 of Supplementary Material C
	- 'Results': outputs from "Comparison_survival_time.R"
	- compare_XXX.png: plots in Section 4 of Supplementary Material C
	- reproduce_plot_comparison_survival_time.R: run this script to reproduce the plots
	  compare_XXX.png
	- child_compare_XXX.Rnw: separate .Rnw files to reproduce the Tables in Section 4 
	  of Supplementary Material C
- Comparison_value_search:
	- comparison_value_search.R: run this script to reproduce comparison_dmsm_dwsurv.RData
	  and the two figures.
	- comparison_dmsm_dwsurv.RData: used to reproduce results in Section 5 of Supplementary Material C
	- XXX.png: plots in Section 5 of Supplementary Material C
- Estimated_optimal_DTR:
	- estimate_optimal_DTR.R: run this script to reproduce the results presented in
	  Section 3 of Supplementary Material C
	- child_ability_XXX.Rnw: separate .Rnw files to reproduce the Tables in Section 3 
	  of Supplementary Material C
    
## Set-up to reproduce the results

The content of this folder should be saved on the user's computer. In all R and .Rnw files, the variable 'mainpath' should be replaced by the location of this folder on the user's 
computer. The R package DTRreg version 1.4 must be installed.

The content of the folders 'Independent30', 'Independent60', 'X1dependent30', 'X1dependent60', 'X1X2dependent30' and 'X1X2dependent60' should be sorted by name.

R version 3.5.1 (2018-07-02) -- "Feather Spray" was used.
