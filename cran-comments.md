## Resubmission
This is a re-submission. In this version I have:

* updated the package's description;
* made the *RiskFactors* plot as an optional feature to the user;
* improved the overall readability of the *Term Premia* and *Fit* graphical functions by removing the legend boxes and the word "Legend". Additionally, I have adapted the colors in the *Term Premia* graph to be color-blind friendly;
* performed several changes to the *Bias_Correc_VAR* function in order to make it more robust;  
* made the *StarFactors* internal, since this function is not called by *GVAR*. Vignette has been adapted accordingly;
* fixed a bug in the *JLL* function for the case $M = 1$;
* fixed a bug in autoplot method to plot the *RiskFactors* type;
* included unit tests to all user facing functions;
* replaced the use of *Convert2JordanForm* and *Reg_K1Q* functions by the *FeedMat_Q* function and its sub-functions in order to enhance the modularity of the code;
* introduced a new helper function, *Load_Excel_Data*, to simplify retrieving data from Excel files. This avoids performing file I/O inside other functions, which now accept data frames or matrices directly as inputs. In this context, I refactored *DatabasePrep*, *Transition_Matrix*, *SpecificMLEInputs*, and *DataForEstimation* so that they no longer take file paths as arguments;
* added a pkgdown site for the package;
* greatly simplified the optimization routine by replacing the use from *Aux_BoundDiag*, *bound2x*, *buildx_vec*, *x2bound*, *pos2x*, *sqrtm_robust*, *True_BoundDiag* and *x2pos* functions integrating them into sevaral other parts of the code;
* made several improvements in the optimization routines. This includes: (i) the development of the *safe_solve* function to help with the inversion of ill-defined matrices; (ii) the implementation of the *check_numeric* function that helps to identify numerical issues along the optimization process; (iii) the replacement of package-based optimization algorithms to base R one, therefore the packages *neldermead* and *fminunc* were removed from *Imports*; (iv) the addition of a new input to the *Optimization* function, so the user can now decide between *L-BFGS-B* and *Nelder-Mead* algorithms or both;
* removed *StarFactor* function since this is no longer needed for the estimation of GVAR setups;
* removed the dataset named *JPSrep*, since this is no longer needed;
* reformatted the whole package using the *styler* package; 
* replaced the binary combination 0/1 for FALSE/TRUE in the package documentation, as well as variable permissible choices;
* included variable dimensions in the documentation of user-facing functions containing matrices;
* improved the data documentation in several dimensions;
* added several files to the data-raw folder in order to replicate the various datasets available in the package;


## Local Checks

I have checked the package locally with:

`devtools::check(args = c("--as-cran"))`

and it passes with **no errors, warnings, or notes**.

I have also verified Rd and documentation integrity using:
`roxygen2::roxygenize()`
`tools::checkRd()`
Both checks completed without issues.

## Additional Information

The only message appearing on CRAN’s automated checks seems to be a **platform-specific HTML rendering warning** under `r-devel-linux-x86_64-debian-gcc`:

*AdjustOptm_BS.html:62:1 (AdjustOptm_BS.Rd:18): Warning: missing </ul> before </div>*
*...*

This message does not appear under other platforms or local CRAN checks.
