## Resubmission
This is a re-submission. In this version I have:

* made the *RiskFactors* plot as an optional feature to the user;
* improved the overall readability of the *Term Premia* and *Fit* graphical functions by removing the legend boxes and the word "Legend". Additionally, I have adapted the colors in the *Term Premia* graph to be color-blind friendly;
* performed several changes to the *Bias_Correc_VAR* function in order to make it more robust;  
* made the *StarFactors* internal, since this function is not called by *GVAR*. Vignette has been adapted accordingly;
* fixed a bug in the *JLL* function for the case $M = 1$;
* fixed a bug in autoplot method to plot the *RiskFactors* type;
* included unit tests to all user facing functions;
* replaced the use of *Convert2JordanForm* and *Reg_K1Q* functions by the *FeedMat_Q* function and its sub-functions in order to enhance the modularity of the code;
* introduced a new helper function, *Load_Excel_Data*, to simplify retrieving data from Excel files. This avoids performing file I/O inside other functions, which now accept data frames or matrices directly as inputs. In this context, I refactored *DatabasePrep*, *Transition_Matrix*, *SpecificMLEInputs*, and *DataForEstimation* so that they no longer take file paths as arguments
