## Resubmission
This is a re-submission. In this version I have:

* made the *RiskFactors* plot as an optional feature to the user;
* improved the overall readability of the *Term Premia* and *Fit* graphical functions by removing the legend boxes and the word "Legend". Additionally, I have adapted the colors in the *Term Premia* graph to be color-blind friendly;
* performed several changes to the *Bias_Correc_VAR* function in order to make it more robust;  
* made the *StarFactors* internal, since this function is not called by *GVAR*. Vignette has been adapted accordingly;
* fixed a bug in the *JLL* function for the case $M = 1$;
* fixed a bug in autoplot method to plot the *RiskFactors* type;
* included unit tests to all user facing functions;
