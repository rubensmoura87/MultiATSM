## Resubmission
This is a re-submission. In this version, I have made substantial changes to the previous version of the package. This includes the:

* update the DOI from Candelon and Moura's paper (forthcoming, JFEC) in the package description, since this paper has now been published;
* change in the labels of the different classes of models to make them more intuitive for the users; 
* inclusion of the *LoadData* function to simplify the process of loading the required dataset used in the several examples;
* development of the *InputsForOpt* function (and its several sub-functions) to better handle the model inputs that need to be used for the estimation of the several model types;
* reshape the *Optimization* function so the optimization routines of the model point estimates, bootstraps and out-of-sample forecast are now unified;
* adaptation of the *Bootstrap* function (and its several sub-functions) so now the confidence intervals of all models can be estimated within the same function set;
* adaptation of the *ForecastYields* function (and its several sub-functions) so now the bond yield forecasts of all model classes can be estimated within the same function set;
* Within the *NumOutputs* function, I have fixed a bug in the *TermPremiaDecompSep* function;
* Finally, the package's vignette and the exportable function documentations have been changed accordingly to reflect all the changes mentioned. Overall, the package interface has been greatly simplified, as it can be noted in Section 4 of the package's vignette.
