## Resubmission
This is a re-submission. In this version I have:


* replaced the references from Candelon and Moura (forthcoming) to Candelon and Moura (2024), since the paper is now officially part of a journal issue;
* included  S3 classes and specific methods for key objects. The package's vignette was incremented accordingly;
* slightly simplified the *LoadData* function to enhance the readability;
* improved the overall efficiency of the model optimization processes by eliminating unnecessary for loops; 
* eliminated the package dependence from the *ceil*, *eig*, *mod*, *isempty*, *lu*, *mrdivide*, *num2str*, *numel*, *rand*, *strcmp*, *strfind*, and *size* functions imported from the *pracma* package. As such, this package has been moved to the to suggests section;
* replaced the use from the *is_empty* function so that the *sjmisc* package no longer appears as a suggested package;
* replaced the use from the *str_length* function so that the *stringr* package no longer appears as a suggested package;
* replaced the use from the *Curry* function so that the *functional* package no longer appears as a suggested package;
* replaced the use from the *readr* function so that the *parse_number* package no longer appears as a suggested package;
* replaced the use from the *demean*, *tic*, *toc* functions so that the *Jmisc* package no longer appears as a suggested package;
* corrected mistakes in the computation of IRF and GIRF for the JLL related models;
* created the *YieldsFitAll* function that encompasses the previous usage of the previously named *YieldsFitAllSep* and *YieldsFitAllJoint* functions;
* reoptimized the estimation routine from the *JLL joint Sigma* model class.

