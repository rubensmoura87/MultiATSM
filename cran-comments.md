## Resubmission
This is a resubmission. In this version I have:

* removed "This package" as the initial words from the description;

* written the package description containing capital letters only in the begining of sentences and names; 

* added the references (doi, https) to the quoted research papers

* specified \value for DataForEstimation.Rd and set InputsForMLEdensity_BS.Rd as non-exportable;

* unwraped the examples for the functions which are executable functions in < 5 seconds. For the others, I have set "\donttest{}";

* removed the commands that were saving files on the user's home filespace. In its current form, the package saves these files on the user's temporary directory; 

* replaced "option(warm=-1)" by "suppressWarnings()" in the Transition_Matrix.Rd.
