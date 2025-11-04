## Resubmission
This is a re-submission. In this version I have:

* replaced the use of the *aes_string()* function from the ggplot2 package, since this function has been deprecated;
* have modified the lines for both the upper and lower confidence interval bounds to be dashed gray lines within the bootstrap setup. Additionally, I have also increased the thickness of the lines for the median and point estimates to make them more easily distinguishable from the confidence bounds;
* adjusted the colors FEVD and GFEVD graphs (point estimate version) to be color-blind friendly;
* modified the FEVD and GFEVD graphs (point estimate version) to start y-axis at zero;
* created the *clean_labels* function to generate legend labels that are automatically cleaned and made more readable, used for generating FEVD and GFEVD graphs.
