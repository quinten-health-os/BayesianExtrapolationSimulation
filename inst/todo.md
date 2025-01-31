## For future work:
 

[] Unit tests

[] Use consistent naming for methods (with Gaussian or Binomial as suffix)

[] Remove target_data$standard_deviation as it is only used for frequentist OCs computation and may introduce some confusion (as instead we use target_data$sample$standard_error)

[] When generating data, remove if (self$summary_measure_likelihood == "normal" && self$sampling_approximation == FALSE), since the data generation should not depend on the summary measure likelihood, but on the sampling approximation (which should be either "Gaussian", or something else)

[] Prior predictive checks for models that use Stan

[] Transform data to handle non zero theta and left null hypothesis space in the Model class

[] Create a print_model_summary method for each model 

[] For consistency across methods, use self$parameters$equivalence_margin instead of self$equivalence_margin for all methods.

[] Generalize the code to allow for different sample sizes in the arms of the target study

[] Implement the NPP in the binary case

[] Implement Test then pool in the binary case, without approximation


[] Make the code more robust using private attributes

[] Use source_data$sample as well, for consistency

[] Make the assumption about the known variance more explicit in the code.

[] Make it clear somewhere (in the vignettes) that the sample standard deviation for the generation of target study data corresponds to the sample std in adults for continuous data.

[] Distinguish public/private/active classes elements

[] Change the TargetData class so that it can be defined without source data as input. 

[] Solve warning of the form : replacing previous import ‘jsonlite::unbox’ by ‘rlang::unbox’ when loading ‘RBExT’

[] Change the drift range computation in the Aprepitant case to constrain the number of drift values (instead of filtering)

 
