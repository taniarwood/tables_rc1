Step zero!

Pick an MCEq table you want to use 


Step one!

Reweight your pickle files with using that flux (and the class to spline the full tables and make zenith integrated splies)

Step two!

A) use a classic version of oscfit (with only one dataloader) to make the systematic parameterizations using the full pickle files (~40 Gigs!).  save them as a pickle file. From now on 
you only need to load the baseline files (with analysis keys only!) and the atmmu files (from data -- > ~3.7 gigs!)

B) use 'remove_weights.py' script to remove the numu weights from the nue files and vis versia in a version of the baseline files that will now only have numu or nue weights. we want this for 
the dataLoaders we will use in the next steps which expect nue only or numu only predictions

Step three!

Run a fit! Use an example file that shows you how to use 'extra_cuts' in oscfit to cut on true energy and load many dataLoaders, each for a true energy bin ! (NuE and NuMu seperate as well!)



Step four!

Plot output with example script!  
