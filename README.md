# tables_rc1


Step zero!

Pick an MCEq table you want to use 


Step one!

Reweight your pickle files using that flux (and the class to spline the full tables and make zenith integrated splies)

Step two!

A) use a classic version of oscfit (with only one dataloader) to make the systematic parameterizations using the full pickle files (~40 Gigs!).  save them as a pickle file. From now on 
you only need to load the baseline files (with analysis keys only!) and the atmmu files (from data -- > ~3.7 gigs!)
The oscfit used for this purpose in this analysis is kept here: https://github.com/trwood84/oscFit_3D_v2.0_andris_daydreamVerson_for_makingSYSFILES  
Check the scripts/analysis_ folder for examples and SYSFILES you can use 

B) use 'remove_weights.py' script to remove the numu weights from the nue files and vis versia in a version of the baseline files that will now only have numu or nue weights. we want this for 
the dataLoaders we will use in the next steps which expect nue only or numu only predictions

Step three!

Run a fit! Use an example file that shows you how to use 'extra_cuts' in oscfit to cut on true energy and load many dataLoaders, each for a true energy bin ! (NuE and NuMu seperate as well!)



Step four!

Plot output with example script!  



NOTES!

Currently in use:  DPMJETIII_h3a_rc1. Look here for latest and greatest, up to date, what was used in final analysis
