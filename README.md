# When can we trust population trends?

This is the code associated with the paper:

Wauchope, H. S., Amano, T., Sutherland, W. J., Johnston, A. (2019) "When can we trust population trends? A method for quantifying the effects of sampling interval and duration" *Methods in Ecology and Evolution*, Early View. 

The code is adaptable to any dataset, but is admittedly fairly long and complex. A description of the components of the code, required data, and how to adapt, see the text at the beginning of the "WauchopeQuantifyingTrends_1.R" file. 

To aid with trouble shooting, there is a Dummy Data file that represents the actual data used with the code. I took the actual data (which requires permission from the Christmas Bird Count to access) and wrote over everything, so they are fake Population names, randomly simulated Negative Binomial counts, Years from 1:30 and Hours set to 1. It’s in an RData format, object name “TYB_Dummy”.
 
It's the same length as the actual data used with the code just to try and keep it as consistent as possible for error checking, but obviously you can reduce it to just a few populations so it doesn’t take a long time to run.
