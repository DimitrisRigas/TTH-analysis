# TTH-analysis

myplot.C plots the gen level diagrams.

reconplot.C plots the recon diagrams

1/12/2025.Version 1.0. //Code cleanup. Gen analysis done. Recon analysis up to sorted b-jets completed

2/12/2025 Version 2.0 //Got a new dataset. Mass of a=12Gev

2-10/12/2025 Trying to adapt the changes and fix all the bugs that come up.

11/12/2025 Version 2.1 //Added cross cleaning,Higgs kinematics, a kinematics, Ht, Met. Added background data. Need more background data because only 4 events remain after cuts. Also updated the plot macros to get the sample size from the corresponding output file depending on the input.

12/12/2025 2.2 // Added cross cleaning on the jet selection choice for b jets and added more kinematic variables

13/12/2025 2.3 // Added more kinematic variables and 2 mroe selection cuts

14/12/2025 2.4 // Fixed the cross cleaning of jet logic. Issue found.  I have way too few events surviving. Solution needs to be found.

7/1/2025 2.5 // Became single leptonic to fix the issue with few events surviving.  Now i have 7 events on signal and 7 events on background of TTto2L2nu. Added a new output file where it saves the kinematic variables as trees . Fixed an issue that came up with the myplot.C.

8/1/2025 Though of a solution in my efforts to save the dileptonic case. Instead of doing boosted and getting 2 double b tagged jets im going do do a semi try where i have atleast 1 double b jet and atleast 3 bjets. Need to think of a way to select the correct b jet though for the reconstruction of Higgs variables because i have the same problem i had with mass 60.  I did atleast one double b jet and although i get more events i lose all of them on the cut of atleast 3 bjets. It could be that im on mass 12 and i have more fully boosted incidents. I will try now getting mass 25 Gev- I restored the dileptonic analysis and saved as a new name Myclass2.C

9/1/2025 Tried going for a selection cut Nb>1. I x5 the signal and my bg survived.
