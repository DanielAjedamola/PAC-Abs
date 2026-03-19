### This repository contains code for the paper titled: 'DATA-DRIVEN CONSTRUCTION OF FINITE ABSTRACTIONS FOR INTERCONNECTED SYSTEMS: A COMPOSITIONAL APPROACH'
### Link: 
### Journal: Submitted to System and Control Letters
### Authors: Daniel Ajeleye*, Anas Makdesi**, and Majid Zamani*
### Contact: daniel.ajeleye@colorado.edu
### Affiliation: University of Colorado Boulder*, USA;  Ludwig Maximilian University of Munich, Germany**

##### Before doing anything with the provided code, please do the following:
##### 1a. Visit here: https://github.com/ericskim/hscc18_modular/tree/hscc18_repeatability for the details on installation and usage of this modified version of tool SCOTS, and ensure the consensus and runningmax examples are working on your device.
##### 1b. Visit here: https://github.com/mkhaled87/scots-ready for the details on installation and usage of SCOTS, and ensure the vehicle example is working on your device.
##### 2. Ensure gurobi is installed on your device, and place the header files gurobi_c.h, gurobi_c++.h, and the license file gurobi.lic in the directory where the .ipynb files can access. 
##### 3. ensure that the directory in the make file (for the .cc files) are properly adjusted as applicable on your local device.
##### 4. 

##### 4. replace in the folder runningmax, the files runningmax.cc and Makefile with the one provided here, and proceed with the same instrcution to compile the files as it was done on the referenced repository in step 1a.
##### 5. replace in the folder consensus, the files consensus.cc, visualize_consensus.py and Makefile with the one provided here, and proceed with the same instrcution to compile the files as it was done on the referenced repository in step 1.
##### 6. a third example was added in folder 'tank', which is a network of fluid-tank system, and it works similarly as those above.

