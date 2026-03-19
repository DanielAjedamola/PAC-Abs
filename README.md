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
### To construct the corresponding abstractions
##### 4a. replace in the folder runningmax, the files runningmax.cc and Makefile with the one provided here in folder Nonlinear Interconnected dt-cs, and proceed with the same instrcution to compile the files as it was done on the referenced repository in step 1a.
##### 4b. replace the files vehicle.cc, vehicle.m with the one provided here in folder vehicle_PAC_guaratee, and proceed with the same instrcution to compile the files as it was done on the referenced repository in step 1b.
##### 4c. replace the files dcdc.cc, dcdc.m with the one provided here in folder DC-DC Boost Converter, and proceed with the same instrcution to compile the files as it was done on the referenced repository in step 1b.
### To construct the corresponding ASFs
##### 5. Run the corresponding .ipynb files in each folder for each example
