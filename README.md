# VaccinationImmunityModeling
Simulations of general respiratory infection model and Colorado specific COVID19 model in MATLAB. 

sir_simulator.m and sir_simulator_masking.m contain the general model ODE system. The associated figure().m files generate the corresponding figures from the paper "Waning Immunity and Partial Vaccination Coverage Lead to Transitions in the Source of Daily Incidence." The CO Datafitting folder contains a copy of the hospitalization data, the CO specific model "CovidVaxMixModel.m", the main file for the data fitting process "fitVaxMixModel_oct2025", and the main file for generating the final figure "processVaxMixModel.m". The simbiology toolbox is required for the CO Datafitting portion of the project. 
