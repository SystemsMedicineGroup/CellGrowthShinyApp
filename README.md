# CellGrowthShinyApp
Modelling Tool for in vitro Cell Growth

# A generic Modelling Shiny App for Trr1 *In Vitro* Data
This Shiny app is designed for the visualization and modelling of growth/viability data obtained from **in vitro** experimental results.  
It allows users to upload experimental data, visualize growth and viability trends, and perform model fitting using MCMC.
  
## Features
- Upload custom CSV datasets for optical density or viability
- Visualize both data types interactively
- Select the model type according to the shapes of the OD plots (Fit experimental data using a generic growth–decay model)
- Set initial values for the model parameters (transition rate, growth rate, carrying capacity and decay rate)
- Adjust the number of MCMC interations to the desired value
- Run the model (subsequent MCMC simulations for parameter estimation are enabled)

# Description of the files
## Data set (csv)
- All In vitro OD data July 2025.csv
- It has experimental ID (exp_id), 
              date of experiment	(date), 
              replicate ID	(rep_id), 
              time points (time), 
              strain type	(strain), 
              concentration of doxycycline added (dox), 
              initial optical density of the broth at time 0 (OD0) and 
              measured optical density	(OD)

## FACS data (csv)
- FACS data.csv
- It has time points (time), 
              concentration of doxycycline added (dox) and 
              proportion of dead cells determined by flow cytometry	dead_cell_prop

## The R script for the ShinyApp
- In vitro Growth Assay_Shinyapp.R

## Dynamic Library Links (for the Fortran subroutines)
- mod1.dll (dll for the model without transition)
- mod1.f90
- mod1.o
- mod1t.f90
- mod2.dll
- mod2.f90 (dll for the model with transition)
- mod2.o
- params_mod.mo

# Note
- Make sure both data types are plotted and the **OD type** is selected before clicking **Run Model**
