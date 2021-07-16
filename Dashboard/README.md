# Dashboard Installation and Running

To run the dashboard you must install R and we recommend installing RStudio.  The dashboard has been developed using R Version 4.1.0 and requires an active internet connection to run.

1. Download R from https://mirror.its.dal.ca/cran/ 
    - Downloads available for Linux, Mac, and Windows

2. Install R by following the installation instructions for your operating system
    - In Windows you will simply need to run the executable downloaded

3. Download RStudio Desktop from https://rstudio.com/products/rstudio/download/#download
    - This should recognize your operating system and provide a link to the correct version of RStudio.
    - RStudio works with Linux/Mac/Windows operating systems

4. Install RStudio Desktop by following the installation instructions for your operating system
    - In Windows you will simply need to run the executable downloaded
  
5. If you do not already have a copy, download the Spatial Closures Dashboard file found at
    - https://github.com/Dave-Keith/Paper_2_SDMs/tree/master/Dashboard/
      - Right click the file **Spatial_closures_dashboard.Rmd** and select *'save link as'* to download the stand-alone file
    - For users familiar with Github you can simply clone/download/fork this repo
        - https://github.com/Dave-Keith/Paper_2_SDMs/
  
6. Open RStudio and then open downloaded Spatial_closures_dashboard.Rmd file
      - File -> Open

7. Press "Run Document" button in the menu bar
      - Note that this step will automatically install a number of the required R packages for this dashboard to run
  
8. The first time you run the Dashboard it should take less than 5 minutes to load.
      - Changing settings within the dashboard can take up to a minute to load when the random fields are updating.

# Using the Interactive Dashboard 

There are several components to the dashboard, a brief overview of the layout

### Model Output

**This includes a comprehensive overview of the covariate and random field results for selected models along with the Intercept only model.**

  - The user can select results from various combinations of
    - *Stock* (cod vs yellowtail)
    - *Season* (Winter, Spring, or Fall)
    - *Era length* (3,5, or 10 years)
    - *Covariates* (Intercept, SST + Dep, or SST + Dep + Sed)
    - *Estimator* (Response, Link, Raw, or Standard Deviation)
      - Response is the most intuitive as it gives the effect size on the proportional scale
  - There are two tabs available
    - *Covariates* tab includes Intercept, SST, Depth, and Sediment type results as appropriate for the chosen model
    - The *Random fields* tab shows the results for each of the random fields used in the model

### Model Diagnostics

**This includes an overview of Model selection results for a subset of the models and shows the results for the 5-fold validation**
  
  - There are two tabs available
    - *Diagnostics* provides the model diagnostic results for a subset of models
      - Statistics include *WAIC*, *CPO*, *DIC*, *RMSE*, and *MAE*
    - *Validation* provides the model validation results for each stock
      - Validation Statistics include *Raw Error*, *RMSE*, *MAE*, and *SD*
      
### Spatial Results 

**This shows the Occurrence Probability (*OP*) results using the prediction grid for each stock using the *final model* **

  - The user can change
    - *Occurrence Probability* all *OP* values > than the maximum value on the slider are used for the figures
    - *Stock* controls the stock of interest
    - *Season* controls the season used 
  - There are four tabs available
    - *Figure* shows a mutli-panel figure of the predictions for each era
    - *Video* is the same information but shown in a looped video
    - *Center of Gravity* shows the COG by era and includes a measure of uncertainty around the COG (+-3 SE)
    - *Time Series* compares the change in area above the selected *OP* threshold (i.e. *core area*) across eras for all of Georges Bank, in Canadian waters, and in U.S. waters on GB

### Closure Time Series

**This compares the change in area above the selected OP threshold (i.e. *core area*) across eras for several area subsets including the various closures on GB**

 - The user can change
    - *Occurrence Probability* all *OP* values > than the maximum value on the slider are used for the figures
    - *Stock* controls the stock of interest
    - *Season* controls the season used 
    - *Closure* controls the closure comparison
      - *None selected* compares the trends in *core area* in U.S. waters with the Canadian Offshore Scallop Fishery (COSF) *core area* trend
      - *Closed Area I* compares the trends in *core area* in U.S. waters with the CA I *core area* trend
      - *Closed Area II* compares the trends in *core area* in U.S. waters with the CA II *core area* trend
      - *Cod Closure* compares the trends in *core area* in COSF with the Cod Closure *core area* trend
      - *Yellowtail Closure* compares the trends in *core area* in COSF with the Yellowtail Closure *core area* trend
      - *All closures* does nothing
      
### Closures Spatial

**This shows the prediction field for a given era with the spatial overlap of the *core area* for the various closures and the survey data points**

 - The user can change
    - *Occurrence Probability* all *OP* values > than the maximum value on the slider are used for the figures
    - *Stock* controls the stock of interest
    - *Season* controls the season used 
    - *Closure* controls which closures to overlay on the map
      - *None selected* shows the base map of *OP*
      - *Closed Area I* displays the location of CA I
      - *Closed Area II* displays the location of CA II
      - *Cod Closure* displays the Cod Closure locations
          - Changing *Survey Observations* slider will change the closed cells displayed
          - For example if the *Survey Observations* slider is set to 2012-2014 only the cells closed between 2012-2014 will be shown
      - *Yellowtail Closure* displays the Yellowtail Closure locations
        - Changing *Survey Observations* slider will change the closed cells displayed
          - For example if the *Survey Observations* slider is set to 2012-2014 only the cells closed between 2012-2014 will be shown
      - *All closures* displays all of the closures
    - *Survey Observations* slider
      1. Controls the era of the prediction field mapped
          - For example, setting the minimum of the slider to 1988 will show the 1987-19XX prediction field
      2. Controls the years of the survey points shown in the figure
    - *Show survey point* allows the user to remove the survey points from the figure
  - There are two tabs available
    - *Spatial Predictions* shows the survey points as presence (yellow) or absence (grey)
    - *Spatial Residuals* shows the survey points as model residuals (data up to 2016) or prediction error (2017-2019 data)
  
  


        
      
