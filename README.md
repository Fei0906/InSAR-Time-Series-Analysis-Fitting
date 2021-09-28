# InSAR Time Series Analysis Fitting

# Introduction:

The shared code and data here is an example used for coseismic deforamtion field reconstruction via InSAR time series fitting. InSAR data was processed by LiCSAR (Lazecky et al., 2020) and StaMPS software (Hooper et al., 2007), using the descending track 006 InSAR time series data from Oct 2014 to Sep 2019 in Iran. We used this data to detect and model the shallow earthquakes happened in Iran, following the Mw 7.3 mainshaock in 2017.

# Files:

    Export: contains some empty files that used to read the SLC date during the programm
    
    parms_aps.mat: StaMPS file, here is used to read the UTC time of the acquistion time
    
    ref*.mat: reference pixels, with and without GACOS tropospheric corretion (Yu et al., 2018). Storing the information of reference area, including the lon, lat, radius of the reference circle area, number of referencing pixels, rms of a linear fitting to the reference pixels, and the LOS of reference pixels.
    
    TSA_ifg_ds.mat: the time series data on deformation area, downsampled by ratio 10 for quick data sharing
    
    Semi.mat: store the semi-variogram fitting (Webster & Oliver, 2007) results using the GBIS software (Bagnardi & Hooper, 2018), later used as the covariance matrix. Here, the sill values of the mainshock interferogram is mannually set as the avarage values of all other interferograms, as the significant coseismic deformations greatly bias the estimation.
    
    runTSA.m: the main function
    
    TSA_EQ_fit.m: the fitting function
    
    TSA_EQ_pixel.m: the plotting function used to plot the results

# How to run: 

  Simply run the runTSA.m and the reconstructed coseismic deformation field will pop out, and you can further click the 'TSA_EQ' button at the bottom left cornor to plot the         time series data for individual pixels.

# Please cite the following paper if you use our code:

Liu, F., Elliott, J. R., Craig, T. J., Hooper,A., & Wright, T. J. (2021). Improving the resolving power of InSAR for earthquakes using time series: A case study in Iran. Geophysical Research Letters, 48, e2021GL093043. https://doi.org/10.1029/2021GL093043

Fei Liu, Sep 2021

