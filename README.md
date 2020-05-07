# seasonal_forecast_scripts

Scripts to download, subset and process seasonal rainfall hindcasts (re-forecasts) from the Copernicus Climate Change Service Climate Data Store (https://cds.climate.copernicus.eu/#!/home)

M. Young, 7/5/2020

List of scripts:

1) download_c3s_seasonal.py (run using submit_download_c3s.sh)

Downloads global gridded seasonal hindcasts of accumulated precipitation from different models for different months.


2) save_subset_sforecast.py (run using submit_subset_sforecast.sh)

Subsets seasonal hindcasts over region of interest, at specified lead time, over a specified accumulation period.


3) save_subset_chirps.py (run using submit_subset_chirps.sh)

Subsets CHIRPSv2.0 daily precipitation data over region of interest and accumulates over period specified from data saved in script 2) and regrids this data to the same horiontal resolution as the hindcasts.


4) save_sforecast_outcomes.py (run using submit_save_sforecast_outcomes.sh)

Reads in seasonal hindcasts saved by script 2) and calculates outcomes (0 or 1) for a given event (e.g. tercile probabilities), specified by the user, for each ensemble member in the forecast. Saves each outcome for each ensemble member in temporary file.


5) save_sforecast_probabilities.py (run using submit_save_sforecast_probabilities.sh)

Reads in temporary outcome ensemble files saved using script 4) and combines these to calculate the event probabilities for the hindcasts. Saves the probabilities.


6) save_chirps_outcomes.py (run using submit_save_chirps_outcomes.sh)

Reads in CHIRPS data saved by script 3) and calculates outcomes (0 or 1) for a given event for evaluation of the hindcast probabilities for the event generated in scripts 4) and 5). Saves outcomes.


7) evaluate_sforecast_chirps.py

Reads in seasonal hindcast data saved by script 2), hindcast probabilities produced by script 5), and equivalent CHIRPS observations produced by scripts 3) and 6). Computes some standard validation statistics - mean precipitation for hindcasts and CHIRPS, bias of hindcast ensemble mean (relative to CHIRPS), hindcast Anomaly Correlation Coefficient (ACC, relative to CHIRPS) and Relative Operating Characteristic (ROC) Area for hindcast probabilities. Plots these statistics as spatial maps over region of interest.


8) average_sforecast_probabilities_region.py

Reads in season hindcast probabilities produce by script 5), subsets these over a region specified by a shapefile (e.g. country border) and computes the average hindcast probabilities over this region. Saves regionally averaged hindcast probabilities in a text file.  
