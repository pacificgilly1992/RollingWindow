# Data averaging using a rolling window
Determining the relationship between two variables can be difficult either because of a lack of data or when the data contains many dependencies. PG can be influenced by many meteorological phenomena, such as rain, snow, fog, charged clouds, conductivity changes, corona and lightning. This PhD analysed the relationship that rainfall had on the PG. Rainfall often occurs simultaneously with other meteorological events causing an increased variability in the PG. Subsetting or acquiring more data was not possible in this situation. Therefore, an algorithm was designed to improve the accuracy of the PG and rainfall relationship by averaging the data using a rolling window technique and allowed a greater number of data points to be grouped together. Overall, the rolling window data averaging method maximises the resolution of a combined data set.

# Code
The development code can be found under EPCC_PGRR_Ensemble.py and can be used to average a 2D dataset using three methods.

1) Normal averaging - This method will bin all the data into equally sized bins. This method is similar to the binning method used for the basic creation of a histogram. The main difference here is that each bin is then averaged.
2) 
