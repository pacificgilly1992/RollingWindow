# Data averaging using a rolling window
Determining the relationship between two variables can be difficult either because of a lack of data or when the data contains many dependencies. PG can be influenced by many meteorological phenomena, such as rain, snow, fog, charged clouds, conductivity changes, corona and lightning. This PhD analysed the relationship that rainfall had on the PG. Rainfall often occurs simultaneously with other meteorological events causing an increased variability in the PG. Subsetting or acquiring more data was not possible in this situation. Therefore, an algorithm was designed to improve the accuracy of the PG and rainfall relationship by averaging the data using a rolling window technique and allowed a greater number of data points to be grouped together. Overall, the rolling window data averaging method maximises the resolution of a combined data set.

# Code
The development code can be found under EPCC_PGRR_Ensemble.py and can be used to average a 2D dataset using three methods.

# Methods
1) Normal averaging - This method will bin all the data into equally sized bins. This method is similar to the binning method used for the basic creation of a histogram. The main difference here is that each bin is then averaged.
2) Variable window averaging - This method will create bins of equal number of data points rather than keeping the bin width constant. This has the advantage of normalising the averaging procedure across the entire dataset and is useful when there is a strong variation in data availablity.
3) Rolling window averaging - This method will create bins based on the previous method but will create a rolling window over the entire dataset. This has the adavantage of maximing the resolution of the dataset you want to average.

For example if we have 100 datapoint and want to have 10 bins, we will select the first 10 datapoints (points 1-10). Then we will roll our bin window by an arbitary amount (say by 1 datapoint). Our next bin will select 10 datapoints starting from point 2 (points 2-11). This process is then repeated (e.g. bin 1 = (1-10), bin 2 = (2-11), bin 3 = (3-12), bin n = (n-(n+m-1)) where n is the bin index and m is the number of datapoints per bin).

A disadvantage of the rolling window method is the boundaries of the dataset will have a poorer resolution. Going back to our example of selecting 10 datapoint bins, for datapoints 5-95 we have a maximum resolution, but for datapoints 1-4 and 96-100 our resolution is degraded or lost as we have a lower number of datapoints to average in comparison.

# Notes
1) In this algorithm we average the data in each bin using either the mean or the median (N.B. the median is good to use when averaging data with large outliers)
2) This algorithm also calculates the standard error of each bin based upon the data being normally distributed. Currently there is no attempt at calulating the standard error for any other data distributions. This could be achieved using bootstrapping for example.
