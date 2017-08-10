# Load in the CBI function
contBoyce <- function(predAtPres, predAtBg, numClasses=10) {
  # contBoyce  Calculates the continuous Boyce index, a measure of model accuracy for presence-only test data (Hirzel, A.H., Le Lay, G., Helfer, V., Randon, C., and Guisan, A.  2006.  Evaluating the ability of habitat suitability models to predict species presences.  Ecological Modeling 199:142-152 ~and~ Boyce, M.S., Vernier, P.R., Nielsen, S.E., and Schmiegelow, F.K.A.  2002.  Evaluating resource selection functions.  Ecological Modeling 157:281-300.  This function uses predictions at a large number of randomly located points (e.g., 10,000) to estimate the proportional coverage of each prediction class rather than actual rasters.  Hence, the index is somewhat approximate, though in practice it is nearly the same as what would be calculated using a prediction raster.
  #
  # ARGUMENTS
  # predAtPres
  # list of predicted values at test presences
  #
  # predAtBg
  # list of predicted values at a large number of randomly located sites
  #
  # numClasses
  # number of classes into which to divide predictions at background sites... Hirzel et al. suggest using 10
  #
  # VALUES
  # Spearman rank correlation coefficient between the proportion of sites in each prediction class and the expected proportion of predictions in each prediction class based on the proportion of the landscape that is in that class (=the continuous Boyce index)
  # 
  # REQUIRED DEPENDANCIES
  #
  #
  # OPTIONAL DEPENDANCIES
  #
  #
  # BAUHAUS
  # 
  #
  # EXAMPLE
  # contBoyce(predAtPres=c(runif(10), predAtBg=runif(10000)^2, numClasses=10) {
  #
  # SOURCE	source('C:/ecology/Cloud Drive/r/SDM/SDM - Calculate Continuous Boyce Index.r')
  #
  # TESTING
  # Tested with data from Steven Callan versus Excel version
  #
  # LICENSE
  # This document is copyright ©2012 by Adam B. Smith.  This document is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3, or (at your option) any later version.  This document is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. Copies of the GNU General Public License versions are available at http://www.R-project.org/Licenses/.
  #
  # AUTHOR	Adam B. Smith | Missouri Botanical Garden, St. Louis, Missouri | adamDOTsmithATmobotDOTorg
  # DATE		2013-03
  # REVISIONS 2013-07-09	Now removes NAs from predictions before calculating CBI.
  #			2014-02-27  Now expands lower end of classes by very small amount to correct for rounding error in floating point precision of predicted values... otherwise can produce NA value if all predictions at presence sites are equal (happens sometimes if predictors are categorical)
  
  #############################
  ## libraries and functions ##
  #############################
  
  ####################
  ## pre-processing ##
  ####################
  
  ## remove NAs
  predAtPres <- predAtPres[!is.na(predAtPres)]
  predAtBg <- predAtBg[!is.na(predAtBg)]
  
  ## calculate right hand side of each bin (expanding just a little to include max prediction in final bin)
  rangeOfPredVals <- (max(c(predAtPres, predAtBg), na.rm=T) - min(c(predAtPres, predAtBg), na.rm=T))  # range
  
  # right hand side of each class (assumes max value is >0)
  classEnd <- seq(
    from=(min(c(predAtPres, predAtBg), na.rm=T) + 2 * rangeOfPredVals/(numClasses + 1)) * 0.999999,
    to=max(c(predAtPres, predAtBg), na.rm=T) * 1.000001,
    length.out=numClasses
  )
  
  ## initiate variables to store predicted/expected (P/E) values
  freqPres <- freqBg <- numeric()
  
  ##########
  ## MAIN ##
  ##########
  
  ## tally proportion of test presences/background sites in each class
  for (countClass in 1:numClasses) {
    
    # number of presence predictions in this class
    freqPres <- c(
      freqPres,
      sum(predAtPres >= classEnd[countClass] - 2 * (rangeOfPredVals / (numClasses + 1)) & predAtPres < classEnd[countClass])
    )
    
    # number of background predictions in this class
    freqBg <- c(
      freqBg,
      sum(predAtBg >= classEnd[countClass] - 2 * (rangeOfPredVals / (numClasses + 1)) & predAtBg < classEnd[countClass])
    )
    
  } # next predicted value class
  
  # add small number to each background frequency class to avoid division by 0 ("small" number is "half" a presence)
  # freqBg <- freqBg + min(0.5, min(freqBg[freqBg > 0]))
  freqBg <- freqBg + 0.5
  
  # calculate PE vector
  PE <- (freqPres / sum(freqPres)) / (freqBg / sum(freqBg))
  
  # calculate continuous Boyce index (CBI)
  CBI <- cor(
    x=1:length(PE),
    y=PE,
    method='spearman'
  )
  
  #####################
  ## post-processing ##
  #####################
  
  return(CBI)
  
}