## Functions for data analysis

wavelengthToWavenumber <- function(wavelength)
  ## Calculates wavenumber (inverse cm) from wavelength (nm)
{ return(10000000/wavelength) }

birgeSponer <- function(procData)
  ## Expects procData to be a dataframe with the following form:
  ##     Wavelength Absorbance v'
  ## 1   ...        ...        ..
  ## 2   ...        ...        ..
  ## .
  ## Returns (v' + 1), delta(nu?) values for a Birge-Sponer plot
{
  nRows <- nrow(procData)
  procData[,4] <- rep(NA, nRows)
  procData[,5] <- rep(NA, nRows)
  procData[,6] <- rep(NA, nRows)
  names(procData)[4] <- "Wavenumber"
  names(procData)[5] <- "deltaNu"
  names(procData)[6] <- "v' + 1"
  for(row in c(1:nRows)) procData[row,4] <- wavelengthToWavenumber(procData[row,1])
  for(row in c(1:nRows-1)) procData[row,5] <- procData[row,4] - procData[row+1,4]
  for(row in c(1:nRows-1)) procData[row,6] <- procData[row,3] + 1
  
  return(procData)
}

secondOrder <- function(procData)
  ## Expects procData to be a dataframe with the following form:
  ##     Wavelength Absorbance v'
  ## 1   ...        ...        ..
  ## 2   ...        ...        ..
  ## .
  ## Returns (v' + 1/2), wavenumber values for a second order plot
{
  nRows <- nrow(procData)
  procData[,4] <- rep(NA, nRows)
  procData[,5] <- rep(NA, nRows)
  names(procData)[4] <- "Wavenumber"
  names(procData)[5] <- "v' + 1/2"
  for(row in c(1:nRows)) {
    procData[row,4] <- wavelengthToWavenumber(procData[row,1])
    procData[row,5] <- procData[row,3] + 0.5
  }
  
  return(procData)
}
