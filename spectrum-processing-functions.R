## This script contains functions that are helpful for processing data related to the 
## 'Electronic Absorption of Molecular Iodine' lab in CHEM 308.

simes <- as.data.frame(matrix(c(541.2, 27, 539.0, 28, 536.9, 29), byrow = TRUE, nrow = 3, ncol = 2))
names(simes)[] <- c("Wavelength", "v'")

## McNaught, I. J., 1980, J. Chem. Educ., Volume 57, Number 2, p. 101
mcnaught <- data.frame(
  Wavelength = c(541.2, 539.0, 536.9, 571.6, 568.6, 565.6, 595.7, 592.0, 588.5),
  XVQN       = c(    0,     0,     0,     1,     1,     1,     2,     2,     2), 
  BVQN       = c(   27,    28,    29,    18,    19,    20,    13,    14,    15)
)

decreaseResolution <- function(dataFrame, n)
  ## Creates a new dataframe by sampling dataFrame for every nth point
{
  return(dataFrame[seq(1, nrow(dataFrame), by = n), ])
}

applyBoxcarWidth <- function(dataFrame, n)
{
  nrdf <- nrow(dataFrame)
  boxcar <- dataFrame[c( (1 + n):(nrdf - n) ), ]
  nrbx <- nrow(boxcar)
  for(i in c(1:nrbx))
  {
    dfrow <- i + n
    ## Calculate mean of dataFrame rows (dfrow - n):(dfrow + n)
    rowMean <- rep(0, ncol(dataFrame)) #initialize
    for(j in c( (dfrow - n):(dfrow + n)) )
      rowMean <- rowMean + dataFrame[j, ] / (2*n + 1)
    boxcar[i, ] <- rowMean
  }
  return(boxcar)
}

peakr <- function(dataFrame, xCol, yCol, findMins = FALSE, returnIndices = FALSE)   ## What an original name
  ## Simple peak finder
{
  indicesOfExtrema <- NA ## initialize
  index <- 1
  nr <- nrow(dataFrame)
  
  for (row in c(2:(nr - 1)))
  {
    cdate <- row ## Candidate index
    cVal <- dataFrame[row, yCol] ## y-val at candidate index
    prevVal <- dataFrame[row - 1, yCol]
    nextVal <- dataFrame[row + 1, yCol]
    
    isMax <- (cVal > prevVal & cVal >= nextVal)
    isMin <- (cVal < prevVal & cVal <= nextVal)
    
    if ( (findMins & isMin) || (!findMins & isMax) )
    {
      ## Found peak from left. Add candidate to indicesOfExtrema and increment index.
      indicesOfExtrema[index] <- cdate
      index <- index + 1
    }
  }
  
  if(returnIndices) return(indicesOfExtrema)
  else return(dataFrame[indicesOfExtrema, c(xCol, yCol)])
}

maxAbsFilter <- function(dataFrame, xCol, yCol, threshold, type = "l") 
  ## Filter indices of extreme values for specific maxima
  ## type = "l", "r", or "both"
{
  minima <- peakr(dataFrame, xCol, yCol, findMins = TRUE, returnIndices = FALSE)
  maxima <- peakr(dataFrame, xCol, yCol, findMins = FALSE, returnIndices = FALSE)
  extrema <- rbind(minima, maxima)
  extrema <- extrema[order(extrema[[1]], decreasing = FALSE), ]
  nr <- nrow(extrema)
  
  indicesOfThresholdMaxima <- NA ## initialize
  i <- 1 ## indexer for indicesOfThresholdMaxima
  
  if(type == "l" || type == "left") {
    for(row in c(2:nr))
    {
      cdateAbs <- extrema[row, 2] ## candidate absorbance
      prevAbs <- extrema[row - 1, 2] ## previous absorbance
      
      if(cdateAbs - prevAbs > threshold) {
        indicesOfThresholdMaxima[i] <- row
        i <- i + 1
      }
    }
  } else if(type == "r" || type == "right") {
    for(row in c(1:(nr-1)))
    {
      cdateAbs <- extrema[row, 2]
      nextAbs <- extrema[row + 1, 2]
      
      if(cdateAbs - nextAbs > threshold) {
        indicesOfThresholdMaxima[i] <- row
        i <- i + 1
      }
    }
  } else if(type == "both") {
    for(row in c(2:(nr-1)))
    {
      cdateAbs <- extrema[row, 2]
      prevAbs <- extrema[row - 1, 2]
      nextAbs <- extrema[row + 1, 2]
      
      if(cdateAbs - prevAbs > threshold && cdateAbs - nextAbs > threshold) {
        indicesOfThresholdMaxima[i] <- row
        i <- i + 1
      }
    }
  }
  
  return(extrema[indicesOfThresholdMaxima, ])
}

indicesToPoints <- function(dataFrame, indices, xCol, yCol) 
  ## Convert vector of indices to dataframe of plottable points
{
  return(dataFrame[indices,c(xCol, yCol)])
}

removePointsByValue <- function(dataFrame, value, ltVal = TRUE, inclusive = FALSE)
  ## This function expects a dataframe of (x,y) coordinate pairs and
  ## will return a subset of the original dataframe with x-values less than 'value' removed. If 
  ## 'ltVal' is set to FALSE, then the function returns a subset with x-values greater than 'value'
  ## removed. 'inclusive' specifies whether or not to remove a point if it's x-value is exactly 'value'.
{
  rowsToRemove <- 0 ## Vector to store row numbers of rows that will be removed from dataFrame
  rIndex <- 1 ## Indexer for rowsToRemove; incremented whenever a value is stored in rowsToRemove
  
  nRows <- nrow(dataFrame)
  rows <- c(1:nRows)
  for(row in rows)
  {
    xVal <- dataFrame[row,1]
    if( (ltVal && xVal < value) || (!ltVal && xVal > value) || (inclusive && xVal == value) ) {
      rowsToRemove[rIndex] <- row
      rIndex <- rIndex + 1
    }
  }
  rowsToReturn <- setdiff(rows, rowsToRemove)
  
  return(dataFrame[rowsToReturn,])
}

removeEveryOtherAfter <- function(dataFrame, value)
  ## Removes every other point after x-value 'value'. In other words, the first x-value larger than 'value'
  ## is kept, the next is removed, the next is kept, and so on. Requires that dataFrame rows are sorted by x-value.
{
  startAtRow <- 0 ## Row at which to start the every-other removal (this row is kept, next is removed, ...)
  nRows <- nrow(dataFrame)
  rows <- c(1:nRows)
  for(row in rows)
  {
    xVal <- dataFrame[row,1]
    if(xVal > value) {
      startAtRow <- row
      break
    }
  }
  rowsToRemove <- seq(row+1, nRows, by = 2)
  rowsToReturn <- setdiff(rows, rowsToRemove)
  
  return(dataFrame[rowsToReturn,])
}

includeRows <- function(pointList, rowIndexVector, dataFrame, xCol, yCol)
  ## Choose particular rows from dataFrame to add to pointList. This function will
  ## add those rows to the end of pointList, then sort pointList by Wavelength.
{
  ## Check if any rowIndexVector elements are already in the pointList; if so,
  ## remove them from rowIndexVector
  removeRows <- c()
  i <- 1 ## indexer for removeRows
  for(row in rowIndexVector)
  {
    wavelengthToAdd <- dataFrame[row,xCol]
    for(point in c(1:nrow(pointList)))
    {
      wavelengthInPointList <- pointList[point, 1]
      if(wavelengthToAdd == wavelengthInPointList)
      {
        removeRows[i] <- row
        i <- i + 1
      }
    }
  }
  rowIndexVector <- setdiff(rowIndexVector, removeRows)
  pointsToAdd <- dataFrame[rowIndexVector, c(xCol,yCol)]
  names(pointsToAdd) <- names(pointList)
  
  ## Add the pointsToAdd to the end of pointList
  pointList <- rbind(pointList, pointsToAdd)
  ## Sort pointList by wavelength
  pointList <- pointList[order(pointList[, 1], decreasing = FALSE), ]
  
  return(pointList)
}

matchRowToLit <- function(dataFrame, litFrame, tolerance = 0.2)
  ## Assign vib. quantum numbers from literature to dataFrame.
  ## dataFrame is expected to be a dataframe that contains the (wavelength, absorbance) coordinates
  ## for peaks that correspond to v'' = 0 to v' = n transitions.
  ## litFrame is expected to be a dataframe that contains (wavelength, v') coordinates.
  ## A v' column will be appended to dataFrame, and v' values will be assigned to this column from literature
  ## if ( abs(dataFrame$Wavelength - litFrame$Wavelength) < tolerance )
{
  nRows <- nrow(dataFrame)
  nCols <- ncol(dataFrame)
  dataFrame[, nCols + 1] <- rep(NA, nRows)
  names(dataFrame)[nCols + 1] <- "BVQN"
  
  for(row in c(1:nRows))
  {
    cwl <- dataFrame[row, 1] ## candidate wavelength for literature match
    for(litRow in c(1:nrow(litFrame)))
    {
      lwl <- litFrame[litRow,1] ## literature wavelength
      if(abs(cwl - lwl) < tolerance) dataFrame[row, nCols + 1] <- litFrame[litRow,2]
    }
  }
  dataFrame <- extrapolateV(dataFrame)
  return(dataFrame)
}

extrapolateV <- function(dataFrame)
{
  ## Sort the dataFrame by values in its first column (wavelength)
  dataFrame <- dataFrame[order(dataFrame[, 1], decreasing = FALSE), ]
  
  nCols <- ncol(dataFrame) ## BVQN should be stored in last column
  maxvRow <- which.max(dataFrame[, nCols]) ## row of max v'
  minvRow <- which.min(dataFrame[, nCols]) ## row of min v'
  vDecreasesByRow <- if(maxvRow < minvRow) TRUE else FALSE
  nRows <- nrow(dataFrame)
  
  if(vDecreasesByRow) {
    ## Start at maxvRow and assign upward, then start at minvRow and assign downward
    if(maxvRow > 1) for(row in c( (maxvRow-1):1 )) {
      vToAssign <- max(dataFrame[, nCols], na.rm = TRUE) + 1
      dataFrame[row, nCols] <- vToAssign
    }
    if(minvRow < nRows) for(row in c( (minvRow+1):nRows )) {
      vToAssign <- min(dataFrame[, nCols], na.rm = TRUE) - 1
      dataFrame[row, nCols] <- vToAssign
    }
  } else {
    ## Start at minvRow and assign upward, then start at maxvRow and assign downward
    if(minvRow > 1) for(row in c( (minvRow-1):1 )) {
      vToAssign <- min(dataFrame[, nCols], na.rm = TRUE) - 1
      dataFrame[row, nCols] <- vToAssign
    }
    if(maxvRow < nRows) for(row in c( (maxvRow+1):nRows )) {
      vToAssign <- max(dataFrame[, nCols], na.rm = TRUE) + 1
      dataFrame[row, nCols] <- vToAssign
    }
  }
  
  return(dataFrame)
}
