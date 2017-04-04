## Create axis to show vibrational quantum numbers (vqn)
createAxis <- function(df, yLow, yHigh, yLabel, startLabelAt = 1, labelEvery = 5, labelBelow = TRUE, mainLabelShift = 2)
  ## head(df)
  ## > Wavelength  bx2Absorbance  XVQN  BVQN
  ## >        ...            ...   ...   ...
{
  # if(startLabelAt > labelEvery){
  #   print("Error: createAxis(): startLabelAt must be less than labelEvery.")  ## To preserve my sanity
  #   return(NA)
  # }

  nRows <- nrow(df) ## Number of x-axis tick marks
  
  if (nRows == 1) {
    print("You don't really need an entire axis for a single point. Have you considered abline()?")
    return(NA)
  } else {  ## (nRows > 1)
    
    # nl <- 0 ## Number of labelled tick marks
    # 
    # if( (nRows - startLabelAt + 1) %% labelEvery == 0) { ## Perfect division
    #   nl <- (nRows - startLabelAt + 1) %/% labelEvery
    # } else { ## Imperfect division
    #   nl <- (nRows - startLabelAt + 1) %/% labelEvery + 1
    # }
    
    distanceUnit <- abs(yHigh - yLow)/3
    r <- (nRows - startLabelAt + 1) %% labelEvery ## remainder
    
    ## Initialize elements of return list
    vertices <- data.frame(AxisX = NA, AxisY = NA)
    labels <- data.frame(LabelX = NA, LabelY = NA, Label = NA, hJust = NA, vJust = NA, stringsAsFactors = FALSE)
    segments <- data.frame(SegmentX = NA, SegmentY = NA)
    
    vIndex <- 1
    lIndex <- 1
    sIndex <- 1
    
    ## Create label for entire axis (this will be the first element in the 'labels' frame)
    if(labelBelow)
      labels[lIndex, ] <- list(df[nRows, 1] + mainLabelShift,  yLow, sprintf("v'' = %d", df$XVQN[1]), 0, 0)
    else
      labels[lIndex, ] <- list(df[nRows, 1] + mainLabelShift, yHigh, sprintf("v'' = %d", df$XVQN[1]), 0, 0)
    
    lIndex <- lIndex + 1
    
    ## Create axis vertices, axis labels, and line segments
    for(row in c(nRows:1)) ## Count down since lowest vqn occurs last in list
    {
      if(row %% labelEvery == r) ## Labeled tick-mark
      { 
        vertices[vIndex + 0, ] <- c(df[row, 1], yLow)
        vertices[vIndex + 1, ] <- c(df[row, 1], yHigh)
        vertices[vIndex + 2, ] <- c(df[row, 1], yLabel)
        vertices[vIndex + 3, ] <- c(df[row, 1], yLow)
        vIndex <- vIndex + 4
        
        ## Creating labels and segments
        if(labelBelow) {
          labels[lIndex + 0, ]   <- list(df$Wavelength[row], yLabel - distanceUnit, as.character(df$BVQN[row]), 0.5, 1)
          lIndex <- lIndex + 1
          
          if(row <= nRows - startLabelAt + 1) {
            segments[sIndex + 0, ] <- c(df[row, 1], yHigh + distanceUnit)
            segments[sIndex + 1, ] <- c(df[row, 1], df[row, 2] - distanceUnit)
            sIndex <- sIndex + 2
          }
        } else {
          labels[lIndex + 0, ]   <- list(df$Wavelength[row], yLabel + 1.5*distanceUnit, as.character(df$BVQN[row]), 0.5, 0)
          lIndex <- lIndex + 1
          
          if(row <= nRows - startLabelAt + 1) {
            segments[sIndex + 0, ] <- c(df[row, 1], yHigh - distanceUnit)
            segments[sIndex + 1, ] <- c(df[row, 1], df[row, 2] + distanceUnit)
            sIndex <- sIndex + 2
          }
        }
      } else { ## Unlabeled tick-mark
        vertices[vIndex + 0, ] <- c(df[row, 1], yLow)
        vertices[vIndex + 1, ] <- c(df[row, 1], yHigh)
        vertices[vIndex + 2, ] <- c(df[row, 1], yLow)
        vIndex <- vIndex + 3
      }
    }
    
    vertices[vIndex + 0, ] <- c(df[1, 1], yHigh)
    vertices[vIndex + 1, ] <- c(df[nRows, 1], yHigh)
    vertices[vIndex + 2, ] <- c(df[nRows, 1], yLow)

    return(list(vertices, labels, segments)) ## Need to return within scope or define in outer scope
  }
}
