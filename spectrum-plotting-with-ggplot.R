## This file is intended to be used in conjunction with "spectrum-processing.R", as it will
## generate all of the objects that will eventually be plotted.


#### BEGIN: Plot spectrum using R's plot() function; export as SVG ####
xmin <- 510
xmax <- 600
ymin <- 0
ymax <- 1
lineType <- 'l'
plotObject <- abs[[y]]~abs[[x]]  ## Defines the x-values and their associated y-values to be plotted
plotColor <- "gray"
pointObject <- fmxPointsR        ## Defines the points to plot on top of the spectrum
pointColor <- "red"
includeSimesAbline <- FALSE

svgFilename <- "unknown.svg"     ## Be sure to use a descriptive name
w <- 20               ## Width
aspectRatio <- 0.618  ## Using the golden ratio (it's genuinely eye-pleasing)
pointSize <- 12
## pointSize defines the thickness (and general size) of certian plot features
##   such as the axes and axis labels. Line thickness of the spectral curve
##   is defined in the plot() function; lwd = 0.4, for example (see ?par).

svg(filename = svgFilename, width = w, height = aspectRatio * w, pointsize = pointSize)

plot(plotObject, lty = lineType, xlim = c(xmin, xmax), ylim = c(ymin, ymax), col = plotColor)
points(pointObject, col = pointColor)

## Simes: "Assign upper-state vibrational quantum number v' with the assistance of
## a few values taken from the literature (Table 31-1 and Fig. 31-4)."
if (includeSimesAbline) abline(v = c(541.2, 539.0, 536.9), col = 'blue') ## Values for I2 from Table 31-1, p. 668
## corresponding vib. quant. wavenum. v' = c(27, 28, 29)

dev.off()
#### END: Plot spectrum using R's plot() function; export as SVG ####

## Plot using ggplot():
install.packages("tidyverse") ## Install tidyverse if necessary
install.packages("svglite")   ## Install svglite to allow us to export ggplot to svg file

library(tidyverse)

## Create axis to show vibrational quantum numbers (vqn)
createAxis <- function(xValues, labelEvery = 5, startLabelAt = 1, yLow = 0.1, yHigh = yLow + 0.0125, yLabel = yLow - 0.005)
{
  nx <- length(xValues) ## Number of x-axis tick marks
  vertices <- 0
  
  if(startLabelAt > labelEvery){
    print("Error: createAxis(): startLabelAt must be less than labelEvery.")  ## To preserve my sanity
    return(NA)
  }
  
  if (nx == 1) {
    print("You don't really need an entire axis for a single point. Have you considered abline()?")
    return(NA)
  } else if (nx <= 0) {
    print("Error: createAxis(): parameter 'nx' must be greater than 1.")
    return(NA)
  } else {  ## (nx > 1)
    nl <- 0 ## Number of labelled tick marks
    
    if( (nx - startLabelAt + 1) %% labelEvery == 0) ## Perfect division
    {
      nl <- (nx - startLabelAt + 1) %/% labelEvery
    } else { ## Imperfect division
      nl <- (nx - startLabelAt + 1) %/% labelEvery + 1
    }

    nv <- 4*nl + 3*(nx - nl)  ## Number of vertices: 3 for each un-labelled tick, 4 for each labelled tick
    r <- (nx - startLabelAt + 1) %% labelEvery ## remainder
    vertices <- matrix(0, nrow = nv, ncol = 2)
    i <- 1
    for(index in c(nx:1)) ## Count down since lowest vqn occurs last in list
    {
      if(index %% labelEvery == r)
      {
        vertices[i+0,] <- c(xValues[index], yLow)
        vertices[i+1,] <- c(xValues[index], yHigh)
        vertices[i+2,] <- c(xValues[index], yLabel)
        vertices[i+3,] <- c(xValues[index], yLow)
        i <- i + 4
      } else {
        vertices[i+0,] <- c(xValues[index], yLow)
        vertices[i+1,] <- c(xValues[index], yHigh)
        vertices[i+2,] <- c(xValues[index], yLow)
        i <- i + 3
      }
    }
    return(as.data.frame(vertices))
  }
}

## Parameters for vqn axis and vertical line segments
peaksOfInterest <- procData ## Renaming for readability (seriously? "filtered max points--some removed"?)
head(peaksOfInterest)
axisTickMarks <- peaksOfInterest[,1]
yLow <- 0.0625
scaleUnit <- 0.0015
vertices <- createAxis(axisTickMarks, labelEvery = 5, startLabelAt = 2,
                       yLow = yLow, yHigh = yLow + 3*scaleUnit, yLabel = yLow - scaleUnit)
nr <- nrow(peaksOfInterest)

## Segment generator:
ySegAbove <- yLow + 4*scaleUnit
ySegBelow <- 3*scaleUnit
genSeg <- function(peaks, index, xCol = 1, yCol = 2, above = ySegAbove, below = ySegBelow)
{
  points <- matrix(0, nrow = 2, ncol = 2)
  points[1,] <- c(peaks[index, xCol], above)
  points[2,] <- c(peaks[index, xCol], peaks[index, yCol] - below)
  return(as.data.frame(points))
}
seg1 <- genSeg(peaksOfInterest, nr-1)
seg2 <- genSeg(peaksOfInterest, nr-6)
seg3 <- genSeg(peaksOfInterest, nr-11)
seg4 <- genSeg(peaksOfInterest, nr-16)
seg5 <- genSeg(peaksOfInterest, nr-21)
seg6 <- genSeg(peaksOfInterest, nr-26)
##segh <- as.data.frame(matrix( c(peaksOfInterest[1,1], peaksOfInterest[nr,1], yLow + 0.0125, yLow + 0.0125), nrow = 2, ncol = 2))

## Parameters for vqn axis labels
yLab <- yLow - 5*scaleUnit
vqnLabels <- as.data.frame(matrix(0, nrow = 7, ncol = 3)) ## 4 total points: axis title + 3 axis labels; each point has x, y, and label
vqnLabels[1,] <- list(peaksOfInterest[1,1] - 2, yLow + 0.0005, "n")
vqnLabels[2,] <- list(peaksOfInterest[nr-1,1], yLab, peaksOfInterest[nr-1,3])
vqnLabels[3,] <- list(peaksOfInterest[nr-6,1], yLab, peaksOfInterest[nr-6,3])
vqnLabels[4,] <- list(peaksOfInterest[nr-11, 1], yLab, peaksOfInterest[nr-11, 3])
vqnLabels[5,] <- list(peaksOfInterest[nr-16, 1], yLab, peaksOfInterest[nr-16, 3])
vqnLabels[6,] <- list(peaksOfInterest[nr-21, 1], yLab, peaksOfInterest[nr-21, 3])
vqnLabels[7,] <- list(peaksOfInterest[nr-26, 1], yLab, peaksOfInterest[nr-26, 3])


## Parameters for final plot
xlim <- c(500, 600)
xby <- 10
ylim <- c(0.05, 0.2)
yby <- 0.025
x_scale <- scale_x_continuous(
  limits = xlim, 
  breaks = seq(xlim[1], xlim[2], by = xby)
)
y_scale <- scale_y_continuous(
  limits = ylim, 
  breaks = seq(ylim[1], ylim[2], by = yby)
)

ggplot(data = absorbance) +
  geom_line(mapping = aes(x = Wavelength, y = roomT5min2), alpha = 0.4) +
  geom_point(
    data = peaksOfInterest, 
    mapping = aes(x = Wavelength, y = Absorbance),
    shape = 1,
    alpha = 1
  ) +
  ## Draw vqn axis
  geom_line(data = vertices, mapping = aes(x = V1, y = V2)) +
  ## Draw vertical lines showing axis alignment
  geom_line(data = seg1, mapping = aes(x = V1, y = V2), linetype = "dashed", alpha = 0.4) +
  geom_line(data = seg2, mapping = aes(x = V1, y = V2), linetype = "dashed", alpha = 0.4) +
  geom_line(data = seg3, mapping = aes(x = V1, y = V2), linetype = "dashed", alpha = 0.4) +
  geom_line(data = seg4, mapping = aes(x = V1, y = V2), linetype = "dashed", alpha = 0.4) +
  geom_line(data = seg5, mapping = aes(x = V1, y = V2), linetype = "dashed", alpha = 0.4) +
  geom_line(data = seg6, mapping = aes(x = V1, y = V2), linetype = "dashed", alpha = 0.4) +
  ##geom_line(data = segh, mapping = aes(x = V1, y = V2)) +
  ## Draw vqn axis labels
  geom_text(data = vqnLabels, mapping = aes(x = V1, y = V2, label = V3, vjust = 0), size = 3.2)+
  x_scale + 
  y_scale +
  labs(
    title = "UV-Vis spectrum of molecular iodine at T = 20 C",
    ## subtitle = "",
    caption = "Circled peaks are due to transitions from v'' = 0 to v' = n where n is given on overlaid axis.",
    x = "Wavelength (nm)", ## use quote() for clean math
    y = "Absorbance"
  ) +
  theme_bw()
## End (Plot using ggplot())

## Save ggplot:
filename = "elec-abs-iodine-roomT-5min2.svg"
width = 8 ## inches
aspectRatio = 0.618 ## golden ratio

ggsave(
  filename,
  device = svg(width = width, height = aspectRatio*width, pointsize = 12)
)
## End (Save ggplot)


## Code for thought:
## Peak finder with threshold parameters
peakIndices <- rep(0,5)
index = 1
tLIndex = 2 ## num indices prior to candiate to determine threshold
tRIndex = 4 ## num indices after candidate to determine threshold
tAbs = 0.01 ## threshold absorbance difference: seems to work well from 0.01 to 0.05 (haven't tested other values)
for (row in c( (1 + tLIndex):(nrow(absorbance) - tRIndex) )) {
  cdate <- row ## Candidate index
  cAbs <- absorbance[row, 2] ## Absorbance at candidate index (col. 2 corresponds to 50.0 C)
  prevAbs <- absorbance[row - 1, 2]
  nextAbs <- absorbance[row + 1, 2]
  farPrevAbs <- absorbance[row - tLIndex, 2]
  farNextAbs <- absorbance[row + tRIndex, 2]
  if (cAbs > prevAbs & cAbs >= nextAbs & (cAbs - farPrevAbs) > tAbs & (cAbs - farNextAbs) > tAbs) {
    ## Found peak from left. Add it to peakIndices and increment index.
    peakIndices[index] <- cdate
    index <- index + 1
  }
}


## Create vqn axis
## Create first tick-mark (a)
vertices[1,] <- c(xValues[1], yHigh)
vertices[2,] <- c(xValues[1], yLow)
## Create midpoints
i <- 3
for(index in c( 2:(length(xValues) - 1) ))
{
  vertices[i+0,] <- c(xValues[index], yLow)
  vertices[i+1,] <- c(xValues[index], yHigh)
  vertices[i+2,] <- c(xValues[index], yLow)
  i <- i + 3
}
## Create final tick-mark
vertices[nv-1,] <- c(xValues[length(xValues)], yLow)
vertices[nv-0,] <- c(xValues[length(xValues)], yHigh)

if( (nx - startLabelAt + 1) %% labelEvery == 0) {
  ## Perfect division requires extra label at end
  nl <- (nx - startLabelAt + 1) %/% labelEvery + 2 ## Include first and last label
} else {
  ## Imperfect division: no extra label required
  
}
## End (Code for thought)