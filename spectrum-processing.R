## This script expects to find <filename> in the current working directory.
## Additionally, this script uses functions that need to be loaded from "spectrum-processing-functions.R".


#### Load spectral data ####
filename <- "boxcar-combined-spectral-data.csv"
absorbance <- read.csv(filename, header = TRUE)
head(absorbance)
tail(absorbance)

#### Compare different boxcar widths to raw data ####
bw1 <- applyBoxcarWidth(absorbance, 1)
bw2 <- applyBoxcarWidth(absorbance, 2)

x <- 1
y <- 2
plot(
  absorbance[[y]] ~ absorbance[[x]],
  xlab = "Wavelength (nm)", ylab = "Absorbance",
  xlim = c(500, 700), ylim = c(0, 0.8),
  lwd = 0.1, col = "black", type = "l"
)
lines(bw1[[y]] ~ bw1[[x]], lwd = 0.1, col = "orange")
lines(bw2[[y]] ~ bw2[[x]], lwd = 1, col = "red")
abline(v = c(541.2, 539.0, 536.9, 571.6, 568.6, 565.6, 595.7, 592.0, 588.5), col = "blue")

absorbanceWithBoxcar <- applyBoxcarWidth(absorbance, 2) ## Really slow...
head(absorbanceWithBoxcar)
tail(absorbanceWithBoxcar)

#### Spectrum processing ####

x <- 1 ## Column of 'absorbance' dataframe that holds x-values, e.g. 1 <-> Wavelength
y <- 2 ## Column of 'absorbance' dataframe that holds y-values, e.g. 2 <-> Absorbance (at some temperature)

## Filename of processed data that will contain peaks and their assigned vqn
## This filename should explain which absorbance is being processed, 
##  e.g. if y == 2, then you're processing the second column of the loaded csv file.
##  Use head(absorbance) if you don't know which column that is!
outputFilename <- "spectraSuite-bxw-2-processed-peaks.csv"

abs <- absorbance[,c(x,y)] ## Selecting and storing the columns of interest from the absorbance df
##abs <- absorbanceWithBoxcar[,c(x,y)]
## abs <- absorbance[19:3500,c(x,y)]
colNames <- c("Wavelength", "Absorbance")    ## Names of columns
names(abs)[] <- colNames
head(abs)
tail(abs)

## Make sure you've loaded the functions from "spectrum-processing-functions.R"

mxInd <- peakr(abs, xCol = 1, yCol = 2, findMins = FALSE, returnIndices = TRUE) ## find indices of maxima
mnInd <- peakr(abs, xCol = 1, yCol = 2, findMins = TRUE, returnIndices = TRUE) ## find indices of minima
extremeInd <- sort(c(mxInd, mnInd), decreasing = FALSE) ## sort extrema indices

ufmx <- indicesToPoints(abs, mxInd, xCol = 1, yCol = 2) ## Unfiltered maximums (for visual analysis)
head(ufmx)

#### Filter parameters ####
absorbanceDeltaThreshold <- 0.01
removeBefore             <- 513
removeAfter              <- 573
getEveryOtherAfter       <- 548

fmxInd <- maxAbsFilter(abs, extremeInd, yCol = 2, threshold = absorbanceDeltaThreshold) ## filter maxima indices
fmxPts <- indicesToPoints(abs, fmxInd, xCol = 1, yCol = 2) ## create dataframe of filtered maxima
fmxPtsR <- removePointsByValue(fmxPts, value = removeBefore, ltVal = TRUE)
fmxPtsR <- removePointsByValue(fmxPtsR, value = removeAfter, ltVal = FALSE)
fmxPtsR <- removeEveryOtherAfter(fmxPtsR, value = getEveryOtherAfter)

names(fmxPtsR)[] <- c("Wavelength", "Absorbance")
head(fmxPtsR)

## Example of adding specific rows to my filtered max points; use only in special circumstances
## fmxPtsR <- includeRows(fmxPtsR, c(100,148,210), absorbance, x, y)
## head(fmxPtsR)

head(simes)
head(mcnaught)

pointsMatchedToLit <- matchRowToLit(fmxPtsR, simes) ## assign v' values to filtered maxima using literature values
pointsMatchedToLit

procData <- extrapolateV(pointsMatchedToLit) ## assign remaining v' values according to previously assigned lit vals
procData ## fully processed data

## Verify your processed peaks by plotting
point.sym <- 1 ## 19 == solid circular point, 1 == open circle
plot(abs[[2]]~abs[[1]], xlab = names(abs)[1], ylab = names(abs)[2], 
     ##xlim = c(545, 550), ## Zoom in to inspect peaks
     ylim = c(0,0.8), 
     type = 'l')
points(ufmx, col = "red", pch = point.sym)
points(fmxPts, col = "blue", pch = point.sym)
##points(fmxPtsR, col = "blue", pch = point.sym, cex = 1.4)
abline(v = 573) ## Draw some lines to help fined sorting values

## Finding the first peak corresponding to v'' = 2
head(ufmx)
pkwl <- 574.71 ## peak wavelength
tol <- 1 ## tolerance
ufmx[which(abs(ufmx[,1] - pkwl) < tol),] ## Filter for rows whose wavelength is within tol of pkwl
fmxPts <- includeRows(fmxPts, c(1281), absorbance, x, y)

## Plotting with ggplot and labeling fmxPts by wavelength
require(ggplot2)
require(ggrepel)
roundedValues <- fmxPts
roundedValues[,1] <- round(fmxPts[,1], digits = 2)
roundedValues <- removePointsByValue(roundedValues, 505)
roundedValues <- removePointsByValue(roundedValues, 620, ltVal = FALSE)
head(roundedValues)
tail(roundedValues)
ggplot(data = abs, mapping = aes(x = Wavelength, y = Absorbance)) +
  scale_y_continuous(limits = c(0, 0.6)) + 
  geom_line() + 
  geom_point(data = ufmx, col = "red", shape = 1) + 
  geom_point(data = fmxPts, col = "blue", shape = 1) + 
  ggrepel::geom_text_repel(
    data = roundedValues, 
    mapping = aes(label = Wavelength, angle = 90),
    size = 1, alpha = 0.6,
    nudge_x = 0, nudge_y = 0.02,
    force = 0.01,
    segment.color = "lightblue")

## Save ggplot:
ggfilename = "unknown-ggplot.svg"
width = 20 ## inches
aspectRatio = 0.618 ## golden ratio

ggsave(
  ggfilename,
  device = svg(width = width, height = aspectRatio*width, pointsize = 12)
)

## Write out the processed data
write.csv(procData, outputFilename, row.names = FALSE)

#### Export filter parameters and associated input/output filenames + hashes ####
require(digest)
filterParametersFilename = "filter-parameters.csv"
filterParameters <- read.csv(filterParametersFilename, header = TRUE, stringsAsFactors = FALSE)
filterParameters <- rbind(filterParameters,
  c(
    filename,                               ## inputDataFilename
    digest(file = filename, algo = "sha1"), ## inputSHA1
    x, y,
    absorbanceDeltaThreshold,
    removeBefore,
    removeAfter,
    getEveryOtherAfter,
    outputFilename,                               ## outputFile
    digest(file = outputFilename, algo = "sha1")  ## outputSHA1
  )
)
head(filterParameters)
write.csv(filterParameters, filterParametersFilename, row.names = FALSE)

#### Simple R plot ####
## Choose two columns to plot
x <- 3    ## Column number of x-values e.g. (Wavelength) in 'absorbance' dataframe
y <- 4    ## Column number of y-values e.g. (Absorbance) in 'absorbance' dataframe

wd <- 20          ## Width
pt.sz <- 12       ## point size
lwd <- 0.2        ## line width; default = 1
lcolor <- "gray"  ## line color; default = "black"

svgfn <- "unknown.svg" ## Make sure to assign a descriptive filename
svg(filename = svgfn, 
    width = wd, 
    height = wd, 
    pointsize = pt.sz)

plot(absorbance[[y]]~absorbance[[x]], xlab = names(absorbance)[x], ylab = names(absorbance)[y],
     type = 'p', col = lcolor, lwd = lwd
     ## , xaxp = c(500, 700, 40)
     , xlim = c(0,0.2), ylim = c(0,0.2)
     , asp = 1
)
abline(a = 0, b = 1)
abline(v = 547)
abline(v = 580)
points(fmxPts, col = 'red')
points(finalData[,1:2], col = 'green')

## Simes: "Assign upper-state vibrational quantum number v' with the assistance of
## a few values taken from the literature (Table 31-1 and Fig. 31-4)."
abline(v = c(541.2, 539.0, 536.9), col = 'blue') ## Values for I2 from Table 31-1, p. 668
## corresponding vib. quant. wavenum. v' = c(27, 28, 29)

dev.off()