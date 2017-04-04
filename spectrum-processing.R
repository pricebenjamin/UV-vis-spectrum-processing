## This script expects to find <filename> in the current working directory.
## Additionally, this script uses functions that need to be loaded from "spectrum-processing-functions.R".


#### Load spectral data ####
filename <- "boxcar-combined-spectral-data.csv"
absorbance <- read.csv(filename, header = TRUE)
head(absorbance)
tail(absorbance)

#### Spectrum processing ####

x <- 1 ## Column of 'absorbance' dataframe that holds x-values, e.g. 1 <-> Wavelength
y <- 2 ## Column of 'absorbance' dataframe that holds y-values, e.g. 2 <-> Absorbance (at some temperature)

## Filename of processed data that will contain peaks and their assigned vqn
## This filename should explain which absorbance is being processed, 
##  e.g. if y == 2, then you're processing the second column of the loaded csv file.
##  Use head(absorbance) if you don't know which column that is!
outputFilename <- "60_C-spectra-suite-boxcar-width-2-processed-peaks.csv"

## Make sure you've loaded the functions from "spectrum-processing-functions.R"

unfilteredMaxima <- peakr(absorbance, xCol = x, yCol = y, findMins = FALSE, returnIndices = FALSE)

filteredMaxima <- maxAbsFilter(absorbance, xCol = x, yCol = y, threshold = 0.01, type = "left")
filteredMaxima <- removePointsByValue(filteredMaxima, value = 505, ltVal = TRUE)
filteredMaxima <- removePointsByValue(filteredMaxima, value = 612, ltVal = FALSE)

## Crudely filter between ground-state vibrational levels
xvqn0Points <- removePointsByValue(filteredMaxima, value = 574.5, ltVal = FALSE)
xvqn0Points <- removeEveryOtherAfter(xvqn0Points, value = 542)

xvqn1Points <- removePointsByValue(filteredMaxima, value = 545, ltVal = TRUE)
xvqn1Points <- removePointsByValue(xvqn1Points, value = 575, ltVal = FALSE)
xvqn1Points <- removeEveryOtherAfter(xvqn1Points, value = 545)

xvqn2Points <- removePointsByValue(filteredMaxima, value = 581.3, ltVal = TRUE)

## Add remaining points that weren't handled by the crude filtering
xvqn0Points <- rbind(xvqn0Points, unfilteredMaxima[which(abs(unfilteredMaxima$Wavelength - 577.46) < 0.1), ])

xvqn1Points <- rbind(xvqn1Points, unfilteredMaxima[which(abs(unfilteredMaxima$Wavelength - 577.91) < 0.1), ])
xvqn1Points <- rbind(xvqn1Points, unfilteredMaxima[which(abs(unfilteredMaxima$Wavelength - 581.22) < 0.1), ])

xvqn2Points <- rbind(xvqn2Points, unfilteredMaxima[which(abs(unfilteredMaxima$Wavelength - 578.48) < 0.1), ])
xvqn2Points <- rbind(xvqn2Points, unfilteredMaxima[which(abs(unfilteredMaxima$Wavelength - 575.34) < 0.1), ])

## Assign XVQN
xvqn0Points <- data.frame(xvqn0Points, XVQN = 0)
xvqn1Points <- data.frame(xvqn1Points, XVQN = 1)
xvqn2Points <- data.frame(xvqn2Points, XVQN = 2)

## Assign BVQN
xvqn0Points <- matchRowToLit(xvqn0Points, mcnaught[which(mcnaught$XVQN == 0), c(1,3)])
xvqn1Points <- matchRowToLit(xvqn1Points, mcnaught[which(mcnaught$XVQN == 1), c(1,3)])
xvqn2Points <- matchRowToLit(xvqn2Points, mcnaught[which(mcnaught$XVQN == 2), c(1,3)])

## Plot to verify
plot(absorbance[[y]]~absorbance[[x]], type = 'l', xlab = "Wavelength (nm)", ylab = "Absorbance",
     xlim = c(500, 700), ylim = c(0, 0.6))
points(unfilteredMaxima, col = "red")
points(filteredMaxima, col = "green")

blueColors <- colorRampPalette(c("blue", "cyan"))(3)
points(xvqn0Points[, c(1,2)], col = blueColors[1], pch = 19)
points(xvqn0Points[, c(1,2)], col = "black", pch = 1)
points(xvqn1Points[, c(1,2)], col = blueColors[2], pch = 19)
points(xvqn1Points[, c(1,2)], col = "black", pch = 1)
points(xvqn2Points[, c(1,2)], col = blueColors[3], pch = 19)
points(xvqn2Points[, c(1,2)], col = "black", pch = 1)

procData <- rbind(xvqn0Points, xvqn1Points, xvqn2Points)
procData

write.csv(procData, outputFilename, row.names = FALSE)

## Plotting with ggplot and labeling filteredMaxima by wavelength
require(ggplot2)
require(ggrepel)
roundedValues <- filteredMaxima
roundedValues[,1] <- round(filteredMaxima[,1], digits = 2)
head(roundedValues)
tail(roundedValues)
ggplot(data = absorbance, mapping = aes(x = Wavelength, y = bxw2Absorbance)) +
  scale_y_continuous(limits = c(0, 0.6)) + 
  geom_line() + 
  geom_point(data = unfilteredMaxima, col = "red", shape = 1) + 
  geom_point(data = filteredMaxima, col = "blue", shape = 1) + 
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

#### Open SVG file for polotting ####
wd <- 20          ## Width
pt.sz <- 12       ## point size
lwd <- 0.2        ## line width; default = 1
lcolor <- "gray"  ## line color; default = "black"

svgfn <- "unknown.svg" ## Make sure to assign a descriptive filename
svg(filename = svgfn, 
    width = wd, 
    height = wd, 
    pointsize = pt.sz)

## Simes: "Assign upper-state vibrational quantum number v' with the assistance of
## a few values taken from the literature (Table 31-1 and Fig. 31-4)."
abline(v = c(541.2, 539.0, 536.9), col = 'blue') ## Values for I2 from Table 31-1, p. 668
## corresponding vib. quant. wavenum. v' = c(27, 28, 29)

dev.off()
