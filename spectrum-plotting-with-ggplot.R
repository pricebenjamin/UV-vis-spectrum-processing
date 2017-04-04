## This file is intended to be used in conjunction with "spectrum-processing.R", as it will
## generate all of the objects that will eventually be plotted.

## Load functions from "spectrum-plotting-functions.R"

head(absorbance, 25)
ggAbsorbance <- absorbance[c(21:nrow(absorbance)), c(1,2)] ## remove undesirable values (spectrometer noise)
names(ggAbsorbance)[2] <- "Absorbance"
head(ggAbsorbance)

head(procData)
ggProcData <- procData
names(ggProcData)[2] <- "Absorbance"
head(ggProcData)

#### Plot using ggplot() ####
require(ggplot2)
require(svglite)

## Parameters for main plot
xlim <- c(500, 625)
xby <- 25
ylim <- c(0, 0.6)
yby <- 0.1
x_scale <- scale_x_continuous(limits = xlim, breaks = seq(xlim[1], xlim[2], by = xby))
y_scale <- scale_y_continuous(limits = ylim, breaks = seq(ylim[1], ylim[2], by = yby))

mainPlot <- ggplot(data = ggAbsorbance) +
  geom_line(mapping = aes(x = Wavelength, y = Absorbance), alpha = 0.6, size = 0.2) +
  x_scale + 
  y_scale +
  labs(title = "UV-Vis spectrum of molecular iodine at T = 60 C", x = "Wavelength [nm]", y = "Absorbance") +
  theme_bw()

mainPlot

## Parameters for xvqn axes and vertical line segments
height <- 0.015
y0 <- 0.05
y1 <- 0.55
y2 <- 0.35
axis0 <- createAxis(ggProcData[which(ggProcData$XVQN == 0), ], y0, y0 + height, y0 - height/3, labelBelow = TRUE,  startLabelAt = 7, labelEvery = 10)
axis1 <- createAxis(ggProcData[which(ggProcData$XVQN == 1), ], y1, y1 - height, y1 + height/3, labelBelow = FALSE, startLabelAt = 6)
axis2 <- createAxis(ggProcData[which(ggProcData$XVQN == 2), ], y2, y2 - height, y2 + height/3, labelBelow = FALSE, startLabelAt = 2)

plotWithAxes <- mainPlot
plotWithAxes <- plotWithAxes + 
  geom_point(
    data = ggProcData[which(ggProcData$XVQN == 0), ], 
    mapping = aes(x = Wavelength, y = Absorbance),
    shape = 1, alpha = 1, size = 1
  ) +
  geom_point(
    data = ggProcData[which(ggProcData$XVQN == 1), ],
    mapping = aes(x = Wavelength, y = Absorbance),
    shape = 2, alpha = 1, size = 1
  ) + 
  geom_point(
    data = ggProcData[which(ggProcData$XVQN == 2), ],
    mapping = aes(x = Wavelength, y = Absorbance),
    shape = 3, alpha = 1, size = 1
  ) + 
  ## Draw v'' = 0, 1, 2 axis
  geom_path(data = data.frame(axis0[1]), mapping = aes(x = AxisX, y = AxisY), size = 0.2) +
  geom_path(data = data.frame(axis1[1]), mapping = aes(x = AxisX, y = AxisY), size = 0.2) + 
  geom_path(data = data.frame(axis2[1]), mapping = aes(x = AxisX, y = AxisY), size = 0.2) +
  ## Draw v'' = 0, 1, 2 axis labels
  geom_text(
    data = data.frame(axis0[2]), ## Labels
    mapping = aes(x = LabelX, y = LabelY, label = Label, hjust = hJust, vjust = vJust),
    size = 2.8
  ) + 
  geom_text(
    data = data.frame(axis1[2]),
    mapping = aes(x = LabelX, y = LabelY, label = Label, hjust = hJust, vjust = vJust),
    size = 2.8
  ) + 
  geom_text(
    data = data.frame(axis2[2]),
    mapping = aes(x = LabelX, y = LabelY, label = Label, hjust = hJust, vjust = vJust),
    size = 2.8
  )

plotWithAxes

## Add vertical line segments showing axis alignment
plotWithSegments <- plotWithAxes
alpha <- 1
size <- 0.2
segments <- as.data.frame(axis0[3])
for(row in seq(1, nrow(segments), by = 2)) {
  plotWithSegments <- plotWithSegments + geom_line(
    data = segments[c(row, row + 1), ], 
    mapping = aes(x = SegmentX, y = SegmentY),
    linetype = "dashed", alpha = alpha, size = size
  )
}
segments <- as.data.frame(axis1[3])
for(row in seq(1, nrow(segments), by = 2)) {
  plotWithSegments <- plotWithSegments + geom_line(
    data = segments[c(row, row + 1), ],
    mapping = aes(x = SegmentX, y = SegmentY),
    linetype = "dashed", alpha = alpha, size = size
  )
}
segments <- as.data.frame(axis2[3])
for(row in seq(1, nrow(segments), by = 2)) {
  plotWithSegments <- plotWithSegments + geom_line(
    data = segments[c(row, row + 1), ],
    mapping = aes(x = SegmentX, y = SegmentY),
    linetype = "dashed", alpha = alpha, size = size
  )
}

plotWithSegments

## End (Plot using ggplot())

## Save ggplot:
filename = "unknown.svg"
width = 8 ## inches
aspectRatio = 0.618 ## golden ratio

ggsave(
  filename,
  device = svg(width = width, height = aspectRatio*width, pointsize = 12)
)
## End (Save ggplot)

#### BEGIN: Plot spectrum using R's plot() function; export as SVG ####
xmin <- 500
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