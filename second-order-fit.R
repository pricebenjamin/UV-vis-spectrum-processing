## Creating Second Order Polynomial Fit for <filename>

filename <- "highT-processed-peaks-with-vqn-45_5C.csv"
procData <- read.csv(filename, header = TRUE)
names(procData)[3] <- "v'"
head(procData)

## Make sure to include functions defined in data-analysis-functions.R

soData <- secondOrder(procData)
head(soData)
tail(soData)

soModel <- lm(Wavenumber ~ poly(`v' + 1/2`, 2, raw = TRUE), soData)
summary(soModel)

## Simple R plot
xVals <- seq(10, 45, by = 0.01)
vFrame <- data.frame(xVals)
names(vFrame)[1] <- "v' + 1/2"

plot(`Wavenumber` ~ `v' + 1/2`, soData)
lines(xVals, predict(soModel, vFrame), lty = "dashed")

#### BEGIN: ggplot ####
library(ggplot2)
library(polyn)

## http://stackoverflow.com/questions/11949331/adding-a-3rd-order-polynomial-and-its-equation-to-a-ggplot-in-r
my.eq <- as.character(signif(as.polynomial(coef(soModel)), 7))
label.text <- paste(paste("italic(y) == "), gsub("x", "~italic(x)", my.eq, fixed = TRUE))

xVals <- seq(14.5, 44.5, by = 0.01)
vFrame <- data.frame(xVals)
names(vFrame)[1] <- "v' + 1/2"
prediction <- data.frame(xVals, predict(soModel, vFrame))
names(prediction)[] <- c("vSamples", "predictedWavenumber")

ggplot(data = soData, mapping = aes(x = `v' + 1/2`, y = Wavenumber)) + 
  geom_point(shape = 1, size = 2) +
  geom_line(data = prediction, 
            mapping = aes(x = vSamples, y = predictedWavenumber), 
            alpha = 0.4, linetype = "dashed") +
  annotate("text", x = 25, y = 19250, 
           label = label.text, 
           parse=TRUE) + 
  theme_bw()


#### Export Second Order fit coefficients and input filename + hash ####
require(digest)
coefficientsFilename = "second-order-coefficients.csv"
soCoef <- read.csv(coefficientsFilename, header = TRUE, stringsAsFactors = FALSE)
soCoef <- rbind(soCoef,
  c(
   inputDataFilename = filename,
   inputSHA1 = digest(file = filename, algo = "sha1"),
   intercept = coef(soModel)[1],
   interceptErr = coef(summary(soModel))[1,2],
   firstOrderCoef = coef(soModel)[2],
   firstOrderCoefErr = coef(summary(soModel))[2,2],
   secondOrderCoef = coef(soModel)[3],
   secondOrderCoefErr = coef(summary(soModel))[3,2],
   rSquared = summary(soModel)$r.squared
  )
)
head(soCoef,10)
write.csv(soCoef, coefficientsFilename, row.names = FALSE)

## Create the first entry and save the file (THIS SHOULD ONLY BE DONE ONCE)
soCoef <- data.frame(
  inputDataFilename = filename,
  inputSHA1 = digest(file = filename, algo = "sha1"),
  intercept = coef(soModel)[1],
  interceptErr = coef(summary(soModel))[1,2],
  firstOrderCoef = coef(soModel)[2],
  firstOrderCoefErr = coef(summary(soModel))[2,2],
  secondOrderCoef = coef(soModel)[3],
  secondOrderCoefErr = coef(summary(soModel))[3,2],
  rSquared = summary(soModel)$r.squared
)
head(soCoef)
write.csv(soCoef, coefficientsFilename, row.names = FALSE)
