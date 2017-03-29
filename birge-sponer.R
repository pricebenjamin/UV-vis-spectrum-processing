## Creating Birge-Sponer plot from <filename>

filename <- "highT-processed-peaks-with-vqn-45_5C.csv"
procData <- read.csv(filename, header = TRUE)
head(procData)
names(procData)[3] <- "v'"
head(procData)

## Make sure to include functions defined in data-analysis-functions.R

#### BEGIN: Birge-Sponer ####
birgeSponerData <- birgeSponer(procData)
head(birgeSponerData)
tail(birgeSponerData)

birgeSponerModel <- lm(`deltaNu` ~ `v' + 1`, birgeSponerData)
summary(birgeSponerModel)

## Simple R plots
plot(`deltaNu`~`v' + 1`, birgeSponerData)
abline(coefficients(birgeSponerModel), lty = "dashed")

## Quantile-Quantile plot
lattice::qqmath( ~ resid(birgeSponerModel),
        xlab = "Theoretical Quantiles",
        ylab = "Residuals"
)

## ggplot
require(ggplot2)
require(polynom)

my.eq <- as.character(signif(as.polynomial(coef(birgeSponerModel)), 6))
label.text <- paste(paste("italic(y) == "), gsub("x", "~italic(x)", my.eq, fixed = TRUE))

ggplot(data = birgeSponerData, mapping = aes(x = `v' + 1`, y = `deltaNu`)) + 
  geom_point(shape = 1) + ## Use hollow circles 
  geom_smooth(method = 'lm', se = FALSE, linetype = "dashed") +
  annotate("text", x = 30, y = 87, 
           label = label.text, 
           parse=TRUE) + 
  theme_bw()


#### Export Birge-Sponer coefficients and input filename + hash ####
require(digest)
coefficientsFilename = "birge-sponer-coefficients.csv"
birgeSponCoef <- read.csv(coefficientsFilename, header = TRUE, stringsAsFactors = FALSE)
birgeSponCoef <- rbind(birgeSponCoef,
  c(
    inputDataFilename = filename,
    inputSHA1 = digest(file = filename, algo = "sha1"),
    intercept = coef(birgeSponerModel)[1],
    interceptErr = coef(summary(birgeSponerModel))[1,2],
    firstOrderCoef = coef(birgeSponerModel)[2],
    firstOrderCoefErr = coef(summary(birgeSponerModel))[2,2],
    rSquared = summary(birgeSponerModel)$r.squared
  )
)
head(birgeSponCoef,10)
write.csv(birgeSponCoef, coefficientsFilename, row.names = FALSE)

## Create the first entry and save the file (THIS SHOULD ONLY BE DONE ONCE)
birgeSponCoef <- data.frame(
  inputDataFilename = filename,
  inputSHA1 = digest(file = filename, algo = "sha1"),
  intercept = coef(birgeSponerModel)[1],
  interceptErr = coef(summary(birgeSponerModel))[1,2],
  firstOrderCoef = coef(birgeSponerModel)[2],
  firstOrderCoefErr = coef(summary(birgeSponerModel))[2,2],
  rSquared = summary(birgeSponerModel)$r.squared
)
head(birgeSponCoef)
write.csv(birgeSponCoef, coefficientsFilename, row.names = FALSE)
