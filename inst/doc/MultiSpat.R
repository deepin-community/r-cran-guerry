## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  warning = FALSE,
  message = FALSE,   # suppress package loading messages
  comment = "#>",
  fig.height = 3,
  fig.width = 3,
  fig.align = "center"
)

## -----------------------------------------------------------------------------
library(Guerry)       # Guerry's data
library(sp)           # management of spatial data
library(ade4)         # multivariate analysis
library(adegraphics)  # graphical representation
library(spdep)        # spatial dependency
library(adespatial)   # multivariate spatial analysis

## -----------------------------------------------------------------------------
names(gfrance85)

## -----------------------------------------------------------------------------
 data(gfrance85)
 df           <- data.frame(gfrance85)[, 7:12]    # the 6 variables
 france.map   <- as(gfrance85, "SpatialPolygons") # the map
 xy           <- coordinates(gfrance85)           # spatial coordinates
 dep.names    <- data.frame(gfrance85)[, 6]       # departement names
 region.names <- data.frame(gfrance85)[, 5]       # region names
 col.region   <- colors()[c(149, 254, 468, 552, 26)] # colors for region

## -----------------------------------------------------------------------------
pca <- dudi.pca(df, scannf = FALSE, nf = 3)

## -----------------------------------------------------------------------------
biplot(pca, plabel.cex = 0.8)

## -----------------------------------------------------------------------------
pca$eig/sum(pca$eig) * 100

## -----------------------------------------------------------------------------
s.corcircle(pca$co)

## ---- fig.width = 4, fig.height  = 4------------------------------------------
s.label(pca$li, ppoint.col = col.region[region.names], plabel.optim = TRUE, plabel.cex = 0.6)
s.Spatial(france.map, col = col.region[region.names], plabel.cex = 0)
s.class(xy, region.names, col = col.region, add = TRUE, ellipseSize = 0, starSize = 0)

## -----------------------------------------------------------------------------
nb <- poly2nb(gfrance85)
lw <- nb2listw(nb, style = "W")

## -----------------------------------------------------------------------------
s.Spatial(france.map, nb = nb, plabel.cex = 0, pSp.border = "white")

## -----------------------------------------------------------------------------
moran.randtest(df, lw)

## ---- fig.width = 4, fig.height=4---------------------------------------------
x <- df[, 3]
x.lag <- lag.listw(lw, df[, 3])
moran.plot(x, lw)
text(x[5], x.lag[5], dep.names[5], pos = 1, cex = 0.8)

## ---- fig.dim = c(6,3)--------------------------------------------------------
moran.randtest(pca$li, lw)
s.value(xy, pca$li[, 1:2], Sp = france.map, pSp.border = "white", symbol = "circle", pgrid.draw = FALSE)

## -----------------------------------------------------------------------------
 bet <- bca(pca, region.names, scannf = FALSE, nf = 2)

## -----------------------------------------------------------------------------
bet$ratio

## ---- fig.dim=c(5,5)----------------------------------------------------------
plot(bet)

## -----------------------------------------------------------------------------
 barplot(bet$eig)
 bet$eig/sum(bet$eig) * 100

## -----------------------------------------------------------------------------
s.arrow(bet$c1, plabel.cex = 0.8)

## ---- fig.dim = c(4,4)--------------------------------------------------------
s.label(bet$ls, as.character(dep.names), ppoint.cex = 0, plabel.optim = TRUE, plabel.col = col.region[region.names], plabel.cex = 0.5)
s.class(bet$ls, fac = region.names, col = col.region, ellipse = 0, add = TRUE)

## ---- fig.dim = c(6,3)--------------------------------------------------------
s.value(xy, bet$ls, symbol = "circle", Sp = france.map, pSp.col = col.region[region.names], pSp.border = "transparent")

## ---- fig.dim = c(6,4)--------------------------------------------------------
poly.xy <- orthobasis.poly(xy, degree = 2)
s.value(xy, poly.xy, Sp = france.map, plegend.drawKey = FALSE)

## -----------------------------------------------------------------------------
pcaiv.xy <- pcaiv(pca, poly.xy, scannf = FALSE, nf = 2)


## -----------------------------------------------------------------------------
sum(pcaiv.xy$eig)/sum(pca$eig) * 100
pcaiv.xy$eig/sum(pcaiv.xy$eig) * 100

## ---- fig.dim=c(5,5)----------------------------------------------------------
plot(pcaiv.xy)

## ---- fig.dim = c(6,6)--------------------------------------------------------
mem1 <- scores.listw(lw)
s.value(xy, mem1[, 1:9], Sp = france.map, plegend.drawKey = FALSE)


## -----------------------------------------------------------------------------
pcaiv.mem <- pcaiv(pca, mem1[,1:10], scannf = FALSE)

## -----------------------------------------------------------------------------
sum(pcaiv.mem$eig)/sum(pca$eig) * 100
pcaiv.mem$eig/sum(pcaiv.mem$eig) * 100

## ---- fig.dim=c(5,5)----------------------------------------------------------
plot(pcaiv.mem)

## -----------------------------------------------------------------------------
 ms <- multispati(pca, lw, scannf = FALSE)

## ---- fig.dim=c(5,5)----------------------------------------------------------
plot(ms)

## -----------------------------------------------------------------------------
summary(ms)

## -----------------------------------------------------------------------------
s.arrow(ms$c1, plabel.cex = 0.8)

## ---- fig.dim = c(4,4)--------------------------------------------------------
s.match(ms$li, ms$ls, plabel.cex = 0)
s.match(ms$li[c(10, 41, 27), ], ms$ls[c(10, 41, 27), ], label = dep.names[c(10, 
     41, 27)], plabel.cex = 0.8, add = TRUE)

## ---- fig.dim = c(6,3)--------------------------------------------------------
s.value(xy, ms$li, Sp = france.map)

## -----------------------------------------------------------------------------
mat <- matrix(NA, 4, 4)
mat.names <- c("PCA", "BCA", "PCAIV-POLY", "PCAIV-MEM", "MULTISPATI")
colnames(mat) <- mat.names[-5]
rownames(mat) <- mat.names[-1]

mat[1, 1] <- procuste.randtest(pca$li[, 1:2], bet$ls[, 1:2])$obs
mat[2, 1] <- procuste.randtest(pca$li[, 1:2], pcaiv.xy$ls[, 1:2])$obs
mat[3, 1] <- procuste.randtest(pca$li[, 1:2], pcaiv.mem$ls[, 1:2])$obs
mat[4, 1] <- procuste.randtest(pca$li[, 1:2], ms$li[, 1:2])$obs
mat[2, 2] <- procuste.randtest(bet$ls[, 1:2], pcaiv.xy$ls[, 1:2])$obs
mat[3, 2] <- procuste.randtest(bet$ls[, 1:2], pcaiv.mem$ls[, 1:2])$obs
mat[4, 2] <- procuste.randtest(bet$ls[, 1:2], ms$li[, 1:2])$obs
mat[3, 3] <- procuste.randtest(pcaiv.xy$ls[, 1:2], pcaiv.mem$ls[, 1:2])$obs
mat[4, 3] <- procuste.randtest(pcaiv.xy$ls[, 1:2], ms$li[, 1:2])$obs
mat[4, 4] <- procuste.randtest(pcaiv.mem$ls[, 1:2], ms$li[, 1:2])$obs

mat

