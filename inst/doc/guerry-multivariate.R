## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  warning = FALSE,
  message = FALSE,   # suppress package loading messages
  comment = "#>",
  fig.height = 5,
  fig.width = 5,
  dpi = 96
)

# packages to be cited here. Code at the end automatically updates packages.bib
to.cite <- c("car", "effects", "ggplot2", "GGally", "ggbiplot")


## ----load---------------------------------------------------------------------
library(Guerry)         # Guerry data
library(car)            # better scatterplots
library(effects)        # Effect Displays for Linear Models
library(ggplot2)        # Elegant Data Visualisations Using the Grammar of Graphics
library(ggrepel)        # better handling of text labels
library(patchwork)      # combine plots
library(heplots)        # Hypothesis-Error plots
library(candisc)        # Visualizing Generalized Canonical Discriminant Analysis
library(dplyr)          # A Grammar of Data Manipulation
library(tidyr)          # Tidy Messy Data
data(Guerry)

## ----guerry-more-vars---------------------------------------------------------
names(Guerry)[-(1:9)]

## ----guerry-long--------------------------------------------------------------
data("Guerry", package="Guerry")
guerry_long <- Guerry |>
  filter(!is.na(Region)) |>
  select(dept:Suicides) |>
  pivot_longer(cols = Crime_pers:Suicides,
               names_to = "variable",
               values_to = "value")
guerry_long

## ----guerry-density1----------------------------------------------------------
ggplot(data = guerry_long,
       aes(x=value, fill=TRUE)) +
  geom_density(alpha=0.2) +
  geom_rug() +
  facet_wrap(~variable, scales="free") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank())

## ----guerry-density2----------------------------------------------------------
col.region   <- colors()[c(149, 254, 468, 552, 26)] # colors for region
ggplot(data = guerry_long,
       aes(x=value, fill=Region)) +
  geom_density(alpha=0.2) +
  geom_rug() +
  facet_wrap(~variable, scales="free") +
  scale_fill_manual(values=col.region) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom",
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank())


## ----lit-pers-scat0-----------------------------------------------------------
ggplot(aes(x=Literacy, y=Crime_pers/1000), data=Guerry) +
  geom_point(size=2) 


## ----lit-pers-scat------------------------------------------------------------
gdf <- Guerry[, c("Literacy", "Crime_pers", "Department")]
gdf$dsq <- mahalanobis(gdf[,1:2], colMeans(gdf[,1:2]), cov(gdf[,1:2]))

ggplot(aes(x=Literacy, y=Crime_pers/1000, label=Department), data=gdf) +
  geom_point(size=2) +
  stat_ellipse(level=0.68, color="blue", size=1.2) +  
  stat_ellipse(level=0.95, color="gray", size=1, linetype=2) + 
  geom_smooth(method="lm", formula=y~x, fill="lightblue") +
  geom_smooth(method="loess", formula=y~x, color="red", se=FALSE) +
  geom_label_repel(data = gdf[gdf$dsq > 4.6,]) +
  theme_bw()

## ----lit-prop-scat------------------------------------------------------------
gdf <- Guerry[, c("Literacy", "Crime_prop", "Department")]
gdf$dsq <- mahalanobis(gdf[,1:2], colMeans(gdf[,1:2]), cov(gdf[,1:2]))

ggplot(aes(x=Literacy, y=Crime_prop/1000, label=Department), data=gdf) +
  geom_point(size=2) +
  stat_ellipse(level=0.68, color="blue", size=1.2) +  
  stat_ellipse(level=0.95, color="gray", size=1, linetype=2) + 
  geom_smooth(method="lm", formula=y~x, fill="lightblue") +
  geom_smooth(method="loess", formula=y~x, color="red", se=FALSE) +
  geom_label_repel(data = gdf[gdf$dsq > 4.6,]) +
  theme_bw()

## ----spm1---------------------------------------------------------------------
par(mar=rep(2,4))
library(car)          # Companion to Applied Regression
scatterplotMatrix(Guerry[,4:9],
                  ellipse=list(levels=0.68), 
                  smooth=FALSE)

## ----corrgram1----------------------------------------------------------------
par(mar=rep(1,4)+.1)
library(corrgram)             # Plot a Correlogram
corrgram(Guerry[,4:9], upper=panel.pie)

## ----corrgram2----------------------------------------------------------------
par(mar=rep(1,4)+.1)
corrgram(Guerry[,4:9], 
         upper=panel.ellipse, 
         order=TRUE,
         lwd=2)

## ----guerry.pca---------------------------------------------------------------
gdata <- Guerry |>
  select(Region, Crime_pers:Suicides) |>   # keep only main variables
  filter(!is.na(Region))                   # delete Corsica (Region==NA)

guerry.pca <- gdata |>
  select(-Region) |>
  prcomp(scale = TRUE)

print(guerry.pca, digits=3)


## ----biplot1------------------------------------------------------------------
#  if(!require(ggbiplot)) remotes::install_github("vqv/ggbiplot")
#  library(ggbiplot) # A ggplot2 based biplot
#  ggbiplot(guerry.pca, groups=gdata$Region,
#           ellipse=TRUE,
#           var.scale = 3, varname.size = 5) +
#    theme_bw() +
#    labs(color="Region") +
#    theme(legend.position = c(0.1, 0.8))

## ----ggbiplot-----------------------------------------------------------------
knitr::include_graphics("figures/ggbiplot.png")

## ----biplot2, out.width="95%"-------------------------------------------------
op <- par(mar=c(5,4,1,1)+.1)
cols = colorspace::rainbow_hcl(5)
covEllipses(guerry.pca$x, 
            group=gdata$Region, 
            pooled=FALSE, 
            fill=TRUE, fill.alpha=0.1,
            col=cols, 
            label.pos=c(3,0,1,1,3), 
            cex=2,
            xlim=c(-4,4), ylim=c(-4,4),
            xlab = "Dimension 1 (35.7 %)", 
            ylab = "Dimension 2 (20.0 %)",
            cex.lab=1.4
            )
points(guerry.pca$x, pch=(15:19)[Guerry$Region], col=cols[Guerry$Region])

candisc::vectors(guerry.pca$rotation, scale=5,  
                 col="black", lwd=3, cex=1.4, 
                 pos = c(4,2,4,2,2,2),
                 xpd=TRUE)
abline(h=0, v=0, col=gray(.70))

## ----mra-models---------------------------------------------------------------
crime.mod1 <- lm(Crime_pers ~  Region + Literacy + Donations +  Infants + Suicides, data=Guerry)
crime.mod2 <- lm(Crime_prop ~  Region + Literacy + Donations +  Infants + Suicides, data=Guerry)

## ----mra-anova----------------------------------------------------------------
Anova(crime.mod1)
Anova(crime.mod2)

## ----mra-effect1-code---------------------------------------------------------
#  plot(predictorEffects(crime.mod1, ~ Region + Literacy + Infants + Suicides),
#       lwd=2, main="")

## ----mra-effect1--------------------------------------------------------------
knitr::include_graphics("figures/mra-effect1.png")

## ----mra-effect2-code---------------------------------------------------------
#  plot(predictorEffects(crime.mod2, ~ Region + Literacy + Infants + Suicides),
#       lwd=2, main="")

## ----mra-effect2--------------------------------------------------------------
knitr::include_graphics("figures/mra-effect2.png")

## ----manova1------------------------------------------------------------------
crime.mod <- lm(cbind(Crime_pers, Crime_prop) ~ 
                Region + Literacy + Donations +  Infants + Suicides, data=Guerry)
Anova(crime.mod)

## ----cqplot-------------------------------------------------------------------
par(mar=c(3,3,1,1)+.1)
labels <- paste0(Guerry$dept,":", Guerry$Department)
cqplot(crime.mod, id.n=4, labels=labels)

## ----heplot-------------------------------------------------------------------
par(mar = c(3,3,1,1)+.1)
heplot(crime.mod, 
       fill=TRUE, fill.alpha=0.05, 
       cex=1.4, cex.lab=1.3 )

## ----candisc------------------------------------------------------------------
crime.can <- candisc(crime.mod)
crime.can

## ----hecan--------------------------------------------------------------------
par(mar = c(3,3,1,1)+.1)
heplot(crime.can, fill=TRUE, fill.alpha=0.1,
       var.col = "black", 
       var.cex = 1.3,
       cex=1.4, cex.lab=1.3)

## ----write-bib, echo = FALSE--------------------------------------------------
# write a packages.bib file of the packages (.packages()) that have been used here
pkgs <- unique(c(to.cite, .packages()))
knitr::write_bib(pkgs, file = here::here("vignettes", "packages.bib"))

