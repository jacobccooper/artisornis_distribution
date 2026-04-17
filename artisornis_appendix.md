Artisornis–Modeling Pipeline
================
Cordeiro, Norbert J., Jacob C. Cooper, Victor J. Mkongewa, Martin Joho,
Jasson John & P. Kariuki Ndang’ang’a
2026-04-17

# Introduction

This is the modeling pipeline for *Artisornis moreaui* niche modeling in
East Africa. Please see the manuscript for detailed information on the
data.

## Packages

*Please note I leave a lot of information here from the original
document by Cobos et al. Please see the package page to learn more.*

This code will also have to be adjusted for your own local `filepath` to
work.

``` r
# to create concave polygons
library(concaveman)

# to create better plots
library(cowplot)

# to better process data
library(data.table)

# to determine geographic distances
library(geosphere)

# to perform niche modeling
library(kuenm2)

# to manipulate spatial data and show countries etc.
library(terra)
```

    ## terra 1.9.11

    ## 
    ## Attaching package: 'terra'

    ## The following object is masked from 'package:data.table':
    ## 
    ##     shift

``` r
library(rnaturalearth)
library(sf)
```

    ## Linking to GEOS 3.13.0, GDAL 3.8.5, PROJ 9.5.1; sf_use_s2() is TRUE

``` r
library(tidyterra)
```

    ## 
    ## Attaching package: 'tidyterra'

    ## The following object is masked from 'package:stats':
    ## 
    ##     filter

``` r
# general data manipulation and plotting
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.2.1     ✔ readr     2.2.0
    ## ✔ forcats   1.0.1     ✔ stringr   1.6.0
    ## ✔ ggplot2   4.0.2     ✔ tibble    3.3.1
    ## ✔ lubridate 1.9.5     ✔ tidyr     1.3.2
    ## ✔ purrr     1.2.1

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::between()          masks data.table::between()
    ## ✖ tidyr::extract()          masks terra::extract()
    ## ✖ dplyr::filter()           masks tidyterra::filter(), stats::filter()
    ## ✖ dplyr::first()            masks data.table::first()
    ## ✖ lubridate::hour()         masks data.table::hour()
    ## ✖ lubridate::isoweek()      masks data.table::isoweek()
    ## ✖ lubridate::isoyear()      masks data.table::isoyear()
    ## ✖ dplyr::lag()              masks stats::lag()
    ## ✖ dplyr::last()             masks data.table::last()
    ## ✖ lubridate::mday()         masks data.table::mday()
    ## ✖ lubridate::minute()       masks data.table::minute()
    ## ✖ lubridate::month()        masks data.table::month()
    ## ✖ lubridate::quarter()      masks data.table::quarter()
    ## ✖ ggplot2::remove_missing() masks kuenm2::remove_missing()
    ## ✖ lubridate::second()       masks data.table::second()
    ## ✖ lubridate::stamp()        masks cowplot::stamp()
    ## ✖ purrr::transpose()        masks data.table::transpose()
    ## ✖ lubridate::wday()         masks data.table::wday()
    ## ✖ lubridate::week()         masks data.table::week()
    ## ✖ lubridate::yday()         masks data.table::yday()
    ## ✖ lubridate::year()         masks data.table::year()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

# Establishing the species’ geographic area

## Training area

Given that this species is restricted to a single mountain range, we
opten to create a rough outline of the range following low-lying
elevations around the mountain range by hand. We did not use a more
complex method–such as `grinnell`–as we know the species does not have
the ability to regularly disperse outside of the mountain range.

## Minimum concave polygons

Given that we have really detailed information on this species’
distribution, we created minimum concave polygons to encompass the
geographic areas in which this species occurs.

``` r
# read occurrence shapefile
# consists of all points used in analyses

occ <- read_sf(paste0(filepath, "occ_joint.gpkg"))

nrow(occ)
```

    ## [1] 589

``` r
# isolate localities

occ_locs <- occ |>
  select(ScientificName,Longitude,Latitude,Population) |>
  unique()

# count of localities
nrow(occ_locs)
```

    ## [1] 265

``` r
# get information on distances between points
# limit to distances under 1 km

x <- occ$Longitude
y <- occ$Latitude

xy <- cbind(x,y) |> 
  as.matrix()

dists <- distm(xy) |> 
  as.numeric()

length(dists[dists < 1000 & dists > 0])
```

    ## [1] 22934

``` r
summary(dists[dists < 1000 & dists > 0])
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##     0.0   389.1   607.9   588.7   829.0   998.8

``` r
# create concave shapefiles around points
poly_pop <- function(shp,n){
  shp |>
    filter(Population == n) |>
    concaveman(concavity = 0.05,length_threshold = 0.01)
}

# create the polygons using information from clusters
# cluster information derived from visual QGIS assessment
poly1 <- poly_pop(occ,1)
poly2 <- poly_pop(occ,2)
poly3 <- poly_pop(occ,3)
poly4 <- poly_pop(occ,4)
poly5 <- poly_pop(occ,5)
poly6 <- poly_pop(occ,6)
poly7 <- poly_pop(occ,7)
poly8 <- poly_pop(occ,8)
```

``` r
# bind polygons together
polyx <- rbind(poly1,poly2,poly3,
               poly4,poly5,poly6,
               poly7,poly8)

# many points ca. 130 m apart, so buffering to 100
# understory nonmigratory bird
# points pretty precise
buff_poly <- st_buffer(polyx,dist = 100)
```

``` r
# get occurrence points
pts <- cbind(occ$Longitude,occ$Latitude)

# not colorblind friendly
# not for publication
# view different polygons
plot(pts, pch=19, col="grey",
     asp=1, xlab = "Longitude", ylab = "Latitude")
plot(buff_poly$polygons[1],add=T,col="red")
plot(buff_poly$polygons[2],add=T,col="blue")
plot(buff_poly$polygons[3],add=T,col="green")
plot(buff_poly$polygons[4],add=T,col="black")
plot(buff_poly$polygons[5],add=T,col="gold")
plot(buff_poly$polygons[6],add=T,col="purple")
plot(buff_poly$polygons[7],add=T,col="pink")
plot(buff_poly$polygons[8],add=T,col="lightblue")
```

![](artisornis_appendix_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

Here, we have a plot of all the individual clusters, as derived from
evaluating the distribution of points in QGIS.

``` r
# write shapefile for future use
st_write(buff_poly, dsn = paste0(filepath, "buffered_polygons.shp"), append = F)
```

## Loading BirdLife Shapefile

``` r
# Load BirdLife shapefile
birdlife <- vect(paste0(filepath, 
                        "arti_birdlife/SppDataRequest.shp"))

# select ONLY area of confidence
birdlife_confidence <- birdlife[1,]

# plot area of confidence
plot(birdlife_confidence[1,])
```

![](artisornis_appendix_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

This will be used for later evaluations.

## Calculating number of points per area

``` r
occ # shapefile with population designation
```

    ## Simple feature collection with 589 features and 4 fields
    ## Geometry type: POINT
    ## Dimension:     XY
    ## Bounding box:  xmin: 38.57426 ymin: -5.177928 xmax: 38.68242 ymax: -4.897577
    ## Geodetic CRS:  WGS 84
    ## # A tibble: 589 × 5
    ##    ScientificName     Longitude Latitude Population                 geom
    ##    <chr>                  <dbl>    <dbl>      <int>          <POINT [°]>
    ##  1 Artisornis_moreaui      38.6    -5.09          3 (38.59999 -5.093741)
    ##  2 Artisornis_moreaui      38.6    -5.10          3  (38.6277 -5.099239)
    ##  3 Artisornis_moreaui      38.6    -5.10          3 (38.62279 -5.097138)
    ##  4 Artisornis_moreaui      38.6    -5.10          3 (38.60003 -5.095966)
    ##  5 Artisornis_moreaui      38.6    -5.10          3 (38.59536 -5.102947)
    ##  6 Artisornis_moreaui      38.6    -5.10          3 (38.60023 -5.096093)
    ##  7 Artisornis_moreaui      38.6    -5.10          3 (38.60095 -5.096093)
    ##  8 Artisornis_moreaui      38.6    -5.10          3  (38.5959 -5.098026)
    ##  9 Artisornis_moreaui      38.6    -5.14          5 (38.62527 -5.135369)
    ## 10 Artisornis_moreaui      38.6    -5.14          5  (38.62532 -5.13593)
    ## # ℹ 579 more rows

``` r
# convert to data frame
occ_dat <- occ |> 
  as.data.frame()

# import CSV of converted coordinates
dat <- read_csv(paste0(filepath, "converted_coords.csv"))
```

    ## Rows: 589 Columns: 6
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## dbl  (5): No_of_inds, UTM_North, UTM_East, Longitude, Latitude
    ## dttm (1): Date
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# join datasets together
by = join_by(Longitude, Latitude)

occ_dat_full <- left_join(dat, occ_dat, by = by, 
                          relationship = "many-to-many")

# plot data
# not colorblind friendly, not for publication
plot(occ_dat_full$Longitude, 
     occ_dat_full$Latitude, 
     col = as.factor(occ_dat_full$Population), 
     asp = 1, pch = 19)
```

![](artisornis_appendix_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
occ_dat_full$year <- year(occ_dat_full$Date)

# get a table of sightings per population
summ_tab <- table(occ_dat_full$Population, occ_dat_full$year) |> 
  # reformat data
  unclass() |> 
  as.data.frame()

# converting data to plot it
plot_tab <- table(occ_dat_full$Population, occ_dat_full$year) |> 
  #unclass() |> 
  as.data.frame() |> 
  rename(Population = Var1, Year = Var2) |> 
  mutate(Year = as.numeric(Year))

# getting the sums
summ_tab$pop <- rownames(summ_tab)

# converting the summation table to a new format
summ_tab |> pivot_longer(cols = !pop, 
                         names_to = "Year", 
                         values_to = "Population") |> 
  mutate(Year = as.numeric(Year)) ->
  ggplot_dat
```

``` r
# create plot for points per population per year
a <- ggplot(data = ggplot_dat, aes(x = Year, 
                              y = Population, 
                              col = pop)) +
  # adjust linewidth
  geom_line(linewidth = 2) +
  # change to viridis colors
  scale_colour_viridis_d() +
  # change y label
  ylab("Observations") +
  # change x scale
  scale_x_continuous(name ="Year", breaks = (2006:2017)) + 
  # change theme
  theme_classic() +
  # change orientation
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # change labels
  labs(color = "Population")

plot(a)
```

![](artisornis_appendix_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

Plot of the number of points per population per year. NOTE that only a
few populations account for the vast majority of records.

    ## quartz_off_screen 
    ##                 2

``` r
map <- ggplot() +
  geom_sf(data = buff_poly, aes(fill = FID)) +
  scale_fill_viridis_c(guide = "none") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

grid_plot <- plot_grid(a, NULL, map, ncol = 3, rel_widths = c(1, -0.25, 1))

grid_plot
```

![](artisornis_appendix_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

Here we have the buffered polygons around all the points that were
shown. We buffer to 100 m to reflect the general area, but keep it
restricted because these birds are territorial and non-migratory.

``` r
# summarize counts per population
ggplot_dat |> 
  group_by(pop) |> 
  dplyr::summarise(pop_count = sum(Population))
```

    ## # A tibble: 8 × 2
    ##   pop   pop_count
    ##   <chr>     <int>
    ## 1 1            18
    ## 2 2            42
    ## 3 3          2229
    ## 4 4            68
    ## 5 5          1254
    ## 6 6             8
    ## 7 7            27
    ## 8 8             3

The vast majority of records come from areas 3 and 5.

``` r
# plot the areas with the most points

plot(buff_poly$polygons,
     xlab = "Longitude", ylab = "Latitude", col = "white")
points(pts, pch=19, col="grey")
#plot(buff_poly$polygons[1],add=T,col="#440154")
#plot(buff_poly$polygons[2],add=T,col="#46327e")
plot(buff_poly$polygons[3],
     add=T,
     col="#365c8d")
#plot(buff_poly$polygons[4],add=T,col="#277f8e")
plot(buff_poly$polygons[5],add=T,col="#1fa187")
```

![](artisornis_appendix_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

``` r
#plot(buff_poly$polygons[6],add=T,col="#4ac16d")
#plot(buff_poly$polygons[7],add=T,col="#a0da39")
#plot(buff_poly$polygons[8],add=T,col="#fde725")
```

Here, we are confirming which areas have the most points across all
years.

``` r
buff_poly <- vect(buff_poly)

areas <- expanse(buff_poly, unit = "km")

print(areas)
```

    ## [1]  0.9795936  0.3408333 21.0642505  2.3041742  5.4184320  0.4660753  3.4996980
    ## [8]  0.2828896

These areas range between ca. 21.1 km$^{2}$ for the largest area and
less than a square kilometer for some of the outlier groups. Keep in
mind it says “square meters”, but we are correcting for this by dividing
by 1,000,000. the areas with the vast majority of records are 21.1
km$^2$ and 5.4 km$^2$, the two largest geographic areas.

``` r
sum(areas) |> round(1)
```

    ## [1] 34.4

In total, the distribution of this species is approximately 34.4 km$^2$.

# Area comparisons

## BirdLife

BirdLife lists this area as 328 km$^2$ (area of occupancy) across the
entire plateau [per their
factsheet](https://datazone.birdlife.org/species/factsheet/long-billed-forest-warbler-artisornis-moreaui).

``` r
birdlife |> 
  expanse(unit = "km") |> sum() |> round(1)
```

    ## [1] 214.8

``` r
birdlife_confidence |> 
  expanse(unit = "km") |> round(1)
```

    ## [1] 179.7

Despite this, the shapefile provided shows an area of 214.8 km$^2$
occupied, with the “confident” area of the species being extant being
179.7 km$^2$.

``` r
expanse(buff_poly, unit = "km") |> sum()
```

    ## [1] 34.35595

``` r
buff_poly |> 
  terra::intersect(birdlife_confidence) |> 
  expanse(unit = "km") |> 
  sum()
```

    ## [1] 30.73753

``` r
30.73753/34.35595
```

    ## [1] 0.8946785

89.5% of the distribution is within the birdlife area of confidence.

``` r
birdlife_confidence |> 
  erase(buff_poly) |> 
  expanse(unit = "km") |> 
  sum()
```

    ## [1] 148.9731

Amount of overprediction.

## Comparing to elevation and protected areas

``` r
# load 
elev <- rast(paste0(filepath, "resampled/all_elev.tif"))
tzparks <- vect("~/Dropbox/Manuscripts/Artisornis/ne_tz_parks.gpkg")

# get Usambara cropping area
shp <- vect("~/Dropbox/Manuscripts/Artisornis/cropping_area.shp")
```

``` r
# inactivating lines for saving data

#png("~/Desktop/elev_plot.png")

par(mfrow = c(1,2))

# c(bottom, left, top, right)
# get plotting area
plot(shp, mar = c(2,0,0.5,0))
# add elevation with correct gradient
plot(elev, add = T, col=map.pal("elevation"),
     plg = list(x=38.83, title="Elevation"))
# add localities
points(pts, col = "#0000004C", pch = 20)
# add a scalebar
sbar(xy = "bottomright", 
     col = "black", fill = "black", 
     below = "km", adj = c(0.5, -0.5))

# get plotting area
plot(shp, mar = c(2, 0, 0.5, 0))
# add elevation with correct gradient
plot(elev, add = T, col=map.pal("elevation"), legend = F)
# add protected areas
plot(tzparks, add = T, lwd = 1.5, col = "grey", alpha = 0.5)
# add population polygons
plot(buff_poly, add = T, col = "white", alpha = 0.75)
# add a scalebar
sbar(xy = "bottomright", 
     col = "black", fill = "black", 
     below = "km", adj = c(0.5, -0.5))
```

![](artisornis_appendix_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

``` r
#dev.off()
```

## Comparing areas of occurrence and protected areas etc.

``` r
buff_poly |> 
  terra::intersect(tzparks) ->
  buff_parks

plot(buff_parks)
```

![](artisornis_appendix_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->

``` r
sum(expanse(buff_parks)) / sum(expanse(buff_poly))
```

    ## [1] 0.7669269

77% of the distribution is in protected areas.

``` r
plot(buff_poly, col = "red")
plot(buff_parks, col = "white", add = T)
```

![](artisornis_appendix_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->
Areas of non-overlap include one small polygon and several edges of the
larger polygons.

# Environmental niche modeling

The following are the necessary steps for performing environmental niche
modeling with the KUENM2 pipeline.

## Download and clip environmental layers

First, we need to download environmental layers for these species.
Because these layers are so large, we are going to be cropping them to a
general East African region immediately.

All data are between 2006 and 2017.

``` r
# get a shapefile of eastern Africa
shp <- vect("~/Dropbox/Manuscripts/Artisornis/east_africa.shp")
```

*Note* that I am also loading a shapefile on the land area from
NaturalEarth data to help crop data, and I have a hidden `savepath` for
locally saving my files.

``` r
# system to run bash
# remove white space at end of URL

files <- trimws(read_lines(paste0(filepath,"envidatS3paths.txt")))

for(i in 1:length(files)){
  print(paste("Starting",i))
  ##NOTE##
  # number at the end may vary #
  fname <- strsplit(files[i],"/")[[1]][11]
  
  # define savepath for saving files
  # wget from here not working right...
  # switching to curl
  curl_download(files[i],
                paste0(savepath,fname))
  # crop file
  temp <- crop(terra::rast(paste0(savepath,fname)),
               shp)
  temp <- mask(temp,
               land)
  # if you want to change extension
  # fname <- gsub(".tif",".asc",fname)
  terra::writeRaster(temp,paste0(savepath,fname),
                     overwrite=T)
}
```

We created these models in a custom made polygon encompassing the East
Usambara, as this is the area where the species is restricted and we can
safely say that the species does not regularly disperse outside of this
area in the current climate.

## KU ENM model creation

We are going to run two iterations of KUENM models—one for the data with
*unique points* only, and one with *all points* to reflect the increased
density at places with repeated observations. For each type of model, we
are using the number of background points as close to the same number of
presence points.

### Coarse scale models

``` r
# get Usambara cropping area
shp <- vect("~/Dropbox/Manuscripts/Artisornis/cropping_area.shp")

# get new data layers
env_files <- list.files(paste0(filepath,"resampled/CHELSA"), pattern = "*.tif$", full.names = T)
```

``` r
# resample elevation and canopy to the same as the environmental datalayers
# ensure we have smaller rasters that will work better

# pull up an example raster for standardizing
samp <- rast(env_files[1])

# pull up elevation raster
rast(paste0(filepath, "resampled/all_elev.tif")) |> 
  # crop to study area
  crop(shp) |> 
  # resample to sample raster
  resample(samp) |> 
  # mask to sample area
  mask(shp) ->
  # create object
  elev_rast
  
# write resampled raster
writeRaster(x = elev_rast, filename = paste0(filepath, "resampled/cropped_elev.tif"))

# create a slope raster from elevation raster
elev_rast |> 
  terrain(v = "slope", neighbors = 8) ->
  slope

# save slope raster
writeRaster(x = slope, 
            filename = paste0(filepath,"resampled/slope.tif"),
            NAflag = -9999,
            overwrite = T)

# load and crop canopy raster
rast(paste0(filepath, "resampled/canopy_filled2x.tif")) |> 
  crop(shp) |> 
  resample(samp) |> 
  mask(shp) ->
  canopy_rast

# save canopy raster
writeRaster(x = canopy_rast, filename = paste0(filepath, "resampled/cropped_canopy.tif"))
```

``` r
# load environmental files for modelling
env_files2 <- c(env_files, 
                paste0(filepath, "resampled/cropped_elev.tif"),
                paste0(filepath, "resampled/cropped_canopy.tif"),
                paste0(filepath, "resampled/slope.tif"))

# load environmental variables and ensure all are of the same extent etc.
variables <- terra::rast(env_files2) |> 
  crop(shp) |> 
  mask(shp)

# rename variables for downstream analysis
# leaving names "as is" is a problem for KUENM
variables |> 
  tidyterra::rename(bio10 = `CHELSA_bio10_1981-2010_V.2.1`,
                    bio11 = `CHELSA_bio11_1981-2010_V.2.1`,
                    bio15 = `CHELSA_bio15_1981-2010_V.2.1`,
                    bio6 = `CHELSA_bio6_1981-2010_V.2.1`,
                    npp = `CHELSA_npp_1981-2010_V.2.1`,
                    pet_penman_max = `CHELSA_pet_penman_max_1981-2010_V.2.1`,
                    pen_penman_min = `CHELSA_pet_penman_min_1981-2010_V.2.1`,
                    elev = "Band 1",
                    canopy = canopy_filled2x) ->
  variables

# fix to vect
occ <- vect(occ)
```

Here, we are running the “prepare data” step. Note that most of the
parameters are left at the defaults.

**We are performing two parallel pipelines for unique localities and all
observations.**

``` r
occ_df_unique <- occ |> 
  as.data.frame() |> 
  unique()

occ_df_all <- occ |> 
  as.data.frame()


d_unique <- prepare_data(algorithm = "maxnet",
                  occ = occ_df_unique,
                  x = "Longitude", y = "Latitude",
                  raster_variables = variables,
                  species = "Artisornis moreaui", 
                  partition_method = "kfolds", 
                  n_partitions = 5, 
                  # round up to nearest hundred from number of occurrences
                  # NOTE that it removes some points from background
                  n_background = 300)

d_all <- prepare_data(algorithm = "maxnet",
                  occ = occ_df_all,
                  x = "Longitude", y = "Latitude",
                  raster_variables = variables,
                  species = "Artisornis moreaui", 
                  partition_method = "kfolds", 
                  n_partitions = 5, 
                  # round up to nearest hundred from number of occurrences
                  # NOTE that it removes some points from background
                  n_background = 600)

saveRDS(d_unique, file.path(filepath, "coarse_prepared_data_unique.rds"))
saveRDS(d_all, file.path(filepath, "coarse_prepared_data_all.rds"))
```

We compare these prepared objects and ensure they ran correctly:

    ## prepared_data object summary
    ## ============================
    ## Species: Artisornis moreaui 
    ## Number of Records: 530 
    ##   - Presence: 265 
    ##   - Background: 265 
    ## Partition Method: kfolds 
    ##   - Number of kfolds: 5 
    ## Continuous Variables:
    ##   - bio10, bio11, bio15, bio6, npp, pet_penman_max, pen_penman_min, elev, canopy, slope 
    ## Categorical Variables: None
    ## PCA Information: PCA not performed
    ## Weights: No weights provided
    ## Calibration Parameters:
    ##   - Algorithm: maxnet 
    ##   - Number of candidate models: 10130 
    ##   - Features classes (responses): lq, lqp 
    ##   - Regularization multipliers: 0.1, 0.5, 2, 3, 1

    ## prepared_data object summary
    ## ============================
    ## Species: Artisornis moreaui 
    ## Number of Records: 1111 
    ##   - Presence: 589 
    ##   - Background: 522 
    ## Partition Method: kfolds 
    ##   - Number of kfolds: 5 
    ## Continuous Variables:
    ##   - bio10, bio11, bio15, bio6, npp, pet_penman_max, pen_penman_min, elev, canopy, slope 
    ## Categorical Variables: None
    ## PCA Information: PCA not performed
    ## Weights: No weights provided
    ## Calibration Parameters:
    ##   - Algorithm: maxnet 
    ##   - Number of candidate models: 10130 
    ##   - Features classes (responses): lq, lqp 
    ##   - Regularization multipliers: 0.1, 0.5, 2, 3, 1

We can look at environmental characteristics using the histogram
function:

``` r
calib_hist <- explore_calibration_hist(data = d_all, 
                                       raster_variables = variables,
                                       include_m = T)

plot_calibration_hist(calib_hist)
```

![](artisornis_appendix_files/figure-gfm/unnamed-chunk-40-1.png)<!-- -->

In the above, gray represents the values across the entire calibration
area, blue are the background values, and green are the presence records
(magnified by 2).

We can also look a the locations of presence and absence points, such as
below for the unique dataset:

``` r
pbg <- explore_partition_geo(data = d_unique, raster_variables = variables[[1]])

terra::plot(pbg)
```

![](artisornis_appendix_files/figure-gfm/unnamed-chunk-41-1.png)<!-- -->
\#### Running the models

For each scenario, we are allowing an “error considered” of 2.5% and 1%
to account for high confidence in our data.

``` r
m_maxnet_2.5_unique <- calibration(data = d_unique,
                        error_considered = 2.5,
                        # omission_rate = 5,
                        parallel = T,
                        ncores = 4)

saveRDS(m_maxnet_2.5_unique, file.path(filepath, "candidate_maxnet_25_unique.rds"))

m_maxnet_1_unique <- calibration(data = d_unique,
                        error_considered = 1,
                        # omission_rate = 5,
                        parallel = T,
                        ncores = 4)

saveRDS(m_maxnet_1_unique, file.path(filepath, "candidate_maxnet_1_unique.rds"))

m_maxnet_2.5_all <- calibration(data = d_all,
                        error_considered = 2.5,
                        # omission_rate = 5,
                        parallel = T,
                        ncores = 4)

saveRDS(m_maxnet_2.5_all, file.path(filepath, "candidate_maxnet_25_all.rds"))

m_maxnet_1_all <- calibration(data = d_all,
                        error_considered = 1,
                        # omission_rate = 5,
                        parallel = T,
                        ncores = 4)

saveRDS(m_maxnet_1_all, file.path(filepath, "candidate_maxnet_1_all.rds"))
```

``` r
m_maxnet_2.5_unique <- readRDS(file.path(filepath, "candidate_maxnet_25_unique.rds"))
m_maxnet_1_unique <- readRDS(file.path(filepath, "candidate_maxnet_1_unique.rds"))
m_maxnet_2.5_all <- readRDS(file.path(filepath, "candidate_maxnet_25_all.rds"))
m_maxnet_1_all <- readRDS(file.path(filepath, "candidate_maxnet_1_all.rds"))
```

#### Fitted Models

``` r
fm_2.5_unique <- fit_selected(calibration_results = m_maxnet_2.5_unique,
                              replicate_method = "kfolds", 
                              n_replicates = 5)

saveRDS(fm_2.5_unique, file.path(filepath, "fitted_models_25_unique.rds"))

fm_1_unique <- fit_selected(calibration_results = m_maxnet_1_unique,
                              replicate_method = "kfolds", 
                              n_replicates = 5)

saveRDS(fm_1_unique, file.path(filepath, "fitted_models_1_unique.rds"))

fm_2.5_all <- fit_selected(calibration_results = m_maxnet_2.5_all,
                              replicate_method = "kfolds", 
                              n_replicates = 5)

saveRDS(fm_2.5_all, file.path(filepath, "fitted_models_25_all.rds"))

fm_1_all <- fit_selected(calibration_results = m_maxnet_1_all,
                              replicate_method = "kfolds", 
                              n_replicates = 5)

saveRDS(fm_1_all, file.path(filepath, "fitted_models_1_all.rds"))
```

##### Response curves for 2.5% error

``` r
all_response_curves(fm_2.5_unique, show_variability = T, show_lines = T)
```

![](artisornis_appendix_files/figure-gfm/unnamed-chunk-47-1.png)<!-- -->![](artisornis_appendix_files/figure-gfm/unnamed-chunk-47-2.png)<!-- -->![](artisornis_appendix_files/figure-gfm/unnamed-chunk-47-3.png)<!-- -->![](artisornis_appendix_files/figure-gfm/unnamed-chunk-47-4.png)<!-- -->![](artisornis_appendix_files/figure-gfm/unnamed-chunk-47-5.png)<!-- -->![](artisornis_appendix_files/figure-gfm/unnamed-chunk-47-6.png)<!-- -->

##### Response curves for 1% error

``` r
all_response_curves(fm_1_unique, show_variability = T, show_lines = T)
```

![](artisornis_appendix_files/figure-gfm/unnamed-chunk-48-1.png)<!-- -->![](artisornis_appendix_files/figure-gfm/unnamed-chunk-48-2.png)<!-- -->![](artisornis_appendix_files/figure-gfm/unnamed-chunk-48-3.png)<!-- -->![](artisornis_appendix_files/figure-gfm/unnamed-chunk-48-4.png)<!-- -->![](artisornis_appendix_files/figure-gfm/unnamed-chunk-48-5.png)<!-- -->![](artisornis_appendix_files/figure-gfm/unnamed-chunk-48-6.png)<!-- -->
\#### Variable importance for 2.5% error

``` r
import_plot <- function(x){
  variable_importance(x, progress_bar = F, verbose = F) |> 
    group_by(predictor) |> 
    summarise(contrib_mean = mean(contribution)) ->
    summ_tab
  print(summ_tab)
  variable_importance(x, progress_bar = F, verbose = F) |> 
    plot_importance()
}
import_plot(fm_2.5_unique)
```

    ## # A tibble: 6 × 2
    ##   predictor      contrib_mean
    ##   <chr>                 <dbl>
    ## 1 bio10              -0.00162
    ## 2 bio15               0.0651 
    ## 3 canopy              0.305  
    ## 4 elev                0.397  
    ## 5 pen_penman_min      0.233  
    ## 6 pet_penman_max      0.00504

![](artisornis_appendix_files/figure-gfm/unnamed-chunk-49-1.png)<!-- -->

``` r
import_plot(fm_2.5_all)
```

    ## # A tibble: 6 × 2
    ##   predictor      contrib_mean
    ##   <chr>                 <dbl>
    ## 1 bio15                0.113 
    ## 2 bio6                 0.0172
    ## 3 canopy               0.264 
    ## 4 elev                 0.269 
    ## 5 pet_penman_max       0.128 
    ## 6 slope                0.209

![](artisornis_appendix_files/figure-gfm/unnamed-chunk-50-1.png)<!-- -->

``` r
import_plot(fm_1_unique)
```

    ## # A tibble: 6 × 2
    ##   predictor      contrib_mean
    ##   <chr>                 <dbl>
    ## 1 bio10              -0.00162
    ## 2 bio15               0.0629 
    ## 3 canopy              0.316  
    ## 4 elev                0.392  
    ## 5 pen_penman_min      0.229  
    ## 6 pet_penman_max      0.00504

![](artisornis_appendix_files/figure-gfm/unnamed-chunk-51-1.png)<!-- -->

``` r
import_plot(fm_1_all)
```

    ## # A tibble: 6 × 2
    ##   predictor      contrib_mean
    ##   <chr>                 <dbl>
    ## 1 bio15                0.113 
    ## 2 bio6                 0.0172
    ## 3 canopy               0.264 
    ## 4 elev                 0.269 
    ## 5 pet_penman_max       0.128 
    ## 6 slope                0.209

![](artisornis_appendix_files/figure-gfm/unnamed-chunk-52-1.png)<!-- -->

Elevation, Canopy, and Slope are extremely important for all models.

#### Model Predictions

``` r
p_2.5_unique <- predict_selected(models = fm_2.5_unique, 
                             new_variables = variables, 
                             progress_bar = F)

p_1_unique <- predict_selected(models = fm_1_unique, 
                             new_variables = variables, 
                             progress_bar = F)

p_2.5_all <- predict_selected(models = fm_2.5_all, 
                             new_variables = variables, 
                             progress_bar = F)

p_1_all <- predict_selected(models = fm_1_all, 
                             new_variables = variables, 
                             progress_bar = F)
```

``` r
u_2.5 <- ggplot() +
  geom_spatraster(data = p_2.5_unique$General_consensus$mean, 
                  aes(fill = mean)) +
  scale_fill_viridis_c() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

u_1 <- ggplot() +
  geom_spatraster(data = p_1_unique$General_consensus$mean, 
                  aes(fill = mean)) +
  scale_fill_viridis_c() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

a_2.5 <- ggplot() +
  geom_spatraster(data = p_2.5_all$General_consensus$mean, 
                  aes(fill = mean)) +
  scale_fill_viridis_c() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

a_1 <- ggplot() +
  geom_spatraster(data = p_1_all$General_consensus$mean, 
                  aes(fill = mean)) +
  scale_fill_viridis_c() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

plot_grid(u_2.5, u_1, a_2.5, a_1, ncol = 2)
```

![](artisornis_appendix_files/figure-gfm/unnamed-chunk-54-1.png)<!-- -->

Maps are very similar across all scenarios.

``` r
# show with clamping
p_clamping_local <- predict_selected(models = fm_2.5_unique, 
                               new_variables = variables,
                               consensus = "mean",
                               extrapolation_type = "EC",
                               progress_bar = F)
```

``` r
terra::plot(p_clamping_local$General_consensus$mean, 
            main = "Clamping - 2.5% unique")
```

![](artisornis_appendix_files/figure-gfm/unnamed-chunk-56-1.png)<!-- -->

``` r
p_noextrap_local <- predict_selected(models = fm_2.5_unique, 
                               new_variables = variables,
                               consensus = "mean",
                               extrapolation_type = "NE",
                               progress_bar = F)

terra::plot(p_noextrap_local$General_consensus$mean,
            main = "No Extrapolation - 2.5% unique")
```

![](artisornis_appendix_files/figure-gfm/unnamed-chunk-57-1.png)<!-- -->

``` r
threshold_mean_2.5_unique <- fm_2.5_unique$thresholds$consensus$mean
threshold_mean_1_unique <- fm_1_unique$thresholds$consensus$mean
threshold_mean_2.5_all <- fm_2.5_all$thresholds$consensus$mean
threshold_mean_1_all <- fm_1_all$thresholds$consensus$mean

thresh_2.5_unique <- 
  (p_2.5_unique$General_consensus$mean >= threshold_mean_2.5_unique) * 1

thresh_1_unique <- 
  (p_1_unique$General_consensus$mean >= threshold_mean_1_unique) * 1

thresh_2.5_all <- 
  (p_2.5_all$General_consensus$mean >= threshold_mean_2.5_all) * 1

thresh_1_all <- 
  (p_1_all$General_consensus$mean >= threshold_mean_1_all) * 1
```

``` r
u_2.5 <- ggplot() +
  geom_spatraster(data = thresh_2.5_unique, aes(fill = mean)) +
  ggtitle("2.5% error - unique") +
  scale_fill_viridis_c(na.value = "#FF000000",
                       breaks = c(0, 1),
                       labels = c("Absent", "Present"),
                       guide = guide_legend()) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90),
        plot.margin = unit(c(5.5, 0, 0, 0), "pt")) +
  labs(fill = NULL)

u_1 <- ggplot() +
  geom_spatraster(data = thresh_1_unique, aes(fill = mean)) +
  ggtitle("1% error - unique") +
  scale_fill_viridis_c(na.value = "#FF000000",
                       breaks = c(0, 1),
                       labels = c("Absent", "Present"),
                       guide = "none") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90),
        plot.margin = unit(c(5.5, 0, 0, 0), "pt")) +
  labs(fill = NULL)

all_2.5 <- ggplot() +
  geom_spatraster(data = thresh_2.5_all, aes(fill = mean)) +
  ggtitle("2.5% error - all") +
  scale_fill_viridis_c(na.value = "#FF000000",
                       breaks = c(0, 1),
                       labels = c("Absent", "Present"),
                       guide = guide_legend()) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90),
        plot.margin = unit(c(5.5, 0, 0, 0), "pt")) +
  labs(fill = NULL)

all_1 <- ggplot() +
  geom_spatraster(data = thresh_1_all, aes(fill = mean)) +
  ggtitle("1% error - all") +
  scale_fill_viridis_c(na.value = "#FF000000",
                       breaks = c(0, 1),
                       labels = c("Absent", "Present"),
                       guide = "none") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90),
        plot.margin = unit(c(5.5, 0, 0, 0), "pt")) +
  labs(fill = NULL)

p1 <- plot_grid(u_2.5, NULL, u_1, all_2.5, NULL, all_1, 
                ncol = 3, rel_widths = c(1, -0.6, 1), align = "hv")
print(p1)
```

![](artisornis_appendix_files/figure-gfm/unnamed-chunk-59-1.png)<!-- -->

``` r
save_plot(plot = p1, filename = paste0(filepath, "ENM_threshold.png"))
```

``` r
# convert to polygons

v_2.5_u <- terra::as.polygons(thresh_2.5_unique)[2]
v_1_u <- terra::as.polygons(thresh_1_unique)[2]
v_2.5_a <- terra::as.polygons(thresh_2.5_all)[2]
v_1_a <- terra::as.polygons(thresh_1_all)[2]

v_2.5_u_rast <- p_2.5_unique$General_consensus$mean
```

``` r
expanse(v_2.5_u, unit = "km") |> round(1)
```

    ## [1] 147.3

``` r
expanse(v_2.5_a, unit = "km") |> round(1)
```

    ## [1] 149.9

``` r
expanse(v_1_u, unit = "km") |> round(1)
```

    ## [1] 157.5

``` r
expanse(v_1_a, unit = "km") |> round(1)
```

    ## [1] 167.8

As expected, the 1% error has a larger area than the 2.5% error. Using
all points is also resulting in a larger area than just using the unique
points.

### High resolution data layers

``` r
# changing to 30 m resolution
# 1 m too fine, too much generalization for elevation

highq_elev <- rast(paste0(filepath,"resampled/all_elev.tif")) |> 
  crop(shp) |> 
  mask(shp)

highq_canopy <- rast("~/Downloads/canopy_filled2x.tif") |> 
  resample(highq_elev) |> 
  crop(shp) |> 
  mask(shp)

highq_slope <- highq_elev |> 
  terrain(v = "slope", neighbors = 8)

stacked <- c(highq_canopy, highq_elev, highq_slope) |> 
  rename(canopy = canopy_filled2x,
         elev = `Band 1`)

plot(stacked)
```

![](artisornis_appendix_files/figure-gfm/unnamed-chunk-62-1.png)<!-- -->

``` r
d_unique <- prepare_data(algorithm = "maxnet",
                  occ = occ_df_unique,
                  x = "Longitude", y = "Latitude",
                  raster_variables = stacked,
                  species = "Artisornis moreaui", 
                  partition_method = "kfolds", 
                  n_partitions = 5, 
                  # round up to nearest hundred from number of occurrences
                  # NOTE that it removes some points from background
                  n_background = 300)

d_all <- prepare_data(algorithm = "maxnet",
                  occ = occ_df_all,
                  x = "Longitude", y = "Latitude",
                  raster_variables = stacked,
                  species = "Artisornis moreaui", 
                  partition_method = "kfolds", 
                  n_partitions = 5, 
                  # round up to nearest hundred from number of occurrences
                  # NOTE that it removes some points from background
                  n_background = 600)

saveRDS(d_unique, file.path(filepath, "fine_prepared_data_unique.rds"))
saveRDS(d_all, file.path(filepath, "fine_prepared_data_all.rds"))
```

We compare these prepared objects and ensure they ran correctly:

    ## prepared_data object summary
    ## ============================
    ## Species: Artisornis moreaui 
    ## Number of Records: 563 
    ##   - Presence: 265 
    ##   - Background: 298 
    ## Partition Method: kfolds 
    ##   - Number of kfolds: 5 
    ## Continuous Variables:
    ##   - canopy, elev, slope 
    ## Categorical Variables: None
    ## PCA Information: PCA not performed
    ## Weights: No weights provided
    ## Calibration Parameters:
    ##   - Algorithm: maxnet 
    ##   - Number of candidate models: 40 
    ##   - Features classes (responses): lq, lqp 
    ##   - Regularization multipliers: 0.1, 3, 1, 0.5, 2

    ## prepared_data object summary
    ## ============================
    ## Species: Artisornis moreaui 
    ## Number of Records: 1186 
    ##   - Presence: 589 
    ##   - Background: 597 
    ## Partition Method: kfolds 
    ##   - Number of kfolds: 5 
    ## Continuous Variables:
    ##   - canopy, elev, slope 
    ## Categorical Variables: None
    ## PCA Information: PCA not performed
    ## Weights: No weights provided
    ## Calibration Parameters:
    ##   - Algorithm: maxnet 
    ##   - Number of candidate models: 40 
    ##   - Features classes (responses): lq, lqp 
    ##   - Regularization multipliers: 0.1, 3, 1, 0.5, 2

We can also look a the locations of presence and absence points, such as
below for the unique dataset:

#### Running the models

For each scenario, we are allowing an “error considered” of 2.5% and 1%
to account for high confidence in our data.

``` r
m_maxnet_2.5_unique <- calibration(data = d_unique,
                        error_considered = 2.5,
                        # omission_rate = 5,
                        parallel = T,
                        ncores = 4)

saveRDS(m_maxnet_2.5_unique, file.path(filepath, "fine_candidate_maxnet_25_unique.rds"))

m_maxnet_1_unique <- calibration(data = d_unique,
                        error_considered = 1,
                        # omission_rate = 5,
                        parallel = T,
                        ncores = 4)

saveRDS(m_maxnet_1_unique, file.path(filepath, "fine_candidate_maxnet_1_unique.rds"))

m_maxnet_2.5_all <- calibration(data = d_all,
                        error_considered = 2.5,
                        # omission_rate = 5,
                        parallel = T,
                        ncores = 4)

saveRDS(m_maxnet_2.5_all, file.path(filepath, "fine_candidate_maxnet_25_all.rds"))

m_maxnet_1_all <- calibration(data = d_all,
                        error_considered = 1,
                        # omission_rate = 5,
                        parallel = T,
                        ncores = 4)

saveRDS(m_maxnet_1_all, file.path(filepath, "fine_candidate_maxnet_1_all.rds"))
```

``` r
m_maxnet_2.5_unique <- readRDS(file.path(filepath, "fine_candidate_maxnet_25_unique.rds"))
m_maxnet_1_unique <- readRDS(file.path(filepath, "fine_candidate_maxnet_1_unique.rds"))
m_maxnet_2.5_all <- readRDS(file.path(filepath, "fine_candidate_maxnet_25_all.rds"))
m_maxnet_1_all <- readRDS(file.path(filepath, "fine_candidate_maxnet_1_all.rds"))
```

#### Fitted Models

``` r
fm_2.5_unique <- fit_selected(calibration_results = m_maxnet_2.5_unique,
                              replicate_method = "kfolds", 
                              n_replicates = 5)

saveRDS(fm_2.5_unique, file.path(filepath, "fine_fitted_models_25_unique.rds"))

fm_1_unique <- fit_selected(calibration_results = m_maxnet_1_unique,
                              replicate_method = "kfolds", 
                              n_replicates = 5)

saveRDS(fm_1_unique, file.path(filepath, "fine_fitted_models_1_unique.rds"))

fm_2.5_all <- fit_selected(calibration_results = m_maxnet_2.5_all,
                              replicate_method = "kfolds", 
                              n_replicates = 5)

saveRDS(fm_2.5_all, file.path(filepath, "fine_fitted_models_25_all.rds"))

fm_1_all <- fit_selected(calibration_results = m_maxnet_1_all,
                              replicate_method = "kfolds", 
                              n_replicates = 5)

saveRDS(fm_1_all, file.path(filepath, "fine_fitted_models_1_all.rds"))
```

##### Response curves for 2.5% error

``` r
all_response_curves(fm_2.5_unique, show_variability = T, show_lines = T)
```

![](artisornis_appendix_files/figure-gfm/unnamed-chunk-70-1.png)<!-- -->![](artisornis_appendix_files/figure-gfm/unnamed-chunk-70-2.png)<!-- -->

##### Response curves for 1% error

``` r
all_response_curves(fm_2.5_all, show_variability = T, show_lines = T)
```

![](artisornis_appendix_files/figure-gfm/unnamed-chunk-71-1.png)<!-- -->![](artisornis_appendix_files/figure-gfm/unnamed-chunk-71-2.png)<!-- -->
\#### Variable importance for 2.5% error

``` r
import_plot(fm_2.5_unique)
```

    ## # A tibble: 2 × 2
    ##   predictor contrib_mean
    ##   <chr>            <dbl>
    ## 1 canopy         0.991  
    ## 2 slope          0.00950

![](artisornis_appendix_files/figure-gfm/unnamed-chunk-72-1.png)<!-- -->

``` r
import_plot(fm_2.5_all)
```

    ## # A tibble: 2 × 2
    ##   predictor contrib_mean
    ##   <chr>            <dbl>
    ## 1 canopy          0.0480
    ## 2 elev            0.952

![](artisornis_appendix_files/figure-gfm/unnamed-chunk-73-1.png)<!-- -->

``` r
import_plot(fm_1_unique)
```

    ## # A tibble: 2 × 2
    ##   predictor contrib_mean
    ##   <chr>            <dbl>
    ## 1 canopy         0.991  
    ## 2 slope          0.00950

![](artisornis_appendix_files/figure-gfm/unnamed-chunk-74-1.png)<!-- -->

``` r
import_plot(fm_1_all)
```

    ## # A tibble: 2 × 2
    ##   predictor contrib_mean
    ##   <chr>            <dbl>
    ## 1 canopy          0.0480
    ## 2 elev            0.952

![](artisornis_appendix_files/figure-gfm/unnamed-chunk-75-1.png)<!-- -->

Elevation, Canopy, and Slope are extremely important for all models.

#### Model Predictions

``` r
p_2.5_unique <- predict_selected(models = fm_2.5_unique, 
                             new_variables = variables, 
                             progress_bar = F)

p_1_unique <- predict_selected(models = fm_1_unique, 
                             new_variables = variables, 
                             progress_bar = F)

p_2.5_all <- predict_selected(models = fm_2.5_all, 
                             new_variables = variables, 
                             progress_bar = F)

p_1_all <- predict_selected(models = fm_1_all, 
                             new_variables = variables, 
                             progress_bar = F)
```

``` r
u_2.5 <- ggplot() +
  geom_spatraster(data = p_2.5_unique$General_consensus$mean, 
                  aes(fill = mean)) +
  scale_fill_viridis_c(na.value = "#FF000000") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

u_1 <- ggplot() +
  geom_spatraster(data = p_1_unique$General_consensus$mean, 
                  aes(fill = mean)) +
  scale_fill_viridis_c(na.value = "#FF000000",
                       guide = "none") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

a_2.5 <- ggplot() +
  geom_spatraster(data = p_2.5_all$General_consensus$mean, 
                  aes(fill = mean)) +
  scale_fill_viridis_c(na.value = "#FF000000") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

a_1 <- ggplot() +
  geom_spatraster(data = p_1_all$General_consensus$mean, 
                  aes(fill = mean)) +
  scale_fill_viridis_c(na.value = "#FF000000",
                       guide = "none") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p1 <- plot_grid(u_2.5, NULL, u_1, a_2.5, NULL, a_1, 
                ncol = 3, rel_widths = c(1, -0.65, 1), align = "hv")

p1
```

![](artisornis_appendix_files/figure-gfm/unnamed-chunk-77-1.png)<!-- -->

Maps are very similar across all scenarios.

``` r
# show with clamping
p_clamping_local <- predict_selected(models = fm_2.5_unique, 
                               new_variables = variables,
                               consensus = "mean",
                               extrapolation_type = "EC",
                               progress_bar = F)
```

``` r
terra::plot(p_clamping_local$General_consensus$mean, 
            main = "Clamping - 2.5% unique")
```

![](artisornis_appendix_files/figure-gfm/unnamed-chunk-79-1.png)<!-- -->

``` r
p_noextrap_local <- predict_selected(models = fm_2.5_unique, 
                               new_variables = variables,
                               consensus = "mean",
                               extrapolation_type = "NE",
                               progress_bar = F)

terra::plot(p_noextrap_local$General_consensus$mean,
            main = "No Extrapolation - 2.5% unique")
```

![](artisornis_appendix_files/figure-gfm/unnamed-chunk-80-1.png)<!-- -->

``` r
threshold_mean_2.5_unique <- fm_2.5_unique$thresholds$consensus$mean
threshold_mean_1_unique <- fm_1_unique$thresholds$consensus$mean
threshold_mean_2.5_all <- fm_2.5_all$thresholds$consensus$mean
threshold_mean_1_all <- fm_1_all$thresholds$consensus$mean

thresh_2.5_unique <- 
  (p_2.5_unique$General_consensus$mean >= threshold_mean_2.5_unique) * 1

thresh_1_unique <- 
  (p_1_unique$General_consensus$mean >= threshold_mean_1_unique) * 1

thresh_2.5_all <- 
  (p_2.5_all$General_consensus$mean >= threshold_mean_2.5_all) * 1

thresh_1_all <- 
  (p_1_all$General_consensus$mean >= threshold_mean_1_all) * 1
```

``` r
u_2.5 <- ggplot() +
  geom_spatraster(data = thresh_2.5_unique, aes(fill = mean)) +
  ggtitle("2.5% error - unique") +
  scale_fill_viridis_c(na.value = "#FF000000",
                       breaks = c(0, 1),
                       labels = c("Absent", "Present"),
                       guide = guide_legend()) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90),
        plot.margin = unit(c(5.5, 0, 0, 0), "pt")) +
  labs(fill = NULL)

u_1 <- ggplot() +
  geom_spatraster(data = thresh_1_unique, aes(fill = mean)) +
  ggtitle("1% error - unique") +
  scale_fill_viridis_c(na.value = "#FF000000",
                       breaks = c(0, 1),
                       labels = c("Absent", "Present"),
                       guide = "none") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90),
        plot.margin = unit(c(5.5, 0, 0, 0), "pt")) +
  labs(fill = NULL)

all_2.5 <- ggplot() +
  geom_spatraster(data = thresh_2.5_all, aes(fill = mean)) +
  ggtitle("2.5% error - all") +
  scale_fill_viridis_c(na.value = "#FF000000",
                       breaks = c(0, 1),
                       labels = c("Absent", "Present"),
                       guide = guide_legend()) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90),
        plot.margin = unit(c(5.5, 0, 0, 0), "pt")) +
  labs(fill = NULL)

all_1 <- ggplot() +
  geom_spatraster(data = thresh_1_all, aes(fill = mean)) +
  ggtitle("1% error - all") +
  scale_fill_viridis_c(na.value = "#FF000000",
                       breaks = c(0, 1),
                       labels = c("Absent", "Present"),
                       guide = "none") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90),
        plot.margin = unit(c(5.5, 0, 0, 0), "pt")) +
  labs(fill = NULL)

p1 <- plot_grid(u_2.5, NULL, u_1, all_2.5, NULL, all_1, 
                ncol = 3, rel_widths = c(1, -0.6, 1), align = "hv")
print(p1)
```

![](artisornis_appendix_files/figure-gfm/unnamed-chunk-82-1.png)<!-- -->

``` r
save_plot(plot = p1, filename = paste0(filepath, "coarse_ENM_threshold.png"))
```

``` r
# convert to polygons

v_2.5_u_f <- terra::as.polygons(thresh_2.5_unique)[2]
v_1_u_f <- terra::as.polygons(thresh_1_unique)[2]
v_2.5_a_f <- terra::as.polygons(thresh_2.5_all)[2]
v_1_a_f <- terra::as.polygons(thresh_1_all)[2]

v_1_a_f_rast <- p_1_all$General_consensus$mean
```

``` r
expanse(v_2.5_u_f, unit = "km") |> round(1)
```

    ## [1] 482.8

``` r
expanse(v_2.5_a_f, unit = "km") |> round(1)
```

    ## [1] 206.1

``` r
expanse(v_1_u_f, unit = "km") |> round(1)
```

    ## [1] 553.5

``` r
expanse(v_1_a_f, unit = "km") |> round(1)
```

    ## [1] 224.8

The area of these plots is much larger; note that the area for the
unique plots is much higher than the “all points” plots.

# Comparing Areas of Extent

We are going to do 1% error for unique localities for coarse scale maps
and 1% error for all localities for the fine scale maps.

``` r
# overlap of polygon of range with model
# coarse scale

buff_poly |> 
  # vect() |> 
  terra::intersect(v_1_u) |> 
  expanse(unit = "km") |> 
  sum()
```

    ## [1] 34.27987

``` r
# overlap of polygon of range with model
# fine scale

buff_poly |> 
  # vect() |> 
  terra::intersect(v_1_a_f) |> 
  expanse(unit = "km") |> 
  sum()
```

    ## [1] 32.97362

``` r
# area not accounted for
# coarse
buff_poly |> 
  # vect() |> 
  erase(v_1_u) |> 
  expanse(unit = "km") |> 
  sum()
```

    ## [1] 0.07607422

``` r
# area not accounted for
# fine
buff_poly |> 
  # vect() |> 
  erase(v_1_a_f) |> 
  expanse(unit = "km") |> 
  sum()
```

    ## [1] 1.382327

``` r
# coarse amount covered
34.27987/sum(expanse(buff_poly, unit = "km"))
```

    ## [1] 0.9977856

``` r
# fine amount covered
32.97362/sum(expanse(buff_poly, unit = "km"))
```

    ## [1] 0.9597646

``` r
#coarse overprediction

missing_max <- v_1_u |> 
  erase(buff_poly)

# model area with no birds
v_1_u |> 
  erase(buff_poly) |> 
  expanse(unit = "km") |> 
  sum()
```

    ## [1] 123.2553

``` r
#fine overprediction

missing_max_f <- v_1_a_f |> 
  erase(buff_poly)

# model area with no birds
v_1_a_f |> 
  erase(buff_poly) |> 
  expanse(unit = "km") |> 
  sum()
```

    ## [1] 191.8442

## Comparing model to parks

``` r
buff_poly |> 
  erase(v_1_u) |> 
  # vect() |> 
  makeValid() ->
  poly_out_coarse

poly_out_coarse |> 
  expanse(unit = "km") |> 
  sum()
```

    ## [1] 0.07607422

``` r
buff_poly |> 
  erase(v_1_a_f) |> 
  # vect() |> 
  makeValid() ->
  poly_out_fine

poly_out_fine |> 
  expanse(unit = "km") |> 
  sum()
```

    ## [1] 1.382327

``` r
buff_poly |> 
  erase(poly_out_coarse) |> 
#  vect() |> 
  makeValid()->
  poly_in_coarse

poly_in_coarse |> 
  expanse(unit = "km") |> 
  sum()
```

``` r
buff_poly |> 
  erase(poly_out_fine) |> 
#  vect() |> 
  makeValid()->
  poly_in_fine

poly_in_fine |> 
  expanse(unit = "km") |> 
  sum()
```

    ## [1] 32.97362

``` r
maxnet_inaccurate <- v_1_u |> 
  erase(buff_poly) |> 
#  vect() |> 
  makeValid()

maxnet_accurate <- v_1_u |> 
  erase(maxnet_inaccurate) |> 
  # vect() |> 
  makeValid()
```

``` r
plot(v_1_u)
plot(maxnet_accurate, col = "red", add = T)
plot(maxnet_inaccurate, col = "grey", add = T)
```

## Clusters and protected areas

``` r
tzparks <- vect("~/Dropbox/Manuscripts/Artisornis/ne_tz_parks.gpkg")

plot(tzparks)
```

![](artisornis_appendix_files/figure-gfm/unnamed-chunk-99-1.png)<!-- -->

``` r
buff_poly |> 
  erase(tzparks) ->
  buff_tz

plot(buff_poly)
plot(buff_tz, add = T, col = "red")
```

![](artisornis_appendix_files/figure-gfm/unnamed-chunk-100-1.png)<!-- -->

``` r
1 - sum(expanse(buff_tz, unit = "km"))/sum(expanse(buff_poly, unit = "km"))
```

    ## [1] 0.7669269

## Protected area concordance

``` r
v_1_u |> 
  erase(buff_poly) |> 
  terra::intersect(tzparks) ->
  maxnet_protect_miss

v_1_u |> 
  # erase(buff_poly) |> 
  terra::intersect(tzparks) ->
  maxnet_protect

plot(shp, col = "#FF000000")
plot(v_1_u, col = "black", add = T)
plot(tzparks, add = T, col = "grey")
plot(buff_poly, col = "blue", add = T)
plot(maxnet_protect_miss, add = T, col = "red")
```

![](artisornis_appendix_files/figure-gfm/unnamed-chunk-102-1.png)<!-- -->

``` r
expanse(maxnet_protect, unit = "km") |> 
  sum()
```

    ## [1] 75.18904

``` r
expanse(v_1_u, unit = "km") |> 
  sum() -> a

expanse(maxnet_protect, unit = "km") |> 
  sum() -> b

b/a
```

    ## [1] 0.477284

``` r
v_1_u |> 
  terra::intersect(tzparks) |> 
  terra::intersect(buff_poly) |> 
  expanse(unit = "km") |> 
  sum()
```

    ## [1] 26.30485

``` r
v_1_u |> 
  terra::intersect(tzparks) |> 
  #terra::intersect(buff_poly) |> 
  expanse(unit = "km") |> 
  sum()
```

    ## [1] 75.18904

``` r
buff_poly |> 
  terra::intersect(tzparks) |> 
  expanse(unit = "km") |> 
  sum()
```

    ## [1] 26.3485

``` r
buff_poly |> 
  terra::intersect(tzparks) |> 
  terra::intersect(v_1_u) |> 
  expanse(unit = "km") |> 
  sum()
```

    ## [1] 26.30485

``` r
24.14594/26.30485
```

    ## [1] 0.9179273

``` r
map1 <- ggplot() +
  geom_spatraster(data = v_2.5_u_rast) +
  labs(fill = "Suitability") +
  scale_fill_viridis_c(na.value = "#FF000000", 
                       alpha = 0.6) +
  geom_spatvector(data = v_2.5_u, 
                  fill = "#FF000000",
                  color = "black",
                  size = 1,
                  na.value = "#FF000000") +
  geom_spatvector(data = buff_poly,
                  fill = "white",
                  alpha = 0.75,
                  size = 0.5,
                  color = "black",
                  na.value = "#FF000000") +
  theme_classic() +
  # change orientation
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
```

    ## Warning in layer_sf(geom = GeomSf, data = data, mapping = mapping, stat = stat, : Ignoring unknown parameters: `na.value`
    ## Ignoring unknown parameters: `na.value`

``` r
map2 <- ggplot() +
  geom_spatraster(data = v_1_a_f_rast) +
  labs(fill = "Suitability") +
  scale_fill_viridis_c(na.value = "#FF000000", 
                       alpha = 0.6, guide = "none") +
  geom_spatvector(data = v_1_a_f, 
                  fill = "#FF000000",
                  color = "black",
                  size = 1,
                  na.value = "#FF000000") +
  geom_spatvector(data = buff_poly,
                  fill = "white",
                  alpha = 0.75,
                  size = 0.5,
                  color = "black",
                  na.value = "#FF000000") +
  theme_classic() +
  # change orientation
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
```

    ## Warning in layer_sf(geom = GeomSf, data = data, mapping = mapping, stat = stat, : Ignoring unknown parameters: `na.value`
    ## Ignoring unknown parameters: `na.value`

``` r
fig3 <- plot_grid(map1, NULL, map2, ncol = 3, rel_widths = c(1, -0.3, 1))

fig3
```

![](artisornis_appendix_files/figure-gfm/unnamed-chunk-109-1.png)<!-- -->

``` r
save_plot(filename = "~/OneDrive - University of Nebraska/UNK/Research/Artisornis/figure3.png",
          plot = fig3)
```

End of comparisons; end of pipeline.
