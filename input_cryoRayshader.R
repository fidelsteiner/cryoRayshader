################################################################################
# visualizing a dynamic crysophere - code to prepare input for the cryoRayshader
# 
# input_cryoRayshader.R
#
# ReadMe:
# The code enables users to visualize spatial data in the cryosphere in 3D animations, including a model of the glacier as well as graphical data visualization.
# Necessary input data is a DEM of the area, spatial data of the glacier surface changing in time (e.g. data retrieved from satellites or models) and data to visualize on top
# including glacier velocities, surface features like lakes etc.
# 
# Note that distributed glacier thickness data exists for all glaciers globally (consult: https://www.glims.org/RGI/rgi60_dl.html). 
# Velocity datasets also become more and more accesible (see for example https://nsidc.org/data/golive/)
#
# Prerequisites are a working installation of R Studio (https://rstudio.com/) as well as working packages as detailed below.
# 
# Details on path settings as well as input parameters are given below.
#
# Created:          2019/11/21
# Latest Revision:  2019/11/21
#
# Jakob F Steiner| PhD candidate | Faculty of Geosciences | Universiteit Utrecht | Princetonlaan 8a, 3584 CB Utrecht 
# Vening Meinesz building, room 4.30 | P.O. Box 80.115, 3508 TC Utrecht | j.f.steiner@uu.nl | www.uu.nl/staff/jfsteiner | www.mountainhydrology.org 
################################################################################

# clear entire workspace (excl. packages)
rm(list=ls())
gc()


######### 
# loading all necessary libraries, use install.packages('') for initial installation)
######### 
library(fields)
library(maptools)
library(raster)
library(rgeos)
library(tripack)
library(rgdal)
library(rayshader)
library(leaflet)
library(tiff)
library(geoviz)
library(magick)
library(gifski)
library(parallel)
library(snowfall)
library(sf)
library(RColorBrewer)
library(pracma)
######### 

######### 
# specify path of the cryoRayshader function
######### 
source("F:\\PhD\\cryoRayshader\\cryoRayshader.R")
#########

######### 
# specify paths for all input data
######### 
# *** required paths for all model versions
path.DEM <- 'F:\\PhD\\MSc_Projects\\GuidedResearch\\FinalReport\\3DglacierRayshader\\Data\\Langtang\\DEM\\langtang_solarDEM.tif'            # DEM of study region (.tif file)
path.ortho <- 'F:\\PhD\\MSc_Projects\\GuidedResearch\\FinalReport\\3DglacierRayshader\\Data\\Langtang\\Satelliteimages\\3m_resolution.tif'  # stallite imagery for visualization (.tif file; RGB three band image)
path.catchoutline <- 'F:\\PhD\\MSc_Projects\\GuidedResearch\\FinalReport\\3DglacierRayshader\\Data\\Langtang\\Outlines\\catch_proj_1.shp'   # outline of study region (.shp file)
path.glacoutline <- 'F:\\PhD\\MSc_Projects\\GuidedResearch\\FinalReport\\3DglacierRayshader\\Data\\Langtang\\Outlines\\Langtang_2015.shp'   # outline of glacier (.shp file)
path.thickness <- 'F:\\PhD\\MSc_Projects\\GuidedResearch\\FinalReport\\3DglacierRayshader\\Data\\Langtang\\Icethickness'                    # folder that contains all rasters (.tif) with ice thickness data (chronological order)
path.output <- 'F:\\PhD\\MSc_Projects\\GuidedResearch\\cryoRayshader\\allCode\\output'                                                      # folder for final figures

# * paths only required for some model versions
path.velocity <- 'F:\\PhD\\MSc_Projects\\GuidedResearch\\Data\\Langtang\\Velocity'                                                          # folder that contains all rasters (.tif) with glacier velocities (chronological order)
path.features <- 'F:\\PhD\\MSc_Projects\\GuidedResearch\\FinalReport\\3DglacierRayshader\\Data\\Langtang\\Features\\ponds'                  # folder that contains all features (e.g. lakes) as .shp files (chronological order)

projecIn <- '+proj=utm +zone=45N +datum=WGS84'        # projection (UTM, default: '+proj=utm +zone=45N +datum=WGS84')

freeCores <- 0                                        # in case model is performed on multiple cores, choose how many of your computer's cores you want to keep open

######### 
# specify 3D Visualization Parameters (also see plot_3D() in https://cran.r-project.org/web/packages/rayshader/rayshader.pdf)
######### 
zscale <<- 30   # ratio between x/y and z spacing. Starting value the resolution of the DEM raster
fov <<- 0       # field of view angle 
theta <<- 340   # rotation around z axis
zoom3D <<- 0.8    # increase to zoom out (max = 1), decrease to zoom in
colemph <<- 1   # set 1 to emphasize area of mass loss with grey colour/shaded

######### 
# choose model version and set model details
######### 
# currently possible models are 
# (a) mass loss of ice ('massLoss'), which only shows the 3D model of the glacier with changing ice mass
# (b) mass loss of ice plus vizualization of data in a graph ('massLossData'), adding a datagraph to the 3D model
# (c) mass loss of ice but with velocities visualized as data ('velocities')
# (d) visualize features on the glacier surface ('features')

# *** required for all model versions
modVersion <- 'velocities'
currentYear <<- 5                                                   # number of glacier thickness data that represents the same year as the catchment DEM (to create DEM without ice)
outputNames <<- c('1974','2006','2009','2010','2014','2015','2016') # output names for the individual timesteps, has to correspond with time steps of input data
glacName <<- 'Langtang'                                             # Glacier Name for image titles
spatialRasterViz <<- 1                                              # set to 1 will project the spatial data to the glacier surface (currently only works for 'velocities')
elevBands <<- c(4500,4800,5000,5350)                                # Choose cutoff elevations for elevation bands, necessary for 'massLossData' and 'velocities' (at least 1 number, has to lie within the glacier elevation range)
finalFileName <<- 'meanAreaPonds'

# *** required only for 'features' model
feature_Viz <- 'meanarea'                                           # visualize area of individual features ('meanarea'), total area of all features ('totalarea') or just number of features ('count')

#########
# Run Model
#########
cryoRayshader(path,projecIn,modVersion)