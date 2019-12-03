################################################################################
# cryoRayshader - visualizing cryospheric change in 3D models using the rayshader package
# 
# cryoRayshader.R
#
# ReadMe: 
# see input details in input_cryoRayshader.R as well as on github (https://github.com/fidelsteiner)
#
# Created:          2019/11/21
# Latest Revision:  2019/11/21
#
# Jakob F Steiner| Faculty of Geosciences | Universiteit Utrecht | Princetonlaan 8a, 3584 CB Utrecht 
# Vening Meinesz building, room 4.30 | P.O. Box 80.115, 3508 TC Utrecht | j.f.steiner@uu.nl | www.uu.nl/staff/jfsteiner | www.mountainhydrology.org 
################################################################################

cryoRayshader <- function(path,modVersion){
 
  ###########################
  # load all necessary data
  ###########################
  DEM <- raster(path.DEM)
  DEM <- projectRaster(DEM,crs=projec)
  ortho_image <- brick(path.ortho)
  glacier_outline <- path.glacoutline
  
  # adding the glacier outline
  ogrInfo(glacier_outline)
  glacier_mask <- readOGR(dsn = glacier_outline)
  glacier_mask <- spTransform(glacier_mask, projec) # Transformation of shp files
  
  #cropping the rasters to the glacier extent
  DEM_glac <- crop(DEM,glacier_mask)
  sat_glac <- crop(ortho_image,glacier_mask)

  # Resampling of ortho imagery
  sat_glac_resampled <- raster::resample(sat_glac,DEM_glac, method='bilinear')
  sat_glac_resampled_base <- sat_glac_resampled

  # Model Cores
  
  if(modVersion == 'massLoss'){
    ice_thickness_all <- list.files(path=path.thickness,pattern='.tif',full.names = T)
    if (length(outputNames) != length(ice_thickness_all))
      stop("The number of output file names (outputNames) does not match the number of input files (files in path.thickness)!")
    
    numSteps <- length(ice_thickness_all)
    
    # creating an ice free surface
    latest_thickness <- raster(ice_thickness_all[currentYear])
    latest_thickness <- crop(latest_thickness, glacier_mask)
    latest_thickness <- raster::resample(latest_thickness,DEM_glac, method='bilinear')
    latest_thickness[is.na(latest_thickness)] <- 0
    actual_surface <- DEM_glac - latest_thickness
    
      multiCoreFunc <- function(i){
      glacierthickness <- raster(ice_thickness_all[i])
      glacierthickness <- crop(glacierthickness, glacier_mask)
      glacierthickness <- raster::resample(glacierthickness,DEM_glac, method='bilinear')
      glacierthickness <- focal(glacierthickness, w=matrix(1,nc=5, nr=5),fun=mean)
      glacierthickness[is.na(glacierthickness)] <- 0
      surfaceDEM <- actual_surface + glacierthickness
      
      if(i>1){
        rgl::clear3d()
      }
      
      png(paste(path.output,'\\',outputNames[i],".png",sep=''), width=ncol(DEM_glac), height=nrow(DEM_glac), units = "px", pointsize = 1)# creating a png image with the dimensions of the cropped DEM
      par(mar = c(0,0,0,0), xaxs = "i", yaxs = "i") #Parameters to create a borderless image
      
      sat_glac_resampled <- sat_glac_resampled_base
      
      # emphasize glacier area with colour/grey
      if(colemph == 1){
      sat_glac_resampled[[1]][which(glacierthickness[] > 50)] <- 211
      sat_glac_resampled[[2]][which(glacierthickness[] > 50)] <- 211
      sat_glac_resampled[[3]][which(glacierthickness[] > 50)] <- 211
      }
      
      plotRGB(sat_glac_resampled, stretch="lin")
      dev.off()
      
      terrain_image <- matrix()
      terrain_image <- png::readPNG(paste(path.output,'\\',outputNames[i],".png",sep=''))
      
      #convert the DEM to a matrix
      elmat = matrix(raster::extract(surfaceDEM,raster::extent(surfaceDEM),buffer=1000),
                     nrow = ncol(surfaceDEM),ncol = nrow(surfaceDEM))
      library(rayshader)
      elmat %>% #the 3d function of Rayshader with labels on the map
        sphere_shade(texture = "imhof1") %>%
        add_overlay(terrain_image, alphalayer = 0.9) %>%
        add_shadow(ray_shade(elmat, anglebreaks = seq(60, 90), sunangle = 270, multicore = TRUE, lambert = FALSE, remove_edges = FALSE)) %>%
        add_shadow(ambient_shade(elmat, multicore = TRUE, remove_edges = FALSE)) %>%
        plot_3d(elmat, zscale = zscale, fov = fov,theta = theta, zoom = zoom3D, phi = 35, windowsize = c(nrow(DEM_glac) * 2,ncol(DEM_glac) * 2), solid = T)
      
      render_snapshot(paste(path.output,'\\',outputNames[i],".png",sep='')) #saving the 3d model in png format
    }
    
    # number of workers (leave one to keep system more responsive)
    nworkers <- detectCores() - freeCores
   
      # run entire function parallelized
      sfInit(parallel=T,cpu=nworkers, type="SOCK")                                                # initialize cluster
      loadlib <- lapply(c('raster','truncnorm'),                                                  # load required packages in cluster
                        function(x) sfLibrary(x, character.only=T))
      sfExportAll()                                                                               # transfer all variables to clusters
      sfClusterSetupRNG()                                                                         # set up random number generator
      
      # identify locations of dhdt measurements that are NA or below the accuracy threshold
      
      ModelRunSeries <- seq(1,numSteps,1);
      tout.multi <- system.time(
        out <- sfLapply(ModelRunSeries, multiCoreFunc)
      )
      sfStop() 
    
    #creating GIF with all the produced PNGs
    pngList <- list.files(path=path.output,pattern='.png',full.names = T)
    magick::image_write_gif(magick::image_read(pngList), path = paste(path.output,"\\",finalFileName,".gif",sep = ""), delay = 10/18)
  }
  
  if(modVersion == 'massLossData'){  
    ice_thickness_all <- list.files(path=path.thickness,pattern='.tif',full.names = T)
  if (length(outputNames) != length(ice_thickness_all))
  stop("The number of output file names (outputNames) does not match the number of input files (files in path.thickness)!")
  
  # Create a baseplot twice the size of the original DEM
  DEM_glac_lineplot <- DEM_glac
  
  x2min <-extent(DEM_glac_lineplot)[1]-(extent(DEM_glac_lineplot)[2]-extent(DEM_glac_lineplot)[1])
  newraster <- raster(ext=extent(x2min, extent(DEM_glac_lineplot)[1], extent(DEM_glac_lineplot)[3], extent(DEM_glac_lineplot)[4]), res=c(res(DEM_glac_lineplot)))
  newraster[] <- cellStats(DEM_glac,'min')
  
  # creating an ice free surface
  latest_thickness <- raster(ice_thickness_all[currentYear])
  latest_thickness <- crop(latest_thickness, glacier_mask)
  latest_thickness <- raster::resample(latest_thickness,DEM_glac, method='bilinear')
  latest_thickness[is.na(latest_thickness)] <- 0
  actual_surface <- DEM_glac - latest_thickness
  
  elevCuts <- length(elevBands) - 1
  massLossTimeLine <- matrix(ncol=length(elevBands),nrow=length(ice_thickness_all)) * NA
  colsMat <- c('#000000',brewer.pal(length(elevBands), 'Reds'))
  
  numSteps <- length(ice_thickness_all)
  
    for(i in 1:numSteps){
    thickBrick <- DEM_glac
    sat_glac_resampled <- sat_glac_resampled_base
    png(paste(path.output,'\\',outputNames[i],".png",sep=''), width=ncol(DEM_glac), height=nrow(DEM_glac), units = "px", pointsize = 1)
    par(mar = c(0,0,0,0), xaxs = "i", yaxs = "i") #Parameters to create a borderless image
    
    for(eC in 1: elevCuts){
      thick <- raster(ice_thickness_all[i])
      thick <- crop(thick, glacier_mask)
      thick <- raster::resample(thick,DEM_glac, method='bilinear')
      thick <- focal(thick, w=matrix(1,nc=5, nr=5),fun=mean)
      thick[DEM_glac<=elevBands[eC]|DEM_glac>elevBands[eC + 1]] <- NA
      thick[is.na(thick)] <- 0
      thickBrick <- stack(thickBrick,thick)

    if(i>1){
      rgl::clear3d()
    }


    if(spatialRasterViz != 1){
    shadeCol <- as.vector(col2rgb(colsMat[eC+1]))
    sat_glac_resampled [[1]][which(thick[] > 50)] <- shadeCol[1]
    sat_glac_resampled [[2]][which(thick[] > 50)] <- shadeCol[2]
    sat_glac_resampled [[3]][which(thick[] > 50)] <- shadeCol[3]
    }
    }
    
    thickAllBands <- raster(ice_thickness_all[i])
    thickAllBands <- crop(thickAllBands, glacier_mask)
    thickAllBands <- raster::resample(thickAllBands,DEM_glac, method='bilinear')
    thickAllBands <- focal(thickAllBands, w=matrix(1,nc=5, nr=5),fun=mean)
    thickAllBands2 <- thickAllBands
    thickAllBands[is.na(thickAllBands)] <- 0
    
    
    if(i == 1){
      initialheight <- thickBrick[[2:length(elevBands)]]
      initialheight2 <- initialheight
      initialheight[initialheight<10] <- NA
    }
    actualheight <- thickBrick[[2:length(elevBands)]]
    actualheight2 <- actualheight
    actualheight[actualheight<10] <- NA
    
    dH <- calc(actualheight2 - initialheight2,sum,na.rm=T)
    
    if(spatialRasterViz == 1){
      shadeCol <- as.vector(col2rgb(colsMat[length(elevBands)]))
      sat_glac_resampled[[1]][which(thickAllBands2[] > 50)] <- 255 - abs(dH[which(thickAllBands2[] > 50)]/max(dH[which(thickAllBands2[] > 50)]))*shadeCol[1]
      sat_glac_resampled[[2]][which(thickAllBands2[] > 50)] <- 255 - abs(dH[which(thickAllBands2[] > 50)]/max(dH[which(thickAllBands2[] > 50)]))*shadeCol[2]
      sat_glac_resampled[[3]][which(thickAllBands2[] > 50)] <- 255 - abs(dH[which(thickAllBands2[] > 50)]/max(dH[which(thickAllBands2[] > 50)]))*shadeCol[3]
      
    }
    
    plotRGB(sat_glac_resampled, stretch="lin")
    dev.off()
    
    massLossTimeLine[i,2:4] <- 100 - cellStats(initialheight - actualheight,sum)/cellStats(initialheight,sum)*100
    massLossTimeLine[i,1] <- 100 - sum(cellStats(initialheight - actualheight,sum))/sum(cellStats(initialheight,sum))*100    

    png(paste(path.output,'\\',outputNames[i],"_graph.png",sep=''), width=nrow(DEM_glac), height=ncol(DEM_glac), res = 100)
    par(mar = c(2,4,2,2), xaxs = "i", yaxs = "i") #Parameters to create a borderless image
    matplot(seq(1,length(ice_thickness_all),1),massLossTimeLine,type='b',lty= 1,lwd = 2, pch=18,ylim=c(50,120),col=colsMat,xlab = '',ylab = 'ice thickness [%]', xaxt="n",main=paste(glacName,outputNames[i]))
    axis(side=1, at=seq(1,length(ice_thickness_all),1), labels = outputNames)
    if(spatialRasterViz == 1){
      legend('topright',c('mean',paste('below',elevBands[2:length(elevBands)],'m')),lty= 1,lwd = 2, pch=18,col=colsMat,bty='n')}
    grid(NULL,NULL)
    dev.off()
    
    terrain_image <- matrix()
    terrain_image <- png::readPNG(paste(path.output,'\\',outputNames[i],".png",sep=''))
    terrain_image1 <- png::readPNG(paste(path.output,'\\',outputNames[i],"_graph.png",sep=''))
    
    terrain_image2 <- array(0, c(nrow(terrain_image), ncol(terrain_image), 3))
    terrain_image2[,,1] <- flipud(fliplr(apply(t(terrain_image1[,,1]), 2, rev)))
    terrain_image2[,,2] <- flipud(fliplr(apply(t(terrain_image1[,,2]), 2, rev)))
    terrain_image2[,,3] <- flipud(fliplr(apply(t(terrain_image1[,,3]), 2, rev)))
    
    terrain_imageFINAL <- array(0, c(nrow(terrain_image), ncol(terrain_image)*2, 3))
    terrain_imageFINAL[,seq(1,dim(terrain_image2)[2],1),] <- terrain_image2
    terrain_imageFINAL[,seq(dim(terrain_image2)[2]+1,dim(terrain_image2)[2]*2,1),] <- terrain_image
    terrain_imageFINAL[,seq(dim(terrain_image2)[2]-2,dim(terrain_image2)[2]+2,1),]<-0
    
    surfaceDEM <- actual_surface + thickAllBands
    finalraster<-merge(newraster,surfaceDEM)
    library(rayshader)
    elmat = matrix(raster::extract(finalraster,raster::extent(finalraster),buffer=1000),
                   nrow = ncol(finalraster),ncol = nrow(finalraster))
    elmat %>% 
      sphere_shade(texture = "imhof1") %>%
      add_overlay(terrain_imageFINAL, alphalayer = 0.9) %>%
      add_shadow(ray_shade(elmat, anglebreaks = seq(60, 90), sunangle = 270, multicore = TRUE, lambert = FALSE, remove_edges = FALSE)) %>%
      add_shadow(ambient_shade(elmat, multicore = TRUE, remove_edges = FALSE)) %>%
      plot_3d(elmat, zscale = zscale, fov = fov,theta = theta, zoom = zoom3D, phi = 45, windowsize = c(nrow(DEM_glac) * 2,ncol(DEM_glac) * 2), solid = T)
    
    
    render_snapshot(paste(path.output,'\\',outputNames[i],"_2Figs.png",sep=''))
  }

  pngList <- list.files(path=paste(path.output,"\\",sep=""),pattern='_2Figs.png',full.names = T)
  magick::image_write_gif(magick::image_read(pngList), 
                          path = paste(path.output,"\\",finalFileName,".gif",sep = ""), 
                          delay = 10/18)
  }
  
  if(modVersion == 'velocities'){
    ice_thickness_all <- list.files(path=path.thickness,pattern='.tif',full.names = T)
    vel_all <- list.files(path=path.velocity,pattern='.tif',full.names = T)
    if (length(outputNames) != length(ice_thickness_all))
      stop("The number of output file names (outputNames) does not match the number of input files (files in path.thickness)!")
    if (length(outputNames) != length(vel_all))
      stop("The number of output file names (outputNames) does not match the number of input files (files in path.velocity)!")
    
    # Create a baseplot twice the size of the original DEM
    DEM_glac_lineplot <- DEM_glac
    
    x2min <-extent(DEM_glac_lineplot)[1]-(extent(DEM_glac_lineplot)[2]-extent(DEM_glac_lineplot)[1])
    newraster <- raster(ext=extent(x2min, extent(DEM_glac_lineplot)[1], extent(DEM_glac_lineplot)[3], extent(DEM_glac_lineplot)[4]), res=c(res(DEM_glac_lineplot)))
    newraster[] <- cellStats(DEM_glac,'min')
    
    # creating an ice free surface
    latest_thickness <- raster(ice_thickness_all[currentYear])
    latest_thickness <- crop(latest_thickness, glacier_mask)
    latest_thickness <- raster::resample(latest_thickness,DEM_glac, method='bilinear')
    latest_thickness[is.na(latest_thickness)] <- 0
    actual_surface <- DEM_glac - latest_thickness
    
    elevCuts <- length(elevBands) - 1
    
    colsMat <- c('#000000',brewer.pal(length(elevBands), 'Reds'))
    
    numSteps <- length(ice_thickness_all)
    velTimeLine <- matrix(ncol=length(elevBands),nrow=length(ice_thickness_all)) * NA
    
    for(i in 1:numSteps){
      velBrick <- DEM_glac
      sat_glac_resampled <- sat_glac_resampled_base
      png(paste(path.output,'\\',outputNames[i],".png",sep=''), width=ncol(DEM_glac), height=nrow(DEM_glac), units = "px", pointsize = 1)
      par(mar = c(0,0,0,0), xaxs = "i", yaxs = "i") #Parameters to create a borderless image
      
      for(eC in 1: elevCuts){
        velSet <- raster(vel_all[i])
        velSet <- crop(velSet, glacier_mask)
        velSet <- raster::resample(velSet,DEM_glac, method='bilinear')
        velSet <- focal(velSet, w=matrix(1,nc=5, nr=5),fun=mean)
        velSet[DEM_glac<=elevBands[eC]|DEM_glac>elevBands[eC + 1]] <- NA
        velSet[velSet[]==0] <- NA
        velBrick <- stack(velBrick,velSet)
        
        if(i>1){
          rgl::clear3d()
        }
        if(spatialRasterViz != 1){
        shadeCol <- as.vector(col2rgb(colsMat[eC+1]))
        sat_glac_resampled [[1]][which(velSet[] > 0)] <- shadeCol[1]
        sat_glac_resampled [[2]][which(velSet[] > 0)] <- shadeCol[2]
        sat_glac_resampled [[3]][which(velSet[] > 0)] <- shadeCol[3]
        }
      }
      
      velAllBands <- raster(vel_all[i])
      velAllBands <- crop(velAllBands, glacier_mask)
      velAllBands <- raster::resample(velAllBands,DEM_glac, method='bilinear')
      velAllBands <- focal(velAllBands, w=matrix(1,nc=5, nr=5),fun=mean)
      velAllBands[velAllBands[]==0] <- NA
      if(spatialRasterViz == 1){
      shadeCol <- as.vector(col2rgb(colsMat[length(elevBands)]))
      
      sat_glac_resampled [[1]][which(velAllBands[] > 0)] <- 255-velAllBands[which(velAllBands[] > 0)]/max(velAllBands[which(velAllBands[] > 0)])*shadeCol[1]
      sat_glac_resampled [[2]][which(velAllBands[] > 0)] <- 255-velAllBands[which(velAllBands[] > 0)]/max(velAllBands[which(velAllBands[] > 0)])*shadeCol[2]
      sat_glac_resampled [[3]][which(velAllBands[] > 0)] <- 255-velAllBands[which(velAllBands[] > 0)]/max(velAllBands[which(velAllBands[] > 0)])*shadeCol[3]
      
      }
      
      plotRGB(sat_glac_resampled, stretch="lin")
      dev.off()
      
      thickAllBands <- raster(ice_thickness_all[i])
      thickAllBands <- crop(thickAllBands, glacier_mask)
      thickAllBands <- raster::resample(thickAllBands,DEM_glac, method='bilinear')
      thickAllBands <- focal(thickAllBands, w=matrix(1,nc=5, nr=5),fun=mean)
      thickAllBands[is.na(thickAllBands)] <- 0
    
      velTimeLine[i,2:4] <- cellStats(velBrick[[2:length(elevBands)]],mean)
      velTimeLine[i,1] <- mean(cellStats(velAllBands,mean))    
      
      png(paste(path.output,'\\',outputNames[i],"_graph.png",sep=''), width=nrow(DEM_glac), height=ncol(DEM_glac), res = 100)
      par(mar = c(2,4,2,2), xaxs = "i", yaxs = "i") #Parameters to create a borderless image
      matplot(seq(1,length(vel_all),1),velTimeLine,type='b',lty= 1,lwd = 2, pch=18,ylim=c(0,20),col=colsMat,xlab = '',ylab = 'velocities m/a', xaxt="n",main=paste(glacName,outputNames[i]))
      axis(side=1, at=seq(1,length(vel_all),1), labels = outputNames)
      grid(NULL,NULL)
      if(spatialRasterViz == 1){
        legend('topright',c('mean',paste('below',elevBands[2:length(elevBands)],'m')),lty= 1,lwd = 2, pch=18,col=colsMat,bty='n')}
      dev.off()
      
      terrain_image <- matrix()
      terrain_image <- png::readPNG(paste(path.output,'\\',outputNames[i],".png",sep=''))
      terrain_image1 <- png::readPNG(paste(path.output,'\\',outputNames[i],"_graph.png",sep=''))
      
      terrain_image2 <- array(0, c(nrow(terrain_image), ncol(terrain_image), 3))
      terrain_image2[,,1] <- flipud(fliplr(apply(t(terrain_image1[,,1]), 2, rev)))
      terrain_image2[,,2] <- flipud(fliplr(apply(t(terrain_image1[,,2]), 2, rev)))
      terrain_image2[,,3] <- flipud(fliplr(apply(t(terrain_image1[,,3]), 2, rev)))
      
      terrain_imageFINAL <- array(0, c(nrow(terrain_image), ncol(terrain_image)*2, 3))
      terrain_imageFINAL[,seq(1,dim(terrain_image2)[2],1),] <- terrain_image2
      terrain_imageFINAL[,seq(dim(terrain_image2)[2]+1,dim(terrain_image2)[2]*2,1),] <- terrain_image
      terrain_imageFINAL[,seq(dim(terrain_image2)[2]-2,dim(terrain_image2)[2]+2,1),]<-0
      
      surfaceDEM <- actual_surface + thickAllBands
      finalraster<-merge(newraster,surfaceDEM)

      elmat = matrix(raster::extract(finalraster,raster::extent(finalraster),buffer=1000),
                     nrow = ncol(finalraster),ncol = nrow(finalraster))
      elmat %>% 
        sphere_shade(texture = "imhof1") %>%
        add_overlay(terrain_imageFINAL, alphalayer = 0.9) %>%
        add_shadow(ray_shade(elmat, anglebreaks = seq(60, 90), sunangle = 270, multicore = TRUE, lambert = FALSE, remove_edges = FALSE)) %>%
        add_shadow(ambient_shade(elmat, multicore = TRUE, remove_edges = FALSE)) %>%
        plot_3d(elmat, zscale = zscale, fov = fov,theta = theta, zoom = zoom3D, phi = 45, windowsize = c(nrow(DEM_glac) * 2,ncol(DEM_glac) * 2), solid = T)
      
      render_snapshot(paste(path.output,'\\',outputNames[i],"_2Figs.png",sep=''))
    }
    
    pngList <- list.files(path=paste(path.output,"\\",sep=""),pattern='_2Figs.png',full.names = T)
    magick::image_write_gif(magick::image_read(pngList[2:length(pngList)]), 
                            path = paste(path.output,"\\",finalFileName,".gif",sep = ""), 
                            delay = 10/18)
  }
  
  if(modVersion == 'features'){
    ice_thickness_all <- list.files(path=path.thickness,pattern='.tif',full.names = T)
    feature_all <- list.files(path=path.features,pattern='.shp',full.names = T)
    if (length(outputNames) != length(ice_thickness_all))
      stop("The number of output file names (outputNames) does not match the number of input files (files in path.thickness)!")
    if (length(outputNames) != length(feature_all))
      stop("The number of output file names (outputNames) does not match the number of input files (files in path.features)!")
    
    # Create a baseplot twice the size of the original DEM
    DEM_glac_lineplot <- DEM_glac
    
    x2min <-extent(DEM_glac_lineplot)[1]-(extent(DEM_glac_lineplot)[2]-extent(DEM_glac_lineplot)[1])
    newraster <- raster(ext=extent(x2min, extent(DEM_glac_lineplot)[1], extent(DEM_glac_lineplot)[3], extent(DEM_glac_lineplot)[4]), res=c(res(DEM_glac_lineplot)))
    newraster[] <- cellStats(DEM_glac,'min')
    
    # creating an ice free surface
    latest_thickness <- raster(ice_thickness_all[currentYear])
    latest_thickness <- crop(latest_thickness, glacier_mask)
    latest_thickness <- resample(latest_thickness,DEM_glac, method='bilinear')
    latest_thickness[is.na(latest_thickness)] <- 0
    actual_surface <- DEM_glac - latest_thickness
    
    elevCuts <- length(elevBands) - 1
    
    colsMat <- c('#000000',brewer.pal(length(elevBands), 'Reds'))
    
    numSteps <- length(ice_thickness_all)
    featureTimeLine <- matrix(ncol=length(elevBands),nrow=length(ice_thickness_all)) * NA
    
    for(i in 1:numSteps){
      featureBrick <- DEM_glac
      sat_glac_resampled <- sat_glac_resampled_base
      png(paste(path.output,'\\',outputNames[i],".png",sep=''), width=ncol(DEM_glac), height=nrow(DEM_glac), units = "px", pointsize = 1)
      par(mar = c(0,0,0,0), xaxs = "i", yaxs = "i") #Parameters to create a borderless image
      
      feature_mask <- readOGR(dsn =feature_all[i])
      feature_mask <- spTransform(feature_mask, projec) # Transformation of shp files
      feature_mask <- crop(feature_mask,glacier_mask)
      
      area_features <- sum(area(feature_mask),na.rm=T)
      meanarea_features <- mean(area(feature_mask),na.rm=T)
      number_features <- length(feature_mask)
      
      coordfeat <- SpatialPointsDataFrame(gCentroid(feature_mask, byid=TRUE), feature_mask@data, match.ID=FALSE)
      elev_features <- extract(DEM_glac, coordfeat)
      feature_raster <- rasterize(feature_mask,DEM_glac)
      
      for(eC in 1: elevCuts){
        if(feature_Viz == 'totalarea'){
          featureTimeLine [i,eC+1] <- sum(area(feature_mask)[which(elev_features>=elevBands[eC]&elev_features<elevBands[eC + 1])],na.rm=T)}
        if(feature_Viz == 'count'){
          featureTimeLine [i,eC+1] <- length(which(elev_features>=elevBands[eC]&elev_features<elevBands[eC + 1]))}
        if(feature_Viz == 'meanarea'){
          featureTimeLine [i,eC+1] <- mean(area(feature_mask)[which(elev_features>=elevBands[eC]&elev_features<elevBands[eC + 1])],na.rm=T)}
        
        if(i>1){
          rgl::clear3d()
        }
      }
        
        sat_glac_resampled [[1]][which(feature_raster[] > 0)] <- 0
        sat_glac_resampled [[2]][which(feature_raster[] > 0)] <- 32
        sat_glac_resampled [[3]][which(feature_raster[] > 0)] <- 255

      
      plotRGB(sat_glac_resampled, stretch="lin")
      dev.off()
      
      thickAllBands <- raster(ice_thickness_all[i])
      thickAllBands <- crop(thickAllBands, glacier_mask)
      thickAllBands <- resample(thickAllBands,DEM_glac, method='bilinear')
      thickAllBands <- focal(thickAllBands, w=matrix(1,nc=5, nr=5),fun=mean)
      thickAllBands[is.na(thickAllBands)] <- 0
      
      if(feature_Viz == 'totalarea'){
        featureTimeLine [i,1] <- mean(featureTimeLine [i,2:dim(featureTimeLine)[2]],na.rm=T)}
      if(feature_Viz == 'count'){
        featureTimeLine [i,1] <- mean(featureTimeLine [i,2:dim(featureTimeLine)[2]],na.rm=T)}
      if(feature_Viz == 'meanarea'){
        featureTimeLine [i,1] <- mean(featureTimeLine [i,2:dim(featureTimeLine)[2]],na.rm=T)}    
      
      png(paste(path.output,'\\',outputNames[i],"_graph.png",sep=''), width=nrow(DEM_glac), height=ncol(DEM_glac), res = 100)
      par(mar = c(2,4,2,2), xaxs = "i", yaxs = "i") #Parameters to create a borderless image
      if(feature_Viz == 'totalarea'){
        matplot(seq(1,length(feature_all),1),featureTimeLine ,type='b',lty= 1,lwd = 2, pch=18,ylim=c(0,round(featureTimeLine[1,1]) * 3),col=colsMat,xlab = '',ylab = 'total area [m2]', xaxt="n",main=paste(glacName,outputNames[i]))
      }
      if(feature_Viz == 'count'){
        matplot(seq(1,length(feature_all),1),featureTimeLine ,type='b',lty= 1,lwd = 2, pch=18,ylim=c(0,round(featureTimeLine[1,1]) * 3),col=colsMat,xlab = '',ylab = 'count [-]', xaxt="n",main=paste(glacName,outputNames[i]))
      }
      if(feature_Viz == 'meanarea'){
        matplot(seq(1,length(feature_all),1),featureTimeLine ,type='b',lty= 1,lwd = 2, pch=18,ylim=c(0,round(featureTimeLine[1,1]) * 3),col=colsMat,xlab = '',ylab = 'mean area [m2]', xaxt="n",main=paste(glacName,outputNames[i]))
      } 
      axis(side=1, at=seq(1,length(feature_all),1), labels = outputNames)
      grid(NULL,NULL)
      if(spatialRasterViz == 1){
        legend('topright',c('mean',paste('below',elevBands[2:length(elevBands)],'m')),lty= 1,lwd = 2, pch=18,col=colsMat,bty='n')}
      dev.off()
      
      terrain_image <- matrix()
      terrain_image <- png::readPNG(paste(path.output,'\\',outputNames[i],".png",sep=''))
      terrain_image1 <- png::readPNG(paste(path.output,'\\',outputNames[i],"_graph.png",sep=''))
      
      terrain_image2 <- array(0, c(nrow(terrain_image), ncol(terrain_image), 3))
      terrain_image2[,,1] <- flipud(fliplr(apply(t(terrain_image1[,,1]), 2, rev)))
      terrain_image2[,,2] <- flipud(fliplr(apply(t(terrain_image1[,,2]), 2, rev)))
      terrain_image2[,,3] <- flipud(fliplr(apply(t(terrain_image1[,,3]), 2, rev)))
      
      terrain_imageFINAL <- array(0, c(nrow(terrain_image), ncol(terrain_image)*2, 3))
      terrain_imageFINAL[,seq(1,dim(terrain_image2)[2],1),] <- terrain_image2
      terrain_imageFINAL[,seq(dim(terrain_image2)[2]+1,dim(terrain_image2)[2]*2,1),] <- terrain_image
      terrain_imageFINAL[,seq(dim(terrain_image2)[2]-2,dim(terrain_image2)[2]+2,1),]<-0
      
      surfaceDEM <- actual_surface + thickAllBands
      finalraster<-merge(newraster,surfaceDEM)
      
      elmat = matrix(raster::extract(finalraster,raster::extent(finalraster),buffer=1000),
                     nrow = ncol(finalraster),ncol = nrow(finalraster))
      elmat %>% 
        sphere_shade(texture = "imhof1") %>%
        add_overlay(terrain_imageFINAL, alphalayer = 0.9) %>%
        add_shadow(ray_shade(elmat, anglebreaks = seq(60, 90), sunangle = 270, multicore = TRUE, lambert = FALSE, remove_edges = FALSE)) %>%
        add_shadow(ambient_shade(elmat, multicore = TRUE, remove_edges = FALSE)) %>%
        plot_3d(elmat, zscale = zscale, fov = fov,theta = theta, zoom = zoom3D, phi = 45, windowsize = c(nrow(DEM_glac) * 2,ncol(DEM_glac) * 2), solid = T)
      #render_label(elmat, x = dim(elmat)[1]/2, y = dim(elmat)[2]/2, z = cellStats(surfaceDEM,max) + 3000, zscale = 28, relativez = T, freetype = F,
      #           text = paste(glacName,outputNames[i]), textcolor = 'red', textsize = 1, linewidth = NULL)      
      
      render_snapshot(paste(path.output,'\\',outputNames[i],"_2Figs.png",sep=''))
    }
    
    pngList <- list.files(path=paste(path.output,"\\",sep=""),pattern='_2Figs.png',full.names = T)
    magick::image_write_gif(magick::image_read(pngList), 
                            path = paste(path.output,"\\",finalFileName,".gif",sep = ""), 
                            delay = 10/18)
  }
}