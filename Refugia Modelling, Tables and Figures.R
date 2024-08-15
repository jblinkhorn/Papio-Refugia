library(terra)
library(magrittr)
library(readxl)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #This assumes you have saved the code and data in a single directory

papio <- read_excel("Supplementary Table 1.xlsx") #load papio presence locations
papio2 <- sf::st_as_sf(papio, coords = c("longitude", "latitude"))
papio_split <- split(papio2, papio2$Species)
papio_split <- c(papio_split, list(papio2))

model_names <- "12345"
baboon_names <- c("P. anubis", "P. cynocephalus", "P. hamadryas", "P. kindae", "P. papio", "P. ursinus", "Papio")

####REFUGIA MODELLING####

biggerlist <- list()
for(i in 1:length(baboon_names)){
  biglist <- list(rast(paste0(baboon_names[[i]],".nc"))) #
  for(j in 1:length(biglist)){values(biglist[[j]])[values(biglist[[j]])>1]=NA}
  biggerlist[[i]] <- biglist}

names(biggerlist) <- baboon_names

big_refugia_list <- list()
big_habitable_list <- list()
big_habitable_list2 <- list()
big_habitable_size <- list()
big_wider_refugia_list <- list()
big_wider_habitable_size <- list()
big_wider_habitable_list <- list()
big_wider_habitable_list2 <- list()
big_cells_masked_all_list <- list()

for(i in 1:length(biggerlist)){
  ####STEPWISE HABITABLE AREA AND REFUGIA####
  refugia_list <- list()
  habitable_list <- list()
  habitable_list2 <- list()
  habitable_size <- list()

    tryCatch({
      prediction_patches_list <- lapply(biggerlist[[i]][[1]], function(x) patches(x, directions = 8, zeroAsNA = TRUE, allowGaps=FALSE )) #identifies patches of habitable cells
      prediction_patches_list2 <- list()
      xxx <- vect(papio_split[[i]]) #selects all presence points
      yyy<- rasterize(xxx, prediction_patches_list[[1]]) #turns them from points to raster cells
      zzz <- cells(yyy) #cell ids
      qqq <- values(prediction_patches_list[[1]]) #values as a list
      aaa <- na.omit(unique(qqq[zzz])) #selects unique 
      prediction_patches_list4 <- list()#new list
      occupied_patches2 <- sum(prediction_patches_list[[1]]==aaa) %>% patches(directions=8,zeroAsNA=T)#which cells match the presence points; to patches
      occupied_patches3 <- as.polygons(occupied_patches2)# to polygons
      qwe <- na.omit(unique(terra::extract(prediction_patches_list[[1]], occupied_patches3)))#select the predict patches overlapping inhabited patches
      prediction_patches_list4[[1]]  <- sum(prediction_patches_list[[1]]==qwe$patches)
      prediction_patches_list4[[1]][prediction_patches_list4[[1]]==0] <- NA # set the raster to na if not one of the inhabited patches so you can use it as a mask
      
      for(k in 2:131){
        ewq <- unique(mask(prediction_patches_list[[k]], prediction_patches_list4[[k-1]])) #mask this raster by previous timestep and identiy unique patches present
        prediction_patches_list4[[k]]  <- sum(prediction_patches_list[[k]]==ewq$patches)
        prediction_patches_list4[[k]][prediction_patches_list4[[k]]==0] <- NA
      }  
      refugia_list[[1]] <- sds(prediction_patches_list4) %>% app(sum, na.rm=F) #this is where all patches intersect - i.e. the refugia
      habitable_list[[1]] <- prediction_past_sds_sum2 <- sds(prediction_patches_list4) %>% app(sum, na.rm=T) #this is the extent of all stepwise contiguous areas of habitability
      habitable_list2 <- prediction_past_sds_sum2 <- sds(prediction_patches_list4)
      habitable_size[[1]] <- unlist(lapply(prediction_patches_list4,function(x) length(which(values(x==1))))) #habitable cell counts through time - this is wrong, as it this is cumulative, rather than the largest ever habitable zone in a timeslice
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

  refugia_list[sapply(refugia_list, is.null)] <- NULL
  habitable_list[sapply(habitable_list, is.null)] <- NULL
  habitable_size[sapply(habitable_size, is.null)] <- NULL
  
  ####SUMMED HABITABLE AREA AND REFUGIA####
  wider_refugia_list <- list()
  wider_habitable_size <- list()
  wider_habitable_list <- list()
  wider_habitable_list2 <- list()
  cells_masked_all_list <- list()

    xprediction_past_list2 <- biggerlist[[i]][[1]]
    fff <- sds(xprediction_past_list2)
    summed_prediction_past <- app(xprediction_past_list2, sum, na.rm=T) #sums across SpatRastDataset;
    all_inhabitable_summed_patches <- patches(summed_prediction_past, directions = 8, zeroAsNA = TRUE, allowGaps=FALSE )#identifies patches
    refugia_patches <- patches(summed_prediction_past==131, directions=8, zeroAsNA = TRUE, allowGaps=FALSE)
    xxx <- vect(papio_split[[i]])
    yyy<- rasterize(xxx, all_inhabitable_summed_patches)
    zzz <- cells(yyy) #cell ids
    qqq <- values(all_inhabitable_summed_patches)
    aaa <- na.omit(unique(qqq[zzz]))
    occupied_patches2 <- sum(all_inhabitable_summed_patches==aaa) %>% patches(directions=8,zeroAsNA=T)
    occupied_patches2_polygons <- as.polygons(occupied_patches2) #creates polygon of occupied patches
    masked_sum_patches <- mask(all_inhabitable_summed_patches, occupied_patches2_polygons)
    masked_refugia <- mask(refugia_patches, occupied_patches2_polygons)
    masked_all <- lapply(xprediction_past_list2, function(y) mask(y, occupied_patches2_polygons))
    length_masked_all <- unlist(lapply(masked_all,function(x) length(which(values(x==1)))))
    cells_masked_all <- lapply(masked_all,function(x) which(values(x==1)))
    wider_habitable_size[[1]] <- length_masked_all
    wider_refugia_list[[1]] <- masked_refugia
    wider_habitable_list[[1]] <- sds(masked_all) %>% app(sum, na.rm=T)
    wider_habitable_list2 <- sds(masked_all)
    cells_masked_all_list[[1]] <- cells_masked_all
  
  big_refugia_list[[i]] <- refugia_list
  big_habitable_list[[i]] <- habitable_list
  big_habitable_list2[[i]] <- habitable_list2
  big_habitable_size[[i]] <- habitable_size
  big_wider_refugia_list[[i]] <- wider_refugia_list
  big_wider_habitable_size[[i]] <- wider_habitable_size
  big_wider_habitable_list[[i]] <- wider_habitable_list
  big_wider_habitable_list2[[i]] <- wider_habitable_list2
  big_cells_masked_all_list[[i]] <- cells_masked_all_list
}

####Table 1 Data####
Table1 <- list()
#Step-Wise
ref_stats_list <- list()
for(k_popn in 1:7){
  tryCatch({
ref_stats <- data.frame()#create dataframe
for(i in c(1:length(big_refugia_list[[k_popn]]))){
  ref_stats[i,1] <- length(cells(big_refugia_list[[k_popn]][[i]]))*3.08025 #cell count for masked refugia #to thousand sq kms
  ref_stats[i,2] <- length(cells(big_habitable_list[[k_popn]][[i]]))*3.08025 #cell count for habitable zones#to thousand sq kms
  ref_stats[i,3] <- round(length(cells(big_refugia_list[[k_popn]][[i]]))/length(cells(big_habitable_list[[k_popn]][[i]])), 3)
  }
names(ref_stats) <- c("Refugia cells", "Habitable cells", "Proportional Refugia")
rownames(ref_stats) <- model_names
ref_stats_list[[k_popn]] <- ref_stats
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
names(ref_stats_list) <- baboon_names
ref_stats_list <- ref_stats_list[!sapply(ref_stats_list, is.null)]
Table1[[1]] <- do.call("rbind", ref_stats_list)

#Summed
for(k_popn in 1:7){
  tryCatch({
    ref_stats <- data.frame()#create dataframe
    for(i in c(1:length(big_refugia_list[[k_popn]]))){
      ref_stats[i,1] <- length(cells(big_wider_refugia_list[[k_popn]][[i]]))*3.08025 #cell count for masked refugia #to thousand sq kms
      ref_stats[i,2] <- length(cells(big_wider_habitable_list[[k_popn]][[i]]))*3.08025 #cell count for habitable zones#to thousand sq kms
      ref_stats[i,3] <- round(length(cells(big_wider_refugia_list[[k_popn]][[i]]))/length(cells(big_wider_habitable_list[[k_popn]][[i]])), 3)
    }
    names(ref_stats) <- c("Refugia cells", "Habitable cells", "Proportional Refugia")
    rownames(ref_stats) <- model_names
    ref_stats_list[[k_popn]] <- ref_stats
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
names(ref_stats_list) <- baboon_names
ref_stats_list <- ref_stats_list[!sapply(ref_stats_list, is.null)]
Table1[[2]] <- do.call("rbind", ref_stats_list)
refugia_method <- c("Step-Wise", "Summed")
for(i in 1:2){xlsx::write.xlsx(Table1[[i]], paste0("Table 2 model_metrics_ref_stats", ".xlsx"), sheetName = refugia_method[[i]], append = T)}

####Figure 2####

big_habitable_stack2 <- lapply(big_habitable_list, rast)
big_refugia_stack <- lapply(big_refugia_list, rast)
big_habitable_stack <- lapply(big_habitable_list, rast)
big_refugia_stack <- lapply(big_refugia_list, rast)
big_wider_habitable_stack <- lapply(big_wider_habitable_list, rast)
big_wider_refugia_stack <- lapply(big_wider_refugia_list, rast)
big_wider_habitable_stack <- lapply(big_wider_habitable_list, rast)
big_wider_refugia_stack <- lapply(big_wider_refugia_list, rast)

####
outline <- vect("africa_arabia.shp")
par(mfrow=c(2,2))

library(ggplot2)

par(mfrow=c(1,2))
for(i in 1:7){
plot(sum(big_habitable_stack[[i]], na.rm=TRUE)/length(big_habitable_list[[i]]), col=rev(map.pal("viridis", 131)), range= c(0, 131))
lines(outline)
points(papio_split[[i]], pch=21, col="grey", bg=alpha("white",0.1), cex=1)
lines(as.polygons(sum(big_refugia_stack[[i]], na.rm=TRUE)/(length(big_refugia_list[[i]]))), col="red", lwd=1.5)

plot(sum(big_wider_habitable_stack[[i]], na.rm=TRUE)/length(big_wider_habitable_list[[i]]), col=rev(map.pal("viridis", 131)), range= c(0, 131))
lines(outline)
points(papio_split[[i]], pch=21, col="grey", bg=alpha("white",0.1), cex=1)
lines(as.polygons(sum(big_wider_refugia_stack[[i]], na.rm=TRUE)/(length(big_wider_refugia_list[[i]]))), col="red", lwd=1.5)
}



####Figure 3####
papio_insol <- readr::read_table("papio_insol.csv", col_names=FALSE)
par(mfcol = c(7,1), mai = c(0.25, .25, 0.1, .75))
y_lims <- list()
for(i in 1:7){
  y_lims[[i]] <- c(min(unlist(as.data.frame(big_wider_habitable_size[[i]]) %>% rowSums()*3.08025/length(big_wider_habitable_size[[i]])),
                       unlist(as.data.frame(big_habitable_size[[i]]) %>% rowSums()*3.08025/length(big_habitable_size[[i]]))),
                   max(unlist(as.data.frame(big_habitable_size[[i]]) %>% rowSums()*3.08025/length(big_habitable_size[[i]])),
                       unlist(as.data.frame(big_wider_habitable_size[[i]]) %>% rowSums()*3.08025/length(big_wider_habitable_size[[i]]))
                       ))
}

time <- papio_insol$X1*-1
for(i in 1:6){
  y1 <- papio_insol$X2
  y2 <- papio_insol$X3
  y3 <- as.data.frame(big_wider_habitable_size[[i]]) %>% rowSums()*3.08025/length(big_wider_habitable_size[[i]])
  y4 <- as.data.frame(big_habitable_size[[i]]) %>% rowSums()*3.08025/length(big_habitable_size[[i]])
  y <- list(y1, y2, y3, y4)
  colors = c("#0077BB", "#33BBEE", "#CC3311", "#EE7733")
  axis_colors <- c("#0077BB", "#33BBEE", "black", "black")
  plot(time, y[[1]], yaxt = "n", xlab = "", main = "", ylab = "", type="l", lwd=2, col=colors[[1]], lty=2, xaxt = "n")
  axis(at = c(0.01, 0.02, 0.03, 0.04, 0.05), side = 4, col=colors[1])
  par(new = TRUE)
  plot(time, y[[2]], axes = FALSE, col = colors[2], xlab = "", ylab = "", type="l", lty=2, lwd=2, xaxt = "n")
  axis(at = c(-0.04, 0, 0.04), side = 4, col = colors[[2]], line=3)
  par(new = TRUE)
  plot(rev(time), y[[3]], col = colors[3], xlab = "", ylab = "", type="l", lty=1, lwd=4, ylim = y_lims[[i]], xaxt = "n")#, yaxt="n")
  axis(1, at=c(0, 20, 40, 60, 80, 100, 120), labels=c("", "", "", "", "", "", ""))
  axis(2, at=pretty(unlist(y_lims[[i]])), labels=F)
  lines(rev(time), y[[4]], col=colors[4], lwd=4)
  par(new = TRUE)
  plot(time, y[[1]], yaxt = "n", xlab = "", main = "", ylab = "", type="l", lwd=0, col=scales::alpha(colors[[1]], 0), lty=2, xaxt = "n")
  axis(at = c(0.01, 0.02, 0.03, 0.04, 0.05), side = 4, col=colors[1])
}

y[[3]] <- as.data.frame(big_wider_habitable_size[[7]]) %>% rowSums()*3.08025/length(big_wider_habitable_size[[7]])
y[[4]] <- as.data.frame(big_habitable_size[[7]]) %>% rowSums()*3.08025/length(big_habitable_size[[7]])
plot(time, y[[1]], yaxt = "n", xlab = "", main = "", ylab = "", type="l", lwd=2, col=colors[[1]], lty=2, xaxt = "n")
axis(at = c(0.01, 0.02, 0.03, 0.04, 0.05), side = 4, col=colors[1])
par(new = TRUE)
plot(time, y[[2]], axes = FALSE, col = colors[2], xlab = "", ylab = "", type="l", lty=2, lwd=2, xaxt = "n")
axis(at = c(-0.04, 0, 0.04), side = 4, col = colors[[2]], line=3)
par(new = TRUE)
plot(rev(time), y[[3]], col = colors[3], xlab = "", ylab = "", type="l", lty=1, lwd=4, ylim = y_lims[[7]])#, yaxt="n")
lines(rev(time), y[[4]], col=colors[4], lwd=4)
axis(at = pretty(unlist(y_lims)), side=2, labels=F)
par(new = TRUE)
plot(time, y[[1]], yaxt = "n", xlab = "", main = "", ylab = "", type="l", lwd=0, col=scales::alpha(colors[[1]], 0), lty=2, xaxt = "n")
axis(at = c(0.01, 0.02, 0.03, 0.04, 0.05), side = 4, col=colors[1])

####Figure 4####
over_list <- list()
under_list <-list()
df <- list()
overlap_list <- list()  
for(k_model in 1:6){
  k_model <- 1
  for(i_popn in 1:7){
    for(m_popn in 1:7){
      for(j_time in 1:131){df[[j_time]] <- length(which(big_cells_masked_all_list[[m_popn]][[k_model]][[j_time]]%in%big_cells_masked_all_list[[i_popn]][[k_model]][[j_time]]))}
      under_list[[m_popn]] <- df
    }
    over_list[[i_popn]] <- under_list}
  overlap_list[[k_model]] <- over_list}


papio_cols <- c(rgb(56, 168, 0, maxColorValue=255, alpha=255),
                rgb(255, 255, 0, maxColorValue=255, alpha=255),
                rgb(204, 204, 204, maxColorValue=255, alpha=255),
                rgb(245, 202, 122, maxColorValue=255, alpha=255),
                rgb(255, 0, 0, maxColorValue=255, alpha=255),
                rgb(168, 56, 0, maxColorValue=255, alpha=255))

a_list <- list()
for(i_popn in 1:6){
  b_list <- list()  
  for(m_popn in 1:6){
    ccc <- list()
    for(j_time in 1:131){
      ccc[[j_time]] <- overlap_list[[1]][[i_popn]][[m_popn]][[j_time]]
    }
    b_list[[m_popn]] <- ccc}  
  a_list[[i_popn]] <- b_list}

ys <- list()
ys[[1]] <-c(min(c(unlist(a_list[[1]][[2]]), unlist(a_list[[1]][[3]]), unlist(a_list[[1]][[4]]), unlist(a_list[[1]][[5]]), unlist(a_list[[1]][[6]]))),
            max(c(unlist(a_list[[1]][[2]]), unlist(a_list[[1]][[3]]), unlist(a_list[[1]][[4]]), unlist(a_list[[1]][[5]]), unlist(a_list[[1]][[6]]))))
ys[[2]] <-c(min(c(unlist(a_list[[2]][[1]]), unlist(a_list[[2]][[3]]), unlist(a_list[[2]][[4]]), unlist(a_list[[2]][[5]]), unlist(a_list[[2]][[6]]))),
            max(c(unlist(a_list[[2]][[1]]), unlist(a_list[[2]][[3]]), unlist(a_list[[2]][[4]]), unlist(a_list[[2]][[5]]), unlist(a_list[[2]][[6]]))))
ys[[3]] <-c(min(c(unlist(a_list[[3]][[2]]), unlist(a_list[[3]][[1]]), unlist(a_list[[3]][[4]]), unlist(a_list[[3]][[5]]), unlist(a_list[[3]][[6]]))),
            max(c(unlist(a_list[[3]][[2]]), unlist(a_list[[3]][[1]]), unlist(a_list[[3]][[4]]), unlist(a_list[[3]][[5]]), unlist(a_list[[3]][[6]]))))
ys[[4]] <-c(min(c(unlist(a_list[[4]][[2]]), unlist(a_list[[4]][[3]]), unlist(a_list[[4]][[1]]), unlist(a_list[[4]][[5]]), unlist(a_list[[4]][[6]]))),
            max(c(unlist(a_list[[4]][[2]]), unlist(a_list[[4]][[3]]), unlist(a_list[[4]][[1]]), unlist(a_list[[4]][[5]]), unlist(a_list[[4]][[6]]))))
ys[[5]] <-c(min(c(unlist(a_list[[5]][[2]]), unlist(a_list[[5]][[3]]), unlist(a_list[[5]][[4]]), unlist(a_list[[5]][[1]]), unlist(a_list[[5]][[6]]))),
            max(c(unlist(a_list[[5]][[2]]), unlist(a_list[[5]][[3]]), unlist(a_list[[5]][[4]]), unlist(a_list[[5]][[1]]), unlist(a_list[[5]][[6]]))))
ys[[6]] <-c(min(c(unlist(a_list[[6]][[2]]), unlist(a_list[[6]][[3]]), unlist(a_list[[6]][[4]]), unlist(a_list[[6]][[5]]), unlist(a_list[[6]][[1]]))),
            max(c(unlist(a_list[[6]][[2]]), unlist(a_list[[6]][[3]]), unlist(a_list[[6]][[4]]), unlist(a_list[[6]][[5]]), unlist(a_list[[6]][[1]]))))

time <- papio_insol$X1*-1
index <- list(c(1,2,3,4,5,6),c(2,1,3,4,5,6),c(3,1,2,4,5,6),c(4,1,2,3,5,6),c(5,1,2,3,4,6),c(6,1,2,3,4,5))
par(mfcol = c(6,1), mai = c(0.25, .25, 0.1, .75))

i <- 1
easy <- index[[i]]
y1 <- papio_insol$X2
y2 <- papio_insol$X3
y3 <- as.data.frame(a_list[[easy[1]]][[easy[2]]])
y4 <- as.data.frame(a_list[[easy[1]]][[easy[3]]])
y5 <- as.data.frame(a_list[[easy[1]]][[easy[4]]])
y6 <- as.data.frame(a_list[[easy[1]]][[easy[5]]])
y7 <- as.data.frame(a_list[[easy[1]]][[easy[6]]])

y <- list(y1, y2, y3, y4, y5, y6, y7)
colors = c("#0077BB", "#33BBEE")
axis_colors <- c("#0077BB", "#33BBEE", "black", "black")

plot(time, y[[1]], yaxt = "n", xlab = "", main = "", ylab = "", type="l", lwd=2, col=colors[[1]], lty=2, xaxt = "n")
axis(at = c(0.01, 0.02, 0.03, 0.04, 0.05), side = 4, col=colors[1])
par(new = TRUE)

plot(time, y[[2]], axes = FALSE, col = colors[2], xlab = "", ylab = "", type="l", lty=2, lwd=2, xaxt = "n")
axis(at = c(-0.04, 0, 0.04), side = 4, col = colors[[2]], line=3)

par(new = TRUE)
plot(rev(time), y[[3]], col = papio_cols[easy[2]], xlab = "", ylab = "", type="l", lty=1, lwd=4, ylim = ys[[i]], xaxt = "n")
axis(1, at=c(0, 20, 40, 60, 80, 100, 120), labels=c("", "", "", "", "", "", ""))
axis(2, at=pretty(unlist(y_lims[[i]])), labels=F)
lines(rev(time), y[[4]], col=papio_cols[easy[3]], lwd=4)
lines(rev(time), y[[5]], col=papio_cols[easy[4]], lwd=4)
lines(rev(time), y[[6]], col=papio_cols[easy[5]], lwd=4)
lines(rev(time), y[[7]], col=papio_cols[easy[6]], lwd=4)

par(new = TRUE)
plot(time, y[[1]], yaxt = "n", xlab = "", main = "", ylab = "", type="l", lwd=0, col=scales::alpha(colors[[1]], 0), lty=2, xaxt = "n")
axis(at = c(0.01, 0.02, 0.03, 0.04, 0.05), side = 4, col=colors[1])

for(i in 2:5){
  easy <- index[[i]]
  y1 <- papio_insol$X2
  y2 <- papio_insol$X3
  y3 <- as.data.frame(a_list[[easy[1]]][[easy[2]]])
  y4 <- as.data.frame(a_list[[easy[1]]][[easy[3]]])
  y5 <- as.data.frame(a_list[[easy[1]]][[easy[4]]])
  y6 <- as.data.frame(a_list[[easy[1]]][[easy[5]]])
  y7 <- as.data.frame(a_list[[easy[1]]][[easy[6]]])
  
  y <- list(y1, y2, y3, y4, y5, y6, y7)
  colors = c("#0077BB", "#33BBEE")
  axis_colors <- c("#0077BB", "#33BBEE", "black", "black")
  
  plot(time, y[[1]], yaxt = "n", xlab = "", main = "", ylab = "", type="l", lwd=2, col=colors[[1]], lty=2, xaxt = "n")
  axis(at = c(0.01, 0.02, 0.03, 0.04, 0.05), side = 4, col=colors[1])
  par(new = TRUE)
  
  plot(time, y[[2]], axes = FALSE, col = colors[2], xlab = "", ylab = "", type="l", lty=2, lwd=2, xaxt = "n")
  axis(at = c(-0.04, 0, 0.04), side = 4, col = colors[[2]], line=3)
  
  par(new = TRUE)
  plot(rev(time), y[[3]], col = papio_cols[easy[2]], xlab = "", ylab = "", type="l", lty=5, lwd=4, ylim = ys[[i]], xaxt = "n")
  axis(1, at=c(0, 20, 40, 60, 80, 100, 120), labels=c("", "", "", "", "", "", ""))
  axis(2, at=pretty(unlist(y_lims[[i]])), labels=F)
  lines(rev(time), y[[4]], col=papio_cols[easy[3]], lwd=4)
  lines(rev(time), y[[5]], col=papio_cols[easy[4]], lwd=4)
  lines(rev(time), y[[6]], col=papio_cols[easy[5]], lwd=4)
  lines(rev(time), y[[7]], col=papio_cols[easy[6]], lwd=4)
  
  par(new = TRUE)
  plot(time, y[[1]], yaxt = "n", xlab = "", main = "", ylab = "", type="l", lwd=0, col=scales::alpha(colors[[1]], 0), lty=2, xaxt = "n")
  axis(at = c(0.01, 0.02, 0.03, 0.04, 0.05), side = 4, col=colors[1])
}

i <- 6
easy <- index[[i]]
y1 <- papio_insol$X2
y2 <- papio_insol$X3
y3 <- as.data.frame(a_list[[easy[1]]][[easy[2]]])
y4 <- as.data.frame(a_list[[easy[1]]][[easy[3]]])
y5 <- as.data.frame(a_list[[easy[1]]][[easy[4]]])
y6 <- as.data.frame(a_list[[easy[1]]][[easy[5]]])
y7 <- as.data.frame(a_list[[easy[1]]][[easy[6]]])

y <- list(y1, y2, y3, y4, y5, y6, y7)
colors = c("#0077BB", "#33BBEE")
axis_colors <- c("#0077BB", "#33BBEE", "black", "black")

plot(time, y[[1]], yaxt = "n", xlab = "", main = "", ylab = "", type="l", lwd=2, col=colors[[1]], lty=2, xaxt = "n")
axis(at = c(0.01, 0.02, 0.03, 0.04, 0.05), side = 4, col=colors[1])
par(new = TRUE)

plot(time, y[[2]], axes = FALSE, col = colors[2], xlab = "", ylab = "", type="l", lty=2, lwd=2, xaxt = "n")
axis(at = c(-0.04, 0, 0.04), side = 4, col = colors[[2]], line=3)

par(new = TRUE)
plot(rev(time), y[[3]], col = papio_cols[easy[2]], xlab = "", ylab = "", type="l", lty=5, lwd=4, ylim = ys[[i]])
axis(2, at=pretty(unlist(y_lims[[i]])), labels=F)
lines(rev(time), y[[4]], col=papio_cols[easy[3]], lwd=4)
lines(rev(time), y[[5]], col=papio_cols[easy[4]], lwd=4)
lines(rev(time), y[[6]], col=papio_cols[easy[5]], lwd=4)
lines(rev(time), y[[7]], col=papio_cols[easy[6]], lwd=4)

par(new = TRUE)
plot(time, y[[1]], yaxt = "n", xlab = "", main = "", ylab = "", type="l", lwd=0, col=scales::alpha(colors[[1]], 0), lty=2, xaxt = "n")
axis(at = c(0.01, 0.02, 0.03, 0.04, 0.05), side = 4, col=colors[1])

####Supplementary Information GIFs####


require("RColorBrewer")
require("rasterVis")
require("gridExtra")
require("scales")
require("latticeExtra")
require("animation")

outline <- vect("africa_arabia.shp")

f <- seq(0, 130000, 1000)
k <- 1
saveGIF({
  for(i in 1:131){
plot(biggerlist[[k]][[1]][[i]], col="grey", legend=F, main = expression(italic("P. anubis")), xlab=paste0(f[[i]], " years ago"))
plot(big_wider_habitable_list2[[k]][[i]], col="#CC3311", add=T, legend=F)
plot(big_habitable_list2[[k]][[i]], col="#EE7733", add=T, legend=F)
lines(outline)
lines(as.polygons(big_wider_refugia_list[[k]][[1]]), col="#33BBEE", lwd=2)
lines(as.polygons(big_refugia_list[[k]][[1]]), col="#0077BB", lwd=2)
legend("bottomleft",
       c("Predicted habitable range","Summed Range","Step-Wise Range"), fill=c("grey", "#CC3311", "#EE7733"), horiz=F, cex=0.8, inset=c(0, 0.032))
legend("bottomright",
       c("Summed Refugia","Step-Wise Refugia"), col=c("#33BBEE", "#0077BB"), lty=1, horiz=F, cex=0.8, inset=c(0, 0.032))
  }
}, interval=0.4, movie.name=paste0(baboon_names[[k]], ".gif"))

k <- 2
saveGIF({
  for(i in 1:131){
    plot(biggerlist[[k]][[1]][[i]], col="grey", legend=F, main = expression(italic("P. cynocephalus")), xlab=paste0(f[[i]], " years ago"))
    plot(big_wider_habitable_list2[[k]][[i]], col="#CC3311", add=T, legend=F)
    plot(big_habitable_list2[[k]][[i]], col="#EE7733", add=T, legend=F)
    lines(outline)
    lines(as.polygons(big_wider_refugia_list[[k]][[1]]), col="#33BBEE", lwd=2)
    lines(as.polygons(big_refugia_list[[k]][[1]]), col="#0077BB", lwd=2)
    legend("bottomleft",
           c("Predicted habitable range","Summed Range","Step-Wise Range"), fill=c("grey", "#CC3311", "#EE7733"), horiz=F, cex=0.8, inset=c(0, 0.032))
    legend("bottomright",
           c("Summed Refugia","Step-Wise Refugia"), col=c("#33BBEE", "#0077BB"), lty=1, horiz=F, cex=0.8, inset=c(0, 0.032))
  }
}, interval=0.4, movie.name=paste0(baboon_names[[k]], ".gif"))

k <- 3
saveGIF({
  for(i in 1:131){
    plot(biggerlist[[k]][[1]][[i]], col="grey", legend=F, main = expression(italic("P. hamadryas")), xlab=paste0(f[[i]], " years ago"))
    plot(big_wider_habitable_list2[[k]][[i]], col="#CC3311", add=T, legend=F)
    plot(big_habitable_list2[[k]][[i]], col="#EE7733", add=T, legend=F)
    lines(outline)
    lines(as.polygons(big_wider_refugia_list[[k]][[1]]), col="#33BBEE", lwd=2)
    lines(as.polygons(big_refugia_list[[k]][[1]]), col="#0077BB", lwd=2)
    legend("bottomleft",
           c("Predicted habitable range","Summed Range","Step-Wise Range"), fill=c("grey", "#CC3311", "#EE7733"), horiz=F, cex=0.8, inset=c(0, 0.032))
    legend("bottomright",
           c("Summed Refugia","Step-Wise Refugia"), col=c("#33BBEE", "#0077BB"), lty=1, horiz=F, cex=0.8, inset=c(0, 0.032))
  }
}, interval=0.4, movie.name=paste0(baboon_names[[k]], ".gif"))

k <- 4
saveGIF({
  for(i in 1:131){
    plot(biggerlist[[k]][[1]][[i]], col="grey", legend=F, main = expression(italic("P. kindae")), xlab=paste0(f[[i]], " years ago"))
    plot(big_wider_habitable_list2[[k]][[i]], col="#CC3311", add=T, legend=F)
    plot(big_habitable_list2[[k]][[i]], col="#EE7733", add=T, legend=F)
    lines(outline)
    lines(as.polygons(big_wider_refugia_list[[k]][[1]]), col="#33BBEE", lwd=2)
    lines(as.polygons(big_refugia_list[[k]][[1]]), col="#0077BB", lwd=2)
    legend("bottomleft",
           c("Predicted habitable range","Summed Range","Step-Wise Range"), fill=c("grey", "#CC3311", "#EE7733"), horiz=F, cex=0.8, inset=c(0, 0.032))
    legend("bottomright",
           c("Summed Refugia","Step-Wise Refugia"), col=c("#33BBEE", "#0077BB"), lty=1, horiz=F, cex=0.8, inset=c(0, 0.032))
  }
}, interval=0.4, movie.name=paste0(baboon_names[[k]], ".gif"))

k <- 5
saveGIF({
  for(i in 1:131){
    plot(biggerlist[[k]][[1]][[i]], col="grey", legend=F, main = expression(italic("P. papio")), xlab=paste0(f[[i]], " years ago"))
    plot(big_wider_habitable_list2[[k]][[i]], col="#CC3311", add=T, legend=F)
    plot(big_habitable_list2[[k]][[i]], col="#EE7733", add=T, legend=F)
    lines(outline)
    lines(as.polygons(big_wider_refugia_list[[k]][[1]]), col="#33BBEE", lwd=2)
    lines(as.polygons(big_refugia_list[[k]][[1]]), col="#0077BB", lwd=2)
    legend("bottomleft",
           c("Predicted habitable range","Summed Range","Step-Wise Range"), fill=c("grey", "#CC3311", "#EE7733"), horiz=F, cex=0.8, inset=c(0, 0.032))
    legend("bottomright",
           c("Summed Refugia","Step-Wise Refugia"), col=c("#33BBEE", "#0077BB"), lty=1, horiz=F, cex=0.8, inset=c(0, 0.032))
  }
}, interval=0.4, movie.name=paste0(baboon_names[[k]], ".gif"))

k <- 6
saveGIF({
  for(i in 1:131){
    plot(biggerlist[[k]][[1]][[i]], col="grey", legend=F, main = expression(italic("P. ursinus")), xlab=paste0(f[[i]], " years ago"))
    plot(big_wider_habitable_list2[[k]][[i]], col="#CC3311", add=T, legend=F)
    plot(big_habitable_list2[[k]][[i]], col="#EE7733", add=T, legend=F)
    lines(outline)
    lines(as.polygons(big_wider_refugia_list[[k]][[1]]), col="#33BBEE", lwd=2)
    lines(as.polygons(big_refugia_list[[k]][[1]]), col="#0077BB", lwd=2)
    legend("bottomleft",
           c("Predicted habitable range","Summed Range","Step-Wise Range"), fill=c("grey", "#CC3311", "#EE7733"), horiz=F, cex=0.8, inset=c(0, 0.032))
    legend("bottomright",
           c("Summed Refugia","Step-Wise Refugia"), col=c("#33BBEE", "#0077BB"), lty=1, horiz=F, cex=0.8, inset=c(0, 0.032))
  }
}, interval=0.4, movie.name=paste0(baboon_names[[k]], ".gif"))

k <- 7
saveGIF({
  for(i in 1:131){
    plot(biggerlist[[k]][[1]][[i]], col="grey", legend=F, main = expression("All"~italic("Papio")), xlab=paste0(f[[i]], " years ago"))
    plot(big_wider_habitable_list2[[k]][[i]], col="#CC3311", add=T, legend=F)
    plot(big_habitable_list2[[k]][[i]], col="#EE7733", add=T, legend=F)
    lines(outline)
    lines(as.polygons(big_wider_refugia_list[[k]][[1]]), col="#33BBEE", lwd=2)
    lines(as.polygons(big_refugia_list[[k]][[1]]), col="#0077BB", lwd=2)
    legend("bottomleft",
           c("Predicted habitable range","Summed Range","Step-Wise Range"), fill=c("grey", "#CC3311", "#EE7733"), horiz=F, cex=0.8, inset=c(0, 0.032))
    legend("bottomright",
           c("Summed Refugia","Step-Wise Refugia"), col=c("#33BBEE", "#0077BB"), lty=1, horiz=F, cex=0.8, inset=c(0, 0.032))
  }
}, interval=0.4, movie.name=paste0(baboon_names[[k]], ".gif"))