#################
### LOAD DATA ###
#################

### Read data
load(file.choose()) ## Choose DATA_TRAINING.RData (about 939.6 Mb)...

### This loads several R objects:
# 'anom.training': matrix (dimension 11315x16703, size 1443.6 Mb). Training set. The (i,j)th entry contains the Red Sea surface temperature anomaly for the i-th day and j-th location in the training set. Missing data are coded as NA's.
# 'loc': matrix (dimension 16703x2, size 0.3 Mb). The i-th row contains the longitude/latitude of the i-th location.
# 'time': vector (length 11315, size 0.6 Mb) of class 'POSIXlt', containing the time of each measurement.
# 'year': vector (length 11315, size 0.1 Mb) of class 'numeric', containing the year of each measurement.
# 'month': vector (length 11315, size 0.1 Mb) of class 'numeric', containing the month of each measurement.
# 'day': vector (length 11315, size 0.1 Mb) of class 'numeric', containing the day of each measurement.
# 'index.training': vector (length 129243945, size 493.0 Mb) of class 'integer', containing the index of each space-time location in the training set, in order to retrieve the training data from the 'anom.training' matrix (of size 11315x16703).
# 'index.validation': vector (length 162000, size 0.6 Mb) of class 'integer', containing the index of each space-time location in the validation set. If 'data' is the true (unknown) data matrix of size 11315x16703 (days x locations), then 'data[index.validation]' retrieves the validation data.

#################
### BENCHMARK ###
#################

### The benchmark prediction is computed as follows:
# 1. We compute the spatio-temporal minimum X(s,t)=min A(s,t) (over some space-time neighborhood)
# 2. We pool all non-missing data values (from X(s,t)) over space and time
# 3. Assuming stationarity, the benchmark prediction at all locations (s,t) is defined as the empirical CDF computed from the pooled data (in step 2). These predictions are the same for all (s,t). 

library(fields) # needed to compute the distance in km from longitude/latitude coordinates 
radius <- 50 # neighborhood radius in kilometers
X.min <- matrix(nrow=nrow(anom.training),ncol=ncol(anom.training)) # this matrix will contain the spatio-temporal minimum X(s,t)=min A(s,t); object size: 1441.9 Mb

### Computation of X(s,t)=min A(s,t). CAREFUL! This takes some time (about 45min on a new laptop)!
for(j in 1:ncol(anom.training)){ # for each location
  print(j) 
  dist.j <- drop(rdist.earth(x1=matrix(loc[j,],nrow=1),x2=loc,miles=FALSE)) # compute the distance between the j-th location and all other locations
  loc.nei <- which(dist.j<=radius) # find its neighbors (within the specified radius of 50km)
  for(i in 1:nrow(anom.training)){ # for each day
    week.i <- i+c(-3:3) # define the week for i-th day (i.e., temporal neighborhood)
    week.i <- week.i[week.i>=1 & week.i<=nrow(anom.training)] # make sure the week contains valid indices
    X.min[i,j] <- min(anom.training[week.i,loc.nei]) # compute the minimum anomaly within the selected space-time neighborhood, and fill up the matrix (if the neighborhood has one or more NA's, the minimum will also be NA)
  }
}

str(X.min)
save(X.min, file = "xmin50.RData")
load(file.choose())

### Computation of benchmark prediction
F.benchmark <- ecdf(X.min[!is.na(X.min)]) # the benchmark is the empirical CDF of non-missing X(s,t)=min A(s,t) values.
plot(F.benchmark)

### Evaluation of benchmark prediction at each design points
xk <- -1+c(1:400)/100 # 'design points' used to evaluate the predicted distribution
F.benchmark.xk <- F.benchmark(xk) # evaluates the predicted distribution at each xk value

### Save prediction
n.validation <- length(index.validation) # number of validation points; this should be equal to 334800
n.xk <- length(xk) # number of 'design points'; this should be equal to 400
prediction <- matrix(F.benchmark.xk,nrow=n.validation,ncol=n.xk,byrow=TRUE) # Resulting prediction matrix for the Benchmark prediction (assuming stationarity over space and time)
name.of.team <- "benchmark" # For your own predictions, specify your team name
save(prediction,file=paste("/my/dropbox/folder/prediction_",name.of.team,".RData",sep="")) # Change path to save in the correct dropbox folder


############################
### COMPUTING the twCRPS ###
############################

### Function calculating the twCRPS
# INPUTS:
# prediction: matrix (dimension 162000x400, size 494.4 Mb). The i-th row contains the predicted distribution for the i-th point in the validation set, evaluated at each of the 400 design points
# true.observations: vector (length 162000, size 1.2 Mb). Observation vector (unknown to the teams), containing the true spatio-temporal minimum X(s,t)=min A(s,t) for each point in the validation set
twCRPS <- function(prediction,true.observations){
  n.validation <- length(true.observations) # number of observations in the validation set (here, equal to 162000)
  xk <- -1+c(1:400)/100 # 'design points' used to evaluate the predicted distribution
  n.xk <- length(xk) # number of design points (here, equal to 400)
  weight <- function(x){ # Define weight function for computing the threshold-weighted CRPS
    return( pnorm((x-1.5)/0.4) )
  }
  wxk <- weight(xk) # weights for each of the 'design points'
  
  true.observations.mat <- matrix(true.observations,nrow=n.validation,ncol=n.xk) # true observations put in matrix form
  xk.mat <- matrix(xk,nrow=n.validation,ncol=n.xk,byrow=TRUE) # design points put in matrix form
  wxk.mat <- matrix(wxk,nrow=n.validation,ncol=n.xk,byrow=TRUE) # weights put in matrix form
  
  twCRPS.res <- 4*mean((prediction-(true.observations.mat<=xk.mat))^2*wxk.mat)
  return(twCRPS.res)
}
true.observations <- X.min.true[index.validation] # this retrieves the 'true spatio-temporal minimum of anomalies', X(s,t), that are in the validation set (and therefore UNKNOWN to the teams!). Here, X.min.true is a matrix of size 11315x16703 (days x locations).
twCRPS.benchmark <- twCRPS(prediction=prediction,true.observations=true.observations)




############################
##### Data Exploration #####
############################

# ------ Libraries -------
library(ggmap)        # pour la carte
# library(viridis)    # pour les couleurs

# library(ggplot2)
library(tibble)       # comme un data.frame (pour les animations)
library(gganimate)    # animations

# -----------------------
# -------- Carte -------- 
# -----------------------

# Coordonnées (pour le cadre de la fenêtre)
sac_borders <- c(bottom = 11, 
                 top = 32,
                 left = 30,
                 right = 45)

# ça peut prendre du temps (baisser le zoom si jamais, mais perte de résolution)
map_background <- get_stamenmap(sac_borders, zoom = 8, maptype = "terrain-background")

# map_background <- get_stamenmap(sac_borders, zoom = 8, maptype = "watercolor")
# ggmap(map_background)      # carte seule


# Fonction pour afficher les anomalies de température sur la carte
#   map :                 carte de type 'ggmap'
#   time_ind :            jour (depuis la date 1985-01-01)
#   limit :               échelle de mesure des anomalies
#   val_set :             TRUE pour afficher les points du validation set, FALSE sinon
#   anom :                TRUE pour afficher les anomalies de temperature, FALSE pour afficher les donnees de X.min (donc min anomalies sur rayon 50 km et 7j)
#   matrice_valid_set :   Matrice avec les coordonnées (dans la matrice anom.training) des points du set de validation
#   return_plot :         TRUE pour retourner le graph (sous forme de variable), FALSE pour simplememt l'afficher 

carte <- function(map, time_ind =  "1985-01-01", limit = c(-2, 2), val_set = FALSE, anom = TRUE, matrice_valid_set = NULL, return_plot = FALSE) {
  time_index <- which(time == time_ind)    # pour avoir l'index associé à la date donnée
  
  if (anom){
    M <- as.data.frame(cbind(loc, anom.training[time_index, ]))
  } else {
    M <- as.data.frame(cbind(loc, X.min[time_index, ]))
  }
  colnames(M) <- c("longitude", "latitude", "anomalies")
  
  map <- ggmap(map) +
    geom_point(data = M, mapping = aes(x = longitude, y = latitude, 
                                       col = anomalies)) +
    # scale_color_distiller(palette = "Spectral", direction = -1, limits = c(-2, 2))    # différents choix pour changer la couleur
    # scale_color_viridis("[°C]", option = "plasma", limits = limit) + 
    scale_colour_gradient2(name = "[°C]", low = "#00AFBB", mid = "white", high = "#FC4E07",
                          space = "Lab", na.value = "grey50", guide = "colourbar",
                          aesthetics = "colour", limits = limit) +
    labs(x = "longitude (°E)", y = "latitude (°N)") + 
    ggtitle(paste("Temperature anomalies (", time[time_index], ")")) + 
    theme(axis.title = element_text(size = 10, face = "bold"), 
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  
  if (val_set) {
    # Pour afficher en noir les points du validation set
  
    # all.equal(is.na(anom.training[index.validation]), rep(TRUE, length(anom.training[index.validation])))   # ok, toutes les valeurs du validation set sont bien inconnues
    
    if (is.null(matrice_valid_set)) {       # si la matrice donnée est nulle, on la calcule
      matrice_valid_set <- get_coord_val_set()
    }

    # print(length(matrice_valid_set[which(matrice_valid_set[,1] == time_index), 2]))      # on vérifie bien qu'on a à chaque fois 500 valeurs 
    Carte <- map + geom_point(data = M[matrice_valid_set[which(matrice_valid_set[,1] == time_index), 2],], mapping = aes(x = longitude, y = latitude), col = "black")
  }else{
    Carte <- map
  }
  if (return_plot){
    return(print(Carte))
  }else{
    return(Carte)
  }
}

# Fonction pour calculer les positions des points de validation dans la matrice 'anom.training'
get_coord_val_set <- function() {
  Location <- (index.validation-1) %/% nrow(anom.training) + 1      
  Day <- (index.validation-1) %% nrow(anom.training) + 1
  
  NewA <- matrix(nrow = length(Location), ncol = 2)
  for (i in 1:length(Location)){
    NewA[i, ] <- c(Day[i],Location[i])
  }
  
  return(NewA)
}

# On calcule une première fois (pour éviter de le refaire à chaque fois qu'on affiche la carte)
Coord_validation_set <- get_coord_val_set()

# les points du validation set commencent à partir de 2007, avec 500 valeurs random 3 fois par mois (le 5, 15 et 25 du mois)
carte(map_background, "2008-01-25", val_set = TRUE, matrice_valid_set = Coord_validation_set)
carte(map_background, "1985-01-01", anom = FALSE)


# Si on veut afficher plusieurs graphiques dans une même fenêtre (utilisation de la librairie 'ggpubr') :
p1 <- carte(map_background, "2008-01-25", val_set = TRUE, matrice_valid_set = Coord_validation_set, return_plot = TRUE)
p2 <- carte(map_background, "2008-01-15", val_set = TRUE, matrice_valid_set = Coord_validation_set, return_plot = TRUE)
ggarrange(p1, p2, ncol = 2)



# ------ Remarques et observations ------ 
# Les données manquantes sont (spatialement) les mêmes sur une durées de 1 mois, et ensuite elles sont changées
# (donc bon point, on garde une certaine continuité dans les données que nous avons, surtout par rapport à la matrice X.min, 
# qui sans cela serait remplie très majoritairement de valeurs Na...) : 
all.equal(which(is.na(anom.training[32,])), which(is.na(anom.training[33, ])))



# ------------------------
# ------ Animations ------
# ------------------------

# Fonction pour avoir une carte animée (/!\ le temps de rendu peut prendre quelques minutes !)
#   start_time :    date de départ, de la forme "année-mois-jour"
#   time_window :   durée (en jours) de l'animation
#   anomalie :      TRUE pour afficher les anomalies de température, FALSE pour afficher les données de X.min (donc min anomalies sur rayon 50 km et 7j)

animated_map <- function(start_time = "1985-01-01", time_window = 30, anomalie = TRUE){
  index_time <- which(time == start_time)
    
  # Calcul de la matrice qui va contenir toutes les données (frame par frame)
  X <- c(); Y <- c(); Time <- c(); anom <- c()
  for (i in 1:time_window) {
    # Comme toutes les données doivent être contenues dans une seule et même matrice, on construit une matrice 'time_window'-fois plus grande : 
    if (anomalie) {
      X <- rbind(X, as.matrix(loc[,1], nrow = ncol(anom.training), ncol = 1))
      Y <- rbind(Y, as.matrix(loc[,2], nrow = ncol(anom.training), ncol = 1))
      Time <- rbind(Time, as.matrix(rep(index_time + i-1, ncol(anom.training))))
      anom <- rbind(anom, as.matrix(anom.training[index_time + i,], nrow = ncol(anom.training), ncol = 1))
    } else {
      X <- rbind(X, as.matrix(loc[,1], nrow = ncol(X.min), ncol = 1))
      Y <- rbind(Y, as.matrix(loc[,2], nrow = ncol(X.min), ncol = 1))
      Time <- rbind(Time, as.matrix(rep(index_time + i-1, ncol(X.min))))
      anom <- rbind(anom, as.matrix(X.min[index_time + i,], nrow = ncol(X.min), ncol = 1))
    }
  }
  
  M <- cbind(X, Y, Time, anom)
  colnames(M) <- c("X", "Y", "Time", "Anomalies")    # on renomme les colonnes
  M <- as_tibble(M)    # on convertit notre matrice en type 'tibble' pour pouvoir gérer l'animation
  
  # Affichage de l'animation (on n'affiche pas la carte (temps de rendu trop long))
  ggplot(M, aes(X, Y, color = Anomalies)) + 
    geom_point(alpha = 0.7, show.legend = TRUE) +
    # scale_color_viridis("[°C]", option = "plasma") +    # ou alors (juste pour changer les couleurs)
    # scale_color_distiller(palette = "Spectral", direction = -1) +
    scale_colour_gradient2(low = "#00AFBB", mid = "white", high = "#FC4E07",
                           space = "Lab", na.value = "grey50", guide = "colourbar",
                           aesthetics = "colour") +
    
    labs(title = paste('Start: ', start_time, "    ","{round(round(frame_time, 2)/time_window * 100)}", "%"), x = 'longitude (°E)', y = 'latitude (°N)') +
    theme(axis.title = element_text(size = 10, face = "bold"), 
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5)) + 
    transition_time(Time) + 
    ease_aes('linear')        # pour avoir une transition progressive, et pas une animation saccadée
}

animated_map(time_window = 30, anomalie = TRUE)
animated_map(time_window = 30, anomalie = FALSE)




