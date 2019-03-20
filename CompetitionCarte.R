# ============================================================================================================================
# --------------------------------------------------------- LOAD DATA --------------------------------------------------------
# ============================================================================================================================


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


# ============================================================================================================================
# --------------------------------------------------------- BENCHMARK --------------------------------------------------------
# ============================================================================================================================

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



# ============================================================================================================================
# --------------------------------------------------- COMPUTING the twCRPS ---------------------------------------------------
# ============================================================================================================================

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




# ============================================================================================================================
# ---------------------------------------------------------- CARTES ----------------------------------------------------------
# ============================================================================================================================

# ------ Libraries 
library(ggmap)        # pour la carte
# library(viridis)    # pour les couleurs

# library(ggplot2)
library(tibble)       # comme un data.frame (pour les animations)
library(gganimate)    # animations

# Coordonnées (pour le cadre de la fenêtre)
sac_borders <- c(bottom = 11, 
                 top = 32,
                 left = 30,
                 right = 45)

# ça peut prendre du temps (baisser le zoom si jamais, mais perte de résolution)
map_background <- get_stamenmap(sac_borders, zoom = 8, maptype = "terrain-background")

# map_background <- get_stamenmap(sac_borders, zoom = 8, maptype = "watercolor")
# ggmap(map_background)      # carte seule


#_#_#_#_#_#_#_#_#_#_#_#_#_________Carte simple__________#_#_#_#_#_#_#_#_#_#_#_#_#

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
carte(map_background, "1985-01-01", val_set = TRUE, matrice_valid_set = Coord_validation_set)
carte(map_background, "1985-01-01", anom = FALSE)


# Si on veut afficher plusieurs graphiques dans une même fenêtre (utilisation de la librairie 'ggpubr') :
p1 <- carte(map_background, "2008-01-15", val_set = TRUE, matrice_valid_set = Coord_validation_set, return_plot = TRUE)
p2 <- carte(map_background, "2008-01-25", val_set = TRUE, matrice_valid_set = Coord_validation_set, return_plot = TRUE)
ggarrange(p1, p2, ncol = 2)



# ------ Remarques et observations 
# Les données manquantes sont (spatialement) les mêmes sur une durées de 1 mois, et ensuite elles sont changées
# (donc bon point, on garde une certaine continuité dans les données que nous avons, surtout par rapport à la matrice X.min, 
# qui sans cela serait remplie très majoritairement de valeurs Na...) : 
all.equal(which(is.na(anom.training[32,])), which(is.na(anom.training[33, ])))


#_#_#_#_#_#_#_#_#_#_#_#_#_________Carte animée__________#_#_#_#_#_#_#_#_#_#_#_#_#

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
    
    # labs(title = paste('Start: ', start_time, "    ","{round(round(frame_time, 2)/time_window * 100)}", "%"), x = 'longitude (°E)', y = 'latitude (°N)') +
    labs(title = paste('Start: ', start_time, "    ","{round(frame_time, 2)}"), x = 'longitude (°E)', y = 'latitude (°N)') +
    
    theme(axis.title = element_text(size = 10, face = "bold"), 
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5)) + 
    transition_time(Time) + 
    ease_aes('linear')        # pour avoir une transition progressive, et pas une animation saccadée
}

animated_map(start_time = "1985-01-01", time_window = 30, anomalie = TRUE)
animated_map(time_window = 30, anomalie = FALSE)


#_#_#_#_#_#_#_#_#_#_#_#_#_________Carte moyenne/variance__________#_#_#_#_#_#_#_#_#_#_#_#_#

# Fonction pour avoir la carte de la moyenne/variance des anomalies de températures selon une période donnée.
#   Which :    'mean' pour la moyenne, 'variance' pour la variance entre la période 'start' et 'end'
#   start :    début de la période
#   end :      fin de la période

Plot_mean_variance <- function(map = map_background, Which = "mean", start = "1985-01-01", end = "1985-01-31"){
  Start <- which(time == start)
  End <- which(time == end)
  
  if(Which == "mean") {
    matr <- apply(anom.training[Start:End,], 2, mean)
    MM <- as.data.frame(cbind(loc, matr))
    colnames(MM) <- c("longitude", "latitude", "autre")
  }else if(Which == "variance") {
    matr <- apply(anom.training[Start:End,], 2, var)
    MM <- as.data.frame(cbind(loc, matr))
    colnames(MM) <- c("longitude", "latitude", "autre")
  }
  
  map1 <- ggmap(map_background) +
    geom_point(data = MM, mapping = aes(x = longitude, y = latitude, 
                                        col = autre)) +
    scale_colour_gradient2(name = "[°C]", low = "#00AFBB", mid = "white", high = "#FC4E07",
                           space = "Lab", na.value = "grey50", guide = "colourbar",
                           aesthetics = "colour") +
    labs(x = "longitude (°E)", y = "latitude (°N)") + 
    ggtitle(paste("(", Which,")", "temperature anomalies (", start, "to", end , ")")) + 
    theme(axis.title = element_text(size = 10, face = "bold"), 
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  print(map1)
}

Plot_mean_variance(Which = "variance", start = "1998-01-01", end = "1998-01-31")


# ============================================================================================================================
# ---------------------------------------------------- Create hidden data ----------------------------------------------------
# ============================================================================================================================

# https://stackoverflow.com/questions/15387328/realistic-simulated-elevation-data-in-r-perlin-noise
perlin_noise <- function( 
  n = 5,   m = 7,    # Size of the grid for the vector field
  N = 100, M = 100   # Dimension of the image
) {
  # For each point on this n*m grid, choose a unit 1 vector
  vector_field <- apply( array( rnorm( 2 * n * m ), dim = c(2,n,m) ), 2:3,    # 3-dimensionnal array
                         function(u) u / sqrt(sum(u^2))         # we rescale values 
  )
  f <- function(x,y) {
    # Find the grid cell in which the point (x,y) is
    i <- floor(x)
    j <- floor(y)
    stopifnot( i >= 1 || j >= 1 || i < n || j < m )     # stop if any of the expressions are not ALL true
    # The 4 vectors, from the vector field, at the vertices of the square
    v1 <- vector_field[,i,j]
    v2 <- vector_field[,i+1,j]
    v3 <- vector_field[,i,j+1]
    v4 <- vector_field[,i+1,j+1]
    # Vectors from the point to the vertices
    u1 <- c(x,y) - c(i,j)
    u2 <- c(x,y) - c(i+1,j)
    u3 <- c(x,y) - c(i,j+1)
    u4 <- c(x,y) - c(i+1,j+1)
    # Scalar products
    a1 <- sum( v1 * u1 )
    a2 <- sum( v2 * u2 )
    a3 <- sum( v3 * u3 )
    a4 <- sum( v4 * u4 )
    # Weighted average of the scalar products
    s <- function(p){ 3 * p^2 - 2 * p^3 }
    p <- s( x - i )
    q <- s( y - j )
    b1 <- (1-p)*a1 + p*a2
    b2 <- (1-p)*a3 + p*a4
    (1-q) * b1 + q * b2
  }
  
  xs <- seq(from = 1, to = n, length = N+1)[-(N+1)]
  ys <- seq(from = 1, to = m, length = M+1)[-(M+1)]
  outer( xs, ys, Vectorize(f) )     # outer product of arrays xs and ys, with the fonction f (vectorized for efficiency) as operation
}


# a :      plus 'a' est proche de 0, moins les petites irrégularités comptent (donc grande zones)
# k :      nombre de couches de fréquence différentes
# quant :  plus c'est proche de 1, moins on a de zones inconnues 
Hide_data_custom <- function(a = 0.6, k = 5, n = 5, m = 7, quant = 0.5){
  
  # Fonction pour avoir les vecteurs des coordonnées longitudes et latitudes (uniques)
  get_size_map <- function(donnees, column = 1){
    result <- c()
    for (i in 2:16703) {
      if (donnees[i,column] > donnees[i-1, column]){
        result <- c(result, donnees[i-1, column])  # update long vector
      }
    }
    result <- c(result, donnees[length(donnees[,1]),column])
    return(result)
  }
  
  loc2 <- loc[order(loc[,2], decreasing = FALSE),]
  long <- get_size_map(loc)
  lat <- get_size_map(loc2, column = 2)
  
  # length(unique(long)) == length(long)      # TRUE
  # length(unique(lat)) == length(lat)        # TRUE
  
  # Fonction pour remplir la matrice moule de 1 si les coordonnées correspondantes sont contenues dans la matrice loc
  fillMoule <- function(mat) {
    for (n in 1:16703){
      x <- which(long == loc[n, 1])
      y <- which(lat  == loc[n, 2])
      mat[length(lat)-y+1, x] <- 1
    }
    return(mat)
  }
  
  # Matrice moule. Va contenir 1 si les coordonnées correspondantes sont contenues dans la matrice loc, et 0 sinon
  moule <- matrix(0, nrow = length(lat), ncol = length(long))
  moule <- fillMoule(moule)
  # image(moule)
  
  # On génère notre carte de bruit :
  perlin_raw <- perlin_noise(n, m, N = length(lat), M = length(long))
  for( i in 2:k ) { 
    perlin_raw <- perlin_raw + a^i * perlin_noise(2^i,2^i,length(lat),length(long))
  }
  
  # image(perlin_raw)
  qq <- quantile(perlin_raw, prob = quant)[[1]]
  
  perlin_raw <- apply(perlin_raw, c(1, 2), function(x){ return(x/max(perlin_raw)) } )    # rescale données de la carte bruit
  # hist(perlin_raw)          # distribution ~ normale
  
  # Si distribution selon loi normale, on sépare au niveau de la moyenne pour avoir autant de zones connues que inconnues
  perlin <- apply(perlin_raw, c(1, 2), function(x){ if (x >= qq){return(1)}else{return(0)} })
  # image(perlin)
  
  # les deux matrices contiennent que des 0 et 1, donc on va faire la multiplication élément par élément pour 'intersecter' les matrices
  result_matrix <- perlin * moule 
  # image(result_matrix)
  # image(result_matrix + moule)
  
  # Fonction pour retrouver notre matrice loc, avec cette fois-ci une troisième colonne pour savoir si la données est manquante ou non (0 ou 1)
  get_back_vector <- function(carte){
    loc_back <- matrix(nrow = length(loc[,1]), ncol = 3)   # pour longitude, latitude, et 3ème pour 1 si donnee connue, 0 sinon
    for (n in 1:length(loc[,1])){
      loc_back[n, 1:2] <- c(loc[n, 1], loc[n, 2])
      x <- which(long == loc[n, 1])
      y <- which(lat  == loc[n, 2])
      loc_back[n, 3] <- carte[length(lat)-y+1, x]
    }
    return(loc_back)
  }
  loc_back <- get_back_vector(result_matrix)
  
  return(loc_back)
}

# Fonction pour afficher la carte (avec les données qu'on a enlevées)
plot_map_missing_data <- function(donnees, colors = c("white", "black")){
  M <- as.data.frame(cbind(donnees[,]))
  M[,3] <- as.factor(M[,3])
  colnames(M) <- c("longitude", "latitude", "state")
  
  ggplot(data = M, mapping = aes(x = longitude, y = latitude), cex = 0.2) +
    geom_point(aes(colour = state), alpha = 0.7) +
    scale_colour_manual(values = colors) + 
    labs(x = "longitude (°E)", y = "latitude (°N)") + 
    ggtitle("state") + 
    theme(axis.title = element_text(size = 10, face = "bold"), 
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
}


loc_back <- Hide_data_custom(a = 1, k = 5, n = 7, m = 7, quant = 0.6)
plot_map_missing_data(loc_back)


# Démarche : On va créer des zones pour le validation set, et dans ces zones on va tirer un certain nombre de points (coordonnées),
# assez centrées (pas proches du bord), qui serviront pour tester nos modèles. Les zones créées sont assez similaires à celles qu'on
# a reçues.
# Pour un mois, les zones cachées restent constantes (et donc changent chaque début de mois), et on va prendre des points de validation
# pour les dates du 5, 15 et 25 du mois (remarque: la date du 5 est éloignée du 1 de plus de trois jours, donc centrée aussi temporellement)

jour <- which(time == "2007-01-01")       # /!\ à partir de 2007, on a 60% de données NaN (20% avant)
loc_backIntersect <- loc_back

for (i in 1:length(anom.training[jour,])){
  if (!(is.na(anom.training[jour,i])) && loc_back[i,3] == 0){       # donnée non manquante
    loc_backIntersect[i,3] <- 0
  }else if (!(is.na(anom.training[jour,i])) && loc_back[i,3] == 1){ # donnée qu'on a retiré mais dont on connait la valeur
    loc_backIntersect[i,3] <- 1
  }else if (is.na(anom.training[jour,i]) && loc_back[i,3] == 0){    # donnée qu'on nous a retiré mais qu'on n'a pas retiré
    loc_backIntersect[i,3] <- 2
  }else{                                                            # donnée qu'on nous a retiré et qu'on a retiré
    loc_backIntersect[i,3] <- 3
  }
}

plot_map_missing_data(loc_backIntersect, colors = c("white", "black", "grey", "red"))



#_#_#_#_#_#_#_#_#_#_#_#_#_________Proportion données NaN__________#_#_#_#_#_#_#_#_#_#_#_#_#

# On va calculer la proportion de données cachées (NaN) en fonction du temps :
proportion <- c()
for (i in 1:nrow(anom.training)) {
  prop <- length(which(is.na(anom.training[i,])))/16703
  proportion <- c(proportion, prop)
}

# On remarque que la proportion est de ~ 20% du début jusqu'au 2007-01-01, et monte à 60% après
# 2007-01-01 correspond à l'année à partir de laquelle les données du validation set sont prises.
plot(1:nrow(anom.training), proportion)
abline(v = which(time == "2007-01-01"), col = "red")
# print(proportion[which(time == "2007-01-01")-1])
# print(proportion[which(time == "2007-01-01")])










