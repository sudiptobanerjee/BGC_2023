##The following program computes the Euclidean distance between two points point1 and point2.
 
distonplane <- function(point1,point2)
{
diff <- point1 - point2
distance <- sqrt(sum(diff^2))
distance
}


#The following program computes the distance on the surface of the earth between two points point1 and point2. Both the points are of the form (Longitude, Latitude)
  
geodetic.distance <- function(point1, point2)
{
R <- 6371
p1rad <- point1 * pi/180
p2rad <- point2 * pi/180
d <- sin(p1rad[2])*sin(p2rad[2])+cos(p1rad[2])*cos(p2rad[2])*cos(abs(p1rad[1]-p2rad[1]))	
d <- acos(d)
R*d
}

#The following function computes the distance on the surface of the earth between two points point1 and point2, by measuring the length of the chord between the two points. This approximation is better than the no.curvature.distance.stretch - both the points are of the form (Longitude, Latitude)

distance.chord <- function(point1,point2)
{
R <- 6371

p1rad <- point1 * pi/180
p2rad <- point2 * pi/180

lat <- p1rad[2]
lon <- p1rad[1]

u1 <- c(cos(lat)*cos(lon), cos(lat)*sin(lon), sin(lat))

lat <- p2rad[2]
lon <- p2rad[1]

u2 <- c(cos(lat)*cos(lon),cos(lat)*sin(lon),sin(lat))

R*sqrt(sum((u1-u2)^2))

}

#Inter-site chordal distance matrix
distance.chord.matrix <- function(long, lat) {

NSITES <- length(long)
latitude <- lat	
longitude <- long	
loc <- cbind(longitude, latitude)	

d <- matrix(nrow=NSITES, ncol=NSITES)	

for(i in 1:(NSITES-1)) {	
d[i,i] <- 0.0		
	for (j in (i+1):NSITES) {					
	d[i,j] <- distance.chord(loc[i,],loc[j,])
	d[j,i] <- d[i,j]
	}	
}
d[NSITES, NSITES] <- 0.0
d
}

#The following program computes the distance on the surface of the earth between two points point1 and point2, ignoring the curvarture, and straightening the arc of latitude and longitude. So it uses Euclidean distance to approximate the geodetic distance, making necessary scale adjustments, but often leads to gross over-estimation, as we move towards the poles. Both the points are of the form (Longitude, Latitude)

distance.stretch <- function(point1,point2)
{
R <- 6371
p1rad <- point1 * pi/180
p2rad <- point2 * pi/180
distance <- R*distonplane(p1rad, p2rad)
distance
}


#Inter-site stretched distance matrix
distance.stretch.matrix <- function(long, lat) {

NSITES <- length(long)
latitude <- lat	
longitude <- long	
loc <- cbind(longitude, latitude)	

d <- matrix(nrow=NSITES, ncol=NSITES)	

for(i in 1:(NSITES-1)) {	
d[i,i] <- 0.0		
	for (j in (i+1):NSITES) {					
	d[i,j] <- distance.stretch(loc[i,],loc[j,])
	d[j,i] <- d[i,j]
	}	
}
d[NSITES, NSITES] <- 0.0
d
}

#The following function transforms spherical coordinates to Euclidean coordinates. Input a Nx2 matrix with each row = (Longitude,Latitude) for N sites.

spherical.to.cartesian <- function(point)
{
R <- 6371

prad <- point * pi / 180

lat <- prad[,2]
lon <- prad[,1]

u <- c(R*cos(lat)*cos(lon), R*cos(lat)*sin(lon), R*sin(lat))

u 
}


#Here is my program to compute geodetic inter-site distance matrix. Input requires *vectors* in degrees. 

geodetic.distance.matrix <- function (long,lat) {		

NSITES <- length(lat)
R <- 6371

latitude <- lat	
longitude <- long	
latlong <- cbind(latitude, longitude)*pi/180	

d <- matrix(nrow=NSITES, ncol=NSITES)	

for(i in 1:(NSITES-1)) {	
d[i,i] <- 1.0		
	for (j in (i+1):NSITES) {					
	d[i,j] <- sin(latlong[i,1]) * sin(latlong[j,1]) + 			cos(latlong[i,1]) * cos(latlong[j,1])	*				cos(abs(latlong[i,2] - latlong[j,2]))
	d[j,i] <- d[i,j]
	}	
}
d[NSITES, NSITES] <- 1.0
d <- R*acos(d)
d

}

#Here is another version of the above program - with a dataframe as input, having lat and long as components.
geodetic.distance.dataframe <- function(dataframe, lat, long) {	
	NSITES <- nrow(dataframe)
		
	latitude <- dataframe$lat	
	longitude <- dataframe$long	

	latlong <- cbind(latitude, longitude)*pi/180	

	d <- matrix(nrow=NSITES, ncol=NSITES)

	for(i in 1:(NSITES-1)) {		
	d[i,i] <- 1.0
		for (j in (i+1):NSITES) {			
			d[i,j] <- sin(latlong[i,1])*sin(latlong[j,1]) + cos(latlong[i,1]) * cos(latlong[j,1])*cos(abs(latlong[i,2] - latlong[j,2]))			
			d[j,i] <- d[i,j]		
		}	
	}

d[NSITES, NSITES] <- 1.0
d <- acos(d)

for (i in 1: NSITES) { 	
	for (j in i:NSITES) { 		
		if (d[i,j] < 0.00000000001) {d[i,j] <- 0}
			d[j,i] <- d[i,j]	 	
		}
	}
d <- ifelse (d < 0.000000000001, 0.0, d)
R <- 6371
d <- R*d
d

}


#Montse Fuentes' program to compute geodetic distance: equivalent to geodetic.distance

rdistearth<-
function(loc1, loc2)
{
        if(missing(loc2))
                loc2 <- loc1
        R <- 6371
        lat <- loc1[, 2]
        lon <- loc1[, 1]
        coslat1 <- cos((lat * pi)/180)
        sinlat1 <- sin((lat * pi)/180)
        coslon1 <- cos((lon * pi)/180)
        sinlon1 <- sin((lon * pi)/180)
        lat <- loc2[, 2]
        lon <- loc2[, 1]
        coslat2 <- cos((lat * pi)/180)
        sinlat2 <- sin((lat * pi)/180)
        coslon2 <- cos((lon * pi)/180)
        sinlon2 <- sin((lon * pi)/180)
        PP1 <- cbind(coslat1 * coslon1, coslat1 * sinlon1, sinlat1)
        PP2 <- cbind(coslat2 * coslon2, coslat2 * sinlon2, sinlat2)
        pp <- (PP1 %*% t(PP2))
        R * acos(ifelse(pp > 1, 1, pp))
}

#Fuentes' style program to compute chordal distance
library(fields)
rdistchord <-
function(loc1, loc2)
{
        if(missing(loc2))
                loc2 <- loc1
        R <- 6371
        lat <- loc1[, 2]
        lon <- loc1[, 1]
        coslat1 <- cos((lat * pi)/180)
        sinlat1 <- sin((lat * pi)/180)
        coslon1 <- cos((lon * pi)/180)
        sinlon1 <- sin((lon * pi)/180)
        lat <- loc2[, 2]
        lon <- loc2[, 1]
        coslat2 <- cos((lat * pi)/180)
        sinlat2 <- sin((lat * pi)/180)
        coslon2 <- cos((lon * pi)/180)
        sinlon2 <- sin((lon * pi)/180)
        PP1 <- cbind(coslat1 * coslon1, coslat1 * sinlon1, sinlat1)
        PP2 <- cbind(coslat2 * coslon2, coslat2 * sinlon2, sinlat2)
        pp <- rdist(PP1,PP2)
        R * pp
}

#Fuentes' style program to compute euclid distance

rdisteuclid <-
function(loc1, loc2)
{
        if(missing(loc2))
                loc2 <- loc1
        R <- 6371
        lat <- loc1[, 2]
        lon <- loc1[, 1]
        x1 <- lat * (pi/180)
        y1 <- lon * (pi/180)
        lat <- loc2[, 2]
        lon <- loc2[, 1]
        x2 <- lat * (pi/180)
        y2 <- lon * (pi/180)
        PP1 <- cbind(x1, y1)
        PP2 <- cbind(x2, y2)
        pp <- rdist(PP1,PP2)
        R * pp
}

#Montse Fuentes' projection on the plane Centered around the center of gravity

lonlat.to.miles<-
function(lon.lat)
{
        x <- lon.lat[, 1]
        y <- lon.lat[, 2]
        mx <- mean(x)
        my <- mean(y)
        temp <- cbind(rep(mx, 2), range(y))
        sy <- rdistearth(temp)[2, 1]
        temp <- cbind(range(x), rep(my, 2))
        sx <- rdistearth(temp)[2, 1]
        temp <- list(x = sx/(max(x) - min(x)), y = sy/(max(y) - min(y)))
        cbind((x - mx) * temp$x, (y - my) * temp$y)
}

