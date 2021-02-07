library(rgdal)
library(wordcloud)
library(reshape)
library(maptools)
library(classInt)
library(FactoMineR)
library(FNN)

########
# READ #
########

f <- readOGR("C:/Users/mathieu/Documents/data/ADMIN-EXPRESS-COG_2-0__SHP__FRA_L93_2019-09-24/ADMIN-EXPRESS-COG_2-0__SHP__FRA_2019-09-24/ADMIN-EXPRESS-COG/1_DONNEES_LIVRAISON_2019-09-24/ADE-COG_2-0_SHP_LAMB93_FR/COMMUNE.shp", "COMMUNE")
deps <- readOGR("C:/Users/mathieu/Documents/data/ADMIN-EXPRESS-COG_2-0__SHP__FRA_L93_2019-09-24/ADMIN-EXPRESS-COG_2-0__SHP__FRA_2019-09-24/ADMIN-EXPRESS-COG/1_DONNEES_LIVRAISON_2019-09-24/ADE-COG_2-0_SHP_LAMB93_FR/DEPARTEMENT.shp", "DEPARTEMENT")


####################################
# DETECTION DE LA DERNIERE SYLLABE #
####################################

vyl <- "[aeiouyîïéè]"
csn <- "[^aeiouyîïéè-]"
 
regles <- paste(paste("(", vyl, "+", csn, "*", "$", ")", sep=""),
                paste("(", vyl, "+", csn, "{1,2}", "(e)[s]*$", ")", sep=""),
                sep="|"
               )

txt <- tolower(f$NOM_COMM)
pos <- regexpr(regles, txt)
sylls <- sapply(1:length(txt), function(i) substr(txt[i], pos[i], nchar(txt[i])))


####################################
# COMPTAGE DES SYLLABES PAR REGION #
####################################

df <- data.frame(region=f$NOM_REGION, lastsyllabus=sylls, value=1)
r <- cast(df, region~lastsyllabus, sum)
rownames(r) <- r$region; r$region <- NULL


#######
# AFC #
#######

afc <- CA(r, ncp=2)

coordsReg <- afc$row$coord[, c(1,2)]
coordsSyl <- afc$col$coord[, c(1,2)]

n <- 50
nn <- get.knnx(coordsSyl, coordsReg, k=n)


#################
# COUCHE REGION #
#################

reg <- unionSpatialPolygons(deps, deps$NOM_REGION)
reg <- spChFIDs(reg, names(reg))
reg <- SpatialPolygonsDataFrame(reg, data=data.frame(nom_reg=names(reg), row.names=names(reg)))


################################
# COUCHE POSITION DES LIBELLES #
################################

out <- vector(mode="list", length=nrow(r))
for (i in seq(along=out)) {
  sylls <- names(r)[nn$nn.index[i, ]]
  weights <- nn$nn.dist[i, ]
  df <- data.frame(sylls=sylls, weights=weights, region=reg$nom_reg[i])
  pts <- spsample(reg[i, ], n=n, type="nonaligned")
  if (length(pts) < n) { # parfois, le compte n'est pas bon, donc on crée des points supplémentaires
    pts <-spRbind(spsample(reg[i, ], n=(n-length(pts)), type="random"), pts)
  }
  out[[i]] <- SpatialPointsDataFrame(pts[1:n, ], data=df)
}

labelspt <- do.call("rbind", out)


#######################
# TAILLE DES LIBELLES #
#######################

nCuts <- 6
cexs <- seq(1.2, 1.5, length.out=nCuts+1)

szl <- lapply(out, function(x) {
  if (length(unique(x$weights))==1) {
    return(rep(1.3, n))
  }
  else {
    ints <- classIntervals(x$weights, nCuts, style="jenks")
    return(cexs[findInterval(x$weights, ints$brk)])
  }
})

szs <- unlist(szl)


#########################
# COULEUR DES LIBELLES #
#########################

cols <- rainbow(nrow(reg))[as.numeric(as.factor(labelspt$region))]


#########
# CARTE #
#########

## WORDCLOUD
png(file="IMG/noms_communes.png", width=2000, height=2000)
coords <- coordinates(labelspt)
plot(coords[,1], coords[,2], type="n", axes=FALSE, xlab=NA, ylab=NA)
nc <- wordlayout(coords[,1], coords[,2], labelspt$sylls, cex=szs)
text(nc[,1]+0.5*nc[,3], nc[,2]+0.5*nc[,4], labelspt$sylls, cex=szs, col=cols)
plot(as(reg, "SpatialLines"), add=T, col=rgb(.1,.1,.1,.1), lty=2)
dev.off()

## MAPTOOLS:POINTLABELS
png(file="IMG/maptools_noms_communes.png", width=2000, height=2000)
plot(coords[, 1], coords[, 2], type="n", axes=FALSE, xlab=NA, ylab=NA)
pointLabel(coords[, 1], coords[, 2], labelspt$sylls, cex=cexs, col=cols, offset=2)
dev.off()