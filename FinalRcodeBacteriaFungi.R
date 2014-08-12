library(vegan)

fungi.species <- read.table("/home/mariya/documents/ITS/ITS.species.txt", header=T,sep="\t")
bacteria.genus <- read.table("~/documents/16S/bacteria.genus.counts", header=T, sep="\t")
row.names(fungi.species) <- fungi.species$taxon
row.names(bacteria.genus) <- bacteria.genus$Taxon_Name
species_comm <- fungi.species[,3:ncol(fungi.species)]
bacteria_comm <- bacteria.genus[,4:ncol(bacteria.genus)]
fungi.t <- t(species_comm)
bacteria.t <- t(bacteria_comm)
nrow(bacteria.t)
#384
fungi.tyson <- fungi.t[grep("...z", rownames(fungi.t)),]
fungi.tyson <- fungi.tyson[,apply(fungi.tyson,2,sum) > 0]
nrow(fungi.tyson)
#384
fungi.tyson <- fungi.tyson[apply(fungi.tyson,1,sum) >= 100,]
nrow(fungi.tyson)
#323
#384-323=61
bacteria.tyson <- bacteria.t[apply(bacteria.t,1,sum) >= 100,]
nrow(bacteria.tyson)
#314
#384-314=70
fungi.pa.tyson <- fungi.tyson
fungi.pa.tyson[fungi.pa.tyson > 0] <- 1
bacteria.pa.tyson <- bacteria.tyson
bacteria.pa.tyson[bacteria.pa.tyson > 0] <- 1
fungi.stand <- decostand(fungi.tyson, method="hellinger")
fungi.pa.stand <- decostand(fungi.pa.tyson, method="hellinger")
bacteria.stand <- decostand(bacteria.tyson, method="hellinger")
bacteria.pa.stand <- decostand(bacteria.pa.tyson, method="hellinger")

meta <- read.table("~/documents/ITS/ScriptJunk/meta_full.txt", header=T, row.names=1)
meta.Wd.Per <- meta[,c(10,12,21,22)]
meta.conduit.D <- meta[,c(12,13)]
meta.Wd.Per <- meta.Wd.Per[rowSums(is.na(meta.Wd.Per)) <= 0,]
meta.conduit.D <- meta.conduit.D[rowSums(is.na(meta.conduit.D)) <= 0,]
tapply(meta.conduit.D$Conduit.D.um., meta.conduit.D$Family, mean)
#Annonaceae    Cannabaceae Caprifoliaceae      Cornaceae   Cupressaceae      Ebenaceae       Fabaceae 
#19.00526       30.34328       17.03445       27.45968        8.28000       25.75913       33.26274 
#Fagaceae   Juglandaceae       Oleaceae       Pinaceae    Platanaceae       Rosaceae    Sapindaceae 
#37.38283       31.15201       18.75476       11.24000       26.77467       19.16279       17.96995 
#Ulmaceae       Vitaceae 
#20.05187       19.07490
tapply(meta.Wd.Per$Wd.Per_N, meta.Wd.Per$Family, mean)
#Annonaceae    Cannabaceae Caprifoliaceae      Cornaceae   Cupressaceae      Ebenaceae       Fabaceae 
#0.3197750      0.3981667             NA      0.1228750      0.1324250      0.2740667      0.2332750 
#Fagaceae   Juglandaceae       Oleaceae       Pinaceae    Platanaceae       Rosaceae    Sapindaceae 
#0.1916400      0.1351000      0.1442000      0.1053833      0.1841667      0.1220833      0.1651167 
#Ulmaceae       Vitaceae 
#0.1619000      0.2193000
tapply(meta.Wd.Per$Wd.Per_C, meta.Wd.Per$Family, mean)
#Annonaceae    Cannabaceae Caprifoliaceae      Cornaceae   Cupressaceae      Ebenaceae       Fabaceae 
#49.09500       47.81667             NA       48.18750       51.49250       48.10000       48.13250 
#Fagaceae   Juglandaceae       Oleaceae       Pinaceae    Platanaceae       Rosaceae    Sapindaceae 
#48.35000       48.05667       47.75000       51.83000       48.79333       48.58667       48.50125 
#Ulmaceae       Vitaceae 
#48.74333       48.29333
meta.Wd.Per.clade <- meta[,c(10,21,22)]
meta.Wd.Per.clade <- meta.Wd.Per.clade[rowSums(is.na(meta.Wd.Per.clade)) <= 0,]
tapply(meta.Wd.Per.clade$Wd.Per_C, meta.Wd.Per.clade$Clade, mean)
#Angiosperms  Gymnosperm 
#48.33299    51.60500
tapply(meta.Wd.Per.clade$Wd.Per_N, meta.Wd.Per.clade$Clade, mean)
#Angiosperms  Gymnosperm 
#0.2087454   0.1234111
meta.rep <- meta
index <- is.na(meta.rep$Conduit.D.um.) & meta.rep$Family == "Pinaceae"
meta.rep$Conduit.D.um.[index] <- 11.24000
index <- is.na(meta.rep$Wd.Per_N) & meta.rep$Family == "Pinaceae"
meta.rep$Wd.Per_N[index] <- 0.1053833 
index <- is.na(meta.rep$Wd.Per_N) & meta.rep$Family == "Juglandaceae"
meta.rep$Wd.Per_N[index] <- 0.1351000
index <- is.na(meta.rep$Wd.Per_N) & meta.rep$Family == "Fagaceae"
meta.rep$Wd.Per_N[index] <- 0.1916400
index <- is.na(meta.rep$Wd.Per_C) & meta.rep$Family == "Pinaceae"
meta.rep$Wd.Per_C[index] <- 51.83000
index <- is.na(meta.rep$Wd.Per_C) & meta.rep$Family == "Juglandaceae"
meta.rep$Wd.Per_C[index] <- 48.05667
index <- is.na(meta.rep$Wd.Per_C) & meta.rep$Family == "Fagaceae"
meta.rep$Wd.Per_C[index] <- 48.35000
index <- is.na(meta.rep$Wd.Per_N) & meta.rep$Family == "Caprifoliaceae"
meta.rep$Wd.Per_N[index] <- 0.2087454
index <- is.na(meta.rep$Wd.Per_C) & meta.rep$Family == "Caprifoliaceae"
meta.rep$Wd.Per_C[index] <- 48.33299

meta.fungi <- meta.rep[row.names(fungi.tyson),]
meta.bacteria <- meta.rep[row.names(bacteria.tyson),]

###***************************####
o.fungi <- subset(fungi.tyson,rownames(fungi.tyson) %in% rownames(meta.fungi[meta.fungi$HarvestYear == 1,]))
t.fungi <- subset(fungi.tyson,rownames(fungi.tyson) %in% rownames(meta.fungi[meta.fungi$HarvestYear == 3,]))
o.meta.fungi <- meta.fungi[row.names(o.fungi),]
t.meta.fungi <- meta.fungi[row.names(t.fungi),]
o.fungi.stand <- decostand(o.fungi, method="hellinger")
o.fungi.dist <- vegdist(o.fungi.stand, method="euclidean")
o.fungi.mds <- metaMDS(o.fungi.stand, dist = "euclidean")
t.fungi.stand <- decostand(t.fungi, method="hellinger")
t.fungi.dist <- vegdist(t.fungi.stand, method="euclidean")
t.fungi.mds <- metaMDS(t.fungi.stand, dist = "euclidean")
#t.fungi.adonis <- adonis(t.fungi.dist ~ Conduit.D.um. * Conduit.L.m. + avg.top.bottom.density + WCProtPerDM * WLigPerDM * WPPerDM * Cellulose.cell * Hemicellulose.hemic + PlotLocation * Pos * factor(PotID), t.meta.fungi) 
#o.fungi.adonis <- adonis(o.fungi.dist ~ Conduit.D.um. * Conduit.L.m. + avg.top.bottom.density + WCProtPerDM * WLigPerDM * WPPerDM * Cellulose.cell * Hemicellulose.hemic + PlotLocation * Pos * factor(PotID), o.meta.fungi) 
t.fungi.ef <- envfit(t.fungi.mds ~ Conduit.D.um. * Conduit.L.m. + avg.top.bottom.density + WCProtPerDM * WLigPerDM * WPPerDM * Cellulose.cell * Hemicellulose.hemic + PlotLocation * Pos * factor(PotID), t.meta.fungi)   
o.fungi.ef <- envfit(o.fungi.mds ~ Conduit.D.um. * Conduit.L.m. + avg.top.bottom.density + WCProtPerDM * WLigPerDM * WPPerDM * Cellulose.cell * Hemicellulose.hemic + PlotLocation * Pos * factor(PotID), o.meta.fungi)   
ordiplot(t.fungi.mds, type="p", display="sites")
plot(t.fungi.ef)
ordiplot(o.fungi.mds, type="p", display="sites")
plot(o.fungi.ef)

o.pa.fungi <- subset(fungi.pa.tyson,rownames(fungi.pa.tyson) %in% rownames(meta.fungi[meta.fungi$HarvestYear == 1,]))
t.pa.fungi <- subset(fungi.pa.tyson,rownames(fungi.pa.tyson) %in% rownames(meta.fungi[meta.fungi$HarvestYear == 3,]))
o.pa.meta.fungi <- meta.fungi[row.names(o.pa.fungi),]
t.pa.meta.fungi <- meta.fungi[row.names(t.pa.fungi),]
o.pa.fungi.stand <- decostand(o.pa.fungi, method="hellinger")
o.pa.fungi.dist <- vegdist(o.pa.fungi.stand, method="euclidean")
o.pa.fungi.mds <- metaMDS(o.pa.fungi.stand, dist = "euclidean")
t.pa.fungi.stand <- decostand(t.pa.fungi, method="hellinger")
t.pa.fungi.dist <- vegdist(t.pa.fungi.stand, method="euclidean")
t.pa.fungi.mds <- metaMDS(t.pa.fungi.stand, dist = "euclidean")
t.pa.fungi.ef <- envfit(t.pa.fungi.mds ~ Conduit.D.um. * Conduit.L.m. + avg.top.bottom.density + WCProtPerDM * WLigPerDM * WPPerDM * Cellulose.cell * Hemicellulose.hemic + PlotLocation * Pos * factor(PotID), t.meta.fungi)   
o.pa.fungi.ef <- envfit(o.pa.fungi.mds ~ Conduit.D.um. * Conduit.L.m. + avg.top.bottom.density + WCProtPerDM * WLigPerDM * WPPerDM * Cellulose.cell * Hemicellulose.hemic + PlotLocation * Pos * factor(PotID), o.meta.fungi)   
#t.pa.fungi.adonis <- adonis(t.pa.fungi.dist ~ Conduit.D.um. * Conduit.L.m. + avg.top.bottom.density + WCProtPerDM * WLigPerDM * WPPerDM * Cellulose.cell * Hemicellulose.hemic + PlotLocation * Pos * factor(PotID), t.meta.fungi) 
#o.pa.fungi.adonis <- adonis(o.pa.fungi.dist ~ Conduit.D.um. * Conduit.L.m. + avg.top.bottom.density + WCProtPerDM * WLigPerDM * WPPerDM * Cellulose.cell * Hemicellulose.hemic + PlotLocation * Pos * factor(PotID), o.meta.fungi) 
ordiplot(t.pa.fungi.mds, type="p", display="sites")
plot(t.pa.fungi.ef)
ordiplot(o.pa.fungi.mds, type="p", display="sites")
plot(o.pa.fungi.ef)

o.bacteria <- subset(bacteria.tyson,rownames(bacteria.tyson) %in% rownames(meta.bacteria[meta.bacteria$HarvestYear == 1,]))
t.bacteria <- subset(bacteria.tyson,rownames(bacteria.tyson) %in% rownames(meta.bacteria[meta.bacteria$HarvestYear == 3,]))
o.meta.bacteria <- meta.bacteria[row.names(o.bacteria),]
t.meta.bacteria <- meta.bacteria[row.names(t.bacteria),]
o.bacteria.stand <- decostand(o.bacteria, method="hellinger")
o.bacteria.dist <- vegdist(o.bacteria.stand, method = "euclidean")
o.bacteria.mds <- metaMDS(o.bacteria.stand, dist = "euclidean")
t.bacteria.stand <- decostand(t.bacteria, method="hellinger")
t.bacteria.dist <- vegdist(t.bacteria.stand, method = "euclidean")
t.bacteria.mds <- metaMDS(t.bacteria.stand, dist = "euclidean")
#t.bacteria.adonis <- adonis(t.bacteria.dist ~ Conduit.D.um. * Conduit.L.m. + avg.top.bottom.density + WCProtPerDM * WLigPerDM * WPPerDM * Cellulose.cell * Hemicellulose.hemic + PlotLocation * Pos * factor(PotID), t.meta.bacteria) 
#o.bacteria.adonis <- adonis(o.bacteria.dist ~ Conduit.D.um. * Conduit.L.m. + avg.top.bottom.density + WCProtPerDM * WLigPerDM * WPPerDM * Cellulose.cell * Hemicellulose.hemic + PlotLocation * Pos * factor(PotID), o.meta.bacteria) 
t.ef.bacteria <- envfit(t.bacteria.mds ~ Conduit.D.um. * Conduit.L.m. + avg.top.bottom.density + WCProtPerDM * WLigPerDM * WPPerDM * Cellulose.cell * Hemicellulose.hemic + PlotLocation * Pos * factor(PotID), t.meta.bacteria)   
o.ef.bacteria <- envfit(o.bacteria.mds ~ Conduit.D.um. * Conduit.L.m. + avg.top.bottom.density + WCProtPerDM * WLigPerDM * WPPerDM * Cellulose.cell * Hemicellulose.hemic + PlotLocation * Pos * factor(PotID), o.meta.bacteria)   

o.pa.bacteria <- subset(bacteria.pa.tyson,rownames(bacteria.pa.tyson) %in% rownames(meta.bacteria[meta.bacteria$HarvestYear == 1,]))
t.pa.bacteria <- subset(bacteria.pa.tyson,rownames(bacteria.pa.tyson) %in% rownames(meta.bacteria[meta.bacteria$HarvestYear == 3,]))
o.pa.meta.bacteria <- meta.bacteria[row.names(o.pa.bacteria),]
t.pa.meta.bacteria <- meta.bacteria[row.names(t.pa.bacteria),]
o.pa.bacteria.stand <- decostand(o.pa.bacteria, method="hellinger")
o.pa.bacteria.dist <- vegdist(o.pa.bacteria.stand, method = "euclidean")
o.pa.bacteria.mds <- metaMDS(o.pa.bacteria.stand, dist = "euclidean")
t.pa.bacteria.stand <- decostand(t.pa.bacteria, method="hellinger")
t.pa.bacteria.dist <- vegdist(t.pa.bacteria.stand, method = "euclidean")
t.pa.bacteria.mds <- metaMDS(t.pa.bacteria.stand, dist = "euclidean")
#t.pa.bacteria.adonis <- adonis(t.pa.bacteria.dist ~ Conduit.D.um. * Conduit.L.m. + avg.top.bottom.density + WCProtPerDM * WLigPerDM * WPPerDM * Cellulose.cell * Hemicellulose.hemic + PlotLocation * Pos * factor(PotID), t.meta.bacteria) 
#o.pa.bacteria.adonis <- adonis(o.pa.bacteria.dist ~ Conduit.D.um. * Conduit.L.m. + avg.top.bottom.density + WCProtPerDM * WLigPerDM * WPPerDM * Cellulose.cell * Hemicellulose.hemic + PlotLocation * Pos * factor(PotID), o.meta.bacteria) 
t.pa.ef.bacteria <- envfit(t.pa.bacteria.mds ~ Conduit.D.um. * Conduit.L.m. + avg.top.bottom.density + WCProtPerDM * WLigPerDM * WPPerDM * Cellulose.cell * Hemicellulose.hemic + PlotLocation * Pos * factor(PotID), t.pa.meta.bacteria)   
o.pa.ef.bacteria <- envfit(o.pa.bacteria.mds ~ Conduit.D.um. * Conduit.L.m. + avg.top.bottom.density + WCProtPerDM * WLigPerDM * WPPerDM * Cellulose.cell * Hemicellulose.hemic + PlotLocation * Pos * factor(PotID), o.pa.meta.bacteria)   

ordiplot(t.bacteria.mds, type="t", display="sites")
ordiplot(t.bacteria.mds, type="p", display="sites")
plot(t.ef.m.bacteria)
plot(t.bacteria.mds, display = "sites", type = "p")
with(t.meta.bacteria, ordiellipse(t.bacteria.mds, Clade, conf = 0.999, label = TRUE))
plot(o.species.mds, display = "sites", type = "p")
with(o.meta.species, ordiellipse(o.species.mds, Clade, conf = 0.95))

#fungi - bacteria mantel test
o.fungi.bac <- o.fungi[row.names(o.bacteria) %in% row.names(o.fungi),]
t.fungi.bac.merge <- merge(t.bacteria, t.fungi, by='row.names')
t.bac.fungi <- t.fungi.bac.merge[,2:1062]
row.names(t.bac.fungi) <- t.fungi.bac.merge$Row.names
t.fungi.bac <- t.fungi.bac.merge[,1063:ncol(t.fungi.bac.merge)]
row.names(t.fungi.bac) <- t.fungi.bac.merge$Row.names
o.fungi.bac.stand <- decostand(o.fungi.bac, method="hellinger")
t.fungi.bac.stand <- decostand(t.fungi.bac, method="hellinger")
t.bac.fungi.stand <- decostand(t.bac.fungi, method="hellinger")
o.fungi.bac.dist <- vegdist(o.fungi.bac.stand, method="euclidean")
t.fungi.bac.dist <- vegdist(t.fungi.bac.stand, method="euclidean")
t.bac.fungi.dist <- vegdist(t.bac.fungi.stand, method="euclidean")
o.fungi.bac.mantel <- mantel(o.fungi.bac.dist, o.bacteria.dist)
t.fungi.bac.mantel <- mantel(t.fungi.bac.dist, t.bac.fungi.dist)


##to do
#envfit with new models *
#try with preseance absence and counts *
#mantel test *
#stacked bar plots
#figures for ESA - AMY
#current vs. initial trait data
#adonis *