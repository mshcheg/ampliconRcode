library(vegan)

species <- read.table("/home/mariya/documents/ITS/ITS.species.txt", header=T,sep="\t")
bacteria.genus <- read.table("~/documents/16S/bacteria.genus.counts", header=T, sep="\t")
#genus <- read.table("/home/mshcheglovitova/AmpliconSeq2/Data/DataParse2/Zanne-ITS/abundance/ITS.genus.txt", header=T, sep="\t")
row.names(species) <- species$taxon
row.names(bacteria.genus) <- bacteria.genus$Taxon_Name
#row.names(genus) <- genus$taxon
#genus_comm <- genus[,3:ncol(genus)]
species_comm <- species[,3:ncol(species)]
bacteria_comm <- bacteria.genus[,4:ncol(bacteria.genus)]
#genus.t <- t(genus_comm)
species.t <- t(species_comm)
bacteria.t <- t(bacteria_comm)
#genus.tyson <- genus.t[grep("...z", rownames(genus.t)),]
species.tyson <- species.t[grep("...z", rownames(species.t)),]
#genus.tyson <- genus.tyson[,apply(genus.tyson,2,sum) > 0]
species.tyson <- species.tyson[,apply(species.tyson,2,sum) > 0]
species.tyson <- species.tyson[apply(species.tyson,1,sum) >= 100,]
bacteria.tyson <- bacteria.t[apply(bacteria.t,1,sum) >= 100,]
#genus.tyson <- genus.tyson[apply(genus.tyson,1,sum) >= 100,]
species.stand <- decostand(species.tyson, method="hellinger")
bacteria.stand <- decostand(bacteria.tyson, method="hellinger")
#genus.stand <- decostand(genus.tyson, method="hellinger")
species.mds <- metaMDS(species.stand, dist = "euclidean")
#genus.mds <- metaMDS(genus.stand, dist = "euclidean")


ordiplot(species.mds, display="sites", type="text")
#ordiplot(genus.mds, display="sites", type="text")

meta <- read.table("~/documents/ITS/ScriptJunk/meta_full.txt", header=T, row.names=1)
meta.species <- meta[row.names(species.tyson),]
#meta.genus <- meta[row.names(genus.tyson),]

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

meta.species <- meta.rep[row.names(species.tyson),]
meta.bacteria <- meta.rep[row.names(bacteria.tyson),]

ordiplot(species.mds)
plot(envfit(species.mds, meta.species[, c(1,2,6,7)]))

#Plot Ordination by Tree Species
tree.sp <- species.tyson
row.names(tree.sp) <- meta.species$Name
tree.sp <- aggregate(tree.sp, by=list(rownames(tree.sp)), FUN=sum)
row.names(tree.sp) <- c("ACRU","AEGL","AMAR","ASTR","CATO","CEOC","CEOC2","COFL","DIVI","FRAM","GLTR","JUNI","JUVI","JUVI2","LOMA","PIEC","PIST","PLOC","PRSE","QUAL","QUVE","QUVE2","ULRU","VIVU")
tree.stand <- decostand(tree.sp[,2:ncol(tree.sp)], method="hellinger")
tree.mds <- metaMDS(tree.stand, dist = "euclidean")
ordiplot(tree.mds, display="sites", type="text")
meta.tree <- meta[!duplicated(meta$Name),]
rownames(meta.tree) <- meta.tree$Name
meta.tree <- meta.tree[row.names(tree.sp),]
tree.dist <- vegdist(tree.stand, method = "euclidean")
tree.clust <- hclust(tree.dist, method = "average")
plot(tree.clust)

tree.fig <- ordiplot(tree.mds, type="none")
points(tree.fig, "sites", pch=19, col="blue", select=meta.tree$Clade == "Angiosperms")
points(tree.fig, "sites", pch=19, col="green", select=meta.tree$Clade == "Gymnosperm")
#ordiellipse(tree.mds, meta.tree$Clade, conf = 0.95, label = TRUE)
plot(envfit(tree.mds, meta.tree[, 14:ncol(meta.tree)], na.rm=T), col="black")

tree.fig <- ordiplot(tree.mds, type="none")
points(tree.fig, "sites", pch=19, col=colors()[375], select=meta.tree$Order == "Cornales")
points(tree.fig, "sites", pch=19, col=colors()[400], select=meta.tree$Order == "Dipsacales")
points(tree.fig, "sites", pch=19, col=colors()[425], select=meta.tree$Order == "Ericales")
points(tree.fig, "sites", pch=19, col=colors()[450], select=meta.tree$Order == "Fabales")
points(tree.fig, "sites", pch=19, col=colors()[475], select=meta.tree$Order == "Fagales")
points(tree.fig, "sites", pch=19, col=colors()[500], select=meta.tree$Order == "Lamiales")
points(tree.fig, "sites", pch=19, col=colors()[525], select=meta.tree$Order == "Magnoliales")
points(tree.fig, "sites", pch=19, col=colors()[550], select=meta.tree$Order == "Pinales")
points(tree.fig, "sites", pch=19, col=colors()[575], select=meta.tree$Order == "Proteales")
points(tree.fig, "sites", pch=19, col=colors()[600], select=meta.tree$Order == "Rosales")
points(tree.fig, "sites", pch=19, col=colors()[625], select=meta.tree$Order == "Sapindales")
points(tree.fig, "sites", pch=19, col=colors()[650], select=meta.tree$Order == "Vitales")
#ordiellipse(tree.mds, meta.tree$Order, conf = 0.95, label = TRUE)
ordicluster(tree.mds, tree.clust, col = "gray")
legend("topright", pch= 19, col= c(colors()[375], colors()[400], colors()[425], colors()[450], colors()[475], colors()[500], colors()[525], colors()[550], colors()[575], colors()[600], colors()[625], colors()[650]),  legend=c("Cornales", "Dipsacales", "Ericales", "Fabales", "Fagales", "Lamiales", "Magnoliales", "Pinales", "Proteales", "Rosales", "Sapindales", "Vitales"))


species.fig <- ordiplot(species.mds, type="none")
points(species.fig, "sites", pch=19, col=colors()[375], select=meta.species$Order == "Cornales")
points(species.fig, "sites", pch=19, col=colors()[400], select=meta.species$Order == "Dipsacales")
points(species.fig, "sites", pch=19, col=colors()[425], select=meta.species$Order == "Ericales")
points(species.fig, "sites", pch=19, col=colors()[450], select=meta.species$Order == "Fabales")
points(species.fig, "sites", pch=19, col=colors()[475], select=meta.species$Order == "Fagales")
points(species.fig, "sites", pch=19, col=colors()[500], select=meta.species$Order == "Lamiales")
points(species.fig, "sites", pch=19, col=colors()[525], select=meta.species$Order == "Magnoliales")
points(species.fig, "sites", pch=19, col=colors()[550], select=meta.species$Order == "Pinales")
points(species.fig, "sites", pch=19, col=colors()[575], select=meta.species$Order == "Proteales")
points(species.fig, "sites", pch=19, col=colors()[600], select=meta.species$Order == "Rosales")
points(species.fig, "sites", pch=19, col=colors()[625], select=meta.species$Order == "Sapindales")
points(species.fig, "sites", pch=19, col=colors()[650], select=meta.species$Order == "Vitales")
#ordiellipse(species.mds, meta.species$Order, conf = 0.95, label = TRUE)
legend("topright", pch= 19, col= c(colors()[375], colors()[400], colors()[425], colors()[450], colors()[475], colors()[500], colors()[525], colors()[550], colors()[575], colors()[600], colors()[625], colors()[650]),  legend=c("Cornales", "Dipsacales", "Ericales", "Fabales", "Fagales", "Lamiales", "Magnoliales", "Pinales", "Proteales", "Rosales", "Sapindales", "Vitales"))



#par(mfrow=c(1,2))
species.fig <- ordiplot(species.mds, type="none", main="Species")
points(species.fig, "sites", pch=19, col="green", select=meta.species$PlotLocation == "L")
points(species.fig, "sites", pch=19, col="blue", select=meta.species$PlotLocation == "H")
#genus.fig <- ordiplot(genus.mds, type="none", main="Genus")
#points(genus.fig, "sites", pch=19, col="green", select=meta.genus$PlotLocation == "L")
#points(genus.fig, "sites", pch=19, col="blue", select=meta.genus$PlotLocation == "H")

**
  par(mfrow=c(1,2))
species.fig <- ordiplot(species.mds, type="none", main="Species")
points(species.fig, "sites", pch=19, col="green", select=meta.species$HarvestYear == "1")
points(species.fig, "sites", pch=19, col="blue", select=meta.species$HarvestYear == "3")
ordiellipse(species.mds, meta.species$HarvestYear, conf = 0.95, label = TRUE)
genus.fig <- ordiplot(genus.mds, type="none", main="Genus")
points(genus.fig, "sites", pch=19, col="green", select=meta.genus$HarvestYear == "1")
points(genus.fig, "sites", pch=19, col="blue", select=meta.genus$HarvestYear == "3")
ordiellipse(genus.mds, meta.genus$HarvestYear, conf = 0.95, label = TRUE)

par(mfrow=c(1,2))
species.fig <- ordiplot(species.mds, type="none", main="Species")
points(species.fig, "sites", pch=19, col="green", select=meta.species$Pos == "t")
points(species.fig, "sites", pch=19, col="blue", select=meta.species$Pos == "b")
genus.fig <- ordiplot(genus.mds, type="none", main="Genus")
points(genus.fig, "sites", pch=19, col="green", select=meta.genus$Pos == "t")
points(genus.fig, "sites", pch=19, col="blue", select=meta.genus$Pos == "b")

**
  par(mfrow=c(1,2))
species.fig <- ordiplot(species.mds, type="none", main="Species")
points(species.fig, "sites", pch=19, col="blue", select=meta.species$Clade == "Angiosperms")
points(species.fig, "sites", pch=19, col="green", select=meta.species$Clade == "Gymnosperm")
ordiellipse(species.mds, meta.species$Clade, conf = 0.95, label = TRUE)
genus.fig <- ordiplot(genus.mds, type="none", main="Genus")
points(genus.fig, "sites", pch=19, col="blue", select=meta.genus$Clade == "Angiosperms")
points(genus.fig, "sites", pch=19, col="green", select=meta.genus$Clade == "Gymnosperm")
ordiellipse(genus.mds, meta.genus$Clade, conf = 0.95, label = TRUE)

#Lo.genus <- subset(genus.tyson,rownames(genus.tyson) %in% rownames(meta[meta$PlotLocation == "L",]))
#Hi.genus <- subset(genus.tyson,rownames(genus.tyson) %in% rownames(meta[meta$PlotLocation == "H",]))
Lo.species <- subset(species.tyson,rownames(species.tyson) %in% rownames(meta[meta$PlotLocation == "L",]))
Hi.species <- subset(species.tyson,rownames(species.tyson) %in% rownames(meta[meta$PlotLocation == "H",]))
Hi.species.stand <- decostand(Hi.species, method="hellinger")
#Hi.genus.stand <- decostand(Hi.genus, method="hellinger")
Hi.species.mds <- metaMDS(Hi.species.stand, dist = "euclidean")
#Hi.genus.mds <- metaMDS(Hi.genus.stand, dist = "euclidean")
Lo.species.stand <- decostand(Lo.species, method="hellinger")
#Lo.genus.stand <- decostand(Lo.genus, method="hellinger")
Lo.species.mds <- metaMDS(Lo.species.stand, dist = "euclidean")
#Lo.genus.mds <- metaMDS(Lo.genus.stand, dist = "euclidean")
Hi.meta.species <- meta.species.rep[row.names(Hi.species),]
#Lo.meta.genus <- meta[row.names(Lo.genus),]
#Hi.meta.genus <- meta[row.names(Hi.genus),]
Lo.meta.species <- meta.species.rep[row.names(Lo.species),]


t.Lo.genus <- subset(Lo.genus,rownames(Lo.genus) %in% rownames(meta[meta$PlotLocation == "L" & meta$Pos == "t",]))
t.Hi.genus <- subset(Hi.genus,rownames(Hi.genus) %in% rownames(meta[meta$PlotLocation == "H" & meta$Pos == "t",]))
b.Lo.genus <- subset(Lo.genus,rownames(Lo.genus) %in% rownames(meta[meta$PlotLocation == "L" & meta$Pos == "b",]))
b.Hi.genus <- subset(Hi.genus,rownames(Hi.genus) %in% rownames(meta[meta$PlotLocation == "H" & meta$Pos == "b",]))
t.Hi.genus.stand <- decostand(t.Hi.genus, method="hellinger")
b.Hi.genus.stand <- decostand(b.Hi.genus, method="hellinger")
t.Lo.genus.stand <- decostand(t.Lo.genus, method="hellinger")
b.Lo.genus.stand <- decostand(b.Lo.genus, method="hellinger")
t.Hi.genus.mds <- metaMDS(t.Hi.genus.stand, dist = "euclidean")
b.Hi.genus.mds <- metaMDS(b.Hi.genus.stand, dist = "euclidean")
t.Lo.genus.mds <- metaMDS(t.Lo.genus.stand, dist = "euclidean")
b.Lo.genus.mds <- metaMDS(b.Lo.genus.stand, dist = "euclidean")
t.Lo.meta.genus <- meta[row.names(t.Lo.genus),]
b.Lo.meta.genus <- meta[row.names(b.Lo.genus),]
t.Hi.meta.genus <- meta[row.names(t.Hi.genus),]
b.Hi.meta.genus <- meta[row.names(b.Hi.genus),]


t.Lo.species <- subset(Lo.species,rownames(Lo.species) %in% rownames(meta[meta$PlotLocation == "L" & meta$Pos == "t",]))
t.Hi.species <- subset(Hi.species,rownames(Hi.species) %in% rownames(meta[meta$PlotLocation == "H" & meta$Pos == "t",]))
b.Lo.species <- subset(Lo.species,rownames(Lo.species) %in% rownames(meta[meta$PlotLocation == "L" & meta$Pos == "b",]))
b.Hi.species <- subset(Hi.species,rownames(Hi.species) %in% rownames(meta[meta$PlotLocation == "H" & meta$Pos == "b",]))
t.Hi.species.stand <- decostand(t.Hi.species, method="hellinger")
b.Hi.species.stand <- decostand(b.Hi.species, method="hellinger")
t.Lo.species.stand <- decostand(t.Lo.species, method="hellinger")
b.Lo.species.stand <- decostand(b.Lo.species, method="hellinger")
t.Hi.species.mds <- metaMDS(t.Hi.species.stand, dist = "euclidean")
b.Hi.species.mds <- metaMDS(b.Hi.species.stand, dist = "euclidean")
t.Lo.species.mds <- metaMDS(t.Lo.species.stand, dist = "euclidean")
b.Lo.species.mds <- metaMDS(b.Lo.species.stand, dist = "euclidean")
t.Lo.meta.species <- meta.species.rep[row.names(t.Lo.species),]
b.Lo.meta.species <- meta.species.rep[row.names(b.Lo.species),]
t.Hi.meta.species <- meta.species.rep[row.names(t.Hi.species),]
b.Hi.meta.species <- meta.species.rep[row.names(b.Hi.species),]

o.t.Lo.genus <- subset(Lo.genus,rownames(Lo.genus) %in% rownames(meta[meta$PlotLocation == "L" & meta$Pos == "t" & meta$HarvestYear == 1,]))
o.t.Hi.genus <- subset(Hi.genus,rownames(Hi.genus) %in% rownames(meta[meta$PlotLocation == "H" & meta$Pos == "t" & meta$HarvestYear == 1,]))
o.b.Lo.genus <- subset(Lo.genus,rownames(Lo.genus) %in% rownames(meta[meta$PlotLocation == "L" & meta$Pos == "b" & meta$HarvestYear == 1,]))
o.b.Hi.genus <- subset(Hi.genus,rownames(Hi.genus) %in% rownames(meta[meta$PlotLocation == "H" & meta$Pos == "b" & meta$HarvestYear == 1,]))
o.t.Hi.genus.stand <- decostand(o.t.Hi.genus, method="hellinger")
o.b.Hi.genus.stand <- decostand(o.b.Hi.genus, method="hellinger")
o.t.Lo.genus.stand <- decostand(o.t.Lo.genus, method="hellinger")
o.b.Lo.genus.stand <- decostand(o.b.Lo.genus, method="hellinger")
o.t.Hi.genus.mds <- metaMDS(o.t.Hi.genus.stand, dist = "euclidean")
o.b.Hi.genus.mds <- metaMDS(o.b.Hi.genus.stand, dist = "euclidean")
o.t.Lo.genus.mds <- metaMDS(o.t.Lo.genus.stand, dist = "euclidean")
o.b.Lo.genus.mds <- metaMDS(o.b.Lo.genus.stand, dist = "euclidean")
o.t.Lo.meta.genus <- meta[row.names(o.t.Lo.genus),]
o.b.Lo.meta.genus <- meta[row.names(o.b.Lo.genus),]
o.t.Hi.meta.genus <- meta[row.names(o.t.Hi.genus),]
o.b.Hi.meta.genus <- meta[row.names(o.b.Hi.genus),]

t.t.Lo.genus <- subset(Lo.genus,rownames(Lo.genus) %in% rownames(meta[meta$PlotLocation == "L" & meta$Pos == "t" & meta$HarvestYear == 3,]))
t.t.Hi.genus <- subset(t.Hi.genus,rownames(t.Hi.genus) %in% rownames(meta[meta$PlotLocation == "H" & meta$Pos == "t" & meta$HarvestYear == 3,]))
t.b.Lo.genus <- subset(b.Lo.genus,rownames(b.Lo.genus) %in% rownames(meta[meta$PlotLocation == "L" & meta$Pos == "b" & meta$HarvestYear == 3,]))
t.b.Hi.genus <- subset(b.Hi.genus,rownames(b.Hi.genus) %in% rownames(meta[meta$PlotLocation == "H" & meta$Pos == "b" & meta$HarvestYear == 3,]))
t.t.Hi.genus.stand <- decostand(t.t.Hi.genus, method="hellinger")
t.b.Hi.genus.stand <- decostand(t.b.Hi.genus, method="hellinger")
t.t.Lo.genus.stand <- decostand(t.t.Lo.genus, method="hellinger")
t.b.Lo.genus.stand <- decostand(t.b.Lo.genus, method="hellinger")
t.t.Hi.genus.mds <- metaMDS(t.t.Hi.genus.stand, dist = "euclidean")
t.b.Hi.genus.mds <- metaMDS(t.b.Hi.genus.stand, dist = "euclidean")
t.t.Lo.genus.mds <- metaMDS(t.t.Lo.genus.stand, dist = "euclidean")
t.b.Lo.genus.mds <- metaMDS(t.b.Lo.genus.stand, dist = "euclidean")
t.t.Lo.meta.genus <- meta[row.names(t.t.Lo.genus),]
t.b.Lo.meta.genus <- meta[row.names(t.b.Lo.genus),]
t.t.Hi.meta.genus <- meta[row.names(t.t.Hi.genus),]
t.b.Hi.meta.genus <- meta[row.names(t.b.Hi.genus),]

o.t.Lo.species <- subset(Lo.species,rownames(Lo.species) %in% rownames(meta[meta$PlotLocation == "L" & meta$Pos == "t" & meta$HarvestYear == 1,]))
o.t.Hi.species <- subset(Hi.species,rownames(Hi.species) %in% rownames(meta[meta$PlotLocation == "H" & meta$Pos == "t" & meta$HarvestYear == 1,]))
o.b.Lo.species <- subset(Lo.species,rownames(Lo.species) %in% rownames(meta[meta$PlotLocation == "L" & meta$Pos == "b" & meta$HarvestYear == 1,]))
o.b.Hi.species <- subset(Hi.species,rownames(Hi.species) %in% rownames(meta[meta$PlotLocation == "H" & meta$Pos == "b" & meta$HarvestYear == 1,]))
o.t.Hi.species.stand <- decostand(o.t.Hi.species, method="hellinger")
o.b.Hi.species.stand <- decostand(o.b.Hi.species, method="hellinger")
o.t.Lo.species.stand <- decostand(o.t.Lo.species, method="hellinger")
o.b.Lo.species.stand <- decostand(o.b.Lo.species, method="hellinger")
o.t.Lo.meta.species <- meta.species.rep[row.names(o.t.Lo.species),]
o.b.Lo.meta.species <- meta.species.rep[row.names(o.b.Lo.species),]
o.t.Hi.meta.species <- meta.species.rep[row.names(o.t.Hi.species),]
o.b.Hi.meta.species <- meta.species.rep[row.names(o.b.Hi.species),]

t.t.Lo.species <- subset(Lo.species,rownames(Lo.species) %in% rownames(meta[meta$PlotLocation == "L" & meta$Pos == "t" & meta$HarvestYear == 3,]))
t.t.Hi.species <- subset(Hi.species,rownames(Hi.species) %in% rownames(meta[meta$PlotLocation == "H" & meta$Pos == "t" & meta$HarvestYear == 3,]))
t.b.Lo.species <- subset(Lo.species,rownames(Lo.species) %in% rownames(meta[meta$PlotLocation == "L" & meta$Pos == "b" & meta$HarvestYear == 3,]))
t.b.Hi.species <- subset(Hi.species,rownames(Hi.species) %in% rownames(meta[meta$PlotLocation == "H" & meta$Pos == "b" & meta$HarvestYear == 3,]))
t.t.Hi.species.stand <- decostand(t.t.Hi.species, method="hellinger")
t.b.Hi.species.stand <- decostand(t.b.Hi.species, method="hellinger")
t.t.Lo.species.stand <- decostand(t.t.Lo.species, method="hellinger")
t.b.Lo.species.stand <- decostand(t.b.Lo.species, method="hellinger")

t.b.Lo.species.dist <- vegdist(t.b.Lo.species.stand, method = "euclidean")
t.b.Lo.species.Cellulose.cell <- vegdist(t.b.Lo.meta.species$Cellulose.cell, method="euclidean")
mantel(t.b.Lo.species.dist, t.b.Lo.species.Cellulose.cell)

t.t.Hi.species.mds <- metaMDS(t.t.Hi.species.stand, dist = "euclidean",4)
t.b.Hi.species.mds <- metaMDS(t.b.Hi.species.stand, dist = "euclidean",4)
t.t.Lo.species.mds <- metaMDS(t.t.Lo.species.stand, dist = "euclidean",4)
t.b.Lo.species.mds <- metaMDS(t.b.Lo.species.stand, dist = "euclidean",4)
t.t.Lo.meta.species <- meta.species.rep[row.names(t.t.Lo.species),]
t.b.Lo.meta.species <- meta.species.rep[row.names(t.b.Lo.species),]
t.t.Hi.meta.species <- meta.species.rep[row.names(t.t.Hi.species),]
t.b.Hi.meta.species <- meta.species.rep[row.names(t.b.Hi.species),]

**
  par(mfrow=c(1,2))
Lo.species.fig <- ordiplot(Lo.species.mds, type="none", main="Low Species")
points(Lo.species.fig, "sites", pch=19, col="blue", select=Lo.meta.species$Clade == "Angiosperms")
points(Lo.species.fig, "sites", pch=19, col="green", select=Lo.meta.species$Clade == "Gymnosperm")
plot(envfit(Lo.species.fig, Lo.meta.species[, 14:ncol(Lo.meta.species)], na.rm=T), add=T)
#ordiellipse(Lo.species.mds, Lo.meta.species$Clade, conf = 0.95, label = TRUE)
Hi.species.fig <- ordiplot(Hi.species.mds, type="none", main="High Species")
points(Hi.species.fig, "sites", pch=19, col="blue", select=Hi.meta.species$Clade == "Angiosperms")
points(Hi.species.fig, "sites", pch=19, col="green", select=Hi.meta.species$Clade == "Gymnosperm")
plot(envfit(Hi.species.fig, Hi.meta.species[, 14:ncol(Hi.meta.species)], na.rm=T), add=T)
#ordiellipse(Hi.species.mds, Hi.meta.species$Clade, conf = 0.95, label = TRUE)

par(mfrow=c(1,2))
Lo.genus.fig <- ordiplot(Lo.genus.mds, type="none", main="Low genus")
points(Lo.genus.fig, "sites", pch=19, col="blue", select=Lo.meta.genus$Clade == "Angiosperms")
points(Lo.genus.fig, "sites", pch=19, col="green", select=Lo.meta.genus$Clade == "Gymnosperm")
ordiellipse(Lo.genus.mds, Lo.meta.genus$Clade, conf = 0.95, label = TRUE)
Hi.genus.fig <- ordiplot(Hi.genus.mds, type="none", main="High genus")
points(Hi.genus.fig, "sites", pch=19, col="blue", select=Hi.meta.genus$Clade == "Angiosperms")
points(Hi.genus.fig, "sites", pch=19, col="green", select=Hi.meta.genus$Clade == "Gymnosperm")
ordiellipse(Hi.genus.mds, Hi.meta.genus$Clade, conf = 0.95, label = TRUE)

**
  par(mfrow=c(1,2))
Lo.species.fig <- ordiplot(Lo.species.mds, type="none", main="Low Species")
points(Lo.species.fig, "sites", pch=19, col="blue", select=Lo.meta.species$HarvestYear == "1")
points(Lo.species.fig, "sites", pch=19, col="green", select=Lo.meta.species$HarvestYear == "3")
ordiellipse(Lo.species.mds, Lo.meta.species$HarvestYear, conf = 0.95, label = TRUE)
Hi.species.fig <- ordiplot(Hi.species.mds, type="none", main="High Species")
points(Hi.species.fig, "sites", pch=19, col="blue", select=Hi.meta.species$HarvestYear == "1")
points(Hi.species.fig, "sites", pch=19, col="green", select=Hi.meta.species$HarvestYear == "3")
ordiellipse(Hi.species.mds, Hi.meta.species$HarvestYear, conf = 0.95, label = TRUE)

**
  par(mfrow=c(1,2))
Lo.genus.fig <- ordiplot(Lo.genus.mds, type="none", main="Low genus")
points(Lo.genus.fig, "sites", pch=19, col="blue", select=Lo.meta.genus$HarvestYear == "1")
points(Lo.genus.fig, "sites", pch=19, col="green", select=Lo.meta.genus$HarvestYear == "3")
ordiellipse(Lo.genus.mds, Lo.meta.genus$HarvestYear, conf = 0.95, label = TRUE)
Hi.genus.fig <- ordiplot(Hi.genus.mds, type="none", main="High genus")
points(Hi.genus.fig, "sites", pch=19, col="blue", select=Hi.meta.genus$HarvestYear == "1")
points(Hi.genus.fig, "sites", pch=19, col="green", select=Hi.meta.genus$HarvestYear == "3")
ordiellipse(Hi.genus.mds, Hi.meta.genus$HarvestYear, conf = 0.95, label = TRUE)

par(mfrow=c(1,2))
Lo.species.fig <- ordiplot(Lo.species.mds, type="none", main="Low Species")
points(Lo.species.fig, "sites", pch=19, col="blue", select=Lo.meta.species$Pos == "t")
points(Lo.species.fig, "sites", pch=19, col="green", select=Lo.meta.species$Pos == "b")
ordiellipse(Lo.species.mds, Lo.meta.species$Pos, conf = 0.95, label = TRUE)
Hi.species.fig <- ordiplot(Hi.species.mds, type="none", main="High Species")
points(Hi.species.fig, "sites", pch=19, col="blue", select=Hi.meta.species$Pos == "t")
points(Hi.species.fig, "sites", pch=19, col="green", select=Hi.meta.species$Pos == "b")
ordiellipse(Hi.species.mds, Hi.meta.species$Pos, conf = 0.95, label = TRUE)

par(mfrow=c(1,2))
Lo.genus.fig <- ordiplot(Lo.genus.mds, type="none", main="Low genus")
points(Lo.genus.fig, "sites", pch=19, col="blue", select=Lo.meta.genus$Pos == "t")
points(Lo.genus.fig, "sites", pch=19, col="green", select=Lo.meta.genus$Pos == "b")
ordiellipse(Lo.genus.mds, Lo.meta.genus$Pos, conf = 0.95, label = TRUE)
Hi.genus.fig <- ordiplot(Hi.genus.mds, type="none", main="High genus")
points(Hi.genus.fig, "sites", pch=19, col="blue", select=Hi.meta.genus$Pos == "t")
points(Hi.genus.fig, "sites", pch=19, col="green", select=Hi.meta.genus$Pos == "b")
ordiellipse(Hi.genus.mds, Hi.meta.genus$Pos, conf = 0.95, label = TRUE)

par(mfrow=c(2,2))
t.Lo.species.fig <- ordiplot(t.Lo.species.mds, type="none", main="Low top")
points(t.Lo.species.fig, "sites", pch=19, col=colors()[375], select=t.Lo.meta.species$Order == "Cornales")
points(t.Lo.species.fig, "sites", pch=19, col=colors()[400], select=t.Lo.meta.species$Order == "Dipsacales")
points(t.Lo.species.fig, "sites", pch=19, col=colors()[425], select=t.Lo.meta.species$Order == "Ericales")
points(t.Lo.species.fig, "sites", pch=19, col=colors()[450], select=t.Lo.meta.species$Order == "Fabales")
points(t.Lo.species.fig, "sites", pch=19, col=colors()[475], select=t.Lo.meta.species$Order == "Fagales")
points(t.Lo.species.fig, "sites", pch=19, col=colors()[500], select=t.Lo.meta.species$Order == "Lamiales")
points(t.Lo.species.fig, "sites", pch=19, col=colors()[525], select=t.Lo.meta.species$Order == "Magnoliales")
points(t.Lo.species.fig, "sites", pch=19, col=colors()[550], select=t.Lo.meta.species$Order == "Pinales")
points(t.Lo.species.fig, "sites", pch=19, col=colors()[575], select=t.Lo.meta.species$Order == "Proteales")
points(t.Lo.species.fig, "sites", pch=19, col=colors()[600], select=t.Lo.meta.species$Order == "Rosales")
points(t.Lo.species.fig, "sites", pch=19, col=colors()[625], select=t.Lo.meta.species$Order == "Sapindales")
points(t.Lo.species.fig, "sites", pch=19, col=colors()[650], select=t.Lo.meta.species$Order == "Vitales")
ordiellipse(t.Lo.species.mds, t.Lo.meta.species$Order, conf = 0.95, label = TRUE, show.groups="Rosales", col=colors()[600])
ordiellipse(t.Lo.species.mds, t.Lo.meta.species$Order, conf = 0.95, label = TRUE, show.groups="Fabales", col=colors()[450])
ordiellipse(t.Lo.species.mds, t.Lo.meta.species$Order, conf = 0.95, label = TRUE, show.groups="Pinales", col=colors()[550])
ordiellipse(t.Lo.species.mds, t.Lo.meta.species$Order, conf = 0.95, label = TRUE, show.groups="Sapindales", col=colors()[625])
#legend("topright", pch= 19, col= c(colors()[375], colors()[400], colors()[425], colors()[450], colors()[475], colors()[500], colors()[525], colors()[550], colors()[575], colors()[600], colors()[625], colors()[650]),  legend=c("Cornales", "Dipsacales", "Ericales", "Fabales", "Fagales", "Lamiales", "Magnoliales", "Pinales", "Proteales", "Rosales", "Sapindales", "Vitales"))
t.Hi.species.fig <- ordiplot(t.Hi.species.mds, type="none", main="High top")
points(t.Hi.species.fig, "sites", pch=19, col=colors()[375], select=t.Hi.meta.species$Order == "Cornales")
points(t.Hi.species.fig, "sites", pch=19, col=colors()[400], select=t.Hi.meta.species$Order == "Dipsacales")
points(t.Hi.species.fig, "sites", pch=19, col=colors()[425], select=t.Hi.meta.species$Order == "Ericales")
points(t.Hi.species.fig, "sites", pch=19, col=colors()[450], select=t.Hi.meta.species$Order == "Fabales")
points(t.Hi.species.fig, "sites", pch=19, col=colors()[475], select=t.Hi.meta.species$Order == "Fagales")
points(t.Hi.species.fig, "sites", pch=19, col=colors()[500], select=t.Hi.meta.species$Order == "Lamiales")
points(t.Hi.species.fig, "sites", pch=19, col=colors()[525], select=t.Hi.meta.species$Order == "Magnoliales")
points(t.Hi.species.fig, "sites", pch=19, col=colors()[550], select=t.Hi.meta.species$Order == "Pinales")
points(t.Hi.species.fig, "sites", pch=19, col=colors()[575], select=t.Hi.meta.species$Order == "Proteales")
points(t.Hi.species.fig, "sites", pch=19, col=colors()[600], select=t.Hi.meta.species$Order == "Rosales")
points(t.Hi.species.fig, "sites", pch=19, col=colors()[625], select=t.Hi.meta.species$Order == "Sapindales")
points(t.Hi.species.fig, "sites", pch=19, col=colors()[650], select=t.Hi.meta.species$Order == "Vitales")
ordiellipse(t.Hi.species.mds, t.Hi.meta.species$Order, conf = 0.95, label = TRUE, show.groups="Rosales", col=colors()[600])
ordiellipse(t.Hi.species.mds, t.Hi.meta.species$Order, conf = 0.95, label = TRUE, show.groups="Fabales", col=colors()[450])
ordiellipse(t.Hi.species.mds, t.Hi.meta.species$Order, conf = 0.95, label = TRUE, show.groups="Pinales", col=colors()[550])
ordiellipse(t.Hi.species.mds, t.Hi.meta.species$Order, conf = 0.95, label = TRUE, show.groups="Sapindales", col=colors()[625])
#legend("topright", pch= 19, col= c(colors()[375], colors()[400], colors()[425], colors()[450], colors()[475], colors()[500], colors()[525], colors()[550], colors()[575], colors()[600], colors()[625], colors()[650]),  legend=c("Cornales", "Dipsacales", "Ericales", "Fabales", "Fagales", "Lamiales", "Magnoliales", "Pinales", "Proteales", "Rosales", "Sapindales", "Vitales"))
b.Lo.species.fig <- ordiplot(b.Lo.species.mds, type="none", main="Low bottom")
points(b.Lo.species.fig, "sites", pch=19, col=colors()[375], select=b.Lo.meta.species$Order == "Cornales")
points(b.Lo.species.fig, "sites", pch=19, col=colors()[400], select=b.Lo.meta.species$Order == "Dipsacales")
points(b.Lo.species.fig, "sites", pch=19, col=colors()[425], select=b.Lo.meta.species$Order == "Ericales")
points(b.Lo.species.fig, "sites", pch=19, col=colors()[450], select=b.Lo.meta.species$Order == "Fabales")
points(b.Lo.species.fig, "sites", pch=19, col=colors()[475], select=b.Lo.meta.species$Order == "Fagales")
points(b.Lo.species.fig, "sites", pch=19, col=colors()[500], select=b.Lo.meta.species$Order == "Lamiales")
points(b.Lo.species.fig, "sites", pch=19, col=colors()[525], select=b.Lo.meta.species$Order == "Magnoliales")
points(b.Lo.species.fig, "sites", pch=19, col=colors()[550], select=b.Lo.meta.species$Order == "Pinales")
points(b.Lo.species.fig, "sites", pch=19, col=colors()[575], select=b.Lo.meta.species$Order == "Proteales")
points(b.Lo.species.fig, "sites", pch=19, col=colors()[600], select=b.Lo.meta.species$Order == "Rosales")
points(b.Lo.species.fig, "sites", pch=19, col=colors()[625], select=b.Lo.meta.species$Order == "Sapindales")
points(b.Lo.species.fig, "sites", pch=19, col=colors()[650], select=b.Lo.meta.species$Order == "Vitales")
ordiellipse(b.Lo.species.mds, b.Lo.meta.species$Order, conf = 0.95, label = TRUE, show.groups="Rosales", col=colors()[600])
ordiellipse(b.Lo.species.mds, b.Lo.meta.species$Order, conf = 0.95, label = TRUE, show.groups="Fabales", col=colors()[450])
ordiellipse(b.Lo.species.mds, b.Lo.meta.species$Order, conf = 0.95, label = TRUE, show.groups="Pinales", col=colors()[550])
ordiellipse(b.Lo.species.mds, b.Lo.meta.species$Order, conf = 0.95, label = TRUE, show.groups="Sapindales", col=colors()[625])
#legend("topright", pch= 19, col= c(colors()[375], colors()[400], colors()[425], colors()[450], colors()[475], colors()[500], colors()[525], colors()[550], colors()[575], colors()[600], colors()[625], colors()[650]),  legend=c("Cornales", "Dipsacales", "Ericales", "Fabales", "Fagales", "Lamiales", "Magnoliales", "Pinales", "Proteales", "Rosales", "Sapindales", "Vitales"))
b.Hi.species.fig <- ordiplot(b.Hi.species.mds, type="none", main="High bottom")
points(b.Hi.species.fig, "sites", pch=19, col=colors()[375], select=b.Hi.meta.species$Order == "Cornales")
points(b.Hi.species.fig, "sites", pch=19, col=colors()[400], select=b.Hi.meta.species$Order == "Dipsacales")
points(b.Hi.species.fig, "sites", pch=19, col=colors()[425], select=b.Hi.meta.species$Order == "Ericales")
points(b.Hi.species.fig, "sites", pch=19, col=colors()[450], select=b.Hi.meta.species$Order == "Fabales")
points(b.Hi.species.fig, "sites", pch=19, col=colors()[475], select=b.Hi.meta.species$Order == "Fagales")
points(b.Hi.species.fig, "sites", pch=19, col=colors()[500], select=b.Hi.meta.species$Order == "Lamiales")
points(b.Hi.species.fig, "sites", pch=19, col=colors()[525], select=b.Hi.meta.species$Order == "Magnoliales")
points(b.Hi.species.fig, "sites", pch=19, col=colors()[550], select=b.Hi.meta.species$Order == "Pinales")
points(b.Hi.species.fig, "sites", pch=19, col=colors()[575], select=b.Hi.meta.species$Order == "Proteales")
points(b.Hi.species.fig, "sites", pch=19, col=colors()[600], select=b.Hi.meta.species$Order == "Rosales")
points(b.Hi.species.fig, "sites", pch=19, col=colors()[625], select=b.Hi.meta.species$Order == "Sapindales")
points(b.Hi.species.fig, "sites", pch=19, col=colors()[650], select=b.Hi.meta.species$Order == "Vitales")
ordiellipse(b.Hi.species.mds, b.Hi.meta.species$Order, conf = 0.95, label = TRUE, show.groups="Rosales", col=colors()[600])
ordiellipse(b.Hi.species.mds, b.Hi.meta.species$Order, conf = 0.95, label = TRUE, show.groups="Fabales", col=colors()[450])
ordiellipse(b.Hi.species.mds, b.Hi.meta.species$Order, conf = 0.95, label = TRUE, show.groups="Pinales", col=colors()[550])
ordiellipse(b.Hi.species.mds, b.Hi.meta.species$Order, conf = 0.95, label = TRUE, show.groups="Sapindales", col=colors()[625])
#legend("topright", pch= 19, col= c(colors()[375], colors()[400], colors()[425], colors()[450], colors()[475], colors()[500], colors()[525], colors()[550], colors()[575], colors()[600], colors()[625], colors()[650]),  legend=c("Cornales", "Dipsacales", "Ericales", "Fabales", "Fagales", "Lamiales", "Magnoliales", "Pinales", "Proteales", "Rosales", "Sapindales", "Vitales"))


pdf("YearOneITSordination.pdf", width=11, height=8.5)
op <- par(mfrow=c(2,2), oma=c(3,0,0,0)) 
#par("usr")
#[1] -2.214074  2.214074 -1.620000  1.620000
o.t.Lo.species.fig <- ordiplot(o.t.Lo.species.mds, type="none", main="Year 1: Low Top", bty = "n", xlab="")
#points(o.t.Lo.species.fig, "sites", pch=16, col=colors()[375], select=o.t.Lo.meta.species$Order == "Cornales")
points(o.t.Lo.species.fig, "sites", pch=16, col=colors()[653], select=o.t.Lo.meta.species$Order == "Dipsacales")
#points(o.t.Lo.species.fig, "sites", pch=16, col=colors()[425], select=o.t.Lo.meta.species$Order == "Ericales")
#points(o.t.Lo.species.fig, "sites", pch=16, col=colors()[450], select=o.t.Lo.meta.species$Order == "Fabales")
points(o.t.Lo.species.fig, "sites", pch=16, col=colors()[475], select=o.t.Lo.meta.species$Order == "Fagales")
#points(o.t.Lo.species.fig, "sites", pch=16, col=colors()[500], select=o.t.Lo.meta.species$Order == "Lamiales")
#points(o.t.Lo.species.fig, "sites", pch=16, col=colors()[525], select=o.t.Lo.meta.species$Order == "Magnoliales")
points(o.t.Lo.species.fig, "sites", pch=16, col=colors()[550], select=o.t.Lo.meta.species$Order == "Pinales")
#points(o.t.Lo.species.fig, "sites", pch=16, col=colors()[575], select=o.t.Lo.meta.species$Order == "Proteales")
points(o.t.Lo.species.fig, "sites", pch=16, col=colors()[600], select=o.t.Lo.meta.species$Order == "Rosales")
#points(o.t.Lo.species.fig, "sites", pch=16, col=colors()[625], select=o.t.Lo.meta.species$Order == "Sapindales")
#points(o.t.Lo.species.fig, "sites", pch=16, col=colors()[650], select=o.t.Lo.meta.species$Order == "Vitales")
ordiellipse(o.t.Lo.species.mds, o.t.Lo.meta.species$Order, conf = 0.95, label = TRUE, show.groups="Rosales", col=colors()[600])
ordiellipse(o.t.Lo.species.mds, o.t.Lo.meta.species$Order, conf = 0.95, label = TRUE, show.groups="Fagales", col=colors()[475])
ordiellipse(o.t.Lo.species.mds, o.t.Lo.meta.species$Order, conf = 0.95, label = TRUE, show.groups="Pinales", col=colors()[550])
o.t.Hi.species.fig <- ordiplot(o.t.Hi.species.mds, type="none", main="Year 1: High Top", bty = "n", xlab="", ylab="")
#points(o.t.Hi.species.fig, "sites", pch=16, col=colors()[375], select=o.t.Hi.meta.species$Order == "Cornales")
points(o.t.Hi.species.fig, "sites", pch=16, col=colors()[653], select=o.t.Hi.meta.species$Order == "Dipsacales")
#points(o.t.Hi.species.fig, "sites", pch=16, col=colors()[425], select=o.t.Hi.meta.species$Order == "Ericales")
#points(o.t.Hi.species.fig, "sites", pch=16, col=colors()[450], select=o.t.Hi.meta.species$Order == "Fabales")
points(o.t.Hi.species.fig, "sites", pch=16, col=colors()[475], select=o.t.Hi.meta.species$Order == "Fagales")
#points(o.t.Hi.species.fig, "sites", pch=16, col=colors()[500], select=o.t.Hi.meta.species$Order == "Lamiales")
#points(o.t.Hi.species.fig, "sites", pch=16, col=colors()[525], select=o.t.Hi.meta.species$Order == "Magnoliales")
points(o.t.Hi.species.fig, "sites", pch=16, col=colors()[550], select=o.t.Hi.meta.species$Order == "Pinales")
#points(o.t.Hi.species.fig, "sites", pch=16, col=colors()[575], select=o.t.Hi.meta.species$Order == "Proteales")
points(o.t.Hi.species.fig, "sites", pch=16, col=colors()[600], select=o.t.Hi.meta.species$Order == "Rosales")
#points(o.t.Hi.species.fig, "sites", pch=16, col=colors()[625], select=o.t.Hi.meta.species$Order == "Sapindales")
#points(o.t.Hi.species.fig, "sites", pch=16, col=colors()[650], select=o.t.Hi.meta.species$Order == "Vitales")
ordiellipse(o.t.Hi.species.mds, o.t.Hi.meta.species$Order, conf = 0.95, label = TRUE, show.groups="Rosales", col=colors()[600])
ordiellipse(o.t.Hi.species.mds, o.t.Hi.meta.species$Order, conf = 0.95, label = TRUE, show.groups="Fagales", col=colors()[475])
ordiellipse(o.t.Hi.species.mds, o.t.Hi.meta.species$Order, conf = 0.95, label = TRUE, show.groups="Pinales", col=colors()[550])
o.b.Lo.species.fig <- ordiplot(o.b.Lo.species.mds, type="none", main="Year 1: Low Bottom", bty = "n")
#points(o.b.Lo.species.fig, "sites", pch=16, col=colors()[375], select=o.b.Lo.meta.species$Order == "Cornales")
points(o.b.Lo.species.fig, "sites", pch=16, col=colors()[653], select=o.b.Lo.meta.species$Order == "Dipsacales")
#points(o.b.Lo.species.fig, "sites", pch=16, col=colors()[425], select=o.b.Lo.meta.species$Order == "Ericales")
#points(o.b.Lo.species.fig, "sites", pch=16, col=colors()[450], select=o.b.Lo.meta.species$Order == "Fabales")
points(o.b.Lo.species.fig, "sites", pch=16, col=colors()[475], select=o.b.Lo.meta.species$Order == "Fagales")
#points(o.b.Lo.species.fig, "sites", pch=16, col=colors()[500], select=o.b.Lo.meta.species$Order == "Lamiales")
#points(o.b.Lo.species.fig, "sites", pch=16, col=colors()[525], select=o.b.Lo.meta.species$Order == "Magnoliales")
points(o.b.Lo.species.fig, "sites", pch=16, col=colors()[550], select=o.b.Lo.meta.species$Order == "Pinales")
#points(o.b.Lo.species.fig, "sites", pch=16, col=colors()[575], select=o.b.Lo.meta.species$Order == "Proteales")
points(o.b.Lo.species.fig, "sites", pch=16, col=colors()[600], select=o.b.Lo.meta.species$Order == "Rosales")
#points(o.b.Lo.species.fig, "sites", pch=16, col=colors()[625], select=o.b.Lo.meta.species$Order == "Sapindales")
#points(o.b.Lo.species.fig, "sites", pch=16, col=colors()[650], select=o.b.Lo.meta.species$Order == "Vitales")
ordiellipse(o.b.Lo.species.mds, o.b.Lo.meta.species$Order, conf = 0.95, label = TRUE, show.groups="Rosales", col=colors()[600])
ordiellipse(o.b.Lo.species.mds, o.b.Lo.meta.species$Order, conf = 0.95, label = TRUE, show.groups="Fagales", col=colors()[475])
ordiellipse(o.b.Lo.species.mds, o.b.Lo.meta.species$Order, conf = 0.95, label = TRUE, show.groups="Pinales", col=colors()[550])
o.b.Hi.species.fig <- ordiplot(o.b.Hi.species.mds, type="none", main="Year 1: High Bottom", bty = "n", ylab="")
#points(o.b.Hi.species.fig, "sites", pch=16, col=colors()[375], select=o.b.Hi.meta.species$Order == "Cornales")
points(o.b.Hi.species.fig, "sites", pch=16, col=colors()[653], select=o.b.Hi.meta.species$Order == "Dipsacales")
#points(o.b.Hi.species.fig, "sites", pch=16, col=colors()[425], select=o.b.Hi.meta.species$Order == "Ericales")
#points(o.b.Hi.species.fig, "sites", pch=16, col=colors()[450], select=o.b.Hi.meta.species$Order == "Fabales")
points(o.b.Hi.species.fig, "sites", pch=16, col=colors()[475], select=o.b.Hi.meta.species$Order == "Fagales")
#points(o.b.Hi.species.fig, "sites", pch=16, col=colors()[500], select=o.b.Hi.meta.species$Order == "Lamiales")
#points(o.b.Hi.species.fig, "sites", pch=16, col=colors()[525], select=o.b.Hi.meta.species$Order == "Magnoliales")
points(o.b.Hi.species.fig, "sites", pch=16, col=colors()[550], select=o.b.Hi.meta.species$Order == "Pinales")
#points(o.b.Hi.species.fig, "sites", pch=16, col=colors()[575], select=o.b.Hi.meta.species$Order == "Proteales")
points(o.b.Hi.species.fig, "sites", pch=16, col=colors()[600], select=o.b.Hi.meta.species$Order == "Rosales")
#points(o.b.Hi.species.fig, "sites", pch=16, col=colors()[625], select=o.b.Hi.meta.species$Order == "Sapindales")
#points(o.b.Hi.species.fig, "sites", pch=16, col=colors()[650], select=o.b.Hi.meta.species$Order == "Vitales")
ordiellipse(o.b.Hi.species.mds, o.b.Hi.meta.species$Order, conf = 0.95, label = TRUE, show.groups="Rosales", col=colors()[600])
ordiellipse(o.b.Hi.species.mds, o.b.Hi.meta.species$Order, conf = 0.95, label = TRUE, show.groups="Fagales", col=colors()[475])
ordiellipse(o.b.Hi.species.mds, o.b.Hi.meta.species$Order, conf = 0.95, label = TRUE, show.groups="Pinales", col=colors()[550])
par(op)
op <- par(usr=c(1,0,1,0), xpd=NA)
legend(0.97,1.1, pch= 16, col= c(colors()[653], colors()[475], colors()[550], colors()[600]),  legend=c("Dipsacales", "Fagales", "Pinales", "Rosales"),xpd = TRUE, horiz = TRUE, inset = c(0, 0), bty = "n")
dev.off()

op <- par(mfrow=c(2,2), oma=c(3,0,0,0))
t.t.Lo.species.fig <- ordiplot(t.t.Lo.species.mds, type="none", main="Year 3: Low Top", bty = "n")
points(t.t.Lo.species.fig, "sites", pch=19, col=colors()[375], select=t.t.Lo.meta.species$Order == "Cornales")
#points(t.t.Lo.species.fig, "sites", pch=19, col=colors()[653], select=t.t.Lo.meta.species$Order == "Dipsacales")
points(t.t.Lo.species.fig, "sites", pch=19, col=colors()[425], select=t.t.Lo.meta.species$Order == "Ericales")
points(t.t.Lo.species.fig, "sites", pch=19, col=colors()[450], select=t.t.Lo.meta.species$Order == "Fabales")
points(t.t.Lo.species.fig, "sites", pch=19, col=colors()[475], select=t.t.Lo.meta.species$Order == "Fagales")
points(t.t.Lo.species.fig, "sites", pch=19, col=colors()[500], select=t.t.Lo.meta.species$Order == "Lamiales")
points(t.t.Lo.species.fig, "sites", pch=19, col=colors()[525], select=t.t.Lo.meta.species$Order == "Magnoliales")
points(t.t.Lo.species.fig, "sites", pch=19, col=colors()[550], select=t.t.Lo.meta.species$Order == "Pinales")
points(t.t.Lo.species.fig, "sites", pch=19, col=colors()[575], select=t.t.Lo.meta.species$Order == "Proteales")
points(t.t.Lo.species.fig, "sites", pch=19, col=colors()[600], select=t.t.Lo.meta.species$Order == "Rosales")
points(t.t.Lo.species.fig, "sites", pch=19, col=colors()[625], select=t.t.Lo.meta.species$Order == "Sapindales")
points(t.t.Lo.species.fig, "sites", pch=19, col=colors()[650], select=t.t.Lo.meta.species$Order == "Vitales")
ordiellipse(t.t.Lo.species.mds, t.t.Lo.meta.species$Order, conf = 0.95, label = TRUE, show.groups="Rosales", col=colors()[600])
ordiellipse(t.t.Lo.species.mds, t.t.Lo.meta.species$Order, conf = 0.95, label = TRUE, show.groups="Fagales", col=colors()[475])
ordiellipse(t.t.Lo.species.mds, t.t.Lo.meta.species$Order, conf = 0.95, label = TRUE, show.groups="Pinales", col=colors()[550])
t.t.Hi.species.fig <- ordiplot(t.t.Hi.species.mds, type="none", main="Year 3: High Top", bty = "n")
points(t.t.Hi.species.fig, "sites", pch=19, col=colors()[375], select=t.t.Hi.meta.species$Order == "Cornales")
#points(t.t.Hi.species.fig, "sites", pch=19, col=colors()[653], select=t.t.Hi.meta.species$Order == "Dipsacales")
points(t.t.Hi.species.fig, "sites", pch=19, col=colors()[425], select=t.t.Hi.meta.species$Order == "Ericales")
points(t.t.Hi.species.fig, "sites", pch=19, col=colors()[450], select=t.t.Hi.meta.species$Order == "Fabales")
points(t.t.Hi.species.fig, "sites", pch=19, col=colors()[475], select=t.t.Hi.meta.species$Order == "Fagales")
points(t.t.Hi.species.fig, "sites", pch=19, col=colors()[500], select=t.t.Hi.meta.species$Order == "Lamiales")
points(t.t.Hi.species.fig, "sites", pch=19, col=colors()[525], select=t.t.Hi.meta.species$Order == "Magnoliales")
points(t.t.Hi.species.fig, "sites", pch=19, col=colors()[550], select=t.t.Hi.meta.species$Order == "Pinales")
points(t.t.Hi.species.fig, "sites", pch=19, col=colors()[575], select=t.t.Hi.meta.species$Order == "Proteales")
points(t.t.Hi.species.fig, "sites", pch=19, col=colors()[600], select=t.t.Hi.meta.species$Order == "Rosales")
points(t.t.Hi.species.fig, "sites", pch=19, col=colors()[625], select=t.t.Hi.meta.species$Order == "Sapindales")
points(t.t.Hi.species.fig, "sites", pch=19, col=colors()[650], select=t.t.Hi.meta.species$Order == "Vitales")
ordiellipse(t.t.Hi.species.mds, t.t.Hi.meta.species$Order, conf = 0.95, label = TRUE, show.groups="Rosales", col=colors()[600])
ordiellipse(t.t.Hi.species.mds, t.t.Hi.meta.species$Order, conf = 0.95, label = TRUE, show.groups="Fagales", col=colors()[475])
ordiellipse(t.t.Hi.species.mds, t.t.Hi.meta.species$Order, conf = 0.95, label = TRUE, show.groups="Pinales", col=colors()[550])
t.b.Lo.species.fig <- ordiplot(t.b.Lo.species.mds, type="none", main="Year 3: Low Bottom", bty = "n")
points(t.b.Lo.species.fig, "sites", pch=19, col=colors()[375], select=t.b.Lo.meta.species$Order == "Cornales")
#points(t.b.Lo.species.fig, "sites", pch=19, col=colors()[653], select=t.b.Lo.meta.species$Order == "Dipsacales")
points(t.b.Lo.species.fig, "sites", pch=19, col=colors()[425], select=t.b.Lo.meta.species$Order == "Ericales")
points(t.b.Lo.species.fig, "sites", pch=19, col=colors()[450], select=t.b.Lo.meta.species$Order == "Fabales")
points(t.b.Lo.species.fig, "sites", pch=19, col=colors()[475], select=t.b.Lo.meta.species$Order == "Fagales")
points(t.b.Lo.species.fig, "sites", pch=19, col=colors()[500], select=t.b.Lo.meta.species$Order == "Lamiales")
points(t.b.Lo.species.fig, "sites", pch=19, col=colors()[525], select=t.b.Lo.meta.species$Order == "Magnoliales")
points(t.b.Lo.species.fig, "sites", pch=19, col=colors()[550], select=t.b.Lo.meta.species$Order == "Pinales")
points(t.b.Lo.species.fig, "sites", pch=19, col=colors()[575], select=t.b.Lo.meta.species$Order == "Proteales")
points(t.b.Lo.species.fig, "sites", pch=19, col=colors()[600], select=t.b.Lo.meta.species$Order == "Rosales")
points(t.b.Lo.species.fig, "sites", pch=19, col=colors()[625], select=t.b.Lo.meta.species$Order == "Sapindales")
points(t.b.Lo.species.fig, "sites", pch=19, col=colors()[650], select=t.b.Lo.meta.species$Order == "Vitales")
ordiellipse(t.b.Lo.species.mds, t.b.Lo.meta.species$Order, conf = 0.95, label = TRUE, show.groups="Rosales", col=colors()[600])
ordiellipse(t.b.Lo.species.mds, t.b.Lo.meta.species$Order, conf = 0.95, label = TRUE, show.groups="Fagales", col=colors()[475])
ordiellipse(t.b.Lo.species.mds, t.b.Lo.meta.species$Order, conf = 0.95, label = TRUE, show.groups="Pinales", col=colors()[550])
t.b.Hi.species.fig <- ordiplot(t.b.Hi.species.mds, type="none", main="Year 3: High Bottom", bty = "n")
points(t.b.Hi.species.fig, "sites", pch=19, col=colors()[375], select=t.b.Hi.meta.species$Order == "Cornales")
#points(t.b.Hi.species.fig, "sites", pch=19, col=colors()[653], select=t.b.Hi.meta.species$Order == "Dipsacales")
points(t.b.Hi.species.fig, "sites", pch=19, col=colors()[425], select=t.b.Hi.meta.species$Order == "Ericales")
points(t.b.Hi.species.fig, "sites", pch=19, col=colors()[450], select=t.b.Hi.meta.species$Order == "Fabales")
points(t.b.Hi.species.fig, "sites", pch=19, col=colors()[475], select=t.b.Hi.meta.species$Order == "Fagales")
points(t.b.Hi.species.fig, "sites", pch=19, col=colors()[500], select=t.b.Hi.meta.species$Order == "Lamiales")
points(t.b.Hi.species.fig, "sites", pch=19, col=colors()[525], select=t.b.Hi.meta.species$Order == "Magnoliales")
points(t.b.Hi.species.fig, "sites", pch=19, col=colors()[550], select=t.b.Hi.meta.species$Order == "Pinales")
points(t.b.Hi.species.fig, "sites", pch=19, col=colors()[575], select=t.b.Hi.meta.species$Order == "Proteales")
points(t.b.Hi.species.fig, "sites", pch=19, col=colors()[600], select=t.b.Hi.meta.species$Order == "Rosales")
points(t.b.Hi.species.fig, "sites", pch=19, col=colors()[625], select=t.b.Hi.meta.species$Order == "Sapindales")
points(t.b.Hi.species.fig, "sites", pch=19, col=colors()[650], select=t.b.Hi.meta.species$Order == "Vitales")
ordiellipse(t.b.Hi.species.mds, t.b.Hi.meta.species$Order, conf = 0.95, label = TRUE, show.groups="Rosales", col=colors()[600])
ordiellipse(t.b.Hi.species.mds, t.b.Hi.meta.species$Order, conf = 0.95, label = TRUE, show.groups="Fagales", col=colors()[475])
ordiellipse(t.b.Hi.species.mds, t.b.Hi.meta.species$Order, conf = 0.95, label = TRUE, show.groups="Pinales", col=colors()[550])
par(op)
op <- par(usr=c(1,0,1,0), xpd=NA)
legend(0.97,1.1, pch= 19, col= c(colors()[375], colors()[425], colors()[450], colors()[475], colors()[500], colors()[525], colors()[550], colors()[575], colors()[600], colors()[625], colors()[650]), legend=c("Cornales", "Ericales", "Fabales", "Fagales", "Lamiales", "Magnoliales", "Pinales", "Proteales", "Rosales", "Sapindales", "Vitales"), xpd = TRUE, horiz = TRUE, inset = c(0, 0), bty = "n")

plot(envfit(o.t.Lo.species.mds, o.t.Lo.meta.species[,13:ncol(o.t.Lo.meta.species)], na.rm=T))

RDA <- rda(o.t.Lo.species ~ Conduit.D.um. + Conduit.L.m. + avg.top.bottom.density + WCProtPerDM + WADFPerDM + WNDFPerDM + WLigPerDM + WPPerDM + Wd.Per_N + Wd.Per_C + Cellulose.cell + Hemicellulose.hemic, o.t.Lo.meta.species, na.action=na.omit)   
anova(RDA)
summary(RDA)
plot(RDA)

CCA1 <- cca(o.t.Lo.species ~ Conduit.D.um. + Conduit.L.m. + avg.top.bottom.density + WCProtPerDM + WADFPerDM + WNDFPerDM + WLigPerDM + WPPerDM + Wd.Per_N + Wd.Per_C + Cellulose.cell + Hemicellulose.hemic, o.t.Lo.meta.species)   
CCA2 <- cca(o.t.Hi.species ~ Conduit.D.um. + Conduit.L.m. + avg.top.bottom.density + WCProtPerDM + WADFPerDM + WNDFPerDM + WLigPerDM + WPPerDM + Wd.Per_N + Wd.Per_C + Cellulose.cell + Hemicellulose.hemic, o.t.Hi.meta.species)   
CCA3 <- cca(o.b.Lo.species ~ Conduit.D.um. + Conduit.L.m. + avg.top.bottom.density + WCProtPerDM + WADFPerDM + WNDFPerDM + WLigPerDM + WPPerDM + Wd.Per_N + Wd.Per_C + Cellulose.cell + Hemicellulose.hemic, o.b.Lo.meta.species)   
CCA4 <- cca(o.b.Hi.species ~ Conduit.D.um. + Conduit.L.m. + avg.top.bottom.density + WCProtPerDM + WADFPerDM + WNDFPerDM + WLigPerDM + WPPerDM + Wd.Per_N + Wd.Per_C + Cellulose.cell + Hemicellulose.hemic, o.b.Hi.meta.species)   

par(mfrow=c(2,2))
plot(CCA1, display=c("wa", "bp"), type="text", choices=c(1,2), main="Year 3: Low top")
plot(CCA2, display=c("wa", "bp"), type="text", choices=c(1,2), main="Year 3: High top")
plot(CCA3, display=c("wa", "bp"), type="text", choices=c(1,2), main="Year 3: Low bottom")
plot(CCA4, display=c("wa", "bp"), type="text", choices=c(1,2), main="Year 3: High bottom")


t.t.Lo.CCA <- cca(t.t.Lo.species ~ Conduit.D.um. + Conduit.L.m. + avg.top.bottom.density + WCProtPerDM + WADFPerDM + WNDFPerDM + WLigPerDM + WPPerDM + Wd.Per_N + Wd.Per_C + Cellulose.cell + Hemicellulose.hemic, t.t.Lo.meta.species)   
t.t.Hi.CCA <- cca(t.t.Hi.species ~ Conduit.D.um. + Conduit.L.m. + avg.top.bottom.density + WCProtPerDM + WADFPerDM + WNDFPerDM + WLigPerDM + WPPerDM + Wd.Per_N + Wd.Per_C + Cellulose.cell + Hemicellulose.hemic, t.t.Hi.meta.species)   
t.b.Lo.CCA <- cca(t.b.Lo.species ~ Conduit.D.um. + Conduit.L.m. + avg.top.bottom.density + WCProtPerDM + WADFPerDM + WNDFPerDM + WLigPerDM + WPPerDM + Wd.Per_N + Wd.Per_C + Cellulose.cell + Hemicellulose.hemic, t.b.Lo.meta.species)   
t.b.Hi.CCA <- cca(t.b.Hi.species ~ Conduit.D.um. + Conduit.L.m. + avg.top.bottom.density + WCProtPerDM + WADFPerDM + WNDFPerDM + WLigPerDM + WPPerDM + Wd.Per_N + Wd.Per_C + Cellulose.cell + Hemicellulose.hemic, t.b.Hi.meta.species)   
o.t.Lo.CCA <- cca(o.t.Lo.species ~ Conduit.D.um. + Conduit.L.m. + avg.top.bottom.density + WCProtPerDM + WADFPerDM + WNDFPerDM + WLigPerDM + WPPerDM + Wd.Per_N + Wd.Per_C + Cellulose.cell + Hemicellulose.hemic, o.t.Lo.meta.species)   
o.t.Hi.CCA <- cca(o.t.Hi.species ~ Conduit.D.um. + Conduit.L.m. + avg.top.bottom.density + WCProtPerDM + WADFPerDM + WNDFPerDM + WLigPerDM + WPPerDM + Wd.Per_N + Wd.Per_C + Cellulose.cell + Hemicellulose.hemic, o.t.Hi.meta.species)   
o.b.Lo.CCA <- cca(o.b.Lo.species ~ Conduit.D.um. + Conduit.L.m. + avg.top.bottom.density + WCProtPerDM + WADFPerDM + WNDFPerDM + WLigPerDM + WPPerDM + Wd.Per_N + Wd.Per_C + Cellulose.cell + Hemicellulose.hemic, o.b.Lo.meta.species)   
o.b.Hi.CCA <- cca(o.b.Hi.species ~ Conduit.D.um. + Conduit.L.m. + avg.top.bottom.density + WCProtPerDM + WADFPerDM + WNDFPerDM + WLigPerDM + WPPerDM + Wd.Per_N + Wd.Per_C + Cellulose.cell + Hemicellulose.hemic, o.b.Hi.meta.species)   

anova(t.t.Lo.CCA, by="term")
anova(t.t.Hi.CCA, by="term")
anova(t.b.Lo.CCA, by="term")
anova(t.b.Hi.CCA, by="term")
anova(o.t.Lo.CCA, by="term")
anova(o.t.Hi.CCA, by="term")
anova(o.b.Lo.CCA, by="term")
anova(o.b.Hi.CCA, by="term")

###***************************####
o.species <- subset(species.tyson,rownames(species.tyson) %in% rownames(meta.species[meta.species$HarvestYear == 1,]))
t.species <- subset(species.tyson,rownames(species.tyson) %in% rownames(meta.species[meta.species$HarvestYear == 3,]))
o.meta.species <- meta.species[row.names(o.species),]
t.meta.species <- meta.species[row.names(t.species),]
t.cca <- cca(t.species ~ Conduit.D.um. * Conduit.L.m. + avg.top.bottom.density * WCProtPerDM * WADFPerDM * WNDFPerDM + WLigPerDM * WPPerDM * Wd.Per_N * Wd.Per_C * Cellulose.cell * Hemicellulose.hemic + PlotLocation * Pos * PotID, t.meta.species)   
o.cca <- cca(o.species ~ Conduit.D.um. * Conduit.L.m. + avg.top.bottom.density * WCProtPerDM * WADFPerDM * WNDFPerDM + WLigPerDM * WPPerDM * Wd.Per_N * Wd.Per_C * Cellulose.cell * Hemicellulose.hemic + PlotLocation * Pos * PotID, o.meta.species)   
o.cca.2 <- cca(o.species ~ Hemicellulose.hemic + Conduit.D.um. + Conduit.L.m. + avg.top.bottom.density + WCProtPerDM + WADFPerDM + WNDFPerDM + WLigPerDM + WPPerDM + Wd.Per_N + Wd.Per_C + Cellulose.cell, o.meta.species)   
o.cca.2 <- cca(o.species ~ Conduit.D.um. * Conduit.L.m. + avg.top.bottom.density * WCProtPerDM * WADFPerDM * WNDFPerDM + WLigPerDM * WPPerDM * Wd.Per_N * Wd.Per_C * Cellulose.cell * Hemicellulose.hemic, o.meta.species)   
o.species.stand <- decostand(o.species, method="hellinger")
o.species.mds <- metaMDS(o.species.stand, dist = "euclidean",4)
t.species.stand <- decostand(t.species, method="hellinger")
t.species.mds <- metaMDS(t.species.stand, dist = "euclidean",k=4, trymax=999, startmax=0.999999)
ordiplot(t.bacteria.mds, type="p", display="sites")
plot(t.ef.m.bacteria)
plot(t.bacteria.mds, display = "sites", type = "p")
with(t.meta.bacteria, ordiellipse(t.bacteria.mds, Clade, conf = 0.999, label = TRUE))

##to do
#envfit with new models
#try with preseance absence and counts
#stacked bar plots
#figures for ESA - AMY
#current vs. initial trait data
#adonis

t.ef <- envfit(t.species.mds, t.meta.species[,c(1,2,7,10,13:ncol(t.meta.species))], permu = 999)
o.ef <- envfit(o.species.mds, o.meta.species[,c(1,2,7,10,13:ncol(t.meta.species))], permu = 999)
t.ef.m <- envfit(t.species.mds ~ Clade + Conduit.D.um. * Conduit.L.m. + avg.top.bottom.density * WCProtPerDM * WADFPerDM * WNDFPerDM + WLigPerDM * WPPerDM * Wd.Per_N * Wd.Per_C * Cellulose.cell * Hemicellulose.hemic + PlotLocation * Pos * factor(PotID), t.meta.species)   
o.ef.m <- envfit(o.species.mds ~ Clade + Conduit.D.um. * Conduit.L.m. + avg.top.bottom.density * WCProtPerDM * WADFPerDM * WNDFPerDM + WLigPerDM * WPPerDM * Wd.Per_N * Wd.Per_C * Cellulose.cell * Hemicellulose.hemic + PlotLocation * Pos * factor(PotID), o.meta.species)   
o.bacteria <- subset(bacteria.tyson,rownames(bacteria.tyson) %in% rownames(meta.bacteria[meta.bacteria$HarvestYear == 1,]))
t.bacteria <- subset(bacteria.tyson,rownames(bacteria.tyson) %in% rownames(meta.bacteria[meta.bacteria$HarvestYear == 3,]))
o.meta.bacteria <- meta.bacteria[row.names(o.bacteria),]
t.meta.bacteria <- meta.bacteria[row.names(t.bacteria),]
o.bacteria.stand <- decostand(o.bacteria, method="hellinger")
o.bacteria.mds <- metaMDS(o.bacteria.stand, dist = "euclidean",4)
t.bacteria.stand <- decostand(t.bacteria, method="hellinger")
t.bacteria.mds <- metaMDS(t.bacteria.stand, dist = "euclidean",k=4, trymax=999, startmax=0.999999)
ordiplot(t.bacteria.mds, type="t", display="sites")
t.ef.bacteria <- envfit(t.bacteria.mds, t.meta.bacteria[,c(1,2,7,10,13:ncol(t.meta.bacteria))], permu = 999)
o.ef.bacteria <- envfit(o.bacteria.mds, o.meta.bacteria[,c(1,2,7,10,13:ncol(t.meta.bacteria))], permu = 999)
t.ef.m.bacteria <- envfit(t.bacteria.mds ~ Clade + Conduit.D.um. * Conduit.L.m. + avg.top.bottom.density * WCProtPerDM * WADFPerDM * WNDFPerDM + WLigPerDM * WPPerDM * Wd.Per_N * Wd.Per_C * Cellulose.cell * Hemicellulose.hemic + PlotLocation * Pos * factor(PotID), t.meta.bacteria)   
o.ef.m.bacteria <- envfit(o.bacteria.mds ~ Clade + Conduit.D.um. * Conduit.L.m. + avg.top.bottom.density * WCProtPerDM * WADFPerDM * WNDFPerDM + WLigPerDM * WPPerDM * Wd.Per_N * Wd.Per_C * Cellulose.cell * Hemicellulose.hemic + PlotLocation * Pos * factor(PotID), o.meta.bacteria)   

o.species.bac <- o.species[row.names(o.bacteria) %in% row.names(o.species),]
o.sp.bac.stand <- decostand(o.species.bac, method="hellinger")
o.sp.bac.mds <- metaMDS(o.sp.bac.stand, dist = "euclidean",4)
o.proc <- procrustes(o.sp.bac.mds, o.bacteria.mds)
plot(o.proc)
protest(o.sp.bac.mds, o.bacteria.mds)

plot(o.species.mds, display = "sites", type = "p")
with(o.meta.species, ordiellipse(o.species.mds, Clade, conf = 0.95))



#Based on BEST model for year 1 and year 3
#Best model for year 1 after permutation test on cca: 
#Model: cca(formula = o.species ~ Conduit.D.um. + WPPerDM + WCProtPerDM + WADFPerDM, data = o.meta.species[, c(1, 2, 7, 13:ncol(o.meta.species))])
#Best Model for year 3:
#Model: cca(formula = t.species ~ Wd.Per_C + PotID + PlotLocation + WCProtPerDM + avg.top.bottom.density + WLigPerDM + WPPerDM, data = t.meta.species[, c(1, 2, 7, 13:ncol(t.meta.species))])
t.cca.best <- cca(t.species ~ Conduit.D.um. + avg.top.bottom.density * WCProtPerDM * WADFPerDM + WLigPerDM * WPPerDM + PlotLocation * PotID, t.meta.species)   
o.cca.best <- cca(o.species ~ Conduit.D.um. + avg.top.bottom.density * WCProtPerDM * WADFPerDM * WLigPerDM * WPPerDM + PlotLocation * PotID, o.meta.species)   




t.cca.A <- cca(t.species ~ Conduit.D.um. * Conduit.L.m., t.meta.species)   
o.cca.A <- cca(o.species ~ Conduit.D.um. * Conduit.L.m., o.meta.species)   
t.cca.F <- cca(t.species ~ avg.top.bottom.density * WCProtPerDM * WADFPerDM * WNDFPerDM, t.meta.species)   
o.cca.F <- cca(o.species ~ avg.top.bottom.density * WCProtPerDM * WADFPerDM * WNDFPerDM, o.meta.species)   
t.cca.C <- cca(t.species ~ WLigPerDM * WPPerDM * Wd.Per_N * Wd.Per_C * Cellulose.cell * Hemicellulose.hemic, t.meta.species)   
o.cca.C <- cca(o.species ~ WLigPerDM * WPPerDM * Wd.Per_N * Wd.Per_C * Cellulose.cell * Hemicellulose.hemic, o.meta.species)   
t.cca.E <- cca(t.species ~ PlotLocation * Pos * PotID, t.meta.species)   
o.cca.E <- cca(o.species ~ PlotLocation * Pos * PotID, o.meta.species)   

anova(t.cca.A, by="term")
anova(o.cca.A, by="term")
anova(t.cca.F, by="term")
anova(o.cca.F, by="term")
anova(t.cca.C, by="term")
anova(o.cca.C, by="term")
anova(t.cca.E, by="term")
anova(o.cca.E, by="term")

#year 1
o.mod1 <- cca(o.species ~ ., o.meta.species[,c(1,2,7,13:ncol(o.meta.species))])
o.mod0 <- cca(o.species ~ 1, o.meta.species[,c(1,2,7,13:ncol(o.meta.species))])
o.mod <- step(o.mod0, scope = formula(o.mod1), test = "perm")

#Year 3
t.mod1 <- cca(t.species ~ ., t.meta.species[,c(1,2,7,13:ncol(t.meta.species))])
t.mod0 <- cca(t.species ~ 1, t.meta.species[,c(1,2,7,13:ncol(t.meta.species))])
t.mod <- step(t.mod0, scope = formula(t.mod1), test = "perm")

#Year 1 Hi
#Year 1 Lo
o.Lo.species <- subset(Lo.species,rownames(Lo.species) %in% rownames(meta[meta$PlotLocation == "L" & meta$HarvestYear == 1,]))
o.Hi.species <- subset(Hi.species,rownames(Hi.species) %in% rownames(meta[meta$PlotLocation == "H" & meta$HarvestYear == 1,]))
o.Lo.meta.species <- meta.species.rep[row.names(o.Lo.species),]
o.Hi.meta.species <- meta.species.rep[row.names(o.Hi.species),]
o.Lo.mod1 <- cca(o.Lo.species ~ ., o.Lo.meta.species[,c(1,2,7,13:ncol(o.Lo.meta.species))])
o.Lo.mod0 <- cca(o.Lo.species ~ 1, o.Lo.meta.species[,c(1,2,7,13:ncol(o.Lo.meta.species))])
o.Lo.mod <- step(o.Lo.mod0, scope = formula(o.Lo.mod1), test = "perm")
o.Hi.mod1 <- cca(o.Hi.species ~ ., o.Hi.meta.species[,c(1,2,7,13:ncol(o.Hi.meta.species))])
o.Hi.mod0 <- cca(o.Hi.species ~ 1, o.Hi.meta.species[,c(1,2,7,13:ncol(o.Hi.meta.species))])
o.Hi.mod <- step(o.Hi.mod0, scope = formula(o.Hi.mod1), test = "perm")

#Year 3 Hi
#Year 3 Lo
t.Lo.species <- subset(Lo.species,rownames(Lo.species) %in% rownames(meta[meta$PlotLocation == "L" & meta$HarvestYear == 3,]))
t.Hi.species <- subset(Hi.species,rownames(Hi.species) %in% rownames(meta[meta$PlotLocation == "H" & meta$HarvestYear == 3,]))
t.Lo.meta.species <- meta.species.rep[row.names(t.Lo.species),]
t.Hi.meta.species <- meta.species.rep[row.names(t.Hi.species),]
t.Lo.mod1 <- cca(t.Lo.species ~ ., t.Lo.meta.species[,c(1,2,7,13:ncol(t.Lo.meta.species))])
t.Lo.mod0 <- cca(t.Lo.species ~ 1, t.Lo.meta.species[,c(1,2,7,13:ncol(t.Lo.meta.species))])
t.Lo.mod <- step(t.Lo.mod0, scope = formula(o.Lo.mod1), test = "perm")
t.Hi.mod1 <- cca(t.Hi.species ~ ., t.Hi.meta.species[,c(1,2,7,13:ncol(t.Hi.meta.species))])
t.Hi.mod0 <- cca(t.Hi.species ~ 1, t.Hi.meta.species[,c(1,2,7,13:ncol(t.Hi.meta.species))])
t.Hi.mod <- step(t.Hi.mod0, scope = formula(t.Hi.mod1), test = "perm")

anova(t.cca, by="term")
anova(o.cca, by="term")



t.t.Lo.rda <- rda(t.t.Lo.species ~ Conduit.D.um. + Conduit.L.m. + avg.top.bottom.density + WCProtPerDM + WADFPerDM + WNDFPerDM + WLigPerDM + WPPerDM + Wd.Per_N + Wd.Per_C + Cellulose.cell + Hemicellulose.hemic, t.t.Lo.meta.species, na.action=na.omit)   
t.t.Hi.rda <- rda(t.t.Hi.species ~ Conduit.D.um. + Conduit.L.m. + avg.top.bottom.density + WCProtPerDM + WADFPerDM + WNDFPerDM + WLigPerDM + WPPerDM + Wd.Per_N + Wd.Per_C + Cellulose.cell + Hemicellulose.hemic, t.t.Hi.meta.species, na.action=na.omit)   
t.b.Lo.rda <- rda(t.b.Lo.species ~ Conduit.D.um. + Conduit.L.m. + avg.top.bottom.density + WCProtPerDM + WADFPerDM + WNDFPerDM + WLigPerDM + WPPerDM + Wd.Per_N + Wd.Per_C + Cellulose.cell + Hemicellulose.hemic, t.b.Lo.meta.species, na.action=na.omit)   
t.b.Hi.rda <- rda(t.b.Hi.species ~ Conduit.D.um. + Conduit.L.m. + avg.top.bottom.density + WCProtPerDM + WADFPerDM + WNDFPerDM + WLigPerDM + WPPerDM + Wd.Per_N + Wd.Per_C + Cellulose.cell + Hemicellulose.hemic, t.b.Hi.meta.species, na.action=na.omit)   
o.t.Lo.rda <- rda(o.t.Lo.species ~ Conduit.D.um. + Conduit.L.m. + avg.top.bottom.density + WCProtPerDM + WADFPerDM + WNDFPerDM + WLigPerDM + WPPerDM + Wd.Per_N + Wd.Per_C + Cellulose.cell + Hemicellulose.hemic, o.t.Lo.meta.species, na.action=na.omit)   
o.t.Hi.rda <- rda(o.t.Hi.species ~ Conduit.D.um. + Conduit.L.m. + avg.top.bottom.density + WCProtPerDM + WADFPerDM + WNDFPerDM + WLigPerDM + WPPerDM + Wd.Per_N + Wd.Per_C + Cellulose.cell + Hemicellulose.hemic, o.t.Hi.meta.species, na.action=na.omit)   
o.b.Lo.rda <- rda(o.b.Lo.species ~ Conduit.D.um. + Conduit.L.m. + avg.top.bottom.density + WCProtPerDM + WADFPerDM + WNDFPerDM + WLigPerDM + WPPerDM + Wd.Per_N + Wd.Per_C + Cellulose.cell + Hemicellulose.hemic, o.b.Lo.meta.species, na.action=na.omit)   
o.b.Hi.rda <- rda(o.b.Hi.species ~ Conduit.D.um. + Conduit.L.m. + avg.top.bottom.density + WCProtPerDM + WADFPerDM + WNDFPerDM + WLigPerDM + WPPerDM + Wd.Per_N + Wd.Per_C + Cellulose.cell + Hemicellulose.hemic, o.b.Hi.meta.species, na.action=na.omit)   

anova(t.t.Lo.rda, by="term")
anova(t.t.Hi.rda, by="term")
anova(t.b.Lo.rda, by="term")
anova(t.b.Hi.rda, by="term")
anova(o.t.Lo.rda, by="term")
anova(o.t.Hi.rda, by="term")
anova(o.b.Lo.rda, by="term")
anova(o.b.Hi.rda, by="term")




t.t.Lo.CCA.s <- summary(t.t.Lo.CCA)  
t.t.Hi.CCA.s <- summary(t.t.Hi.CCA)
t.b.Lo.CCA.s <- summary(t.b.Lo.CCA) 
t.b.Hi.CCA.s <- summary(t.b.Hi.CCA)  
o.t.Lo.CCA.s <- summary(o.t.Lo.CCA) 
o.t.Hi.CCA.s <- summary(o.t.Hi.CCA) 
o.b.Lo.CCA.s <- summary(o.b.Lo.CCA) 
o.b.Hi.CCA.s <- summary(o.b.Hi.CCA)

anova(CCA1, by="term")


RDA <- rda(species.tyson ~ Conduit.D.um. + Conduit.L.m. + avg.top.bottom.density + WCProtPerDM + WADFPerDM + WNDFPerDM + WLigPerDM + WPPerDM + Wd.Per_N + Wd.Per_C + Cellulose.cell + Hemicellulose.hemic, meta.species, na.action=na.omit)   
CCA <- cca(species.tyson ~ Conduit.D.um. * Conduit.L.m. + avg.top.bottom.density * WCProtPerDM * WADFPerDM * WNDFPerDM + WLigPerDM * WPPerDM * Wd.Per_N * Wd.Per_C * Cellulose.cell * Hemicellulose.hemic, meta.species, na.action=na.omit)   
plot(CCA, display=c("wa", "bp"), type="text", choices=c(1,2))




