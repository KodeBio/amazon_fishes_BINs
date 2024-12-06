# Load utility functions
source("utils.R")

#### ----------------------------
#### 1. Load and Clean Data
#### ----------------------------

amz_fish_bin <- read.csv(file = "data/amz_BOLD_bins.csv")
amz_seqs <- read.fasta(dir = "data/amz_BOLD_alig.fasta")

#removing sequences that were manually filtered in BioEdit
amz_fish_bin <- amz_fish_bin[which(amz_fish_bin$sequenceID%in%names(amz_seqs)),]

#CChanging “~” to “-”
amz_seqs <- lapply(amz_seqs, function(i) gsub("~", "-", i))

#basins
basins <- c("Upper Amazon", "Negro", "Madeira", "Trombetas", "Paru/Jari", "Tapajos", "Xingu", "Tocantins")
basins_initials <- c("UA", "Ne", "Ma", "Tr", "PJ", "Ta", "Xi", "To")
amz_fish_bin[amz_fish_bin$basin=="Alto Amazonas", "basin"] <- "Upper Amazon"

amz_colors <- c("#badca7", "#859e94", "#d3e6c9", "#bcc8c2",
                "#cae3bc", "#88aa9b", "#c2d08b", "#badca7")


#### ----------------------------
#### 2. Basin Data Evaluation
#### ----------------------------

# Graph 1: Source of information (coords or articles) ====
#sequences found manually and by coordinates
coord_seq <- sum(!is.na(amz_fish_bin$lat))
manual_seq <- sum(is.na(amz_fish_bin$lat))

#Bins encontrados manualmente y por coordenadas
pre_BIN <- unique(amz_fish_bin$bin_uri[!is.na(amz_fish_bin$lat)])
post_BIN <- unique(amz_fish_bin$bin_uri[is.na(amz_fish_bin$lat)])

coord_BINs <- length(setdiff(pre_BIN, post_BIN))
mixed_BINs <- length(intersect(pre_BIN, post_BIN))
manual_BINs <- length(setdiff(post_BIN, pre_BIN))

#of mixed BINs, how many seqs are by coordinates and how many are manual
coord_mixed_seq <- length(amz_fish_bin$sequenceID[amz_fish_bin$bin_uri%in%intersect(pre_BIN, post_BIN) 
                                                  & !is.na(amz_fish_bin$lat)])
manual_mixed_seq <- length(amz_fish_bin$sequenceID[amz_fish_bin$bin_uri%in%intersect(pre_BIN, post_BIN) 
                                                   & is.na(amz_fish_bin$lat)])

#obtaining the proportions
S1 <- (coord_seq-coord_mixed_seq)/(coord_seq+manual_seq)
S2 <- (coord_mixed_seq)/(coord_seq+manual_seq)
S3 <- (manual_mixed_seq)/(coord_seq+manual_seq)
S4 <- (manual_seq-manual_mixed_seq)/(coord_seq+manual_seq)

B1 <- coord_BINs/(coord_BINs+mixed_BINs+manual_BINs)
B2 <- mixed_BINs/(coord_BINs+mixed_BINs+manual_BINs)
B3 <- mixed_BINs/(coord_BINs+mixed_BINs+manual_BINs)
B4 <- manual_BINs/(coord_BINs+mixed_BINs+manual_BINs)

##graphing
plot(0, 0, type = "n", xlab = "Number of sequences", ylab = "Number of BINs", 
     xlim = c(0, 1), ylim = c(0, 0.65), 
     xaxs="i", yaxs="i", xaxs = "i", axes = FALSE)

G1TS <- 0.8

#altura bin, seq wide
rect(S1, 0, S1+S2, B2, col="#859e94", border ="#859e94")
rect(S1+S2, 0, S1+S2+S3, B3, col="#d3e6c9",  border ="#d3e6c9")
rect(0, 0, S1, B1, col="#859e94")
rect(S1+S2+S3, 0, 1, B4, col="#d3e6c9")
segments(S1, B2, S1+S2+S3, B3, col="black")

#colocando los valores en el gráfico
text(S1/2, 0.03, labels = paste0("(", coord_seq-coord_mixed_seq, ")"), cex = G1TS)
text(S1+S2/2, 0.03, labels = paste0("(", coord_mixed_seq, ")"), cex = G1TS)
text(S1+S2+S3/2, 0.03, labels = paste0("(", manual_mixed_seq, ")"), cex = G1TS)
text((S1+S2+S3+1)/2, 0.03, labels = paste0("(", manual_seq-manual_mixed_seq, ")"), cex = 0.8)


text(S1/2, B1+0.06, cex = G1TS, 
     labels = paste0("BINs obtained through"))
text(S1/2, B1+0.03, cex = G1TS,
     labels = paste0("coordinates (", coord_BINs, ")"))

text(S1+(S2+S3)/2, B2+0.03, cex = G1TS,
     labels = paste0("Mixed BINs (", mixed_BINs, ")"))

text((S1+S2+S3+1)/2, B4+0.06, cex = G1TS, xpd = T,
     labels = paste0("BINs obtained through"))
text((S1+S2+S3+1)/2, B4+0.03, cex = G1TS, xpd = T,
     labels = paste0("bibliography (", manual_BINs, ")"))


#axis
axis(side=1, at=c(0, S2+S1, 1), labels = FALSE)
axis(side=1, at = c((S2+S1)/2, (S2+S1+1)/2), col.ticks = "white", cex.axis=G1TS,
     labels = c(paste0("By coordinates (", coord_seq, ")"), 
                paste0("By bibliography (", manual_seq, ")")))

axis(side=2, las=2, at=(seq(0, 360, 40)*0.58)/360, labels = seq(0, 360, 40))


# Number of BINs and sequences per basin ====
bin_x_basin <- unique(paste(amz_fish_bin$bin_uri, amz_fish_bin$basin, sep = "_"))
bin_x_basin <- sub(".*_(.*)", "\\1", bin_x_basin)

bin_num <- aggregate(amz_fish_bin$bin_uri, by = list(Basin = amz_fish_bin$basin), 
                     FUN = function(x) c(Sequences = length(x), BINs = length(unique(x))))
bin_num <- do.call(data.frame, bin_num)
colnames(bin_num) <- c("Basin", "Sequences", "BINs")

bin_num$Sequences <- as.numeric(bin_num$Sequences)
bin_num$BINs <- as.numeric(bin_num$BINs)
bin_num <- bin_num[match(basins, bin_num$Basin), ]
bin_num

#graphing
plot(0, 0, type = "n", xlab = "Sub-basin", ylab = "Total", xlim = c(0.5, nrow(bin_num)*1.5), 
     ylim = c(0,1000), axes = FALSE, yaxs = "i")
abline(h = seq(0,1000,100), lty = "dotted", col = "grey70")
bp <- barplot(bin_num$Sequences, names.arg = basins_initials, 
              col = amz_colors, border = FALSE, cex.names = 1,
              ylim = c(0, 1000), yaxs = "i", las = 1, space = 0.5, axes = FALSE, add = TRUE)

lines(bp, bin_num$BINs, type = "o", col = "black", pch = 20)
text(x = bp, y = bin_num$Sequences, labels = bin_num$Sequences, pos = 3, cex = 1)
text(x = bp, y = bin_num$BINs, labels = bin_num$BINs, pos = 3, col = "black", cex = 0.8)
axis(side = 2, at = seq(0,1000,100), labels = seq(0,1000,100), las=1, cex.axis = 0.8)
box()

legend("topright", legend = c("Sequences", "BINs"), 
       fill = c(amz_colors[1], NA), border = c(NA, NA),
       col = c(NA, "black"), 
       pch = c(NA, 19), lty = c(NA, 1),
       bty = "n") # bty = "n" removes the border of the caption box


#### ----------------------------
#### 4. Genetic Distance Within Basins
#### ----------------------------
basin_homogeneity <- function(basin, database = amz_fish_bin, seq_database = amz_seqs, 
                              BINs = NULL, alg="max", dgp="k2p", average=FALSE, minseqs=3){
  minibase <- database[database$basin%in%basin,]
  if(is.null(BINs)) BINs <- unique(minibase$bin_uri)

  distances <- vector()
  for(BIN in BINs){
    seqIDs <- minibase$sequenceID[minibase$bin_uri==BIN]
    if(length(seqIDs) < minseqs) next
    seqs <- unlist(seq_database[as.character(seqIDs)])
    
    dist <-  polidistance(seqs, FUN = dgp)
    distances <- c(distances, match.fun(alg)(dist))
    names(distances)[length(distances)] <- BIN
  }
  distances <- 1/(distances+1)
  
  if(!average) return(distances)
  else return(mean(distances))
}

intraBasin_max <- lapply(basins, function(i){
  basin_homogeneity(i, database=amz_fish_bin, seq_database=amz_seqs, 
                    dgp="k2p", alg="max", minseqs=3)
})
names(intraBasin_max) <- basins

# graphic
intraBasin_max_xd <- intraBasin_max
par(mar=c(6.1, 4.1, 0.1, 2.1))
vioplot::vioplot(intraBasin_max_xd, cex.axis = 0.8, col = amz_colors, border = amz_colors,
                 yaxs="i", las=2, ylim = c(min(unlist(intraBasin_max_xd))-0.05, max(unlist(intraBasin_max_xd))+0.05), axes = F)
stripchart(intraBasin_max_xd, vertical=TRUE, method="jitter", 
           pch=20, col="black", add=TRUE, cex = 0.8, jitter=0.2)
points(1:length(intraBasin_max_xd), sapply(intraBasin_max_xd, median), pch = 16, col = "white")
# Adding the mean
intraBasin_max_mean <- sapply(intraBasin_max_xd, mean)
for(i in seq_along(intraBasin_max_mean)){
  segments(i-0.5, intraBasin_max_mean[i], i+0.5, intraBasin_max_mean[i], "black", lwd=2.5)
  if(i==length(intraBasin_max_mean)) break
  segments(i+0.5, intraBasin_max_mean[i], i+0.5, intraBasin_max_mean[i+1], "black", lwd=2.5, lty = 3)
}

#statistical test
kruskal.test(intraBasin_max_xd)

#### ----------------------------
#### 5. Genetic Distance Between Basins
#### ----------------------------
# Functions
dg_comparator1x1 <- function(basin1, basin2, database = amz_fish_bin, seq_database = amz_seqs,
                             dgp="k2p", alg="min", average=TRUE){
  bas1_BINs <- unique(database$bin_uri[which(database$basin%in%basin1)])
  bas2_BINs <- unique(database$bin_uri[which(database$basin%in%basin2)])
  BINscomun <- intersect(bas1_BINs, bas2_BINs)
  
  distances <- vector()
  for(BIN in BINscomun){
    seqIDs1 <- database$sequenceID[database$basin%in%basin1 & database$bin_uri==BIN]
    seqIDs2 <- database$sequenceID[database$basin%in%basin2 & database$bin_uri==BIN]
    
    seqs1 <- unlist(seq_database[as.character(seqIDs1)])
    seqs2 <- unlist(seq_database[as.character(seqIDs2)])
    
    dist <- polidistance(seqs1, seqs2, FUN = dgp)
    dist <- match.fun(alg)(dist, na.rm = TRUE)
    dist <- 1/(dist+1)
    
    distances <- c(distances, dist)
  }
  names(distances) <- BINscomun
  
  if(average) return(round(mean(distances), 2))
  else return(distances)
}

basin_similarity <- function(basins, database = amz_fish_bin, seq_database = amz_seqs,
                             dgp = "k2p", alg = "min", average = FALSE) {
  comparisons <- combn(basins, 2)
  
  inter_BIN <- lapply(1:ncol(comparisons), function(i) {
    dg_comparator1x1(
      basin1 = comparisons[1, i],
      basin2 = comparisons[2, i],
      database = database,
      seq_database = seq_database,
      dgp = dgp,
      alg = alg,
      average = average
    )
  })
  names(inter_BIN) <- apply(comparisons, 2, function(pair) paste(pair[1], "x", pair[2]))

  if (!average) return(inter_BIN)

  matriz <- matrix(NA, ncol = length(basins), nrow = length(basins),
                   dimnames = list(basins, basins))
  for (i in seq_len(ncol(comparisons))) {
    basin1 <- comparisons[1, i]
    basin2 <- comparisons[2, i]

    matriz[basin1, basin2] <- inter_BIN[[i]]
    matriz[basin2, basin1] <- inter_BIN[[i]]
  }
  
  return(matriz)
}
# dg_comprator:
# performs different comparisons between the basins of the matrices depending on 
# the algorithm used. All are performed between the common BINs of the basins.
# average=TRUE, provide the average; average=FALSE, provide all values

k2dg_intraBIN_min <- basin_similarity(basins, dgp = "k2p", alg = "min", average = TRUE)

green_palette <- colorRampPalette(c("#658a6e", "floralwhite"))(25)

breaks <- seq(0.5, 1, length.out = length(green_palette) + 1)

gplots::heatmap.2(
  k2dg_intraBIN_min_med, 
  tracecol="purple",
  na.rm = TRUE, 
  trace = "none",
  cexRow = 0.8, 
  cexCol = 0.8,
  col = green_palette,
  key.title = "Color Key",
  key.xlab = "Genetic distance",
  key.ylab = "Measure",
  margins = c(5.6, 5.6),
  # breaks = breaks,
  hclustfun = function(x) hclust(x, method = "average")
)

dev.off()

gen_dist_plot(k2dg_intraBIN_min, mincol = "#658a6e",  maxcol="floralwhite")

#### ----------------------------
#### 6. East-West Boundary Analysis
#### ----------------------------
#Functions
region_homogeneity <- function(region1, region2, database=amz_fish_bin, seq_database=amz_seqs, 
                               alg="max", dgp="k2p", minseqs=3){
  BINs1 <- unique(database$bin_uri[database$basin %in% region1])
  BINs2 <- unique(database$bin_uri[database$basin %in% region2])
  BINs <- intersect(BINs1, BINs2)
  
  sufficient_BINs <- sapply(BINs, function(BIN) {
    R1_count <- sum(database$basin %in% region1 & database$bin_uri == BIN)
    R2_count <- sum(database$basin %in% region2 & database$bin_uri == BIN)
    R1_count >= 2 & R2_count >= 2
  })
  BINs <- BINs[sufficient_BINs]
  
  R1_max_dg <- basin_homogeneity(region1, database = database, seq_database = seq_database, 
                                 BINs = BINs, alg = alg, dgp = dgp,
                                 average=FALSE, minseqs=minseqs)
  R2_max_dg <- basin_homogeneity(region2, database = database, seq_database = seq_database, 
                                 BINs = BINs, alg = alg, dgp = dgp,
                                 average=FALSE, minseqs=minseqs)
  RT_max_dg <- basin_homogeneity(basins, database = database, seq_database = seq_database, 
                                 BINs = BINs, alg = alg, dgp = dgp,
                                 average=FALSE, minseqs=minseqs)
  
  return(list("Región 1"=R1_max_dg, 
              "Región 2"=R2_max_dg,
              "Región total"=RT_max_dg,
              "BINs comparados"=BINs
  ))
}

regions <- list(
  list(
    "West"=c("Upper Amazon", "Negro", "Madeira", "Trombetas", "Paru/Jari", "Tapajos"),
    "East"=c("Xingu", "Tocantins")
  ),
  list(
    "West"=c("Upper Amazon", "Negro", "Madeira", "Trombetas", "Paru/Jari"),
    "East"=c("Tapajos", "Xingu", "Tocantins")
  ),
  list(
    "West"=c("Upper Amazon", "Negro", "Madeira", "Trombetas", "Tapajos"),
    "East"=c("Paru/Jari", "Xingu", "Tocantins")
  ),
  list(
    "West"=c("Upper Amazon", "Negro", "Madeira", "Trombetas"),
    "East"=c("Paru/Jari", "Tapajos", "Xingu", "Tocantins")
  ),
  list(
    "West"=c("Upper Amazon", "Negro", "Madeira"),
    "East"=c("Trombetas", "Paru/Jari", "Tapajos", "Xingu", "Tocantins")
  ),
  list(
    "West"=c("Upper Amazon", "Negro"),
    "East"=c("Madeira", "Trombetas", "Paru/Jari", "Tapajos", "Xingu", "Tocantins")
  ),
  list(
    "West"=c("Upper Amazon", "Madeira"),
    "East"=c("Negro", "Trombetas", "Paru/Jari", "Tapajos", "Xingu", "Tocantins")
  )
)

limit_max <- list()
for(i in 1:length(regiones)){
  limit_max[[i]] <- region_homogeneity(region1=regions[[i]][[1]], 
                                       region2=regions[[i]][[2]], 
                                       database = amz_fish_bin, 
                                       seq_database = amz_seqs,
                                       dgp = "k2p", 
                                       alg = "max",
                                       minseqs=2)
}

#Graphic
limit_compare_plot(regions, limit_max, xsize=0.80, ysize=0.79, BINsize = T,
                   col1 = "#bad9a6", col2 = "#658a6e", tracecol = "black")
kruskal.test(lapply(limit_max, "[[", 1))
kruskal.test(lapply(limit_max, "[[", 2))


#### ----------------------------
#### 7. Supplementary Analysis: Common BINs
#### ----------------------------
## Matrix of the number of BINs in common per basin ##
basin_BINs <- lapply(basins, function(basin) {
  unique(amz_fish_bin$bin_uri[amz_fish_bin$basin == basin])
})
names(basin_BINs) <- basins

commonBINs <- matrix(NA, ncol = length(basins), nrow = length(basins),
                     dimnames = list(basins, basins))

for (i in seq_along(basins)) {
  for (j in i:length(basins)) {
    commonBINs[i, j] <- length(intersect(basin_BINs[[i]], basin_BINs[[j]]))
  }
}
commonBINs

#graph
#The upset graph plots only the unique intersections, so a data set is 
#constructed with that structure.
USdata <- matrix(NA, nrow = 0, ncol = length(basins))
colnames(USdata) <- basins
a <- 1
for(i in 1:nrow(commonBINs)){
  for(j in a:ncol(commonBINs)){
    if(i==j) next
    colpos1 <- which(rownames(commonBINs)[i]==colnames(USdata))
    colpos2 <- which(colnames(commonBINs)[j]==colnames(USdata))
    miniUSdata <- matrix(0, nrow = commonBINs[i,j], ncol = length(basins))
    miniUSdata[,c(colpos1, colpos2)] <- 1
    USdata <- rbind(USdata, miniUSdata)
  }
  a <- a+1
}
colnames(USdata)[c(1, 5)] <- c("Upper_Amazon", "Paru_Jari")
combin <- combn(colnames(USdata),2)
intersections <- list()
for(i in 1:ncol(combin)){
  intersections[[i]] <- list(combin[1,i], combin[2,i])
}

USdata <- as.data.frame(USdata)
UpSetR::upset(USdata, 
              intersections = intersections, 
              text.scale = 1.2,
              matrix.color = "#668a6e", 
              main.bar.color = "#c1cf8b", 
              sets.bar.color = "#c1cf8b", 
              order.by = "freq",
              shade.alpha = 1, 
              matrix.dot.alpha = 1,
              sets.x.label = "", 
              mainbar.y.label = "Number of BINs in common")











