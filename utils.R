flattenlist <- function(x){
  morelists <- sapply(x, function(k) class(k)[1]=="list")
  out <- c(x[!morelists], unlist(x[morelists], recursive=FALSE))
  if(sum(morelists)){ 
    Recall(out)
  }else{
    return(out)
  }
}#removes the lists from lists and places them in a single list

intermean <- function(x){
  media <- vector()
  for(i in 2:length(x)){
    media <- c(media, mean(c(x[i-1], x[i])))
  }
  return(media)
}#returns the average of every second element of a vector

limit_compare_plot <- function(regions, results=NULL, BINsize=FALSE, basins=NULL,
                               col1="darkolivegreen1", col2="#0c5329", col3="gray30",
                               border1=FALSE, border2=FALSE, border3=FALSE,
                               tracecol="purple", pch.cex = 2, ysize=0.75, xsize=0.75,
                               sort=FALSE){
  if(is.null(basins)){
    cuencas <- c("Upper Amazon","Negro","Madeira","Trombetas",
                 "Paru/Jari", "Tapajos","Xingu","Tocantins")
  }else cuencas <- basins
    
  par_bu <- par()
  
  values <- list()
  totBINs <- vector()
  for(i in 1:length(regions)){
    values[[i]] <- t(data.frame(Region.1 = mean(results[[i]][[1]]), 
                                Region.2 = mean(results[[i]][[2]]), 
                                Region.total = mean(results[[i]][[3]])))
    totBINs <- c(totBINs, length(results[[i]][[4]]))
  }
 
  if(sort){
    orden <- order(unlist(lapply(results, "[[", 3)))
    values <- values[orden]
    regions <- regions[orden]
  }
  
  #defining the limits of the graph
  par(mar = c(3,0.5,0.5,1))
  par(mgp = c(1.5,0.5,0))
  nombres_posy <- c(length(regions), length(regions)/ysize)
  nombres_posx <- c(max(unlist(values))*(xsize-1), 0)
  
  rangox <- c(nombres_posx[1],1.05*max(unlist(values)))
  rangoy <- c(0, nombres_posy[2])
  plot(0, 0, type="n", xlab="value", ylab="", ylim=rangoy, cex.lab=0.8,
       xlim=rangox, xaxs="i", yaxs="i", axes=F)
  axis(side = 1, at = pretty(rangox)[pretty(rangox)>=0], 
       labels = pretty(rangox)[pretty(rangox)>=0], las=1, cex.axis=0.8)
  segments(0,0,1.05*max(unlist(values)),0)
  
  # abline(h=nombres_posy[1], col="red")
  # abline(v=nombres_posx[2], col="red")
  
  #Placing the name of the basins
  texto_poslim <- seq(nombres_posx[1], nombres_posx[2], length.out = length(cuencas)+1)
  
  # abline(v=texto_poslim, col="gray")
  
  for(i in 1:length(texto_pos)){
    text(x=intermean(texto_poslim)[i], y=nombres_posy[1]+0.2, labels=cuencas[i], srt=90, adj=0, cex = 0.8)
  }
  
  # abline(h=1:(nombres_posy[1]-1), col="gray30")
  
  # abline(h=0.5:(nombres_posy[1]-0.5), col="gray50")
  
  #placing rectangles below the dots as decoration
  texto_poslim <- rev(texto_poslim)
  rectsel <- 1:length(cuencas)
  rectsel <- rectsel[1:length(cuencas)%%2 !=0]
  for(i in rectsel){
    rect(texto_poslim[i+1], 0, texto_poslim[i], nombres_posy[1], col="gray90", border = F)
  }
  
  #Marking the points that correspond to the different boundaries
  texto_poslim <- rev(texto_poslim)
  for(i in 0.5:(length(regions)-0.5)){
    points(intermean(texto_poslim), rep(i, length(cuencas)), 
           pch=19, cex = pch.cex, col="gray80")
  }
  for(i in 0.5:(length(regions)-0.5)){
    oeste <- cuencas %in% regions[[i+0.5]][[1]]
    este <- cuencas %in% regions[[i+0.5]][[2]]
    points(intermean(texto_poslim)[oeste], rep(i, sum(oeste)), 
           pch = 21, cex = pch.cex, bg = col1, col = border1)
    points(intermean(texto_poslim)[este], rep(i, sum(este)), 
           pch = 21, cex = pch.cex, bg = col2, col = border2)
  }
  
  #Placing the bars
  for(i in 0:(length(regions)-1)){
    posy <- seq(i+0.15, i+0.85, length.out = 4)
    posx <- rev(values[[i+1]][,1])
    colores <- c(col3, col2, col1)
    borders <- list(border3, border2, border1)
    for(j in 1:(length(posy)-1)){
      rect(0, posy[j], posx[j], posy[j+1], col=colores[j], border=borders[[j]])
    }
  }
  segments(0,0,0,nombres_posy[1], "black")
  
  #placing the number of BINs by comparison
  if(BINsize){
    totBINs_pos <- (totBINs-0)/(max(totBINs)-0) * (0.8*max(unlist(values))-0)+0
    # totBINs_pos <- rev(totBINs_pos)
    for(i in 1:length(totBINs)){
      segments(totBINs_pos[i], i-1, totBINs_pos[i], i, col=tracecol, lwd = 2)
      if(i==length(totBINs)) next
      segments(totBINs_pos[i], i, totBINs_pos[i+1], i, col=tracecol, lty = 2)
    }
    
    ec <- lm(unique(sort(totBINs))~unique(sort(totBINs_pos)))
    labels <- round(ec$coefficients[2]*pretty(rangox)[pretty(rangox)>=0]+ec$coefficients[1], 0)
    
    labels[1] <- "" #placing the 0 further to the right
    label1pos <- ((as.numeric(labels[length(labels)])*0.5/100)-ec$coefficients[1])/ec$coefficients[2]
    axis(3, at = label1pos, labels = 0, tck = 0,
         las=1, cex.axis=0.8, pos=nombres_posy[1])
    
    axis(3, at = pretty(rangox)[pretty(rangox)>=0], labels = labels, 
         las=1, cex.axis=0.8, pos=nombres_posy[1])
    text(x=mean(rangox), cex = 0.8, labels="Number of compared BINs",
         y=nombres_posy[1]+(nombres_posy[2]-nombres_posy[1])*0.35)
  }
  
  on.exit(par(par_bu), add = T)
}#returns a graph to compare metrics at different east-west boundaries

gen_dist_plot <- function(gd, arrange = NULL, num_lett = NULL, matrixmin=TRUE,
                          mincol="darkolivegreen1", maxcol="#0c5329"){
  if(is.null(arrange)) arrange <- rownames(gd)
  
  #obtaining color scale
  allgd <- vector()
  for(i in 1:(nrow(gd)-1)){
    allgd <- c(allgd, gd[i,(i+1):nrow(gd)])
  }
  allgd <- as.vector(allgd)
  if(matrixmin) allgd <- c(allgd, min(gd, na.rm=T))
  
  scaled_data <- scale(allgd, center = min(allgd), scale = max(allgd) - min(allgd))
  
  # color_palette <- colorRampPalette(c("#dbf3ff", "black"))(100)
  # color_palette <- colorRampPalette(c("greenyellow", "#0c5329"))(100)
  color_palette <- colorRampPalette(c(mincol, maxcol))(100)
  
  rango <- ceiling(scaled_data*100)
  rango[which(rango==0)] <- 1
  
  colors <- color_palette[rango]
  
  col_df <- data.frame("gd"=allgd,
                       "color"=colors)
  
  #graph
  par(mar=c(0, 0, 0, 0))
  plot(0,type="n", axes=F, xlab="", ylab="", xaxs="i", yaxs="i",
       xlim = c(0, 1), ylim = c(0, 1))
  
  vdiv <- seq(0, 1, length.out=nrow(gd)+1)
  hdiv <- seq(1/10, 1, length.out=nrow(gd))
  
  for(i in 1:nrow(gd)){
    text((vdiv[i]+vdiv[i+1])/2, 1/20, rownames(gd)[i], adj = c(0.5, 0.5), cex = 0.8)
  }
  
  for(i in 1:nrow(gd)){
    subCue <- setdiff(arrange, rownames(gd)[i])
    
    coldf_basin <- data.frame("cuenca"=subCue,
                              "gd"=NA,
                              "color"=NA)
    for(j in 1:(nrow(gd)-1)){
      coldf_basin[j,2] <- ifelse(is.na(gd[rownames(gd)[i], subCue[j]]), yes = gd[subCue[j], rownames(gd)[i]],
                                 no = gd[rownames(gd)[i], subCue[j]])
      coldf_basin[j,3] <- unique(col_df$color[which(coldf_basin[j,2]==col_df$gd)])
    }
    coldf_basin <- coldf_basin[order(coldf_basin$gd),]
    
    x1 <- (vdiv[i+1]+vdiv[i])/2
    x2 <- (7*vdiv[i+1]+vdiv[i])/8
    for(j in 1:(nrow(gd)-1)){
      y1 <- hdiv[j]
      y2 <- hdiv[j+1]
      rect(x1, y1, x2, y2, col = coldf_basin$color[j], border = F)
    }
    
    if(!is.null(num_lett)) cuenca <- substring(coldf_basin$cuenca, 1, num_lett)
    else cuenca <- coldf_basin$cuenca
    
    for(j in 1:nrow(gd)){
      #In the middle
      # text((11*vdiv[i]+5*vdiv[i+1])/16, (hdiv[j+1]+hdiv[j])/2, cuenca[j], adj = c(1, 0.5), cex = 0.8)
      #to the left
      x <- (vdiv[i]+vdiv[i+1])/2 - (vdiv[i+1]-vdiv[i])*0.05
      y <- (hdiv[j+1]+hdiv[j])/2
      text(x, y, cuenca[j], adj = c(1, 0.5), cex = 0.8)
    }
  }
  abline(h=1/10, col="black")
}# elaborates a graph based on a matrix with genetic distances between watersheds

consensus <- function(sequences, threshold=0.5){
  if(length(unique(nchar(sequences)))!=1) stop("Unaligned sequences")

  sequences <- unlist(sequences)
  seqs_DNAbin <- Biostrings::DNAStringSet(sequences)
  consensus <- DECIPHER::ConsensusSequence(seqs_DNAbin, threshold=threshold)
  consensus <- as.character(consensus)
  
  return(consensus)
} #obtains a consensus sequence from the most common nucleotide

representative <- function(sequences){
  selected <- which.max(sapply(sequences, function(k){
    N <- c("A","T","C","G")
    sum(sapply(N, function(x) sum(unlist(gregexpr(x, k)) != -1)))
  }))
  sequences <- unlist(sequences[selected])
  return(sequences)
} #of a group of sequences, returns the one with the most nucleotides.

`%Lin%` <- function(x, y){
  bool <- sapply(y, function(i) any(x%in%i))
  return(bool)
}#x vector, y list. If any of the elements of x is in an element of y: T 

read.fasta <- function(file, dir = NULL){
  
  if(!is.null(dir)) crudo <- readLines(dir)
  else crudo <- file
  
  nom.pos <- which(startsWith(crudo, ">"))
  nombres <- lapply(crudo[nom.pos], function(k){
    sub(".", "", k)
  })
  
  seqs <- crudo
  seqs[nom.pos] <- 0
  seqs <- paste(seqs, collapse = "")
  seqs <- unlist(strsplit(seqs, split = 0))
  seqs <- seqs[-1]
  
  secuencias <- as.list(seqs)
  names(secuencias) <- unlist(nombres)
  
  return(secuencias)
}

write.fasta <- function(x, file){
  for(i in 1:length(x)){
    nombre <- paste0(">", names(x)[i])
    write(nombre, file = file, append = T)
    
    secuencia <- x[[i]]
    for(j in 1:ceiling(nchar(secuencia)/70)){
      sec <- substr(secuencia, (j-1)*70+1, (j)*70)
      write(sec, file = file, append = T)
    }
  }
}#from a list with sequences creates a fasta

coord_dist <- function(coord1, coord2){
  d <- sqrt((coord1[1] - coord2[1])^2 + (coord1[2] - coord2[2])^2)
  return(d)
}

pdistance <- function(seq1, seq2){
  seq1 <- unlist(strsplit(seq1, split = ""))
  seq1[which(!seq1%in%c("A","T","C", "G"))] <- NA
  seq1[which(seq1=="A")] <- 1
  seq1[which(seq1=="T")] <- 2
  seq1[which(seq1=="C")] <- 3
  seq1[which(seq1=="G")] <- 4
  
  seq2 <- unlist(strsplit(seq2, split = ""))
  seq2[which(!seq2%in%c("A","T","C", "G"))] <- NA
  seq2[which(seq2=="A")] <- -1
  seq2[which(seq2=="T")] <- -2
  seq2[which(seq2=="C")] <- -3
  seq2[which(seq2=="G")] <- -4
  
  coinc <- as.numeric(seq1) + as.numeric(seq2)
  coinc <- coinc[which(!is.na(coinc))]
  p <- sum(coinc!=0)
  seq_len <- length(coinc)
  return(p*100/seq_len)
}

k2pdistance <- function(seq1, seq2){
  seq1 <- unlist(strsplit(seq1, split = ""))
  seq1[which(!seq1%in%c("A","T","C", "G"))] <- NA
  seq1[which(seq1=="A")] <- 1
  seq1[which(seq1=="G")] <- 2
  seq1[which(seq1=="T")] <- -1
  seq1[which(seq1=="C")] <- -2
  
  seq2 <- unlist(strsplit(seq2, split = ""))
  seq2[which(!seq2%in%c("A","T","C", "G"))] <- NA
  seq2[which(seq2=="A")] <- 1
  seq2[which(seq2=="G")] <- 2
  seq2[which(seq2=="T")] <- -1
  seq2[which(seq2=="C")] <- -2
  
  coinc <- as.numeric(seq1) + as.numeric(seq2)
  coinc <- coinc[which(!is.na(coinc))]
  p <- length(which(coinc%in%c(-3, 3)))/length(coinc)
  q <- length(which(coinc%in%c(-1, 0, 1)))/length(coinc)
  d <- -0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q)) * 100
  return(d)
}

polidistance <- function(seqs1, seqs2=NULL, FUN="k2p", names=FALSE){
  if(is.null(seqs2)){
    comb <- combn(1:length(seqs1),  2)
    distances <- vector()
    for(i in 1:ncol(comb)){
      if(FUN=="p"){
        dist <- pdistance(seqs1[comb[1,i]], seqs1[comb[2,i]])
      }else if(FUN=="k2p"){
        dist <- k2pdistance(seqs1[comb[1,i]], seqs1[comb[2,i]])
      }
      if(names) names(dist) <- paste0(names(seqs1[comb[1,i]]), "x", names(seqs1[comb[2,i]]))
      distances <- c(distances, dist)
    }
    return(distances)
  }else{
    comb <- expand.grid(1:length(seqs1), 1:length(seqs2))
    distances <- vector()
    for(i in 1:nrow(comb)){
      if(FUN=="p"){
        dist <- pdistance(seqs1[comb[i,1]], seqs2[comb[i,2]])
      }else if(FUN=="k2p"){
        dist <- k2pdistance(seqs1[comb[i,1]], seqs2[comb[i,2]])
      }
      if(names) names(dist) <- paste0(names(seqs1[comb[i,1]]), "x", names(seqs2[comb[i,2]]))
      distances <- c(distances, dist)
    }
    return(distances)
  }
}#from two sets of sequences, we obtain a vector with the distances of all the sequences from each other.

nuc_div <- function(sequences){
  seqs_matrix <- t(sapply(sequences, function(x) strsplit(x, "")[[1]]))
  seqs_DNAbin <- ape::as.DNAbin(seqs_matrix)
  Ndiv <- pegas::nuc.div(seqs_DNAbin)
  return(Ndiv*100)
}#returns the nucleotide diversity of a sequencing vector



