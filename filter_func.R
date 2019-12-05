
## ------ divide each taxonomic level ------##


## ------ noise removal --------------------##
## ------first, for species level ----------##
noise.removal <- function(data, percent = 0.01, percent2 = 0.001, top=NULL){
  ## data are relative abundances
  ## rows are species, columns are samples ##
  
  ## species with max relative abundance <0.1% in all samples and 
  ## species whose average abundance across all samples below 0.01% were removed.
  Matrix <- data
  bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent 
  Matrix_1 <- Matrix[bigones,]
  print(percent)
  rowmax <- apply(Matrix_1, 1, max)
  bigones2 <- rowmax >= percent2 
  print(percent2)
  Matrix_2 <- Matrix_1[bigones2,]
  return(Matrix_2)
}

noise.removal.path <- function(data, percent = 1e-06, top=NULL){
	## data are relative abundances
	##1. filter pathways with zero abundance in >15% of samples
	##2. filter pathways whose max abundance are lower than 1e-06

  Matrix <- data
  zero.num <- apply(Matrix, 1, function(data){y <- length(data[data ==0])})
  col.num <- ncol(data)
  bigones <- zero.num<=0.15*col.num
  Matrix_1 <- Matrix[bigones,]
  
  rowmax <- apply(Matrix_1, 1, max)
  bigones2 <- rowmax >= percent 
  Matrix_2 <- Matrix_1[bigones2,]
  
  return(Matrix_2)
}

