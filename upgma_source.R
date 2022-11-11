library(ape)


init_gds <- function(seqq){
  gds.data <- Table(seqq)
  gds.meta <- Meta(seqq)
  gds.columns <- Columns(seqq)
  gds.exp <-gds.data[, -2:-1]
  gds.exp<- as.matrix(log(gds.exp))
  gds.exp.dist <- dst(gds.exp)
  
  return (gds.exp.dist)
}

#makes the distance matrix 
dst <- function(m){
  apply(m, 2, function(x){
    apply(m, 2, function(y){
      sqrt( sum((x-y)^2, na.rm=TRUE) )
    })
  })
}

##function to get a vector containing the smallest element of the matrix 
## and its coordinates (col & row)
min.dist_mat <-function(tabll){
  min_v <- max(tabll)
  min_v.c <- 0
  min_v.r <- 0
  for(row in 2:nrow(tabll)){
    for(col in 1:(row-1)){
      if(tabll[row,col] < min_v){
        min_v <- tabll[row,col]
        min_v.r <- row
        min_v.c <- col
      }
    }
  }
  coord<-c(min_v, min_v.r,min_v.c)
  return(coord)
}

val_to_change <-function(tab, values){
  x1 <-c( )
  x2 <- x1
  x3 <- x1
  c_row <- x1
  for (col in c(1:ncol(tab))){ # parsing through the columns (could also be rows)
    if (col == values[2] | col == values[3]){ # only in cols with the minimum 
      if (col == values[2]){ #separating the values per column for ease of calculus
        for (row in c(1:nrow(tab))){ 
          if (tab[row,col] > values[1]){ 
            x1 <-append(x1, c(tab[row,col]))
            c_row <- append(c_row, row) # keeping track of the rows indexes 
          }
        }
      }
      if (col == values[3]){
        for (row in c(1:nrow(tab))){ 
          if (tab[row, col]> values[1]){ 
            x2 <- append(x2, c(tab[row , col]))
          }
        }
      }
    }else{
      for (row in c(1:nrow(tab))){ # storing the values of the unafected columns
        if (row != values[2] & row != values[3]){
          x3 <- append (x3, c(tab[row, col]))
        }
      }
    } 
  }
  return(list(x1, x2, x3, c_row))
}

upgma <- function(tmp){
  col.ids <- c(1:ncol(tmp))
  colnames(tmp)<- as.character(col.ids)
  
  family_mat <- matrix(NA, nrow=(ncol(tmp)*2)-4, ncol=3) #initialization of output mat
  node_count <- 0
  temp_v =1
  
  names_count <- c((nrow(tmp)+1))
  #beginning of the "correct" names reference for the next matrix, here the name 
  #of the first merged node (always placed in the beginning)
  
  for (i in (1:c(ncol(tmp)-2))){
    
    coords<- min.dist_mat(tmp)
    #output the minimal value of the matrix and it's coordinates 
    #(min, row=x, col=y)
    
    new_val <- c(val_to_change(tmp, coords))
    #outputs the values to recalculate ( in columns x(1) & y(2) )
    # the values to keep as is (3) and the row/ col number of the latter (4) 
    
    trnas_mat <- matrix(unlist(new_val[1:2]), ncol = 2, nrow = nrow(tmp)-2)
    #values to recalculate as matrix (1 & 2) 
    
    calculus <-c( ) #reinitialized at each iteration
    
    for (i in c(1:nrow(trnas_mat))){
      calculus <- append(calculus, (((trnas_mat[i,1]) + (trnas_mat[i,2]))/2), 
                         after=length(calculus)) 
    } #calculating the new values of the sequences' distances to the merged node (=mean) 
    
    calculus <- append(calculus, 0, after = 0)
    
    phase2_mat <- matrix ( NA, nrow =(nrow(tmp)-1), ncol =(ncol(tmp)-1))
    #creation of the new matrix 
    
    phase2_mat[,1] <- calculus
    phase2_mat[1,] <- calculus
    
    rest_mat <- matrix(unlist(new_val[3]), ncol = ncol(tmp)-2, nrow = nrow(tmp)-2)
    
    z <- 1
    for (x in( 2: c(nrow(phase2_mat)))){
      for (y in (2: c(ncol(phase2_mat)))){
        phase2_mat[x,y] <-rest_mat[z]
        z <- z+1
      }
    } # fills in the unchanged values into the new matrix
    
    if (temp_v >= 2){
      names_count <- c(as.integer(colnames(tmp)[1])+1)
      col.ids <- c(as.integer(colnames(tmp)))
    }
    #only skipped in the first iteration of the loop where the names are already correct
    
    coord_as_name <- c(as.integer(colnames(tmp)[coords[2]]),as.integer
                       (colnames(tmp)[coords[3]]))
    # 'translates' the matrix's minimum's coordinates into the "correct" names
    #by "correct" i mean the names that will be used in the output matrix that 
    #stay coherent over every iteration instead of being overwritten. 
    
    for (num in c(1: length(col.ids))){
      if (col.ids[num] != coord_as_name[1] & col.ids[num] != coord_as_name[2]){
        names_count <- append(names_count, col.ids[num] )
      }
    }#appends the names of the columns that won't be merged to the next matrix's 
    
    colnames(phase2_mat)<- as.character(names_count)
    #setting the new names to the new matrix
    
    for (j in c(1:2)){
      if (temp_v == 1 | temp_v ==2){
        family_mat[temp_v,] <- c(as.integer(colnames(phase2_mat)[1]), 
                                 as.integer( coord_as_name[j]), coords[1] / 2)
        temp_v <- temp_v + 1
      }else{
        family_mat[temp_v,] <- c(as.integer(colnames(phase2_mat)[1]), 
                                 as.integer(coord_as_name[j]), coords[1] / 2) 
        temp_v <- temp_v + 1
      }
      node_count <- node_count + 0.5
    }#making the descendants matrix (parent, child , distance) using the "correct" names 
    #set previously 
    
    tmp <- phase2_mat
  }
  return(list(family_mat, node_count))
}
