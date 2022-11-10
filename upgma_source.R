library(ape)

##function to get a vector containing the smallest element of the matrix 
## and its coordinates (col & row)

#makes the distance matrix 
dst <- function(m){
  apply(m, 2, function(x){
    apply(m, 2, function(y){
      sqrt( sum((x-y)^2, na.rm=TRUE) )
    })
  })
}

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
  for (col in c(1:ncol(tab))){
    if (col == values[2] | col == values[3]){
      if (col == values[2]){
      for (row in c(1:nrow(tab))){ #or as a condition do if val > values[1] (min val)
        if (tab[row,col] > values[1]){ #that should work!
          x1 <-append(x1, c(tab[row,col]))
          c_row <- append(c_row, row)
         }
        }
      }
      if (col == values[3]){
      for (row in c(1:nrow(tab))){ #or as a condition do if val > values[1] (min val)
        if (tab[row, col]> values[1]){ #that should work!
          x2 <- append(x2, c(tab[row , col]))
         }
        }
      }
    }else{
      for (row in c(1:nrow(tab))){
        if (row != values[2] & row != values[3]){
          x3 <- append (x3, c(tab[row, col]))
          }
        }
      } 
    }
  return(list(x1, x2, x3, c_row))
}


#upgma <- function(tmp){
  #for (i in (1:c(ncol(tmp)-2))){
    
  #}
 # return(list(family_mat, node_count))
#}

