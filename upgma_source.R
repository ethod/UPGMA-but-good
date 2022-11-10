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


upgma <- function(tmp){
  col.ids <- c(1:ncol(tmp))
  colnames(tmp)<- as.character(col.ids)
  old <- tmp
  #stores the original table for reference since it gets replaced at each iteration
  
  family_mat <- matrix(NA, nrow=(length(tmp)*2)-1, ncol=3)
  node_count <- 0
  temp_v =1
  
  names_count <- c((nrow(old)+1))
  #"corect" names reference for filling output matrix 
  
  for (i in (1:c(ncol(tmp)-2))){
    #initialization of the output matrix that is appended to each run  ends here &&&&&&&&&&&&&&&&&&&&&
    
    coords<- min.dist_mat(tmp)
    #output the minimal value of the matrix and it's coordinates (min, row=x, col=y)
    
    new_val <- c(val_to_change(tmp, coords))
    #outputs the values to recalculate ( in columns x(1) & y(2) ), the values to 
    #keep as is (3) and the row/ col number of the latter (4) 
    
    trnas_mat <- matrix(unlist(new_val[1:2]), ncol = 2, nrow = nrow(tmp)-2)
    #values to recalculate as matrix (1 & 2) 
    
    
    calculus <-c( )
    
    for (i in c(1:nrow(trnas_mat))){
      calculus <- append(calculus, (((trnas_mat[i,1])+(trnas_mat[i,2]))/2), 
                         after=length(calculus)) # so calc mean of each row
    }
    
    calculus <- append(calculus, 0, after = 0)
    
    
    phase2_mat <- matrix ( NA, nrow=(nrow(tmp)-1), ncol=(ncol(tmp)-1))
    
    
    phase2_mat[,1]=calculus
    phase2_mat[1,]=calculus
    
    rest_mat <- matrix(unlist(new_val[3]), ncol=ncol(tmp)-2, nrow=nrow(tmp)-2)
    
    z <- 1
    for (x in( 2: c(nrow(phase2_mat)))){
      for (y in (2: c(ncol(phase2_mat)))){
        phase2_mat[x,y] <-rest_mat[z]
        z <- z+1
      }
    }
    
    #should contain the 'real names' equivalent to current coordinates & leftover val
    #use debug function? 
    #give row & colnames 
    #get rownames/ col of the amtrix -> as.integer(rowname (tmp)[2])
    #col.ids if we merge then col.ids <- c(max(col.ids +1), col.ids[c(-a;-b)])
    #parent id separated? 
    #to acces the data using 'fake names ' -> do as character (safer version of 
    #things though not very pretty)
    #actually claing them nodes id
    
    if (temp_v >= 2){
      names_count <- c(as.integer(colnames(tmp)[1])+1)
      col.ids <- c(as.integer(colnames(tmp)))
    }
    
    coord_as_name <- c(as.integer(colnames(tmp)[coords[2]]),as.integer
                       (colnames(tmp)[coords[3]]))
    
    
    for (num in c(1: length(col.ids))){
      if (col.ids[num] != coord_as_name[1] & col.ids[num] != coord_as_name[2]){
        names_count <- append(names_count, col.ids[num] )
      }
    }
    
    #= 7 3 4 
    #naming the new mat
    
    colnames(phase2_mat)<- as.character(names_count)
    
    
    #temporarily store parent, child & distance value during upgma iteration
    #good for parent, not for children names idkkkkkk need to code loop of 'oh you 
    #were here before, that means you were called 'lqjgh' idk how probably hard 
    
    
    for (j in c(1:2)){
      if (temp_v == 1 | temp_v ==2){
        family_mat[temp_v,] <- c(as.integer(colnames(phase2_mat)[1]), 
                                 as.integer( coord_as_name[j]), coords[1] / 2 )
        temp_v <- temp_v + 1
      }else{
        family_mat[temp_v,] <- c(as.integer(colnames(phase2_mat)[1]), 
                                 as.integer(coord_as_name[j]), coords[1] / 2 ) 
        temp_v <- temp_v + 1
      }
      node_count <- node_count + 0.5
    }
    #inconsistency issue of naming the child since we rewrite the matrix, dunno 
    #how to fiw, also parent name gets fucked up 
    #format
    tmp <- phase2_mat
    
  }
  return(list(family_mat, node_count))
}

