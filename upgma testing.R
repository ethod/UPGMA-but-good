source('upgma_source.R')

library(ape)
library(Biostrings)

require("GEOquery")

#initialization

gds <- getGEO('GDS3850')

gds.data <- Table(gds)
gds.meta <- Meta(gds)
gds.columns <- Columns(gds)

gds.exp <-gds.data[, -2:-1]

## rownames(gds.exp)<- gds.data[,1]

gds.exp<- as.matrix(log(gds.exp))

gds.exp.dist <- dst(gds.exp)

image(x=1:nrow(gds.exp.dist), y=1:nrow(gds.exp.dist), gds.exp.dist)

#ref matrix parent child etc

#see source file

##iterative test of vague upgma starts here really &&&&&&&&&&&&&&&&&&&&&&&&&&&&&

for (i in (1:c(ncol(tmp)-2))){
  
  coords<- min.dist_mat(tmp)
  
  new_val <- c( )
  new_val <- val_to_change(tmp, coords)
  
  trnas_mat <- matrix(unlist(new_val[1:2]), ncol = 2, nrow = nrow(tmp)-2)
  
  calculus <-c( )
  
  for (i in c(1:nrow(trnas_mat))){
    calculus <- append(calculus, (((trnas_mat[i,1])+(trnas_mat[i,2]))/2), 
                       after=length(calculus)) # so calc mean of each row
  }
  
  names(calculus) <- names(unlist(new_val))[1:length(calculus)]
  
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
  
  #format (put upper side tri to MAX)
  
  phase2_mat[!lower.tri(phase2_mat)] <-max(phase2_mat+1)
  tmp <- phase2_mat
  
  #temporarily store parent, child & distance value during upgma iteration
  family_mat <- matrix(NA, nrow=length(tmp)*2, ncol=3)
  node_count <- 0
  temp_v =1
  for (j in c(2:3)){
    if (temp_v == 1 | temp_v ==2){
      family_mat[temp_v,] <- c(ncol(tmp) + 1, coords[j], coords[1] / 2 )
      temp_v <- temp_v + 1
    }else{
      family_mat[temp_v,] <- c(family_mat[temp_v-1, 1], coords[j], coords[1] / 2)
      temp_v <- temp_v + 1
    }
    node_count <- node_count + 0.5
  }
  
}


####
#first setup lines (not in upgma loop)

tmp <- read.table("mat_01.txt", stringsAsFactors=FALSE, header=TRUE)
col.ids <- c(1:ncol(tmp))
colnames(tmp)<- as.character(col.ids)
old <- tmp
#stores the original table for reference since it gets replaced at each iteration

family_mat <- matrix(NA, nrow=length(tmp)*2, ncol=3)
node_count <- 0
temp_v =1

names_count <- c((nrow(old)+1))
#"corect" names reference for filling output matrix 

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

names(calculus) <- names(unlist(new_val))[1:length(calculus)]

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
}

coord_as_name <- as.integer(c(colnames(tmp)[coords[2]], colnames(tmp)[coords[3]]))



names_count <- append( names_count, col.ids[c(-(coord_as_name[1]), -(coord_as_name[2]))])
#naming the new mat

colnames(phase2_mat)<- as.character(names_count)


#temporarily store parent, child & distance value during upgma iteration
#good for parent, not for children names idkkkkkk need to code loop of 'oh you 
#were here before, that means you were called 'lqjgh' idk how probably hard 


for (j in c(1:2)){
  if (temp_v == 1 | temp_v ==2){
    family_mat[temp_v,] <- as.integer(c(colnames(phase2_mat)[1], coord_as_name[j], coords[1] / 2 ))
    temp_v <- temp_v + 1
  }else{
    family_mat[temp_v,] <- as.integer(c(colnames(phase2_mat)[1], coord_as_name[j], coords[1] / 2))
    temp_v <- temp_v + 1
  }
  node_count <- node_count + 0.5
}
#inconsistency issue of naming the child since we rewrite the matrix, dunno 
#how to fiw, also parent name gets fucked up 
#format
tmp <- phase2_mat



##upgma(tmp)

##wooo

descendance <- list()

class(descendance) <- 'phylo'

descendance$edge <- matrix (c(family_mat[,1], family_mat[,2]))

descendance$Nnode <- node_count

#this is where i'd put the <- labels IF I HAD ONE

plot.phylo(descendance)



##################### tests

#not 100% sure about this one 
mat_v2 <- tmp
mat_v2[!lower.tri(mat_v2)] <-max(mat_v2+1)

#other version by lars
tmp2 <- tmp
tmp2[ !lower.tri(tmp2) ] <- max(tmp + 1) ##sets values in upper right tri much 
##higher so they won't be detected when searching for min 
tmp2.mi <- which.min(as.matrix(tmp2)) 
tmp2.mc <- 1 + (tmp2.mi-1) %/% ncol(tmp) ##-1 bc of R counting from 1
tmp2.mr <- 1 + (tmp2.mi-1) %% ncol(tmp)

#ref matrix parent child etc
descendance <- matrix(nrow=0, ncol=3)
colnames(descendance) <- c("parent","child", "distance")

#joining algorithm
#function needs to take (input matrix and either output matrix or code it in)

#upgma <- function(df){} 




for (col in c(1:ncol(tmp))){
  x <- append(x, col, after= length(x))
}
x<-c()
x1<- c( )   
       #calc new val

       #don't do shit? 
       # this is going great
  
#to make new matrix -> store data as logical vector then fill in by row of size (og matrix -1 col et 1 row)
#maybe doing new mat [col AB ]<- c(new data vector of 44 x 6 + new data)



coords <- min.dist_mat(tmp)
dist_child <- c(coords[1]/2)
for (val in (nrow(tmp)-2)){
  
  for (i in tmp[coords[2]]){ #so column y (coords: min_v; y;x)
    for (j in (col!=i[coords[1]] & i[coords[2]])){
      
      }
    }
  }
  

#i am thinking

## using ape, use table with 'paren' & 'child' columns (using rbind as the 
##algorithm goes w/ matrix where nrow=0 & ncol=2, with additional column for distance)
##parent of A et B is AB. giving A & B indix (1:ncol of distance matrix)
##parent; child
##6 ; 1
##6; 2

## AB is now indix 6 still 3, 4, 5 for the others, maybe using named vector = 1 to ncol(table) 
##Names (named vector) <- cosl names (table) 



#the error basically means 'you're trying to index a list with a 0 but i'm not??
distmatrix[[2]]
min_value

## it's not doing anything i am giving up for now 


coord<-which(distmatrix== 2 ,arr.ind=TRUE)

#but replacing 2 by the smallest val variable 
#so function that takes coord[[1]], coord[[2]]

#i am going to murder someone



