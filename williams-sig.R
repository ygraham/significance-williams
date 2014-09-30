# Copyright 2014 Yvette Graham 
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Lesser General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/> 

args <- commandArgs(T);
library(psych)

LP <- args[1]
HUMAN.FN <- args[2]
METRICS.FN <- args[3]
WRITE.DIR <- args[4]

if( (((LP=="") | (!file.exists(HUMAN.FN))) | (!file.exists(METRICS.FN))) | (WRITE.DIR=="") ){
  cat( paste( "usage: <lang pair> <human scores filename> <metric scores filename> <write directory>\n"))
  quit("no")
}

dir.create(WRITE.DIR,showWarnings=F)

CORR.FN <- paste(WRITE.DIR,"/pearson-corr.",LP,sep="")
MATRIX.FN <- paste(WRITE.DIR,"/williams-results.",LP,sep="")

h <- read.table( HUMAN.FN, header=T, colClasses=c('factor','factor','factor','factor', 'numeric'))
h <- h[ which( (h$METRIC=="HUMAN") & (h$LP==LP) ), ]

if( length( h$LP ) < 2 ){
  cat(paste("Error, too few human scores found for systems\n"))
  cat(paste("You should provide more MT systems scored by humans, did you specify the correct language pair? (see README file)\n"))
  quit("no")
} 

a <- read.table( METRICS.FN, header=T, colClasses=c('factor','factor','factor','factor','numeric'))
a <- a[ which( a$LP==LP ), ]

if( length( h$LP ) < 4 ){
  cat(paste("Error, too few metric scores found for systems\n"))
  cat(paste("You should provide more MT systems scored by at least two metrics\n"))
  quit("no")
} 

h.scrs <- c()

for(s in sort(unique(unlist( a$SYSTEM ))) ){
  sys <- h[ which( h$SYSTEM==s), ]
  if( length(sys$SCORE) != 1 ){
    cat( paste( "Error with human MT system scores, too few:",  length(sys$SCORE) ,"\n"))
    quit("no")
  }
  h.scrs <- c(h.scrs, sys$SCORE)
}

N <- length(h.scrs)

sig.c <- 0

sink(CORR.FN)
cat(paste("# --------------------------------------------------------------------\n"))
cat(paste("# Pearson Correlation with Human Scores\n"))
cat(paste("# Language Pair:",LP,"\n"))
cat(paste("# --------------------------------------------------------------------\n"))

# print matrix of results
for( m1 in sort( unique( unlist( a$METRIC ) ) ) ){
    
  cat(paste(m1),"") 
  
  m1.scrs <- c() 
  
  for( s in sort( unique( unlist( a$SYSTEM ) ) ) ){
    sys <- a[ which( a$METRIC==m1 & a$SYSTEM==s ), ]

    if( length(sys$SCORE) != 1 ){
      sink()
      cat( paste( "Error with metric scores, each MT system must be scored exactly once by each metric.\n"))
      cat( paste( "Number of scores for MT system", s," by metric ",m1,":", length(sys$SCORE),"\n"))
      quit("no")
     }

    m1.scrs <- c( m1.scrs, sys$SCORE)

  }

  cat(paste( abs( cor( h.scrs,  m1.scrs, method="pearson")),"\n"))
}
sink()


sink(MATRIX.FN)

cat(paste("# -----------------------------------------------------------------------------------\n"))
cat(paste("# William's Test Results \n"))
cat(paste("# -----------------------------------------------------------------------------------\n"))
cat(paste("# Results read as follows:\n"))
cat(paste("#    If the p-value in any given cell is lower than, for example, 0.05,              \n"))
cat(paste("#    the correlation with human scores of the metric named leftmost in that ROW \n"))
cat(paste("#                           is significantly higher than\n"))
cat(paste("#    the correlation with human scores of the metric named at the top of that COLUMN.\n"))
cat(paste("# -----------------------------------------------------------------------------------\n"))

# print matrix of results
for( m1 in sort( unique( unlist( a$METRIC ) ) ) ){
  cat(paste( "\t", m1, sep="" ))
}
cat(paste("\n"))

for( m1 in sort( unique( unlist( a$METRIC ) ) ) ){
    
  cat(paste(m1)) 
  
  m1.scrs <- c() 
  
  for( s in sort( unique( unlist( a$SYSTEM ) ) ) ){
    sys <- a[ which( a$METRIC==m1 & a$SYSTEM==s ), ]

    if( length(sys$SCORE) != 1 ){
      sink()
      cat( paste( "Error with metric scores, each MT system must be scored exactly once by each metric.\n"))
      cat( paste( "Number of scores for MT system", s," by metric ",m1,":", length(sys$SCORE),"\n"))
      quit("no")
     }

    m1.scrs <- c( m1.scrs, sys$SCORE)

  }

  for( m2 in sort( unique( unlist( a$METRIC ) ) ) ){
        
    cat(paste("\t"))

    if( m1==m2 ){

      cat( paste( "-" ) ) 

    }else{

      m2.scrs <- c()
        
      for( s in sort( unique( unlist( a$SYSTEM ) ) ) ){
        sys <- a[ which( a$METRIC==m2 & a$SYSTEM==s ), ]

        if( length(sys$SCORE) != 1 ){
          sink()
          cat( paste( "Error with metric scores, each MT system must be scored exactly once by each metric.\n"))
          cat( paste( "Number of scores for MT system", s," by metric ",m2,":", length(sys$SCORE),"\n"))
          quit("no")
        }
        
        m2.scrs <- c( m2.scrs, sys$SCORE)
      }

      # Perform a one-tailed test for this pair of metrics and human scores
      # Tests if the increase in correlation between m1 (and human judgment) and m2 (and human judgment) is significant

      if( abs( cor( h.scrs,  m1.scrs, method="pearson")) > abs( cor( h.scrs,  m2.scrs, method="pearson")) ){

        p <- r.test( 
          n = N, 
  	  r12 = abs( cor( h.scrs,  m1.scrs, method="pearson")), 
 	  r13 = abs( cor( h.scrs,  m2.scrs, method="pearson")), 
	  r23 = abs( cor( m1.scrs, m2.scrs, method="pearson")), 
 	  pooled=T, twotailed=F)$p

          cat( paste( p, sep="" ) )
      }else{
        cat( paste( "-", sep="" ) )
      }     
    }
  }
  cat(paste("\n"))
}

sink()
cat( paste( "Pearson Correlations of Metrics with Human Scores written to: ", CORR.FN, "\n"))
cat( paste( "Significance Test Results written to: ", MATRIX.FN, "\n"))
