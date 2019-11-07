#to be able to run this script you need to have told the Maxent model to produce background predictions. If you are running MaxEnt in R this means putting the argument (after "args") "writebackgroundpredictions=true" as true not false. 

#put the working directory as the folder where the maxent results are 

#read in the background predictions

list.files()

#read in the file 
backgroundpredictions <- read.csv("species_backgroundPredictions.csv")

#we need the last column so will set the number as x
x <- length(backgroundpredictions)

#extract the cloglog/logistic results
backgroundclog <- backgroundpredictions[,x]

#now read in the sample predictions for testing
samplepredictions <- read.csv("species_samplePredictions.csv")

#we need the last column again of logistic or cloglog predictions so set a second x
x2 <- length(samplepredictions)

#extract the cloglog/logistic results for sample
sampleclog <- samplepredictions[,x2]

#set n the number of pseuabsences used for backgroudn predictions by MaxEnt
n <- 10000

#set threshold value
th <- 0.0972

TSS_calculations <- function (sample_clog, prediction_clog, n, th) {
  
  xx <- sum(sample_clog > th)
  yy <- sum(prediction_clog > th)
  xxx <- sum(sample_clog < th)
  yyy <- sum(prediction_clog < th)
  
  ncount <- sum(xx,yy,xxx,yyy)
  
  overallaccuracy <- (xx + yyy)/ncount 
  sensitivity <- xx / (xx + xxx)
  specificity <- yyy / (yy + yyy)
  tss <- sensitivity + specificity - 1
  
  #kappa calculations
  a <- xx + xxx
  b <- xx + yy
  c <- yy + yyy
  d <- xxx + yyy
  e <- a * b
  f <- c * d
  g <- e + f
  h <- g / (ncount * ncount)
  hup <- overallaccuracy - h
  hdown <- 1 - h
  
  kappa <- hup/hdown
  Po <- (xx + yyy) / ncount
  Pe <- ((b/ncount) * (a/ncount)) + ((d/ncount) * (c/ncount))
  Px1 <- Po - Pe
  Px2 <- 1 - Pe
  Px3 <- Px1/Px2
  
  tx1 <- xx + yyy
  tx2 <- 2 * a * c
  tx3 <- a - c
  tx4 <- xx - yyy
  tx5 <- ncount * (( 2 * a ) - tx4)
  tx6 <- ncount * tx1
  
  kappamax <- (tx6 - tx2 - (tx3 * tx4)) / ((tx5 - tx3) - (tx3 * tx4))
  
  cat(" Maxent results for model with\n",a,"test sample predictions\n",c ,"background predicitons\n\n TSS value:        ", tss,"\n Overall accuracy: ",overallaccuracy,"\n Sensitivity:      ",sensitivity,"\n Specificity:      ",specificity,"\n Kappa:            ",kappa,"\n Kappa max:        ",kappamax)
  

}


#run the function, the input values are the sampleclog values, then the background clog values, the sample number for the pseudo absences and then threshold value
TSS_calculations(sampleclog,backgroundclog,n,th)








