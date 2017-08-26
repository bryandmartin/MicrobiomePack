## Set up environment
rm(list=ls())
setwd("/Users/Bryan/Google Drive/Daniela/Microbiome")

## Packages
library(readr)
library(tibble)
library(dplyr)
library(compositions)
## Read in data
covars <- as_data_frame(read_csv("~/Google Drive/Daniela/Microbiome/whitman-covariateinfo.csv"))
data <- as_data_frame(read_csv("~/Google Drive/Daniela/Microbiome/whitman-data.csv"))

### I know this is inefficient, practicing with tibble because I want to be more like Hadley Wickham.
## Get data with rows named as samples, columns named as OTUs
# OTU names
colns <- c("S",as.character(data$X1))
# S names
S <- colnames(data)[-1]
# (rows=S,cols=OTU)
data <- as_data_frame((t(data[,-1])))
# Add column of samples
data <- add_column(data,S,.before=TRUE)
colnames(data) <- colns
inherits(data,"tbl_df")
######## FUNCTION ########
#### Helper function to make data compositional
## Assumes that compositions (samples) are across rows
makeComp <- function(X) {
  return(X/rowSums(X))
}

######## FUNCTION ########
#### Get pertubation (by p) of a row Wi for apply
getPurt <- function(Wi,base,p=0.05) {
  zeros <- which(Wi==0)
  nz <- length(zeros)
  # nonzeros <- which(Wi != 0)
  Q <- length(Wi)
  Z.purt <- rep(0,Q)
  if(nz != 0) {
    Z.purt[zeros] <- (Wi[zeros]+p) / (sum(Wi)+p*nz)
    Z.purt[-zeros] <- (Wi[-zeros]) / (sum(Wi)+p*nz)
  } else {
    Z.purt <- Wi / sum(Wi)
  }

  # Return Y.purt
  return(log(Z.purt[-base]/Z.purt[base]))
}

######## FUNCTION ########
#### Get logratios, WITH PERTURBATIONS
# Takes in: Raw counts W, base index D
logratios <- function(W,base,p=0.05) {
  W <- as.matrix(W)
  # get purturbed Y, apply returns arguments as columns
  # CHECK: Apply forces transpose. Worth it?
  Y.purt <- t(apply(W,1,getPurt,base=base,p=p))
  attr(Y.purt,"center")=apply(Y.purt,2,mean)
  return(Y.purt)
}


#### Recall: W is raw data
#### Y is logratios

W <- as.matrix(data[,-1])

whichpart <- function(x, n=75) {
  nx <- length(x)
  p <- nx-n
  xp <- sort(x, partial=p)[p]
  return(which(x > xp))
}

W <- W[,whichpart(colSums(W),n=75)]

## Fornow, debugging, lets take top 75


#Y.p <- logratios(W,base=100,p=0.05)


