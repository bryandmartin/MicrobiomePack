##### Start with Day 0, ID A
W <- as.matrix(data[,-1])
# filtering
covars0 <- covars %>% 
  filter(Day == 0) %>%
  filter(ID == "A")
# get the samples of interest
samp0 <- covars0$X1

# extract subset of interest
#Z.0A <- Z.comp[which(S %in% samp0),]
# turn into acomp for compositions package
#Z.0A <- acomp(Z.0A)


W <- W[which(S %in% samp0),]
