source("C:/Users/uchem/Dropbox/Projects/R/Think Bayes/thinkbayes.r")


Die <- setRefClass("Die",
                   fields = list(sides = "numeric"),
                   contains = "Pmf",
                   methods = list(
                     initialize = function(sides){
                       callSuper(1:sides)
                     }
                   ))
d6 <- Die(6)
d6$Print()
dice <- sapply(1:3, function(i) d6)
three <- SampleSum(dice, 1000)
three$Print()
three_exact <- d6$add(d6$add(d6))
three_exact$Print()
plot(three_exact$Render(), type = 'l', col="blue", lty="dashed", lwd=2, ylim = c(0,0.15))
lines(three$Render(), type = 'l',col="skyblue", lwd=2)
legend("topright", legend = c("sample","exact"), fill = c("blue", "skyblue"))


RandomMax <- function(dists){
  total <- max(sapply(dists, function(dist) dist$Random()))
}
SampleMax <- function(dists, n){
  pmf <- MakePmfFromList(sapply(1:n, function(i) RandomMax(dists)))
}


PmfMax <- function(pmf1, pmf2){
  res <- Pmf()
  for(i in names(pmf1$Items())){
    v1 <- as.numeric(i)
    p1 <- pmf1$Items()[[i]]
    for(j in names(pmf2$Items())){
      v2 <- as.numeric(j)
      p2 <- pmf2$Items()[[j]]
      res$Incr(max(v1,v2), p1*p2)
    }
  }
  
  return(res)
}

dice <- sapply(1:3, function(i) three_exact)
three.exact.max <- PmfMax(three_exact, PmfMax(three_exact,three_exact))
three.exact.max$Print()
plot(three.exact.max$Render(),
     type = 'l', 
     col="blue", 
     lwd=2, 
     ylim = c(0,0.20), xlim = c(3,19))

best_attr_cdf <- three_exact$Max(6)
best_attr_pmf <- best_attr_cdf$MakePmf()



#Mixture Model

d6 <- Die(6)
d8 <- Die(8)

mix <- Pmf()

for (die in c(d6, d8)){
  for(outcome in die$Items()){
    mix$Incr(outcome, die$Items()[[outcome]])
  }
}

mix$Normalize()
mix$Print()
d8
pmf.dice <- Pmf()
pmf.dice$Set("Die(4)",5)
pmf.dice$Set("Die(6)",4)
pmf.dice$Set("Die(8)",3)
pmf.dice$Set("Die(12)",2)
pmf.dice$Set("Die(20)",1)

mix <- MakeMixture(pmf.dice)
mix$Normalize()
mix$Print()
library(ggplot2)
ggplot(mix$Render(), aes(x = V1, y = V2)) +
  geom_bar(stat="identity",fill = "skyblue")+
  labs(x  = "outcomes", y = "Probability")

