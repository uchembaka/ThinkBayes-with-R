source("......./thinkbayes.r")

Euro <- setRefClass("Euro",
                    contains = "Suite",
                    methods = list(
                      Likelihood = function(data, hypo){
                        x <- as.numeric(hypo)
                        if(data == 'H') return (x/100)
                        else return(1-x/100)
                      },
                      Plot = function(color = "blue"){
                        plot(.self$Render(), type = 'l',
                             xlab = "x",
                             ylab = "Probability",
                             col = color,
                             lwd = 2,
                             ylim = c(0,max(.self$Render()[,2]))
                             )
                      }
                    )
)#Euro class

euro <- Euro(0:101)
dataset <- append(rep('H', 140), rep('T', 110))

for(data in dataset){
  euro$Update(data)
}


euro$Plot()
legend("topright", legend = c("Uniform Prior"), fill = c("blue"))

euro$MaximumLikelihood()#56
cat('Mean:', euro$Mean())
cat('Median:', Percentile(euro, 50))
cat("CI: ", CredibleInterval(euro, 90))

TrianglePrior <- function(){
  suite <- Euro()
  for(x in 0:50){
    suite$Set(x,x)
  }
  for(x in 51:101){
    suite$Set(x,100-x)
  }
  suite$Normalize()
  return(suite)
}


euro2 <- TrianglePrior()
euro2$Plot("skyblue")
euro <- Euro(0:101)
lines(euro$Render(), col = "blue", lwd =2)
legend("topright", legend = c("Uniform","triangle"), fill = c("blue", "skyblue"))

for(data in dataset){
  euro2$Update(data)
  euro$Update(data)
}
euro$Plot()
lines(euro2$Render(), col = "skyblue")
legend("topright", legend = c("Uniform","triangle"), fill = c("blue", "skyblue"))


Euro <- setRefClass("Euro",
                    contains = "Suite",
                    methods = list(
                      Likelihood = function(data, hypo){
                        x <- as.numeric(hypo)/100
                        heads <- data[1]
                        tails <- data[2]
                        like <- x^heads * (1-x)^tails
                        return(like)
                      },
                      Plot = function(color = "blue"){
                        plot(.self$Render(), type = 'l',
                             xlab = "x",
                             ylab = "Probability",
                             col = color,
                             lwd = 2,
                             ylim = c(0,max(.self$Render()[,2]))
                        )
                      }
                    )
)#Euro class: with updated likelihood

data <- c(140, 110)
euro <- Euro(0:101)
euro$Update(data)
euro$Plot()


#using beta class form thinkbayes.r
beta <- Beta()
beta$Update(c(140, 110))
beta$Mean()

pmf <- beta$MakePmf()
pmf$Print()
