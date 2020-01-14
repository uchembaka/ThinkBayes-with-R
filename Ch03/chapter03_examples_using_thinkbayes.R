source("...../thinkbayes.r")
#Chapter 3 exercise using source code from thinkbayes.r


#3.1 Dice Problem

Dice <- setRefClass("Dice",
                    contains = "Suite",
                    methods = list(
                      Likelihood = function(data, hypo){
                        if(as.numeric((hypo)) < as.numeric(data)) return(0)
                        else return(1/as.numeric((hypo)))
                      }
                    )
)#Dice Class
dice <- Dice$new(c('4','6','8','12','20'))
dice$Update('6')
dice$Print()

for (i in c('6','8','7','7','5','4')){
  dice$Update(i)
}
dice$Print()


Train <- setRefClass("Train",
                     fields = list(hypotheses = "numeric"),
                     contains = "Dice",
                     methods = list(
                       initialize = function(hypos){
                         hypotheses <<- as.numeric(hypos)
                         callSuper(.self)
                       },
                       Plot = function(){
                         plot(hypotheses, unlist(.self$ls), type = 'l', 
                              xlab = "Number of trains",
                              ylab = "Probability",
                              col = "blue",
                              lwd = 2)
                       }
                     )
)#Train class

train <- Train$new(1:3)
train$Update('60')
train$Plot()

print(train$Mean())

for(i in c(500, 1000, 2000)){
  train <- Train$new(1:i)
  for (data in c('60', '30', '90')){
    train$Update(data)
  }
  cat(i,": ",train$Mean(),"\n")
}


#Train example using power law prior
Train <- setRefClass("Train",
                     fields = list(hypotheses = "numeric"),
                     contains = "Dice",
                     methods = list(
                       initialize = function(hypos,alpha = 1){
                         hypotheses <<- as.numeric(hypos)
                         callSuper(hypos)
                         .self$PowerLaw(alpha)
                       },
                       PowerLaw = function(alpha = 1){
                         for (val in names(.self$ls)) {
                           .self$Set(val,as.numeric(val)^(-alpha))
                         }
                         .self$Normalize()
                       },
                       Plot = function(){
                         plot(hypotheses, unlist(.self$ls), type = 'l', 
                              xlab = "Number of trains",
                              ylab = "Probability",
                              col = "blue",
                              lwd = 2)
                       }
                     )
)#Train class

#updated example
for(i in c(500, 1000, 2000)){
  train <- Train$new(1:i)
  for (data in c('60', '30', '90')){
    train$Update(data)
  }
  cat(i,": ",train$Mean(),"\n")
}

train <- Train$new(1:1000)
train$Update('60')
train$Plot()
train <- Train$new(1:1000, alpha = 0)
train$Update('60')
lines(unlist(train$ls), col="skyblue", lwd = 2)
legend("topright", legend = c('Uniform Prior', "Power law prior"), fill = c("skyblue","blue"))

#Remember the train used is that of the for loop above with three updates
cdf <- train$MakeCdf()
interval <- c(cdf$Percentile(5),cdf$Percentile(95))
print(interval)
