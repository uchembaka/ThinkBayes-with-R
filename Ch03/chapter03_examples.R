source("....../Ch02/chapter02_examples.R")#link code to chapter 02

#3.1 Dice Problem

Dice <- setRefClass("Dice",
                    contains = "Suite",
                    methods = list(
                      initialize = function(hypos){
                        for (i in hypos) {
                          .self$Set(i,1)
                        }
                        .self$Normalize()
                      },
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
                     contains = "Suite",
                     methods = list(
                       initialize = function(hypos){
                         hypotheses <<- as.numeric(hypos)
                         for (i in hypos) {
                           .self$Set(i,1)
                         }
                         .self$Normalize()
                       },
                       Likelihood = function(data, hypo){
                         if(as.numeric((hypo)) < as.numeric(data)) return(0)
                         else return(1/as.numeric((hypo)))
                       },
                       Plot = function(){
                         plot(hypotheses, unlist(.self$x), type = 'l', 
                              xlab = "Number of trains",
                              ylab = "Probability",
                              col = "blue",
                              lwd = 2)
                       }
                     )
)#Train class

#Example
train <- Train$new(as.character(1:1000))
train$Update('60')
train$Plot()

#Mean method for Pmf Class
#NOTE: For this method to be avialable to Train class you will have to run this method,
#then re-run the Suite class from chapter02, then run the Train class again.
Pmf$methods(Mean = function(){
  hypos <- as.numeric(names(x))
  if(any(is.na(hypos)) == TRUE) return ("This method only works for numeric hypotheses")
  else{
    total = 0
    for (i in 1:length(x)){
      total = total + hypos[i]*x[[i]]
    }
    return (total)
  }
})

print(train$Mean())

for(i in c(500, 1000, 2000)){
  train <- Train$new(as.character(1:i))
  for (data in c('60', '30', '90')){
    train$Update(data)
  }
  cat(i,": ",train$Mean(),"\n")
}

#Train class with power law prior
Train <- setRefClass("Train",
                     fields = list(hypotheses = "numeric"),
                     contains = "Suite",
                     methods = list(
                       initialize = function(hypos, alpha=1.0){
                         hypotheses <<- as.numeric(hypos)
                         for (i in hypos) {
                           .self$Set(i,as.numeric(i)^(-alpha))
                         }
                         .self$Normalize()
                       },
                       Likelihood = function(data, hypo){
                         if(as.numeric((hypo)) < as.numeric(data)) return(0)
                         else return(1/as.numeric((hypo)))
                       },
                       Plot = function(){
                         plot(hypotheses, unlist(.self$x), type = 'l', 
                              xlab = "Number of trains",
                              ylab = "Probability",
                              col = "blue",
                              lwd = 2)
                       }
                     )
)#Train class

#updated example
for(i in c(500, 1000, 2000)){
  train <- Train$new(as.character(1:i))
  for (data in c('60', '30', '90')){
    train$Update(data)
  }
  cat(i,": ",train$Mean(),"\n")
}

train <- Train$new(as.character(1:1000))
train$Update('60')
train$Plot()
train <- Train$new(as.character(1:1000), alpha = 0)
train$Update('60')
lines(unlist(train$x), col="skyblue", lwd = 2)
legend("topright", legend = c('Uniform Prior', "Power law prior"), fill = c("skyblue","blue"))

#Percentile function

Percentile <- function (pmf, percentage){
  p <- percentage/100.0
  total <- 0
  for (i in 1:length(pmf$x)){
    total = total + pmf$x[[i]]
    if(total >= p){
      return (as.numeric(names(pmf$x[i])))
    }
  }
}

#Example
train <- Train$new(as.character(1:1000))
for (data in c('60', '30', '90')){
  train$Update(data)
}

interval <- c(Percentile(train, 5), Percentile(train, 95)) 
print(interval)

#This Cdf class will contain only the necessary methods need for the example as described in chapter 02
#It will be modelled as described in thinkbayes.py
Cdf <- setRefClass("Cdf",
                   fields = list(xs = "vector", ps = "vector", name = "character"),
                   methods = list(
                     initialize = function(xs = NA, ps = NA, name = ""){
                       if(all(is.na(xs)) == TRUE){
                         .self$xs <<- vector()
                       }else{.self$xs <<- xs}
                       if(all(is.na(ps)) == TRUE){
                         .self$ps <<- vector()
                       }else{.self$ps <<- ps}
                       .self$name <<- name
                     },
                     Value = function(p){
                       if(p < 0 || p > 1){
                         stop("Probability p must be in range[0,1]")
                       }
                       if(p == 0) return(.self$xs[1])
                       if(p == 1) return(.self$xs[length(.self$xs)])
                       index <- findInterval(p, ps)#replacing bisect 
                       if(p == .self$ps[index]) return (.self$xs[index])
                       else return (.self$xs[index+1])

                     },
                     Percentile = function(p){
                       return (.self$Value(p/100))
                     }
                   )
)#CDF class

#add new method to PmF class (MakeCdf)
Pmf$methods(MakeCdf = 
              function(name = NA){
  return (MakeCdfFromPmf(pmf = .self, name = name))
})

MakeCdfFromPmf = function(pmf, name = NA){
  return (MakeCdfFromItems(pmf$Items()))
}

MakeCdfFromItems = function(items, name=""){
  runsum <- 0;
  xs <- vector()
  cs <- vector()
  
  items <- sort(items)
  
  for (i in 1:length(items)){
    runsum <- runsum + items[[i]]
    xs <- append(xs, as.numeric(names(items[i])))
    cs <- append(cs, runsum)
  }
  total <- runsum
  ps <- vector()
  for(c in cs){
    ps <- append(ps, c/total)
  }
  cdf <- Cdf$new(xs, ps, name)
  return (cdf)
}
#NOTE: To make the MakeCdf() method available to Suite class, re-run the Suite class from chapter 02,
#then run the second version of Train class again. This is so that Suite class inherits the new Pmf 
#method MakeCdf() and Train inherits the new method from Suite class

cdf <- train$MakeCdf()
interval <- c(cdf$Percentile(5),cdf$Percentile(95))
print(interval)
