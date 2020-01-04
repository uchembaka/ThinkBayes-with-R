#Code written as defined in the book chapter 2 (Computational Statistics)
#NOT as in the thinkbayes.py
Pmf <- setRefClass("Pmf", 
            fields = list(x = "list"),
            methods = list(
              Set = function(rv, pr){
                x[[rv]] <<- pr
              },
              Incr = function(word){
                if (word %in% names(x)){
                  x[[word]] <<- x[[word]]+1
                }else{
                  x[[word]] <<- 1
                }
              },
              Normalize = function(){
                sum <- Reduce("+",x)
                for(i in 1:length(x)){
                  x[i] <<- x[[i]]/sum
                }
              },
              Prob = function(rv){
                return(x[[rv]])
              },
              Mult = function(rv, pr){
                x[[rv]] <<- x[[rv]] * pr
              }
            )
)#Pmf Class

#item method to Pmf
Pmf$methods(Items = function(){
  items <- c()
  hypos <- names(x)
  for (name in 1:length(hypos)){
    value <- .self$Prob(hypos[name]);
    items <- append(items, value)
  }
  names(items) <- hypos
  return (items)
})

#Six-sided die example
die <- Pmf()
for(i in 1:6){
  die$Set(i,1/6.0)
}
#print prob
for(i in 1:6){
  print(die$Prob(i))
}

#Word list example
word_list <- c("Think","Bayes","By", "Allen", "Downy", "Rep", "Rep")
word_count <- Pmf()
for(i in 1:length(word_list)){
  word_count$Incr(word_list[i])
}
print(word_count$Prob("Rep"))#before normalizing
word_count$Normalize()
print(word_count$Prob("Rep"))#after normalizing

#Cookie Problem
cookie <- Pmf()
cookie$Set("Bowl2", 0.5)
cookie$Set("Bowl1", 0.5)
cookie$Mult("Bowl2", 0.5)
cookie$Mult("Bowl1", 0.75)
cookie$Normalize()
print(cookie$Prob("Bowl1"))
print(cookie$Prob("Bowl2"))

#The Bayesian framework

Cookie <- setRefClass("Cookie",
                      fields = list(mixes = "list"),
                      contains = "Pmf",
                      methods = list(
                        initialize = function(hypos){
                          mixes <<- list("Bowl1" = list("vanilla" = 0.75, "chocolate" = 0.25),
                                       "Bowl2" = list("vanilla" = 0.5, "chocolate" = 0.5))
                          for (i in hypos) {
                            .self$Set(i,1)
                          }
                          .self$Normalize()
                        }
                      )
                      )
#Cookie class first example
hypos <- c("Bowl1", "Bowl2")
cookie2 <- Cookie$new(hypos)
print(cookie2$Prob("Bowl1"))# should be 0.5

Cookie$methods(Likelihood = function(data, hypo){
  mix <- .self$mixes[[hypo]]
  like <- mix[[data]]
  return (like)
})
Cookie$methods(Update = function(data){
  for(hypo in names(x)){
    like = .self$Likelihood(data, hypo)
    .self$Mult(hypo, like)
  }
  .self$Normalize()
})

cookie2 <- Cookie$new(hypos)#Override previous def
cookie2$Update("vanilla")
cookie2$Items()

#generalization to case with more than one draw with replacement
dataset = c("vanilla", "chocolate", "vanilla")
for (i in dataset){
  cookie2$Update(i)  
}
cookie2$Items()



#Monty Hall Problem
Monty <- setRefClass("Monty",
                     contains = "Pmf",
                     methods = list(
                       initialize = function(hypos){
                         for (i in hypos) {
                           .self$Set(i,1)
                         }
                         .self$Normalize()
                       },
                       Likelihood = function(data, hypo){
                         if (hypo == data){
                           return (0)
                         }else if(hypo == 'A'){
                           return (0.5)
                         }else return (1)
                       },
                       Update = function(data){
                         for(hypo in names(x)){
                           like = .self$Likelihood(data, hypo)
                           .self$Mult(hypo, like)
                         }
                         .self$Normalize()
                       }
                     )
)#Monty Class

#Example
hypos <- c('A','B','C')
monty <- Monty$new(hypos) 
monty$Update('B')
monty$Items()


#Encapsulating the framework
Suite <- setRefClass("Suite",
                     contains = "Pmf",
                     methods = list(
                       initialize = function(){
                         stop("This class is acting as an abstract class. Use any of the subclasses")
                       },
                       Update = function(data){
                         for(hypo in names(x)){
                           like = .self$Likelihood(data, hypo)
                           .self$Mult(hypo, like)
                         }
                         .self$Normalize()
                       },
                       Print = function(){
                         items.name <- c(names(.self$Items()))
                         i <- 1
                         for(item in .self$Items()){
                           cat(items.name[i],": ", item,"\n")
                           i = i+1
                         }
                       }
                     )
)#Suite class

#re-write Cookie class and MOnty Class
Cookie <- setRefClass("Cookie",
                      fields = list(mixes = "list"),
                      contains = "Suite",
                      methods = list(
                        initialize = function(hypos){
                          mixes <<- list("Bowl1" = list("vanilla" = 0.75, "chocolate" = 0.25),
                                         "Bowl2" = list("vanilla" = 0.5, "chocolate" = 0.5))
                          for (i in hypos) {
                            .self$Set(i,1)
                          }
                          .self$Normalize()
                        },
                        Likelihood = function(data, hypo){
                          mix <- .self$mixes[[hypo]]
                          like <- mix[[data]]
                          return (like)
                        }
                      )
)#Cookie Class updated

Monty <- setRefClass("Monty",
                     contains = "Suite",
                     methods = list(
                       initialize = function(hypos){
                         for (i in hypos) {
                           .self$Set(i,1)
                         }
                         .self$Normalize()
                       },
                       Likelihood = function(data, hypo){
                         if (hypo == data){
                           return (0)
                         }else if(hypo == 'A'){
                           return (0.5)
                         }else return (1)
                       }
                     )
)#Monty Class updated

#M&M Problem
MandM <- setRefClass("MandM",
                     fields = list(mix94 = "list", mix96 = "list", hypoA = "list", hypoB = "list", hypotheses = "list"),
                     contains = "Suite",
                     methods = list(
                       initialize = function(hypos){
                         mix94 <<- list('brown' = 30,
                                       'yellow' = 20,
                                       'red' = 20,
                                       'green' = 10,
                                       'orange' = 10,
                                       'tan' = 10)
                         
                         mix96 <<- list('blue' = 24,
                                        'green' = 20,
                                        'orange' = 16,
                                        'yellow' = 14,
                                        'red' = 13,
                                        'brown' = 13)
                         
                         hypoA <<- list('bag1' = mix94, 'bag2' = mix96)
                         hypoB <<- list('bag1' = mix96, 'bag2' = mix94)
                         
                         hypotheses <<- list('A' = hypoA, 'B' = hypoB)
                         
                         for (i in hypos) {
                           .self$Set(i,1)
                         }
                         .self$Normalize()
                       },
                       Likelihood = function(data, hypo){
                         bag <- data[1]
                         color <- data[2]
                         mix <- hypotheses[[hypo]][[bag]]
                         like <- mix[[color]]
                         return(like)
                       }
                     )
)#MandM class

#Examples after encapsulation 
hypos <- c("Bowl1", "Bowl2")
cookie <- Cookie$new(hypos)
cookie$Update("vanilla")
cookie$Print()

hypos <- c('A','B','C')
monty <- Monty$new(hypos) 
monty$Update('B')
monty$Print()

mandm <- MandM$new(c('A','B'))
mandm$Update(c('bag1', 'yellow'))
mandm$Update(c('bag2', 'green'))
mandm$Print()
