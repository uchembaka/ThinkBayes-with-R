source(".....thinkbayes.r")
#Chapter 2 exercise using source code from thinkbayes.r


#6 sided die example
die <- Pmf$new()
for(i in 1:6){
  die$Set(i, 1/6)
}

die$Print()

word.list <- c("Hello", "Hello", "world")
word.count <- Pmf()
for(word in word.list){
  word.count$Incr(word, 1)
}

word.count$Print()

word.count$Normalize()

word.count$Print()#after normalization

#Cookie Problem
cookie <- Pmf()
cookie$Set("Bowl2", 0.5)
cookie$Set("Bowl1", 0.5)
cookie$Print()
cookie$Mult("Bowl2", 0.5)
cookie$Mult("Bowl1", 0.75)
cookie$Print()
cookie$Normalize()
cookie$Print()

#The Bayesian framework

Cookie <- setRefClass("Cookie",
                      fields = list(mixes = "list"),
                      contains = "Suite",
                      methods = list(
                        initialize = function(hypos){
                          mixes <<- list("Bowl1" = list("vanilla" = 0.75, "chocolate" = 0.25),
                                         "Bowl2" = list("vanilla" = 0.5, "chocolate" = 0.5))
                          callSuper(hypos)
                        },
                        Likelihood = function(data, hypo){
                          mix <- .self$mixes[[hypo]]
                          like <- mix[[data]]
                          return (like)
                        })
)#Cookie Class

Monty <- setRefClass("Monty",
                     contains = "Suite",
                     methods = list(
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
                         
                         callSuper(hypos)
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

