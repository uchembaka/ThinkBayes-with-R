#  Classes and functions as defined in thinkbayes.py
#  Pmf: represents a probability mass function (map from values to probs).
#  DictWrapper: private parent class for Hist and Pmf.

Odds <- function(p){
  #' Compute odds for a given probability
  #' 
  #' Example: p = 0.75 means 75 for and 25 agains, or 3:1 odds in favor.
  #' 
  #' Note: when p = 1, Odds(1) = inifinity
  #' 
  #' @param p: float 0-1
  #' 
  #' @return float odds
  #' 
  
  if(p == 1) inivisible(Inf)
  
  inivisible(p/(p-1))
}#Odds()

Probability <- function(o){
  #'Computes the probability corresponding to given odds.
  #'
  #'@param o: float odds, strictly positive
  #'
  #'@return float probability
  #'
  #' @example o = 2 means 2:1 odds in favour, or 2/3 probability
  
  invisible(o(o+1))
}#Probability

Probability2 <- function(yes, no){
  #' Computes the probability corresponding to given odds.
  #' 
  #' @param yes, no: int or float odds in favor
  #' 
  #' @example yes=2, no=1 means 2:1 odds in favor, or 2/3 probability.
  
  invisible(yes/(yes+no))
}#Probability2


#' Interpolator class
#' 
#' Represents a mapping between sorted sequences; performs linear interp.
#' 
#' @field xs: sorted list
#' @field yx: sorted list
#' 
Interpolator <- setRefClass("Interpolator",
                            fields = list(xs = "vector", ys = "vector"),
                            methods = list(
                              initialize = function(xs, ys){
                                .self$xs <<- xs
                                .self$ys <<- ys
                              },#init
                              
                              Lookup = function(x){
                                "Looks up x and returns the corresponding value of y."
                                return(.self$Bisect(x, .self$xs, .self$ys))
                              },#Lookup()
                              
                              Reverse = function(y){
                                "Looks up y and returns the corresponding value of x."
                                
                                return(y, .self$ys, .self$xs)
                              },#Reverse()
                              
                              Bisect = function(x, xs, ys){
                                "Helper function."
                                
                                if(x <= xs[1]) return (.self$ys[1])
                                if(x >= xs[length(.self$xs)]) return (.self$ys[length(.self$ys)])
                                
                                i = findInterval(x, .self$xs)+1
                                
                                y = ys[i-1] + ((ys[i] - ys[i-1])/(xs[i] - xs[i-1]))*(x - xs[i-1])
                                
                                invisible(y)
                              }#Bisect
                            )
)#Interpolator class

#' DictWrapper Class
DictWrapper <- setRefClass("DictWrapper",
                           fields = list(ls = "list", name = "character", log = "logical"),
                           methods = list(
                             initialize = function(values = NULL, name = ''){
                               "Initializes the distribution.
                               
                               hypos: sequence of hypotheses"
                               
                               .self$name <<- name
                               .self$ls <<- list()
                               
                               #flag whether the distribution is under a log transform
                               .self$log <<- FALSE
                               
                               if(all(is.null(values))) return ()
                               
                               init_methods <- c(.self$InitPmf,
                                                 .self$InitMapping,
                                                 .self$InitSequence,
                                                 .self$InitFailure)
                               
                               for (method in init_methods){
                                 tryCatch({
                                   method(values)
                                   break
                                 },
                                 error = function(e) {""}
                                 )#trycatch()
                               }#for()
                               
                               if(length(.self) > 0){
                                 .self$Normalize()
                               }
                             },#initialize()
                             
                             InitSequence = function(values){
                              "Initializes with a sequence of equally-likely values.
                               
                               values: sequence of values"
                               
                               for(value in values){
                                 .self$Set(as.character(value), 1)
                               }
                             },#InitSequence()
                             
                             InitMapping = function(values){
                               "Initializes with a map from value to probability.
                               
                               values: map from value to probability"

                               if(is.null(names(values))) return(stop(""))
                               for(value in 1:length(values)){
                                 .self$Set(names(values[value]), value)
                               }
                               
                             },#InitMapping()
                             
                             InitPmf = function (values){
                               "Initializes with a Pmf.
                                
                               values: Pmf object"
                               for(name in names(values$Items())){
                                 .self$Set(name, values$Items()[[name]])
                               }
                             },#InitPmf()
                             
                             InitFailure = function(values){
                               "Raises an error."
                               stop("None of the initialization methods worked.")
                             },#InitFailure()
                             
                             len = function(){
                               return(length(.self$ls))
                             },
                             
                             contains = function(value){
                               return (.self$ls[as.character(value)])
                             },
                             
                             Copy = function(name = NA){
                               new <- .self$copy()
                               if(is.na(name)){
                                 invisible (new)
                               }
                               new$name <- name
                               invisible(new)
                             },#Copy()
                             
                             Scale = function(factor){
                               "Multiplies the values by a factor.
                                
                               factor: what to multiply by
                                
                               Returns: new object"
                               
                               if(any(is.na(as.numeric(names(.self$ls))))){
                                 stop("You can only scale numeric values")
                               }
                               new <- .self$Copy()
                               new$ls <- list()
                               
                               for(name in names(.self$Items())){
                                 newValues <- as.numeric(name)*factor
                                 new$Set(as.character(newValues), .self$Items()[[name]])
                               }
                               invisible (new)
                             },#Scale()
                             
                             Log = function(m = NA){
                               "Log transforms the probabilities.
                               
                               Removes values with probability 0.
                               
                               Normalizes so that the largest logprob is 0."
                               
                               if(.self$log){
                                 stop("PmF/Hist already under a log transform")
                               }
                               .self$log <- TRUE
                               
                               if(is.na(m)) m <- .self$MaxLike()
                               
                               for(name in names(.self$Items())){
                                 if(.self$Items()[[name]] == 0) .self$Remove(name)
                                 else .self$Set(name, log(.self$Items()[[name]]/m))
                               }
                             },#log()
                             
                             Exp = function(m = NA){
                               "Exponentiates the probabilities.
                               
                               m: how much to shift the ps before exponentiating
                               
                               If m is None, normalizes so that the largest prob is 1."
                               
                               if(!.self$log){
                                 stop("Pmf/Hist  not under a log transform")
                               }
                               .self$log <- FALSE
                               
                               if(is.na(m)) m <- .self$MaxLike()
                               
                               items <- .self$Items()
                               for(i in 1:length(items)){
                                  .self$Set(names(items[i]), exp(items[[i]]-m))
                               }
                             },#Exp()
                             
                             GetList = function(){
                               invisible(.self$ls)
                             },#GetList()
                             
                             SetList = function(ls){
                               .self$ls <<- ls
                             }, #SetList()
                             
                             Values = function(){
                               "Gets an unsorted sequence of values.
                               
                               Note: one source of confusion is that the keys of this
                               dictionary are the values of the Hist/Pmf, and the
                               values of the dictionary are frequencies/probabilities."
                               
                               if(any(is.na(as.numeric(names(.self$ls))))) invisible(names(.self$ls))
                               else invisible(as.numeric(names(.self$ls)))
                             },#Values()
                             
                             Items = function(){
                               invisible (unlist(.self$ls))
                             }, #ITems()
                             
                             Render = function(){
                               "Generates a sequence of points suitable for plotting.
                               
                               Returns:
                                  vector of (sorted value sequence, freq/prob sequence)
                               "
                              invisible(cbind(as.numeric(names(.self$ls)), unlist(.self$ls, use.names = F)))
                             },#Render()
                             
                             Print = function(){
                               "Prints the values and freqs/probs in ascending order."
                               items = .self$Items()[order(unlist(.self$Items(), use.names = F))]
                               items.name <- c(names(items))
                               i <- 1
                               for(item in items){
                                 cat(sprintf("%*s",5,items.name[i]),": ", item,"\n")
                                 i = i+1
                               }
                             },#Print()
                             
                             Set = function(x, y=0){
                               "Sets the freq/prob associated with the value x.
                               
                               Args:
                                  x: number value
                                  y: number freq or prob"
                                  
                               .self$ls[as.character(x)] <- y
                             },#Set()
                             
                             Incr = function(x, term = 1){
                              "Increments the freq/prob associated with the value x.
                               
                               Args:
                                  x: number value
                                  term: how much to increment by"
                                
                               if (x %in% names(.self$ls)){
                                 .self$ls[as.character(x)] <- .self$ls[[as.character(x)]]+1
                               }else{
                                 .self$ls[as.character(x)] <- 1
                               }
                             },#Incr()
                             
                             Mult = function(x, factor){
                               "Scales the freq/prob associated with the value x.
                               
                               Args:
                                 x: number value
                               factor: how much to multiply by"
                               
                               if (x %in% names(.self$ls)){
                                 .self$ls[as.character(x)] <- .self$ls[[as.character(x)]] * factor
                               }else{
                                 .self$ls[as.character(x)] <- 0 * factor
                               }
                               
                             },#Mult()
                             
                             Remove = function(x){
                               "Removes a value.
                                
                               Throws an exception if the value is not there.
                               
                               Args:
                                 x: value to remove"
                              
                               if(as.character(x) %in% names(.self$ls)){
                                 .self$ls <- .self$ls[names(.self$ls) != as.character(x)]
                               }else{
                                 stop("no such value")
                               }
                             }, #Remove()
                             
                             Total = function(){
                               "Returns the total of the frequencies/probabilities in the map."
                               
                               total <- sum(unlist(.self$ls, use.names = F))
                               return(total)
                             }, #Total()
                             
                             MaxLike = function(){
                               "Returns the largest frequency/probability in the map."
                               
                               return(max(unlist(.self$ls, use.names = F)))
                             }#MaxLike()
                           )#methods()
)#DictWrapper Class


#' Histogram class
#' 
#' Represents a histogram, which is a map from values to frequencies.
#' frequencies are integer counters.
Hist <- setRefClass("Hist",
                    contains = "DictWrapper",
                    methods = list(
                      Freq = function(x){
                        "Gets the frequency associated with the value x.
                        
                        Args:
                            x: number value
                            
                        Returns:
                            frequency
                        "
                        
                        if (x %in% names(.self$ls)){
                          return(.self$ls[[as.character(x)]])
                        }else{
                          return(0)
                        }
                      },#Freq()
                      
                      Freqs = function(xs){
                        "Gets frequencies for a sequence of values."
                        return(sapply(xs, function(x) .self$Freq(x)))
                      },#Freqs()
                      
                      IsSubset = function(other){
                        "Checks whether the values in this histogram are a subset of
                        the values in the given histogram."
                        
                        for(name in names(.self$Items())){
                          val <- name
                          freq <- .self$Items()[[name]]
                          if (freq > other$Freq(val)) return(FALSE)
                        }
                        
                        return(TRUE)
                      },#IsSubset
                      
                      Subtract = function(other){
                        "Subtracts the values in the given histogram from this histogram."
                        for(name in names(.self$Items())){
                          val <- name
                          freq <- .self$Items()[[name]]
                          .self$Incr(val, -freq)
                        }
                      }#Subtract
                    )
)#Hist class


#' Represents a probability mass function.
#' 
#' probabilities are floating-point.
#' Pmfs are not necessarily normalized.
Pmf <- setRefClass("Pmf",
                   contains = "DictWrapper",
                   methods = list(
                     Prob = function(x, default=0){
                       "Gets the probability associated with the value x.
                        
                       Args:
                          x: number value
                          default: value to return if the key is not there
                       
                       Returns:
                          float probability"
                       
                       if (x %in% names(.self$ls)){
                         return(.self$ls[[as.character(x)]])
                       }else{
                         return(default)
                       }
                     },#Prob()
                     
                     MakeCdf = function(name = NA){
                       "Makes a Cdf."
                       
                       return (MakeCdfFromPmf(pmf = .self, name = name))
                     },#MakeCdf()
                     
                     Probs = function(xs){
                       " Gets probabilities for a sequence of values"
                       return(sapply(xs, function(x) .self$Prob(x)))
                     },#Probs()
                     
                     ProbGreater = function(x){
                       temp <- .self$ls[which(as.numeric(names(.self$ls)) > as.numeric(x))]
                       temp <- as.vector(unlist(temp, use.names = F))
                       return (sum(temp))
                     },#ProbGreater
                     
                     ProbLesser = function(x){
                       temp <- .self$ls[which(as.numeric(names(.self$ls)) < as.numeric(x))]
                       temp <- as.vector(unlist(temp, use.names = F))
                       return (sum(temp))
                     },#ProbLesser
                     
                     Normalize = function(fraction = 1.0){
                       
                       "Normalizes this PMF so the sum of all probs is fraction.
                        
                       Args:
                          fraction: what the total should be after normalization
                       
                       Returns: the total probability before normalizing
                       "
                       
                       if(.self$log == TRUE) stop("Pmf is under a log transform")
                       
                       
                       total = .self$Total()
                       if(total == 0.0){
                         warning("Normalize: total probability is zero")
                         invisible (total)
                       }
                       
                       factor <- fraction/total
                       for(x in 1:length(.self$ls)){
                         .self$ls[x] <- .self$ls[[x]] * factor
                       }
                       invisible(total)
                     },#Normalize()
                     
                     Random = function(){
                       "Chooses a random element from this PMF.
                        
                       Returns:
                         float value from the Pmf"
                      
                       if(.self$len() == 0) stop("Pmf contains no value")
                       
                       target <- runif(1)
                       total <- 0.0
                       for(x in 1:length(.self$ls)){
                         total <- total+.self$ls[[x]]
                         if(total >= target) return(as.numeric(names(.self$ls[x])))
                       }
                     },#Random()
                     
                     Mean = function(){
                       "Computes the mean of a PMF.
                       
                       Returns:
                         float mean"
                       
                       x <- as.numeric(names(.self$ls))
                       mu <- 0.0
                       for(i in 1:length(.self$ls)){
                         mu <- mu + .self$ls[[i]]* x[i]
                       }
                       return (mu)
                     },#Mean()
                     
                     Var = function(mu = NA){
                       "Computes the variance of a PMF.
                       
                       Args:
                         mu: the point around which the variance is computed;
                         if omitted, computes the mean
                       
                       Returns:
                        float variance"
                      
                       if(is.na(mu)) mu <- .self$Mean()
                       
                       x <- as.numeric(names(.self$ls))
                       var <- 0.0
                       for(i in 1:length(.self$ls)){
                         var <- var + .self$ls[[i]]*(x[i]-mu)^2 
                       }
                       return (var)
                     },#Var()
                     
                     MaximumLikelihood = function(){
                       "Returns the value with the highest probability.
                       
                       Returns: float probability"
                       
                       return(max(as.numeric(names(.self$ls[which(.self$Items() == max(.self$Items()))]))))
                     },#MaximumLikelihood
                     
                     CredibleInterval = function(percentage = 90){
                       "Computes the central credible interval.
                       
                       If percentage=90, computes the 90% CI.
                       
                       Args:
                         percentage: float between 0 and 100
                       
                       Returns:
                         sequence of two floats, low and high"
                      
                       cdf <- .selfMakeCdf()
                       return (cdf.CredibleInterval(percentage))
                     },#CredibleInterval()
                     
                     add = function(other){
                       "Computes the Pmf of the sum of values drawn from self and other.
                       
                       other: another Pmf
                       
                       returns: new Pmf"
                       tryCatch(return(.self$AddPmf(other)), 
                                finally = return(.self$AddConstant(other)))
                     },#add()
                     
                     AddPmf = function(other){
                       "Computes the Pmf of the sum of values drawn from self and other.
                       
                       other: another Pmf
                       
                       returns: new Pmf
                       "
                       pmf <- Pmf$new()
                       for(i in 1:length(.self$Items())){
                         v1 <- as.numeric(names(.self$Items()[i]))
                         p1 <- .self$Items()[[i]]
                         for(j in 1:length(other$Items())){
                           v2 <- as.numeric(names(other$Items()[j]))
                           p2 <- other$Items()[[j]]
                           print(v1+v2)
                           pmf$Incr(v1+v2, p1*p2)
                         }
                       }
                       invisible(pmf)
                     },#AddPmf()
                     
                     AddConstant = function(other){
                       "Computes the Pmf of the sum a constant and  values from self.
                       
                       other: a number
                       
                       returns: new Pmf
                       "
                       pmf <- Pmf$new()
                       for(i in 1:length(.self$Items())){
                         v1 <- as.numeric(names(.self$Items()[i]))
                         p1 <- .self$Items()[[i]]
                         pmf$Set(v1+other, p1)
                       }
                       
                       invisible(pmf)
                     },#AddConstant()
                     
                     sub = function(other){
                       "Computes the Pmf of the diff of values drawn from self and other.
                       
                       other: another Pmf
                       
                       returns: new Pmf
                       "
                       pmf <- Pmf$new()
                       for(i in 1:length(.self$Items())){
                         v1 <- as.numeric(names(.self$Items()[i]))
                         p1 <- .self$Items()[[i]]
                         for(j in 1:length(other$Items())){
                           v2 <- as.numeric(names(other$Items()[j]))
                           p2 <- other$Items()[[j]]
                           pmf$Incr(v1-v2, p1*p2)
                         }
                       }
                       invisible(pmf)
                     },#sub()
                     
                     Max = function(k){
                       "Computes the CDF of the maximum of k selections from this dist.
                       
                       k: int
                       
                       returns: new Cdf
                       "
                       cdf <- .self$MakeCdf()
                       cdf$ps <- sapply(cdf$ps, function(x) x^k)
                       invisible(cdf)
                     }#max()
                   )
)#Pmf Class


#' Represents a joint distribution.
#' 

Joint <- setRefClass("Joint",
                     contains = "Pmf",
                     methods = list(
                       
                       unpack = function(v){
                         "helper function"
                         invisible (sapply(unlist(strsplit(v,","), use.names = F), function(x) as.numeric(x)))
                       },#unpack
                       
                       pack = function(v){
                         "helper function"
                         invisible(paste(v, collapse = ","))
                       },#pack
                       
                       Marginal = function(i , name=''){
                         "Gets the marginal distribution of the indicated variable.
                         
                         i: index of the variable we want
                         
                         Returns: Pmf
                         "
                         pmf <- Pmf(name = name)
                         for(val in names(.self$Items())){
                           pmf$Incr(.self$unpack(val)[i], .self$Items()[[val]])
                         }
                         
                         invisible(pmf)
                       },#Marginal()
                       
                       Conditional = function(i, j, val, name=''){
                         "Gets the conditional distribution of the indicated variable.
                         
                         Distribution of vs[i], conditioned on vs[j] = val.
                         
                         i: index of the variable we want
                         j: which variable is conditioned on
                         val: the value the jth variable has to have
                         
                         Returns: Pmf
                         "
                         pmf <- Pmf$new(name = name)
                         for(vs in names(.self$Items())){
                           if(.self$unpack(vs)[j] != val) next
                           pmf$Incr(.self$unpack(vs)[i], .self$Items()[[vs]])
                         }
                         pmf$Normalize()
                         
                         invisible(pmf)
                       },#Conditional()
                       
                       MaxLikeInterval = function(percentage = 90){
                         "Returns the maximum-likelihood credible interval.
                         
                         If percentage=90, computes a 90% CI containing the values
                         with the highest likelihoods.
                         
                         percentage: float between 0 and 100
                         
                         Returns: list of values from the suite
                         "
                         
                         interval <- c()
                         total <- 0
                         
                         t <- sort(.self$Items(), decreasing = T)
                         
                         for(val in names(t)){
                           interval <- append(interval, val)
                           total <- total+t[[val]]
                           if(total >= percentage/100) break
                         }
                         return(interval)
                         
                       }#MaxLikeInterval()
                     )
)#Joint pmf class


#' Make Joint Distribution
#' 
#' Joint distribution of values from pmf1 and pmf2
#' 
#' @param pmf1: Pmf Object
#' @param pmf2: Pmf Object
#' 
#' @return Joint pmf of values pairs
#' 
MakeJoint = function(pmf1, pmf2){
  joint <- Joint()
  for(v1 in names(pmf1$Items())){
    for(v2 in names(pmf2$Items())){
      joint$Set(joint$pack(c(v1,v2)), pmf1$Items()[[v1]]*pmf2$Items()[[v2]])
    }
  }
  invisible(joint)
}#MakeJoint

#' MakeHistFromList
#' 
#' Makes a histogram from an unsorted sequence of values
#' 
#' @param t: sequence of numbers
#' @param name: name of histogram
#' 
#' @return Hist object
#' 
MakeHistFromList <- function(t, name=""){
  hist <- Hist(name = name)
  for(x in t){
    hist$Incr(x)
  }
  
  invisible(hist)
}#MakeHistFromList()

#' MakeHistFromList
#' 
#' Makes a histogram from a map from values to frequencies.
#' 
#' @param d: list that maps value to frequencies
#' @param name: name of this histogram
#' 
#' @return Hist object
#' 
MakeHistFromDict = function(d, name=""){
  invisible(Hist(d, name))
}#MakeHistFromDict

#' MakePmfFromList
#' 
#' Makes a PMF from an unsorted sequence of values.
#' 
#' @param t: sequence of number
#' @param name: name of this PMF
#' 
#' @return pmf object
#' 
MakePmfFromList = function(t, name=""){
  hist <- MakeHistFromList(t)
  d <- hist$GetList()
  pmf <- Pmf(d, name)
  pmf$Normalize()
  invisible(pmf)
}#MakePmfFromList()

#' MakePmfFromDict
#' 
#' Makes a PMF from a map from values to probabilities
#' 
#' @param d: list thta maps value to probabilities
#' @param name: name of this pmf
#' 
#' @return pmf object
#' 
MakePmfFromDict =function(d, name=""){
  pmf <- Pmf(d, name)
  pmf$Normalize()
  invisible(pmf)
}#MakePmfFromDict()

#' Represents a cumulative distribution function.
#' 
#'  @field xs: sequence of values
#'  @field ps: sequence of probabilities
#'  @field name: string used as a graph label.

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
                     
                     Copy = function(name = NA){
                       "Returns a copy of this Cdf.
                       
                       Args:
                          name: string name for the new Cdf
                       "
                       new <- .self$copy()
                       if(is.na(name)){
                         invisible (new)
                       }
                       new$name <- name
                       invisible(new)
                     },#Copy()
                     
                     MakePmf = function(name=NA){
                       "Makes a Pmf"
                       
                       invisible(MakePmfFromCdf(name = name))
                     },#MakePmf
                     
                     Values = function(){
                       invisible(.self$xs)
                     },#Values
                     
                     Items = function(){
                       "Returns a sorted sequence of (value, probability) pairs."
                       temp <- .self$ps
                       names(temp) <- .self$xs
                       invisible(temp)
                     },#Items()
                     
                     Append = function(x, p){
                       "Add an (x, p) pair to the end of this CDF.
                       
                       Note: this us normally used to build a CDF from scratch, not
                       to modify existing CDFs.  It is up to the caller to make sure
                       that the result is a legal CDF.
                       "
                       .self$xs <<- append(.self$xs, x)
                       .self$ps <<- append(.self$ps, p)
                     },#Append()
                     
                     Shift = function(term){
                       "Adds a term to the xs.
                       
                       term: how much to add
                       "
                       new <- .self$Copy()
                       new$xs <- sapply(.self$xs, function(x) x+term)
                       invisible(new)
                     },#Shift()
                     
                     Scale = function(factor){
                       "Multiplies the xs by a factor
                       
                       factor: what to multiply by
                       "
                       new <- .self$Copy()
                       new$xs <- sapply(.self$xs, function(x) x*factor)
                       invisible(new)
                     },#Scale()
                     
                     Prob = function(x){
                       "Returns CDF(x), the probability that corresponds to value x
                       
                       Args:
                          x: number
                          
                      Returns:
                          float probability
                       "
                       
                       if(x < .self$xs[1]) return(0.0)
                       
                       index <- findInterval(x, .self$xs)+1
                       return(.self$ps[index]-1)
                     },
                     
                     Value = function(p){
                       "Returns InverseCDF(p), the value that corresponds to probability p.
                       
                       Args:
                          p: number in the range [0, 1]
                          
                        Returns:
                          number value
                       "
                       if(p < 0 || p > 1){
                         stop("Probability p must be in range[0,1]")
                       }
                       if(p == 0) return(.self$xs[1])
                       if(p == 1) return(.self$xs[length(.self$xs)])
                       index <- findInterval(p, .self$ps)+1
                       if((index-1) == 0) return(.self$xs[1])
                       else{
                         if(p == .self$ps[index-1]) return (.self$xs[index-1])
                         else return (.self$xs[index])
                       }
                     },#Value()
                     
                     Percentile = function(p){
                       "Returns the value that corresponds to percentile p.
                       
                       Args:
                          p: number in the range [0, 100]
                          
                        Returns:
                          number value
                       "
                       return (.self$Value(p/100))
                     },#Percentile()
                     
                     Random = function(){
                       "Chooses a random value from this distribution."
                       return(.self$Value(runif(1)))
                     },#Random()
                     
                     Sample = function(n){
                       "Generates a random sample from this distribution.
                       
                       Args:
                          n: int length of the sample
                       "
                       sample <- vector("numeric", length = n)
                       return(sample, function(x) .self$Random())
                     },#Sample()
                     
                     Mean = function(){
                       "Computes the mean of a CDF
                       
                       Returns:
                          mean
                       "
                       old_p = 0
                       total = 0
                       
                     }
                   )
)#CDF class


#' MakeCdfFromPmf
#' 
#' Makes a CDF from a Pmf object.
#' 
#' @param pmf: Pmf.Pmf object
#' name: string name for the data.
#' 
#' @return Cdf object
MakeCdfFromPmf = function(pmf, name = NA){
  return (MakeCdfFromItems(pmf$Items()))
}

#' MakeCdfFromItems
#' 
#' Makes a cdf from an unsorted sequence of (value, frequency) pairs.
#' 
#' @param items: unsorted sequence of (value, frequency) pairs
#' name: string name for this CDF
#' 
#' @return  cdf: list of (value, fraction) pairs
MakeCdfFromItems = function(items, name=""){
  runsum <- 0;
  xs <- vector()
  cs <- vector()
  
  #items <- sort(items)
  
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


#'Represents a suite of hypotheses and their probabilities.
Suite <- setRefClass("Suite",
                     contains = "Pmf",
                     methods = list(
                       Update = function(data){
                         "Updates each hypothesis based on the data.
                         data: any representation of the data
                         returns: the normalizing constant
                         "
                         
                         for(hypo in names(.self$ls)){
                           like <- .self$Likelihood(data, hypo)
                           .self$Mult(hypo, like)
                         }
                         invisible(.self$Normalize())
                       },#Update()
                       
                       LogUpdate = function(data){
                         "Updates a suite of hypotheses based on new data.
                         
                         Modifies the suite directly; if you want to keep the original, make a copy.
                         
                         Note: unlike Update, LogUpdate does not normalize.
                         
                         Args:
                            data: any representation of the data
                         "
                         
                         for(hypo in names(x)){
                           like <- .self$LogLikelihood(data, hypo)
                           .self$Incr(hypo, like)
                         }
                       },#LogUpdate()
                       
                       UpdateSet = function(dataset){
                         "Updates each hypothesis based on the dataset.
                         
                         This is more efficient than calling Update repeatedly because
                         it waits until the end to Normalize.
                         
                         Modifies the suite directly; if you want to keep the original, make a copy.
                         
                         dataset: a sequence of data
                         
                         returns: the normalizing constant
                         "
                         
                         for (data in dataset){
                           for(hypo in names(x)){
                             like <- .self$Likelihood(data, hypo)
                             .self$Mult(hypo, like)
                           }
                         }
                         invisible(.self$Normalize())
                       },#UpdateSet()
                       
                       LogUpdateSet = function(dataset){
                         "Updates each hypothesis based on the dataset.
                         
                         Modifies the suite directly; if you want to keep the original, make a copy.
                         
                         dataset: a sequence of data
                         
                         returns: None
                         "
                         for(data in dataset){
                           .self$LogUpdate(data)
                         }
                       },#LogUpdateSet()
                       
                       Likelihood = function(data, hypo){
                         "Computes the likelihood of the data under the hypothesis.
                         
                         hypo: some representation of the hypothesis
                         data: some representation of the data
                         "
                         stop("Unimplemented Method")
                       },#Likelihood()
                       
                       LogLikelihood = function(data, hypo){
                         "Computes the log likelihood of the data under the hypothesis.
                         
                         hypo: some representation of the hypothesis
                         data: some representation of the data
                         "
                         stop("Unimplemented Method")
                       },#LogLikelihood(),
                       
                       Print = function(){
                         items.name <- c(names(.self$Items()))
                         i <- 1
                         for(item in .self$Items()){
                           cat(items.name[i],": ", item,"\n")
                           i = i+1
                         }
                       },#Print()
                       
                       MakeOdds = function(){
                         "Transforms from probabilities to odds.
                         
                         Values with prob=0 are removed.
                         "
                         for(i in 1:length(.self$Items())){
                           hypo <- .self$Items()[i]
                           prob <- .self$Items()[[i]]
                           if (prob == 0){
                             .self$Remove(hypo)
                           }else{
                             .self$Set(hypo, Odds(prob))
                           }
                         }
                       },#MakeOdds()
                       
                       MakeProbs = function(){
                         "Transforms from odds to probabilities."
                         for(i in 1:length(.self$Items())){
                           hypo <- .self$Items()[i]
                           odds <- .self$Items()[[i]]
                           .self$Set(hypo, Probability(odds))
                         }
                       }#MakeProbs()
                     )
)#Suite class

#' Percentile
#' 
#' Computes a percentile of a given pmf
#' 
#' @param pmf: a pmf object
#' @param percentage: value 0 - 100
#' 

Percentile <- function(pmf, percentage){
  p <- percentage/100
  total <- 0
  for(val in names(pmf$Items())){
    total <- total + pmf$Items()[[val]]
    if(total >= p) return(val)
  }
}#Percentile

#' CredibleInterval
#' 
#' Computes a credible interval for a given distribution.
#' If percentage=90, computes the 90% CI.
#' 
#' @param pmf: Pmf object representing a posterior distribution
#' @param percentage: float between 0 and 100
#' 
#' @return sequence of two floats, low and high
#' 
CredibleInterval <- function(pmf, percentage = 90){
  cdf <- pmf$MakeCdf()
  prob <- (1-percentage/100)/2
  interval <- c(cdf$Value(prob), cdf$Value(1-prob))
  return(interval)
}#CredibleInterval


#' Beta Distribution
#' 
#' Represents a Beta distribution
#' See http://en.wikipedia.org/wiki/Beta_distribution
#' 
Beta <- setRefClass("Beta",
                    fields = list(alpha = "numeric", beta = "numeric", name = "character"),
                    methods = list(
                      initialize = function(alpha =1, beta = 1, name = ""){
                        .self$alpha <<- alpha
                        .self$beta <<- beta
                        .self$name <<- name
                      },#initialize
                      
                      Update = function(data){
                        "Updates a Beta distribution
                        
                        data: numeric vector (heads, tails)
                        "
                        heads <- data[1]
                        tails <- data[2]
                        .self$alpha <- .self$alpha + heads
                        .self$beta <- .self$beta + tails
                      },#Update()
                      
                      Mean = function(){
                        "Computes the mean of this distribution."
                        return(.self$alpha/(.self$alpha+.self$beta))
                      },#Mean()
                      
                      Random = function(){
                        "Generates a random variate from this distribution."
                        return(rbeta(1,.self$alpha, .self$beta))
                      },#Random()
                      
                      EvalPdf = function(x){
                        "Evaluates the PDF at x."
                        x^(.self$alpha-1)*(1-x)^(.self$beta-1)
                      },#EvalPdf()
                      
                      MakePmf = function(steps = 100, name=""){
                        "Returns a Pmf of this distribution.
                        
                        Note: Normally, we just evaluate the PDF at a sequence
                        of points and treat the probability density as a probability mass
                        
                        But if alpha or beta is less than one, we have to be
                        more careful because the PDF goes to infinity at x=0
                        and x=1.  In that case we evaluate the CDF and compute
                        differences.
                        "
                        if(.self$alpha < 1 | .self$beta < 1){
                          cdf <- .self$MakeCdf()
                          pmf <- cdf$MakePmf()
                          invisible(pmf)
                        }
                        
                        xs <- sapply(0:steps, function(x) x/steps)
                        probs <- sapply(xs, function(x) .self$EvalPdf(x))
                        d <- probs
                        names(d) <- xs
                        pmf <- MakePmfFromDict(d, name)
                        invisible(pmf)
                      },#MakePmf()
                      
                      MakeCdf = function(steps = 100){
                        "Returns the CDF of this distribution."
                        xs <- sapply(0:steps, function(x) x/steps)
                        ps <- sapply(xs, function(x) pbeta(x, .self$alpha, .self$beta))
                        cdf <- Cdf(xs, ps)
                        invisible(cdf)
                      }#Make
                    )
)#Beta class
