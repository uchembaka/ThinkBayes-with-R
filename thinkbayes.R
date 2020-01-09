#  Class and functions as defined in thinkbayes.py
#  Pmf: represents a probability mass function (map from values to probs).
#  DictWrapper: private parent class for Hist and Pmf.


#' DictWrapper Class
DictWrapper <- setRefClass("DictWrapper",
                           fields = list(ls = "list", name = "character", log = "logical"),
                           methods = list(
                             initialize = function(values = NA, name = ''){
                               "Initializes the distribution.
                               
                               hypos: sequence of hypotheses"
                               
                               .self$name <<- name
                               .self$ls <<- list()
                               
                               #flag whether the distribution is under a log transform
                               .self$log <<- FALSE
                               
                               if(all(is.na(values))) return ()
                               
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
                               
                               items <- values$Item()
                               for(value in 1:length(items)){
                                 .self$Set(names(items[value]), value)
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
                               return (ls[as.character(value)])
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
                               
                               items <- .self$Items()
                               for(value in 1:length(items)){
                                 newValues <- as.numeric(names(items[value]))*factor
                                 new$Set(as.character(newValues), items[[value]])
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
                               
                               items <- .self$Items()
                               for(i in 1:length(items)){
                                 if(items[[i]] == 0) .self$Remove(names(item[i]))
                                 else .self$Set(names(items[i]), log(items[[i]]/m))
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
                               items <- c()
                               items <- unlist(.self$ls, use.names = F)
                               names(items) <- names(.self$ls)
                               invisible (items)
                             }, #ITems()
                             
                             Render = function(){
                               "Generates a sequence of points suitable for plotting.
                               
                               Returns:
                                  vector of (sorted value sequence, freq/prob sequence)"
                                  
                              pl <- .self$ls[order(unlist(.self$ls, use.names = F))]
                              invisible(cbind(as.numeric(names(pl)), unlist(pl, use.names = F)))
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

#' Represents a probability mass function.
#' 
#' Values can be any hashable type; probabilities are floating-point.
#' Pmfs are not necessarily normalized.
Pmf <- setRefClass("PmF",
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
                         invisible(.self$ls[[as.character(x)]])
                       }else{
                         invisible(default)
                       }
                     },#Prob()
                     
                     MakeCdf = function(name = NA){
                       "Makes a Cdf."
                       
                       return (MakeCdfFromPmf(pmf = .self, name = name))
                     },#MakeCdf()
                     
                     Probs = function(xs){
                       " Gets probabilities for a sequence of values"
                       
                       probs = c()
                       for(x in xs){
                         probs <- append(probs, .self$Prob(x))
                       }
                       invisible (probs)
                     },#Probs()
                     
                     ProbGreater = function(x){
                       temp <- .self$ls[which(as.numeric(names(.self$ls)) > as.numeric(x))]
                       temp <- as.vector(unlist(temp, use.names = F))
                       invisible (sum(temp))
                     },#ProbGreater
                     
                     ProbLesser = function(x){
                       temp <- .self$ls[which(as.numeric(names(.self$ls)) < as.numeric(x))]
                       temp <- as.vector(unlist(temp, use.names = F))
                       invisible (sum(temp))
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
                         if(total >= target) invisible (as.numeric(names(.self$ls[x])))
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
                       invisible (mu)
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
                       invisible (var)
                     },#Var()
                     
                     MaximumLikelihood = function(){
                       "Returns the value with the highest probability.
                       
                       Returns: float probability"
                       
                       invisible(max(as.numeric(names(.self$ls[which(.self$Items() == max(.self$Items()))]))))
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


#' Represents a cumulative distribution function.
#' 
#'  @field xs: sequence of values
#'  ps: sequence of probabilities
#'  name: string used as a graph label.

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
                       index <- findInterval(p, .self$ps)
                       if(index == 0){
                         tmp <- .self$ps
                         tmp <- append(tmp, p)
                         tmp <- sort(tmp)
                         index <- findInterval(p, tmp)
                       }#replacing bisect 
                       
                       if((index-1) == 0) return(.self$xs[1])
                       else{
                         if(p == .self$ps[index-1]) return (.self$xs[index-1])
                         else return (.self$xs[index])
                       }
                     },
                     Percentile = function(p){
                       "Returns the value that corresponds to percentile p.
                       
                       Args:
                          p: number in the range [0, 100]
                          
                        Returns:
                          number value
                       "
                       return (.self$Value(p/100))
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



