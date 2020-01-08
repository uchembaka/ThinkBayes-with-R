#' Class and functions as defined in thinkbayes.py
#' Pmf: represents a probability mass function (map from values to probs).
#' DictWrapper: private parent class for Hist and Pmf.


DictWrapper <- setRefClass("DictWrapper",
                           fields = list(ls = "list", name = "character", log = "logical"),
                           methods = list(
                             initialize = function(values = NA, name = ''){
                               #' Initializes the distribution.
                               #' 
                               #' hypos: sequence of hypotheses
                               
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
                               #'Initializes with a sequence of equally-likely values.
                               #'
                               #'values: sequence of values
                               
                               for(value in values){
                                 .self$Set(as.character(value), 1)
                               }
                             },#InitSequence()
                             
                             InitMapping = function(values){
                               #'Initializes with a map from value to probability.
                               #'
                               #'values: map from value to probability

                               if(is.null(names(values))) return(stop(""))
                               for(value in 1:length(values)){
                                 .self$Set(names(values[value]), value)
                               }
                               
                             },#InitMapping()
                             
                             InitPmf = function (values){
                               #' Initializes with a Pmf.
                               #' 
                               #'values: Pmf object
                               
                               items <- values$Item()
                               for(value in 1:length(items)){
                                 .self$Set(names(items[value]), value)
                               }
                             },#InitPmf()
                             
                             InitFailure = function(values){
                               #' Raises an error.
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
                                 return (new)
                               }
                               new$name <- name
                               return(new)
                             },#Copy()
                             
                             Scale = function(factor){
                               #' Multiplies the values by a factor.
                               #' 
                               #' factor: what to multiply by
                               #' 
                               #' Returns: new object
                               
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
                               return (new)
                             },#Scale()
                             
                             Log = function(m = NA){
                               #'Log transforms the probabilities.
                               #'
                               #'Removes values with probability 0.
                               #'
                               #'Normalizes so that the largest logprob is 0.
                               
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
                               #' Exponentiates the probabilities.
                               #' 
                               #'m: how much to shift the ps before exponentiating
                               #'
                               #'If m is None, normalizes so that the largest prob is 1.
                               
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
                               return(.self$ls)
                             },#GetList()
                             
                             SetList = function(ls){
                               .self$ls <<- ls
                             }, #SetList()
                             
                             Values = function(){
                               #' Gets an unsorted sequence of values.
                               #'
                               #' Note: one source of confusion is that the keys of this
                               #' dictionary are the values of the Hist/Pmf, and the
                               #' values of the dictionary are frequencies/probabilities.
                               
                               if(any(is.na(as.numeric(names(.self$ls))))) return(names(.self$ls))
                               else return(as.numeric(names(.self$ls)))
                             },#Values()
                             
                             Items = function(){
                               return (.self$ls)
                             }, #ITems()
                             
                             Render = function(){
                               #' Generates a sequence of points suitable for plotting.
                               
                               #' Returns:
                                  #'vector of (sorted value sequence, freq/prob sequence)
                              pl <- .self$ls[order(unlist(.self$ls, use.names = F))]
                              return(cbind(as.numeric(names(pl)), unlist(pl, use.names = F)))
                             },#Render()
                             
                             Print = function(){
                               #' Prints the values and freqs/probs in ascending order.
                               items = .self$Items()[order(unlist(.self$Items(), use.names = F))]
                               items.name <- c(names(items))
                               i <- 1
                               for(item in items){
                                 cat(items.name[i],": ", item,"\n")
                                 i = i+1
                               }
                             },#Print()
                             
                             Set = function(x, y=0){
                               #' Sets the freq/prob associated with the value x.
                               #'
                               #' Args:
                                  #' x: number value
                                  #' y: number freq or prob
                                  
                               .self$ls[as.character(x)] <- y
                             },#Set()
                             
                             Incr = function(x, term = 1){
                               #' Increments the freq/prob associated with the value x.
                               #'
                               #' Args:
                                  #' x: number value
                                  #' term: how much to increment by
                                
                               if (x %in% names(.self$ls)){
                                 .self$ls[as.character(x)] <- .self$ls[[as.character(x)]]+1
                               }else{
                                 .self$ls[as.character(x)] <- 1
                               }
                             },#Incr()
                             
                             Mult = function(x, factor){
                               #' Scales the freq/prob associated with the value x.
                               
                               #' Args:
                                 #" x: number value
                               #' factor: how much to multiply by
                               
                               .self$ls[as.character(x)] <- .self$ls[[as.character(x)]] * factor
                             },#Mult()
                             
                             Remove = function(x){
                               #' Removes a value.
                               #' 
                               #' Throws an exception if the value is not there.
                               
                               #' Args:
                                 #' x: value to remove
                              
                               if(as.character(x) %in% names(.self$ls)){
                                 .self$ls <- .self$ls[names(.self$ls) != as.character(x)]
                               }else{
                                 stop("no such value")
                               }
                             }, #Remove()
                             
                             Total = function(){
                               #' Returns the total of the frequencies/probabilities in the map.
                               
                               total <- sum(unlist(.self$ls, use.names = F))
                               return(total)
                             }, #Total()
                             
                             MaxLike = function(){
                               #' Returns the largest frequency/probability in the map.
                               
                               return(max(unlist(.self$ls, use.names = F)))
                             }#MaxLike()
                           )#methods()
)#DictWrapper Class

