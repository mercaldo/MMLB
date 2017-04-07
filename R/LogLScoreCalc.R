LogLScoreCalc <-
function( params, subjectData, CondLike = FALSE, ParamLengths,
                         EmpiricalCheeseCalc = FALSE, Q, W, Z, AdaptiveQuad=FALSE){#, offset){ #, SampProbs=sampprobs){
  
  betaM <- params[1:ParamLengths[1]]
  gamma <- params[ (ParamLengths[1]+1): (ParamLengths[1]+ParamLengths[2]) ]
  sigma <- exp(params[(ParamLengths[1]+ParamLengths[2]+1): (ParamLengths[1]+ParamLengths[2]+ParamLengths[3])])
  
  .Call("LogLScoreCalc_CALL", betaM, gamma, sigma, subjectData, CondLike, #SampProbs,
        EmpiricalCheeseCalc, Q, W, Z,  AdaptiveQuad) #offset,
}
