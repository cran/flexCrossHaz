print.flexCrossHaz <-
function(x, digits = max(3, getOption("digits") - 3),...){
  cat("***Crossing Hazards in Cox model***\n")
  cat("Variable: ", names(x$pterms), "\n")
  cat("Estimate:", round(x$est,2), "( SE = ", round(x$s.e.,3),")\n")
    invisible(x)
}
  
  

