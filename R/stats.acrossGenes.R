stats.acrossGenes <-
function (diff, alpha, tail) {
  # Descr:    calculates statistics across genes
  # Deps:     statshelpers.qntls
  #           statshelpers.cv
  # I/p:      diff
  #           alpha
  #           tail

  dbgBool = get("P2C2M.flg.dbgBool", envir=p2c2m.globalVars)
  if (dbgBool) {
    cat("\n", xtermStyle::style("DEBUGMODE> stats.acrossGenes", fg="red"), sep="")
  }

  # Perform calculations per row (i.e per MCMC generation)
  acrossG_Sum = rowSums(diff)
  acrossG_Mean = rowMeans(diff)
  acrossG_Median = rowMedians(diff)
  acrossG_Mode = rowModes(diff)
  acrossG_CV = statshelpers.cv(diff)

  # "acrossG" is a matrix with four columns hereafter
  acrossG = cbind(acrossG_Sum, 
                  acrossG_Mean, 
                  acrossG_Median,
                  acrossG_Mode,
                  acrossG_CV)

  qntls = statshelpers.qntls(acrossG, alpha, tail)
  sigSgns = statshelpers.sigsgn(qntls, tail)

  subCol1 = c(apply(acrossG, MARGIN=2, mean))
  subCol2 = c(apply(acrossG, MARGIN=2, sd))
  subCol3 = c(sigSgns)

  outD = paste(round(subCol1, 2), 
               paste("(", "\u00B1", round(subCol2, 2), ")", sep=""), 
               subCol3)

  return(outD)
}
