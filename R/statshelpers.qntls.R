statshelpers.qntls <-
function (ind, alpha, tail) {
  # Descr:  generates a distribution of quantiles;
  #         quantile levels are hardcoded in variable "qntlLevels"
  # Deps:   (none)
  #         ind = a data frame
  #         alpha
  #         tail

  dbgBool = get("P2C2M.flg.dbgBool", envir=p2c2m.globalVars)
  if (dbgBool) {
    cat("\n", xtermStyle::style("DEBUGMODE> statshelpers.qntls",
        fg="red"), sep="")
  }

  #################################
  # Set quantile levels via alpha #
  #################################
  if (tail=="1l" | tail=="1r") {
    qntlLevels = c(alpha, 1-alpha)
  }
  # Adjusting alpha value for two-tailed test
  if (tail=="2") {
    qntlLevels = c((alpha/2), 1-(alpha/2))
  }

  qntls = t(apply(ind, MARGIN=2, quantile, qntlLevels, na.rm=T))
  #qntls = rbind(qntls, quantile(acrossGenes, qntlLevels, na.rm=T))

  return(qntls)
}
