calc.ndc <-
function(gTree, sTree, assoc) {
  # Descr:  returns the number of deep coalescences for an entire tree
  # Deps:   calc.parse
  #         calchelpers.gtreeparse
  # I/p:    sTree
  #         gTree
  #         assoc
  # Note:   ndc = "number of deep coalescences"

  #dbgBool = get("P2C2M.flg.dbgBool", envir=p2c2m.globalVars)
  #if (dbgBool) {
  #  cat("\n",xtermStyle::style("DEBUGMODE> calc.ndc",fg="red"),sep="")
  #}

  handle = calc.parse(sTree, assoc)
  nodes = handle$nodes
  dmvD = handle$dmvD
  ndc = c()
  for (node in nodes) {
    tempData = calchelpers.gtreeparse(sTree, gTree, assoc, dmvD, node)
    ndc = c(ndc, tempData$lNd-1)
  }

  return(sum(ndc))
}
