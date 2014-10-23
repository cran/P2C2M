calc.parse <-
function(sTree, assoc) {
  # Descr:  parses species tree nodes for metric calculation
  # Deps:   calchelpers.dmvparse
  # I/p:    sTree
  #         assoc

#  dbgBool = get("P2C2M.flg.dbgBool", envir=p2c2m.globalVars)
#  if (dbgBool) {
#    cat("\n",xtermStyle::style("DEBUGMODE> calc.parse",fg="red"),sep="")
#  }

  spNames = sTree$tip.label
  nSp = length(spNames)
  nBr = (2*nSp)-1
  tiplist = list()
  for(i in 1:nSp) {
    tiplist[[spNames[i]]]=assoc[which(assoc[,1]==spNames[i]),2]
  }
  dmvD = calchelpers.dmvparse(sTree, nSp)
  # returns demographic info of the species tree
  nTips = lapply(tiplist, length)
  # evaluate number of tips per species tree
  if(any(nTips==1)) {
    tmp = dmvD[-which(nTips==1),]
  } else {
    tmp = dmvD
  }
  nodes = tmp[,"node"]
  names(nodes) = NULL

  outd = list()
  outd$nodes = nodes
  outd$dmvD = dmvD

  return(outd)
}
