calc.coal_reid <-
function(gTree, sTree, assoc, ploidy) {
  # Descr:    calculates the probability of a whole gene tree
  #           by calculating the probability of subtrees defined 
  #           by their MRCA (i.e. their node); done for all 
  #           nodes in a tree, values are summed up thereafter
  # Deps:     calc.parse
  #           calchelpers.brprob
  # I/p:      sTree
  #           gTree
  #           assoc
  #           ploidy
  # Note:     gtp = "gene tree probability" (a form of coalescent likelihood)

  #dbgBool = get("P2C2M.flg.dbgBool", envir=p2c2m.globalVars)
  #if (dbgBool) {
  #  cat("\n",xtermStyle::style("DEBUGMODE> calc.coal_reid",fg="red"),sep="")
  #}

  handle = calc.parse(sTree, assoc)
  nodes = handle$nodes
  dmvD = handle$dmvD

  lnP = c()
  for(node in nodes) {
    lnP = c(lnP, log(calchelpers.brprob(sTree, gTree, assoc, 
                                        ploidy, dmvD, node)))
  }

  return(sum(lnP))
}
