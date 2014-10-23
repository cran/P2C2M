corehelpers.repl <-
function (simTree, assoc) {
  # Descr:    corehelpers.repl
  # Deps:     -
  # I/p:      simTree
  #           assoc

  #dbgBool = get("P2C2M.flg.dbgBool", envir=p2c2m.globalVars)
  #if (dbgBool) {
  #  cat("\n", xtermStyle::style("DEBUGMODE> corehelpers.repl", fg="red"), 
  #      sep="")
  #}

  nwLbls = assoc[,1][match(simTree$tip.label, assoc[,2])]
  simTree$tip.label = sapply(as.list(as.data.frame(nwLbls))[[1]], 
                             function(x){toString(x)})

  return(simTree)
}
