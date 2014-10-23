readtree.phybase <-
function(inFn) {
  # Descr:  read tree of format phybase
  # Deps:   (none)
  # I/p:    inFn

  beastVers = get("p2c2m.flg.beastVers", envir=p2c2m.globalVars)
  dbgBool = get("P2C2M.flg.dbgBool", envir=p2c2m.globalVars)

  if (dbgBool) {
    cat("\n", xtermStyle::style("DEBUGMODE> readtree.phybase", fg="red"), 
        sep="")
  }

  # Specify name of parsing script
  pySc = system.file("exec", paste("BEAST2phybase_", beastVers, ".py", sep=""), 
                     package="P2C2M")
  system(paste("python2", pySc, inFn))
  # Read tree via PHYBASE
  pTrees = phybase::read.tree.string(paste(rmext(inFn), ".P2C2M.phyb", sep=""))
  # Remove leading white space from tree specs
  pTrees$tree = gsub("^ .", "", pTrees$tree)

  return(pTrees)
}
