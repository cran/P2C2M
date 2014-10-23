readhelpers.parse <-
function(inFn) {
  # Descr:  dealing with metadata
  # Deps: -
  # I/p:  inFn

  dbgBool = get("P2C2M.flg.dbgBool", envir=p2c2m.globalVars)
  if (dbgBool) {
    cat("\n",xtermStyle::style("DEBUGMODE> readhelpers.parse",fg="red"),sep="")
  }

  tree = ape::read.nexus(inFn)
  # Load inFn line by line
  lines = scan(file = inFn, what = "", sep = "\n", quiet=T)
  # Extract all those lines that contain keyword "tree STATE"
  outD = c()
  outD$"info" = tree
  outD$"data" = lines[grep("tree STATE", lines)]

  #outd[["info"]] = tree
  #outd[["data"]] = lines[grep("tree STATE", lines)]

  return(outD)
}
