readtree.gtree <-
function(inFn, locus) {
  # Descr:  read gene tree    
  # Deps:   ?     
  # I/p:    inFn
  #         locus
  # Note:   little change compared to Noah's scripts

  dbgBool = get("P2C2M.flg.dbgBool", envir=p2c2m.globalVars)
  vrbBool = get("P2C2M.flg.vrbBool", envir=p2c2m.globalVars)

  if (dbgBool) {
    cat("\n", xtermStyle::style("DEBUGMODE> readtree.gtree", fg="red"), sep="")
  }

  treeD = readhelpers.parse(inFn)
  # Remove everything in lines except tree definition				  
  treeD$data = gsub("tree STATE_.*\\[&R\\] ", "", treeD$data)

  if (vrbBool) {
    cat("\tN of tree specs loaded: ", length(treeD$data), "\n", sep="")
  }
  brD = lapply(treeD$data, readhelpers.extract, locus, sTreeFlg=F)

  tree = treeD$info

  # Append branch rate stats to trees
  if(class(tree)=="multiPhylo") {
    for(i in 1:length(tree)) {
      mode(brD[[i]]) = "numeric"
      tree[[i]][["rate"]] = as.matrix(brD[[i]][,c(-1,-2)])
    }
  }
  if(class(tree)=="phylo") {
    mode(brD[[1]]) = "numeric"
    tree[["rate"]] = as.numeric(brD[[1]][,c(-1,-2)])
  }

  return(tree)
}
