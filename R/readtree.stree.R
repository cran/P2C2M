readtree.stree <-
function(inFn) {
  # Descr:  read species tree
  # Deps:   ?
  # I/p:    inFn
  # Note:   little change compared to Noah's scripts

  dbgBool = get("P2C2M.flg.dbgBool", envir=p2c2m.globalVars)
  vrbBool = get("P2C2M.flg.vrbBool", envir=p2c2m.globalVars)

  if (dbgBool) {
    cat("\n", xtermStyle::style("DEBUGMODE> readtree.stree", fg="red"), sep="")
  }

  treeD = readhelpers.parse(inFn)
  # extract only the tree specification
  treeD$data = gsub("tree STATE_[[:digit:]]+[[:space:]]\\[[^]]*\\] = \\[&R\\] ",
                    "", treeD$data)
  if (vrbBool) {
    cat("\tN of tree specs loaded: ", length(treeD$data), "\n", sep="")
  }
  brD = lapply(treeD$data, readhelpers.extract, locus=NA, sTreeFlg=T)

  tree = treeD$info
  # Notes by Noah: Appending branch width stats to trees; vector "check" is 
  # included as a verification and should exactly match vector "edge.length"
  if(length(brD[[1]][1,])==6) {
    if(class(tree)=="multiPhylo") {
      for(i in 1:length(tree)) {
        # tree[[i]][["check"]]<-as.numeric(brD[[i]][,2])
        tree[[i]][["dmt"]] = as.numeric(brD[[i]][,3])
        tree[[i]][["dmv_start"]] = as.numeric(brD[[i]][,4])
        tree[[i]][["dmv_end"]] = as.numeric(brD[[i]][,5])
        tree[[i]][["dmv"]] = as.numeric(brD[[i]][,6])
      }
    }
    if(class(tree)=="phylo") {
      #	tree[["check"]]<-as.numeric(brD[[1]][,2])
      tree[["dmt"]] = as.numeric(brD[[1]][,3])
      tree[["dmv_start"]] = as.numeric(brD[[1]][,4])
      tree[["dmv_end"]] = as.numeric(brD[[1]][,5])
      tree[[i]][["dmv"]] = as.numeric(brD[[i]][,6])
    }
  }

  if(length(brD[[1]][1,])==3) {	
    if(class(tree)=="multiPhylo") {
      for(i in 1:length(tree)) {
        # tree[[i]][["check"]]<-as.numeric(brD[[i]][,2])
        tree[[i]][["dmv"]]<-as.numeric(brD[[i]][,3])
      }
    }
    if(class(tree)=="phylo") {
      #	tree[["check"]]<-as.numeric(brD[[1]][,2])
      tree[["dmv"]] = as.numeric(brD[[1]][,3])
    }
  }

  return(tree)
}
