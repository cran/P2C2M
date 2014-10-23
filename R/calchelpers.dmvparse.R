calchelpers.dmvparse <-
function(sTree, nSp) {
  # Descr:    parsing demographic data (dmv and dmt) of a branch
  # Deps:     -
  # I/p:      sTree
  #           nSp

#  dbgBool = get("P2C2M.flg.dbgBool", envir=p2c2m.globalVars)
#  if (dbgBool) {
#    cat("\n", xtermStyle::style("DEBUGMODE> calchelpers.dmvparse", fg="red"),
#        sep="")
#  }

    # generating rows of node number - branch length pairs
    dmvD = cbind(sTree$edge[,2], sTree$edge.length)
    # adding another branch of length Inf
    dmvD = rbind(dmvD, c((nSp+1), Inf))
    # adding dmv values
    dmvD = cbind(dmvD, sTree$dmv)
    # order the matrix by the first column
    dmvD = dmvD[order(dmvD[,1]),]
    # calc branchling times of the species tree (via ape-function)
    stBt = ape::branching.times(sTree)

    # sort the branching times by their names (which are numbers)
    # TFL may not be necessary, as stBt may already be sorted
    stBt = stBt[order(as.numeric(names(stBt)))]

    # add x zeros to the beginning of branching times list, 
    # where x = number of species in species tree
    pre = structure(rep(0, nSp),.Names=c(1:nSp))
    stBt = c(pre, stBt)

    # add the column "stBt" to the matrix "dmvD"
    dmvD = cbind(dmvD, stBt)
    colnames(dmvD) = c("node", "length", "dmv", "sbt")
    rownames(dmvD) = c(1:length(stBt))

return(dmvD)
}
