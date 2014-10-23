calchelpers.nodetips <-
function(tree, node) {
  # Descr:    returns the tip numbers for a given node
  # Deps:     calchelpers.nodetips
  # I/p:      tree = a phylog. tree with numbered internal nodes and numbered tips
  #           node = a node in that tree, specified as an integer

#  dbgBool = get("P2C2M.flg.dbgBool", envir=p2c2m.globalVars)
#  if (dbgBool) {
#    cat("\n", xtermStyle::style("DEBUGMODE> calchelpers.nodetips", fg="red"),
#        sep="")
#  }

    nTips = length(tree$tip.label)
    if (node <= nTips) {node}
    else 
        {
     	outd = numeric()
     	k = tree$edge[which(tree$edge[,1] == node),2]
     	for (j in k)
            {
            if (j <= nTips) {outd = c(outd, j)}
            else {outd = c(outd, calchelpers.nodetips(tree, j))}
            }
        outd        
        }
}
