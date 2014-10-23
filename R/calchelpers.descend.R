calchelpers.descend <-
function (sTree, gTree, assoc, cNode) {
##################################
# Function "calchelpers.descend" #
##################################
# Descr:    returns a tree containing descendants of all gene lineages 
#           that pass through a node; here this function is used to 
#           obtain the gene tree that is part of a species tree
# Deps:     calchelpers.nodetips
# I/p:      sTree
#           gTree
#           assoc
#           cNode = coalescence node

#  dbgBool = get("P2C2M.flg.dbgBool", envir=p2c2m.globalVars)
#  if (dbgBool) {
#    cat("\n", xtermStyle::style("DEBUGMODE> calchelpers.descend", fg="red"),
#        sep="")
#  }

    spNames = sTree$tip.label
    # TFL means "Unless cNode is 'NA' AND unless the coalescence node is 
    # the root, do ..."
    if (!is.na(cNode) && cNode != (length(spNames)+1))
        {
        subtreeTaxa = spNames[calchelpers.nodetips(sTree, cNode)]
        gTreeDscr = c()
        for(taxon in subtreeTaxa)
            {
            nodesRaw = assoc[which(assoc[,1] == taxon),2]
            if (is.character(nodesRaw[1]) == T)
                {nodes = nodesRaw}
            if (is.character(nodesRaw[1]) == F)
                {nodes = as.numeric(nodesRaw)}
            gTreeDscr = c(gTreeDscr, nodes)
            }
        taxaIndex = match(gTreeDscr, gTree$tip.label)
        # TFL removes all terminal branches of the gene tree 
        # that are not in the subtree; gTree must be in DNAbin-format
        subtree = ape::drop.tip(gTree, gTree$tip.label[-taxaIndex])
        }
    else {subtree = gTree}

return(subtree)
}
