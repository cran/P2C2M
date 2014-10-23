ms.exec <-
function(assoc, tree, nTips, ploidy, pFHndl) {
  ######################
  # Function "ms.exec" #
  ######################
  # Descr:  conduction simulations via Hudson's ms   
  # Deps:   calchelpers.nodetips
  # I/p:    tree = the species tree, inferred via the multispecies 
  #                coalescent model
  #         nTips = the total number of tips to be simulated
  #         assoc = a dataframe with two columns and n rows
  #         pFHndle
  #
  # Note:   "assoc": The first column contains the tip labels of the 
  #                  population tree and the second contains the number of 
  #                  alleles to be sampled from that population.
  #
  #         "ploidy": refers to *BEAST's "ploidy" in the xml files
  #                   and modifies the DMV values. When all loci have the same 
  #                   ploidy, it should be left as 1. When ploidy varies, it 
  #                   should be 0.5 for mitoch. and 2 for diploid nuclear.
  #
  #         "tree": all trees must be species trees with associated 
  #                 'dmv' values.
  #
  #         dmv-values: "tree$dmv" = Nu haploid alleles, or 2Nu individuals; 
  #                     hence, multiply by 2 to get the standard 4Nu 
  #                     population scaled mutation rate. BEAST's.dmvparse is 
  #                     equal to Nu, where N is the pop size in alleles. 
  #                     Hence: 2Nu in diploid individuals.

  dbgBool = get("P2C2M.flg.dbgBool", envir=p2c2m.globalVars)
  if (dbgBool) {
    cat("\n", xtermStyle::style("DEBUGMODE> ms.exec", fg="red"), sep="")
  }
###########################################
# 1. Preparatory actions, adjusting theta #
###########################################
  nSp = length(tree$tip.label)

  theta = tree$dmv*2*ploidy
  # relTheta-values are theta-values devided by the root dmv (which is the 
  # last one in the list, hence: tail(x,1)
  relTheta = theta/tail(theta, 1)

  # Branch lengths of the species tree also must be divided by root theta
  tree$edge.length = tree$edge.length/tail(theta, 1)
  # Set up a data frame of branching times, where the row names (which equal
  # the node number) are their own column
  brTms = as.data.frame(sort(ape::branching.times(tree)))
  brTms[,2] = rownames(brTms)
  # nodeNames is a list of node names plus the root
  nodeNames = c(tree$edge[,2], (nSp+1))

##############################
# 2. Converting associations #
##############################
  # FROM: alpinus,m206395a
  #       alpinus,m206395b
  #       alpinus,m206400a
  #       alpinus,m206400b
  # TO:   alpinus, 4
  assoc = as.matrix(as.data.frame(table(assoc[,1])))

###################################
# 3. Set initial pop sizes for ms #
###################################
  popSizes = ""
  for(i in 1:nSp) {
    spName = which(tree$tip.label==assoc[i,1])
    # tree$dmv and tree$edge are in the same order
    size_subpop = relTheta[which(tree$edge[,2]==spName)]
   popSizes = paste(popSizes, "-n", i, size_subpop)
  }

##############################
# 4. Set island model for ms #
##############################
  islModel = paste("-I", length(tree$tip.label), paste(assoc[,2], collapse=' '))

####################################
# 5. Set tree specification for ms #
####################################
  treeSpec = ""
  for (row in 1:nrow(brTms)) {

    # "brTms" are the ages of all nodes, sorted from smallest to largest
    # "nodeNames" is a list of node names plus the root
    brTmName = brTms[row,2]
    brTmVal = brTms[row,1]

    # Calculate a scaling factor
    scaledBrTme = relTheta[which(nodeNames==brTmName)]

    # Take those branches that have the shortes length
    children = sort(tree$edge[tree$edge[,1] == brTmName, 2])
    child1 = sort(calchelpers.nodetips(tree, children[1]))[1]
    child2 = sort(calchelpers.nodetips(tree, children[2]))[1]
    tips = sort(c(child1, child2))

    # Check which of the elements in 
    # "tree$tip.label[tips[2]]" are also in "assoc[,1]"
    popI = which(assoc[,1]==tree$tip.label[tips[2]])
    popJ = which(assoc[,1]==tree$tip.label[tips[1]])

    treeSpec = paste(treeSpec, "-ej", brTmVal, popI, popJ, 
                               "-en", brTmVal, popJ, scaledBrTme)	
  }

###################################################
# 6. Combining specifications for full ms command #
###################################################
  # Path to Hudson's ms (Hudson 2002)
  pathToMS = system.file("msdir", "ms", package="P2C2M")
  nReps = get("p2c2m.flg.nReps", envir=p2c2m.globalVars)

  fullCmd = paste(pathToMS, nTips, nReps, islModel, popSizes, treeSpec, "-T")
  # Logging ms command to file
  loghelpers.wrt1(fullCmd, pFHndl, colFlg=F)

#################
# 7. Execute ms #
#################
  # The grep command catches the tree spec via the trailing semicolon;
  # the first backslash necessary due to R, the second due to bash shell
  msOutp = system(paste(fullCmd, "| grep \\;"), intern=T)

##############################
# 8. Read in simulated trees #
##############################
  if (nReps == 1) {
    trees = list()
    trees[[1]] = ape::read.tree(text=msOutp)
  }
  else {
    trees = ape::read.tree(text=msOutp)
  }

  # Undoing the adjustment of species tree brlens 
  # as was necessary for ms-simulation
  if (length(trees) > 1) {
    for (t in 1:length(trees)) {
      trees[[t]]$edge.length = trees[[t]]$edge.length*tail(theta, 1)
    }
  }
  else {
    trees$edge.length = trees$edge.length*tail(theta, 1)
  }

  return(trees)
}
