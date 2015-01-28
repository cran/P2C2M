p2c2m.analyze <-
function (inFn, inData, pFHndl) {
  # Descr:  core of script, coordinates the metric calculation
  # Deps:   (various)
  # I/p:    inFn
  #         inData
  #         pFHndl

  dbgBool = get("P2C2M.flg.dbgBool", envir=p2c2m.globalVars)
  mpiBool = get("P2C2M.flg.mpiBool", envir=p2c2m.globalVars)
  vrbBool = get("P2C2M.flg.vrbBool", envir=p2c2m.globalVars)

  if (dbgBool) {
    cat("\n", xtermStyle::style("DEBUGMODE> core.analyze", fg="red"), sep="")
  }

##########################
# 1. Preparatory actions #
##########################
  descrStats = get("P2C2M.flg.descrStats", envir=p2c2m.globalVars)
  nLoci = length(inData$loci)
  NsTrees = length(inData$sTrees)
  singleAllele = get("p2c2m.flg.singleAllele", envir=p2c2m.globalVars)
  # generate empty matrices of dimensions 'matrix'
  post_dist = post_pred_dist = matrix(nrow=NsTrees, ncol=nLoci)
  # label the columns of the empty matrices by the locus names
  colnames(post_dist) = colnames(post_pred_dist) = inData$loci
  # Initilizing data 
  data = list()

##########################################################
# 2. Metrics for the trees of the posterior distribution #
##########################################################
  loghelpers.prnt1("Metrics for the trees of the posterior distribution", pFHndl, tofileFlg=T, 
                   color="green")

 ##########################################################
 # 2.a. Looping through genes, calc. descript. statistics #
 ##########################################################
  # Loop through all gene tree names
  for (j in 1:nLoci) {
    loghelpers.prnt2(paste("Locus '", inData$loci[j], "'", sep=""), pFHndl, 
                     tofileFlg=F, color="blue")
    if (vrbBool) {cat("Analyzing gene trees:", "\n", sep="")}
    # Loop through all species trees
    for (i in 1:NsTrees) {
      # Display tree number on screen
      if (vrbBool) {cat(" ", i, sep="")}
      # Save tree number as replicate ID
      assign("p2c2m.flg.repID", paste("post_dist", i, sep="."), 
             envir = p2c2m.globalVars)

      post_dist[i,j] = corehelpers.metrics(inData$gTrees[[j]][[i]],
                                     inData$pTrees$tree[[i]],
                                     inData$pTrees$names,
                                     inData$sTrees[[i]],
                                     inData$assoc,
                                     inData$ploidy[j],
                                     descrStats,
                                     singleAllele)
    }
    if (vrbBool) {cat("\n")}
  }

 ###########################
 # 2.b. Assembling results #
 ###########################
  for (s in descrStats) {
    # Populating list "data" with the results of the trees of the posterior distribution
    data[[s]]$post_dist$unsort = as.data.frame(rtnlc(post_dist, which(s==descrStats)))
    # Assigning colnames
    colnames(data[[s]]$post_dist$unsort) = inData$loci
    # Generating a sorted result clone
    data[[s]]$post_dist$sorted = apply(data[[s]]$post_dist$unsort, MARGIN=2, sort)
    # Writing results to parameter file
    loghelpers.wrt2(cat(paste(s, "[post_dist][unsort]\n", sep="")), pFHndl)
    loghelpers.wrt2(data[[s]]$post_dist$unsort, pFHndl)
    loghelpers.wrt2(cat(paste(s, "[post_dist][sorted]\n", sep="")), pFHndl)
    loghelpers.wrt2(data[[s]]$post_dist$sorted, pFHndl)
  }

#####################################################################
# 3. Metrics for the trees of the posterior predictive distribution #
#####################################################################
  loghelpers.prnt1("Metrics for the trees of the posterior predictive distribution",
                   pFHndl, tofileFlg=T, color="green")

 ##################################
 # 3.a. Initiate MPI (if applic.) #
 ##################################
  # Determining number of processors (i.e. CPUs) available
  if (Sys.info()['sysname'] == "Linux") {
    nproc = as.numeric(system("nproc", intern=T))
  }
  # MacOS option
  if (Sys.info()['sysname'] == "Darwin") {
    nproc = as.numeric(system("sysctl -n hw.ncpu", intern=T))
  }

  if (mpiBool && nproc >= 3) {
    loghelpers.prnt2(paste("\n", "Initialize parallelization", "|", 
                           "N of processors:", toString(nproc), sep=" "),
                     pFHndl, tofileFlg=F, color="yellow")
    # Specifying one master and nproc-1 slaves
    #Rmpi:: mpi.spawn.Rslaves(nslaves=nproc-1, needlog=F)
    Rmpi:: mpi.spawn.Rslaves(nslaves=nproc-1)
    # Passing separate functions to MPI environment
    mpiEnvir = c(calc.coal_reid, calc.gsi, calc.ndc, calc.parse, calc.coal_liu, 
                 calchelpers.brprob, calchelpers.descend, calchelpers.dmvparse,
                 calchelpers.gtreeparse, calchelpers.nodetips, 
                 calchelpers.probcalc, corehelpers.metrics, corehelpers.repl, 
                 frmtMntse, apTreeshape::as.treeshape.phylo,  
                 ape::branching.times, ape::drop.tip, genealogicalSorting::gsi, 
                 phybase::loglikeSP)
    #lapply(mpiEnvir, Rmpi::mpi.bcast.Robj2slave, all=T)
    lapply(mpiEnvir, Rmpi::mpi.bcast.Robj2slave)
  }

 ##############################
 # 3.b. Looping through genes #
 ##############################
  for (j in 1:nLoci) {
    loghelpers.prnt2(paste("Locus '", inData$loci[j], "'", sep=""), pFHndl, 
                     tofileFlg=T, color="blue")
    if (vrbBool) {cat("Analyzing the posterior predictive trees:", "\n", sep="")}
    nTips = length(inData$gTrees[[j]][[1]]$tip.label)
    for (i in 1:NsTrees) {
      # Display tree number on screen
      if (vrbBool) {cat(" ", i, sep="")}
      # Save tree number as replicate ID
      assign("p2c2m.flg.repID", paste("post_pred_dist", j, i, sep="."), 
             envir = p2c2m.globalVars)


 ############################
 # 3.c. Simulation of trees #
 ############################
      # print(inData$assoc)
      # speciesName   alleleName
      # [1,] "A"      "a1"
      # [2,] "A"      "a10"
      # [3,] "A"      "a11"
      # [4,] "A"      "a12"
      # [5,] "A"      "a2"

      # A constituent-frequency table is required for the simulations in ms,
      # because ms saves the alleles as numbers such as "1", "2", "3", ...
      CFTable = readhelpers.constit(inData$assoc)

      # print(CFTable)
      #       aList  V2
      # Var1  "A"    "1"
      # Var1  "A"    "2"
      # Var1  "A"    "3"

      # Simulate gene trees
      simTrees = ms.exec(CFTable, inData$sTrees[[i]], nTips, 
                         inData$ploidy[j], pFHndl)


      nReps = get("p2c2m.flg.nReps", envir=p2c2m.globalVars)
      ## Possible improvement in TFL: why not "length(simTrees)" instead of 
      ## "as.numeric(nReps)"
      simReps = matrix(nrow=as.numeric(nReps), ncol=1)

      # Loop over simulated replicates
      ## Possible improvement: Use apply function instead of loop
      for (k in 1:nReps) {
        # Generating a simTree version where tips have been replaced with 
        # corresponding species names
        simTreePrnt = corehelpers.repl(simTrees[[k]], CFTable)
        # Write the simulated tree to file
        loghelpers.wrt1(ape::write.tree(simTreePrnt), pFHndl, colFlg=F)
      }

 ##################################################################
 # 3.d. Calculation of descript. statistics depend. on MPI status #
 ##################################################################
      if (mpiBool && nproc >= 3) {
        # Applying function "corehelpers.metrics" in parallel
        #  The second element in mpi.parLapply has to be the function that
        #  the data is to be applied
        simReps = Rmpi::mpi.parLapply(simTrees,
                                      corehelpers.metrics,
                                      inData$pTrees$tree[[i]],
                                      inData$pTrees$names,
                                      inData$sTrees[[i]],
                                      CFTable,
                                      inData$ploidy[j],
                                      descrStats,
                                      singleAllele)
        simReps = as.matrix(as.character(simReps))
      }

      else {
        for (k in 1:nReps) {
        simReps[k,1] = corehelpers.metrics(simTrees[[k]],
                                           inData$pTrees$tree[[i]],
                                           inData$pTrees$names,
                                           inData$sTrees[[i]],
                                           CFTable,
                                           inData$ploidy[j],
                                           descrStats,
                                           singleAllele)
        }
      }

 #########################################
 # 3.e. Parsing of average metric values #
 #########################################
      # Extracting the listcolumns of the matrix simReps via GSO.rtnlc
      sumD = list()
      for (s in descrStats) {
        sumD[[s]] = Mode(rtnlc(simReps, which(s==descrStats)))
      }
      sumD = unlist(sumD)
      names(sumD) = NULL

      # Save the replicate values to variable "post_pred_dist"
      post_pred_dist[i,j] = toString(sumD)
      
    }
  if (vrbBool) {cat("\n")}
  }

 ##################
 # 3.f. Close MPI #
 ##################
  if (mpiBool && nproc >= 3) {
    # Closing slaves
    loghelpers.prnt2(paste("\n", "Close parallelization", "\n", sep=""), 
                     pFHndl, tofileFlg=F, color="yellow")
    Rmpi::mpi.close.Rslaves(dellog=T)
  }

 ###########################
 # 3.g. Assembling results #
 ###########################
  for (s in descrStats) {
    # Populating list "data" with the simulated results
    # Both "as.matrix" and "as.data.frame" are critical below
    data[[s]]$post_pred_dist$unsort = as.matrix(as.data.frame(rtnlc(post_pred_dist, which(s==descrStats))))
    # Assigning colnames
    colnames(data[[s]]$post_pred_dist$unsort) = inData$loci
    # Generating a sorted result clone
    data[[s]]$post_pred_dist$sorted = apply(data[[s]]$post_pred_dist$unsort, MARGIN=2, sort)
    # Writing results to parameter file
    loghelpers.wrt2(cat(paste(s, "[post_pred_dist][unsort]\n", sep="")), pFHndl)
    loghelpers.wrt2(data[[s]]$post_pred_dist$unsort, pFHndl)
    loghelpers.wrt2(cat(paste(s, "[post_pred_dist][sorted]\n", sep="")), pFHndl)
    loghelpers.wrt2(data[[s]]$post_pred_dist$sorted, pFHndl)
  }

#####################################
# 4. Calculating metric differences #
#####################################
  loghelpers.prnt1("Calculating differences btw. the posterior and the posterior predictive distribution",
                   pFHndl, tofileFlg=T, color="green")
  
  for (s in descrStats) {
    # Populating list "data" with the differences
    srtBool = get("P2C2M.flg.srtBool", envir=p2c2m.globalVars)
    if (srtBool) {
      data[[s]]$dif = statshelpers.diffrnce(data[[s]]$post_dist$sorted, 
                                            data[[s]]$post_pred_dist$sorted)
      # Remove the data set not chosen
      data[[s]]$post_dist$unsort = NULL
      data[[s]]$post_pred_dist$unsort = NULL
    }
    else {
      data[[s]]$dif = statshelpers.diffrnce(data[[s]]$post_dist$unsort, 
                                            data[[s]]$post_pred_dist$unsort)
      # Remove the data set not chosen
      data[[s]]$post_dist$sorted = NULL
      data[[s]]$post_pred_dist$sorted = NULL
    }
    # Writing differences to file
    loghelpers.wrt2(cat(paste(s, "[dif]\n", sep="")), pFHndl)
    loghelpers.wrt2(data[[s]]$dif, pFHndl)
  }

  return(data)
}
