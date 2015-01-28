p2c2m.readstarb <-
function (apWd, inFn, pFHndl) {
  # Descr:  coordinates the reading of all input      
  # Deps:   readtree.stree
  #         readtree.phybase
  # I/p:    apWd = absolute path to working directory
  #         inFn = infile name (name of Beauti XML outd file)
  #         inLogFn
  #         pFHndl

  descrStats = get("P2C2M.flg.descrStats", envir=p2c2m.globalVars)
  dbgBool = get("P2C2M.flg.dbgBool", envir=p2c2m.globalVars)
  vrbBool = get("P2C2M.flg.vrbBool", envir=p2c2m.globalVars)

  if (dbgBool) {
    cat("\n", xtermStyle::style("DEBUGMODE> readstarb.read", fg="red"), sep="")
  }

#####################################
# 1. Parsing data via Python script #
#####################################
  loghelpers.prnt1("Parsing data from xml inputfile", pFHndl, tofileFlg=T, 
                   color="green")

  # Specify name of parsing script
  beastVers = get("p2c2m.flg.beastVers", envir=p2c2m.globalVars)
  pySc = system.file("exec", 
               paste("parseBEAUTiXML_", beastVers, ".py", sep=""), 
               package="P2C2M")
  # Load infile in Python
  rPython::python.assign("inFn", paste(apWd, inFn, sep=""))
  rPython::python.load(pySc)

######################################
# 2. Saving parsed info to variables #
######################################
  # Loading data of Species-Allele matrix into variable "SAmtrx"
  SAmtrx = do.call(rbind, rPython::python.get("SAmtrx"))
  colnames(SAmtrx) = c("speciesName", "alleleName")

  # Loading data of Locus-Ploidy-Matrix into variable "LPmtrx"
  LPmtrx = do.call(rbind, rPython::python.get("LPmtrx"))
  colnames(LPmtrx) = c("locus", "ploidy")
  loci = LPmtrx[,1]
  ploidy = as.numeric(LPmtrx[,2])

  # CURRENTLY PUT ON HOLD:
  #    # Loading data of Gene-Partition-matrix into variable "GPmtrx"        
  #    GPmtrx = do.call(rbind, rPython::python.get("GPmtrx"))
  #    colnames(GPmtrx) = c("genes","partitions")

  # CURRENTLY PUT ON HOLD:
  #    # Loading alignment-dictionary into variable "ALNdict"
  #    alignm = list()
  #    ALNdict = rPython::python.get("ALNdict")
  #    for (i in 1:length(loci))
  #        {
  #        # "order(names(ALNdict[[i]]))" assures that genes and loci are 
  #        # arranged in the same order as in LPmtrx
  #        listorder = order(names(ALNdict[[i]]))
  #        alignm[[loci[i]]] = ALNdict[[i]][listorder]
  #        }

  if (dbgBool) {
    writeLines(paste("\n","SPECIES ALLELE MATRIX (col.names=T)", sep=""), 
               pFHndl)
    loghelpers.wrt2(SAmtrx, pFHndl)
    writeLines(paste("\n", "LOCI", sep=""), pFHndl)
    loghelpers.wrt2(loci, pFHndl)
    writeLines(paste("\n", "PLOIDY", sep=""), pFHndl)
    loghelpers.wrt2(ploidy, pFHndl)
  }

##################################################
# 3. Loading gene trees via function "readgtree" #
##################################################
  loghelpers.prnt1("Loading gene trees", pFHndl, tofileFlg=T, color="green")

  gTrees = list()
  for (locus in loci) {
    gTreeName = paste(apWd, locus, ".trees",sep="")
    loghelpers.prnt2(paste("Empirical gene trees of locus '", locus,"'",sep=""),
                     pFHndl, tofileFlg=T, color="blue")
    # Loading gene trees via function "readgtree"
    gTrees[[locus]] = readtree.gtree(gTreeName, locus)
    loghelpers.wrt1(ape::write.tree(gTrees[[locus]]), pFHndl, colFlg=F)

    if (vrbBool) {
      cat("\t", "N of gene trees loaded: ", length(gTrees[[locus]]), "\n", 
          sep="")
    }
  }

####################################
# 4. Loading regular species trees #
####################################
  loghelpers.prnt1("Loading species trees", pFHndl, tofileFlg=T, color="green")

  streeName = paste(apWd, "species.trees", sep="")
  loghelpers.prnt2(paste("Regular species trees", sep=""), pFHndl, tofileFlg=T, 
                   color="blue")
  # Load sTrees into string (not list!)
  sTrees = readtree.stree(streeName)
  loghelpers.wrt1(ape::write.tree(sTrees), pFHndl, colFlg=F)

####################################
# 5. Loading phybase species trees #
####################################

  # See if flag descrStats contains "COAL_LIU", for which phybase species trees are
  # needed
  if ("COAL_LIU" %in% descrStats) {
    loghelpers.prnt2(paste("Phybase species trees",sep=""), pFHndl, tofileFlg=T, 
                       color="blue")

    # Load pTrees as list
    pTrees = readtree.phybase(streeName)
    loghelpers.wrt1(pTrees$tree, pFHndl, colFlg=F)

    if (vrbBool | dbgBool) {  
      if (length(sTrees) != length(pTrees$tree)) {
        cat("\n", xtermStyle::style("ERROR: N(regular species trees) != \
N(phybase species trees)", fg="red"), "\n", sep="")
      }
      if (length(sTrees) != length(pTrees$tree)) {
        cat("\t","N of species trees loaded: ",length(sTrees),"\n",sep="")
      }
    }
  }

  if (!"COAL_LIU" %in% descrStats) {
    pTrees = NULL
  }

# CURRENTLY PUT ON HOLD:
#     #######################
#     # 6. Loading log file #
#     #######################
#       loghelpers.prnt1("Loading logfile", pFHndl, tofileFlg=F, color="green")
#       log = read.table(paste(apWd, inLogFn, sep=""), header = T)

#####################################################
# 7. Setting up ordered allele-species associations #
#####################################################
  loghelpers.prnt1("Generating ordered allele-species associations", pFHndl, 
                     tofileFlg=T, color="green")

  # FOR DEVELOPER: Entire section could potentially be replace with:
  #                assoc = SAmtrx

  assoc = c()
  # Get tips from first tree of first locus in variable "gTrees"
  tips = gTrees[[1]][[1]]$tip.label
  for (tip in tips) {
    assoc = rbind(assoc, SAmtrx[which(SAmtrx[,2]==tip),])
  }
  assoc = as.matrix(as.data.frame(assoc))
  loghelpers.wrt2(assoc, pFHndl)

###################################
# 8. Saving data to variable outD #
###################################

  if (dbgBool) {
    # With 20 input trees, 5% are 2 trees
    if (length(sTrees) > 39) {
      cat("\n", xtermStyle::style("DEBUGMODE> *** REDUCING DATA SET TO 5% ***",
          fg="red"), "\n", sep="")
      reduced = length(sTrees)*0.05
      pTrees = head(pTrees, reduced)
      sTrees = head(sTrees, reduced)
      gTrees = lapply(gTrees, head, reduced)
    }
  }

  outD = list()
  # CURRENTLY PUT ON HOLD: #class(outD) = "starbeast.data"

  # CURRENTLY PUT ON HOLD: # outd$alignm = alignm
  outD$assoc = assoc
  outD$loci = loci
  outD$gTrees = gTrees
  # CURRENTLY PUT ON HOLD: # outD$log = log
  # CURRENTLY PUT ON HOLD: # outD$parts = GPmtrx
  outD$ploidy = ploidy
  outD$pTrees = pTrees
  outD$sTrees = sTrees

  return(outD)
}
