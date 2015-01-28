p2c2m.complete <-
function (path="/home/user/Desktop/", xml.file="beast.xml", 
                           descr.stats="COAL_REID,NDC", beast.vers="1.8", 
                           single.allele=c("O"), num.reps=100, 
                           use.sorted=FALSE, use.mpi=FALSE, verbose=FALSE, 
                           dbg=FALSE) {
  # Descr:  initiates the entire package by starting the wrapper
  # Deps:   p2c2m.init, p2c2m.readstarb, p2c2M.analyze, p2c2M.statcalc
  # I/p:    path = absolute path to Working directory
  #         xml.file = name of Beauti XML infile
  #         descr.stats = list of descriptive statistics
  #         beast.vers = version fo BEAST
  #         single.allele = name of outgroup
  #         num.reps = number of replicates
  #         use.sorted = decision if to use sorting
  #         use.mpi = decision if to use MPI
  #         verbose = flag if status output printed on screen
  #         dbg = debugging flag

  startTime = Sys.time()

  # Check if selected descriptive stats are valid
  descrStats = toupper(sort(unlist(strsplit(descr.stats, split=","))))
  selection = which(descrStats %in% c("COAL_REID","COAL_LIU","GSI","NDC"))

  # Error handling
  if (length(descrStats) != length(selection)) {
      stop(cat("\nERROR: Incorrect specification of descript. statistic(s)\n"))
  }

  # Set up variables and options
  p2c2m.init(xml.file, descrStats, beast.vers, single.allele, num.reps, 
             use.sorted, use.mpi, verbose, dbg)

####################################
# 1. Opening parameter file handle #
####################################
  dbgBool = get("P2C2M.flg.dbgBool", envir=p2c2m.globalVars)

  if (dbgBool) {
    # Opening parameter file handle ("prmFile" = "parameter file handle")
    # In command "file", "a" stands for append
    prmFile = file(paste(rmext(xml.file), ".P2C2M.prmts.txt", sep=""), "a")
    # Write name of infile to parameter file handle
    writeLines(xml.file, prmFile)
    cat("\n")
  }

##########################
# 2. Conducting analyses #
##########################
  # Generate indata
  treeData = p2c2m.readstarb(path, xml.file, prmFile)

  if (dbgBool) {
    cat("\n", xtermStyle::style(
        "DEBUGMODE> *** SAVING ALL LOADED TREES ***", 
        fg="red"), "\n", sep="")
    p2c2m.dbg.parsedInput = treeData
    save(p2c2m.dbg.parsedInput, file = paste(rmext(xml.file), 
         ".p2c2m.dbg.parsedInput.rda", sep=""))
  }

  # Generate results data
  descrData = p2c2m.analyze(xml.file, treeData, prmFile)
  # Caclulate statistics
  resultData = p2c2M.statcalc(path, xml.file, treeData$loci, 
                              descrData, prmFile)

######################################################
# 3. Printing runtime, closing parameter file handle #
######################################################
  endTime = Sys.time()

  if (dbgBool) {
  # Write runtime to parameter file
    runTime = paste("\n", "Analysis started at:", startTime,
                    "\n", "Analysis ended at:", endTime)
    loghelpers.wrt2(cat(runTime), prmFile)

    close(prmFile)
  }

##########################
# 4. Summarizing outdata #
##########################
  outD = list()
  outD$rawStats = descrData
  outD$results = resultData

return(outD)
}
