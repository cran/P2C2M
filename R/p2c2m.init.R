p2c2m.init <-
function (xml.file, descr.stats = c("GTP","NDC"), 
                       beast.vers="1.8", single.allele=c("O"), num.reps=100, 
                       use.sorted=F, use.mpi=F, verbose=F, dbg=F) {
  # Descr:  initiates the global variables
  # Deps:   (none)
  # I/p:    xml.file = xml infile
  #         beast.vers = version fo BEAST
  #         out.grp = name of outgroup
  #         num.reps = number of replicates
  #         use.sorted = decision if to use sorting
  #         use.mpi = decision if to use MPI
  #         dbg = debugging flag

#########################
# 1. Set global options #
#########################
  # TFL avoids line-wrapping in command "print()"
  options(width=1000)
  # TFL avoids numbers in scientific e-notation (i.e. "1.2345e+03" = 1234.5)
  options(scipen=1000)

################################
# 2. Set environment variables #
################################
  # Note: I am not using true global variables (i.e. envir = .GlobalEnv), 
  #       because such could have interfered with user-space.
  assign("p2c2m.flg.xmlFile", xml.file, envir = p2c2m.globalVars)
  assign("P2C2M.flg.descrStats", descr.stats, envir = p2c2m.globalVars)
  assign("p2c2m.flg.beastVers", beast.vers, envir = p2c2m.globalVars) 
  assign("p2c2m.flg.singleAllele", single.allele, envir = p2c2m.globalVars)
  assign("p2c2m.flg.nReps", num.reps, envir = p2c2m.globalVars)
  assign("P2C2M.flg.srtBool", use.sorted, envir = p2c2m.globalVars)
  assign("P2C2M.flg.mpiBool", use.mpi, envir = p2c2m.globalVars)
  assign("P2C2M.flg.vrbBool", verbose, envir = p2c2m.globalVars)
  assign("P2C2M.flg.dbgBool", dbg, envir = p2c2m.globalVars)
}
