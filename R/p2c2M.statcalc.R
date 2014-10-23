p2c2M.statcalc <-
function (path, xml.file, loci, resultData, prmFile) {
  # Descr:    coordinates executing of modules
  # Deps:     (various)
  # I/p:      path = absolute path to Working directory
  #           xml.file = name of infile
  #           loci = names of loci
  #           prmFile

  dbgBool = get("P2C2M.flg.dbgBool", envir=p2c2m.globalVars)

  if (dbgBool) {
    cat("\n", xtermStyle::style("DEBUGMODE> p2c2M.statcalc", fg="red"),
        "\n", sep="")
  }

#####################################
# 1. Setting different alpha-values #
#####################################
  loghelpers.prnt1("Summarizing results", prmFile, tofileFlg=T,
                   color="green")

  if (dbgBool) {
    # Initiating results file
    outFn = paste(rmext(xml.file), ".P2C2M.rslts.csv", sep="")
    write(paste(xml.file, "\n", sep=""), outFn, append=T)
  }

  # Setting outdata
  outD = list()

  # Setting alpha values
  alphaValues = c(0.1, 0.05, 0.01)
  for (val in alphaValues) {

#############################
# 2. Calculating statistics #
#############################
    results = stats.main(resultData, loci, val)
    valStr = paste("alpha", as.character(val), sep="")
    outD[[valStr]] = results

    if (dbgBool) {
      # Writing to PRMT file
      # "\u03B1" stands for unicode character "alpha"
      loghelpers.wrt2(cat(paste("\n", "\u03B1 = ", val, sep="")), prmFile)
      loghelpers.wrt2(results, prmFile)
      # Writing to RSLT file
      write(paste("\n", "\u03B1 = ", val, sep=""), outFn, append=T)
      # Note: Without suppressWarnings(), the follow. error msg. would be 
      # printed: "In write.table(outD, outFn, append = T, row.names = T, 
      # col.names = T: appending column names to file)"
      # Note: "col.names=NA" is important, because the column name for a column 
      # of row names is filled with ""
      suppressWarnings(write.table(results, outFn, append=T, row.names=T,
                                   col.names=NA, sep=","))
    }
  }

#####################
# 3. Writing legend #
#####################
  legend = "Differences are calculated as `empirical-simulated` (i.e. the stronger the number deviates from zero, the stronger the deviation of the empirical tree from the coalescent model). Each cell contains the following information in said order: mean, standard deviation and significance level. Codes in square brackets indicate the number of tails. Alpha values are automatically adjusted for the number of tails."
  outD$legend = legend

  if (dbgBool) {
    write(paste("\n", legend, "\n", sep=""), outFn, append=T)
  }

return(outD)
}
