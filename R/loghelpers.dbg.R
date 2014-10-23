loghelpers.dbg <-
function(inList, varName, title) {
  #############################
  # Function "loghelpers.dbg" #
  #############################
  # Descr:  write debug data to file
  # Deps:   -
  # I/p:    inList
  #         title

  dbgBool = get("P2C2M.flg.dbgBool", envir=p2c2m.globalVars)

  if (dbgBool) {
    cat("\n", xtermStyle::style(paste(
        "DEBUGMODE> *** SAVING", title, "TO FILE ***"), 
        fg="red"), "\n", sep="")

    # Readying the data to be saved
    outVarName = paste("p2c2m.dbg.", varName, sep="")
    assign(outVarName, inList)

    # Specifying the outFileName
    xmlFile = get("p2c2m.flg.xmlFile", envir=p2c2m.globalVars)
    outFileName = paste(rmext(xmlFile), ".", outVarName, ".rda", sep="")

    # Saving initial R object
    if (!file.exists(outFileName)) {
      save(outVarName, file=outFileName)
    }

    # Appending to existing R object
    if (file.exists(outFileName)) {
      old.objects = load(outFileName)
      new.objects = c(old.objects, outVarName)
      save(new.objects, file = outFileName)
    }
  #  if (file.exists(outFileName)) {
  #    old.objects = load(outFileName)
  #    save(list=c(old.objects, outVarName), file = outFileName)
  #  }
  }
}
