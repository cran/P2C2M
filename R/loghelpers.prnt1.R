loghelpers.prnt1 <-
function(inD, pFHndl, tofileFlg, color) {
  ###############################
  # Function "loghelpers.prnt1" #
  ###############################
  # Descr:  printing log; with canvas of pound signs
  # Deps:   -
  # I/p:    inD
  #         pFHndl
  #         tofileFlg
  #         color

  dbgBool = get("P2C2M.flg.dbgBool", envir=p2c2m.globalVars)
  vrbBool = get("P2C2M.flg.vrbBool", envir=p2c2m.globalVars)

  if (vrbBool) {
    # Logging to screen
    cat("\n", xtermStyle::style(toupper(inD), fg=color), "\n", sep="")

    # Writing to parameter file
    if (dbgBool) {
      if (tofileFlg) {
        tmpD = paste(rep("#", length(unlist(strsplit(inD, split="")))), collapse="")
        writeLines(paste("\n", "##", tmpD, "##", sep=""), pFHndl)
        writeLines(paste("#", toupper(inD), "#"), pFHndl)
        writeLines(paste("##", tmpD, "##", sep=""), pFHndl)
      }
    }
  }
}
