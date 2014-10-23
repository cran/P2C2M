loghelpers.prnt2 <-
function(inD, pFHndl, tofileFlg, color) {
  ###############################
  # Function "loghelpers.prnt2" #
  ###############################
  # Descr:  printing log; regular
  # Deps:   -
  # I/p:    inD
  #         pFHndl
  #         tofileFlg
  #         color

  dbgBool = get("P2C2M.flg.dbgBool", envir=p2c2m.globalVars)
  vrbBool = get("P2C2M.flg.vrbBool", envir=p2c2m.globalVars)

  if (vrbBool) {
    # Logging to screen
    cat("\n", xtermStyle::style(inD, fg=color), "\n", sep="")

    # Writing to parameter file
    dbgBool = get("P2C2M.flg.dbgBool", envir=p2c2m.globalVars)
    if (dbgBool) {
      if (tofileFlg) {
        writeLines(paste("\n",toupper(inD), sep=""), pFHndl)
      }
    }
  }
}
