loghelpers.wrt2 <-
function(inD, pFHndl) {
  ##############################
  # Function "loghelpers.wrt2" #
  ##############################
  # Descr:  write table under formatting scheme 2
  # Deps:   -
  # I/p:    inD
  #         pFHndl
  #         colFlg

  dbgBool = get("P2C2M.flg.dbgBool", envir=p2c2m.globalVars)

  if (dbgBool) {
    sink(pFHndl, append=T)
    cat("\n")
    print(inD)
    sink()
  }
}
