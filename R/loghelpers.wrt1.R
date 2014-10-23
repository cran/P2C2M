loghelpers.wrt1 <-
function(inD, pFHndl, colFlg) {
  ##############################
  # Function "loghelpers.wrt1" #
  ##############################
  # Descr:  write table under formatting scheme 1
  # Deps:   -
  # I/p:    inD
  #         pFHndl
  #         colFlg

  dbgBool = get("P2C2M.flg.dbgBool", envir=p2c2m.globalVars)

  if (dbgBool) {
    write.table(inD, pFHndl, append=T, row.names=F, col.names=colFlg, 
                quote=F, sep=",")
  }
}
