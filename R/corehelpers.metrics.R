corehelpers.metrics <-
function(gTree, pTree, pTreeNames, sTree, 
                               assoc, ploidy, descrStats, singleAllele) {
  # Descr:    corehelpers.metrics
  # Deps:     calc.coal_reid
  #           calc.ndc
  #           calc.coal_liu
  # I/p:      gTree
  #           pTree
  #           pTreeNames
  #           sTree
  #           assoc
  #           ploidy
  #           descrStats
  #           singleAllele

  #dbgBool = get("P2C2M.flg.dbgBool", envir=p2c2m.globalVars)
  #if (dbgBool) {
  #  cat("\n", xtermStyle::style("DEBUGMODE> corehelpers.metrics", fg="red"), sep="")
  #}

  outL = list()

  if ("GSI" %in% descrStats) {
    gsi = calc.gsi(gTree, assoc, singleAllele)
    # Setting a specific number of significands
    outL$GSI = frmtMntse(gsi, 4)
  }

  if ("COAL_REID" %in% descrStats) {
    gtp = calc.coal_reid(gTree, sTree, assoc, ploidy)
    outL$COAL_REID = frmtMntse(gtp, 4)
  }

  if ("NDC" %in% descrStats) {
    ndc = calc.ndc(gTree, sTree, assoc)
    outL$NDC = frmtMntse(ndc, 4)
  }

  if ("COAL_LIU" %in% descrStats) {
    ray = calc.coal_liu(gTree, pTree, pTreeNames, assoc)
    outL$COAL_LIU = frmtMntse(ray, 4)
  }

  outD = unlist(outL)
  names(outD) = NULL
  outD = toString(outD)

  return(outD)
}
