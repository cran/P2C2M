corehelpers.metrics <-
function(gTree, pTree, pTreeNames, sTree, 
                               assoc, ploidy, descrStats, singleAllele) {
  # Descr:    corehelpers.metrics
  # Deps:     calc.gtp
  #           calc.ndc
  #           calc.ray
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
  #  cat("\n", xtermStyle::style("DEBUGMODE> corehelpers.emp", fg="red"), sep="")
  #}

  outL = list()

  if ("GSI" %in% descrStats) {
    gsi = calc.gsi(gTree, assoc, singleAllele)
    # Setting a specific number of significands
    outL$GSI = frmtMntse(gsi, 4)
  }

  if ("GTP" %in% descrStats) {
    gtp = calc.gtp(gTree, sTree, assoc, ploidy)
    outL$GTP = frmtMntse(gtp, 4)
  }

  if ("NDC" %in% descrStats) {
    ndc = calc.ndc(gTree, sTree, assoc)
    outL$NDC = frmtMntse(ndc, 4)
  }

  if ("RAY" %in% descrStats) {
    ray = calc.ray(gTree, pTree, pTreeNames, assoc)
    outL$RAY = frmtMntse(ray, 4)
  }

  outD = unlist(outL)
  names(outD) = NULL
  outD = toString(outD)

  return(outD)
}
