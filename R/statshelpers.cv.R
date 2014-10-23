statshelpers.cv <-
function(diff) {
  # Descr:    Calculates the coefficient of variance
  # Deps:     (none)
  # I/p:      inD
  dbgBool = get("P2C2M.flg.dbgBool", envir=p2c2m.globalVars)
  if (dbgBool) {
    cat("\n", xtermStyle::style("DEBUGMODE> statshelpers.cv", fg="red"), sep="")
  }

  # Structure of "diff" at this point:
  #        [gene1]   [gene2]  [gene3]
  # [1,]  -2.72395  -9.48870 45.79565
  # [2,]  14.97560  38.06380 88.60285
  # [3,]   8.76280  21.09895 50.51950

#  inD = as.data.frame(diff)
#  inD = ifelse(is.nan(diff), NA, inD)
  Stdv = apply(diff, MARGIN=1, sd, na.rm=T)
  Mean = rowMeans(diff)
  cv = Stdv/Mean

  return(cv)
}
