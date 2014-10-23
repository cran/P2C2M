statshelpers.diffrnce <-
function (emp, sim) {
  ####################################
  # Function "statshelpers.diffrnce" #
  ####################################
  # Descr:    calculates the difference between empirical and simulated
  # Deps:     -
  # I/p:      emp
  #           sim

  dbgBool = get("P2C2M.flg.dbgBool", envir=p2c2m.globalVars)
  if (dbgBool) {
    cat("\n", xtermStyle::style("DEBUGMODE> statshelpers.diffrnce", 
        fg="red"), sep="")
  }

  # Note: The difference is "empirical-simulated". Since "empirical" is worse 
  # than "simulated" whenever it does not conform to the coalescent model (i.e.
  # has larger values), significant differences will be more positive than 
  # non-significant differences.
  diff = emp - sim
  # TFL converts from type "list" to type "double"; is important, because
  # is.infinite and is.nan can only work on type "double"
  diff = as.matrix(diff)

  # Removing diff. values that are infinite ("Inf")
  diff = ifelse(is.infinite(diff), NA, diff)
  # Removing diff. values that are "NaN" ("not a number", i.e. 0/0)
  diff = ifelse(is.nan(diff), NA, diff)

  return(diff)
}
