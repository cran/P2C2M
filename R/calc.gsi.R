calc.gsi <-
function(gTree, assoc, singleAllele) {
  # Descr:  returns the Genealogical Sorting Index for a given tree
  # Deps:   -
  # I/p:    gTree
  #         assoc
  #         singleAllele
  # Note:   gsi = "genealogical sorting index"

  #dbgBool = get("P2C2M.flg.dbgBool", envir=p2c2m.globalVars)
  #if (dbgBool) {
  #  cat("\n", xtermStyle::style("DEBUGMODE> calc.gsi", fg="red"), sep="")
  #}

  # matching the correct associations
  speciesAssoc = assoc[,1][match(gTree$tip.label, assoc[,2])]

## 1. Generate list of all species (i.e. unique elements of "speciesAssoc")
  species = unique(speciesAssoc)
  # Remove single allele species
  species = suppressWarnings(species[species != singleAllele])

## 2. Enclosing data in a list
  data = list()
  data$species = species
  data$gTree = gTree
  data$speciesAssoc = speciesAssoc

### 3. DEBUG mode
#    logdata = list(list(GSI=data))
#    names(logdata) = get("p2c2m.flg.repID", envir=p2c2m.globalVars)
#    loghelpers.dbg(logdata, "DescrStatsRawInput", "INPUT OF 'GSI'")

## 4. Calculating descriptive statistic
  # Calculating the gsi for every "group" (i.e. species other than 
  # single allele species)
  outL = c()
  for (sp in data$species) {
    outL = c(outL, genealogicalSorting::gsi(data$gTree, sp, data$speciesAssoc))
  }

  # Ensuring that the gsi values together add up to a number between 0 and 1
  outd = lapply(outL, function(x){x/length(outL)})
  outd = sum(unlist(outd))

  return(outd)
}
