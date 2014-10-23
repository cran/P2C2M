stats.main <-
function(ind, loci, alpha) {
  # Descr:  generate significance tables
  # Deps:   stats.outmaker
  # I/p:    ind
  #         loci
  #         alpha
  # Note:   CV = coefficient of variance
  #         avgs = arithmetic means
  #         qntls = quantiles

  dbgBool = get("P2C2M.flg.dbgBool", envir=p2c2m.globalVars)
  if (dbgBool) {
    cat("\n", xtermStyle::style("DEBUGMODE> stats.main", fg="red"), "\n", sep="")
  }

##############################
# 1. Setting number of tails #
##############################

  # T-tests for all descriptive statistics are two-tailed, because there is no
  # a priori reason in which direction they should differ.
  tailL = list()
  tailL$GTP = "2"
  tailL$RAY = "2"
  tailL$NDC = "2"
  tailL$GSI = "2"
  ## T-tests for NDC should be left one-tailed, because trees not compliant 
  ## with the coalescent model have a higher number of deep coalescences.
  #NDCtail = "1l"
  ## T-tests for GSI should be right one-tailed, because trees not compliant 
  ## with the coalescent model have lower values.
  #GSItail = "1r"

#############################
# 2. Inferring significance #
#############################

  stats = get("P2C2M.flg.descrStats", envir=p2c2m.globalVars)
  perGene = acrGenes = list()
  for (stat in stats) {
    perGene[[stat]] = stats.perGene(ind[[stat]]$dif, alpha, tailL[[stat]])
    acrGenes[[stat]] = stats.acrossGenes(ind[[stat]]$dif, alpha, tailL[[stat]])
  }

  #RAY_perG = stats.perGene(ind$RAY$dif, alpha, tail=RAYtail)
  #NDC_perG = stats.perGene(ind$NDC$dif, alpha, tail=NDCtail)
  #GSI_perG = stats.perGene(ind$GSI$dif, alpha, tail=GSItail)

  #RAY_acrossG = stats.acrossGenes(ind$RAY$dif, alpha, tail=RAYtail)
  #NDC_acrossG = stats.acrossGenes(ind$NDC$dif, alpha, tail=NDCtail)
  #GSI_acrossG = stats.acrossGenes(ind$GSI$dif, alpha, tail=GSItail)

#######################
# 3. Combining output #
#######################
  outL = list()

  # perGene output
  outL$perGene = sapply(perGene, cbind)
  #outL$perGene = cbind(GTP_perG, RAY_perG, NDC_perG, GSI_perG)
  rownames(outL$perGene) = c(loci)

  # acrossGene output
  outL$acrGenes = sapply(acrGenes, cbind)
  #outL$acrossGenes = cbind(GTP_acrossG, RAY_acrossG, NDC_acrossG, GSI_acrossG)
  rownames(outL$acrGenes) = c("Sum", "Mean", "Median", "Mode", "CV")

  # Naming rows of output
  names = c()
  for (stat in stats) {
    names = c(names, paste(stat, "[", tailL[[stat]], "]", sep=""))
  }
  colnames(outL$perGene) = colnames(outL$acrGenes) = names

  return(outL)
}
