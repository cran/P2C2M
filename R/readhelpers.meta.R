readhelpers.meta <-
function(metaData, metaCols, treePart, sTreeFlg) {
  # Descr:  dealing with metadata
  # Deps:   -
  # I/p:    metaData
  #         metaCols
  #         treePart
  #         sTreeFlg

  dbgBool = get("P2C2M.flg.dbgBool", envir=p2c2m.globalVars)
  if (dbgBool) {
    cat("\n", xtermStyle::style("DEBUGMODE> readhelpers.meta", fg="red"), 
        sep="")
  }

  # 1. Setting up empty matrix (called "brData") of dimensions: 
  #    rows = length(meta), columns = metaCols+3 or +2
  if(metaCols==1) {
    brData = array(dim=c(length(metaData), metaCols+2))
  }
  if(sTreeFlg && metaCols==3) {
    brData = array(dim=c(length(metaData), metaCols+3))
  }

  # 2. Save node name, branch length and dmvaphic info to matrix

  for(i in 1:length(metaData)) {
    # brInfo: "br", "length"
    brInfo = unlist(strsplit(treePart[i], ":"))
    metaInfo = unlist(strsplit(metaData[i], ","))

    if(metaCols==1) {
      brData[i,] = c(brInfo, metaInfo)
    }
    if(sTreeFlg && metaCols==3) {
      # Calculation of dmvValue as mean between dmv_start and dmv_end
      dmvValue = mean(as.numeric(c(metaInfo[2], metaInfo[3])))
      # Fill the first two columns of the matrix "brData" with node 
      # name and branch length info,and fill the remaining three columns 
      # with the metadata info
      brData[i,] = c(brInfo, metaInfo, dmvValue)
    }
  }

  # 3. Remove all those brData entries that do not contain 
  # a branch name (i.e. that represent nodes)
  brData = brData[which(brData[,1]!=""),]

  # 4. Naming brData
  # name the columns of matrix "brData" depending on 
  # the number of metadata infos
  if(sTreeFlg && metaCols==1) {
    colnames(brData) = c("br", "length", "dmv")
  }

  if(!sTreeFlg && metaCols==1) {
    colnames(brData) = c("br", "length", "rate")
  }

  if(sTreeFlg && metaCols==3) {
    colnames(brData) = c("br", "length", "dmt", "dmv_start", "dmv_end", "dmv")
  }
  # name the rows of of matrix "brData" by the node name
  rownames(brData) = brData[,1]

  return(brData)
}
