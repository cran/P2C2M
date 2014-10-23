readhelpers.extract <-
function(tree, locus, sTreeFlg) {
  # Descr:  extracting tree and branch rate vectors, 
  #         and order vectors appropriately
  # Deps:   -
  # I/p:    tree
  #         locus
  #         sTreeFlg

  dbgBool = get("P2C2M.flg.dbgBool", envir=p2c2m.globalVars)
  if (dbgBool) {
    cat("\n", xtermStyle::style("DEBUGMODE> readhelpers.extract", fg="red"),
        sep="")
  }

  # 1. Adding internal node labels
  # Count the number of taxa by counting number of square brackets, 
  # b/c each taxon is followed by metadata
  if (sTreeFlg==F) {
    ntax = (length(unlist(strsplit(tree, "\\[")))+1)/2
  }
  if (sTreeFlg==T) {
    ntax = length(unlist(strsplit(tree, "\\[")))/2
  }

  # Calculate number of nodes
  if (sTreeFlg==F) {
    nodes = (ntax+2):((2*ntax)-1)
  }
  if (sTreeFlg==T) {
    nodes = (ntax+1):((2*ntax)-1)
  }

  # Add internal node labels
  for (node in nodes) {
    if (sTreeFlg==F) {
      repl = paste(")", node, ":[", sep="")
      tree = sub("):\\[", repl, tree)
    }

    if (sTreeFlg==T) {
      repl = paste(")", node, "[", sep="")
      tree = sub(")\\[", repl, tree)
    }
  }

  # 2. Extract all lines with "dmv/dmt" | "rate" info
  # Split tree spec at every square bracket 
  wntdList = strsplit(tree, "\\[|\\]")
  # Extract those list elements that contain keyword "dmv/dmt" | "rate"
  if (sTreeFlg==F) {wntdLines = grep("rate", unlist(wntdList))}
  if (sTreeFlg==T) {wntdLines = grep("dmv", unlist(wntdList))}
  metaRaw = unlist(wntdList)[wntdLines]

  # Count seperate pieces of info (i.e. columns) 
  # related to "dmv/dmt" | "rate" values
  metaCols = length(unlist(strsplit(metaRaw[1], ",")))
  # Parse "dmv/dmt" | "rate" values as raw numbers
  if (sTreeFlg==F) {metaData = gsub("&|rate=|\\{|\\}", "", metaRaw)}
  if (sTreeFlg==T) {metaData = gsub("&|dm.=|\\{|\\}", "", metaRaw)}

  # Necessary for starBeast.v.1.8. and higher: Remove locus name plus trailing 
  # period from metadata element
  metaData = gsub(paste(locus, ".", sep=""), "", metaData)

  # 3. Extract node names and ultrametric branch length info
  # Remove any metadata from tree spec
  treePart = gsub("\\[[^]]*\\]", "\\[\\]", tree)
  treePart = unlist(strsplit(treePart, ",|)"))
  treePart = gsub("\\(|\\)|;|\\[|\\]", "", treePart)
  if (sTreeFlg==T) {
    treePart[length(treePart)] = paste(treePart[length(treePart)], NA, sep=":")
  }

  # 4. Metadata extraction
  brData = readhelpers.meta(metaData, metaCols, treePart, sTreeFlg)

  # 5. Generate a transl. matrix btw. tree and "brData"
  # Remove any metadata from tree spec
  aString = gsub("\\[[^]]*\\]", "", tree)

  # Convert tree from string to APE object
  aTree = ape::read.tree(text=aString)

  if (sTreeFlg==F) {
    translate = cbind(aTree$node.label[-1], (ntax+2):nodes[length(nodes)])
  }
  if (sTreeFlg==T) {
    translate = cbind(aTree$node.label, (ntax+1):nodes[length(nodes)])
  }

  # Generate preliminary tree spec matrix with placeholders in second column
  translate = rbind(translate, cbind(aTree$tip.label, 1:ntax))

  # Name the rows by the second column
  rownames(translate) = translate[,2]
  # Extract items in "translate" that are named as in "aTree$edge[,2]" 
  # and save them to "ordTransl". This step orders the items in 
  # "translate" to the order specified in "aTree".
  ordTransl = translate[as.character(aTree$edge[,2]),]
  # Remove brData entries that do not contain 
  # a branch name (i.e. that represent nodes)
  if (sTreeFlg==T) {
    ordTransl = ordTransl[which(ordTransl[,1]!=""),]
  }

  # Obtain these entries in "brData" that are labelled as given in 
  # "ordTransl[,1]". This step obtains a subset of entries from 
  # brData and saves them as ordbrData
  outd = brData[ordTransl[,1],]

  # Add the flast entry of the matrix "brData" to the matrix "final_brData"
  if (sTreeFlg==T) {
    outd = rbind(outd, brData[length(brData[,2]),])
  }
  # remove all rownames from "final_brData"
  rownames(outd) = NULL

  return(outd)
}
