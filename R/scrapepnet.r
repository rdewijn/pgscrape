#' \code{getPNetPredictions} returns a data frame with predicted upstream kinases obtained from the phosphoNET website
#' @param urlBase url base string for phospoNET html query. Default: "http://www.phosphonet.ca/kinasepredictor.aspx"
#' @param uni a uniprot ID
#' @param ps phs]ophorylated residue
#' @return a data frame with the phosphoNET prediction table. Returns NULL if the uni ps combination is not found.
#' @examples
#' # For obtaining predictions for Tyrosine 486 form protein with uniprot ID Q14247
#' getPNetPredictions(uni = "Q14247", ps = "Y486")
#' @import XML dplyr
#' @export
getPNetPredictions = function(urlBase = "http://www.phosphonet.ca/kinasepredictor.aspx", uni,ps){
  url = paste(urlBase,"?uni=",uni,"&ps=",ps, sep = "")
  aDoc = htmlTreeParse(url, useInternalNode = TRUE)

  PredictedKinaseXPath = "//tr[td[@class='td-SideHeader td-RowHeader' and contains(., 'Kinase')] ]"
  rawTable = as.data.frame(t(xpathSApply(aDoc, PredictedKinaseXPath, getChildrenStrings)))
  if( dim(rawTable)[2] > 0){
    predTable = data.frame(Kinase_rank_str = rawTable[,1],
                         Kinase_rank = 1:dim(rawTable)[1],
                         Kinase_name = rawTable[,2], Kinase_UniprotID = rawTable[,3],
                         Kinase_PredictorVersion2Score = rawTable[,4])
    #also get the p-site sequence
    pSiteNode = getNodeSet(aDoc, "//tr[contains(., 'P-Site Sequence')]")[[1]]
    aMatch = regexpr("P-Site Sequence:(?<sequence>[[:upper:]]+)\r", xmlValue(pSiteNode), perl = TRUE)
    ini = attr(aMatch, "capture.start")
    fini = ini + attr(aMatch, "capture.length")
    pSeq = substr(xmlValue(pSiteNode), ini, fini-1)
    predTable = data.frame(predTable, phosphoNET_PSiteSequence = pSeq,
                                      Substrate_UniprotID = uni,
                                      Substrate_pSite = ps,
                                      Substrate_pSitePosition = as.numeric(substring(ps,2)),
                                      Substrate_phosphoLink = paste(uni, substring(ps,2), sep = "_")
                           )

  } else {
    predTable = NULL
  }
  return(predTable)
}
#' \code{getPNetPredictions} returns a data frame with predicted upstream kinases obtained from the phosphoNET website for multiple phosphosites
#' @param urlBase url base string for phospoNET html query. Default: "http://www.phosphonet.ca/kinasepredictor.aspx"
#' @param uni a vectyor containing uniprot IDs
#' @param ps a vector containing phosphorylated residues, same length as uni
#' @return a data frame with the phosphoNET predictions for the uni ps combination that were found found.
#' @examples
#' # For obtaining predictions for Tyrosine 486, Threonine 24, ans Serine 11 of protein with uniprot ID Q14247
#' getPNetPredictionTable(uni = c("Q14247", "Q14247","Q14247"), ps = c("Y486", "T24","S11"))
#'@export
getPNetPredictionTable = function(urlBase = "http://www.phosphonet.ca/kinasepredictor.aspx", uni, ps){
  if (length(uni) != length(ps)) stop ("uni and ps must be vectors of the same length")
  qdf = data.frame(uni = uni, ps = ps)
  result = qdf %>% group_by(1:nrow(qdf)) %>% do({
    aPred = getPNetPredictions(urlBase, .$uni, .$ps)
  })
  return(result[,-1])
}


