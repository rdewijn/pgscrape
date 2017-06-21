library(pgscrape)
library(dplyr)

p2p = read.table("Pep2Protein-PTK-Full-86312-UniProt-HUMAN-2014-11-05.txt", header = TRUE, sep = "\t")

scrape = p2p %>% group_by(PepProtein_PhosLink, ID, PepProtein_SeqHomology, pident) %>% do(
  {
    aProt = .$PepProtein_UniprotID
    aSite =  paste(.$PepProtein_Residue, .$PepProtein_PhosSite, sep = "")
    print(paste(aProt, aSite))
    res = getPNetPredictions(uni = aProt,
                             ps  = aSite)
    if(is.null(res)){
      return(data.frame())
    }
    return(res)

  }
)
save(file = "170620 scrape.RData", scrape)
