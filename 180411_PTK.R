library(pgscrape)
library(dplyr)

p2p = read.table("Pep2Protein-Short-PTK-86312-86402-swiss_prot-HUMAN-2018-04-11.txt", header = TRUE, sep = "\t")

p2p = p2p %>% distinct( PepProtein_PhosLink, ID, PepProtein_SeqHomology, .keep_all = TRUE)
scrape = p2p %>% group_by(PepProtein_PhosLink, ID, PepProtein_SeqHomology) %>% do(
  {
    aProt = .$PepProtein_UniprotID
    aSite =  paste(.$PepProtein_Residue, .$PepProtein_PhosSite, sep = "")
    print(paste(aProt, aSite))
    res = getPNetPredictions(uni = aProt,
                             ps  = aSite)
    if(is.null(res)){
      res = data.frame()
    }
    res
  }
)
save(file = "180411 scrape.RData", scrape)
