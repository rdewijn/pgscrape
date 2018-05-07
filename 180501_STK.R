library(pgscrape)
library(dplyr)

p2p = read.table("./87102/Pep2Protein-Short-STK-87102-swiss_prot-HUMAN-2018-05-01.txt" , header = TRUE, sep = "\t")
p2p = p2p %>% distinct( PepProtein_PhosLink, ID, PepProtein_SeqHomology, .keep_all = TRUE)
p2p = p2p %>% arrange(PepProtein_PhosLink, ID, PepProtein_SeqHomology)
scrape = p2p %>% group_by(PepProtein_PhosLink, ID, PepProtein_SeqHomology) %>% do(
  {
    aProt = .$PepProtein_UniprotID
    aSite =  paste(.$PepProtein_Residue, .$PepProtein_PhosSite, sep = "")

    res = try(getPNetPredictions(uni = aProt,
                             ps  = aSite))

    if(is.null(res) | class(res) == "try-error"){
      print(paste(aProt, aSite))
      res = data.frame(Kinase_name = NA)
    }
    res
  }
)
save(file = "180503_87102_scrape.RData", scrape)
