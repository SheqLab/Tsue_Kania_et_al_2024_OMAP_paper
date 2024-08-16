#make sure you have downloaded packages installed: tidyverse, pheatmap, cluster, factoextra
#at CoreQuant, export with renamed column and conditions
#install command: install.packages("tidyverse")
#going to use some fo devin's code from ODB_rc_pq_040921.R

#import packages
library(tidyverse)
library(pheatmap)
library(cluster)
library(factoextra)
library(reshape2)
library(corrplot)
library(readxl)

#wd working directory
setwd("~/Documents/UW/Shechner Lab/Biotinylated Oligo exp/Mass Spec/MS analysis/ODB paper 2021 analysis")

#read data into tibble
odb_exp5 = read_tsv("protein_quant_35.tsv")
uniqueNucleolar = read_xlsx("unique_nucleolar_genes_listOnly.xlsx")
uniqueNucleolarList <- uniqueNucleolar %>% select("Uniprot")
uniqueNuclear = read_xlsx("unique_nuclear_genes_listOnly.xlsx")
uniqueNuclearList <-uniqueNuclear %>% select("Uniprot")

#filter contaminants and deal missing values?
pb.scaled = odb_exp5 %>%
  filter(!grepl("##",`Protein Id`) & !grepl("contaminant",`Protein Id`)) %>%
  select(contains("Protein Id") | contains("scaled"))

#make correlation matrices
#make correlation heatmap of all replicates except blank channel
scaledOnly <- pb.scaled %>% select(contains("scaled") & !contains("blank"))
corrAll = cor(scaledOnly)

pheatmap(corrAll, name = "Exp 5 correlation matrix")



#fast test correlation plots all work with library corrplot and reshape2
corrplot(corrAll)
melt_corrAll = melt(corrAll)
ggplot(data = melt_corrAll, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()

#heat maps for entire exp 5
pheatmap(pb.scaled %>% select(contains("scaled")), 
         name = "ODB", #title of legend
         column_title = "Scaled", row_title = "Condition",
         cluster_cols = FALSE, labels_row = pb.scaled$`Protein Id`,
)

pheatmap(pb.scaled %>% select(`its1_1m_1 scaled`,`its1_1m_2 scaled`,`its1_1m_3 scaled`,
                              `7sk-1m_1 scaled`,`7sk-1m_2 scaled`,`7sk-1m_3 scaled`, 
                              '7sk-10m_1 scaled', '7sk-10m_2 scaled', '7sk-10m_3 scaled'), 
         name = "ODB", #title of legend
         column_title = "Scaled", row_title = "Condition",
         cluster_cols = FALSE, labels_row = pb.scaled$`Protein Id`
)

#mean, t test and sd
pb.sm = pb.scaled %>%
  rownames_to_column(var = "Protein") %>%
  rowwise() %>%
  mutate(
    mean.scr_1m = mean(c_across("scrmbl-1m_1 scaled":"scrmbl-1m_3 scaled")),
    mean.scr_10m = mean(c_across("scrmbl-10m_1 scaled":"scrmbl-10m_3 scaled")),
    mean.7sk_1m = mean(c_across("7sk-1m_1 scaled":"7sk-1m_3 scaled")),
    mean.7sk_10m = mean(c_across("7sk-10m_1 scaled":"7sk-10m_3 scaled")),
    mean.its1_1m = mean(c_across("its1_1m_1 scaled":"its1_1m_3 scaled")),
    sd.scr_1m = sd(c_across("scrmbl-1m_1 scaled":"scrmbl-1m_3 scaled")),
    sd.scr_10m = sd(c_across("scrmbl-10m_1 scaled":"scrmbl-10m_3 scaled")),
    sd.7sk_1m = sd(c_across("7sk-1m_1 scaled":"7sk-1m_3 scaled")),
    sd.7sk_10m = sd(c_across("7sk-10m_1 scaled":"7sk-10m_3 scaled")),
    sd.its1_1m = sd(c_across("its1_1m_1 scaled":"its1_1m_3 scaled")),
    t.7sk_1m = t.test(c_across("7sk-1m_1 scaled":"7sk-1m_3 scaled"),
                      c_across("scrmbl-10m_1 scaled":"scrmbl-10m_3 scaled"))$p.value,
    t.7sk_1mS1 = t.test(c_across("7sk-1m_1 scaled":"7sk-1m_3 scaled"),
                        c_across("scrmbl-1m_1 scaled":"scrmbl-1m_3 scaled"))$p.value,
    t.7sk_10m = t.test(c_across("7sk-10m_1 scaled":"7sk-10m_3 scaled"),
                       c_across("scrmbl-10m_1 scaled":"scrmbl-10m_3 scaled"))$p.value,
    t.its1_1m = t.test(c_across("its1_1m_1 scaled":"its1_1m_3 scaled"),
                       c_across("scrmbl-10m_1 scaled":"scrmbl-10m_3 scaled"))$p.value,
    t.its1_1mS1 = t.test(c_across("its1_1m_1 scaled":"its1_1m_3 scaled"),
                         c_across("scrmbl-1m_1 scaled":"scrmbl-1m_3 scaled"))$p.value,
    t.itsVS7sk1m = t.test(c_across("its1_1m_1 scaled":"its1_1m_3 scaled"), 
                          c_across("7sk-1m_1 scaled":"7sk-1m_3 scaled"))$p.value,
    t.7sk10mVS7sk1m = t.test(c_across("7sk-10m_1 scaled":"7sk-10m_3 scaled"),
                             c_across("7sk-1m_1 scaled":"7sk-1m_3 scaled"))$p.value,
    l2r.toS.71m = log(mean.7sk_1m/mean.scr_10m,2),
    l2r.toS.71mS1 = log(mean.7sk_1m/mean.scr_1m,2),
    l2r.toS.710m = log(mean.7sk_10m/mean.scr_10m,2),
    l2r.toS.i1m = log(mean.its1_1m/mean.scr_10m,2),
    l2r.toS.i1mS1 = log(mean.its1_1m/mean.scr_1m,2),
    l2r.toS.itsVS7sk = log(mean.its1_1m/mean.7sk_1m,2),
    l2r.toS.7sk10mVS7sk1m = log(mean.7sk_10m/mean.7sk_1m,2)
  )
#p adjust with BH
pb.sm = pb.sm %>%
  add_column(
    q.7sk1m = p.adjust(pb.sm$t.7sk_1m,method = "BH"),
    q.7sk10m = p.adjust(pb.sm$t.7sk_10m,method = "BH"),
    q.its1m = p.adjust(pb.sm$t.its1_1m,method = "BH"),
    q.itsVS7sk = p.adjust(pb.sm$t.itsVS7sk1m, method = "BH"),
    q.7sk10mVS7sk1m = p.adjust(pb.sm$t.7sk10mVS7sk1m, method = "BH"),
    q.7sk1mS1 = p.adjust(pb.sm$t.7sk_1mS1, method = "BH"),
    q.its1mS1 = p.adjust(pb.sm$t.its1_1mS1, method = "BH")
  )

#add column for log pvalue
pb.sm = pb.sm %>%
  add_column(
    log10BHpval.7sk_1m = -1*log10(pb.sm$q.7sk1m),
    log10BHpval.7sk_10m = -1*log10(pb.sm$q.7sk10m),
    log10BHpval.its_1m = -1*log10(pb.sm$q.its1m), 
    log10BHpval.itsVS7sk = -1*log10(pb.sm$q.itsVS7sk),
    log10BHpval.7sk10mVS7sk1m = -1*log10(pb.sm$q.7sk10mVS7sk1m),
    log10BHpval.7sk_1mS1 = -1*log10(pb.sm$q.7sk1mS1),
    log10BHpval.its_1mS1 = -1*log10(pb.sm$q.its1mS1),
    log10pval.7sk_1m = -1*log10(pb.sm$t.7sk_1m), 
    log10pval.7sk_10m = -1*log10(pb.sm$t.7sk_10m),
    log10pval.its_1m = -1*log10(pb.sm$t.its1_1m),
    log10pval.itsVS7sk = -1*log10(pb.sm$t.itsVS7sk1m),
    log10pval.7sk10mVS7sk1m = -1*log10(pb.sm$t.7sk10mVS7sk1m),
    log10pval.7sk_1mS1 = -1*log10(pb.sm$t.7sk_1mS1),
    log10pval.its_1mS1 = -1*log10(pb.sm$t.its1_1mS1)
  )

pb.plc.scaled = pb.sm %>%
  select(contains("scaled") | contains("Protein")) %>%
  filter_if(~is.numeric(.), all_vars(!is.infinite(.))) %>%
  pivot_longer(contains("scaled"), names_to = "Condition", values_to = "Scaled Ratio") %>%
  mutate(simpleCond = unlist(
    lapply(Condition, FUN = function(x){ return(unlist(str_split(x,'\\_'))[1]) }))
  )


#collapse all of the l2r into one tibble
pb.plc.lr = pb.sm %>%
  select(contains("l2r") | contains("Protein")) %>%
  filter_if(~is.numeric(.), all_vars(!is.infinite(.))) %>%
  pivot_longer(contains("l2r"), names_to = "Condition", values_to = "Log Ratio")

#match lists add +1 and -1 if in HPA nucleolar or nuclear protein list respectively
#make a uniprot accession column in all data
allProt = pb.sm %>% add_column(Uniprot.wIso = unlist(
  lapply(pb.sm$`Protein Id`, FUN = function(x){ return(unlist(str_split(x,'\\|'))[2]) }))
) %>%
  mutate(Uniprot = unlist(
    lapply(Uniprot.wIso, FUN = function(x){ return(unlist(str_split(x,'\\-'))[1]) }))
  )
#match the unique nucleolar and unique nuclear protein lists from HPA to our data
nucleolarProt_unique = semi_join(allProt, uniqueNucleolarList)
nucleolarProt_unique = nucleolarProt_unique %>% add_column(1)
nuclearProt_unique = semi_join(allProt, uniqueNuclearList)
nuclearProt_unique = nuclearProt_unique %>% add_column(-1)


#plot histogram for l2r ITS vs Scram and percentage within these lists
nucleolarProt_uniqueL <- nucleolarProt_unique %>% select("Protein Id", "l2r.toS.i1mS1")
nucleolarProt_unique %>% ggplot(aes(x=l2r.toS.i1mS1)) + geom_histogram() +
  xlab("log2 FC (ITS/scram)") + ylab("Count") + ggtitle("ITS1 min vs Scram 1 min Nucleolar List")
nuclearProt_unique %>% ggplot(aes(x=l2r.toS.i1mS1)) + geom_histogram() +
  xlab("log2 FC (ITS/scram)") + ylab("Count") + ggtitle("ITS1 min vs Scram 1 min Nuclear List")
nuclearProt_unique %>% ggplot(aes(x=l2r.toS.71mS1)) + geom_histogram() +
  xlab("log2 FC (7SK/scram)") + ylab("Count") + ggtitle("7SK1 min vs Scram 1 min Nuclear List")
nucleolarProt_unique %>% ggplot(aes(x=l2r.toS.71mS1)) + geom_histogram() +
  xlab("log2 FC (7SK/scram)") + ylab("Count") + ggtitle("7SK1 min vs Scram 1 min Nucleolar List")
#try to plot histograms for l2r 7sk vs scram


#generate extent of correlation plots of bio reps 
#try l2r of experimental over negative control for each replicate. pause.
#test = pb.sm %>%
#  add_column(
#    l2rits_rep1 = log(pb.sm$`its1_1m_1 scaled`/pb.sm$`scrmbl-1m_1 scaled`,2)
#  )
#  mutate(l2rits_rep1 = log(`its1_1m_1 scaled`/`scrmbl-1m_1 scaled`,2))


#  l2rits_rep1 = log(pb.sm$`its1_1m_1 scaled`/pb.sm$`scrmbl-1m_1 scaled`,2)
#  mutate(
#    l2rits_rep1 = log(pb.sm$`its1_1m_1 scaled`/pb.sm$`scrmbl-1m_1 scaled`,2)
#  ) 


#simple volcano plot transform data first so easy to plot from tibble
vcITS1m <- pb.sm %>% select("Protein Id", "l2r.toS.i1m", "log10BHpval.its_1m", "log10pval.its_1m")
vcITS1mS1 <- pb.sm %>% select("Protein Id", "l2r.toS.i1mS1", "log10BHpval.its_1mS1", "log10pval.its_1mS1")
vc7sk1m <- pb.sm %>% select("Protein Id", "l2r.toS.71m", "log10BHpval.7sk_1m", "log10pval.7sk_1m")
vc7sk1mS1 <- pb.sm %>% select("Protein Id", "l2r.toS.71mS1", "log10BHpval.7sk_1mS1", "log10pval.7sk_1mS1")
vc7sk10m <- pb.sm %>% select("Protein Id", "l2r.toS.710m", "log10BHpval.7sk_10m", "log10pval.7sk_10m")
vc7skVSits <- pb.sm %>% select("Protein Id", "l2r.toS.itsVS7sk", "log10BHpval.itsVS7sk", "log10pval.itsVS7sk")
vc7sk10mVS7sk1m <- pb.sm %>% select("Protein Id", "l2r.toS.7sk10mVS7sk1m", "log10BHpval.7sk10mVS7sk1m", "log10pval.7sk10mVS7sk1m")

#plot BH adjust pval and pval 
vcITS1m %>% ggplot(aes(l2r.toS.i1m,log10BHpval.its_1m)) + geom_point() +
  xlab("log2 FC (ITS/scram)") + ylab("-log10 adjust p-value")  + 
  geom_hline(yintercept = 2, linetype = 2, alpha = 0.5) + ggtitle("ITS1 min vs Scram 10 min BH adjust")
vcITS1mS1 %>% ggplot(aes(l2r.toS.i1mS1, log10BHpval.its_1mS1)) + geom_point() +
  geom_hline(yintercept = 2, linetype = 2, alpha = 0.5) + 
  xlab("log2 FC (ITS/scram 1 min)") + ylab("-log10 adjust pval") + ggtitle("ITS1 min vs Scram 1 min BH adjust")
vc7sk1mS1 %>% ggplot(aes(l2r.toS.71mS1, log10BHpval.7sk_1mS1)) + geom_point() +
  xlab("log2 FC (7sk/scram)") + ylab("-log10 adjust p-value")  + 
  geom_hline(yintercept = 2, linetype = 2, alpha = 0.5) + ggtitle("7SK 1 min vs Scram 1 min BH adjust")

#ITS vs scram not adjusted pval
vcITS1m %>% ggplot(aes(l2r.toS.i1m, log10pval.its_1m)) + geom_point() +
  xlab("log2 FC (ITS/scram)") + ylab("-log10 p-value")  + 
  geom_hline(yintercept = 2, linetype = 2, alpha = 0.5) + ggtitle("ITS1 min vs Scram 10 min")
vcITS1mS1 %>% ggplot(aes(l2r.toS.i1mS1, log10pval.its_1mS1)) + geom_point() +
  geom_hline(yintercept = 2, linetype = 2, alpha = 0.5) + 
  xlab("log2 FC (ITS/scram 1 min)") + ylab("-log10 pval") + ggtitle("ITS1 min vs Scram 1 min")

vc7sk1m %>% ggplot(aes(l2r.toS.71m,log10BHpval.7sk_1m)) + geom_point()
vc7sk1m %>% ggplot(aes(l2r.toS.71m, log10pval.7sk_1m)) + geom_point()

vc7sk10m %>% ggplot(aes(l2r.toS.710m,log10BHpval.7sk_10m)) + geom_point()
vc7sk10m %>% ggplot(aes(l2r.toS.710m, log10pval.7sk_10m)) + geom_point()

vc7skVSits %>% ggplot(aes(l2r.toS.itsVS7sk,log10BHpval.itsVS7sk)) + geom_point() +
  xlab("log2 FC (ITS/7SK)") + ylab("-log10 adjusted p-value") + 
  geom_hline(yintercept = 1.3, linetype = 2.5, alpha = 0.5)
vc7skVSits %>% ggplot(aes(l2r.toS.itsVS7sk,log10pval.itsVS7sk)) + geom_point() +
  xlab("log2 FC (ITS/7SK)") + ylab("-log10 p-value")  + 
  geom_hline(yintercept = 2, linetype = 2, alpha = 0.5)

vc7sk10mVS7sk1m %>% ggplot(aes(l2r.toS.7sk10mVS7sk1m,log10BHpval.7sk10mVS7sk1m)) + geom_point()
vc7sk10mVS7sk1m %>% ggplot(aes(l2r.toS.7sk10mVS7sk1m,log10pval.7sk10mVS7sk1m)) + geom_point()

#generate lists of proteins based on thresholding for pval 0.01
nucProt <- vcITS1mS1 %>%
  # Filter for significant observations 
  filter(log10BHpval.its_1mS1 >= 2 & (l2r.toS.i1mS1 >= 2.5 | l2r.toS.i1mS1 <= -2.5)) 

##GO term analysis
#for all nucleolar proteins in cluster 1,4,7,12
#change directory
test = read.delim("./All Nucleolar Proteins/allNucProt_DAVID_indProteinAssignemnts_CC.txt")

#heatmaps

#looking at 7sk conditions
pheatmap(pb.scaled %>% select(`its1_1m_1 scaled`,`its1_1m_2 scaled`,`its1_1m_3 scaled`,
                               `7sk-1m_1 scaled`,`7sk-1m_2 scaled`,`7sk-1m_3 scaled`, '7sk-10m_1 scaled',
                              '7sk-10m_2 scaled', '7sk-10m_3 scaled'), 
                              name = "ODB", #title of legend
                              column_title = "Scaled", row_title = "Condition",
                              cluster_cols = FALSE, labels_row = pb.scaled$`Protein Id`,
                              clustering_method = "median"
)
#test that
pheatmap(pb.scaled %>% select(`its1_1m_1 scaled`,`its1_1m_2 scaled`,`its1_1m_3 scaled`,
                              `7sk-1m_1 scaled`,`7sk-1m_2 scaled`, '7sk-10m_1 scaled',`7sk-1m_3 scaled`,
                              '7sk-10m_2 scaled', '7sk-10m_3 scaled'), 
         name = "ODB", #title of legend
         column_title = "Scaled", row_title = "Condition",
         cluster_cols = FALSE, labels_row = pb.scaled$`Protein Id`,
         clustering_method = "median"
)


#look at 7sk known protein interactors within exp 35

ggplot(pb.plc.scaled %>% filter(grepl("LARP7",`Protein Id`) |
                                  grepl("MEPCE",`Protein Id`)|
                                  grepl("HEXI1",`Protein Id`)|
                                  grepl("CDK9",`Protein Id`)|
                                  grepl("CCNT1",`Protein Id`)),
       aes(factor(simpleCond,levels = c("blank scaled","scrmbl-1m","scrmbl-10m","7sk-1m","7sk-10m","its1")),
           `Scaled Ratio`,fill = simpleCond)) +
  geom_boxplot() +
  facet_wrap(~`Protein Id`,nrow = 1,scales = "free") +
  theme_classic()

ggplot(pb.plc.scaled %>% filter(grepl("LARP7",`Protein Id`) |
                                  grepl("MEPCE",`Protein Id`)|
                                  grepl("HEXI1",`Protein Id`)|
                                  grepl("CCNT1",`Protein Id`)|
                                  grepl("CDK9",`Protein Id`)),
       aes(factor(simpleCond,levels = c("blank scaled","scrmbl-1m","scrmbl-10m","7sk-1m","7sk-10m","its1")),
           `Scaled Ratio`,fill = Condition)) +
  geom_bar(stat="identity",position="dodge") +
  facet_wrap(~`Protein Id`,nrow = 1,scales = "free") +
  theme_classic()



#looking at exp 6 dataset of ITS time titration
odb = read_tsv("protein_quant_48.tsv")

scaledExp6 = odb %>%
  filter(!grepl("##",`Protein Id`) & !grepl("contaminant",`Protein Id`)) %>%
  select(contains("Protein Id") | contains("scaled"))

scaled.meanExp6 = scaledExp6 %>%
  rownames_to_column(var = "Protein") %>%
  rowwise() %>%
  mutate(
    mean.7sk_1m = mean(c_across("7sk_ctrl_1 scaled":"7sk_ctrl_3 scaled")),
    mean.its1_1s = mean(c_across("its1_1s_1 scaled":"its1_1s_3 scaled")),
    mean.its1_1m = mean(c_across("its1_1min_1 scaled":"its1_1min_3 scaled")),
    mean.its1_10m = mean(c_across("its1_10min_1 scaled":"its1_10min_3 scaled")),
    mean.its1_100m = mean(c_across("its1_100min_1 scaled":"its1_100min_3 scaled")),
    t.its1_1s = t.test(c_across("its1_1s_1 scaled":"its1_1s_3 scaled"),
                       c_across("7sk_ctrl_1 scaled":"7sk_ctrl_3 scaled"))$p.value,
    t.its1_1m = t.test(c_across("its1_1min_1 scaled":"its1_1min_3 scaled"),
                       c_across("7sk_ctrl_1 scaled":"7sk_ctrl_3 scaled"))$p.value,
    t.its1_10m = t.test(c_across("its1_10min_1 scaled":"its1_10min_3 scaled"),
                        c_across("7sk_ctrl_1 scaled":"7sk_ctrl_3 scaled"))$p.value,
    t.its1_100m = t.test(c_across("its1_100min_1 scaled":"its1_100min_3 scaled"),
                         c_across("7sk_ctrl_1 scaled":"7sk_ctrl_3 scaled"))$p.value,
    l2r.to7.1s = log(mean.its1_1s/mean.7sk_1m,2),
    l2r.to7.1m = log(mean.its1_1m/mean.7sk_1m,2),
    l2r.to7.10m = log(mean.its1_10m/mean.7sk_1m,2),
    l2r.to7.100m = log(mean.its1_100m/mean.7sk_1m,2)
  )

scaled.meanExp6 = scaled.meanExp6 %>%
  add_column(
    q.1s = p.adjust(scaled.meanExp6$t.its1_1s,method = "BH"),
    q.1m = p.adjust(scaled.meanExp6$t.its1_1m,method = "BH"),
    q.10m = p.adjust(scaled.meanExp6$t.its1_10m,method = "BH"),
    q.100m = p.adjust(scaled.meanExp6$t.its1_100m,method = "BH")
  )
#look at 7sk interactors in this data set> not done. hav to fix the columns fo plc scaled so that has time point in it
pb.plc.scaled.itsLabel = scaled.meanExp6 %>%
  select(contains("scaled") | contains("Protein")) %>%
  filter_if(~is.numeric(.), all_vars(!is.infinite(.))) %>%
  pivot_longer(contains("scaled"), names_to = "Condition", values_to = "Scaled Ratio") %>%
  mutate(simpleCond = unlist(
    lapply(Condition, FUN = function(x){ return(unlist(str_split(x,'\\_'))[1]) }))
  )
ggplot(pb.plc.scaled.itsLabel %>% filter(grepl("LARP7",`Protein Id`) |
                                  grepl("MEPCE",`Protein Id`)|
                                  grepl("HEXI1",`Protein Id`)|
                                  grepl("CDK9",`Protein Id`)|
                                  grepl("CCNT1",`Protein Id`)),
       aes(factor(simpleCond,levels = c("blank scaled","its1_1s","its1_1m","its1_10m","its1_100m","7sk_ctrl")),
           `Scaled Ratio`,fill = simpleCond)) +
  geom_boxplot() +
  facet_wrap(~`Protein Id`,nrow = 1,scales = "free") +
  theme_classic()

ggplot(pb.plc.scaled.itsLabel %>% filter(grepl("LARP7",`Protein Id`) |
                                  grepl("MEPCE",`Protein Id`)|
                                  grepl("HEXI1",`Protein Id`)|
                                  grepl("CCNT1",`Protein Id`)|
                                  grepl("CDK9",`Protein Id`)),
       aes(factor(simpleCond,levels = c("blank scaled","7sk","its1")),
           `Scaled Ratio`,fill = Condition)) +
  geom_bar(stat="identity",position="dodge") +
  facet_wrap(~`Protein Id`,nrow = 1,scales = "free") +
  theme_classic()

#kmedoid clustering optimization
fviz_nbclust(scaled.meanExp6 %>%
               select(contains("mean")), pam, method = "wss")
#wss method: for total within sum of square
k3_kmedoid = pam(scaled.meanExp6 %>%
                   select(contains("mean")),k=3,nstart=25)
k12_kmedoid = pam(scaled.meanExp6 %>%
                    select(contains("mean")), k=12, nstart=25)

fviz_cluster(k12_kmedoid, data = scaled.mean %>%
               select(contains("mean")),
             geom = "point")
pre.line.kmedcluster = scaled.mean %>%
  rownames_to_column(var = "Protein") %>%
  select(contains("mean.") | contains("Protein")) %>%
  add_column(clust = k7_kmedoid$cluster) %>%
  pivot_longer(contains("mean."), names_to = "Condition", values_to = "Scaled")

ggplot(pre.line.kmedcluster, aes(x = Condition, y = Scaled, group = Protein, color = factor(clust))) +
  geom_line(alpha = 0.4) +
  geom_line(data = pre.line.kmedcluster %>% filter(Protein == "sp|P06748|NPM_HUMAN"), color = "black") +
  geom_label(data = pre.line.kmedcluster %>% filter(Protein == "sp|P06748|NPM_HUMAN",Condition == "mean.7SK_1m"),
             aes(label = Protein), nudge_y = 2, color = "black") +
  geom_line(data = pre.line.kmedcluster %>% filter(Protein == "sp|O94992|HEXI1_HUMAN"), color = "black") +
  geom_label(data = pre.line.kmedcluster %>% filter(Protein == "sp|O94992|HEXI1_HUMAN",Condition == "mean.7SK_10m"),
             aes(label = Protein), nudge_y = 2, color = "black") +
  geom_line(data = pre.line.kmedcluster %>% filter(Protein == "sp|Q9H0S4|DDX47_HUMAN"), color = "red") +
  geom_label(data = pre.line.kmedcluster %>% filter(Protein == "sp|Q9H0S4|DDX47_HUMAN",Condition == "mean.ITS_1m"),
             aes(label = Protein), nudge_y = 2, color = "black") +
  geom_line(data = pre.line.kmedcluster %>% filter(Protein == "sp|Q7L2J0|MEPCE_HUMAN"), color = "red") +
  geom_label(data = pre.line.kmedcluster %>% filter(Protein == "sp|Q7L2J0|MEPCE_HUMAN",Condition == "mean.ITS_1m"),
             aes(label = Protein), nudge_y = 2, color = "black") +
  facet_wrap(vars(clust)) +
  theme_bw()




