# process pq 46
# based on https://uc-r.github.io/kmeans_clustering
library(tidyverse)
library(pheatmap)
#library(ComplexHeatmap)
library(cluster)
library(factoextra)

#set your working directory. previously from Schweppe lab
#setwd("C:\\Users\\Mintaka\\Dropbox (Personal)\\UW\\Collaborations\\Triforce\\Shechner_ODB\\rna_odb_paper")
#setwd("C:\\Users\\Mammmals\\Dropbox\\UW\\Collaborations\\Triforce\\Shechner_ODB\\rna_odb_paper")
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/UW/Shechner Lab/Biotinylated Oligo exp/Mass Spec/MS analysis/ODB paper 2021 analysis")

#odb = read_tsv("protein_quant_46.tsv")
#ITS labeling titration 
odb = read_tsv("protein_quant_48.tsv")

scaled = odb %>%
  filter(!grepl("##",`Protein Id`) & !grepl("contaminant",`Protein Id`)) %>%
  select(contains("Protein Id") | contains("scaled"))

scaled.piv = scaled %>%
  pivot_longer(!"Protein Id", names_to = "Condition", values_to = "Scaled")

colnames(scaled.piv) = c("Protein","Condition","Scaled")

pheatmap(scaled %>% select(contains("scaled")), 
         name = "ODB", #title of legend
         column_title = "Scaled", row_title = "Condition",
         cluster_cols = FALSE, labels_row = scaled$`Protein Id`,
)

#nucleolar
ggplot(scaled.piv %>% filter(str_detect(Protein,"NPM")), aes(x = Condition, y = Scaled)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw()
#excluded
ggplot(scaled.piv %>% filter(str_detect(Protein,"TRA2B")), aes(x = Condition, y = Scaled)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw()


scaled.mean = scaled %>%
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

scaled.mean = scaled.mean %>%
  add_column(
    q.1s = p.adjust(scaled.mean$t.its1_1s,method = "BH"),
    q.1m = p.adjust(scaled.mean$t.its1_1m,method = "BH"),
    q.10m = p.adjust(scaled.mean$t.its1_10m,method = "BH"),
    q.100m = p.adjust(scaled.mean$t.its1_100m,method = "BH")
  )

scaled.mean.p = scaled.mean %>%
  select(contains("t.") | contains("Protein")) %>%
  pivot_longer(-contains("Protein"))

#ggplot(scaled.mean.q,aes(x= value)) +
 # geom_histogram()

scaled.mean.q = scaled.mean %>%
  select(contains("q.") | contains("Protein")) %>%
  pivot_longer(-contains("Protein"))

ggplot(scaled.mean.q,aes(x= value, fill = name)) +
  geom_histogram() +
  facet_wrap(~name)

pheatmap(scaled.mean %>%
           select(mean.its1_1m,mean.7sk_1m), 
         name = "ODB", #title of legend
         column_title = "Scaled", row_title = "Condition",
         cluster_cols = FALSE, labels_row = scaled.mean$Protein,
)

pheatmap(scaled.mean %>%
           select(contains("l2r.")) %>%
           filter_if(~is.numeric(.), all_vars(!is.infinite(.))), 
         name = "ODB", #title of legend
         column_title = "Log Ratio", row_title = "Condition",
         cluster_cols = FALSE, labels_row = scaled.mean$Protein,
)

pre.line = scaled.mean %>%
  select(contains("mean.") | contains("Protein Id")) %>%
  pivot_longer(!"Protein Id", names_to = "Condition", values_to = "Scaled")

ggplot(pre.line , aes(x = Condition, y = Scaled, group = `Protein Id`)) +
  geom_line(alpha = 0.4) +
  theme_bw()

scaled.mean.z = scaled.mean %>%
  select(contains("mean")) %>% scale()

### clustering
k1 = kmeans(scaled.mean %>%
              select(contains("mean")),centers=1,nstart=25)
k2 = kmeans(scaled.mean %>%
              select(contains("mean")),centers=2,nstart=25)
k3 = kmeans(scaled.mean %>%
              select(contains("mean")),centers=3,nstart=25)
k4 = kmeans(scaled.mean %>%
              select(contains("mean")),centers=4,nstart=25)
k5 = kmeans(scaled.mean %>%
              select(contains("mean")),centers=5,nstart=25)
k12 = kmeans(scaled.mean %>%
               select(contains("mean")),centers=12,nstart=25)
k12.lr = kmeans(scaled.mean %>%
               select(contains("l2r."))%>%
                 filter_if(~is.numeric(.), all_vars(!is.infinite(.)))
               ,centers=12,nstart=25)

#pam of scaled.mean is kmedoid analysis
pam12 = pam(scaled.mean %>%
               select(contains("mean")),12,metric = "euclidean", stand = FALSE)
pam12.lr = pam(scaled.mean %>%
              select(contains("l2r."))%>%
              filter_if(~is.numeric(.), all_vars(!is.infinite(.))),
            12,metric = "euclidean", stand = FALSE)

pam6 = pam(scaled.mean %>%
              select(contains("mean")),6,metric = "euclidean", stand = FALSE)
#kmedoid plot of 12 clusetr
fviz_cluster(pam12, data = scaled.mean %>%
               select(contains("mean")),
             geom = "point")

fviz_cluster(k1,data = scaled.mean %>%
               select(contains("mean")),
             geom = "point")
fviz_cluster(k2,data = scaled.mean %>%
               select(contains("mean")),
             geom = "point")
fviz_cluster(k3,data = scaled.mean %>%
               select(contains("mean")),
             geom = "point")
fviz_cluster(k4,data = scaled.mean %>%
               select(contains("mean")),
             geom = "point")
fviz_cluster(k5,data = scaled.mean %>%
               select(contains("mean")),
             geom = "point")
fviz_cluster(k12,data = scaled.mean %>%
               select(contains("mean")),
             geom = "point")

fviz_nbclust(scaled.mean %>%
               select(contains("mean")),kmeans,method="wss")

fviz_nbclust(scaled.mean %>%
               select(contains("mean")),kmeans,method="silhouette")


gap_stat <- clusGap(scaled.mean %>%
                      select(contains("mean")), FUN = kmeans, nstart = 25,
                    K.max = 20, B = 50)
fviz_gap_stat(gap_stat)
#kmedoid version of gapstat_put in 7/30/24 not sure where is in code
gap_stat_kmed <- clusGap(scaled.mean %>% select(contains("mean")),
                    FUN = pam,
                    nstart = 25,
                    K.max = 20, # max clusters to consider
                    B = 50)     # total bootstrapped iterations
fviz_gap_stat(gap_stat_kmed)
fviz_gap_stat(gap_stat_kmed, linecolor = "steelblue", maxSE = list(method =
                                                                "globalSEmax", SE.factor = 1))
gap_stat_kmed_1 <- clusGap(scaled.mean %>% select(contains("mean")),
                         FUN = pam,
                         nstart = 25,
                         K.max = 15, # max clusters to consider
                         B = 50)     # total bootstrapped iterations
fviz_gap_stat(gap_stat_kmed_1, linecolor = "steelblue", maxSE = list(method =
                                                                     "globalSEmax", SE.factor = 2))

#make new for 12 cluster with l2r and cluster assignment, also includes k means and kmedoid analysis
pre.line.cluster = scaled.mean %>%
  select(contains("mean.") | contains("Protein")) %>%
  add_column(clust = k12$cluster) %>%
  add_column(clust.pam = pam12$cluster) %>%
  pivot_longer(contains("mean."), names_to = "Condition", values_to = "Scaled")

pre.line.cluster.lr = scaled.mean %>%
  select(contains("l2r.") | contains("Protein")) %>%
  filter_if(~is.numeric(.), all_vars(!is.infinite(.))) %>%
  add_column(clust = k12.lr$cluster) %>%
  add_column(clust.pam = pam12.lr$cluster) %>%
  pivot_longer(contains("l2r."), names_to = "Condition", values_to = "Log Ratio")

# ggplot(pre.line.cluster , aes(x = Condition, y = Scaled, group = Protein, color = factor(clust))) +
#   geom_line(alpha = 0.4) +
#   theme_bw()

ggplot(pre.line.cluster , aes(x = Condition, y = Scaled, group =  `Protein Id`, color = factor(clust))) +
  geom_line(alpha = 0.4) +
  geom_line(data = pre.line.cluster %>% filter( `Protein Id` == "sp|P06748|NPM_HUMAN"), color = "black") +
  geom_label(data = pre.line.cluster %>% filter( `Protein Id` == "sp|P06748|NPM_HUMAN",Condition == "mean.its1_100m"),
             aes(label = `Protein Id`), nudge_y = 2, color = "black") +
  facet_wrap(vars(clust)) +
  theme_bw()

ggplot(pre.line.cluster.lr , aes(x = Condition, y = `Log Ratio`, group =  `Protein Id`, color = factor(clust))) +
  geom_line(alpha = 0.4) +
  geom_line(data = pre.line.cluster.lr %>% filter( `Protein Id` == "sp|P06748|NPM_HUMAN"), color = "black") +
  geom_label(data = pre.line.cluster.lr %>% filter( `Protein Id` == "sp|P06748|NPM_HUMAN",Condition == "mean.its1_100m"),
             aes(label = `Protein Id`), nudge_y = 2, color = "black") +
  facet_wrap(vars(clust),scales = "free") +
  theme_bw()

plc.k1 = pre.line.cluster %>%
  filter(clust == 1)

plc.k5 = pre.line.cluster %>%
  filter(clust == 5)

#kmediods
ggplot(pre.line.cluster , aes(x = Condition, y = Scaled, group =  `Protein Id`, color = factor(clust.pam))) +
  geom_line(alpha = 0.4) +
  geom_line(data = pre.line.cluster %>% filter( `Protein Id` == "sp|P06748|NPM_HUMAN"), color = "black") +
  geom_label(data = pre.line.cluster %>% filter( `Protein Id` == "sp|P06748|NPM_HUMAN",Condition == "mean.its1_100m"),
             aes(label = `Protein Id`), nudge_y = 2, color = "black") +
  facet_wrap(vars(clust.pam)) +
  theme_bw()

ggplot(pre.line.cluster.lr , aes(x = Condition, y = `Log Ratio`, group =  `Protein Id`, color = factor(clust.pam))) +
  geom_line(alpha = 0.4) +
  geom_line(data = pre.line.cluster.lr %>% filter( `Protein Id` == "sp|P06748|NPM_HUMAN"), color = "black") +
  geom_label(data = pre.line.cluster.lr %>% filter( `Protein Id` == "sp|P06748|NPM_HUMAN",Condition == "mean.its1_100m"),
             aes(label = `Protein Id`), nudge_y = 2, color = "black") +
  facet_wrap(vars(clust.pam),scales = "free") +
  theme_bw()

#pull out the individual cluster of proteins
plc.pam4 = pre.line.cluster %>%
  filter(clust.pam == 4)
plc.pam8 = pre.line.cluster %>%
  filter(clust.pam == 8)
plc.pam6 = pre.line.cluster %>%
  filter(clust.pam == 6)
plc.pam11 = pre.line.cluster %>%
  filter(clust.pam == 11)

plc.pam7.lr = pre.line.cluster.lr %>% filter(clust.pam == 7)

fviz_nbclust(scaled.mean %>%
               select(contains("mean")),pam,method="wss")

#plot individual cluster protein profiles
ggplot(pre.line.cluster %>% filter(clust.pam == 4), aes(x = Condition, y = Scaled, group =  `Protein Id`, color = factor(clust.pam))) +
  geom_line(alpha = 0.4) +
  theme_bw()

ggplot(plc.pam7.lr, aes(x = Condition, y = `Log Ratio`, group =  `Protein Id`, color = factor(clust.pam))) +
  geom_line(alpha = 0.4) +
  theme_bw()

## HPA

##nuc = read_tsv("C:\\Users\\Mammmals\\Dropbox\\UW\\Collaborations\\Triforce\\Shechner_ODB\\HPA\\subcell_location_Nucleoli.tsv")
#nuc = read_tsv("C:\\Users\\Mammmals\\Dropbox\\UW\\Collaborations\\Triforce\\Shechner_ODB\\HPA\\subcell_location_Nucleoli_All.tsv")
#nuc = read_tsv("C:\\Users\\Mintaka\\Dropbox (Personal)\\UW\\Collaborations\\Triforce\\Shechner_ODB\\HPA\\subcell_location_Nucleoli_All.tsv")
nuc = read_tsv("subcell_location_Nucleoli_All.tsv")

#https://www.proteinatlas.org/humanproteome/subcellular/nucleoli
sum.nuc = tibble(
  Nucleoli = 1410,
  Other.Organelles = 11631,
  Secreted = 774,
  Not.Mapped = 6275
)

#add column for uniprot identifier to scaled means of each protein for all conditions
sm.nuc = scaled.mean %>%
  add_column(Uniprot.wIso = unlist(
    lapply(scaled.mean$`Protein Id`, FUN = function(x){ return(unlist(str_split(x,'\\|'))[2]) }))
  ) %>%
  mutate(Uniprot = unlist(
    lapply(Uniprot.wIso, FUN = function(x){ return(unlist(str_split(x,'\\-'))[1]) }))
    )
#add Reliability classifier to our data sorted by uniprot identifier
sm.nuc.merge = sm.nuc %>%
  left_join(nuc %>% select(Uniprot,`Reliability (IF)`),by="Uniprot")

#add k12 medoid analysis and log ratios then make a long DF
plc.lr.nuc = sm.nuc.merge %>%
  select(contains("l2r.") | contains("Protein") | `Reliability (IF)`) %>%
  filter_if(~is.numeric(.), all_vars(!is.infinite(.))) %>%
  add_column(clust = k12.lr$cluster) %>%
  add_column(clust.pam = pam12.lr$cluster) %>%
  pivot_longer(contains("l2r."), names_to = "Condition", values_to = "Log Ratio")

#plot all supported and enhanced proteins present in our exp
ggplot(plc.lr.nuc %>% filter(`Reliability (IF)` == "Supported" | `Reliability (IF)` == "Enhanced" ),
       aes(x = Condition, y = `Log Ratio`, group =  `Protein Id`, color = factor(clust.pam))) +
  geom_line(alpha = 0.4) + labs(x = "ITS Biotinylation time titration", y="Log 2 Ratio of ITS time/7SK 1 min")
  theme_bw()

#plot for each protein in our dataset that is labeled by HPA as approved, enhanced etc classifier
ggplot(plc.lr.nuc , aes(x = Condition, y = `Log Ratio`, group =  `Protein Id`, color = factor(`Reliability (IF)`))) +
  geom_line(alpha = 0.4) +
  geom_line(data = plc.lr.nuc %>% filter( `Protein Id` == "sp|P06748|NPM_HUMAN"), color = "black") +
  geom_label(data = plc.lr.nuc %>% filter( `Protein Id` == "sp|P06748|NPM_HUMAN",Condition == "l2r.to7.100m"),
             aes(label = `Protein Id`), color = "black") +
  facet_wrap(vars(`Reliability (IF)`),scales = "free") +
  theme_classic() +
  ylim(-10,30)


#for each kmedoid cluster, highlight the HPA classifier
ggplot(plc.lr.nuc , aes(x = Condition, y = `Log Ratio`, group =  `Protein Id`, color = factor(`Reliability (IF)`))) +
  geom_line(alpha = 0.4, aes(color = factor(`Reliability (IF)`))) +
  geom_line(data = plc.lr.nuc %>% filter( `Protein Id` == "sp|P06748|NPM_HUMAN"), color = "black") +
  geom_label(data = plc.lr.nuc %>% filter( `Protein Id` == "sp|P06748|NPM_HUMAN",Condition == "l2r.to7.100m"),
             aes(label = `Protein Id`), color = "black") +
  facet_wrap(vars(clust.pam),scales = "free") +
  theme_classic() +
  ylim(-10,30)

## establish clusters for fishers exact testing
#matrix of count from each classified protein in our dataset
es.plc.lr.nuc = plc.lr.nuc %>%
  distinct(`Protein Id`,clust.pam,`Reliability (IF)`) %>%
  group_by(clust.pam) %>% 
  summarise(
    e.only = sum(grepl("Enhanced",`Reliability (IF)`)),
    a.only = sum(grepl("Approved",`Reliability (IF)`)),
    s.only = sum(grepl("Supported",`Reliability (IF)`)),
    u.only = sum(grepl("Uncertain",`Reliability (IF)`)),
    na.only = sum(is.na(`Reliability (IF)`)),
    n = n(),
    es = sum(grepl("Enhanced",`Reliability (IF)`) | grepl("Supported",`Reliability (IF)`)),
            n.es = n() - es) %>%
  add_column(nuc = 1410, all = sum(sum.nuc) - nuc)

#number of each classified protein in our dataset
intersList = es.plc.lr.nuc %>%
  summarise(
    e.only.sum = sum(e.only),
    a.only.sum = sum(a.only),
    s.only.sum = sum(s.only),
    u.only.sum = sum(u.only),
    na.only.sum = sum(na.only),
    es.sum = sum(es)
  )
  

plc.lr.nuc %>%
  distinct(`Protein Id`,clust.pam,`Reliability (IF)`) %>%
  group_by(clust.pam) %>% 
  summarise(es = sum(grepl("Enhanced",`Reliability (IF)`) |
                       grepl("Supported",`Reliability (IF)`) |
                       grepl("Approved",`Reliability (IF)`)
                       ),
            n.es = n() - es) %>%
  add_column(nuc = 1410, all = sum(sum.nuc) - nuc)


#fisher exact test and BH adjusted p val called fish.q
es.plc.lr.nuc = es.plc.lr.nuc %>%
  rowwise() %>%
  mutate(fish = fisher.test(matrix(unlist(c_across("es":"all")),nrow = 2))$p.value)
  

es.plc.lr.nuc = es.plc.lr.nuc %>%
  add_column(fish.q = p.adjust(es.plc.lr.nuc$fish, method = "BH"))

mosaicplot(matrix(unlist(es.plc.lr.nuc[1,-1]),nrow = 2))

ft.t = fisher.test(matrix(unlist(es.plc.lr.nuc[1,-1]),nrow = 2))

#output clusters taht are <0.05 adjust p vale from fisher exact test
es.plc.lr.nuc %>% filter(fish.q < 0.05) %>% select(clust.pam)

# determining previous unannotated (in HPA) nucleolar proteins, first pull out cluster 1,4,7,12
plc.lr.nuc.specific = plc.lr.nuc %>%
  filter(clust.pam %in% unlist(es.plc.lr.nuc %>% filter(fish.q < 0.05) %>% select(clust.pam)))

#generate list of proteins in these clusters 1,4,7,12 
plc.lr.nuc.specific.UNIprot = plc.lr.nuc.specific %>%
  add_column(Uniprot.wIso = unlist(
    lapply(plc.lr.nuc.specific$`Protein Id`, FUN = function(x){ return(unlist(str_split(x,'\\|'))[2]) }))
  ) %>%
  mutate(Uniprot = unlist(
    lapply(Uniprot.wIso, FUN = function(x){ return(unlist(str_split(x,'\\-'))[1]) }))
  )
NucleoProteinDistinct = plc.lr.nuc.specific.UNIprot %>% select(`Protein Id`, clust.pam, `Reliability (IF)`, Uniprot) %>% unique()


write.csv(NucleoProteinDistinct, file = "AllNucleolarProteins_in_cluster1,4,7,12.csv")

#plot set of 1,4,7,12 clusters
ggplot(plc.lr.nuc.specific , aes(x = Condition, y = `Log Ratio`, group =  `Protein Id`, color = clust.pam)) +
  geom_line(alpha = 0.4) +
  geom_hline(yintercept = 0) +
  theme_classic()

#plot the proteins that are classified as NA in the nucleolar protein clusters 1,4,7,12
ggplot(plc.lr.nuc.specific %>%
         filter(is.na(`Reliability (IF)`)), aes(x = Condition, y = `Log Ratio`, group =  `Protein Id`, color = clust.pam)) +
  geom_line(alpha = 0.4) +
  geom_hline(yintercept = 0) +
  theme_classic()

#list of nucleolar proteins that are HPA non assigned
newNucProt.sum = plc.lr.nuc.specific %>%
  summarise(
    e.only.sum = sum(e.only),
    a.only.sum = sum(a.only),
    s.only.sum = sum(s.only),
    u.only.sum = sum(u.only),
    na.only.sum = sum(na.only),
    es.sum = sum(es)
  )

#generate list of 118 proteins that are HPA NA
newNucProt = plc.lr.nuc.specific %>% filter(is.na(`Reliability (IF)`)) %>% select(`Protein Id`, clust.pam, Condition, `Log Ratio`)
newNucProtDistinct = plc.lr.nuc.specific %>% filter(is.na(`Reliability (IF)`)) %>% select(`Protein Id`, clust.pam) %>% unique()
newNucProtDistinct.UNIprot = newNucProtDistinct %>%
  add_column(Uniprot.wIso = unlist(
    lapply(newNucProtDistinct$`Protein Id`, FUN = function(x){ return(unlist(str_split(x,'\\|'))[2]) }))
  ) %>%
  mutate(Uniprot = unlist(
    lapply(Uniprot.wIso, FUN = function(x){ return(unlist(str_split(x,'\\-'))[1]) }))
  )
write.csv(newNucProtDistinct.UNIprot, file = "HPA NA nucleolar proteins.csv")

#Figure ODB panel B
odb.panel.b = read_tsv("protein_quant_35.tsv")

pb.scaled = odb.panel.b %>%
  filter(!grepl("##",`Protein Id`) & !grepl("contaminant",`Protein Id`)) %>%
  select(contains("Protein Id") | contains("scaled"))

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
    t.7sk_10m = t.test(c_across("7sk-10m_1 scaled":"7sk-10m_3 scaled"),
                       c_across("scrmbl-10m_1 scaled":"scrmbl-10m_3 scaled"))$p.value,
    t.its1_1m = t.test(c_across("its1_1m_1 scaled":"its1_1m_3 scaled"),
                        c_across("scrmbl-10m_1 scaled":"scrmbl-10m_3 scaled"))$p.value,
    l2r.toS.71m = log(mean.7sk_1m/mean.scr_10m,2),
    l2r.toS.710m = log(mean.7sk_10m/mean.scr_10m,2),
    l2r.toS.i1m = log(mean.its1_1m/mean.scr_10m,2)
  )

pb.sm = pb.sm %>%
  add_column(
    q.1s = p.adjust(pb.sm$t.7sk_1m,method = "BH"),
    q.1m = p.adjust(pb.sm$t.7sk_10m,method = "BH"),
    q.10m = p.adjust(pb.sm$t.its1_1m,method = "BH")
  )

pb.plc.scaled = pb.sm %>%
  select(contains("scaled") | contains("Protein")) %>%
  filter_if(~is.numeric(.), all_vars(!is.infinite(.))) %>%
  pivot_longer(contains("scaled"), names_to = "Condition", values_to = "Scaled Ratio") %>%
  mutate(simpleCond = unlist(
    lapply(Condition, FUN = function(x){ return(unlist(str_split(x,'\\_'))[1]) }))
  )

pb.plc.lr = pb.sm %>%
  select(contains("l2r") | contains("Protein")) %>%
  filter_if(~is.numeric(.), all_vars(!is.infinite(.))) %>%
  pivot_longer(contains("l2r"), names_to = "Condition", values_to = "Log Ratio")


pb.plc.lr %>% filter(grepl("DDX47",`Protein Id`) |
                   grepl("RPF1",`Protein Id`)|
                   grepl("UTP6",`Protein Id`)|
                   grepl("FTSJ3",`Protein Id`)|
                   grepl("UBTF",`Protein Id`)|
                   grepl("NOL10",`Protein Id`)
)

ggplot(pb.plc.scaled %>% filter(grepl("DDX47",`Protein Id`) |
                             grepl("NPM",`Protein Id`)|
                             grepl("UTP6",`Protein Id`)|
                             grepl("NUCL",`Protein Id`)|
                               grepl("TLN1",`Protein Id`)|
                               grepl("CCAR2",`Protein Id`)|
                             grepl("NOL10",`Protein Id`)),
       aes(factor(simpleCond,levels = c("blank scaled","scrmbl-1m","scrmbl-10m","7sk-1m","7sk-10m","its1")),
           `Scaled Ratio`,fill = simpleCond)) +
  geom_boxplot() +
  facet_wrap(~`Protein Id`,nrow = 1,scales = "free") +
  theme_classic()

ggplot(pb.plc.scaled %>% filter(grepl("DDX47",`Protein Id`) |
                                  grepl("NPM",`Protein Id`)|
                                  grepl("NUCL",`Protein Id`)|
                                  grepl("NOL10",`Protein Id`)),
       aes(factor(simpleCond,levels = c("blank scaled","scrmbl-1m","scrmbl-10m","7sk-1m","7sk-10m","its1")),
           `Scaled Ratio`,fill = Condition)) +
  geom_bar(stat="identity",position="dodge") +
  facet_wrap(~`Protein Id`,nrow = 1,scales = "free") +
  theme_classic()

pheatmap(pb.scaled %>% select(contains("scaled")), 
         name = "ODB", #title of legend
         column_title = "Scaled", row_title = "Condition",
         cluster_cols = FALSE, labels_row = pb.scaled$`Protein Id`,
)

pheatmap(pb.scaled %>% select(contains("scaled")), 
         name = "ODB", #title of legend
         column_title = "Scaled", row_title = "Condition",
         cluster_cols = T, labels_row = pb.scaled$`Protein Id`,
)

pheatmap(pb.scaled %>% select(`its1_1m_1 scaled`,`its1_1m_2 scaled`,`its1_1m_3 scaled`,
                              `7sk-1m_1 scaled`,`7sk-1m_2 scaled`,`7sk-1m_3 scaled`), 
         name = "ODB", #title of legend
         column_title = "Scaled", row_title = "Condition",
         cluster_cols = FALSE, labels_row = pb.scaled$`Protein Id`,
         clustering_method = "median"
)

#Figure ODB panel D
library(corrr)
odb.panel.d = read_tsv("protein_quant_93.tsv")

pd.scaled = odb.panel.d %>%
  filter(!grepl("##",`Protein Id`) & !grepl("contaminant",`Protein Id`)) %>%
  select(contains("Protein Id") | contains("scaled"))

pd.sm = pd.scaled %>%
  rownames_to_column(var = "Protein") %>%
  rowwise() %>%
  mutate(
    mean.7_panc = mean(c_across("7sk_panc_1 scaled":"7sk_panc_3 scaled")),
    mean.i.panc = mean(c_across("its1_panc_1 scaled":"its1_panc_3 scaled")),
    mean.i.8988t = mean(c_across("its1_8988t_1 scaled":"its1_8988t_3 scaled")),
    mean.i.aspc = mean(c_across("its1_aspc_1 scaled":"its1_aspc_3 scaled")),
    mean.i.suit3 = mean(c_across("its1_suit_1 scaled":"its1_suit_3 scaled")),
    sd.7_panc = sd(c_across("7sk_panc_1 scaled":"7sk_panc_3 scaled")),
    sd.i.panc = sd(c_across("its1_panc_1 scaled":"its1_panc_3 scaled")),
    sd.i.8988t = sd(c_across("its1_8988t_1 scaled":"its1_8988t_3 scaled")),
    sd.i.aspc = sd(c_across("its1_aspc_1 scaled":"its1_aspc_3 scaled")),
    sd.i.suit3 = sd(c_across("its1_suit_1 scaled":"its1_suit_3 scaled")),
    t.i.panc = t.test(c_across("its1_panc_1 scaled":"its1_panc_3 scaled"),
                      c_across("7sk_panc_1 scaled":"7sk_panc_3 scaled"))$p.value,
    t.i.8988 = t.test(c_across("its1_8988t_1 scaled":"its1_8988t_3 scaled"),
                       c_across("7sk_panc_1 scaled":"7sk_panc_3 scaled"))$p.value,
    t.i.aspc = t.test(c_across("its1_aspc_1 scaled":"its1_aspc_3 scaled"),
                       c_across("7sk_panc_1 scaled":"7sk_panc_3 scaled"))$p.value,
    t.i.suit3 = t.test(c_across("its1_suit_1 scaled":"its1_suit_3 scaled"),
                       c_across("7sk_panc_1 scaled":"7sk_panc_3 scaled"))$p.value,
    l2r.toIp.ip = log(mean.i.panc/mean.7_panc,2),
    l2r.toIp.i8 = log(mean.i.8988t/mean.7_panc,2),
    l2r.toIp.ia = log(mean.i.aspc/mean.7_panc,2),
    l2r.toIp.is = log(mean.i.suit3/mean.7_panc,2)
  )

pd.sm = pd.sm %>%
  add_column(
    q.toIp.ip = p.adjust(pd.sm$t.i.panc,method = "BH"),
    q.i.8988 = p.adjust(pd.sm$t.i.8988,method = "BH"),
    q.i.aspc = p.adjust(pd.sm$t.i.aspc,method = "BH"),
    q.i.suit3 = p.adjust(pd.sm$t.i.suit3,method = "BH")
  )

cor.pd = pd.sm %>% select(contains("l2r")) %>% correlate()

rplot(cor.pd)

cor.pd = pd.sm %>% select(contains("scaled")) %>% correlate()

rplot(cor.pd)

plc.lr.pd = pd.sm %>%
  select(contains("l2r.") | contains("Protein")) %>%
  filter_if(~is.numeric(.), all_vars(!is.infinite(.))) %>%
  pivot_longer(contains("l2r."), names_to = "Condition", values_to = "Log Ratio")

ggplot(plc.lr.pd,
       aes(x = Condition, y = `Log Ratio`, group =  `Protein Id`)) +
  geom_line(alpha = 0.4) +
  theme_bw()

pheatmap(pd.sm %>%
           select(contains("l2r.")) %>%
           filter_if(~is.numeric(.), all_vars(!is.infinite(.))), 
         name = "ODB Pancreatic", #title of legend
         column_title = "Log Ratio", row_title = "Condition",
         #cluster_cols = FALSE,
         labels_row = pd.sm$Protein,
)

pd.plc.scaled = pd.sm %>%
  select(contains("scaled") | contains("Protein")) %>%
  filter_if(~is.numeric(.), all_vars(!is.infinite(.))) %>%
  pivot_longer(contains("scaled"), names_to = "Condition", values_to = "Scaled Ratio") %>%
  mutate(simpleCond = unlist(
    lapply(Condition, FUN = function(x){ 
      return(paste0(unlist(str_split(x,'\\_'))[1],
                    "_",
                    unlist(str_split(x,'\\_'))[2]))
      }))
  )

pd.plc.lr = pd.sm %>%
  select(contains("l2r") | contains("Protein")) %>%
  filter_if(~is.numeric(.), all_vars(!is.infinite(.))) %>%
  pivot_longer(contains("l2r"), names_to = "Condition", values_to = "Log Ratio")

ggplot(pd.plc.scaled %>% filter(grepl("DDX47",`Protein Id`) |
                                  grepl("NPM",`Protein Id`)|
                                  grepl("UTP6",`Protein Id`)|
                                  grepl("NUCL",`Protein Id`)|
                                  grepl("TLN1",`Protein Id`)|
                                  grepl("CCAR2",`Protein Id`)|
                                  grepl("NOL10",`Protein Id`)),
       aes(factor(simpleCond,levels = c("blank scaled","scrmbl-1m","scrmbl-10m","7sk-1m","7sk-10m","its1")),
           `Scaled Ratio`,fill = simpleCond)) +
  geom_boxplot() +
  facet_wrap(~`Protein Id`,nrow = 1) +
  theme_classic()

pd.pam12.lr = pam(pd.sm %>%
                 select(contains("l2r."))%>%
                 filter_if(~is.numeric(.), all_vars(!is.infinite(.))),
               12,metric = "euclidean", stand = FALSE)

fviz_cluster(pd.pam12.lr,data = pd.sm %>%
               select(contains("l2r.")),
             geom = "point")

pd.pr = pd.sm %>%
  select(contains("scaled")) %>% prcomp()

pd.pr.out = as.data.frame(pd.pr[[2]]) %>%
  rownames_to_column() %>%
  mutate(simpleCond = unlist(
    lapply(rowname, FUN = function(x){ 
      return(paste0(unlist(str_split(x,'\\_'))[1],
                    "_",
                    unlist(str_split(x,'\\_'))[2]))
    }))
  )

ggplot(pd.pr.out, aes(PC1, PC2, color = simpleCond)) +
  geom_point() +
  theme_classic()

