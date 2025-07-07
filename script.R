#Installation von PopdictR
#devtools::install_github("jogrue/regexhelpeR")
#devtools::install_github("jogrue/multidictR")
#devtools::install_github("jogrue/popdictR")

#downgrade quanteda
#remotes::install_version("quanteda", version = "2.1.2", repos = "http://cran.us.r-project.org")
#pak::pkg_install("quanteda@2.1.2") #Alternative to the other command

#Library

library(readxl)
library(dplyr)
library(quanteda)
library(regexhelpeR)
library(multidictR)
library(popdictR)
library(tidyverse)
library(writexl)
library(stringr)
library(MLmetrics)
library(binom)
library(DescTools)
library(lme4)
library(influence.ME)
library(broom)
library(sandwich)
library(lmtest)

######################################################################################################################
########################################   Data loading and prep    ##################################################
######################################################################################################################

#set wd
setwd("")

#loading of the transcripts
blattnamen <- excel_sheets("raw_data_zentrum.xlsx")
path <- "raw_data_zentrum.xlsx"
data_list_zentrum <- lapply(blattnamen, function(blatt) {
  read_excel(path, sheet = blatt)
})
names(data_list_zentrum) <- blattnamen
data_zentrum <- bind_rows(data_list_zentrum, .id = "Blattname")


blattnamen <- excel_sheets("raw_data_pro.xlsx")
path <- "raw_data_pro.xlsx"
data_list_pro <- lapply(blattnamen, function(blatt) {
  read_excel(path, sheet = blatt)
})
names(data_list_pro) <- blattnamen
data_pro <- bind_rows(data_list_pro, .id = "Blattname")

data_zentrum <- data_zentrum %>%
  rename(show = Blattname)
data_pro <- data_pro %>%
  rename(show = Blattname)

data_pro <- data_pro %>%
  mutate(
    start_time = format(start_time, format = "%H:%M:%S"),
    end_time   = format(end_time, format = "%H:%M:%S")
  )

data_zentrum <- data_zentrum %>%
  mutate(
    start_time = format(start_time, format = "%H:%M:%S"),
    end_time   = format(end_time, format = "%H:%M:%S")
  )

#data prep for dictionary analysis -> remove symbols, punctions and numbers
data_pro <- data_pro %>%
  mutate(text = text %>%
           str_replace_all("[[:punct:]]", " ") %>%
           str_replace_all("[[:digit:]]", " ") %>%           
           str_replace_all("[^A-Za-zÄÖÜäöüß//s]", " ") %>%    
           str_squish())                                    

data_zentrum <- data_zentrum %>%
  mutate(text = text %>%
           str_replace_all("[[:punct:]]", " ") %>%           
           str_replace_all("[[:digit:]]", " ") %>%           
           str_replace_all("[^A-Za-zÄÖÜäöüß//s]", " ") %>%     
           str_squish())

#combinding the dataframes
data_pro <- data_pro %>%
  mutate(sender = "Puls4")

data_zentrum <- data_zentrum %>%
  mutate(sender = "ORF2")

data_all <- bind_rows(data_pro, data_zentrum)

#creation of the corpus
corpus_all <-corpus(data_all, text_field = "text")



######################################################################################################################
########################################   Dictionary Analysis    ####################################################
######################################################################################################################

#run the popdictR dictionary on the corpus
results_all <- run_popdict(corpus_all)

#convert the results to data frames
results_all <- convert(results_all, to = "data.frame")




######################################################################################################################
########################################   Dictionary Validiation    #################################################
######################################################################################################################

#creating turn_id
results_all <- results_all %>%
  group_by(show) %>%
  mutate(
    speaker_change = (speaker != lag(speaker, default = first(speaker))),
    turn_id = cumsum(replace_na(speaker_change, FALSE))
  ) %>%
  ungroup()

 turns <- results_all %>%
   group_by(show, turn_id, speaker) %>%
   summarize(
     full_text = paste(text, collapse = ". "),
     start_doc_id = min(doc_id),
     end_doc_id = max(doc_id),
     n_sentences = sum(n_sentences, na.rm = TRUE),
     dict_gruendl_2020 = sum(dict_gruendl_2020, na.rm = TRUE),
     .groups = "drop"
   )


#random sample hits
set.seed(321)
treffer_gruendl_sample <- turns %>%
  filter(dict_gruendl_2020 > 0) %>%
  sample_n(25)

#write_xlsx(treffer_gruendl_sample, path = "stichprobe_hits.xlsx")

#random sample non-hits
set.seed(321)
sample <- turns %>%
  filter(dict_gruendl_2020 == 0) %>%
  dplyr::sample_n(100)
#write_xlsx(sample, path = "stichprobe_nonhits.xlsx")

#loading of the manual coded sample
sample_coded <- read_excel("sample_coded_combined.xlsx")

TP <- sum(sample_coded$dict_gruendl_2020 == 1 & sample_coded$manual_coded ==1)
FP <- sum(sample_coded$dict_gruendl_2020 == 1 & sample_coded$manual_coded == 0)
FN <- sum(sample_coded$dict_gruendl_2020 == 0 & sample_coded$manual_coded == 1)

precision <- TP / (TP + FP)
print(precision) 
recall <- TP / (TP + FN)
print(recall) 
f1 <- 2 * (precision * recall) / (precision + recall)
print(f1) 

#amount of false negatives in the sample
print(FN)

#percentage of all speechparts that were coded populist
mean(turns$dict_gruendl_2020 == 1, na.rm = TRUE) * 100

#calculte pearson corr
sample_coded_neg <-  sample_coded %>%
  filter(dict_gruendl_2020 == 0)
sample_coded_pos <- sample_coded %>%
  filter(dict_gruendl_2020 == 1)

set.seed(321)
sample_coded_neg <- sample_coded_neg %>% 
  slice_sample(n = 25)

sample_coded_small <- bind_rows(sample_coded_neg, sample_coded_pos)
cor(sample_coded_small$manual_coded, sample_coded_small$dict_gruendl_2020, method = "pearson")


######################################################################################################################
#############################################   Data Prep    #########################################################
######################################################################################################################

#### creation of dummy variables for hypotheses testing ####

#Dummy variables for H1.1

results_all$role_politician <- dplyr::case_when(
  is.na(results_all$role) ~ NA_character_,
  results_all$role == "Politiker/in" ~ "politician",
  TRUE ~ "non-politician"
)

#Dummy variable for H1.2
results_all$party_populist <- dplyr::case_when(
  is.na(results_all$party) ~ NA_character_,
  results_all$party == "FPÖ" ~ "Populist",
  results_all$party == "NA" ~ "NA",
  TRUE ~ "not Populist"
)
results_all$party_populist <- na_if(results_all$party_populist, "NA")

#Dummy variable for H1.3
results_all$party_pole <- dplyr::case_when(
  is.na(results_all$party) ~ NA_character_,
  results_all$party == "FPÖ" ~ "pole",
  results_all$party == "GRÜNE" ~ "pole",
  results_all$party == "NA" ~ "NA",
  TRUE ~ "not pole"
)
results_all$party_pole <- na_if(results_all$party_pole, "NA")


#### Data prep for regression analysis ####

results_all$show <- as.factor(results_all$show)

results_all$sender <- as.factor(results_all$sender)
results_all$sender <- relevel(results_all$sender, ref = "ORF2")

results_all$role_politician <- as.factor(results_all$role_politician)
results_all$role_politician <- relevel(results_all$role_politician, ref = "non-politician")

results_all$party_populist <- as.factor(results_all$party_populist)
results_all$party_populist <- relevel(results_all$party_populist, ref = "not Populist")


results_all$party_pole <- as.factor(results_all$party_pole)
results_all$party_pole <- relevel(results_all$party_pole, ref = "not pole")



######################################################################################################################
#############################################   Exporting Dataset    #################################################
######################################################################################################################


#saving the dataset
saveRDS(results_all, file = "dataset.rds")
#loading the datafile
#results_all <- readRDS("dataset.rds")

#creation of a new subset with only politicians
results_politicians <- subset(results_all, role_politician == "politician")

gruendl_dict <- popdictR::gruendl_terms
print(gruendl_dict)

######################################################################################################################
########################################    simple analysis    #######################################################
######################################################################################################################

########################################   general descriptive analysis  ##############################################

#Number of Different Speakers
count(results_all)
results_all %>%
  filter(sender == "ORF2") %>%
  count()
results_all %>%
  filter(sender == "Puls4") %>%
  count()

#sentences by role
table_1<- results_all %>%
  group_by(role) %>%
  summarise(
    unique_speakers = n_distinct(speaker),  
    total_sentences = n()                   
  ) %>%
  arrange(desc(total_sentences))
print(table_1)
#write_xlsx(table_1, path = "Datenanalyse/table_1.xlsx")

#sentences by party
table_2<- results_all %>%
  group_by(party) %>%
  summarise(
    unique_speakers = n_distinct(speaker),  
    total_sentences = n()                   
  ) %>%
  arrange(desc(total_sentences))
print(table_2)
#write_xlsx(table_2, path = "Datenanalyse/table_2.xlsx")

#table to show how many populist sentences were coded
table(results_all$dict_gruendl_2020, useNA = "always")
table(results_politicians$dict_gruendl_2020, useNA = "always")

#percentage of all sentences that were coded populist
mean(results_all$dict_gruendl_2020 == 1, na.rm = TRUE) * 100


########################################   descriptive and simple inferential analysis  ##############################

#populist ratio for every role
table_3 <- results_all %>%
  group_by(role) %>%
  summarize(
    total_sentences = sum(n_sentences),
    populist_hits = sum(dict_gruendl_2020),
    populist_ratio = populist_hits / total_sentences * 100
  ) %>%
  filter(total_sentences >= 10) %>% 
  arrange(desc(populist_ratio))
print(table_3)
#write_xlsx(table_3, path = "Datenanalyse/table_3.xlsx")

#populist ratio for every party
table_4 <- results_all %>%
  group_by(party) %>%
  summarize(
    total_sentences = sum(n_sentences),
    populist_hits = sum(dict_gruendl_2020),
    populist_ratio = populist_hits / total_sentences * 100
  ) %>%
  filter(total_sentences >= 10) %>% 
  arrange(desc(populist_ratio))
print(table_4)
#(table_4, path = "Datenanalyse/table_4.xlsx")

#######################################################
################## H1.1-H1.3 ##########################
#######################################################


################################################## H1.1 Populism - Pole Parties ###################################

#Ratio
party_results_pole <- results_all %>%
  filter(party_pole != "NA") %>%
  group_by(party_pole) %>%
  summarize(
    total_sentences = sum(n_sentences, na.rm = TRUE),
    populist_hits = sum(dict_gruendl_2020, na.rm = TRUE),
    populist_ratio = populist_hits / total_sentences * 100
  )
print(party_results_pole)

#simple graphical plot
ggplot(party_results_pole, aes(x = reorder(party_pole, -populist_ratio), y = populist_ratio, fill = party_pole)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(
    "pole" = "black",
    "not pole" = "darkgrey"
  )) +
  labs(
    title = "",
    x = "Parteityp",
    y = "Sätze mit populistischen Ausdrücken"
  ) +
  scale_y_continuous(
    labels = scales::percent_format(scale = 1),
    limits = c(0, 1)
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(),
    legend.position = "none"
  )

#significance Tests
party_results_pole <- party_results_pole %>%
  mutate(non_populist_hits = total_sentences - populist_hits)

chi_matrix <- party_results_pole %>%
  select(party_pole, populist_hits, non_populist_hits) %>%
  column_to_rownames("party_pole") %>%
  as.matrix()
print(chi_matrix)

chisq.test(chi_matrix)
chisq.test(chi_matrix, correct = FALSE)
fisher.test(chi_matrix)
CramerV(chi_matrix)

#################################################### H1.2 Populism - Populist Parties ########################################

#ratio
party_results_pop <- results_all %>%
  group_by(party_populist) %>%
  summarize(
    total_sentences = sum(n_sentences, na.rm = TRUE),
    populist_hits = sum(dict_gruendl_2020, na.rm = TRUE),
    populist_ratio = populist_hits / total_sentences * 100
  )
print(party_results_pop)

#drop all non-politicians
party_results_pop <- party_results_pop %>% 
  filter(party_populist != "NA")
print(party_results_pop)

#simple graphical plot
ggplot(party_results_pop, aes(x = reorder(party_populist, -populist_ratio), y = populist_ratio, fill = party_populist)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(
    "Populist" = "black",
    "not Populist" = "darkgrey"
  )) +
  labs(
    title = "",
    x = "Parteityp",
    y = "Sätze mit populistischen Ausdrücken"
  ) +
  scale_y_continuous(
    labels = scales::percent_format(scale = 1),
    limits = c(0, 1)
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(),
    legend.position = "none"
  )

#significance Tests
party_results_pop <- party_results_pop %>%
  mutate(non_populist_hits = total_sentences - populist_hits)

chi_matrix <- party_results_pop %>%
  select(party_populist, populist_hits, non_populist_hits) %>%
  column_to_rownames("party_populist") %>%
  as.matrix()
print(chi_matrix)

chisq.test(chi_matrix)
CramerV(chi_matrix)
fisher.test(chi_matrix)


################################### H1.3 Politicians / non-politicians ##########################################

#ratio
role_politician <- results_all %>%
  group_by(role_politician) %>%
  summarize(
    total_sentences = sum(n_sentences, na.rm = TRUE),
    populist_hits = sum(dict_gruendl_2020, na.rm = TRUE),
    populist_ratio = populist_hits / total_sentences * 100
  ) %>%
  arrange(desc(populist_ratio))
print(role_politician)

#simple graphical plot
ggplot(role_politician, aes(x = reorder(role_politician, -populist_ratio), 
                            y = populist_ratio, 
                            fill = role_politician)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(
    "politician" = "black",
    "non-politician" = "darkgrey"
  )) +
  labs(
    title = "",
    x = "Rolle",
    y = "Sätze mit populistischen Ausdrücken"
  ) +
  scale_y_continuous(
    labels = scales::percent_format(scale = 1),
    limits = c(0, 1)
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text()
  )


#significance tests
role_politician <- role_politician %>%
  mutate(non_populist_hits = total_sentences - populist_hits)

chi_matrix <- role_politician %>%
  select(role_politician, populist_hits, non_populist_hits) %>%
  column_to_rownames("role_politician") %>%
  as.matrix()
print(chi_matrix)

chisq.test(chi_matrix)
CramerV(chi_matrix)
fisher.test(chi_matrix)


#######################################################
################## media ##############################
#######################################################

########################################### H2.1 public/private broadcasters  ##########################################

#ratio
channel_results <- results_all %>%
  group_by(sender) %>%
  summarize(
    total_sentences = sum(n_sentences),
    populist_hits = sum(dict_gruendl_2020),
    populist_ratio = populist_hits / total_sentences * 100
  )
print(channel_results)

#simple plot
ggplot(channel_results, aes(x = reorder(sender, -populist_ratio), y = populist_ratio, fill = sender)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(
    "ORF2" = "darkgrey",
    "Puls4" = "black"
  )) +
  labs(
    title = "",
    x = "Sender",
    y = "Sätze mit populistischen Ausdrücken"
  ) +
  scale_y_continuous(labels = scales::percent_format(scale = 1), 
                     limits = c(0, 1)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(),
    legend.position = "none"
  )

#significance Tests
channel_results <- channel_results %>%
  mutate(non_populist_hits = total_sentences - populist_hits)

chi_matrix <- channel_results %>%
  select(sender, populist_hits, non_populist_hits) %>%
  column_to_rownames("sender") %>%
  as.matrix()
print(chi_matrix)

chisq.test(chi_matrix)
CramerV(chi_matrix)
fisher.test(chi_matrix)

###################################### H2.2 / public/private + actors ####################################


########## H1.1 ###########

#ratio
channel_results_pole <- results_all %>%
  filter(!is.na(party_pole)) %>%
  mutate(role_sender = paste(sender, party_pole, sep = " – ")) %>%
  group_by(role_sender) %>%
  summarize(
    populist_hits = sum(dict_gruendl_2020),
    total = n(),
    non_populist_hits = total - populist_hits,
    populist_ratio = populist_hits / total * 100
  )
print(channel_results_pole)

#simple plot
channel_results_pole %>%
  ggplot(aes(
    x = reorder(role_sender, -populist_ratio),
    y = populist_ratio,
    fill = substr(role_sender, 1, 4)  
  )) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(
    "ORF2" = "grey",
    "Puls" = "black"
  )) +
  labs(
    title = "",
    x = "Sender – Parteipol",
    y = "Sätze mit populistischen Ausdrücken"
  ) +
  scale_y_continuous(
    limits = c(0, 1.5),
    labels = scales::percent_format(scale = 1)
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

#significance tests
chi_matrix <- channel_results_pole %>%
  select(role_sender, populist_hits, non_populist_hits) %>%
  column_to_rownames("role_sender") %>%
  as.matrix()
print(chi_matrix)

chisq.test(chi_matrix)
fisher.test((chi_matrix))
CramerV(chi_matrix)




########## H1.2 ###########

#ratio
channel_results_pop <- results_all %>%
  filter(!is.na(party_populist)) %>%
  mutate(role_sender = paste(sender, party_populist, sep = " – ")) %>%
  group_by(role_sender) %>%
  summarize(
    total_sentences = n(),  
    populist_hits = sum(dict_gruendl_2020),
    populist_ratio = populist_hits / total_sentences * 100,
    .groups = "drop"
  )
print(channel_results_pop)

#simple plot
channel_results_pop %>%
  ggplot(aes(
    x = reorder(role_sender, -populist_ratio),
    y = populist_ratio,
    fill = substr(role_sender, 1, 4)  # sender extrahieren für Farbe
  )) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(
    "ORF2" = "grey",
    "Puls" = "black"
  )) +
  labs(
    title = "",
    x = "Sender – Parteityp",
    y = "Sätze mit populistischen Ausdrücken"
  ) +
  scale_y_continuous(
    limits = c(0, 1.5),
    labels = scales::percent_format(scale = 1)
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

#sigificance tests
channel_results_pop <- channel_results_pop %>%
  mutate(non_populist_hits = total_sentences - populist_hits)
chi_matrix <- channel_results_pop %>%
  select(role_sender, populist_hits, non_populist_hits) %>%
  column_to_rownames("role_sender") %>%
  as.matrix()
print(chi_matrix)

chisq.test(chi_matrix)
fisher.test(chi_matrix)
CramerV(chi_matrix)


########## H1.3 ###########

#ratio
channel_results_role <- results_all %>%
  group_by(sender, role_politician) %>%
  summarize(
    total_sentences = n(),  
    populist_hits = sum(dict_gruendl_2020),
    populist_ratio = populist_hits / total_sentences * 100,
    .groups = "drop"
  )
print(channel_results_role)

#simple plot
ggplot(channel_results_role, aes(
  x = reorder(paste(sender, role_politician, sep = " – "), -populist_ratio),
  y = populist_ratio,
  fill = sender
)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(
    "ORF2" = "grey",
    "Puls4" = "black"
  )) +
  labs(
    title = "",
    x = "Sender – Rolle",
    y = "Sätze mit populistischen Ausdrücken"
  ) +
  scale_y_continuous(
    labels = scales::percent_format(scale = 1),
    limits = c(0, 1.5)
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

#sigificance tests
channel_results_role_2 <- results_all %>%
  mutate(role_sender = paste(sender, role_politician, sep = " – ")) %>%
  group_by(role_sender) %>%
  summarize(
    populist_hits = sum(dict_gruendl_2020),
    total = n(),
    non_populist_hits = total - populist_hits
  )

chi_matrix <- channel_results_role_2 %>%
  select(role_sender, populist_hits, non_populist_hits) %>%
  column_to_rownames("role_sender") %>%
  as.matrix()
print(chi_matrix)
chisq.test(chi_matrix)
CramerV(chi_matrix)
fisher.test(chi_matrix)



######################################################################################################################
########################################    Regression models    #####################################################
######################################################################################################################


#####################   Logit_politicians - H1.3   ##########################

model_0.4 <- glm(dict_gruendl_2020 ~ role_politician + sender,
                         data = results_all,
                         family = binomial)
summary(model_0.4)

exp(model_0.4$coef[2])
exp(model_0.4$coef[3])

#### Testing Model Fit ###

model_null <- glm(dict_gruendl_2020 ~ 1, data = results_all, family = binomial)

ll_full <- as.numeric(logLik(model_0.4))
ll_null <- as.numeric(logLik(model_null))

G <- 2 * (ll_full - ll_null)
pchisq(G, df = 2, lower.tail = FALSE)

pseudo_r2 <- 1 - (model_0.4$deviance/model_0.4$null.deviance)
pseudo_r2

### testing model assumptions ###

#testing outliers
plot(model_0.4, which = 4, id.n = 3) 
abline(h = 4 / (nobs(model_0.4) - length(coef(model_0.4))), col = "red", lty = 2)

influence.measures(model_0.4)
cd <- cooks.distance(model_0.4)
plot(cd) 


cutoff <- 4 / (nobs(model_0.4) - length(coef(model_0.4)))
influential_points <- which(cd > cutoff)
influential_df_1 <- results_all[cd > cutoff, ] 

model.data <- augment(model_0.4) %>% 
  mutate(index = 1:n()) 
inf_out_1 <- model.data %>% 
  filter(abs(.std.resid) > 3)

#testing sample size
table(results_all$role_politician, results_all$dict_gruendl_2020)
table(results_all$sender, results_all$dict_gruendl_2020)

#testing Multicolinearity
car::vif(model_0.4)



#####################   Logit_pole_parties H1.1   ##########################

results_politicians <- results_politicians %>%
  filter(!is.na(party_pole), !is.na(party_populist))


model_0.2 <- glm(dict_gruendl_2020 ~ party_pole + sender,
                         data = results_politicians,
                         family = binomial)
summary(model_0.2)
exp(model_0.2$coef[2])
exp(model_0.2$coef[3])

### Testing Model Fit ###
model_null_pol <- glm(dict_gruendl_2020 ~ 1, data = results_politicians, family = binomial)

ll_full <- as.numeric(logLik(model_0.2))
ll_null <- as.numeric(logLik(model_null_pol))

G <- 2 * (ll_full - ll_null)
pchisq(G, df = 2, lower.tail = FALSE)

pseudo_r2 <- 1 - (model_0.2$deviance/model_0.2$null.deviance)
pseudo_r2

### testing model assumptions ###

#testing outliers
plot(model_0.2, which = 4, id.n = 3) 
abline(h = 4 / (nobs(model_0.2) - length(coef(model_0.2))), col = "red", lty = 2)

influence.measures(model_0.2)
cd <- cooks.distance(model_0.2)
plot(cd) 

cutoff <- 4 / (nobs(model_0.2) - length(coef(model_0.2)))
influential_points <- which(cd > cutoff)
influential_df_2 <- results_politicians[cd > cutoff, ] 

model.data <- augment(model_0.2) %>% 
  mutate(index = 1:n()) 
inf_out_2 <- model.data %>% 
  filter(abs(.std.resid) > 3)

#testing sample size
table(results_politicians$party_pole, results_politicians$dict_gruendl_2020)
table(results_politicians$sender, results_politicians$dict_gruendl_2020)

#testing Multicolinearity
car::vif(model_0.2)



#####################   Logit_populist_parties H1.2   ##########################

model_0.3 <- glm(dict_gruendl_2020 ~ party_populist + sender,
                         data = results_politicians,
                         family = binomial)
summary(model_0.3)

exp(model_0.3$coef[2])
exp(model_0.3$coef[3])


#Testing Model Fit

ll_full <- as.numeric(logLik(model_0.3))
ll_null <- as.numeric(logLik(model_null_pol))

G <- 2 * (ll_full - ll_null)
pchisq(G, df = 2, lower.tail = FALSE)

pseudo_r2 <- 1 - (model_0.3$deviance/model_0.3$null.deviance)
pseudo_r2


### testing model assumptions ###

#testing outliers
plot(model_0.3, which = 4, id.n = 3) 
abline(h = 4 / (nobs(model_0.3) - length(coef(model_0.3))), col = "red", lty = 2)

influence.measures(model_0.3)
cd <- cooks.distance(model_0.3)
plot(cd)


cutoff <- 4 / (nobs(model_0.3) - length(coef(model_0.3)))

influential_points <- which(cd > cutoff)
influential_df_3 <- results_politicians[cd > cutoff, ] 

model.data <- augment(model_0.2) %>% 
  mutate(index = 1:n()) 
inf_out_3 <- model.data %>% 
  filter(abs(.std.resid) > 3)

#testing sample size
table(results_politicians$party_populist, results_politicians$dict_gruendl_2020)
table(results_politicians$sender, results_politicians$dict_gruendl_2020)

#testing Multicolinearity
car::vif(model_0.3)


#Export all Influential Outliers for the regression models
saveRDS(influential_df_1, file = "inf_out_polit.rds")
saveRDS(influential_df_2, file = "inf_out_pol.rds")
saveRDS(influential_df_3, file = "inf_out_pop.rds")



