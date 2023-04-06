library(tidyverse)
library(glmnet)
library(ungeviz) # For plots only. From C. Wilke https://github.com/wilkelab/ungeviz
######
## load data
######
# read all phenotypic data
datafile <- "RFG raw data 12_10_2018.xlsx"

rawdata <-  readxl::read_xlsx(datafile)
sigvalue <- 0.01 # significance threshold for ctrl vs drug t-tests

## Ampicillin/Sulbactam 64 and Ampicillin 2048 are duplicated, will have to figure that out

rawdata <- rawdata %>% filter(DrugType!="Ampicillin/Sulbactam 64")
procdata <- rawdata %>%
  separate(col = 'DrugType', into = c('drug', 'dose'), sep = ' ', remove = F)%>%
  mutate(dose = as.numeric(dose))%>%
  group_by(drug)%>%
  mutate(dose = dose/min(dose))%>%
  ungroup()

drug_data <- procdata %>%
  gather(key = 'haplotype', value='pheno', -c('drug', 'dose', 'DrugType')) %>%
  separate(col = 'haplotype', into=c('treat', 's1', 's2', 's3', 's4'), sep=c(1,2,3,4))%>%
  mutate(haplotype = str_c(s1,s2,s3,s4),
         batch = DrugType,
         dStrength = dose,
         dose = if_else(treat=='c', 0, dose), 
         treat = if_else(treat =='c', 'control', 'drug'),
         pheno = if_else(pheno <0, 0, pheno),
         DrugType = NULL)%>%
  group_by(drug)%>%
  mutate(dStrength= dense_rank(dStrength))%>%
  ungroup()%>%
  mutate(hapfact = factor(haplotype, 
                          levels = c("0000", "0001", "0010", "0100", "1000", 
                                                "0011", "0101", "1001", "0110", "1010", "1100", 
                                                "0111", "1011", "1101", "1110", "1111"), 
                          labels=c("MEGN","MEGD","MESN","MKGN","LEGN","MESD","MKGD","LEGD",
                                   "MKSN","LESN","LKGN","MKSD","LKSN","LKGD","LESD","LKSD"),
                          ordered =T),
         alleleCount =as.numeric(s1) +as.numeric(s2) +as.numeric(s3) + as.numeric(s4))

man_labels <- tibble(drug =c("Cefprozil", "Ceftazidime", "Piperacillin/Tazobactam", "Amoxicillin/ClavulanicAcid",
                             "Ampicillin/Sulbactam", "Ampicillin", "Amoxicillin", "Cefepime", "Cefotaxime", "Cefotetan"),
                     drug_form =c("Cefprozil", "Ceftazidime", "Piperacillin/\nTazobactam", "Amoxicillin/\nClavulanic Acid",
                                  "Ampicillin/\n Sulbactam", "Ampicillin", "Amoxicillin", "Cefepime", "Cefotaxime", "Cefotetan"),
                     lab =c("czl", "czd", "p+t", "a+c", "a+s","amp", "amx", "cfp", "cfx", "cft"))

treat_summary <- drug_data %>% select(batch, drug, dStrength) %>% 
  unique() %>%
  mutate(druglab = man_labels$lab[match(drug, man_labels$drug)],
         batchlab = paste(druglab, dStrength, sep="_"))

include_maxdose <- drug_data %>% 
  group_by(drug) %>%
  arrange(desc(dStrength)) %>% 
  filter(row_number()==1)%>%
  pull(batch)

##### Raw data explore ----

# Let's describe the effect on WT
wt_test <- procdata %>% group_by(DrugType) %>% 
  nest()%>%
  mutate(ttest = map_dbl(data, function(d) t.test(d$h0000, d$c0000, alternative = "less")$p.value))%>%
  mutate(sig = ttest < sigvalue)

wt_test|> filter(sig)
include_filter1 <- wt_test$DrugType[wt_test$sig]
##### Control subset (batch effect) ----

zeroset <- filter(drug_data, treat=='control')

normfactors <-filter(zeroset, haplotype=='0000', treat == 'control')%>%
  group_by(batch)%>%
  summarize(avg = mean(pheno), std = sd(pheno))%>%
  ungroup()


##### Drugabbility calcs ----

# max dose set, removing drugs that didn't have effect on WT ("include_filter1")
# use relY, growth relative to WT-ctrl

all_maxdose <- filter(drug_data, batch %in% include_maxdose, batch %in% include_filter1)%>%
  mutate(mu = normfactors$avg[match(batch, normfactors$batch)],
         std = normfactors$std[match(batch, normfactors$batch)],
         relY = pheno/mu,
         relY2 = (pheno-mu)/std) |>
  select(-c(mu, std)) %>%
  mutate(drug = man_labels$drug_form[match(drug, man_labels$drug)])

calc_resistance <- function(dp, hp, tp, df){
  x <- df %>% filter(drug == dp, haplotype==hp, treat==tp)
  y <- df %>% filter(drug == dp, haplotype=='0000', treat=='control')
  t.test(x$relY, y$relY, alternative = "less")$p.value
}

all_maxdose_drug <- all_maxdose %>% filter(treat =='drug')
all_maxdose_summary <- all_maxdose_drug%>%
  select(c(drug, haplotype, hapfact, treat, relY))%>%
  group_by(drug, haplotype, hapfact, treat)%>%
  summarize(relY=mean(relY))%>%
  ungroup()%>%
  unique()%>%
  mutate(resistance_p = pmap_dbl(list(dp=drug, hp=haplotype, tp=treat), calc_resistance, df = all_maxdose))%>%
  mutate(sig = resistance_p < sigvalue)

#### Drugabbility panel plots ----

drugg_df <-  all_maxdose_summary |>
  filter(treat=='drug')|>
  mutate(relY =  if_else(sig, relY, 1), # relative growth is 1 if there's no significant diff in growth
         susceptibility = 1 - relY)

vdrug <- drugg_df |>
  group_by(haplotype, hapfact) |>
  summarise(mu = mean(susceptibility), s = sd(susceptibility)/sqrt(n())) |>
  ungroup()

ddrug <- drugg_df |>
  group_by(drug) |>
  summarise(mu = mean(susceptibility), s = sd(susceptibility)/sqrt(n())) |>
  ungroup()

describe_drugability_plot <- ggplot(all_maxdose_drug, aes(y=relY, x=hapfact))+
  ungeviz::geom_hpline(data=all_maxdose_summary, aes(y = relY, color = sig & treat =='drug'), width=0.5, alpha=0.8)+
  geom_jitter(width=0.1, alpha=0.4, size=0.9) +
  geom_hline(yintercept=1, lty=2)+
  facet_grid(drug~.)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, vjust = 0.5), 
        legend.position = 'none', 
        strip.text.y = element_text(angle=0, size=8),
        strip.background = element_blank(),
        panel.border = element_rect(color='black', fill=NA),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.ticks.y = element_line())+
  scale_y_continuous(breaks =c(0, 1.0))+
  labs(x="", y="Relative Growth")+
  scale_color_manual(values=c('grey70', 'darkorange'))

vd_plot <- ggplot(vdrug)+
  geom_errorbar(aes(ymin=mu-0.01, ymax=mu+s, x=hapfact), width=0.15, color='grey40')+
  geom_col(aes(x=hapfact, y = mu), fill='grey70', width=0.5)+
  scale_y_continuous(limits = c(0,1), expand = c(0,0), breaks =c(0,0.5, 1), position = "right")+
  coord_cartesian(ylim = c(0,1))+
  theme_bw()+
  theme(strip.text.x = element_blank(),
        strip.placement = 'none',
        axis.line.y = element_line(),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.ticks.x = element_blank())+
  labs(x="", y="Variant vulnerability")


dd_plot <- ggplot(ddrug)+
  geom_errorbar(aes(ymin=mu-0.01, ymax=mu+s, x=drug), width=0.15, color='grey40')+
  geom_col(aes(x=drug, y = mu), fill='grey70', width=0.5)+
  facet_wrap(drug~., ncol=1, scales='free_y')+
  scale_y_continuous(limits = c(0,1), expand = c(0,0), breaks =c(0,0.5, 1),)+
  coord_flip(expand=0.5)+
  theme_bw()+
  theme(strip.text.x = element_blank(),
        strip.placement = 'none',
        axis.line.y = element_blank(),
        axis.line.x = element_line(),
        axis.text.y = element_blank(),
        axis.ticks.y.left = element_blank(),
        strip.text.y = element_text(angle=0, size=8),
        strip.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.ticks.y = element_line())+
  labs(x="", y="Drug applicability")

bltmp <- ggplot()+theme_void()

drug_panelplot <- cowplot::plot_grid(describe_drugability_plot, dd_plot, vd_plot, bltmp,
                                     ncol = 2, nrow = 2,  
                                     align = "hv", 
                                     axis= 'lrb', 
                                     rel_widths = c(2, 1), rel_heights = c(3, 1))


ggsave("fig_2_paneled_2drugabilities.png", drug_panelplot, width =8, height=8)
ggsave("fig_2_paneled_2drugabilities.pdf", drug_panelplot, width =8, height=8, useDingbats=FALSE)

#### V-druggability vs. One-step V-druggability. ----
#### One-step V-druggability is the average of all one-mutant neighbors
bin_to_label <- function(x){
  out <-  switch(x, 
                 "0000" = "MEGN",
                 "0001" = "MEGD",
                 "0010" = "MESN",
                 "0100" = "MKGN",
                 "1000" = "LEGN",
                 "0011" = "MESD",
                 "0101" = "MKGD",
                 "1001" = "LEGD",
                 "0110" = "MKSN",
                 "1010" = "LESN",
                 "1100" = "LKGN",
                 "0111" = "MKSD",
                 "1110" = "LKSN",
                 "1101" = "LKGD",
                 "1011" = "LESD",
                 "1111" = "LKSD",
                 character(0)  )
  out
}

hammdist <- function(x,y){
  xvec <- unlist(strsplit(x, split=""))
  yvec <- unlist(strsplit(y, split=""))
  sum(xvec != yvec) 
}

onesteplist <- tibble(x =vdrug$haplotype, y = vdrug$haplotype) |> 
  complete(x,y)|>
  mutate(ham = map2_dbl(x,y, hammdist)) |>
  filter(ham==1)

calc_onestep_drugg <- function(hap, druggabilities, dists){
  onedists <- dists |> filter(x == hap)
  onedruggs <- druggabilities |> filter(haplotype %in% onedists$y)
  mean(onedruggs$mu)
}

onestep_drugg <- vdrug |>
  mutate(mu_1step = map_dbl(haplotype, calc_onestep_drugg, druggabilities = vdrug, dists = onesteplist))|>
  mutate(deg = map_dbl(haplotype, function(x) sum(as.numeric(unlist(strsplit(x, split="")) ))))|>
  mutate(letts = map_chr(haplotype, bin_to_label))

onetep_lm <- summary(glm(mu_1step ~ mu, data =onestep_drugg, family = quasibinomial('logit')))
# Correlation value:
with(onestep_drugg, cor(mu_1step, mu))

onestep_plot <- ggplot(onestep_drugg,aes(x=mu, y=mu_1step))+ 
  geom_text(aes(label=letts), size=2.5)+
  #stat_smooth(method = "glm")+
  theme_bw()+
  labs(x="Variant vulnerability", y="One-step variant vulnerability")
ggsave("fig_3a_onestep_v_druggability.png", onestep_plot, width =5, height=3)
ggsave("fig_3a_onestep_v_druggability.pdf", onestep_plot, width =5, height=3, useDingbats=FALSE)


#### mutational path plots ----
library(igraph, pos='package:base', quietly = T)
library(ggraph, pos='package:base', quietly = T)
library(ggimage)

paths <- function(x){
  if(is.character(x)){
    out <-  switch(x, 
                   "0000" = c("0001", "0010", "0100", "1000"),
                   "0001" = c("0011", "0101", "1001"),
                   "0010" = c("0011", "0110", "1010"),
                   "0100" = c("0101", "0110", "1100"),
                   "1000" = c("1001", "1010", "1100"),
                   "0011" = c("0111", "1011"),
                   "0101" = c("0111", "1101"),
                   "1001" = c("1011", "1101"),
                   "0110" = c("0111", "1110"),
                   "1010" = c("1011", "1110"),
                   "1100" = c("1101", "1110"),
                   "0111" = c("1111"),
                   "1110" = c("1111"),
                   "1101" = c("1111"),
                   "1011" = c("1111"),
                   character(0)  )
  }else{
    out<- switch(as.character(x), 
                 "0" = c(1, 10, 100, 1000),
                 "1" = c(11, 101, 1001),
                 "10" = c(11, 110, 1010),
                 "100" = c(101, 110, 1100),
                 "1000" = c(1001, 1010, 1100),
                 "11" = c(111, 1011),
                 "101" = c(111, 1101),
                 "1001" = c(1011, 1101),
                 "110" = c(111, 1110),
                 "1010" = c(1011, 1110),
                 "1100" = c(1101, 1110),
                 "111" = c(1111),
                 "1110" = c(1111),
                 "1101" = c(1111),
                 "1011" = c(1111),
                 numeric(0)  )
  }
  return(out)
}
bin_allele_count <- function(bin){
  nchar(bin) - nchar(gsub('1', '', bin))
}
bin_int_to_char <- function(x){
  out<- switch(as.character(x), 
               "0" = "0000",
               "1" = "0001",
               "10" = "0010",
               "100" = "0100",
               "11" =  "0011",
               "101" = "0101",
               "110" = "0110",
               "111" = "0111",
               as.character(x) )
  return(out)
}

mutant_plot <- function(df,xvar='mutcount', yvar='relY'){
  
  nodelist <- tibble(binary = sort(map_chr(c(0, 1, 10, 100, 1000, 11, 101,110, 1001, 1010,1100, 111, 1011, 1101,1110, 1111), bin_int_to_char))) %>%
    left_join(., df, by='binary') |>
    mutate(letts = map_chr(binary, bin_to_label))
    #mutate(img = paste0("img/", binary, ".png"))
  
  xpos <- with(nodelist, get(xvar))
  ypos <- with(nodelist, get(yvar))

  templateNet <- tibble(from = nodelist$binary)%>%
    mutate(to = map(from, paths) )%>%unnest(cols=c(to)) 
  
  graf <- igraph::graph_from_data_frame(templateNet, vertices = nodelist$binary)

  out <- ggraph::ggraph(graf, layout = data.frame(x = xpos, y = ypos)) +
    ggraph::geom_edge_link0(edge_width=0.2, color='darkgreen', alpha=0.9)+ 
    geom_node_label(data=nodelist, aes(x=mutcount, y=relY, label = letts), label.padding = unit(0.09, "lines"), size=1.9)+
    scale_y_continuous(breaks = c(0,0.5, 1), expand = c(0, 0))+
    coord_cartesian(clip='off')+
    #theme(legend.position = 'none',
    #      plot.background = element_rect(fill='grey95'),
    #      plot.margin=unit(c(0,0,0,0), "null"),
    #      axis.ticks.x = element_blank(), 
    #      panel.grid.minor = element_blank(), 
    #      axis.text.x = element_blank(), 
    #      panel.grid.major = element_blank())+
    theme_minimal()+
    theme(axis.ticks.x = element_line(color='grey90'), panel.grid.minor = element_blank(), panel.grid.major.x  = element_blank() )+
    expand_limits(x=c(-0.1, 4.1), y = c(0,1))+
    labs(x="", y="")
  out
}

origami <- all_maxdose_summary |> 
  filter(treat=='drug')|>
  rename(binary = "haplotype")|>
  mutate(mutcount = map_int(binary,  bin_allele_count)) |>
  group_by(drug) |> nest() |>
  mutate(mplot = map(data, mutant_plot))|>
  mutate(mplot = map2(mplot,drug, function(p, d) p+labs(subtitle=d)))

origami_panel <- cowplot::plot_grid(
          origami$mplot[[1]]+labs(y="Relative growth")+theme(axis.text.x = element_blank()) , 
          ggplot()+theme_nothing(),
          origami$mplot[[2]]+labs(y="Relative growth")+theme(axis.text.x = element_blank()), 
          origami$mplot[[3]]+theme(axis.text.x = element_blank()), 
          origami$mplot[[4]]+labs(y="Relative growth")+theme(axis.text.x = element_blank()), 
          origami$mplot[[5]]+theme(axis.text.x = element_blank()), 
          origami$mplot[[6]]+labs(y="Relative growth", x="Number of mutations"), 
          origami$mplot[[7]]+labs(x="Number of mutations"), ncol=2)

ggsave(filename = 'fig_S1_byMutNumber.pdf', origami_panel, width = 7, height = 9, useDingbats=FALSE)
ggsave(filename = 'fig_S1_byMutNumber.png', origami_panel, width = 7, height = 9)

vdrug_origami <- mutant_plot(vdrug|>
                               rename(relY="mu", binary="haplotype")|>
                               mutate(mutcount = map_int(binary,  bin_allele_count)) )+
  labs(x="Number of mutations", y="Variant vulnerability")

ggsave(filename = 'fig_3b_vdrug_byMutNumber.pdf', vdrug_origami, width = 5, height = 3, useDingbats=FALSE)
ggsave(filename = 'fig_3b_vdrug_byMutNumber.png', vdrug_origami, width = 5, height = 3)

  

#### LASSO per drug ----
one_regreg <- function(df, my_alpha = 1){
  my_alpha <- 1 #1 is LASSO and 0 is ridge
  
  mat <- model.matrix(~ s1*s2*s3*s4, data = df)
  mat <- mat[,-which(colnames(mat) == "(Intercept)")]
  
  reg <- cv.glmnet(x = mat, 
                   y = df$relY2, 
                   nfolds = 100,
                   family='gaussian', 
                   alpha = my_alpha, 
                   standardize = FALSE,
                   keep = TRUE, 
                   grouped = FALSE,
                   intercept = TRUE)
  
  coefs <- coef(reg, s = "lambda.1se")
  
  tibble(term = rownames(coefs), value = unname(coefs[,1]) )
}

lasso_by_drug <- all_maxdose %>% filter(treat == 'drug')%>%
  group_by(drug)%>%
  nest()%>%
  mutate(terms = map(data, one_regreg) )

nice_terms <- rbind(
  c("s41", "A", "1"),            
  c("s31", "B", "1"),          
  c("s21", "C", "1"),           
  c("s11", "D", "1"),           
  c("s31:s41", "AxB", "2"),        
  c("s21:s41", "AxC", "2"),       
  c("s21:s31", "BxC", "2"),        
  c("s11:s41", "AxD", "2"),        
  c("s11:s31", "BxD", "2"),        
  c("s11:s21", "CxD", "2"),        
  c("s21:s31:s41","AxBxC", "3"),   
  c("s11:s31:s41","AxBxD", "3"),   
  c("s11:s21:s41","AxCxD", "3"),    
  c("s11:s21:s31","BxCxD", "3"),   
  c("s11:s21:s31:s41","AxBxCxD", "4"))


nice_terms <- rbind(
  c("s41", "(M69L) A", "1"),            
  c("s31", "(E104K) B", "1"),          
  c("s21", "(G238S) C", "1"),           
  c("s11", "(N276D) D", "1"),           
  c("s31:s41", "AxB", "2"),        
  c("s21:s41", "AxC", "2"),       
  c("s21:s31", "BxC", "2"),        
  c("s11:s41", "AxD", "2"),        
  c("s11:s31", "BxD", "2"),        
  c("s11:s21", "CxD", "2"),        
  c("s21:s31:s41","AxBxC", "3"),   
  c("s11:s31:s41","AxBxD", "3"),   
  c("s11:s21:s41","AxCxD", "3"),    
  c("s11:s21:s31","BxCxD", "3"),   
  c("s11:s21:s31:s41","AxBxCxD", "4"))

 

lasso_drug_nice <- lasso_by_drug %>% 
  unnest(terms) %>%
  filter(term !="(Intercept)")%>%
  mutate(term = factor(term, levels = nice_terms[,1], labels = nice_terms[,2]),
         term_cat = nice_terms[match(term, nice_terms[,2]),3])

paleta <- c("grey80", wesanderson::wes_palette("Zissou1", n = 5)[c(1,3,4,5)])

lasso_plot <- ggplot(lasso_drug_nice)+
  geom_hline(yintercept=0, linewidth=0.2)+
  geom_col(aes(x = term, y = value, fill = term_cat)) +
  facet_grid(drug~., scales = 'free_y')+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, vjust = 0.5), 
        legend.position = 'none', 
        strip.text.y = element_text(angle=0, size=8),
        strip.background = element_blank(),
        panel.border = element_rect(color='black', fill=NA),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.ticks.y = element_line())+
  labs(x="", y="Effect on growth (SD)")+
  scale_y_continuous(breaks = scales::trans_breaks(identity, identity, n = 3))+
  scale_fill_manual(values = paleta)

ggsave("fig_4_relative_lasso_by_drug.pdf", lasso_plot, width=5, height=6, useDingbats=FALSE)
ggsave("fig_4_relative_lasso_by_drug.png", lasso_plot, width=5, height=6)

lasso_plotb <- ggplot(lasso_drug_nice)+
  geom_vline(xintercept=0, linewidth=0.2)+
  geom_col(aes(y = term, x = value, fill = term_cat)) +
  facet_grid(.~drug, scales = 'free_x')+
  theme_bw()+
  theme(axis.text.x = element_text(angle=0, size=8), 
        legend.position = 'none', 
        strip.text.y = element_text(),
        strip.background = element_blank(),
        panel.border = element_rect(color='black', fill=NA),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.ticks.x = element_line())+
  labs(y="", x="Effect on growth (SD)")+
  scale_x_continuous(breaks = scales::trans_breaks(identity, identity, n = 3))+
  scale_fill_manual(values = paleta)+
  scale_y_discrete(limits = rev(levels(lasso_drug_nice$term)))
ggsave("fig_4alt_relative_lasso_by_drug.png", lasso_plotb, width=8, height=3)

# Effect sizes:
lasso_drug_nice |> arrange(desc(value))

