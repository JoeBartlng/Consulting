###
### Key Driver Analysis
### Date: 26Mar2021
### Author: Joe Bartling
###

rm(list=ls())
dev.off(dev.list()["RStudioGD"])

library(tidyverse)
library(misty)
library(relaimpo)
library(ggpubr)
library(labelled)


setwd("/home/joe/Documents/Consulting/Aaron Zaber - RSM")

colorBlind1  <- c("#000000", "#E69F00")

colorBlind2  <- c("#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

colorBlind3  <- c("#F0E442", "#0072B2", "#D55E00", "#CC79A7")

colorBlind4  <- c("#D55E00", "#CC79A7")

colorBlind5 <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")


colorBlind6 <- c( "#DDCC77", "#117733", "#332288", "#AA4499", 
                 "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")


relweights <- function(fit,LOB,DV,pal,...){   
  # fit <- Audit_fit1
  # LOB <- 'Audit'
  # DV <- 'Overall Experience'
  # pal <- colorBlind1 
  R <- cor(fit$model)   
  nvar <- ncol(R)          
  rxx <- R[2:nvar, 2:nvar] 
  rxy <- R[2:nvar, 1]      
  svd <- eigen(rxx)        
  evec <- svd$vectors                           
  ev <- svd$values         
  delta <- diag(sqrt(ev))  
  lambda <- evec %*% delta %*% t(evec)        
  lambdasq <- lambda ^ 2   
  beta <- solve(lambda) %*% rxy           
  rsquare <- colSums(beta ^ 2)                   
  rawwgt <- lambdasq %*% beta ^ 2    
  import <- (rawwgt / rsquare) * 100 
  lbls <- names(fit$model[2:nvar])   
  rownames(import) <- lbls
  colnames(import) <- "Weights"
  n <- nobs(fit)
  
  import <- as.data.frame(import) %>% 
    rownames_to_column(var='IV') %>% 
    mutate(labls=paste0(round(Weights,1),"%"), method='Johnson') %>% 
    dplyr::select(IV, Weights,labls,method)
  
  Shapley <- calc.relimp(fit , method = c("lmg"),rela=TRUE)

  
  
  lmg <- as.data.frame(Shapley$lmg) %>% 
    rownames_to_column(var = "IV") %>% 
    mutate(Weights=`Shapley$lmg`*100, labls=paste0(round(`Shapley$lmg`*100,1),"%"),method='Shapley') %>% 
    dplyr::select(IV, Weights,labls,method)
  
  imptnce <- import %>% 
    bind_rows(lmg) %>% 
    mutate(method=as.factor(method),
           IV=substr(IV,1,75) ,
           IV=gsub('[.]',' ',IV))
  
  imptnce$method <- relevel(imptnce$method, "Shapley")
  
  top <- round(max(imptnce$Weights),0)+5
  
  plt_johnson <- ggplot(data=imptnce, aes(x=reorder(IV, Weights), y=Weights, fill=method)) +
    geom_bar(stat="identity",, position=position_dodge())+
    scale_fill_manual(values=pal) +
    geom_text(aes(label=labls), position=position_dodge(width=0.9), hjust=-.25) +
    ylab("% of R-Square") +
    scale_y_continuous(limits = c(0, top)) +
    xlab("Predictor Variables") +
    labs(title = paste0("Relative Importance of IVs"),
         subtitle = paste0("LOB: ",LOB, " -- DV: ",DV),
         caption = paste0("Model Rsq=",round(rsquare,2)," n=",n)) +
    theme(axis.text.y = element_text(angle = 0, hjust=1,size=10),legend.position = c(0.9, 0.1)) +
    guides(fill = guide_legend(reverse = TRUE)) +
    coord_flip()

  return(plt_johnson)
}

keydrv1 <- read.sav("Post-Engagement Responses_OPEN.sav") 

vec <- NULL

for (i in seq_along(names(keydrv1))) {
  vec <- c(vec,var_label(keydrv1)[[i]])
  output <- vec
}
output

keydrv2 <- as_tibble(keydrv1)

names(keydrv2) <- make.names(vec)

str(keydrv2)

keydrv2 <- keydrv2 %>% 
  mutate(LOB_P=if_else(`Line.of.Business` ==1,'Audit',
                       if_else(`Line.of.Business` ==2,'Tax','Consulting')))

table(keydrv2$Line.of.Business)
## Data Setup 

Audit <- keydrv2 %>% 
  filter(LOB_P=='Audit')

Audit2 <- Audit[,c(23,4:12,21:22)]

Audit2 %>%
  dplyr::select(-LOB_P) %>%  # replace to your needs
  summarise_all(funs(sum(is.na(.))))

sum(complete.cases(Audit2))

str(Audit2)

Audit_no_newteam <- Audit2 %>% 
  dplyr::select(-All.new.team.members.demonstrated.an.appropriate.understanding.of.my.organization.and.the.engagement.details)

Audit_no_newteam  %>%
  dplyr::select(-LOB_P) %>%  # replace to your needs
  summarise_all(funs(sum(is.na(.))))

sum(complete.cases(Audit_no_newteam ))


Tax <- keydrv2 %>% 
  filter(LOB_P=='Tax')

Tax2 <- Tax[,c(23,4:6,10,13,16:22)]

Tax2 %>%
  dplyr::select(-LOB_P) %>%  # replace to your needs
  summarise_all(funs(sum(is.na(.))))

sum(complete.cases(Tax2))

Consulting <- keydrv2 %>% 
  filter(LOB_P=='Consulting')

Consulting2 <- Consulting [,c(23,3:6,10:15,21:22)]

Consulting2 %>%
  dplyr::select(-LOB_P) %>%  # replace to your needs
  summarise_all(funs(sum(is.na(.))))

str(Consulting2)

table(Consulting2$Consulting.Service.Line)

sum(complete.cases(Consulting2))

Consulting_complete <- Consulting2[complete.cases(Consulting2),]

table(Consulting_complete$Consulting.Service.Line)

Consulting_Risk <- Consulting_complete %>% 
  filter(Consulting.Service.Line=="Risk Consulting")

Consulting_Technology <- Consulting_complete %>% 
  filter(Consulting.Service.Line=="Technology Consulting")

### Model Function

reg_it <- function(df,dv,ivs,ivf){   

  df_names <- names(df)
  
  df_f1 <- as.formula(paste0(df_names[dv],'~',paste0(df_names[ivs:ivf],collapse = '+')))
  
  df_fit1 <- lm(df_f1,df)
  
  summary(df_fit1 )
  return(df_fit1)
}

##Audit

Audit_fit1 <- reg_it(Audit2,2,5,12)
Audit_fit2 <- reg_it(Audit2,3,5,12)
Audit_fit3 <- reg_it(Audit2,4,5,12)

Audit_results1 <- relweights(Audit_fit1,'Audit','Overall Experience', colorBlind1 )
Audit_results1 

Audit_results2 <- relweights(Audit_fit2,'Audit', 'Final Deliverables', colorBlind1 )
Audit_results2 

Audit_results3 <- relweights(Audit_fit3,'Audit','LTR', colorBlind1)
Audit_results3 

Audit_no_newteam_fit1 <- reg_it(Audit_no_newteam ,2,5,11)
Audit_no_newteam_fit2 <- reg_it(Audit_no_newteam ,3,5,11)
Audit_no_newteam_fit3 <- reg_it(Audit_no_newteam ,4,5,11)

Audit_resultsNNT1 <- relweights(Audit_no_newteam_fit1,'Audit-w/o New Team Members','Overall Experience', colorBlind6)
Audit_resultsNNT1 

Audit_resultsNNT2 <- relweights(Audit_no_newteam_fit2,'Audit-w/o New Team Members','Final Deliverables', colorBlind6)
Audit_resultsNNT2 

Audit_resultsNNT3 <- relweights(Audit_no_newteam_fit3,'Audit-w/o New Team Members','LTR', colorBlind6)
Audit_resultsNNT3 

### Tax

Tax_fit1 <- reg_it(Tax2,2,5,13)
Tax_fit2 <- reg_it(Tax2,3,5,13)
Tax_fit3 <- reg_it(Tax2,4,5,13)

Tax_results1 <- relweights(Tax_fit1,'Tax','Overall Experience', colorBlind2 )
Tax_results1 

Tax_results2 <- relweights(Tax_fit2,'Tax', 'Final Deliverables', colorBlind2 )
Tax_results2 

Tax_results3 <- relweights(Tax_fit3,'Tax','LTR', colorBlind2)
Tax_results3 

### Consulting

Consulting_fit1 <- reg_it(Consulting2,3,6,13)
Consulting_fit2 <- reg_it(Consulting2,4,6,13)
Consulting_fit3 <- reg_it(Consulting2,5,6,13)

Consulting_results1 <- relweights(Consulting_fit1,'Consulting','Overall Experience', colorBlind3 )
Consulting_results1 

Consulting_results2 <- relweights(Consulting_fit2,'Consulting', 'Final Deliverables', colorBlind3 )
Consulting_results2 

Consulting_results3 <- relweights(Consulting_fit3,'Consulting','LTR', colorBlind3)
Consulting_results3 

### Consulting Risk

Consulting_Risk_fit1 <- reg_it(Consulting_Risk,3,6,13)
Consulting_Risk_fit2 <- reg_it(Consulting_Risk,4,6,13)
Consulting_Risk_fit3 <- reg_it(Consulting_Risk,5,6,13)

Consulting_Risk_results1 <- relweights(Consulting_Risk_fit1,'Consulting-Risk','Overall Experience', colorBlind4 )
Consulting_Risk_results1 

Consulting_Risk_results2 <- relweights(Consulting_Risk_fit2,'Consulting-Risk', 'Final Deliverables', colorBlind4 )
Consulting_Risk_results2 

Consulting_Risk_results3 <- relweights(Consulting_Risk_fit3,'Consulting-Risk','LTR', colorBlind4)
Consulting_Risk_results3 

### Consulting Technology

Consulting_Technology_fit1 <- reg_it(Consulting_Technology,3,6,13)
Consulting_Technology_fit2 <- reg_it(Consulting_Technology,4,6,13)
Consulting_Technology_fit3 <- reg_it(Consulting_Technology,5,6,13)

Consulting_Technology_results1 <- relweights(Consulting_Technology_fit1,'Consulting-Technology','Overall Experience', colorBlind5 )
Consulting_Technology_results1 

Consulting_Technology_results2 <- relweights(Consulting_Technology_fit2,'Consulting-Technology', 'Final Deliverables', colorBlind5 )
Consulting_Technology_results2 

Consulting_Technology_results3 <- relweights(Consulting_Technology_fit3,'Consulting-Technology','LTR', colorBlind5)
Consulting_Technology_results3 


## Without Satisfaction.of.frequency.of.conversations.with.RSM.professionals.outside.engagement


## Set up data

Audit3 <- Audit2[,1:11]
Audit_no_newteam3 <- Audit_no_newteam[,1:11]
Tax3 <- Tax2[,1:12]
Consulting3 <- Consulting2[,1:12]
Consulting_Risk3 <- Consulting_Risk[,1:12]
Consulting_Technology3 <- Consulting_Technology[,1:12]

##Audit w/o Sat

Audit_fit1b <- reg_it(Audit3,2,5,11)
Audit_fit2b <- reg_it(Audit3,3,5,11)
Audit_fit3b <- reg_it(Audit3,4,5,11)

Audit_results1b <- relweights(Audit_fit1b,'Audit','Overall Experience', colorBlind1 )
Audit_results1b 

Audit_results2b <- relweights(Audit_fit2b,'Audit', 'Final Deliverables', colorBlind1 )
Audit_results2b 

Audit_results3b <- relweights(Audit_fit3b,'Audit','LTR', colorBlind1)
Audit_results3b 

Audit_no_newteam_fit1b <- reg_it(Audit_no_newteam3 ,2,5,11)
Audit_no_newteam_fit2b <- reg_it(Audit_no_newteam3 ,3,5,11)
Audit_no_newteam_fit3b <- reg_it(Audit_no_newteam3 ,4,5,11)

Audit_resultsNNT1b <- relweights(Audit_no_newteam_fit1b,'Audit-w/o New Team Members','Overall Experience', colorBlind6)
Audit_resultsNNT1b 

Audit_resultsNNT2b <- relweights(Audit_no_newteam_fit2b,'Audit-w/o New Team Members','Final Deliverables', colorBlind6)
Audit_resultsNNT2b 

Audit_resultsNNT3b <- relweights(Audit_no_newteam_fit3b,'Audit-w/o New Team Members','LTR', colorBlind6)
Audit_resultsNNT3b 

### Tax

Tax_fit1b <- reg_it(Tax3,2,5,12)
Tax_fit2b <- reg_it(Tax3,3,5,12)
Tax_fit3b <- reg_it(Tax3,4,5,12)

Tax_results1b <- relweights(Tax_fit1b,'Tax','Overall Experience', colorBlind2 )
Tax_results1b 

Tax_results2b <- relweights(Tax_fit2b,'Tax', 'Final Deliverables', colorBlind2 )
Tax_results2b 

Tax_results3b <- relweights(Tax_fit3b,'Tax','LTR', colorBlind2)
Tax_results3b 

### Consulting

Consulting_fit1b <- reg_it(Consulting3,3,6,12)
Consulting_fit2b <- reg_it(Consulting3,4,6,12)
Consulting_fit3b <- reg_it(Consulting3,5,6,12)

Consulting_results1b <- relweights(Consulting_fit1b,'Consulting','Overall Experience', colorBlind3 )
Consulting_results1b 

Consulting_results2b <- relweights(Consulting_fit2b,'Consulting', 'Final Deliverables', colorBlind3 )
Consulting_results2b 

Consulting_results3b <- relweights(Consulting_fit3b,'Consulting','LTR', colorBlind3)
Consulting_results3b 

### Consulting Risk

Consulting_Risk_fit1b <- reg_it(Consulting_Risk3,3,6,12)
Consulting_Risk_fit2b <- reg_it(Consulting_Risk3,4,6,12)
Consulting_Risk_fit3b <- reg_it(Consulting_Risk3,5,6,12)

Consulting_Risk_results1b <- relweights(Consulting_Risk_fit1b,'Consulting-Risk','Overall Experience', colorBlind4 )
Consulting_Risk_results1b 

Consulting_Risk_results2b <- relweights(Consulting_Risk_fit2b,'Consulting-Risk', 'Final Deliverables', colorBlind4 )
Consulting_Risk_results2b 

Consulting_Risk_results3b <- relweights(Consulting_Risk_fit3b,'Consulting-Risk','LTR', colorBlind4)
Consulting_Risk_results3b 

### Consulting Technology

Consulting_Technology_fit1b <- reg_it(Consulting_Technology3,3,6,12)
Consulting_Technology_fit2b <- reg_it(Consulting_Technology3,4,6,12)
Consulting_Technology_fit3b <- reg_it(Consulting_Technology3,5,6,12)

Consulting_Technology_results1b <- relweights(Consulting_Technology_fit1b,'Consulting-Technology','Overall Experience', colorBlind5 )
Consulting_Technology_results1b 

Consulting_Technology_results2b <- relweights(Consulting_Technology_fit2b,'Consulting-Technology', 'Final Deliverables', colorBlind5 )
Consulting_Technology_results2b 

Consulting_Technology_results3b <- relweights(Consulting_Technology_fit3b,'Consulting-Technology','LTR', colorBlind5)
Consulting_Technology_results3b 



## Without Approx..frequency.of.conversations.with.RSM.professionals.outside.engagement 


## Set up data

Audit4 <- Audit2[,c(1:10,12)]
Audit_no_newteam4 <- Audit_no_newteam[,,c(1:10,12)]
Tax4 <- Tax2[,c(1:11,13)]
Consulting4 <- Consulting2[,c(1:11,13)]
Consulting_Risk4 <- Consulting_Risk[,c(1:11,13)]
Consulting_Technology4 <- Consulting_Technology[,c(1:11,13)]

##Audit w/o Sat other SAT

Audit_fit1c <- reg_it(Audit4,2,5,11)
Audit_fit2c <- reg_it(Audit4,3,5,11)
Audit_fit3c <- reg_it(Audit4,4,5,11)

Audit_results1c <- relweights(Audit_fit1c,'Audit','Overall Experience', colorBlind1 )
Audit_results1c 

Audit_results2c <- relweights(Audit_fit2c,'Audit', 'Final Deliverables', colorBlind1 )
Audit_results2c 

Audit_results3c <- relweights(Audit_fit3c,'Audit','LTR', colorBlind1)
Audit_results3c 

Audit_no_newteam_fit1c <- reg_it(Audit_no_newteam4 ,2,5,11)
Audit_no_newteam_fit2c <- reg_it(Audit_no_newteam4 ,3,5,11)
Audit_no_newteam_fit3c <- reg_it(Audit_no_newteam4 ,4,5,11)

Audit_resultsNNT1c <- relweights(Audit_no_newteam_fit1c,'Audit-w/o New Team Members','Overall Experience', colorBlind6)
Audit_resultsNNT1c 

Audit_resultsNNT2c <- relweights(Audit_no_newteam_fit2c,'Audit-w/o New Team Members','Final Deliverables', colorBlind6)
Audit_resultsNNT2c 

Audit_resultsNNT3c <- relweights(Audit_no_newteam_fit3c,'Audit-w/o New Team Members','LTR', colorBlind6)
Audit_resultsNNT3c 

### Tax

Tax_fit1c <- reg_it(Tax4,2,5,12)
Tax_fit2c <- reg_it(Tax4,3,5,12)
Tax_fit3c <- reg_it(Tax4,4,5,12)

Tax_results1c <- relweights(Tax_fit1c,'Tax','Overall Experience', colorBlind2 )
Tax_results1c 

Tax_results2c <- relweights(Tax_fit2c,'Tax', 'Final Deliverables', colorBlind2 )
Tax_results2c 

Tax_results3c <- relweights(Tax_fit3c,'Tax','LTR', colorBlind2)
Tax_results3c 

### Consulting

Consulting_fit1c <- reg_it(Consulting4,3,6,12)
Consulting_fit2c <- reg_it(Consulting4,4,6,12)
Consulting_fit3c <- reg_it(Consulting4,5,6,12)

Consulting_results1c <- relweights(Consulting_fit1c,'Consulting','Overall Experience', colorBlind3 )
Consulting_results1c 

Consulting_results2c <- relweights(Consulting_fit2c,'Consulting', 'Final Deliverables', colorBlind3 )
Consulting_results2c 

Consulting_results3c <- relweights(Consulting_fit3c,'Consulting','LTR', colorBlind3)
Consulting_results3c 

### Consulting Risk

Consulting_Risk_fit1c <- reg_it(Consulting_Risk4,3,6,12)
Consulting_Risk_fit2c <- reg_it(Consulting_Risk4,4,6,12)
Consulting_Risk_fit3c <- reg_it(Consulting_Risk4,5,6,12)

Consulting_Risk_results1c <- relweights(Consulting_Risk_fit1c,'Consulting-Risk','Overall Experience', colorBlind4 )
Consulting_Risk_results1c 

Consulting_Risk_results2c <- relweights(Consulting_Risk_fit2c,'Consulting-Risk', 'Final Deliverables', colorBlind4 )
Consulting_Risk_results2c 

Consulting_Risk_results3c <- relweights(Consulting_Risk_fit3c,'Consulting-Risk','LTR', colorBlind4)
Consulting_Risk_results3c 

### Consulting Technology

Consulting_Technology_fit1c <- reg_it(Consulting_Technology4,3,6,12)
Consulting_Technology_fit2c <- reg_it(Consulting_Technology4,4,6,12)
Consulting_Technology_fit3c <- reg_it(Consulting_Technology4,5,6,12)

Consulting_Technology_results1c <- relweights(Consulting_Technology_fit1c,'Consulting-Technology','Overall Experience', colorBlind5 )
Consulting_Technology_results1c 

Consulting_Technology_results2c <- relweights(Consulting_Technology_fit2c,'Consulting-Technology', 'Final Deliverables', colorBlind5 )
Consulting_Technology_results2c 

Consulting_Technology_results3c <- relweights(Consulting_Technology_fit3c,'Consulting-Technology','LTR', colorBlind5)
Consulting_Technology_results3c 

## end

library(officer)
library(rvg)

plots <- list(Audit_results1,Audit_results2,Audit_results3,
              Audit_resultsNNT1,Audit_resultsNNT2, Audit_resultsNNT3,
              Tax_results1, Tax_results2, Tax_results3,
              Consulting_results1, Consulting_results2, Consulting_results3,
              Consulting_Risk_results1,Consulting_Risk_results2,Consulting_Risk_results3,
              Consulting_Technology_results1,Consulting_Technology_results2,Consulting_Technology_results3,
              Audit_results1b,Audit_results2b,Audit_results3b,
              Audit_resultsNNT1b,Audit_resultsNNT2b, Audit_resultsNNT3b,
              Tax_results1b, Tax_results2b, Tax_results3b,
              Consulting_results1b, Consulting_results2b, Consulting_results3b,
              Consulting_Risk_results1b,Consulting_Risk_results2b,Consulting_Risk_results3b,
              Consulting_Technology_results1b,Consulting_Technology_results2b,Consulting_Technology_results3b,
              Audit_results1c,Audit_results2c,Audit_results3c,
              Audit_resultsNNT1c,Audit_resultsNNT2c, Audit_resultsNNT3c,
              Tax_results1c, Tax_results2c, Tax_results3c,
              Consulting_results1c, Consulting_results2c, Consulting_results3c,
              Consulting_Risk_results1c,Consulting_Risk_results2c,Consulting_Risk_results3c,
              Consulting_Technology_results1c,Consulting_Technology_results2c,Consulting_Technology_results3c 
              )

LOB <- c('Audit -', 'Audit w/o New Team -', 'Tax', 'Consulting -','Consulting-Risk -','Consulting-Technology -')
DV <- c("Overall Sat -","Deliverables Sat -","LTR -")
type <- c("Full Model", 'w/o Approx Conversations Freq', 'w/o Sat Freq Conversations')

labl <- NULL
ctr <- 0

for (k in 1:3) 
{
  for (i in 1:6)
  {
    for (j in 1:3)
    {
        ctr <- ctr+1
        labl[ctr] <- as.vector(paste(LOB[i],DV[j],type[k]))
        print(labl[ctr])
  
    }
  }
}
labl[2]
plots[2]



doc <- read_pptx()
ctr <- 0
for (p in plots) {
  ctr <- ctr+1
  plot_dml <- rvg::dml(ggobj = p )
  doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
  doc <- ph_with(x = doc, value = labl[ctr],
                 location = ph_location_type(type = "title") )
  doc <- ph_with(doc, plot_dml, location = ph_location_type(type = "body") )
}
print(doc, target = "my_plots3.pptx")

