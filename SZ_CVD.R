library(data.table)
library(dplyr)
setwd("/Users/zhoujiayi/Library/CloudStorage/OneDrive-TheUniversityofHongKong-Connect")

###### 0. load dataset ########
dx <- readRDS("raw_data/schizophrenia/Scz_1993_202103_DX_all.rds")  #diagnosis,sex,birth.date
death <- readRDS("raw_data/schizophrenia/Scz_1993_202103_death_all.rds")  #death date, death cause
rx <- readRDS("raw_data/schizophrenia/Scz_1999_202103_RX_all.rds") #prescription

###### 1. prevalent schizophrenia diagnosed between 2005-2015 #######
dx <- setDT(readRDS("raw_data/schizophrenia/Scz_1993_202103_DX_all.rds"))
schizophrenia <- dx[as.Date(Reference.Date)>=as.Date("2005-01-01")&as.Date(Reference.Date)<=as.Date("2015-12-31")&
                      grepl("^295",All.Diagnosis.Code..ICD9.,ignore.case=T)&!grepl("^295.5",All.Diagnosis.Code..ICD9.,ignore.case=T)]
schizophrenia <- schizophrenia[,.SD[which.min(as.Date(Reference.Date))],by="Reference.Key"] ##49,540
schizophrenia <- schizophrenia[,.(Reference.Key,index.date=Reference.Date)]
saveRDS(schizophrenia,"schizophrenia_CVD/data/schizophrenia_2005.2015.rds")

###### 2.flowchart ########
###### exclusion1: CVD/ASCVD on or before the index date #####
schizophrenia <- setDT(readRDS("schizophrenia_CVD/data/schizophrenia_2005.2015.rds"))
schizophrenia[,index.date:=index.date+365] #### index date as 1 year after first diagnosis

source("schizophrenia_CVD/Coding_Book.R") # Load coding book
dx <- setDT(readRDS("raw_data/schizophrenia/Scz_1993_202103_DX_all.rds"))
dx <- merge(dx,schizophrenia,all.x=T,by="Reference.Key")[!is.na(index.date)]
covar_dz <- function(ct, name,icd) {
  hx <- dx[(as.Date(Reference.Date)<=as.Date(index.date))&grepl(icd,All.Diagnosis.Code..ICD9.,ignore.case=T),unique(Reference.Key)]
  ct[, c(paste0("dx.b4.",name)) := as.numeric(Reference.Key %in% hx)]
  cat(length(hx),"\n")
  return(ct)
}
for(i in 1:5) {
  cat(paste0(out.list[i]), "...")
  schizophrenia <- covar_dz(schizophrenia, out.list[i], cat.out[i])
}
schizophrenia <- covar_dz(schizophrenia, "cvd", CVD.code)

schizophrenia <- schizophrenia[dx.b4.cvd==0] ###46311
saveRDS(schizophrenia,"schizophrenia_CVD/data/schizophrenia_exclude.cvd.rds")

###### exclusion2: no birth date or sex #####
schizophrenia <- setDT(readRDS("schizophrenia_CVD/data/schizophrenia_exclude.cvd.rds"))
dx <- setDT(readRDS("raw_data/schizophrenia/Scz_1993_202103_DX_all.rds"))
birth <- merge(schizophrenia,dx[,.(Reference.Key,Sex,Date.of.Birth..yyyy.mm.dd.,Exact.date.of.birth)],all.x=T,by="Reference.Key")
schi.date.birth <- birth %>% group_by(Reference.Key) %>% 
  mutate(dob=case_when(Exact.date.of.birth=="Y" ~ last(Date.of.Birth..yyyy.mm.dd.),
                       Exact.date.of.birth=="N" ~ last(Date.of.Birth..yyyy.mm.dd.)))
setorder(schi.date.birth,Reference.Key,-Exact.date.of.birth)
schizophrenia <- setDT(schi.date.birth)[,.SD[1],by="Reference.Key"][,.(Reference.Key,index.date,Sex,dob)] 
schizophrenia <- schizophrenia[(!is.na(dob))] #46311-46286=25
schizophrenia <- schizophrenia[Sex!="U"] #46286-1=46285
###### exclusion3: age<=25 #####
schizophrenia[,age:=as.numeric(round((index.date-dob)/365.25,0))]
schizophrenia <- schizophrenia[age>=26] #46285-4586=41699
###### exclusion4: die on or before the index date or dead person with missing death date#####
death <- setDT(readRDS("raw_data/schizophrenia/Scz_1993_202103_death_all.rds"))
death <- death[!is.na(`Date of Registered Death`)|`Registered Death Identifier`!="Alive"][,.(Reference.Key=`Reference Key`,death.date=`Date of Registered Death`,`Exact date of death`,`Registered Death Identifier`)]
excl <- death[!is.na(`Registered Death Identifier`)&is.na(death.date)]
schizophrenia.death <- merge(schizophrenia,death,all.x = T,by="Reference.Key")
schizophrenia.death <- schizophrenia.death[is.na(death.date)|(!is.na(death.date)&index.date<as.Date(death.date))]  #41699-1414=40285
schizophrenia.death <- schizophrenia.death[!Reference.Key %in% excl$Reference.Key]  #40285-2=40283
saveRDS(schizophrenia.death,"schizophrenia_CVD/data/schizophrenia_age.death.rds")

###### exclusion5: anticoagulant or antiplatelet within 6 months before the index date #####
schizophrenia <- setDT(readRDS("schizophrenia_CVD/data/schizophrenia_age.death.rds"))
rx <- setDT(readRDS("raw_data/schizophrenia/Scz_1999_202103_RX_all.rds"))
rx <- rx %>% filter((!is.na(Drug.Frequency))&(!is.na(Dispensing.Date))) %>% 
  filter(!(Drug.Frequency == "ONCE"|Drug.Frequency == "AT ONCE"|grepl("HOURS",Drug.Frequency))) 
source("schizophrenia_CVD/Coding_Book.R") # Load coding book
anticoagulant <- rx[grepl(Anticoagulant.antiplatelet.code,Therapeutic.Classification.BNF,ignore.case=T)][,Dispensing.Date:=as.Date(Dispensing.Date)]
anticoagulant.drug <- merge(schizophrenia,anticoagulant[,.(Reference.Key,Dispensing.Date)],all.x = T,by="Reference.Key")
anticoagulant.drug <- anticoagulant.drug[(index.date-Dispensing.Date)<=180&(index.date-Dispensing.Date)>=0,unique(Reference.Key)]
schizophrenia <- schizophrenia[!Reference.Key %in% anticoagulant.drug]  #40283-856=39427
saveRDS(schizophrenia,"schizophrenia_CVD/data/schizophrenia_patient.list.rds")

###### 3.predictors #########
schizophrenia <- setDT(readRDS("schizophrenia_CVD/data/schizophrenia_patient.list.rds"))
####3.1 Duration of schizophrenia####
dx <- setDT(readRDS("raw_data/schizophrenia/Scz_1993_202103_DX_all.rds"))
dur.sz <- dx[grepl("^295",All.Diagnosis.Code..ICD9.,ignore.case=T)&!grepl("^295.5",All.Diagnosis.Code..ICD9.,ignore.case=T)]
dur.sz <- merge(schizophrenia,dur.sz[,.(Reference.Key,first.sz=Reference.Date)],all.x = T,by="Reference.Key")
dur.sz <- dur.sz[,.SD[which.min(as.Date(first.sz))],by="Reference.Key"]
schizophrenia <- dur.sz[,dur.sz:=as.numeric((index.date-first.sz)/365.25)]
####3.2 comorbidities####
dx <- merge(dx,schizophrenia[,.(Reference.Key,index.date)],all.x=T,by="Reference.Key")
source("schizophrenia_CVD/Coding_Book.R") # Load coding book
covar_dz <- function(ct, name,icd) {
  hx <- dx[(as.Date(Reference.Date)<=as.Date(index.date))&grepl(icd,All.Diagnosis.Code..ICD9.,ignore.case=T),unique(Reference.Key)]
  ct[, c(paste0("dx.b4.",name)) := as.numeric(Reference.Key %in% hx)]
  cat(length(hx),"\n")
  return(ct)
}
for(i in 1:13) {
  cat(paste0(dx.list[i]), "...")
  schizophrenia <- covar_dz(schizophrenia, dx.list[i], cat.dx[i])
}
saveRDS(schizophrenia,"schizophrenia_CVD/data/schizophrenia_comorbidities.rds")

####3.3 prescription####
schizophrenia <- setDT(readRDS("schizophrenia_CVD/data/schizophrenia_comorbidities.rds"))
source("schizophrenia_CVD/Coding_Book.R")
rx <- setDT(readRDS("raw_data/schizophrenia/Scz_1999_202103_RX_all.rds"))
rx <- rx %>% filter((!is.na(Drug.Frequency))&(!is.na(Dispensing.Date))) %>% 
  filter(!(Drug.Frequency == "ONCE"|Drug.Frequency == "AT ONCE"|grepl("HOURS",Drug.Frequency))) 
rx.sz <- merge(schizophrenia[,.(Reference.Key,index.date)],rx[,.(Reference.Key,Dispensing.Date=as.Date(Dispensing.Date),Dispensing.Duration,Dispensing.Duration.Unit,bnf=Therapeutic.Classification.BNF,drug.code=Drug.Item.Code)],all.x=T,by="Reference.Key")
table(rx.sz$Dispensing.Duration.Unit)
table(is.na(rx.sz$Dispensing.Duration))
rx.sz[Dispensing.Duration.Unit=="Day(s)",dispensing.end:=as.Date(Dispensing.Date)+Dispensing.Duration]
rx.sz[Dispensing.Duration.Unit=="Week(s)",dispensing.end:=as.Date(Dispensing.Date)+Dispensing.Duration*7]
rx.sz[Dispensing.Duration.Unit=="Month(s)",dispensing.end:=as.Date(Dispensing.Date)+Dispensing.Duration*30]
rx.sz[is.na(Dispensing.Duration),dispensing.end:=Dispensing.Date]

# rx.sz <- rx.sz[index.date-Dispensing.Date>=0]
rx.sz.6m <- rx.sz[index.date-Dispensing.Date>=0&index.date-Dispensing.Date<=180]
table(rx.sz$Dispensing.Duration.Unit)

filterPrescription <- function(data,data_6m,d.name,d.code,target,code.type="dcode"){
  # Prescribing of antipsychotics was considered if patients filled at least two consecutive prescriptions of 
  # antipsychotic medication in the same class (i.e., FGA or SGA) 
  # or a particular antipsychotic agent (e.g., haloperidol) within 6 months preceding the index date.
  tic <- Sys.time()
  consec.df <- data_6m %>%
    { if (code.type == "bnf") filter(.,grepl(d.code,bnf))
      else if (code.type == "dcode") filter(.,grepl(d.code,drug.code))
      else print("The code.type assignment is wrong.") }%>%
    arrange(Reference.Key, Dispensing.Date) %>%
    group_by(Reference.Key) %>%
    mutate(D.end.previous = lag(dispensing.end,n=1L)) %>% # Shift the start date up by 1 ob
    ungroup() %>%
    mutate(interval = Dispensing.Date - D.end.previous) %>% # Calculate the interval between two obs (Start date(n) - End date(n-1))
    filter((!is.na(interval))&(interval<=30)) %>% # Interval between the end date (previous prescription) and the start date (latter prescription) no longer than 30 days
    distinct(Reference.Key)
  target[[d.name]] <- with(target, ifelse(Reference.Key %in% consec.df$Reference.Key, 1, 0))
  toc <- Sys.time()
  print (toc-tic)
  
  # Duration (Years)
  tic <- Sys.time()
  years.df <- data %>%
    { if (code.type == "bnf") filter(.,grepl(d.code,bnf))
      else if (code.type == "dcode") filter(.,grepl(d.code,drug.code))
      else print("The code.type assignment is wrong.") }%>%
    filter(Dispensing.Date <= index.date)
  
  years.df <- years.df %>%
    arrange(Reference.Key, Dispensing.Date) %>%
    mutate(D.Date.y = as.numeric(format(as.Date(Dispensing.Date),"%Y"))) %>%
    distinct(Reference.Key, D.Date.y, .keep_all = T) %>% 
    group_by(Reference.Key) %>%
    mutate(Years.sum = n()) %>%
    ungroup() %>%
    distinct(Reference.Key,.keep_all = T) %>% 
    select(Reference.Key,Years.sum) %>%
    setnames("Years.sum", paste(d.name,".years.sum",sep=""))
  target = merge(target,years.df,by="Reference.Key",all.x = T)
  toc <- Sys.time()
  print (toc-tic)
  
  return(target)
}

schizophrenia <- filterPrescription(rx.sz,rx.sz.6m,"Antihyper",Antihyper.code,schizophrenia,code.type="bnf")
schizophrenia <- filterPrescription(rx.sz,rx.sz.6m,"Antipsy",Antipsy.code,schizophrenia,code.type="bnf")
schizophrenia <- filterPrescription(rx.sz,rx.sz.6m,"Antidepre",Antidepre.code,schizophrenia,code.type="bnf")

for (i in 1:16){
  print (i)
  schizophrenia <- filterPrescription(rx.sz,rx.sz.6m,cat.fga.list[i],cat.fga[i],schizophrenia)
}
schizophrenia <- filterPrescription(rx.sz,rx.sz.6m,"FGA",cat.fga.all,schizophrenia)
for (i in 1:11){
  print (i)
  schizophrenia <- filterPrescription(rx.sz,rx.sz.6m,cat.sga.list[i],cat.sga[i],schizophrenia)
}
schizophrenia <- filterPrescription(rx.sz,rx.sz.6m,"SGA",cat.sga.all,schizophrenia)

schizophrenia <- as.data.frame(schizophrenia)
continous.cols <- grep("sz.duration|years.sum",names(schizophrenia),value=T)
schizophrenia[continous.cols] <- lapply(schizophrenia[continous.cols],as.numeric)
schizophrenia[continous.cols] <- mutate_all(schizophrenia[continous.cols], ~coalesce(.,0))
saveRDS(schizophrenia,"schizophrenia_CVD/data/schizophrenia_drug.rds")

##### 4. outcome ######
schizophrenia <- setDT(readRDS("schizophrenia_CVD/data/schizophrenia_drug.rds"))
dx <- setDT(readRDS("raw_data/schizophrenia/Scz_1993_202103_DX_all.rds"))
dx <- merge(dx,schizophrenia,all.y=T,by="Reference.Key")[!is.na(index.date)]
source("schizophrenia_CVD/Coding_Book.R") # Load coding book
covar_dz <- function(ct, name,icd) {
  hx <- dx[(as.Date(Reference.Date)>as.Date(index.date))&grepl(icd,All.Diagnosis.Code..ICD9.,ignore.case=T)]
  hx <- hx[,.SD[which.min(Reference.Date)],by="Reference.Key"]
  ct[, c(paste0("oc.",name)) := as.numeric(Reference.Key %in% hx$Reference.Key)]
  ct <- merge(ct,hx[,.(Reference.Key,Reference.Date)],all.x=T,by="Reference.Key")
  colnames(ct)[colnames(ct) == "Reference.Date"] <- paste0("oc.date.",name)
  cat(nrow(hx),"\n")
  return(ct)
}
for(i in 1:5) {
  cat(paste0(out.list[i]), "...")
  schizophrenia <- covar_dz(schizophrenia, out.list[i], cat.out[i])
}
schizophrenia <- covar_dz(schizophrenia, "cvd", CVD.code)

death <- readRDS("raw_data/schizophrenia/Scz_1993_202103_death_all.rds")  #death date, death cause
death <- death[!is.na(`Date of Registered Death`)|`Registered Death Identifier`!="Alive"]
death <- death[,.(Reference.Key=`Reference Key`,death.date=`Date of Registered Death`,death.icd=`Death Cause (Main Cause)`,death.icd2=`Death Cause (Supplementary Cause)`)]

schizophrenia <- merge(schizophrenia,death[,.(Reference.Key,death.icd)],all.x = TRUE,by="Reference.Key")
schizophrenia$death.date <- as.Date(schizophrenia$death.date)
death_dz <- function(ct, name,icd) {
  hx <- schizophrenia[(!is.na(death.date))&(as.Date(death.date)<=as.Date("2021-03-31"))&grepl(icd,death.icd,ignore.case=T),unique(Reference.Key)]
  ct[, c(paste0("death.",name)) := as.numeric(Reference.Key %in% hx)]
  ct[get(paste0("death.",name))==1,c(paste0("death.date.",name)) :=death.date]
  cat(length(hx),"\n")
  return(ct)
}
for(i in 1:4) {
  cat(paste0(death.list[i]), "...")
  schizophrenia <- death_dz(schizophrenia, death.list[i], death.code[i])
}

schizophrenia[death.ascvd==1,date_ASCVD:=death.date.ascvd]
schizophrenia[oc.MI==1|oc.Stroke==1,date_ASCVD:=pmin(oc.date.Stroke,oc.date.MI,death.date.ascvd,death.date,na.rm=T)]
summary(schizophrenia$date_ASCVD)

schizophrenia[death.stroke==1,date_stroke:=death.date.ascvd]
schizophrenia[oc.Stroke==1,date_stroke:=pmin(oc.date.Stroke,death.date.stroke,death.date,na.rm=T)]
schizophrenia[oc.MI==1,date_nonfatal.mi:=oc.date.MI]
schizophrenia[death.chd==1,date_chd:=pmin(death.date.chd,death.date,na.rm=T)]

check <- schizophrenia[outcome.stroke==1&outcome.nonfatal.mi==1] %>% select(grep("date|death|oc|outcome",names(schizophrenia),value=T))

schizophrenia[date_chd!=date_ASCVD,date_chd:=NA]
schizophrenia[date_stroke!=date_ASCVD,date_stroke:=NA]
schizophrenia[date_nonfatal.mi!=date_ASCVD,date_nonfatal.mi:=NA]

schizophrenia[,`:=`(date_CVD=pmin(oc.date.cvd,death.date.cvd,na.rm=T))]
schizophrenia[oc.cvd==1,date_CVD:=pmin(date_CVD,death.date,na.rm=T)]

library(lubridate)
schizophrenia$Final.date.10yr <- as.Date(schizophrenia$index.date) %m+% years(10)
OUTCOMES <- c("CVD","ASCVD","stroke","nonfatal.mi","chd")
for(oc in OUTCOMES) {
  schizophrenia[, c(paste0("censor.date.", oc)) := pmin(get(paste0("date_",oc)), death.date, as.Date("2021-03-31"),Final.date.10yr, na.rm=T)]
  schizophrenia[, c(paste0("time.to.censor.", oc)) := as.numeric(get(paste0("censor.date.", oc)) - index.date)]
  schizophrenia[, c(paste0("outcome.", oc)) := as.numeric(!is.na(get(paste0("date_",oc))) &!is.na(get(paste0("censor.date.",oc))) & 
                                                            get(paste0("censor.date.", oc))==get(paste0("date_",oc)) & 
                                                            get(paste0("date_",oc))>=index.date)]
}

saveRDS(schizophrenia,"schizophrenia_CVD/data/final.dataset.rds")


#### 5. table one ########
cohort <- setDT(readRDS("schizophrenia_CVD/data/final.dataset.rds"))
cohort <- cohort %>%
  mutate(hyper.combo = dx.b4.hyper + Antihyper*2) %>% # 0: 30813; 1: 504; 2: 6400; 3: 1710
  mutate(FGA.SGA = FGA + SGA*2) # 0: 8303; 1: 15655; 2: 12134; 3: 3158
cohort$hyper.combo[cohort$hyper.combo==3] = 2
var <- names(cohort)
var <- var[!var %in% c("Reference.Key","dob","`Exact date of death`","`Registered Death Identifier`","death.icd")]
var.factor <- var[!var %in% c("age","index.date","death.date","first.sz","dur.sz",grep("sum|date|time.to.censor",names(cohort),value=T))]
cohort <- as.data.frame(cohort)
cohort[var.factor] <- lapply(cohort[var.factor],factor)
library(tableone)
table <- CreateTableOne(vars = var,data = setDT(cohort),test=T)
table1 <- print(table)
write.csv(table1, file = "schizophrenia_CVD/data/table1_all.characteristics.csv")

table.sex <- CreateTableOne(vars = var,strata = "Sex",data = setDT(cohort),test=T)
table1.sex <- print(table.sex,test=T, smd=T, dropEqual=T, noSpaces=T)
write.csv(table1.sex, file = "schizophrenia_CVD/data/table1_by.sex.csv")

table.oc <- CreateTableOne(vars = var,strata = "outcome.ASCVD",data = setDT(cohort),test=T)
table1.oc <- print(table.oc,test=T, smd=T, dropEqual=T, noSpaces=T)
write.csv(table1.oc, file = "schizophrenia_CVD/data/table1_by.outcome.csv")

##### 6. Cox ######
library(survival)
library(survminer)
library(StepReg)
library(sjPlot)
cohort <- setDT(readRDS("schizophrenia_CVD/data/final.dataset.rds"))
cohort <- cohort %>%
  mutate(hyper.combo = dx.b4.hyper + Antihyper*2) %>% # 0: 30703; 1: 478; 2: 6389; 3: 1680
  mutate(FGA.SGA = FGA + SGA*2) # 0: 8303; 1: 15655; 2: 12134; 3: 3158
cohort$hyper.combo[cohort$hyper.combo==3] = 2
cohort$Calendar.year <- as.numeric(format(cohort$index.date, "%Y"))

var <- names(cohort)
var <- var[!var %in% c("Reference.Key","dob","`Exact date of death`","`Registered Death Identifier`","death.icd")]
var.factor <- var[!var %in% c("age","index.date","death.date","first.sz","dur.sz","Calendar.year",grep("sum|date|time.to.censor",names(cohort),value=T))]
cohort <- as.data.frame(cohort)
cohort[var.factor] <- lapply(cohort[var.factor],factor)
table(cohort$outcome.CVD)
cohort$outcome.CVD <- with(cohort, ifelse(outcome.CVD==0, 1, 2))
table(cohort$outcome.CVD)
cohort$outcome.ASCVD <- with(cohort, ifelse(outcome.ASCVD==0, 1, 2))

set.seed(321)
sample <- sample(c(TRUE, FALSE), nrow(cohort), replace=TRUE, prob=c(0.7,0.3))
train <- cohort[sample,]
test <- cohort[!sample,]

predictors <- c("hyper.combo","dx.b4.af","dx.b4.diabetes","dx.b4.cancer","dx.b4.dyslipidemia","dx.b4.COPD","dx.b4.depression",
                "dx.b4.substance.misuse","dx.b4.anxiety.disorder","dx.b4.personality","dx.b4.alcohol.misuse","Antidepre",
                "FGA.SGA","Calendar.year","Antidepre.years.sum","FGA.years.sum","SGA.years.sum")

### correlation matrix ###
cor.m <- train[,c("Sex","age",predictors)]
cor.m$Sex <- as.numeric(cor.m$Sex)
library(Hmisc)
cor <- rcorr(as.matrix(cor.m))$r
cor.p <- rcorr(as.matrix(cor.m))$P
write.csv(cor, file = "schizophrenia_CVD/data/correlation.train.csv")

### age and sex-adjusted univariate model
options(scipen = 200)

hr_function <- function(survival.time,outcome){
  base <- coxph(as.formula(paste0("Surv(time =", survival.time,", event =", outcome,") ~ Sex+age")), data = train)
  hr <- cbind(exp(cbind(coef(base),confint(base))),melt(c(round(summary(base)[['coefficients']][,'Pr(>|z|)'],3))))
for(p in predictors) {
  m <- coxph(as.formula(paste0("Surv(time =", survival.time,", event =", outcome,") ~ Sex+age+",p)), data = train)
  hr <- rbind(hr,cbind(exp(cbind(coef(m),confint(m))),melt(c(round(summary(m)[['coefficients']][,'Pr(>|z|)'],3))))[-(1:2),])
}
  colnames(hr) <- c("HR","lower","upper","p")
  return(hr)
}
hr_ascvd <- hr_function("time.to.censor.ASCVD","outcome.ASCVD")
hr_cvd <- hr_function("time.to.censor.CVD","outcome.CVD")

## backward elimination by BIC ##
###### ascvd
back.ascvd <- stepwiseCox(formula= as.formula(paste0("Surv(time = time.to.censor.ASCVD, event = outcome.ASCVD) ~ Sex+age+",paste(predictors,collapse = "+"))),
            data=train,
            selection="backward",
            select="SBC",
            method="efron")
ascvd.train <- coxph(as.formula(paste0("Surv(time = time.to.censor.ASCVD, event = outcome.ASCVD) ~", 
                                       paste(c(t(back.ascvd$`Selected Varaibles`)),collapse = "+"))), data = train)
summary(ascvd.train)
cbind(exp(cbind(coef(ascvd.train),confint(ascvd.train))),melt(c(round(summary(ascvd.train)[['coefficients']][,'Pr(>|z|)'],3))))
cox.zph(ascvd.train)
par(mfrow=c(2,3))
plot(cox.zph(ascvd.train))

###### cvd 
back.CVD <- stepwiseCox(formula= as.formula(paste0("Surv(time = time.to.censor.CVD, event = outcome.CVD) ~ Sex+age+",paste(predictors,collapse = "+"))),
                        data=train,
                        selection="backward",
                        select="SBC",
                        method="efron")
CVD.train <- coxph(as.formula(paste0("Surv(time = time.to.censor.CVD, event = outcome.CVD) ~", 
                                       paste(c(t(back.CVD$`Selected Varaibles`)),collapse = "+"))), data = train)
cbind(exp(cbind(coef(CVD.train),confint(CVD.train))),melt(c(round(summary(CVD.train)[['coefficients']][,'Pr(>|z|)'],3))))
# tab_model(CVD.train, collapse.ci = F, show.ci = 0.95)  
cox.zph(CVD.train)
par(mfrow=c(2,4))
plot(cox.zph(CVD.train))

## backward elimination by AIC (sensitivity analysis) ##
###### ascvd
back.ascvd <- stepwiseCox(formula= as.formula(paste0("Surv(time = time.to.censor.ASCVD, event = outcome.ASCVD) ~ Sex+age+",paste(predictors,collapse = "+"))),
                          data=train,
                          selection="backward",
                          select="AIC",
                          method="efron")
ascvd.train <- coxph(as.formula(paste0("Surv(time = time.to.censor.ASCVD, event = outcome.ASCVD) ~", 
                                       paste(c(t(back.ascvd$`Selected Varaibles`)),collapse = "+"))), data = train)
summary(ascvd.train)
cbind(exp(cbind(coef(ascvd.train),confint(ascvd.train))),melt(c(round(summary(ascvd.train)[['coefficients']][,'Pr(>|z|)'],3))))
cox.zph(ascvd.train)
par(mfrow=c(2,3))
plot(cox.zph(ascvd.train))

###### cvd 
back.CVD <- stepwiseCox(formula= as.formula(paste0("Surv(time = time.to.censor.CVD, event = outcome.CVD) ~ Sex+age+",paste(predictors,collapse = "+"))),
                        data=train,
                        selection="backward",
                        select="AIC",
                        method="efron")
CVD.train <- coxph(as.formula(paste0("Surv(time = time.to.censor.CVD, event = outcome.CVD) ~", 
                                     paste(c(t(back.CVD$`Selected Varaibles`)),collapse = "+"))), data = train)
cbind(exp(cbind(coef(CVD.train),confint(CVD.train))),melt(c(round(summary(CVD.train)[['coefficients']][,'Pr(>|z|)'],3))))
# tab_model(CVD.train, collapse.ci = F, show.ci = 0.95)  
cox.zph(CVD.train)
par(mfrow=c(2,4))
plot(cox.zph(CVD.train))


#### 7. performance ######
## Format the output by (C-statistics, D-statistics, High-risk/ Low-risk patients developing `OUTCOME`) 
Cuno <- function(data, indices){
  new.data <- data[indices,]
  cuno <- survival::concordance(ascvd.train, timewt="n/G2",newdata=new.data)
  return(cuno$concordance)
}
format_output <- function(fit.model,new.data){
  h0 = basehaz(fit.model, centered=T) 
  S0 = exp(-tail(h0$hazard,1)) # C
  preds <- predict(fit.model, newdata = new.data, type ="risk", se.fit=T)
  cvd.result <- new.data[,c("Sex","event")]
  cvd.result$Risk <- 1-S0^preds$fit
  pred.event <- predict(fit.model, newdata = new.data, type ="expected")
  cvd.result$Expected <- pred.event
  OE_ratio <- sum(new.data$event)/sum(cvd.result$Expected)
  SE <- sqrt(1/sum(new.data$event) + 1/sum(cvd.result$Expected))
  lower_CI <- OE_ratio - 1.96*SE
  upper_CI <- OE_ratio + 1.96*SE
  cat("Total O/E ratio", round(OE_ratio,2),"(",round(lower_CI,2),"–",round(upper_CI,2),")\n")
  C.stats <- concordance(fit.model,newdata=new.data,timewt="n/G2")
  boot <- boot(data=new.data, statistic=Cuno, R = 1000)
  C.stats.ci <- boot.ci(boot, type="norm")
  cat("C: ",round(C.stats$concordance,2),"(",round(C.stats.ci$normal[2],2),"-",round(C.stats.ci$normal[3],2),")\n")
  
  D.stats <- royston(fit.model,newdata=new.data)
  cat("D: ",round(D.stats[1],2),"(",round(D.stats[1]-D.stats[2]*1.96,2),"–",round(D.stats[1]+D.stats[2]*1.96,2),")\n")
  
  cvd.result.high <- cvd.result[(cvd.result$Risk>0.2),]
  ratio.high <- round(sum(cvd.result.high$event)*100/dim(cvd.result.high)[1],1)
  
  cvd.result.low <- cvd.result[(cvd.result$Risk<=0.2),]
  ratio.low <- round(sum(cvd.result.low$event)*100/dim(cvd.result.low)[1],1)
  cat("High risk", sum(cvd.result.high$event),"/",dim(cvd.result.high)[1],"(",ratio.high,")\n")
  cat("Low risk", sum(cvd.result.low$event),"/",dim(cvd.result.low)[1],"(",ratio.low,")\n")
  
  for(sex in c("M","F")){
    cat("===================",sex,"===================\n")
    
    C.stats <- concordance(fit.model,newdata=new.data[new.data$Sex==sex,],timewt="n/G2")
    boot <- boot(data=new.data[new.data$Sex==sex,], statistic=Cuno, R = 1000)
    C.stats.ci <- boot.ci(boot, type="norm")
    cat("C: ",round(C.stats$concordance,2),"(",round(C.stats.ci$normal[2],2),"-",round(C.stats.ci$normal[3],2),")\n")
    
    D.stats <- royston(fit.model,newdata=new.data[new.data$Sex==sex,])
    cat("D: ",round(D.stats[1],2),"(",round(D.stats[1]-D.stats[2]*1.96,2),"–",round(D.stats[1]+D.stats[2]*1.96,2),")\n")
    
    cvd.result.high <- cvd.result[(cvd.result$Risk>0.2)&(cvd.result$Sex==sex),]
    ratio.high <- round(sum(cvd.result.high$event)*100/dim(cvd.result.high)[1],1)
    
    cvd.result.low <- cvd.result[(cvd.result$Risk<=0.2)&(cvd.result$Sex==sex),]
    ratio.low <- round(sum(cvd.result.low$event)*100/dim(cvd.result.low)[1],1)
    cat("High risk", sum(cvd.result.high$event),"/",dim(cvd.result.high)[1],"(",ratio.high,")\n")
    cat("Low risk", sum(cvd.result.low$event),"/",dim(cvd.result.low)[1],"(",ratio.low,")\n")
    
    OE_ratio <- sum(new.data[new.data$Sex==sex,"event"])/sum(cvd.result[cvd.result$Sex==sex,"Expected"])
    SE <- sqrt(1/sum(new.data[new.data$Sex==sex,"event"]) + 1/sum(cvd.result[cvd.result$Sex==sex,"Expected"]))
    # lower_CI <- exp(log(OE_ratio) - 1.96*SE)
    # upper_CI <- exp(log(OE_ratio) + 1.96*SE)
    lower_CI <- OE_ratio - 1.96*SE
    upper_CI <- OE_ratio + 1.96*SE
    cat("O/E ratio", round(OE_ratio,2),"(",round(lower_CI,2),"–",round(upper_CI,2),")\n")
  }
  return(cvd.result)
}
library(boot)
test.ascvd <- test
test.ascvd$outcome.ASCVD <- test.ascvd$outcome.ASCVD-1
test.ascvd$event <- test.ascvd$outcome.ASCVD
ascvd.result = format_output(ascvd.train,test.ascvd)
h0 = basehaz(ascvd.train, centered=T) 
S0 = exp(-tail(h0$hazard,1))

test.cvd <- test
test.cvd$outcome.CVD <- test.cvd$outcome.CVD-1
test.cvd$event <- test.cvd$outcome.CVD
CVD.result = format_output(CVD.train,test.cvd)


#### calibration plot #######
library(predtools)
p <- calibration_plot(data=cvd.result, obs="event", pred="Risk", group = "Sex", y_lim = c(0, 0.4), x_lim = c(0, 0.4)) 
p$calibration_plot + scale_x_continuous(name="10-y Predicted CVD Risk, %",
                                        breaks=c(0.0,0.1,0.2,0.3,0.4),
                                        labels=c(0,10,20,30,40)) +
  scale_y_continuous(name="Observed CVD Risk, %",
                     breaks=c(0.0,0.1,0.2),
                     labels=c(0,10,20)) +
  geom_point(size=5) +
  theme(axis.text=element_text(face="plain",size=16),text=element_text(face="bold", size=18))

ggsave("./schizophrenia_CVD/data/HR.cvd.jpeg",dpi = 1200)

p <- calibration_plot(data=ascvd.result, obs="event", pred="Risk", group = "Sex", y_lim = c(0, 0.4), x_lim = c(0, 0.4)) 
p$calibration_plot + scale_x_continuous(name="10-y Predicted ASCVD Risk, %",
                                        breaks=c(0.0,0.1,0.2,0.3,0.4),
                                        labels=c(0,10,20,30,40)) +
  scale_y_continuous(name="Observed ASCVD Risk, %",
                     breaks=c(0.0,0.1,0.2),
                     labels=c(0,10,20)) +
  geom_point(size=5) +
  theme(axis.text=element_text(face="plain",size=16),text=element_text(face="bold", size=18))

ggsave("./schizophrenia_CVD/data/HR.ascvd.jpeg",dpi = 1200)


#### secondary analysis ##########
### model of predictors selected by BIC + added predictors
## ascvd
names(train)
add <- c("Chlorpromazine","Chlorpromazine.years.sum","Deanxit","Deanxit.years.sum","Flupenthixol","Flupenthixol.years.sum",
         "Haloperidol","Haloperidol.years.sum","Pericyazine","Pericyazine.years.sum","Perphenazine","Perphenazine.years.sum",
         "Pimozide","Pimozide.years.sum","Sulpiride","Sulpiride.years.sum","Thioridazine","Thioridazine.years.sum",
         "Thiothixene","Thiothixene.years.sum","Trifluoperazine","Trifluoperazine.years.sum","Zuclopenthixol","Zuclopenthixol.years.sum",
         "Clozapine","Clozapine.years.sum","Amisulpride","Amisulpride.years.sum","Aripiprazole","Aripiprazole.years.sum",
         "Olanzapine","Olanzapine.years.sum","Paliperidone","Paliperidone.years.sum","Quetiapine","Quetiapine.years.sum",
         "Risperidone","Risperidone.years.sum","Sertindole","Sertindole.years.sum","Ziprasidone","Ziprasidone.years.sum")
back.ascvd.bic <- stepwiseCox(formula= as.formula(paste0("Surv(time = time.to.censor.ASCVD, event = outcome.ASCVD) ~ Sex+age+",paste(predictors,collapse = "+"))),
                          data=train,
                          selection="backward",
                          select="SBC",
                          method="efron")
ascvd.train.bic <- coxph(as.formula(paste0("Surv(time = time.to.censor.ASCVD, event = outcome.ASCVD) ~", 
                                       paste(c(t(back.ascvd.bic$`Selected Varaibles`),add),collapse = "+"))), data = train)
cbind(exp(cbind(coef(ascvd.train.bic),confint(ascvd.train.bic))),melt(c(round(summary(ascvd.train.bic)[['coefficients']][,'Pr(>|z|)'],3))))

## cvd 
back.CVD.bic <- stepwiseCox(formula= as.formula(paste0("Surv(time = time.to.censor.CVD, event = outcome.CVD) ~ Sex+age+",paste(predictors,collapse = "+"))),
                        data=train,
                        selection="backward",
                        select="SBC",
                        method="efron")
CVD.train.bic <- coxph(as.formula(paste0("Surv(time = time.to.censor.CVD, event = outcome.CVD) ~", 
                                     paste(c(t(back.CVD.bic$`Selected Varaibles`),add),collapse = "+"))), data = train)
cbind(exp(cbind(coef(CVD.train.bic),confint(CVD.train.bic))),melt(c(round(summary(CVD.train.bic)[['coefficients']][,'Pr(>|z|)'],3))))

### model of predictors selected by AIC + added predictors
back.ascvd.aic <- stepwiseCox(formula= as.formula(paste0("Surv(time = time.to.censor.ASCVD, event = outcome.ASCVD) ~ Sex+age+",paste(predictors,collapse = "+"))),
                              data=train,
                              selection="backward",
                              select="AIC",
                              method="efron")
ascvd.train.aic <- coxph(as.formula(paste0("Surv(time = time.to.censor.ASCVD, event = outcome.ASCVD) ~", 
                                           paste(c(t(back.ascvd.aic$`Selected Varaibles`),add),collapse = "+"))), data = train)
cbind(exp(cbind(coef(ascvd.train.aic),confint(ascvd.train.aic))),melt(c(round(summary(ascvd.train.aic)[['coefficients']][,'Pr(>|z|)'],3))))

## cvd 
back.CVD.aic <- stepwiseCox(formula= as.formula(paste0("Surv(time = time.to.censor.CVD, event = outcome.CVD) ~ Sex+age+",paste(predictors,collapse = "+"))),
                            data=train,
                            selection="backward",
                            select="AIC",
                            method="efron")
CVD.train.aic <- coxph(as.formula(paste0("Surv(time = time.to.censor.CVD, event = outcome.CVD) ~", 
                                         paste(c(t(back.CVD.aic$`Selected Varaibles`),add),collapse = "+"))), data = train)
cbind(exp(cbind(coef(CVD.train.aic),confint(CVD.train.aic))),melt(c(round(summary(CVD.train.aic)[['coefficients']][,'Pr(>|z|)'],3))))
