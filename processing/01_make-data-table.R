# -------------------------------------------------------------------------
#' Make a data.frame containing the data for analysis
#' Combine data from different sources to make one data.frame with all data
#' 
# -------------------------------------------------------------------------


# MIR 3D7 -----------------------------------------------------------------

# data extracted from Table S2:
inoculum_MIR <- c(1e6,1e7,1e8,1e9)
n_wells_MIR <- c(3,3,3,3)
n_pos_wells_MIR <- c(0,2,2,3)

# combine data:
data.MIR.3D7 <- data.frame(inoculum=inoculum_MIR,n_wells=n_wells_MIR,n_pos_wells=n_pos_wells_MIR,
                           group=rep("MIR_3D7",length(inoculum_MIR)))

# remove unused variables:
rm(inoculum_MIR,n_wells_MIR,n_pos_wells_MIR)


# MIR Dd2 -----------------------------------------------------------------

# data extracted from Table S2:
inoculum_MIR_table_Dd2 <- c(1e4,1e5,2e5,1e6,1e7,2e7,1e8,1e9)
n_wells_MIR_table_Dd2 <- c(3,6,96,9,9,5,3,3)
n_pos_wells_MIR_table_Dd2 <- c(0,0,1,3,9,3,3,3)

# combine data:
data.MIR.table.Dd2 <- data.frame(inoculum=inoculum_MIR_table_Dd2,n_wells=n_wells_MIR_table_Dd2,n_pos_wells=n_pos_wells_MIR_table_Dd2,
                                 group=rep("MIR_Dd2",length(inoculum_MIR_table_Dd2)))

# remove unused variables:
rm(inoculum_MIR_table_Dd2,n_wells_MIR_table_Dd2,n_pos_wells_MIR_table_Dd2)


# NSG mouse data ----------------------------------------------------------

# parameters for mouse experiments:
lod <- 0.02 # limit of detection (in percent parasitemia)
volume <- 2 # mouse blood volume in ml
hematocrit <- 70/100 # mouse hematocrit 70%

# load data from 2018 (use all data of treated mice):
data_pd2018 <- read_excel("raw-data/NSG-data.xlsx",sheet="PD2018",col_names = FALSE)
data_pd2018 <- t(data_pd2018)
colnames(data_pd2018) <- data_pd2018[1,]
data_pd2018 <- apply(data_pd2018,2,function(x) as.numeric(as.character(x)))
data_pd2018 <- data_pd2018[,c(19,which(grepl("M5717",colnames(data_pd2018))))]
data_pd2018 <- data_pd2018[rowSums(is.na(data_pd2018))<ncol(data_pd2018)-1,] # remove NA rows
data_pd2018 <- as.data.frame(data_pd2018)
colnames(data_pd2018) <- c(colnames(data_pd2018)[c(1:6)],"6 mpk M5717_M2")
data_NSG_2018 <- as.data.frame(lapply(data_pd2018,function(x) ifelse(x>lod,x,lod)))

# load data from 2019 (use data for the 5 mice treated with the 12mg/kg dose):
data_pd2019 <- read_excel("raw-data/NSG-data.xlsx",sheet="PD2019",col_names = FALSE)
data_pd2019 <- t(data_pd2019)
colnames(data_pd2019) <- data_pd2019[1,]
data_pd2019 <- data_pd2019[2:dim(data_pd2019)[1],2:dim(data_pd2019)[2]]
colnames(data_pd2019)[9] <- "time"
data_pd2019 <- apply(data_pd2019,2,function(x) as.numeric(as.character(x)))
data_pd2019 <- as.data.frame(data_pd2019)
data_pd2019 <- data_pd2019[rowSums(is.na(data_pd2019[,1:8]))<8,] # discard rows without any percent parasitemia measurements
data_pd2019$time[22] <- 41 # add missing time information
time <- sort(unique(data_pd2019$time))
M_2019_12mpk_1_1 <- rep(NA,length(time))
M_2019_12mpk_1_1[which(time%in%data_pd2019$time)] <- data_pd2019$M5717_12mpk_cage1_1
M_2019_12mpk_1_2 <- rep(NA,length(time))
M_2019_12mpk_1_2[which(time%in%data_pd2019$time)] <- data_pd2019$M5717_12mpk_cage1_2
M_2019_12mpk_2_3 <- rep(NA,length(time))
M_2019_12mpk_2_3[which(time%in%data_pd2019$time)] <- data_pd2019$M5717_12mpk_cage2_3
M_2019_12mpk_2_4 <- rep(NA,length(time))
M_2019_12mpk_2_4[which(time%in%data_pd2019$time)] <- data_pd2019$M5717_12mpk_cage2_4
M_2019_12mpk_2_5 <- rep(NA,length(time))
M_2019_12mpk_2_5[which(time%in%data_pd2019$time)] <- data_pd2019$M5717_12mpk_cage2_5
data_NSG_2019 <- as.data.frame(cbind(time,M_2019_12mpk_1_1,M_2019_12mpk_1_2,
                                 M_2019_12mpk_2_3,M_2019_12mpk_2_4,M_2019_12mpk_2_5))
data_NSG_2019 <- as.data.frame(lapply(data_NSG_2019,function(x) ifelse(x>lod,x,lod)))

# load data from 2020 (use all data of treated mice):
data_pd2020 <- read_excel("raw-data/NSG-data.xlsx",sheet="PD2020",col_names = FALSE)
data_pd2020 <- data_pd2020[c(1:16),c(2:4,6,7)] # omit control mice, only treated mice of interest
colnames(data_pd2020) <- c("Days post infection","hematocrit_M1","hematocrit_M2",
                           "12mpk_M5717_M1","12mpk_M5717_M2")
data_pd2020 <- apply(data_pd2020,2,function(x) as.numeric(as.character(x)))
data_pd2020 <- as.data.frame(data_pd2020[rowSums(is.na(data_pd2020))<ncol(data_pd2020)-1,]) # remove NA rows
data_NSG_2020 <- as.data.frame(lapply(data_pd2020[,c(1,4,5)],function(x) ifelse(x>lod,x,lod)))

# data from all years for mice treated with 12 mg/kg and mice treated with 30mg/kg from 2018:
time <- c(min(c(data_NSG_2018$Days.post.infection,data_NSG_2019$time,data_NSG_2020$Days.post.infection)):
            max(c(data_NSG_2018$Days.post.infection,data_NSG_2019$time,data_NSG_2020$Days.post.infection)))
data_NSG <- data.frame(time=time)
data_NSG <- dplyr::left_join(data_NSG,data_NSG_2018[,c(1,which(grepl("12",names(data_NSG_2018))))],by=c('time'='Days.post.infection'))
data_NSG <- dplyr::left_join(data_NSG,data_NSG_2018[,c(1,which(grepl("30",names(data_NSG_2018))))],by=c('time'='Days.post.infection'))
names(data_NSG)[which(names(data_NSG)%in%names(data_NSG_2018))] <- paste(names(data_NSG)[which(names(data_NSG)%in%names(data_NSG_2018))],"_2018",sep="")
data_NSG <- dplyr::left_join(data_NSG,data_NSG_2019[,c(1,which(grepl("12",names(data_NSG_2019))))],by=c('time'='time'))
data_NSG <- dplyr::left_join(data_NSG,data_NSG_2020[,c(1,which(grepl("12",names(data_NSG_2020))))],by=c('time'='Days.post.infection'))
names(data_NSG)[which(names(data_NSG)%in%names(data_NSG_2020))] <- paste(names(data_NSG)[which(names(data_NSG)%in%names(data_NSG_2020))],"_2020",sep="")
data_NSG <- as.data.frame(data_NSG[rowSums(is.na(data_NSG))<ncol(data_NSG)-1,]) # remove NA rows

# LDA-style data:
inoculum_NSG <- t(data_NSG[1,c(2:ncol(data_NSG))]*(hematocrit*volume)/(100*hRBC_vol)) 
# for the mice with known hematocrit, use the known hematocrit:
inoculum_NSG[which(grepl("2020",rownames(inoculum_NSG)) & grepl("M1",rownames(inoculum_NSG)))] <- 
  data_NSG[1,which(grepl("2020",names(data_NSG)) & grepl("M1",names(data_NSG)))]*(data_pd2020$hematocrit_M1[1]/100*volume)/(100*hRBC_vol)
inoculum_NSG[which(grepl("2020",rownames(inoculum_NSG)) & grepl("M2",rownames(inoculum_NSG)))] <- 
  data_NSG[1,which(grepl("2020",names(data_NSG)) & grepl("M2",names(data_NSG)))]*(data_pd2020$hematocrit_M2[1]/100*volume)/(100*hRBC_vol)

n_wells_NSG <- rep(1,length(inoculum_NSG)) # one well (mouse) each
n_pos_wells_NSG <- rep(1,length(inoculum_NSG)) # all positive (i.e. with resistant parasites) but "M_2019_12mpk_2_4", "M_2019_12mpk_2_5" and "X12mpk_M5717_M2_2020"
n_pos_wells_NSG[rownames(inoculum_NSG)%in%c("M_2019_12mpk_2_4","M_2019_12mpk_2_5","X12mpk_M5717_M2_2020")] <- 0

# combine data for all mice:
data.NSG <- data.frame(inoculum=as.numeric(inoculum_NSG),n_wells=n_wells_NSG,n_pos_wells=n_pos_wells_NSG,
                                 group=rep("NSG",length(inoculum_NSG)))

if(nrow(data.NSG)>length(unique(data.NSG$inoculum))){
  inoculum <- sort(unique(data.NSG$inoculum))
  data.NSG <- data.frame(inoculum=inoculum,
                         n_wells=apply(as.data.frame(inoculum),1,function(x){sum(data.NSG$n_wells[data.NSG$inoculum==x])}),
                         n_pos_wells=apply(as.data.frame(inoculum),1,function(x){sum(data.NSG$n_pos_wells[data.NSG$inoculum==x])}),
                         group=rep("NSG",length(inoculum)))
}

# remove unused variables:
rm(lod,volume,hematocrit,data_pd2018,data_NSG_2018,data_pd2019,time,M_2019_12mpk_1_1,M_2019_12mpk_1_2,M_2019_12mpk_2_3,M_2019_12mpk_2_4,M_2019_12mpk_2_5,
   data_NSG_2019,data_pd2020,data_NSG_2020,data_NSG,inoculum_NSG,n_wells_NSG,n_pos_wells_NSG,inoculum)


# VIS data ----------------------------------------------------------------

# parameters for VIS data:
blood_vol <- 5000 # blood volume in ml

# VIS data with treatment parasite number for all treatment cohorts:
data_coh1 <- read_excel("raw-data/VIS-data.xlsx",sheet="Cohort 1",col_names = TRUE)
data_coh2 <- read_excel("raw-data/VIS-data.xlsx",sheet="Cohort 2",col_names = TRUE)
data_coh3 <- read_excel("raw-data/VIS-data.xlsx",sheet="Cohort 3",col_names = TRUE)
data_all_coh <- rbind(data_coh1,data_coh2,data_coh3)
recr <- as.data.frame(cbind(c(unique(data_coh1$ID),unique(data_coh2$ID),unique(data_coh3$ID)),
                            c(0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,0,0,0,0)))
names(recr) <- c("ID","recrudescence")
dosis <- c(rep(400,8),rep(150,6),rep(800,8))

# compute the parasite number at treatment:
tr_par <- vector(mode = "numeric", length = dim(recr)[1]) # parasite number for each individual at treatment
for(i in 1:dim(recr)[1]){
  id_tmp <- recr$ID[i]
  tr_par[i] <- round(data_all_coh$`Geometric mean parasitaemia (parasites/ml)`[data_all_coh$ID==id_tmp & data_all_coh$Time==unique(data_all_coh$Treatment.Day)])*blood_vol
}
recr <- cbind(recr,tr_par,dosis)

# make LDA data using the parasite number at treatment:
inoculum_VIS <- recr$tr_par
n_wells_VIS <- matrix(1,dim(recr)[1],1)
n_pos_wells_VIS <- as.numeric(as.character(recr$recrudescence))

# combine data:
data.VIS <- data.frame(inoculum=as.numeric(inoculum_VIS),n_wells=n_wells_VIS,n_pos_wells=n_pos_wells_VIS,
                       group=rep("VIS",length(inoculum_VIS)))

# remove unused variables:
rm(blood_vol,data_coh1,data_coh2,data_coh3,data_all_coh,recr,dosis,tr_par,i,id_tmp,inoculum_VIS,n_wells_VIS,n_pos_wells_VIS)


# field isolates regrowth data --------------------------------------------

# culture data:
volume <- 6 # culture volume 6mL
hematocrit <- 0.04 # 4% hematocrit

# load data:
field.data <- read.csv("raw-data/Regrowth-field-data.csv")

RBCs_num_per_Culture <- (volume*hematocrit)/hRBC_vol # cellvolume (=volume*hematocrit) / volume of a human RBC
field.data.paras.num <- field.data$Parasitemia[field.data$Time==field.data$TreatTime]*RBCs_num_per_Culture/100 # parasite number

inoculum_field <- sort(unique(field.data.paras.num))
n_wells_field <- sapply(c(1:length(inoculum_field)),function(x){sum(field.data.paras.num==inoculum_field[x])})
n_pos_wells_field <- sapply(c(1:length(inoculum_field)),function(x){sum(field.data$Resistance[which(field.data.paras.num==inoculum_field[x])])})

# combine data:
data.field <- data.frame(inoculum=inoculum_field,n_wells=n_wells_field,n_pos_wells=n_pos_wells_field,group=rep("field",length(inoculum_field)))

# PMR between first parasitemia measurement and treatment (to use for stochastic simulations of experiments):
# par.init <- unlist(lapply(unique(field.data$TreatTime),function(x){mean(field.data$Parasitemia[field.data$TreatTime==x & field.data$Time==2]*RBCs_num_per_Culture/100)}))
# time.treatment <- unique(field.data$TreatTime)-2
# gen.treatment <- round(time.treatment*24/40)
# par.treatment <- unlist(lapply(unique(field.data$TreatTime),function(x){mean(field.data$Parasitemia[field.data$TreatTime==x & field.data$Time==x]*RBCs_num_per_Culture/100)}))
# distr <- unlist(lapply(unique(field.data$TreatTime),function(x){length(field.data$Parasitemia[field.data$TreatTime==x & field.data$Time==2]*RBCs_num_per_Culture/100)}))
# mean.pmr <- round(mean(rep(round((par.treatment/par.init)^(1/(round((time.treatment-2)*24/40)))),unlist(lapply(unique(field.data$TreatTime),function(x){length(unique(field.data$isolates[field.data$TreatTime==x]))})))))

# remove unused variables:
rm(volume,hematocrit,field.data,RBCs_num_per_Culture,field.data.paras.num,inoculum_field,n_wells_field,n_pos_wells_field)


# field data summary ------------------------------------------------------

# create field summary data:
data.field.sum <- data.frame(inoculum=mean(data.field$inoculum[data.field$group=="field"]),n_wells=sum(data.field$n_wells[data.field$group=="field"]),
                             n_pos_wells=sum(data.field$n_pos_wells[data.field$group=="field"]),group="field_sum")


# Regrowth 3D7 data -------------------------------------------------------

# 3D7 experiments by Laurent Dembele and his group:
volume <- 6 # culture volume 6mL
hematocrit <- 0.04 # 4% hematocrit

# load data:
ld_3D7 <- read_excel("raw-data/Regrowth-3D7-data.xlsx")
ld_3D7 <- ld_3D7[rowSums(is.na(ld_3D7))<ncol(ld_3D7),] # remove NA rows
ld_3D7.paras <- ld_3D7[c(2:nrow(ld_3D7)),6] # parasitemias in the data
ld_3D7.paras <- as.numeric(ld_3D7.paras$...6)
ld_3D7.mutants <- c(0,0,0,0,1,1,0,1,0)

RBCs_num_per_Culture <- (volume*hematocrit)/hRBC_vol # cellvolume (=volume*hematocrit) / volume of a human RBC
ld_3D7.paras.num <- ld_3D7.paras*RBCs_num_per_Culture/100 # parasite number
# mean number of parasites at inoculation on day 0: mean(as.numeric(as.matrix(ld_3D7[c(2:10),2]))*RBCs_num_per_Culture/100)
# mean number of parasites at treatment on day 4: mean(ld_3D7.paras.num)

inoculum_ld_3D7 <- sort(unique(ld_3D7.paras.num))
n_wells_ld_3D7 <- sapply(c(1:length(inoculum_ld_3D7)),function(x){sum(ld_3D7.paras.num==inoculum_ld_3D7[x])})
n_pos_wells_ld_3D7 <- sapply(c(1:length(inoculum_ld_3D7)),function(x){sum(ld_3D7.mutants[ld_3D7.paras.num==inoculum_ld_3D7[x]])})

# combine data:
data.ld_3d7 <- data.frame(inoculum=inoculum_ld_3D7,n_wells=n_wells_ld_3D7,n_pos_wells=n_pos_wells_ld_3D7,group=rep("ld_3D7",length(inoculum_ld_3D7)))

# remove unused variables:
rm(volume,hematocrit,ld_3D7,ld_3D7.paras,ld_3D7.mutants,RBCs_num_per_Culture,inoculum_ld_3D7,n_wells_ld_3D7,n_pos_wells_ld_3D7)

# LD's 52x 3D7 data:
inoculum_ld_3D7_x52 <- mean(ld_3D7.paras.num) # assume the parasitemia at treatment is the mean parasitemia of the 9 culutures with known parasitemia
n_wells_ld_3D7_x52 <- 52 # 52 cultures
n_pos_wells_ld_3D7_x52 <- 3 # 3 cultures for 3D7 had resistant parasites

# combine data:
data.ld_3D7_x52 <- data.frame(inoculum=inoculum_ld_3D7_x52,n_wells=n_wells_ld_3D7_x52,n_pos_wells=n_pos_wells_ld_3D7_x52,
                              group=rep("ld_3D7_x52",length(inoculum_ld_3D7_x52)))

# remove unused variables:
rm(inoculum_ld_3D7_x52,ld_3D7.paras.num,n_wells_ld_3D7_x52,n_pos_wells_ld_3D7_x52)


# combine all data into one data.frame ------------------------------------

data.all <- rbind(data.MIR.3D7,data.MIR.table.Dd2,data.NSG,data.VIS,data.field,data.field.sum,data.ld_3d7,data.ld_3D7_x52)

# remove unused variables:
rm(data.MIR.3D7,data.MIR.table.Dd2,data.NSG,data.VIS,data.field,data.field.sum,data.ld_3d7,data.ld_3D7_x52)



