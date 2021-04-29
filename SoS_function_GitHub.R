#Function to obtain the Sentinel species########
#sp_data is a matrix with species abundances as columns and samples as rows
#group_vect is a vector with with the pressure levels for each sample, keeping the same order than sp_data and with the same length than the number of row of sp_data
#Pressure levels has to be introduced as numbers. #The level 1 is the reference conditions, ideally has to be the No Effort level. If this is not possible then it can be Low effort. 
#Besito is the sesnitivity value of the species. Is a data frame with the species name and sensitivity value (see table 2 Serrano et al., unpublished)
#The lines from 10 to 89 compute SIMPER (adapted from Farriols et al 2015) 

SoS <- function(sp_data, group_vect, Besito, table){
  
  ####combinations within group using SIMPER
  
  table_list <- list()
  
  group <- as.character(unique(group_vect))
  cutoff <- 90
  
  for (g in group){
    
    dbg<- na.omit(sp_data[group_vect== g, ])
    take <- t(combn(1:nrow(dbg), 2))
    take
    md<-numeric(ncol(dbg))
    
    
    me<-numeric(ncol(dbg))
    contr <- matrix(nrow = nrow(take), ncol = ncol(dbg+2))
    for (j in 1:nrow(take)) {
      for (i in 1:ncol(dbg)) {	
        md[i] <- 2*min(dbg[take[j , 1],i], dbg[take[j, 2],i])
        me <- dbg[take[j, 1], ] + dbg[take[j, 2], ]
        contr[j, ] <- 100*(md / sum(me))
      }
    }
    
    contr<-data.frame(contr)
    colnames(contr)<-c(names(dbg))
    ndbg<-ncol(dbg)
    x1<-colMeans(contr)
    x2<-apply(contr, 2, sd)
    df1<-data.frame(x1)
    df2<-data.frame(x2)
    x<-c(names(dbg))
    y1<-df1$x1
    y2<-df2$x2
    z2<-data.frame(Species=rep(NA, ndbg), Similarity=rep(NA, ndbg), SD=rep(NA, ndbg))
    z2$Species<-x
    z2$Similarity<-y1
    z2$SD<-y2
    
    ##ordenation of species by similarity
    z2$Similarity <- as.numeric(as.character(z2$Similarity))
    db<-z2[order(z2$Similarity,decreasing = TRUE),]
    b<-sum(db$Similarity)
    ##contribution to mean similarity
    db$contribution<-((100*db$Similarity)/b)
    
    ##calculus of accumulated mean similarity \% and cut according to 90\%
    
    d<-0
    for (i in 1:ncol(dbg)){
      d<-(db$contribution[i]+d)
      db$acum[i]<-d
    }
    
    db$Similarity <- round(db$Similarity,2)
    db$SD <- round(db$SD,2)
    db$contribution <- round(db$contribution,2)
    db$acum <- round(db$acum,2)
    db$acum[1]<- ifelse(db$acum[1]>90, 90, db$acum[1])
    
    r<-db[which(db$acum<=cutoff & db$contribution>0),]
    
    Abundance0<-apply(dbg,2,mean)
    Ab_SD0<-apply(dbg,2,sd)
    Ab_SD10<-data.frame(Ab_SD0)
    Ab0<-data.frame(Abundance0)
    x0<-c(names(dbg))
    y0<-round(Ab0$Abundance0,2)
    w0<-round(Ab_SD10$Ab_SD0,2)	
    z0<-cbind(x0,y0,w0)
    rb0<-merge(z0, r, by.x = "x0", by.y = "Species")
    table<-rb0[order(rb0$Similarity,decreasing = TRUE),]
    names(table)<-c("Species", "Av.Abund", "SD.Abund", "Av.Si", "SD.Si", "Contr", "Cum")
    table$Group <- g
    table_list[[g]] <- table
    
  }
  
  table_final <- do.call(rbind, table_list) #END OF SIMPER
  
  #We select species from Group 1 
  Sp_Group1 <- table_final[table_final$Group=="1",]
  Sp_Group1_vect <- Sp_Group1$Species
  names(Besito) <- c("Species", "BESITO")
  Besito_SIMPER <- unique(Besito[Besito$Species%in%Sp_Group1_vect,])
  
  #We compute frecuencyes
  SpNames <- colnames(sp_data)
  sp_data$Filter <- group_vect 
  LowEffort_Data <- sp_data[sp_data$Filter=="1",]
  Freq <- vector()
  Per <- vector()
  for( i in 1:(ncol(LowEffort_Data)-1)){
    ColOnly <- LowEffort_Data[,i]
    ColOnly_Bi <- ifelse(ColOnly>0,1,0)
    Freq[i] <- sum(ColOnly_Bi)
    Per[i] <- (sum(ColOnly_Bi)/nrow(LowEffort_Data))*100
  }
  
  MyFreqMatrix <- cbind.data.frame(SpNames, Freq, Per)
  Th <- ifelse(round(length(group_vect[group_vect==1])/10)<=2,2,round(length(group_vect[group_vect==1])/10))
  MyFreqMatrix <- MyFreqMatrix[MyFreqMatrix$Freq>=Th,]
  names(MyFreqMatrix) <- c("Species", "Freq", "Per")
  
  MyFreqMatrixWithBesito <- unique(merge(MyFreqMatrix, Besito, by="Species"))
  MyFreqMatrixWithBesito <- MyFreqMatrixWithBesito[order(MyFreqMatrixWithBesito$Freq, decreasing=TRUE),]
  
  #We select species with a Sensitivity value of 5
  
  Sens5_Simper <-  Besito_SIMPER[Besito_SIMPER$BESITO==5,]
  SIM_Freq_Sens <- Sens5_Simper
  
  Mns <- "10 Species reached after include species with a sensitive of 5 from SIMPER"
  
  if (length(SIM_Freq_Sens)>=10) {
    print(Mns)
    return(SIM_Freq_Sens)
  }
  
  Sens5_Freq <-  MyFreqMatrixWithBesito[MyFreqMatrixWithBesito$BESITO==5,]
  Sp_Vect5 <- unique(c(as.character(Sens5_Simper$Species), as.character(Sens5_Freq$Species)))
  MyCut<- ifelse(length(Sp_Vect5)>=10,10,length(Sp_Vect5))
  MyCutSp <- Sp_Vect5[MyCut]
  MyCut_df <- Sens5_Freq[Sens5_Freq$Species==MyCutSp,]
  MyCut_Frq <- MyCut_df$Freq
  MyMat <-   Sens5_Freq[Sens5_Freq$Freq>=MyCut_Frq,]
  SIM_Freq_Sens <- unique(c(as.character(Sens5_Simper$Species), as.character(MyMat$Species)))
  SIM_Freq_Sens <- SIM_Freq_Sens[SIM_Freq_Sens!="integer(0)"]
  
  Mns <- "10 Species reached after include species with a sensitive of 5 ordered by Frecuency"
  
  if (length(SIM_Freq_Sens)>=10) {
    print(Mns)
    return(SIM_Freq_Sens)
  }
  
  #We select species with a Sensitivity value of 4
  
  Sens4_Simper <-  Besito_SIMPER[Besito_SIMPER$BESITO==4,]
  SIM_Freq_Sens <- unique(c(SIM_Freq_Sens,Sens4_Simper$Species))
  
  Mns <- "10 Species reached after include species with a sensitive of 4 from SIMPER"
  
  if (length(SIM_Freq_Sens)>=10) {
    print(Mns)
    return(SIM_Freq_Sens)
  }
  
  Sens4_Freq <-  MyFreqMatrixWithBesito[MyFreqMatrixWithBesito$BESITO==4,]
  Sp_Vect4 <- unique(c(as.character(SIM_Freq_Sens), as.character(Sens4_Simper$Species),as.character(Sens4_Freq$Species)))
  MyCut <- ifelse(length(Sp_Vect4)>=10,10,length(Sp_Vect4))
  MyCutSp <- Sp_Vect4[MyCut]
  MyCut_df <- Sens4_Freq[Sens4_Freq$Species==MyCutSp,]
  MyCut_Frq <- MyCut_df$Freq
  MyMat <-   Sens4_Freq[Sens4_Freq$Freq>=MyCut_Frq,]
  SIM_Freq_Sens <- unique(c(as.character(SIM_Freq_Sens),as.character(Sens4_Simper$Species), as.character(MyMat$Species)))
  SIM_Freq_Sens <- SIM_Freq_Sens[SIM_Freq_Sens!="integer(0)"]
  
  Mns <- "10 Species reached after include species with a sensitive of 4 ordered by Frecuency"
  
  if (length(SIM_Freq_Sens)>=10) {
    print(Mns)
    return(SIM_Freq_Sens)
  }
  
  #We select species with a Sensitivity value of 3
  
  Sens3_Simper <-  Besito_SIMPER[Besito_SIMPER$BESITO==3,]
  SIM_Freq_Sens <- unique(c(SIM_Freq_Sens,Sens3_Simper$Species))
  
  Mns <- "10 Species reached after include species with a sensitive of 3 from SIMPER"
  
  if (length(SIM_Freq_Sens)>=10) {
    print(Mns)
    return(SIM_Freq_Sens)
  }
  
  Sens3_Freq <-  MyFreqMatrixWithBesito[MyFreqMatrixWithBesito$BESITO==3,]
  Sp_Vect3 <- unique(c(as.character(SIM_Freq_Sens), as.character(Sens3_Simper$Species),as.character(Sens3_Freq$Species)))
  MyCut<- ifelse(length(Sp_Vect3)>=10,10,length(Sp_Vect3))
  MyCutSp <- Sp_Vect3[MyCut]
  MyCut_df <- Sens3_Freq[Sens3_Freq$Species==MyCutSp,]
  MyCut_Frq <- MyCut_df$Freq
  MyMat <-   Sens3_Freq[Sens3_Freq$Freq>=MyCut_Frq,]
  SIM_Freq_Sens <- unique(c(as.character(SIM_Freq_Sens),as.character(Sens3_Simper$Species), as.character(MyMat$Species)))
  SIM_Freq_Sens <- SIM_Freq_Sens[SIM_Freq_Sens!="integer(0)"]
  
  Mns <- "5 Species reached after include species with a sensitive of 3 ordered by Frecuency"
  
  if (length(SIM_Freq_Sens)>=5) {
    print(Mns)
    return(SIM_Freq_Sens)
  }
  
  Sens2_Simper <-  Besito_SIMPER[Besito_SIMPER$BESITO==2,]
  SIM_Freq_Sens <- unique(c(SIM_Freq_Sens,Sens2_Simper$Species))
  
  Mns <- "5 Species reached after include species with a sensitive of 2 from SIMPER"
  
  if (length(SIM_Freq_Sens)>=5) {
    print(Mns)
    return(SIM_Freq_Sens)
  }
  
  Sens2_Freq <-  MyFreqMatrixWithBesito[MyFreqMatrixWithBesito$BESITO==2,]
  Sp_Vect2 <- unique(c(as.character(SIM_Freq_Sens), as.character(Sens2_Simper$Species),as.character(Sens2_Freq$Species)))
  MyCut<- ifelse(length(Sp_Vect2)>=5,5,length(Sp_Vect2))
  MyCutSp <- Sp_Vect2[MyCut]
  MyCut_df <- Sens2_Freq[Sens2_Freq$Species==MyCutSp,]
  MyCut_Frq <- MyCut_df$Freq
  MyMat <-   Sens2_Freq[Sens2_Freq$Freq>=MyCut_Frq,]
  SIM_Freq_Sens2 <- unique(c(as.character(SIM_Freq_Sens),as.character(Sens3_Simper$Species), as.character(MyMat$Species)))
  SIM_Freq_Sens2 <- SIM_Freq_Sens2[SIM_Freq_Sens2!="integer(0)"]
  
  Mns <- "5 Species reached after include species with a sensitive of 2 ordered by Frecuency"
  if (length(SIM_Freq_Sens2)>=5) {
    print(Mns)
    return(SIM_Freq_Sens2)
  }
  
  
  Mns <- "Careful: Loop finished without reach minimum level of species"
  print(Mns)
  return(SIM_Freq_Sens2)
  
  
}
