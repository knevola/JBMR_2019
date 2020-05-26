# Z score function
z_score <- function(age,bone,data, age_block){
  data_0 <- data
  ages<-seq(from = min(data_0[[age]]), to = max(data_0[[age]]), by = age_block)
  for (i in 1:(length(ages)-1)){
    data_0$Age_Group[data_0[[age]] >= ages[i] & data_0[[age]] < ages[i+1]] <- paste("Group",i)
  }
  data_0$Age_Group[data_0[[age]] >= ages[i+1]] <- paste("Group",(i+1))
  means<-aggregate(data_0[[bone]],by = list(data_0$Age_Group),FUN = mean)
  names(means)<- c("Age_Group", "mean_BMD")
  sds<-aggregate(data_0[[bone]],by = list(data_0$Age_Group),FUN = sd)
  names(sds)<-c("Age_Group", "sd_BMD")
  
  zscore_tools <- merge(means, sds, by = "Age_Group")
  
  data_1 <- merge(data_0, zscore_tools, by = "Age_Group")
  
  data_1$zscore <- (data_1[[bone]] - data_1$mean_BMD)/data_1$sd_BMD
  return(data_1)
}


