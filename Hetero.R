setwd('~/Desktop/data')
library(broom)


#IC: zz500 id = 1  IF: hs300 id = 2 IH: sz50 id = 3
{
  y <- read.csv("~/Desktop/data/data_new.csv",header = TRUE)
  rf <- read.csv("~/Desktop/data/rf.csv",header = TRUE)
  index_sz50 <- read.csv("~/Desktop/data/index_sz50.csv",header = TRUE)
  index_zz500 <- read.csv("~/Desktop/data/index_zz500.csv",header = TRUE)
  index_hs300 <- read.csv("~/Desktop/data/index_hs300.csv",header = TRUE)
  divd <- read.csv("~/Desktop/data/divd.csv",header = TRUE)
}


# data sourcing
{
  #0 data sourcing
  
  data_IC <- y[grep("IC",y$Agmtcd),]  #IC: zz500
  data_IF <- y[grep("IF",y$Agmtcd),]  #IF: hs300
  data_IH <- y[grep("IH",y$Agmtcd),]  #IH: sz50
  
  
  
  
  #1. calculate mean of spot return, risk free return
  
  #1.1 risk free return
  summary(rf)
  mean_rf1m <- mean (rf_1m)
  mean_rf2m <- mean (rf_2m)
  mean_rf3m <- mean (rf_3m)
  mean_rf6m <- mean (rf_6m)
  var_rf1m <- var(rf_1m)
  sd_rf1m <- var_rf1m ^ 0.5
  var_rf2m <- var(rf_2m)
  sd_rf2m <- var_rf2m ^ 0.5
  var_rf3m <- var(rf_3m)
  sd_rf3m <- var_rf3m ^ 0.5
  var_rf6m <- var(rf_6m)
  sd_rf6m <- var_rf6m ^ 0.5
  paste("[mean]rf_1m:", mean_rf1m, "rf_2m:", mean_rf2m, "rf_3m:", mean_rf3m, "rf_6m:", mean_rf6m)
  paste("[var]rf_1m:", var_rf1m, "rf_2m:", var_rf2m, "rf_3m:", var_rf3m, "rf_6m:", var_rf6m)
  paste("[s.d.]rf_1m:", sd_rf1m, "rf_2m:", sd_rf2m, "rf_3m:", sd_rf3m, "rf_6m:", sd_rf6m)
  
}


# Using moving average return method
{
  # zz500 part
  {
    #2. estimates of expected return (1m, 2m, 3m, 6m)
    
    
    #2.1 convert rf into continuously compounding
    #rf$rf_1m <- log(1+rf$rf_1m/252/100)
    #rf$rf_2m <- log(1+rf$rf_2m/252/100)
    #rf$rf_3m <- log(1+rf$rf_3m/252/100)
    #转化到其他dataset
    
    #2.1.1 calculate dividend
    divd_reg_IF <- divd[match(data_IF$Trddt,divd[,1]),2]
    divd_reg_IC <- divd[match(data_IC$Trddt,divd[,1]),2]
    divd_reg_IH <- divd[match(data_IH$Trddt,divd[,1]),2]
    
    
    #2.2 moving average return - zz500 [completed]
    #2.2.1 1m
    
    sink("output_zz500_30.txt")
    rmu_zz500_30 <- filter(index_zz500$index_zz500, rep(1,30)/30,method="convolution",sides =1)
    rmu_zz500_30 <- na.omit(rmu_zz500_30)
    rmu_zz500_30 <- data.frame(index_zz500$date[30:1884],rmu_zz500_30)
    r_muzz500_reg_30 <- rmu_zz500_30[match(data_IC$Trddt,rmu_zz500_30[,1]),2]
    
    
    p_bs <- data_IC$Spot*exp((data_IC$rf_1m-divd_reg_IC)*data_IC$dt)   #p_bs = s0*e^(rf-q)dt
    p_mu <- data_IC$Spot * exp((r_muzz500_reg_30-divd_reg_IC)*data_IC$dt) #p_mu = s0*e^(rmu-q)dt
    # p_mu <- data_IC$Spot * exp((mean(r_muzz500_reg_30)-divd_reg_IC)*data_IC$dt) #p_mu = s0*e^(rmu-q)dt
    p_spot <- data_IC$Clsprc
    
    df_zz500_30 <- p_spot - p_bs
    df_muzz500_30 <- p_mu - p_bs
    summary(df_zz500_30)
    summary(df_muzz500_30)
    # 
    # df_mu_perc <- df_zz500_30 / p_spot
    # df_bs_perc <- df_muzz500_30 / p_spot
    
    
    hist(df_zz500_30)
    hist(df_muzz500_30)
    
    # x <- -as.numeric(df_muzz500_30 >= 0)
    
    cor.test(p_spot,p_bs,method='pearson')
    cor.test(p_spot,p_mu,method='pearson')
    cor.test(df_zz500_30,df_muzz500_30,method = 'pearson')
    
    
    reg_IC_30 <- lm(df_zz500_30~data_IC$X+df_muzz500_30)  #regression
    summary(reg_IC_30)
    reg_IC_30_noiv <- lm(df_zz500_30~df_muzz500_30)
    summary(reg_IC_30_noiv)
    
    tidy(reg_IC_30)
    tidy(reg_IC_30_noiv)
    
    # reg_IC_perc <- lm(df_bs_perc~data_IC$X+df_mu_perc)  #regression
    # summary(reg_IC_perc)
    
    sink()
    
    #2.2.2 2m
    
    sink("output_zz500_60.txt")
    rmu_zz500_60 <- filter(index_zz500$index_zz500, rep(1,60)/60,method="convolution",sides =1)
    rmu_zz500_60 <- na.omit(rmu_zz500_60)
    rmu_zz500_60 <- data.frame(index_zz500$date[60:1884],rmu_zz500_60)
    r_muzz500_reg_60 <- rmu_zz500_60[match(data_IC$Trddt,rmu_zz500_60[,1]),2]
    
    
    p_bs <- data_IC$Spot*exp((data_IC$rf_2m-divd_reg_IC)*data_IC$dt)   #p_bs = s0*e^(rf-q)dt
    p_mu <- data_IC$Spot * exp((r_muzz500_reg_60-divd_reg_IC)*data_IC$dt) #p_mu = s0*e^(rmu-q)dt
    p_spot <- data_IC$Clsprc
    
    df_zz500_60 <- p_spot - p_bs
    df_muzz500_60 <- p_mu - p_bs
    summary(df_zz500_60)
    summary(df_muzz500_60)
    
    cor.test(p_spot,p_bs,method='pearson')
    cor.test(p_spot,p_mu,method='pearson')
    cor.test(df_zz500_60,df_muzz500_60,method = 'pearson')
    
    reg_IC_60 <- lm(df_zz500_60~data_IC$X+df_muzz500_60)  #regression
    summary(reg_IC_60)
    reg_IC_60_noiv <- lm(df_zz500_60~df_muzz500_60)
    summary(reg_IC_60_noiv)
    
    tidy(reg_IC_60)
    tidy(reg_IC_60_noiv)
    
    result_60 <- data.frame(tidy(reg_IC_60))
    result_60_noiv <- data.frame(tidy(reg_IC_60_noiv))
    
    sink()
    
    #2.2.3 3m
    
    sink("output_zz500_90.txt")
    rmu_zz500_90 <- filter(index_zz500$index_zz500, rep(1,90)/90,method="convolution",sides =1)
    rmu_zz500_90 <- na.omit(rmu_zz500_90)
    rmu_zz500_90 <- data.frame(index_zz500$date[90:1884],rmu_zz500_90)
    r_muzz500_reg_90 <- rmu_zz500_90[match(data_IC$Trddt,rmu_zz500_90[,1]),2]
    
    
    p_bs <- data_IC$Spot*exp((data_IC$rf_3m-divd_reg_IC)*data_IC$dt)   #p_bs = s0*e^(rf-q)dt
    p_mu <- data_IC$Spot * exp((r_muzz500_reg_90-divd_reg_IC)*data_IC$dt) #p_mu = s0*e^(rmu-q)dt
    p_spot <- data_IC$Clsprc
    
    df_zz500_90 <- p_spot - p_bs
    df_muzz500_90 <- p_mu - p_bs
    summary(df_zz500_90)
    summary(df_muzz500_90)
    
    cor.test(p_spot,p_bs,method='pearson')
    cor.test(p_spot,p_mu,method='pearson')
    cor.test(df_zz500_90,df_muzz500_90,method = 'pearson')
    
    reg_IC_90 <- lm(df_zz500_90~data_IC$X+df_muzz500_90)  #regression
    summary(reg_IC_90)
    reg_IC_90_noiv <- lm(df_zz500_90~df_muzz500_90)
    summary(reg_IC_90_noiv)
    
    tidy(reg_IC_90)
    tidy(reg_IC_90_noiv)
    
    result_90 <- data.frame(tidy(reg_IC_90))
    result_90_noiv <- data.frame(tidy(reg_IC_90_noiv))
    
    sink()
  }
  
  # sz50 part
  {
    #2.3 moving average return - sz50
    #2.3.1 1m
    
    sink("output_sz50_30.txt")
    rmu_sz50_30 <- filter(index_sz50$index_sz50, rep(1,30)/30,method="convolution",sides =1)
    rmu_sz50_30 <- na.omit(rmu_sz50_30)
    rmu_sz50_30 <- data.frame(index_sz50$Date[30:1884],rmu_sz50_30)
    r_musz50_reg_30 <- rmu_sz50_30[match(data_IH$Trddt,rmu_sz50_30[,1]),2]
    
    
    p_bs <- data_IH$Spot*exp((data_IH$rf_1m-divd_reg_IH)*data_IH$dt)   #p_bs = s0*e^(rf-q)dt
    p_mu <- data_IH$Spot * exp((r_musz50_reg_30-divd_reg_IH)*data_IH$dt) #p_mu = s0*e^(rmu-q)dt
    p_spot <- data_IH$Clsprc
    
    df_sz50_30 <- p_spot - p_bs
    df_musz50_30 <- p_mu - p_bs
    summary(df_sz50_30)
    summary(df_musz50_30)
    
    cor.test(p_spot,p_bs,method='pearson')
    cor.test(p_spot,p_mu,method='pearson')
    cor.test(df_sz50_30,df_musz50_30,method = 'pearson')
    
    reg_IH_30 <- lm(df_sz50_30~data_IH$X+df_musz50_30)  #regression
    summary(reg_IH_30)
    reg_IH_30_noiv <- lm(df_sz50_30~df_musz50_30)
    summary(reg_IH_30_noiv)
    
    tidy(reg_IH_30)
    tidy(reg_IH_30_noiv)
    
    result_30 <- data.frame(tidy(reg_IH_30))
    result_30_noiv <- data.frame(tidy(reg_IH_30_noiv))
    
    sink()
    
    
    
    #2.3.2 2m
    
    sink("output_sz50_60.txt")
    rmu_sz50_60 <- filter(index_sz50$index_sz50, rep(1,60)/60,method="convolution",sides =1)
    rmu_sz50_60 <- na.omit(rmu_sz50_60)
    rmu_sz50_60 <- data.frame(index_sz50$Date[60:1884],rmu_sz50_60)
    r_musz50_reg_60 <- rmu_sz50_60[match(data_IH$Trddt,rmu_sz50_60[,1]),2]
    
    
    p_bs <- data_IH$Spot*exp((data_IH$rf_2m-divd_reg_IH)*data_IH$dt)   #p_bs = s0*e^(rf-q)dt
    p_mu <- data_IH$Spot * exp((r_musz50_reg_60-divd_reg_IH)*data_IH$dt) #p_mu = s0*e^(rmu-q)dt
    p_spot <- data_IH$Clsprc
    
    df_sz50_60 <- p_spot - p_bs
    df_musz50_60 <- p_mu - p_bs
    summary(df_sz50_60)
    summary(df_musz50_60)
    
    cor.test(p_spot,p_bs,method='pearson')
    cor.test(p_spot,p_mu,method='pearson')
    cor.test(df_sz50_60,df_musz50_60,method = 'pearson')
    
    reg_IH_60 <- lm(df_sz50_60~data_IH$X+df_musz50_60)  #regression
    summary(reg_IH_60)
    reg_IH_60_noiv <- lm(df_sz50_60~df_musz50_60)
    summary(reg_IH_60_noiv)
    
    tidy(reg_IH_60)
    tidy(reg_IH_60_noiv)
    
    result_60 <- data.frame(tidy(reg_IH_60))
    result_60_noiv <- data.frame(tidy(reg_IH_60_noiv))
    
    sink()
    
    
    #2.3.3 3m
    
    sink("output_sz50_90.txt")
    rmu_sz50_90 <- filter(index_sz50$index_sz50, rep(1,90)/90,method="convolution",sides =1)
    rmu_sz50_90 <- na.omit(rmu_sz50_90)
    rmu_sz50_90 <- data.frame(index_sz50$Date[90:1884],rmu_sz50_90)
    r_musz50_reg_90 <- rmu_sz50_90[match(data_IH$Trddt,rmu_sz50_90[,1]),2]
    
    
    p_bs <- data_IH$Spot*exp((data_IH$rf_3m-divd_reg_IH)*data_IH$dt)   #p_bs = s0*e^(rf-q)dt
    p_mu <- data_IH$Spot * exp((r_musz50_reg_90-divd_reg_IH)*data_IH$dt) #p_mu = s0*e^(rmu-q)dt
    p_spot <- data_IH$Clsprc
    
    df_sz50_90 <- p_spot - p_bs
    df_musz50_90 <- p_mu - p_bs
    summary(df_sz50_90)
    summary(df_musz50_90)
    
    cor.test(p_spot,p_bs,method='pearson')
    cor.test(p_spot,p_mu,method='pearson')
    cor.test(df_sz50_90,df_musz50_90,method = 'pearson')
    
    reg_IH_90 <- lm(df_sz50_90~data_IH$X+df_musz50_90)  #regression
    summary(reg_IH_90)
    reg_IH_90_noiv <- lm(df_sz50_90~df_musz50_90)
    summary(reg_IH_90_noiv)
    
    tidy(reg_IH_90)
    tidy(reg_IH_90_noiv)
    
    result_90 <- data.frame(tidy(reg_IH_90))
    result_90_noiv <- data.frame(tidy(reg_IH_90_noiv))
    
    sink()
    
    
  }
  
  # hs300 part
  {
    #2.4 moving average return - hs300
    #2.4.1 1m
    
    sink("output_hs300_30.txt")
    rmu_hs300_30 <- filter(index_hs300$rmu_hs300, rep(1,30)/30,method="convolution",sides =1)
    rmu_hs300_30 <- na.omit(rmu_hs300_30)
    rmu_hs300_30 <- data.frame(index_hs300$Date[30:2010],rmu_hs300_30)
    r_muhs300_reg_30 <- rmu_hs300_30[match(data_IF$Trddt,rmu_hs300_30[,1]),2]
    
    
    p_bs <- data_IF$Spot*exp((data_IF$rf_1m-divd_reg_IF)*data_IF$dt)   #p_bs = s0*e^(rf-q)dt
    p_mu <- data_IF$Spot * exp((r_muhs300_reg_30-divd_reg_IF)*data_IF$dt) #p_mu = s0*e^(rmu-q)dt
    p_spot <- data_IF$Clsprc
    
    df_hs300_30 <- p_spot - p_bs
    df_muhs300_30 <- p_mu - p_bs
    summary(df_hs300_30)
    summary(df_muhs300_30)
    
    cor.test(p_spot,p_bs,method='pearson')
    cor.test(p_spot,p_mu,method='pearson')
    cor.test(df_hs300_30,df_muhs300_30,method = 'pearson')
    
    reg_IF_30 <- lm(df_hs300_30~data_IF$X+df_muhs300_30)  #regression
    summary(reg_IF_30)
    reg_IF_30_noiv <- lm(df_hs300_30~df_muhs300_30)
    summary(reg_IF_30_noiv)
    
    tidy(reg_IF_30)
    tidy(reg_IF_30_noiv)
    
    result_30 <- data.frame(tidy(reg_IF_30))
    result_30_noiv <- data.frame(tidy(reg_IF_30_noiv))
    
    sink()
    
    
    
    
    #2.4.2 2m
    
    sink("output_hs300_60.txt")
    rmu_hs300_60 <- filter(index_hs300$rmu_hs300, rep(1,60)/60,method="convolution",sides =1)
    rmu_hs300_60 <- na.omit(rmu_hs300_60)
    rmu_hs300_60 <- data.frame(index_hs300$Date[60:2010],rmu_hs300_60)
    r_muhs300_reg_60 <- rmu_hs300_60[match(data_IF$Trddt,rmu_hs300_60[,1]),2]
    
    
    p_bs <- data_IF$Spot*exp((data_IF$rf_2m-divd_reg_IF)*data_IF$dt)   #p_bs = s0*e^(rf-q)dt
    p_mu <- data_IF$Spot * exp((r_muhs300_reg_60-divd_reg_IF)*data_IF$dt) #p_mu = s0*e^(rmu-q)dt
    p_spot <- data_IF$Clsprc
    
    df_hs300_60 <- p_spot - p_bs
    df_muhs300_60 <- p_mu - p_bs
    summary(df_hs300_60)
    summary(df_muhs300_60)
    
    cor.test(p_spot,p_bs,method='pearson')
    cor.test(p_spot,p_mu,method='pearson')
    cor.test(df_hs300_60,df_muhs300_60,method = 'pearson')
    
    reg_IF_60 <- lm(df_hs300_60~data_IF$X+df_muhs300_60)  #regression
    summary(reg_IF_60)
    reg_IF_60_noiv <- lm(df_hs300_60~df_muhs300_60)
    summary(reg_IF_60_noiv)
    
    tidy(reg_IF_60)
    tidy(reg_IF_60_noiv)
    
    result_60 <- data.frame(tidy(reg_IF_60))
    result_60_noiv <- data.frame(tidy(reg_IF_60_noiv))
    
    sink()
    
    #2.4.3 3m
    
    sink("output_hs300_90.txt")
    rmu_hs300_90 <- filter(index_hs300$rmu_hs300, rep(1,90)/90,method="convolution",sides =1)
    rmu_hs300_90 <- na.omit(rmu_hs300_90)
    rmu_hs300_90 <- data.frame(index_hs300$Date[90:2010],rmu_hs300_90)
    r_muhs300_reg_90 <- rmu_hs300_90[match(data_IF$Trddt,rmu_hs300_90[,1]),2]
    
    
    p_bs <- data_IF$Spot*exp((data_IF$rf_3m-divd_reg_IF)*data_IF$dt)   #p_bs = s0*e^(rf-q)dt
    p_mu <- data_IF$Spot * exp((r_muhs300_reg_90-divd_reg_IF)*data_IF$dt) #p_mu = s0*e^(rmu-q)dt
    p_spot <- data_IF$Clsprc
    
    df_hs300_90 <- p_spot - p_bs
    df_muhs300_90 <- p_mu - p_bs
    summary(df_hs300_90)
    summary(df_muhs300_90)
    
    cor.test(p_spot,p_bs,method='pearson')
    cor.test(p_spot,p_mu,method='pearson')
    cor.test(df_hs300_90,df_muhs300_90,method = 'pearson')
    
    reg_IF_90 <- lm(df_hs300_90~data_IF$X+df_muhs300_90)  #regression
    summary(reg_IF_90)
    reg_IF_90_noiv <- lm(df_hs300_90~df_muhs300_90)
    summary(reg_IF_90_noiv)
    
    tidy(reg_IF_90)
    tidy(reg_IF_90_noiv)
    
    result_90 <- data.frame(tidy(reg_IF_90))
    result_90_noiv <- data.frame(tidy(reg_IF_90_noiv))
    
    sink()
    
    
  }
  
  # all contracts combined
  {
    #All contracts combined
    
    data_IC$r_mu_all <- NULL
    data_IF$r_mu_all <- NULL
    data_IH$r_mu_all <- NULL
    
    
    #2.4 moving average return - hs300
    #2.4.1 1m
    
    data_IC <- data.frame(data_IC, r_muzz500_reg_30)
    data_IF <- data.frame(data_IF, r_muhs300_reg_30)
    data_IH <- data.frame(data_IH, r_musz50_reg_30)
    
    # change colname 
    names(data_IC)[17] <- "r_mu_all"
    names(data_IF)[17] <- "r_mu_all"
    names(data_IH)[17] <- "r_mu_all"
    
    
    data_all <- data.frame(rbind(data_IC,data_IF,data_IH))
    divd_reg <- c(divd_reg_IC, divd_reg_IF, divd_reg_IH)
    
    sink("output_all_30.txt")
    
    
    p_bs <- data_all$Spot*exp((data_all$rf_1m-divd_reg)*data_all$dt)   #p_bs = s0*e^(rf-q)dt
    p_mu <- data_all$Spot * exp((data_all$r_mu_all-divd_reg)*data_all$dt) #p_mu = s0*e^(rmu-q)dt
    p_spot <- data_all$Clsprc
    
    df_all_30 <- p_spot - p_bs
    df_muall_30 <- p_mu - p_bs
    summary(df_all_30)
    summary(df_muall_30)
    
    cor.test(p_spot,p_bs,method='pearson')
    cor.test(p_spot,p_mu,method='pearson')
    cor.test(df_all_30,df_muall_30,method = 'pearson')
    
    reg_all_30 <- lm(df_all_30~data_all$X+df_muall_30)  #regression
    summary(reg_all_30)
    reg_all_30_noiv <- lm(df_all_30~df_muall_30)
    summary(reg_all_30_noiv)
    
    tidy(reg_all_30)
    tidy(reg_all_30_noiv)
    
    result_30 <- data.frame(tidy(reg_all_30))
    result_30_noiv <- data.frame(tidy(reg_all_30_noiv))
    
    sink()
    
    
    
    
    #2.4.2 2m
    
    sink("output_all_60.txt")
    
    data_IC$r_mu_all <- NULL
    data_IF$r_mu_all <- NULL
    data_IH$r_mu_all <- NULL
    
    data_IC <- data.frame(data_IC, r_muzz500_reg_60)
    data_IF <- data.frame(data_IF, r_muhs300_reg_60)
    data_IH <- data.frame(data_IH, r_musz50_reg_60)
    
    names(data_IC)[17] <- "r_mu_all"
    names(data_IF)[17] <- "r_mu_all"
    names(data_IH)[17] <- "r_mu_all"
    
    data_all <- data.frame(rbind(data_IC,data_IF,data_IH))
    
    p_bs <- data_all$Spot*exp((data_all$rf_2m-divd_reg)*data_all$dt)   #p_bs = s0*e^(rf-q)dt
    p_mu <- data_all$Spot * exp((data_all$r_mu_all-divd_reg)*data_all$dt) #p_mu = s0*e^(rmu-q)dt
    p_spot <- data_all$Clsprc
    
    df_all_60 <- p_spot - p_bs
    df_muall_60 <- p_mu - p_bs
    summary(df_all_60)
    summary(df_muall_60)
    
    cor.test(p_spot,p_bs,method='pearson')
    cor.test(p_spot,p_mu,method='pearson')
    cor.test(df_all_60,df_muall_60,method = 'pearson')
    
    reg_all_60 <- lm(df_all_60~data_all$X+df_muall_60)  #regression
    summary(reg_all_60)
    reg_all_60_noiv <- lm(df_all_60~df_muall_60)
    summary(reg_all_60_noiv)
    
    tidy(reg_all_60)
    tidy(reg_all_60_noiv)
    
    result_60 <- data.frame(tidy(reg_all_60))
    result_60_noiv <- data.frame(tidy(reg_all_60_noiv))
    
    sink()
    
    
    #2.4.3 3m
    
    sink("output_all_90.txt")
    
    data_IC$r_mu_all <- NULL
    data_IF$r_mu_all <- NULL
    data_IH$r_mu_all <- NULL
    
    data_IC <- data.frame(data_IC, r_muzz500_reg_60)
    data_IF <- data.frame(data_IF, r_muhs300_reg_60)
    data_IH <- data.frame(data_IH, r_musz50_reg_60)
    
    names(data_IC)[17] <- "r_mu_all"
    names(data_IF)[17] <- "r_mu_all"
    names(data_IH)[17] <- "r_mu_all"
    
    data_all <- data.frame(rbind(data_IC,data_IF,data_IH))
    
    p_bs <- data_all$Spot*exp((data_all$rf_3m-divd_reg)*data_all$dt)   #p_bs = s0*e^(rf-q)dt
    p_mu <- data_all$Spot * exp((data_all$r_mu_all-divd_reg)*data_all$dt) #p_mu = s0*e^(rmu-q)dt
    p_spot <- data_all$Clsprc
    
    df_all_90 <- p_spot - p_bs
    df_muall_90 <- p_mu - p_bs
    summary(df_all_90)
    summary(df_muall_90)
    
    cor.test(p_spot,p_bs,method='pearson')
    cor.test(p_spot,p_mu,method='pearson')
    cor.test(df_all_90,df_muall_90,method = 'pearson')
    
    reg_all_90 <- lm(df_all_90~data_all$X+df_muall_90)  #regression
    summary(reg_all_90)
    reg_all_90_noiv <- lm(df_all_90~df_muall_90)
    summary(reg_all_90_noiv)
    
    tidy(reg_all_90)
    tidy(reg_all_90_noiv)
    
    result_90 <- data.frame(tidy(reg_all_90))
    result_90_noiv <- data.frame(tidy(reg_all_90_noiv))
    
    sink()
    
  }
  
  
}


# Using risk premium method
{
  # zz500 part
  {
    sink("output_riskp_zz500.txt")
    
    r_muzz500_raw <- index_zz500[match(data_IC$Trddt,index_zz500[,1]),2]
    # riskp_IC <- data_IC$r_mu_all - data_IC$rf_1m
    riskp_IC <- r_muzz500_raw - data_IC$rf_1m
    riskp_IC_mean <- mean(riskp_IC)
    divd_IC_mean <- mean(divd_reg_IC)
    riskp_IC_reg <- riskp_IC_mean + divd_IC_mean
    
    p_bs <- data_IC$Spot*exp((data_IC$rf_1m-divd_IC_mean)*data_IC$dt) 
    p_mu <- data_IC$Spot * exp((riskp_IC_reg+data_IC$rf_1m-divd_reg_IC)*data_IC$dt)
    p_spot <- data_IC$Clsprc
    
    df_zz500_riskp <- p_spot - p_bs
    df_muzz500_riskp <- p_mu - p_bs
    
    summary(df_zz500_riskp)
    summary(df_muzz500_riskp)
    
    cor.test(p_spot,p_bs,method='pearson')
    cor.test(p_spot,p_mu,method='pearson')
    cor.test(df_zz500_riskp,df_muzz500_riskp,method = 'pearson')
    
    reg_IC_riskp <- lm(df_zz500_riskp~data_IC$X+df_muzz500_riskp)  #regression
    summary(reg_IC_riskp)
    
    reg_IC_riskp_noiv <- lm(df_zz500_riskp~df_muzz500_riskp)
    summary(reg_IC_riskp_noiv)
    
    tidy(reg_IC_riskp)
    tidy(reg_IC_riskp_noiv)
    
    sink()
    
  }
  
  # sz50 part
  {
    sink("output_riskp_sz50.txt")
    
    r_musz50_raw <- index_sz50[match(data_IH$Trddt,index_sz50[,1]),2]
    riskp_IH <- r_musz50_raw - data_IH$rf_1m
    riskp_IH_mean <- mean(riskp_IH)
    divd_IH_mean <- mean(divd_reg_IH)
    riskp_IH_reg <- riskp_IH_mean + divd_IH_mean
    
    p_bs <- data_IH$Spot*exp((data_IH$rf_1m-divd_IH_mean)*data_IH$dt) 
    p_mu <- data_IH$Spot * exp((riskp_IH_reg+data_IH$rf_1m-divd_reg_IH)*data_IH$dt)
    p_spot <- data_IH$Clsprc
    
    df_sz50_riskp <- p_spot - p_bs
    df_musz50_riskp <- p_mu - p_bs
    
    summary(df_sz50_riskp)
    summary(df_musz50_riskp)
    
    cor.test(p_spot,p_bs,method='pearson')
    cor.test(p_spot,p_mu,method='pearson')
    cor.test(df_sz50_riskp,df_musz50_riskp,method = 'pearson')
    
    reg_IH_riskp <- lm(df_sz50_riskp~data_IH$X+df_musz50_riskp)  #regression
    summary(reg_IH_riskp)
    
    reg_IH_riskp_noiv <- lm(df_sz50_riskp~df_musz50_riskp)
    summary(reg_IH_riskp_noiv)
    
    tidy(reg_IH_riskp)
    tidy(reg_IH_riskp_noiv)
    
    sink()
    
    
  }
  
  # hs300 part
  {
    
    sink("output_riskp_hs300.txt")
    
    r_muhs300_raw <- index_hs300[match(data_IF$Trddt,index_hs300[,1]),2]
    riskp_IF <- r_muhs300_raw - data_IF$rf_1m
    riskp_IF_mean <- mean(riskp_IF)
    divd_IF_mean <- mean(divd_reg_IF)
    riskp_IF_reg <- riskp_IF_mean + divd_IF_mean
    
    p_bs <- data_IF$Spot*exp((data_IF$rf_1m-divd_IF_mean)*data_IF$dt) 
    p_mu <- data_IF$Spot * exp((riskp_IF_reg+data_IF$rf_1m-divd_reg_IF)*data_IF$dt)
    p_spot <- data_IF$Clsprc
    
    df_hs300_riskp <- p_spot - p_bs
    df_muhs300_riskp <- p_mu - p_bs
    
    summary(df_hs300_riskp)
    summary(df_muhs300_riskp)
    
    cor.test(p_spot,p_bs,method='pearson')
    cor.test(p_spot,p_mu,method='pearson')
    cor.test(df_hs300_riskp,df_muhs300_riskp,method = 'pearson')
    
    reg_IF_riskp <- lm(df_hs300_riskp~data_IF$X+df_muhs300_riskp)  #regression
    summary(reg_IF_riskp)
    
    reg_IF_riskp_noiv <- lm(df_hs300_riskp~df_muhs300_riskp)
    summary(reg_IF_riskp_noiv)
    
    tidy(reg_IF_riskp)
    tidy(reg_IF_riskp_noiv)
    
    sink()
    
    
  }
  
  # All contracts combined
  {
    
    sink("output_riskp_all.txt")
    
    riskp_all <- r_mu_raw- data_all$rf_1m
    riskp_all_mean <- mean(riskp_all)
    divd_all_mean <- mean(divd_reg)
    riskp_all_reg <- riskp_all_mean + divd_all_mean
    
    p_bs <- data_all$Spot*exp((data_all$rf_1m-divd_all_mean)*data_all$dt) 
    p_mu <- data_all$Spot * exp((riskp_all_reg+data_all$rf_1m-divd_reg)*data_all$dt)
    # p_mu <- data_all$Spot * exp(riskp_all_reg*data_all$dt)
    p_spot <- data_all$Clsprc
    
    df_all_riskp <- p_spot - p_bs
    df_muall_riskp <- p_mu - p_bs
    
    summary(df_all_riskp)
    summary(df_muall_riskp)
    
    cor.test(p_spot,p_bs,method='pearson')
    cor.test(p_spot,p_mu,method='pearson')
    cor.test(df_all_riskp,df_muall_riskp,method = 'pearson')
    
    reg_all_riskp <- lm(df_all_riskp~data_all$X+df_muall_riskp)  #regression
    summary(reg_all_riskp)
    
    reg_all_riskp_noiv <- lm(df_all_riskp~df_muall_riskp)
    summary(reg_all_riskp_noiv)
    
    tidy(reg_all_riskp)
    tidy(reg_all_riskp_noiv)
    
    sink()
    
    
    
  }

}

# Using both methods

alpha <- 0.7

#All contracts combined


data_IC$r_mu_all <- NULL
data_IF$r_mu_all <- NULL
data_IH$r_mu_all <- NULL


data_IC <- data.frame(data_IC, r_muzz500_reg_90)
data_IF <- data.frame(data_IF, r_muhs300_reg_90)
data_IH <- data.frame(data_IH, r_musz50_reg_90)

names(data_IC)[17] <- "r_mu_all"
names(data_IF)[17] <- "r_mu_all"
names(data_IH)[17] <- "r_mu_all"

data_all <- data.frame(rbind(data_IC,data_IF,data_IH))

{
  
  sink("output_all_com.txt")
  
  # for (i in seq(0.1,1,0.1)){
  
  # print(i)
  # 
  # alpha <- i
  
  alpha <- 0.2
  rmu_com_reg <- alpha * mean(data_all$r_mu_all) + (1-alpha) * (riskp_all_mean+mean(data_all$rf_3m))
  
  p_bs <- data_all$Spot*exp((data_all$rf_3m-divd_reg)*data_all$dt)   #p_bs = s0*e^(rf-q)dt
  p_mu <- data_all$Spot * exp((rmu_com_reg-divd_reg)*data_all$dt) #p_mu = s0*e^(rmu-q)dt
  p_spot <- data_all$Clsprc
  
  df_all_com <- p_spot - p_bs
  df_muall_com <- p_mu - p_bs
  
  reg_all_com <- lm(df_all_com~data_all$X+df_muall_com)  #regression
  summary(reg_all_com)
  reg_all_com_noiv <- lm(df_all_com~df_muall_com)
  summary(reg_all_com_noiv)
  
  # r_square <- c()
  # r_square_noiv <- c()
  r_square[2] <- summary(reg_all_com)$adj.r.squared
  r_square_noiv[2] <- summary(reg_all_com_noiv)$adj.r.squared
  
  tidy(reg_all_com)
  tidy(reg_all_com_noiv)
  
  
  # }
  
  
  sink()
}



#4. analysis of basic results

  #4.1 plots for regression analysis

# fitted vs true values
# IC
data_IC$predregIC30 <- predict(reg_IC_30,newdata = data_IC)
plot1 <- ggplot(data = data_IC, mapping = aes(x= predregIC30, y= df_zz500_30))+
  geom_point(alpha = 0.2, color = "black")+
  geom_smooth(aes(x = predregIC30,y = df_zz500_30),color = "black")+
  labs(x='fitted',y='real values',title="fitted vs real values (1m)")

data_IC$predregIC60 <- predict(reg_IC_60,newdata = data_IC)
plot2 <- ggplot(data = data_IC, mapping = aes(x= predregIC60, y= df_zz500_60))+
  geom_point(alpha = 0.2, color = "black")+
  geom_smooth(aes(x = predregIC60,y = df_zz500_60),color = "black")+
  labs(x='fitted',y='real values',title="fitted vs real values (2m)")

data_IC$predregIC90 <- predict(reg_IC_90,newdata = data_IC)
plot3 <- ggplot(data = data_IC, mapping = aes(x= predregIC90, y= df_zz500_90))+
  geom_point(alpha = 0.2, color = "black")+
  geom_smooth(aes(x = predregIC90,y = df_zz500_90),color = "black")+
  labs(x='fitted',y='real values',title="fitted vs real values (3m)")

# IF
data_IF$predregIF30 <- predict(reg_IF_30,newdata = data_IF)
plot4 <- ggplot(data = data_IF, mapping = aes(x= predregIF30, y= df_hs300_30))+
  geom_point(alpha = 0.2, color = "black")+
  geom_smooth(aes(x = predregIF30,y = df_hs300_30),color = "black")+
  labs(x='fitted',y='real values',title="fitted vs real values (1m)")

data_IF$predregIF60 <- predict(reg_IF_60,newdata = data_IF)
plot5 <- ggplot(data = data_IF, mapping = aes(x= predregIF60, y= df_hs300_60))+
  geom_point(alpha = 0.2, color = "black")+
  geom_smooth(aes(x = predregIF60,y = df_hs300_60),color = "black")+
  labs(x='fitted',y='real values',title="fitted vs real values (2m)")

data_IF$predregIF90 <- predict(reg_IF_90,newdata = data_IF)
plot6 <- ggplot(data = data_IF, mapping = aes(x= predregIF90, y= df_hs300_90))+
  geom_point(alpha = 0.2, color = "black")+
  geom_smooth(aes(x = predregIF90,y = df_hs300_90),color = "black")+
  labs(x='fitted',y='real values',title="fitted vs real values (3m)")


# IH
data_IH$predregIH30 <- predict(reg_IH_30,newdata = data_IH)
plot7 <- ggplot(data = data_IH, mapping = aes(x= predregIH30, y= df_sz50_30))+
  geom_point(alpha = 0.2, color = "black")+
  geom_smooth(aes(x = predregIH30,y = df_sz50_30),color = "black")+
  labs(x='fitted',y='real values',title="fitted vs real values (1m)")

data_IH$predregIH60 <- predict(reg_IH_60,newdata = data_IH)
plot8 <- ggplot(data = data_IH, mapping = aes(x= predregIH60, y= df_sz50_60))+
  geom_point(alpha = 0.2, color = "black")+
  geom_smooth(aes(x = predregIH60,y = df_sz50_60),color = "black")+
  labs(x='fitted',y='real values',title="fitted vs real values (2m)")

data_IH$predregIH90 <- predict(reg_IH_90,newdata = data_IH)
plot9 <- ggplot(data = data_IH, mapping = aes(x= predregIH90, y= df_sz50_90))+
  geom_point(alpha = 0.2, color = "black")+
  geom_smooth(aes(x = predregIH90,y = df_sz50_90),color = "black")+
  labs(x='fitted',y='real values',title="fitted vs real values (3m)")

# residuals
# IC
plot10 <- ggplot(data = data_IC, mapping = aes(x= predregIC30, y= predregIC30 - df_zz500_30))+
  geom_point(alpha = 0.2, color = "black")+
  geom_smooth(aes(x = predregIC30,y =  predregIC30 - df_zz500_30),color = "black")+
  labs(x='fitted',y='real values',title="fitted vs residuals (1m)")

plot11 <- ggplot(data = data_IC, mapping = aes(x= predregIC60, y= predregIC60 - df_zz500_60))+
  geom_point(alpha = 0.2, color = "black")+
  geom_smooth(aes(x = predregIC60,y =  predregIC60 - df_zz500_60),color = "black")+
  labs(x='fitted',y='real values',title="fitted vs residuals (2m)")

plot12 <- ggplot(data = data_IC, mapping = aes(x= predregIC90, y= predregIC90 - df_zz500_90))+
  geom_point(alpha = 0.2, color = "black")+
  geom_smooth(aes(x = predregIC90,y = predregIC90 - df_zz500_90),color = "black")+
  labs(x='fitted',y='real values',title="fitted vs residuals (3m)")
  
# IF
plot13 <- ggplot(data = data_IF, mapping = aes(x= predregIF30, y= predregIF30 - df_hs300_30))+
  geom_point(alpha = 0.2, color = "black")+
  geom_smooth(aes(x = predregIF30,y =  predregIF30 - df_hs300_30),color = "black")+
  labs(x='fitted',y='real values',title="fitted vs residuals (1m)")

plot14 <- ggplot(data = data_IF, mapping = aes(x= predregIF60, y= predregIF60 - df_hs300_60))+
  geom_point(alpha = 0.2, color = "black")+
  geom_smooth(aes(x = predregIF60,y =  predregIF60 - df_hs300_60),color = "black")+
  labs(x='fitted',y='real values',title="fitted vs residuals (2m)")

plot15 <- ggplot(data = data_IF, mapping = aes(x= predregIF90, y= predregIF90 - df_hs300_90))+
  geom_point(alpha = 0.2, color = "black")+
  geom_smooth(aes(x = predregIF90,y = predregIF90 - df_hs300_90),color = "black")+
  labs(x='fitted',y='real values',title="fitted vs residuals (3m)")


# IH
plot16 <- ggplot(data = data_IH, mapping = aes(x= predregIH30, y= predregIH30 - df_sz50_30))+
  geom_point(alpha = 0.2, color = "black")+
  geom_smooth(aes(x = predregIH30,y =  predregIH30 - df_sz50_30),color = "black")+
  labs(x='fitted',y='real values',title="fitted vs residuals (1m)")

plot17 <- ggplot(data = data_IH, mapping = aes(x= predregIH60, y= predregIH60 - df_sz50_60))+
  geom_point(alpha = 0.2, color = "black")+
  geom_smooth(aes(x = predregIH60,y =  predregIH60 - df_sz50_60),color = "black")+
  labs(x='fitted',y='real values',title="fitted vs residuals (2m)")

plot18 <- ggplot(data = data_IH, mapping = aes(x= predregIH90, y= predregIH90 - df_sz50_90))+
  geom_point(alpha = 0.2, color = "black")+
  geom_smooth(aes(x = predregIH90,y = predregIH90 - df_sz50_90),color = "black")+
  labs(x='fitted',y='real values',title="fitted vs residuals (3m)")

# organize plots


subvp_1 <- viewport(width=1/2,height=1/3,x=1/4,y=5/6)
subvp_2 <- viewport(width=1/2,height=1/3,x=1/4,y=3/6)
subvp_3 <- viewport(width=1/2,height=1/3,x=1/4,y=1/6)
subvp_4 <- viewport(width=1/2,height=1/3,x=3/4,y=5/6)
subvp_5 <- viewport(width=1/2,height=1/3,x=3/4,y=3/6)
subvp_6 <- viewport(width=1/2,height=1/3,x=3/4,y=1/6)

dev.new()
print(plot1,vp = subvp_1)
print(plot2,vp = subvp_2)
print(plot3,vp = subvp_3)
print(plot10,vp = subvp_4)
print(plot11,vp = subvp_5)
print(plot12,vp = subvp_6)
# jpeg(file="zz500.jpeg")

dev.new()
print(plot4,vp = subvp_1)
print(plot5,vp = subvp_2)
print(plot6,vp = subvp_3)
print(plot13,vp = subvp_4)
print(plot14,vp = subvp_5)
print(plot15,vp = subvp_6)

dev.new()
print(plot7,vp = subvp_1)
print(plot8,vp = subvp_2)
print(plot9,vp = subvp_3)
print(plot16,vp = subvp_4)
print(plot17,vp = subvp_5)
print(plot18,vp = subvp_6)

    #4.2 rse
rmse <- function(y, f) {
  sqrt(mean((y - f)^2))
}

data_IC$predregIC30_riskp <- predict(reg_IC_riskp,newdata = data_IC)
data_IF$predregIF30_riskp <- predict(reg_IF_riskp,newdata = data_IF)
data_IH$predregIH30_riskp <- predict(reg_IH_riskp,newdata = data_IH)

sink("rmse.csv")
rmse(df_zz500_30,data_IC$predregIC30 )
rmse(df_zz500_60,data_IC$predregIC60 )
rmse(df_zz500_90,data_IC$predregIC90 )
rmse(df_hs300_30,data_IF$predregIF30 )
rmse(df_hs300_60,data_IF$predregIF60 )
rmse(df_hs300_90,data_IF$predregIF90 )
rmse(df_sz50_30,data_IH$predregIH30 )
rmse(df_sz50_60,data_IH$predregIH60 )
rmse(df_sz50_90,data_IH$predregIH90 )
rmse(df_zz500_riskp,data_IC$predregIC30_riskp)
rmse(df_hs300_riskp,data_IF$predregIF30_riskp)
rmse(df_sz50_riskp,data_IH$predregIH30_riskp)

sink()

  #4.2 statistics for c_mu1m, c_mu2m, c_mu3m, c_mu6m


  #5.1 different maturity


#6. regress results for allowing parameters to float


