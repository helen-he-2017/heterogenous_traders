setwd('~/Desktop/data')
library(broom)
library(grid)

# data sourcing
{
  y <- read.csv("~/Desktop/data/data_new.csv",header = TRUE)
  rf <- read.csv("~/Desktop/data/rf.csv",header = TRUE)
  index_sz50 <- read.csv("~/Desktop/data/rmu_sz50.csv",header = TRUE)
  index_zz500 <- read.csv("~/Desktop/data/rmu_zz500.csv",header = TRUE)
  index_hs300 <- read.csv("~/Desktop/data/rmu_hs300.csv",header = TRUE)
  divd <- read.csv("~/Desktop/data/divd.csv",header = TRUE)
  
  data_IC <- y[grep("IC",y$Agmtcd),]  #IC: zz500
  data_IF <- y[grep("IF",y$Agmtcd),]  #IF: hs300
  data_IH <- y[grep("IH",y$Agmtcd),]  #IH: sz50
  
  divd_reg_IF <- divd[match(data_IF$Trddt,divd[,1]),2]
  divd_reg_IC <- divd[match(data_IC$Trddt,divd[,1]),2]
  divd_reg_IH <- divd[match(data_IH$Trddt,divd[,1]),2]
  
  
  coefficient_IC_iv_dfmu <- c()
  coefficient_IC_iv_X <- c()
  coefficient_IC_noiv_dfmu <- c()
  coefficient_IC_iv_intercept <- c()
  coefficient_IC_noiv_intercept <- c()
  R_squared_IC_iv<- c()
  R_squared_IC_noiv <- c()
  
  coefficient_IF_iv_dfmu <- c()
  coefficient_IF_iv_X <- c()
  coefficient_IF_noiv_dfmu <- c()
  coefficient_IF_iv_intercept <- c()
  coefficient_IF_noiv_intercept <- c()
  R_squared_IF_iv <- c()
  R_squared_IF_noiv <- c()
  
  coefficient_IH_iv_dfmu <- c()
  coefficient_IH_iv_X <- c()
  coefficient_IH_noiv_dfmu <- c()
  coefficient_IH_iv_intercept <- c()
  coefficient_IH_noiv_intercept <- c()
  R_squared_IH_iv<- c()
  R_squared_IH_noiv<- c()
  
  coefficient_all_iv_dfmu <- c()
  coefficient_all_iv_X <- c()
  coefficient_all_noiv_dfmu <- c()
  coefficient_all_iv_intercept <- c()
  coefficient_all_noiv_intercept <- c()
  R_squared_all_iv <- c()
  R_squared_all_noiv<- c()
  
}

# IC
for (rolldate in seq(10,360,10))
{
  rmu_zz500_30 <- filter(index_zz500$index_zz500, rep(1,rolldate)/rolldate,method="convolution",sides =1)
  rmu_zz500_30 <- na.omit(rmu_zz500_30)
  rmu_zz500_30 <- data.frame(index_zz500$date[rolldate:2403],rmu_zz500_30)
  r_muzz500_reg_30 <- rmu_zz500_30[match(data_IC$Trddt,rmu_zz500_30[,1]),2]
  
  
  p_bs <- data_IC$Spot*exp((data_IC$rf_1m-divd_reg_IC)*data_IC$dt)   #p_bs = s0*e^(rf-q)dt
  # p_mu <- data_IC$Spot * exp((mean(r_muzz500_reg_30)-divd_reg_IC)*data_IC$dt) #p_mu = s0*e^(rmu-q)dt
  p_mu <- data_IC$Spot * exp((r_muzz500_reg_30-divd_reg_IC)*data_IC$dt) #p_mu = s0*e^(rmu-q)dt
  p_spot <- data_IC$Clsprc
  
  df_zz500_30 <- p_spot - p_bs
  df_muzz500_30 <- p_mu - p_bs
  # summary(df_zz500_30)
  # summary(df_muzz500_30)
  # 
  # df_mu_perc <- df_zz500_30 / p_spot
  # df_bs_perc <- df_muzz500_30 / p_spot
  
  
  # hist(df_zz500_30)
  # hist(df_muzz500_30)
  # 
  # cor.test(p_spot,p_bs,method='pearson')
  # cor.test(p_spot,p_mu,method='pearson')
  # cor.test(df_zz500_30,df_muzz500_30,method = 'pearson')
  
  
  reg_IC_30 <- lm(df_zz500_30~data_IC$X+df_muzz500_30)  #regression
  # summary(reg_IC_30)
  reg_IC_30_noiv <- lm(df_zz500_30~df_muzz500_30)
  # summary(reg_IC_30_noiv)
  
  # tidy(reg_IC_30)
  # tidy(reg_IC_30_noiv)
  
  coefficient_IC_iv_dfmu[rolldate/10] <- as.numeric(reg_IC_30$coefficients[3])
  coefficient_IC_noiv_dfmu[rolldate/10] <- as.numeric(reg_IC_30_noiv$coefficients[2])
  coefficient_IC_iv_X[rolldate/10] <- as.numeric(reg_IC_30$coefficients[2])
  coefficient_IC_iv_intercept[rolldate/10] <- as.numeric(reg_IC_30$coefficients[1])
  coefficient_IC_noiv_intercept[rolldate/10] <- as.numeric(reg_IC_30_noiv$coefficients[1])
  R_squared_IC_iv[rolldate/10] <- as.numeric(summary(reg_IC_30)$adj.r.squared)
  R_squared_IC_noiv[rolldate/10] <- as.numeric(summary(reg_IC_30_noiv)$adj.r.squared)
  
  # reg_IC_perc <- lm(df_bs_perc~data_IC$X+df_mu_perc)  #regression
  # summary(reg_IC_perc)
  
  # sink()
}

# IF
for (rolldate in seq(10,360,10))
{
  rmu_hs300_30 <- filter(index_hs300$index_hs300, rep(1,rolldate)/rolldate,method="convolution",sides =1)
  rmu_hs300_30 <- na.omit(rmu_hs300_30)
  rmu_hs300_30 <- data.frame(index_hs300$date[rolldate:2403],rmu_hs300_30)
  r_muhs300_reg_30 <- rmu_hs300_30[match(data_IF$Trddt,rmu_hs300_30[,1]),2]
  
  
  p_bs <- data_IF$Spot*exp((data_IF$rf_1m-divd_reg_IF)*data_IF$dt)   #p_bs = s0*e^(rf-q)dt
  # p_mu <- data_IF$Spot * exp((mean(r_muhs300_reg_30)-divd_reg_IF)*data_IF$dt) #p_mu = s0*e^(rmu-q)dt
  p_mu <- data_IF$Spot * exp((r_muhs300_reg_30-divd_reg_IF)*data_IF$dt) #p_mu = s0*e^(rmu-q)dt
  p_spot <- data_IF$Clsprc
  
  df_hs300_30 <- p_spot - p_bs
  df_muhs300_30 <- p_mu - p_bs
  # summary(df_hs300_30)
  # summary(df_muhs300_30)
  # 
  # df_mu_perc <- df_hs300_30 / p_spot
  # df_bs_perc <- df_muhs300_30 / p_spot
  
  
  # hist(df_hs300_30)
  # hist(df_muhs300_30)
  # 
  # cor.test(p_spot,p_bs,method='pearson')
  # cor.test(p_spot,p_mu,method='pearson')
  # cor.test(df_hs300_30,df_muhs300_30,method = 'pearson')
  
  
  reg_IF_30 <- lm(df_hs300_30~data_IF$X+df_muhs300_30)  #regression
  # summary(reg_IF_30)
  reg_IF_30_noiv <- lm(df_hs300_30~df_muhs300_30)
  # summary(reg_IF_30_noiv)
  
  # tidy(reg_IF_30)
  # tidy(reg_IF_30_noiv)
  
  coefficient_IF_iv_dfmu[rolldate/10] <- as.numeric(reg_IF_30$coefficients[3])
  coefficient_IF_noiv_dfmu[rolldate/10] <- as.numeric(reg_IF_30_noiv$coefficients[2])
  coefficient_IF_iv_X[rolldate/10] <- as.numeric(reg_IF_30$coefficients[2])
  coefficient_IF_iv_intercept[rolldate/10] <- as.numeric(reg_IF_30$coefficients[1])
  coefficient_IF_noiv_intercept[rolldate/10] <- as.numeric(reg_IF_30_noiv$coefficients[1])
  R_squared_IF_iv[rolldate/10] <- as.numeric(summary(reg_IF_30)$adj.r.squared)
  R_squared_IF_noiv[rolldate/10] <- as.numeric(summary(reg_IF_30_noiv)$adj.r.squared)
  
  # reg_IF_perc <- lm(df_bs_perc~data_IF$X+df_mu_perc)  #regression
  # summary(reg_IF_perc)
  
  # sink()
}

# IH
for (rolldate in seq(10,360,10))
{
  rmu_sz50_30 <- filter(index_sz50$index_sz50, rep(1,rolldate)/rolldate,method="convolution",sides =1)
  rmu_sz50_30 <- na.omit(rmu_sz50_30)
  rmu_sz50_30 <- data.frame(index_sz50$date[rolldate:2403],rmu_sz50_30)
  r_musz50_reg_30 <- rmu_sz50_30[match(data_IH$Trddt,rmu_sz50_30[,1]),2]
  
  
  p_bs <- data_IH$Spot*exp((data_IH$rf_1m-divd_reg_IH)*data_IH$dt)   #p_bs = s0*e^(rf-q)dt
  # p_mu <- data_IH$Spot * exp((mean(r_musz50_reg_30)-divd_reg_IH)*data_IH$dt) #p_mu = s0*e^(rmu-q)dt
  p_mu <- data_IH$Spot * exp((r_musz50_reg_30-divd_reg_IH)*data_IH$dt) #p_mu = s0*e^(rmu-q)dt
  p_spot <- data_IH$Clsprc
  
  df_sz50_30 <- p_spot - p_bs
  df_musz50_30 <- p_mu - p_bs
  # summary(df_sz50_30)
  # summary(df_musz50_30)
  # 
  # df_mu_perc <- df_sz50_30 / p_spot
  # df_bs_perc <- df_musz50_30 / p_spot
  
  
  # hist(df_sz50_30)
  # hist(df_musz50_30)
  # 
  # cor.test(p_spot,p_bs,method='pearson')
  # cor.test(p_spot,p_mu,method='pearson')
  # cor.test(df_sz50_30,df_musz50_30,method = 'pearson')
  
  
  reg_IH_30 <- lm(df_sz50_30~data_IH$X+df_musz50_30)  #regression
  # summary(reg_IH_30)
  reg_IH_30_noiv <- lm(df_sz50_30~df_musz50_30)
  # summary(reg_IH_30_noiv)
  
  # tidy(reg_IH_30)
  # tidy(reg_IH_30_noiv)
  
  coefficient_IH_iv_dfmu[rolldate/10] <- as.numeric(reg_IH_30$coefficients[3])
  coefficient_IH_noiv_dfmu[rolldate/10] <- as.numeric(reg_IH_30_noiv$coefficients[2])
  coefficient_IH_iv_X[rolldate/10] <- as.numeric(reg_IH_30$coefficients[2])
  coefficient_IH_iv_intercept[rolldate/10] <- as.numeric(reg_IH_30$coefficients[1])
  coefficient_IH_noiv_intercept[rolldate/10] <- as.numeric(reg_IH_30_noiv$coefficients[1])
  R_squared_IH_iv[rolldate/10] <- as.numeric(summary(reg_IH_30)$adj.r.squared)
  R_squared_IH_noiv[rolldate/10] <- as.numeric(summary(reg_IH_30_noiv)$adj.r.squared)
  
  # reg_IH_perc <- lm(df_bs_perc~data_IH$X+df_mu_perc)  #regression
  # summary(reg_IH_perc)
  
  # sink()
}

# all
divd_reg_all <- c(divd_reg_IC, divd_reg_IF, divd_reg_IH)

for (rolldate in seq(10,360,10))
{
  data_IC$r_mu_all <- NULL
  data_IF$r_mu_all <- NULL
  data_IH$r_mu_all <- NULL
  
  rmu_zz500_30 <- filter(index_zz500$index_zz500, rep(1,rolldate)/rolldate,method="convolution",sides =1)
  rmu_zz500_30 <- na.omit(rmu_zz500_30)
  rmu_zz500_30 <- data.frame(index_zz500$date[rolldate:2403],rmu_zz500_30)
  r_muzz500_reg_30 <- rmu_zz500_30[match(data_IC$Trddt,rmu_zz500_30[,1]),2]
  
  rmu_hs300_30 <- filter(index_hs300$index_hs300, rep(1,rolldate)/rolldate,method="convolution",sides =1)
  rmu_hs300_30 <- na.omit(rmu_hs300_30)
  rmu_hs300_30 <- data.frame(index_hs300$date[rolldate:2403],rmu_hs300_30)
  r_muhs300_reg_30 <- rmu_hs300_30[match(data_IF$Trddt,rmu_hs300_30[,1]),2]
  
  
  rmu_sz50_30 <- filter(index_sz50$index_sz50, rep(1,rolldate)/rolldate,method="convolution",sides =1)
  rmu_sz50_30 <- na.omit(rmu_sz50_30)
  rmu_sz50_30 <- data.frame(index_sz50$date[rolldate:2403],rmu_sz50_30)
  r_musz50_reg_30 <- rmu_sz50_30[match(data_IH$Trddt,rmu_sz50_30[,1]),2]
  
  
  data_IC <- data.frame(data_IC, r_muzz500_reg_30)
  data_IF <- data.frame(data_IF, r_muhs300_reg_30)
  data_IH <- data.frame(data_IH, r_musz50_reg_30)
  
  names(data_IC)[15] <- "r_mu_all"
  names(data_IF)[15] <- "r_mu_all"
  names(data_IH)[15] <- "r_mu_all"
  
  data_all <- data.frame(rbind(data_IC,data_IF,data_IH))

  p_bs <- data_all$Spot*exp((data_all$rf_1m-divd_reg_all)*data_all$dt)   #p_bs = s0*e^(rf-q)dt
  # p_mu <- data_all$Spot * exp((mean(r_muall_reg_30)-divd_reg_all)*data_all$dt) #p_mu = s0*e^(rmu-q)dt
  p_mu <- data_all$Spot * exp((data_all$r_mu_all-divd_reg_all)*data_all$dt) #p_mu = s0*e^(rmu-q)dt
  p_spot <- data_all$Clsprc
  
  df_all_30 <- p_spot - p_bs
  df_muall_30 <- p_mu - p_bs
  # summary(df_all_30)
  # summary(df_muall_30)
  # 
  # df_mu_perc <- df_all_30 / p_spot
  # df_bs_perc <- df_muall_30 / p_spot
  
  
  # hist(df_all_30)
  # hist(df_muall_30)
  # 
  # cor.test(p_spot,p_bs,method='pearson')
  # cor.test(p_spot,p_mu,method='pearson')
  # cor.test(df_all_30,df_muall_30,method = 'pearson')
  
  
  reg_all_30 <- lm(df_all_30~data_all$X+df_muall_30)  #regression
  # summary(reg_all_30)
  reg_all_30_noiv <- lm(df_all_30~df_muall_30)
  # summary(reg_all_30_noiv)
  
  # tidy(reg_all_30)
  # tidy(reg_all_30_noiv)
  
  coefficient_all_iv_dfmu[rolldate/10] <- as.numeric(reg_all_30$coefficients[3])
  coefficient_all_noiv_dfmu[rolldate/10] <- as.numeric(reg_all_30_noiv$coefficients[2])
  coefficient_all_iv_X[rolldate/10] <- as.numeric(reg_all_30$coefficients[2])
  coefficient_all_iv_intercept[rolldate/10] <- as.numeric(reg_all_30$coefficients[1])
  coefficient_all_noiv_intercept[rolldate/10] <- as.numeric(reg_all_30_noiv$coefficients[1])
  R_squared_all_iv[rolldate/10] <- as.numeric(summary(reg_all_30)$adj.r.squared)
  R_squared_all_noiv[rolldate/10] <- as.numeric(summary(reg_all_30_noiv)$adj.r.squared)
  
  # reg_all_perc <- lm(df_bs_perc~data_all$X+df_mu_perc)  #regression
  # summary(reg_all_perc)
  
  # sink()
  
  
}

# plot
{plot1 <- ggplot(data = data.frame(a,coefficient_IC_iv_dfmu), mapping = aes(x= a, y= coefficient_IC_iv_dfmu))+
  geom_point(alpha = 0.2, color = "black")+
  geom_smooth(aes(x = a,y = coefficient_IC_iv_dfmu),color = "black")+
  labs(x='time scope',y='coefficient (dfmu)',title="coefficient vs time (zz500)")

plot2 <- ggplot(data = data.frame(a,coefficient_IF_iv_dfmu), mapping = aes(x= a, y= coefficient_IF_iv_dfmu))+
  geom_point(alpha = 0.2, color = "black")+
  geom_smooth(aes(x = a,y = coefficient_IF_iv_dfmu),color = "black")+
  labs(x='time scope',y='coefficient (dfmu)',title="coefficient vs time (hs300)")

plot3 <- ggplot(data = data.frame(a,coefficient_IH_iv_dfmu), mapping = aes(x= a, y= coefficient_IH_iv_dfmu))+
  geom_point(alpha = 0.2, color = "black")+
  geom_smooth(aes(x = a,y = coefficient_IH_iv_dfmu),color = "black")+
  labs(x='time scope',y='coefficient (dfmu)',title="coefficient vs time (sz50)")

plot4 <- ggplot(data = data.frame(a,R_squared_IC_iv), mapping = aes(x= a, y= R_squared_IC_iv))+
  geom_point(alpha = 0.2, color = "black")+
  geom_smooth(aes(x = a,y = R_squared_IC_iv),color = "black")+
  labs(x='time scope',y='R^2',title="R^2 vs time (zz500)")

plot5 <- ggplot(data = data.frame(a,R_squared_IF_iv), mapping = aes(x= a, y= R_squared_IF_iv))+
  geom_point(alpha = 0.2, color = "black")+
  geom_smooth(aes(x = a,y = R_squared_IF_iv),color = "black")+
  labs(x='time scope',y='R^2',title="R^2 vs time (hs300)")

plot6 <- ggplot(data = data.frame(a,R_squared_IH_iv), mapping = aes(x= a, y= R_squared_IH_iv))+
  geom_point(alpha = 0.2, color = "black")+
  geom_smooth(aes(x = a,y = R_squared_IH_iv),color = "black")+
  labs(x='time scope',y='R^2',title="R^2 vs time (sz50)")

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
print(plot4,vp = subvp_4)
print(plot5,vp = subvp_5)
print(plot6,vp = subvp_6)


plot7 <- ggplot(data = data.frame(a,coefficient_all_iv_dfmu), mapping = aes(x= a, y= coefficient_all_iv_dfmu))+
  geom_point(alpha = 0.2, color = "black")+
  geom_smooth(aes(x = a,y = coefficient_all_iv_dfmu),color = "black")+
  labs(x='time scope',y='coefficient (dfmu)',title="coefficient vs time (all contracts)")

plot8 <- ggplot(data = data.frame(a,R_squared_all_iv), mapping = aes(x= a, y= R_squared_all_iv))+
  geom_point(alpha = 0.2, color = "black")+
  geom_smooth(aes(x = a,y = R_squared_all_iv),color = "black")+
  labs(x='time scope',y='R^2',title="R^2 vs time (all contracts)")

print(plot7)
print(plot8)
}

reg_results <- data.frame(coefficient_IC_iv_dfmu,coefficient_IF_iv_dfmu,coefficient_IH_iv_dfmu,coefficient_IC_iv_X,coefficient_IF_iv_X,coefficient_IH_iv_X,coefficient_IC_iv_intercept,coefficient_IF_iv_intercept,coefficient_IH_iv_intercept,coefficient_all_iv_dfmu,coefficient_all_iv_X,coefficient_all_iv_intercept,coefficient_IC_noiv_dfmu,coefficient_IF_noiv_dfmu,coefficient_IH_noiv_dfmu,coefficient_IC_noiv_intercept,coefficient_IF_noiv_intercept,coefficient_IH_noiv_intercept,coefficient_all_noiv_dfmu,coefficient_all_noiv_intercept,R_squared_IC_iv,R_squared_IF_iv,R_squared_IH_iv,R_squared_all_iv,R_squared_IC_noiv,R_squared_IF_noiv,R_squared_IH_noiv,R_squared_all_noiv)
write.table(reg_results, file = "reg results", append = FALSE, quote = TRUE, sep = " ", eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = c("escape", "double"))
