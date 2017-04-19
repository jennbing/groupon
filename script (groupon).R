# this script computes the numerical results to compare the profitability
# of Groupon and Groupon Now model using a time - continuum generalization of  
# two-period groupon model developed by [Edelman et al. 2014, To groupon or not to 
# groupon - the profitability of deep discounts].

# install library packages
install.packages("rootSolve")
install.packages("abind")
install.packages("pracma")

# library
library(rootSolve)
library(abind)
library(pracma)

# function R codes (Please change the filename accordingly)
source("../groupons.R")

# output files for plots (Please change the filenames accordingly)
groupon_file="../groupon-deals.csv"
livingsocial_file="../ls-deals.csv"
pc_file="../groupon plots/profit_contour.png"
pdvk_file="../groupon plots/profit_disc_vouc_kappa.png"
pdv2p_file="../groupon plots/profit_disc_vouc_two_popu.png"
pdv2_file="../groupon plots/profit_disc_vouc_two.png"
rfc_file="../groupon plots/r_F_combined.png"
gaoc_file="../groupon plots/gamma_alpha_op_combined.png"
pdvc_file="../groupon plots/profit_disc_vouc_combined.png"
hs_file="../groupon plots/hist_sales.png"

# model parameters
V_G=10 ## decay constant of population G exponential distribution
V_F=20 ## ... population F ...
lambda=0.7 ## proportion of consumer from G population with access to voucher 
		## service
marg_cost=2 ## firm product or service marginal cost
beta=0.5 ## beta, firm fraction of revenue received using the voucher service
# omega=0.2 ## omega, fraction of consumers not redeem the discount vouchers
r_G0=0.5 ## max. fraction of increase in return probability of population G 
		## thanks to discount voucher
r_F0=0.2 ## ... population F ...
T_k=24*7 ## time scale of information spread due to long-term advertising and 
		## social effect in hours
T=24*364 ## duration of study (hours in a day * days in a year)
r_G=0.3 ## return probability of G population

# phase space of profit against gamma and alpha
r_F=matrix(c(0.6,0.3),ncol=1) ## F population return probability
r_G0=matrix(c(0.5,0.8),ncol=1) ## max. fraction of increase in return probability of the G population thanks to discount voucher
r_F0=matrix(c(0.2,0.2),ncol=1) ## max. fraction of increase in return probability of the F population thanks to discount voucher
profit_disc_vouc_all=matrix(0,nrow=length(seq(0.025,0.6,0.025)),ncol=length(seq(0.025,1,0.025)))
png(pc_file,width=5,height=5,units="in",res=300)
par(mfrow=c(2,2),oma = c(0,0,2,0),mar=c(1,1,3,1),mgp=c(2,1,0))
for (i in (1:2)){
	firm_norm_op_day=firm_norm_price_op2(lambda,r_G,r_F[i],V_G,V_F,marg_cost)
	price_op_day=firm_norm_op_day$price
for (j in (1:2)){
	index2=1
	for (k in (seq(0.025,1,0.025))){
		index1=1
	for (l in (seq(0.025,0.6,0.025))){
		gamma1=matrix(k,1,1)
		alpha1=matrix(l,1,1)
		if (alpha1<0) beta=1
		if (alpha1>=0) beta=0.8
		profit_disc_vouc_all[index1,index2]=profit_disc_vouc(lambda,r_G0[j],r_F0[j],r_G,r_F[i],V_G,V_F,price_op_day,marg_cost,gamma1,alpha1,beta)
		index1=index1+1}
		index2=index2+1}
	profit_norm=profit_disc_vouc(lambda,r_G0[j],r_F0[j],r_G,r_F[i],V_G,V_F,price_op_day,marg_cost,matrix(0,1,1),matrix(1,1,1),beta)
	profit_disc_vouc_all=profit_disc_vouc_all/profit_norm
	contour(seq(0.025,0.6,0.025),seq(0.025,1,0.025),profit_disc_vouc_all,levels=seq(1,30,1),main=bquote(r^F ~ "=" ~ .(r_F[i]) ~ ", " ~ r[0]^G ~ "=" ~ .(r_G0[j]) ~ ", " ~ r[0]^F ~ "=" ~ .(r_F0[j])))
	dummy=which(profit_disc_vouc_all==max(profit_disc_vouc_all),arr.ind=TRUE)
	dummy[1]=seq(0.025,0.6,0.025)[dummy[1]]
	dummy[2]=seq(0.025,1,0.025)[dummy[2]]
	points(dummy[1],dummy[2],pch=3,col="red",cex=2)}}
mtext(text=expression("Firm profit " ~ pi[gamma~","~alpha] ~ " contour, " ~ "y-axis: " ~ gamma ~ ", x-axis: " ~ alpha),outer = TRUE)
dev.off()

# phase space for firms with different popularity 
kappa=matrix(c(0.8,0.2),ncol=1) ## F population return probability
r_F=matrix(c(0.3,0.6),ncol=1)
r_G0=matrix(0.5,1,1) ## max. fraction of increase in return probability of the G population thanks to discount voucher
r_F0=matrix(0.2,1,1) ## max. fraction of increase in return probability of the F population thanks to discount voucher
profit_disc_vouc_all=matrix(0,nrow=length(seq(0.025,0.6,0.025)),ncol=length(seq(0.025,1,0.025)))
png(pdvk_file,width=5,height=5,units="in",res=300)
par(mfrow=c(2,2),oma = c(0,0,2,0),mar=c(1,1,3,1),mgp=c(2,1,0))
for (i in (1:2)){
for (j in (1:2)){
	firm_norm_op_day=firm_norm_price_op2(lambda,r_G,r_F[j],V_G,V_F,marg_cost)
	price_op_day=firm_norm_op_day$price
	index2=1
	for (k in (seq(0.025,1,0.025))){
		index1=1
	for (l in (seq(0.025,0.6,0.025))){
		gamma1=matrix(k,1,1)
		alpha1=matrix(l,1,1)
		kappa1=matrix(kappa[i],1,1)
		if (alpha1<0) beta=1
		if (alpha1>=0) beta=0.8
		profit_disc_vouc_all[index1,index2]=profit_disc_vouc_kappa(lambda,r_G0,r_F0,r_G,r_F[j],V_G,V_F,price_op_day,marg_cost,gamma1,alpha1,beta,kappa1)
		index1=index1+1}
		index2=index2+1}
	profit_norm=profit_disc_vouc_kappa(lambda,r_G0,r_F0,r_G,r_F[j],V_G,V_F,price_op_day,marg_cost,matrix(0,1,1),matrix(1,1,1),beta,kappa1)
	profit_disc_vouc_all=profit_disc_vouc_all/matrix(profit_norm,nrow=nrow(profit_disc_vouc_all),ncol=ncol(profit_disc_vouc_all))
	contour(seq(0.025,0.6,0.025),seq(0.025,1,0.025),profit_disc_vouc_all,levels=seq(1,30,1),main=bquote(kappa ~ "=" ~ .(kappa1) ~ ", " ~ r^F ~ "=" ~ .(r_F[j])))
	dummy=which(profit_disc_vouc_all==max(profit_disc_vouc_all),arr.ind=TRUE)
	dummy[1]=seq(0.025,0.6,0.025)[dummy[1]]
	dummy[2]=seq(0.025,1,0.025)[dummy[2]]
	points(dummy[1],dummy[2],pch=3,col="red",cex=2)}}
mtext(text=expression("Firm profit " ~ pi[kappa] ~ " contour, " ~ "y-axis: " ~ gamma ~ ", x-axis: " ~ alpha),outer = TRUE)
dev.off()

# phase space for two sets of discount rates (taking account of popularity)
kappa=matrix(c(0.8,0.2),ncol=1) ## F population return probability
r_F=matrix(c(0.3,0.6),ncol=1)
r_G0=matrix(0.5,1,1) ## max. fraction of increase in return probability of the G population thanks to discount voucher
r_F0=matrix(0.2,1,1) ## max. fraction of increase in return probability of the F population thanks to discount voucher
alpha1=matrix(c(0.6,0.3),ncol=1)
profit_disc_vouc_all=matrix(0,nrow=length(seq(0.025,1,0.025)),ncol=length(seq(0.025,1,0.025)))
png(pdv2p_file,width=5,height=5,units="in",res=300)
par(mfrow=c(2,2),oma = c(0,0,2,0),mar=c(1,1,3,1),mgp=c(2,1,0))
for (i in (1:2)){
for (j in (1:2)){
	firm_norm_op_day=firm_norm_price_op2(lambda,r_G,r_F[j],V_G,V_F,marg_cost)
	price_op_day=firm_norm_op_day$price
	index2=1
	for (k in (seq(0.025,1,0.025))){
		index1=1
	for (l in (seq(0.025,1,0.025))){
		if ((k+l)>1) next
		gamma1=matrix(c(k,l),ncol=1)
		kappa1=matrix(kappa[i],1,1)
		profit_disc_vouc_all[index1,index2]=profit_disc_vouc_kappa(lambda,r_G0,r_F0,r_G,r_F[j],V_G,V_F,price_op_day,marg_cost,gamma1,alpha1,beta,kappa1)
		index1=index1+1}
		index2=index2+1}
	profit_norm=profit_disc_vouc_kappa(lambda,r_G0,r_F0,r_G,r_F[j],V_G,V_F,price_op_day,marg_cost,matrix(0,1,1),matrix(1,1,1),beta,kappa1)
	profit_disc_vouc_all=profit_disc_vouc_all/matrix(profit_norm,nrow=nrow(profit_disc_vouc_all),ncol=ncol(profit_disc_vouc_all))
	profit_disc_vouc_all[profit_disc_vouc_all==0]=NaN
	contour(seq(0.025,1,0.025),seq(0.025,1,0.025),profit_disc_vouc_all,levels=seq(1,30,1),main=bquote(kappa ~ "=" ~ .(kappa1) ~ ", " ~ r^F ~ "=" ~ .(r_F[j])))
	dummy=which(profit_disc_vouc_all==max(profit_disc_vouc_all,na.rm=TRUE),arr.ind=TRUE)
	dummy[1]=seq(0.025,1,0.025)[dummy[1]]
	dummy[2]=seq(0.025,1,0.025)[dummy[2]]
	points(dummy[1],dummy[2],pch=3,col="red",cex=2)
	points(c(0,1),c(1,0),type="l",col="blue")}}
text=bquote("Firm profit " ~ pi[kappa] ~ " contour, " ~ "y-axis: " ~ gamma[1] ~"("~ alpha[1] ~ "=" ~ .(alpha1[1]) ~"), x-axis: " ~ gamma[2] ~"("~ alpha[2] ~ "=" ~ .(alpha1[2]) ~ ")")
mtext(text=text,outer = TRUE)
dev.off()

# phase space for two sets of discount rates
r_F=matrix(c(0.6,0.3),ncol=1) ## F population return probability
r_G0=matrix(c(0.5,0.8),ncol=1) ## max. fraction of increase in return probability of the G population thanks to discount voucher
r_F0=matrix(c(0.2,0.2),ncol=1) ## max. fraction of increase in return probability of the F population thanks to discount voucher
alpha1=matrix(c(0.6,0.3),ncol=1) ## discount rates
profit_disc_vouc_all=matrix(0,nrow=length(seq(0.025,1,0.025)),ncol=length(seq(0.025,1,0.025)))
png(pdv2_file,width=5.5,height=5,units="in",res=300)
par(mfrow=c(2,2),oma = c(0,0,2,0),mar=c(1,1,3,1),mgp=c(2,1,0))
for (i in (1:2)){
	firm_norm_op_day=firm_norm_price_op2(lambda,r_G,r_F[i],V_G,V_F,marg_cost)
	price_op_day=firm_norm_op_day$price
for (j in (1:2)){
	index2=1
	for (k in (seq(0.025,1,0.025))){
		index1=1
	for (l in (seq(0.025,1,0.025))){
		if ((k+l)>1) next
		gamma1=matrix(c(k,l),ncol=1)
		profit_disc_vouc_all[index1,index2]=profit_disc_vouc(lambda,r_G0[j],r_F0[j],r_G,r_F[i],V_G,V_F,price_op_day,marg_cost,gamma1,alpha1,beta)
		index1=index1+1}
		index2=index2+1}
	profit_norm=profit_disc_vouc(lambda,r_G0[j],r_F0[j],r_G,r_F[i],V_G,V_F,price_op_day,marg_cost,matrix(0,1,1),matrix(1,1,1),beta)
	profit_disc_vouc_all=profit_disc_vouc_all/matrix(profit_norm,nrow=nrow(profit_disc_vouc_all),ncol=ncol(profit_disc_vouc_all))
	profit_disc_vouc_all[profit_disc_vouc_all==0]=NaN
	contour(seq(0.025,1,0.025),seq(0.025,1,0.025),profit_disc_vouc_all,levels=seq(1,30,1),main=bquote(r^F ~ "=" ~ .(r_F[i]) ~ ", " ~ r[0]^G ~ "=" ~ .(r_G0[j]) ~ ", " ~ r[0]^F ~ "=" ~ .(r_F0[j])))
	dummy=which(profit_disc_vouc_all==max(profit_disc_vouc_all,na.rm=TRUE),arr.ind=TRUE)
	dummy[1]=seq(0.025,1,0.025)[dummy[1]]
	dummy[2]=seq(0.025,1,0.025)[dummy[2]]
	points(dummy[1],dummy[2],pch=3,col="red",cex=2)
	points(c(0,1),c(1,0),type="l",col="blue")}}
text=bquote("Firm profit " ~ pi[gamma~","~alpha] ~ " contour, " ~ "y-axis: " ~ gamma[1] ~"("~ alpha[1] ~ "=" ~ .(alpha1[1]) ~"), x-axis: " ~ gamma[2] ~"("~ alpha[2] ~ "=" ~ .(alpha1[2]) ~ ")")
mtext(text=text,outer=TRUE)
dev.off()

## new return probability of F population 
r_F_ave=matrix(c(0.41,0.01),ncol=1)
r_F_day=abind(matrix(0,ncol=3*60),0.095+0.095*cos(8*pi*matrix(seq(1,18*60,1),
	nrow=1)/(24*60)-pi),matrix(0,ncol=3*60),along=2)
r_F_week=0.2*cos(pi*matrix(seq(1,24*60,1),nrow=1)/(24*60)-pi/2)
r_F_day=matrix(r_F_day,nrow=1)
r_F_week=matrix(r_F_week,nrow=1)

## plots
png(rfc_file,width=3.3,height=3.3,units="in",res=300)
plot(t(r_F_ave[2]+r_F_day+r_F_week),type="l",col="black",main=expression("F population return probability, "~ r^{F}),ylab=expression(r^{F}),xlab="time (min) ",ylim=c(0,0.8))
points(t(r_F_ave[1]+r_F_day+r_F_week),type="l",col="black")
#points(t(r_F_ave[1]+r_F_day+r_F_week),type="l",col="black")
dev.off()

# Firm optimal alpha and optimal
r_G0=matrix(c(0.5,0.8),ncol=1) ## max. fraction of increase in return probability of the G population thanks to discount voucher
r_F0=matrix(c(0.2,0.2),ncol=1) ## max. fraction of increase in return probability of the F population thanks to discount voucher
png(gaoc_file,width=5.5,height=5,units="in",res=300)
par(mfrow=c(2,2),oma = c(0,0,2,0),mar=c(2.5,2,3,2.5),mgp=c(2,1,0))
gamma=matrix(0.2,1,1)
for (i in (1:2)){
for (j in (1:2)){
	firm_norm_op_day=firm_norm_price_op2(lambda,r_G,r_F_ave[i]+r_F_day+r_F_week,V_G,V_F,marg_cost)
	price_op_day=firm_norm_op_day$price
	alpha_op_day=firm_alpha_op(lambda,r_G0[j],r_F0[j],r_G,r_F_ave[i]+r_F_day+r_F_week,V_G,V_F,price_op_day,marg_cost,gamma,beta,2)
	alpha_op_day1=firm_alpha_op(lambda,r_G0[j],r_F0[j],r_G,r_F_ave[i]+r_F_day+r_F_week,V_G,V_F,price_op_day,marg_cost,gamma,beta,1)
	gamma_op_day=firm_gamma_op(lambda,r_G0[j],r_F0[j],r_G,r_F_ave[i]+r_F_day+r_F_week,V_G,V_F,price_op_day,marg_cost,alpha_op_day1,beta,2)	
	peak=max(r_F_ave[i]+r_F_day+r_F_week)
	plot(t(alpha_op_day),type="l",lty=2,col="red",ylab="",ylim=c(0.2,0.35),xlab="",main=bquote(r[peak]^F ~ "=" ~ .(peak) ~ ", " ~ r[0]^G ~ "=" ~ .(r_G0[j]) ~ ", " ~ r[0]^F ~ "=" ~ .(r_F0[j]))) 
	points(matrix(alpha_op_day1,nrow=length(alpha_op_day),ncol=1),type="l",col="red")
	par(new=TRUE)
	plot(t(gamma_op_day),type="l",col="blue",xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(0,0.5))
	axis(side=4)}} #,at=c(seq(0,max(gamma_op_day),max(gamma_op_day)/3)))}}#,labels=as.character(c(6e-5,9e-5,12e-5)))}}
mtext(text=expression("Optimal " ~ alpha ~ "* (red, left y-axis) and" ~ gamma ~ "* (blue, right y-axis) over time (min)"), outer=TRUE)
dev.off()

# Firm normal profit and profit with the issuance of discount vouchers (optimal alpha and gamma)
r_G0=matrix(c(0.5,0.8),ncol=1) ## max. fraction of increase in return probability of the G population thanks to discount voucher
r_F0=matrix(c(0.2,0.2),ncol=1) ## max. fraction of increase in return probability of the F population thanks to discount voucher
png(pdvc_file,width=5,height=5,units="in",res=300)
par(mfrow=c(2,2),oma = c(0,0,2,0),mar=c(1,1,3,1),mgp=c(2,1,0))
gamma=matrix(0.2,1,1)
for (i in (1:2)){
for (j in (1:2)){
	firm_norm_op_day=firm_norm_price_op2(lambda,r_G,r_F_ave[i]+r_F_day+r_F_week,V_G,V_F,marg_cost)
	price_op_day=firm_norm_op_day$price
	alpha_op_day=firm_alpha_op(lambda,r_G0[j],r_F0[j],r_G,r_F_ave[i]+r_F_day+r_F_week,V_G,V_F,price_op_day,marg_cost,gamma,beta,2)
	alpha_op_day1=firm_alpha_op(lambda,r_G0[j],r_F0[j],r_G,r_F_ave[i]+r_F_day+r_F_week,V_G,V_F,price_op_day,marg_cost,gamma,beta,1)
	gamma_op_day=firm_gamma_op(lambda,r_G0[j],r_F0[j],r_G,r_F_ave[i]+r_F_day+r_F_week,V_G,V_F,price_op_day,marg_cost,alpha_op_day1,beta,2)	
	profit_disc_vouc_day=profit_disc_vouc(lambda,r_G0[j],r_F0[j],r_G,r_F_ave[i]+r_F_day+r_F_week,V_G,V_F,price_op_day,marg_cost,matrix(0,1,1),matrix(1,1,1),beta)
	profit_disc_vouc_day1=profit_disc_vouc(lambda,r_G0[j],r_F0[j],r_G,r_F_ave[i]+r_F_day+r_F_week,V_G,V_F,price_op_day,marg_cost,gamma,alpha_op_day1,beta)
	profit_disc_vouc_day2=profit_disc_vouc(lambda,r_G0[j],r_F0[j],r_G,r_F_ave[i]+r_F_day+r_F_week,V_G,V_F,price_op_day,marg_cost,gamma_op_day,alpha_op_day1,beta)
	peak=max(r_F_ave[i]+r_F_day+r_F_week)
	plot(t(profit_disc_vouc_day),type="l",col="black",main=bquote(r[peak]^F ~ "=" ~ .(peak) ~ ", " ~ r[0]^G ~ "=" ~ .(r_G0[j]) ~ ", " ~ r[0]^F ~ "=" ~ .(r_F0[j])),ylab="",xlab="",ylim=c(0,2)) 
	points(t(profit_disc_vouc_day1),type="l",col="red")
	points(t(profit_disc_vouc_day2),type="l",col="blue")}}
mtext(text=expression("Firm profit " ~ pi ~ ", "~ pi[gamma~","~alpha] ~ "(" ~ alpha ~ "*), " ~ pi[gamma~","~alpha] ~ "(" ~ gamma ~ "*)" ~ " over time"), outer=TRUE)
dev.off()


data=read.csv(groupon_file,header=TRUE)
dist_price=data[,3]-data[,4]
num_sold=data[,6]
freq=matrix(0,36)
count=0
price=5+((0:35)*10)
for (i in 1:35) {
	count=count+1
	freq[count]=sum(num_sold[(dist_price-price[i])<5 & (price[i]-dist_price)<=5])
}
count=count+1
freq[count]=sum(num_sold[dist_price>=350])
freq=freq/sum(num_sold)
png(hs_file,width=8,height=5,units="in",res=300)
par(mfrow=c(1,2))
plot(price,freq,log="y",col="blue",pch=19,xlab="discounted price",ylab="probability",ylim=c(1e-5,1))
log_est=summary(lm(log(freq[freq!=0])~price[freq!=0]))
log_est
log_est$coefficients[1]
log_est$coefficients[2]
lines(price,exp(price*log_est$coefficients[2]+ log_est$coefficients[1]),col="red",yaxt="n")
legend("topright",c("data","logarithmic fit"),col=c("blue","red"),lty=c(0,1),pch=c(19,NA_integer_),cex=0.8)
title("Groupon Histogram of Sales")

data=read.csv(livingsocial_file,header=TRUE)
dist_price=data[,4]-data[,5]
num_sold=data[,7]
freq=matrix(0,36)
count=0
price=5+((0:35)*10)
for (i in 1:35) {
	count=count+1
	freq[count]=sum(num_sold[(dist_price-price[i])<5 & (price[i]-dist_price)<=5])
}
count=count+1
freq[count]=sum(num_sold[dist_price>=350])
freq=freq/sum(num_sold)
plot(price,freq,log="y",col="blue",pch=19,xlab="discounted price",ylab="probability",ylim=c(1e-5,1))
log_est=summary(lm(log(freq[freq!=0])~price[freq!=0]))
log_est
log_est$coefficients[1]
log_est$coefficients[2]
lines(price,exp(price*log_est$coefficients[2]+ log_est$coefficients[1]),col="red")
legend("topright",c("data","logarithmic fit"),col=c("blue","red"),lty=c(0,1),pch=c(19,NA_integer_),cex=0.8)
title("LivingSocial Histogram of Sales")
dev.off()