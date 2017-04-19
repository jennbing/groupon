## functions here are used to compute numerical studies using the time - 
## continuum model based on the groupon model developed by 
## [Edelman et al. 2014, To groupon or not to groupon - 
## the profitability of deep discounts].

## ------------ firm normal demand, profit, and optimal price --------------

## firm normal operation demand and profit
# input:
# lambda - proportion of consumers from G population
# r_G - return probability of G population
# r_F - ... of F population
# V_G - decay constant in exponential distribution of G population
# V_F - ... of F population
# price - firm product or service price 
# marg_cost - firm product or service marginal cost
# output:
# demand - firm normal operation demand
# profit - firm normal operation profit
firm_norm=function(lambda,r_G,r_F,V_G,V_F,price,marg_cost){ 
	demand=lambda*r_G*(1-exp_cdf(price/r_G,V_G))+(1-lambda)*r_F*
		(1-exp_cdf(price/r_F,V_F))
	profit=demand*(price-marg_cost)
	out=list("demand"=demand,"profit"=profit)
return(out)}

# price optimization (algorithm 1)
# input:
# lambda - proportion of consumers from G population
# r_G - return probability of G population
# r_F - ... of F population
# V_G - decay constant in exponential distribution of G population
# V_F - ... of F population
# price - firm product or service price 
# marg_cost - firm product or service marginal cost
# output:
# profit_op - firm normal operation optimized profit
# price_op - firm normal operation optimal price
firm_norm_price_op1=function(lambda,r_G,r_F,V_G,V_F,marg_cost){
	price_op=matrix(0,nrow=1,ncol=length(r_F))
	profit_op=0
	for (i in 1:length(r_F)){	
		r_F2=r_F[i]
		firm_norm_deri=function(price,lambda1=lambda,r_G1=r_G,r_F1=r_F2,
			V_G1=V_G,V_F1=V_F,marg_cost1=marg_cost){
			out=(r_G1*lambda1*(1-exp_cdf(price/r_G1,V_G1))+
				r_F1*(1-lambda1)*(1-exp_cdf(price/r_F1,V_F1))+
				(price-marg_cost1)*(-r_G1*lambda1*exp_pdf(price/r_G1,
				V_G1)/r_G1-r_F1*(1-lambda1)*
				exp_pdf(price/r_F1,V_F1)/r_F1))
			return(out)}
		price_op[i]=uniroot.all(firm_norm_deri,c(0,100*marg_cost))[1]}
	for (i in unique(price_op)){
		profit_sum1=sum(firm_norm(lambda,r_G,r_F,V_G,V_F,i,marg_cost)
				$profit)
		if (profit_sum1>profit_op){
			profit_op=profit_sum1
			price_op1=i}}
	firm_norm_op=list("profit"=profit_op,"price"=price_op1)
	return(firm_norm_op)}

# price optimization (algorithm 2 - faster)
# input:
# lambda - proportion of consumers from G population
# r_G - return probability of G population
# r_F - ... of F population
# V_G - decay constant in exponential distribution of G population
# V_F - ... of F population
# price - firm product or service price 
# marg_cost - firm product or service marginal cost
# output:
# profit_op - firm normal operation optimized profit
# price_op - firm normal operation optimal price
firm_norm_price_op2=function(lambda,r_G,r_F,V_G,V_F,marg_cost){
	price_op=matrix(0,nrow=1,ncol=length(unique(r_F)))
	profit_op=0
	index=1
	for (i in unique(r_F)){	
		firm_norm_deri=function(price,lambda1=lambda,r_G1=r_G,r_F1=i,
			V_G1=V_G,V_F1=V_F,marg_cost1=marg_cost){
			out=(r_G1*lambda1*(1-exp_cdf(price/r_G1,V_G1))+
				r_F1*(1-lambda1)*(1-exp_cdf(price/r_F1,V_F1))+
				(price-marg_cost1)*(-r_G1*lambda1*exp_pdf(price/r_G1,
				V_G1)/r_G1-r_F1*(1-lambda1)*
				exp_pdf(price/r_F1,V_F1)/r_F1))
			return(out)}
		price_op[index]=uniroot.all(firm_norm_deri,c(0,100*marg_cost))[1]
		index=index+1}
	for (i in unique(price_op)){
		profit_sum1=sum(firm_norm(lambda,r_G,r_F,V_G,V_F,i,marg_cost)
				$profit)
		if (profit_sum1>profit_op){
			profit_op=profit_sum1
			price_op1=i}}
	firm_norm_op=list("profit"=profit_op,"price"=price_op1)
	return(firm_norm_op)}

## ------------ firm profit taking account of discount vouchers -------------

## firm profit taking account of discount voucher 
## (this function takes in time series of gamma and alpha)
# input:
# lambda - proportion of consumer from G population 
# r_G0 - max. fractional increase in return probability of G population
# r_F0 - max. fractional increase in return probability of F population
# r_G - return probability of G population
# r_F - ... F population
# V_G - decay constant in exponential distribution of G population
# V_F - ... F population
# price_op - firm normal operation optimal price
# marg_cost - firm product or service marginal cost
# gamma - fraction of consumer taking up the discount vouchers
# alpha - 1 - voucher discount rate
# beta - fraction received by firm using the voucher service
# output:
# profit - firm profit taking account of discount vouchers
# profit_disc_vouc1=function(lambda,r_G0,r_F0,r_G,r_F,V_G,V_F,price_op,marg_cost,
	# gamma,alpha,beta){
	# gamma_deri=diff(abind(gamma[length(gamma)],gamma,along=2),lag=1)
	# alpha_deri=diff(abind(alpha[length(gamma)],alpha,along=2),lag=1)
	# inc_r_G=r_G0*(1-r_G)*(1-exp(-(gamma+1-alpha+abs(gamma_deri-alpha_deri)*
	#	T_k2)))
	# inc_r_F=r_F0*(1-r_F)*(1-exp(-(gamma+1-alpha+abs(gamma_deri-alpha_deri)*
	# 	T_k2)))
	#inc_r_G=r_G0*(1-r_G)*(1-exp(-(gamma[1,]+1-alpha[1,])))
	#inc_r_F=r_F0*(1-r_F)*(1-exp(-(gamma[1,]+1-alpha[1,])))
	#profit=gamma[1,]*lambda*(r_G+inc_r_G)*(alpha[1,]*beta*price_op-marg_cost)*
	#	(1-exp_cdf(alpha[1,]*price_op/(r_G+inc_r_G),V_G))+
	#	(1-gamma[1,])*lambda*(r_G+inc_r_G)*(price_op-marg_cost)*
	#	(1-exp_cdf(price_op/(r_G+inc_r_G),V_G))+(1-lambda)*(r_F+inc_r_F)*
	#	(price_op-marg_cost)*(1-exp_cdf(price_op/(r_F+inc_r_F),V_F))
	#return(profit)}

## firm profit taking account of discount voucher 
## (this function takes in set of time series of gamma and alpha)
# input:
# lambda - proportion of consumer from G population 
# r_G0 - max. fractional increase in return probability of G population
# r_F0 - max. fractional increase in return probability of F population
# r_G - return probability of G population
# r_F - ... F population
# V_G - decay constant in exponential distribution of G population
# V_F - ... F population
# price_op - firm normal operation optimal price
# marg_cost - firm product or service marginal cost
# gamma - fraction of consumer taking up the discount vouchers
# alpha - 1 - voucher discount rate
# beta - fraction received by firm using the voucher service
# output:
# profit - firm profit taking account of discount vouchers
profit_disc_vouc=function(lambda,r_G0,r_F0,r_G,r_F,V_G,V_F,price_op,marg_cost,
	gamma,alpha,beta){
		if (ncol(gamma)==1) gamma=repmat(gamma,1,ncol(alpha))
		if (ncol(alpha)==1) alpha=repmat(alpha,1,ncol(gamma))
		gamma[alpha==1]=0
		alpha[gamma==0]=1
		inc_r_G=r_G0*(1-r_G)*(1-exp(-(apply(gamma,2,sum)+
				apply(1-alpha,2,sum)/nrow(alpha))))
		inc_r_F=r_F0*(1-r_F)*(1-exp(-(apply(gamma,2,sum)+
				apply(1-alpha,2,sum)/nrow(alpha))))
	profit=0
	for (i in (1:nrow(gamma))){
		profit=profit+gamma[i,]*lambda*(r_G+inc_r_G)*(alpha[i,]*beta*
			price_op-marg_cost)*(1-exp_cdf(alpha[i,]*price_op/(r_G+
			inc_r_G),V_G))}
		profit=profit+(1-apply(gamma,2,sum))*lambda*(r_G+inc_r_G)*
			(price_op-marg_cost)*(1-exp_cdf(price_op/
			(r_G+inc_r_G),V_G))+(1-lambda)*(r_F+inc_r_F)*(price_op-
			marg_cost)*(1-exp_cdf(price_op/(r_F+inc_r_F),V_F))
	return(profit)}

## firm profit taking account of non-redeemed discount voucher
## (this function takes in set of time series of gamma and alpha)
# input:
# lambda - proportion of consumer from G population 
# r_G0 - max. fractional increase in return probability of G population
# r_F0 - max. fractional increase in return probability of F population
# r_G - return probability of G population
# r_F - ... F population
# V_G - decay constant in exponential distribution of G population
# V_F - ... F population
# price_op - firm normal operation optimal price
# marg_cost - firm product or service marginal cost
# gamma - fraction of consumer taking up the discount vouchers
# alpha - 1 - voucher discount rate
# beta - fraction received by firm using the voucher service
# omega - fraction of non-redeemed discount vouchers
# output:
# profit - firm profit taking account of discount vouchers
profit_disc_vouc_omega=function(lambda,r_G0,r_F0,r_G,r_F,V_G,V_F,price_op,marg_cost,
	gamma,alpha,beta,omega){
		if (ncol(gamma)==1) gamma=repmat(gamma,1,ncol(alpha))
		if (ncol(alpha)==1) alpha=repmat(alpha,1,ncol(gamma))
		gamma[alpha==1]=0
		alpha[gamma==0]=1
		inc_r_G=r_G0*(1-r_G)*(1-exp(-(apply(gamma,2,sum)+
				apply(1-alpha,2,sum)/nrow(alpha))))
		inc_r_F=r_F0*(1-r_F)*(1-exp(-(apply(gamma,2,sum)+
				apply(1-alpha,2,sum)/nrow(alpha))))
	profit=0
	for (i in (1:nrow(gamma))){
		profit=profit+gamma[i,]*lambda*(r_G+inc_r_G)*(alpha[i,]*beta*
			price_op-marg_cost)*(1-exp_cdf(alpha[i,]*price_op/(r_G+
			inc_r_G),V_G))+omega*gamma[i,]*lambda*(r_G+inc_r_G)*
			marg_cost*(1-exp_cdf(alpha[i,]*price_op/(r_G+
			inc_r_G),V_G))}
		profit=profit+(1-apply(gamma,2,sum))*lambda*(r_G+inc_r_G)*
			(price_op-marg_cost)*(1-exp_cdf(price_op/
			(r_G+inc_r_G),V_G))+(1-lambda)*(r_F+inc_r_F)*(price_op-
			marg_cost)*(1-exp_cdf(price_op/(r_F+inc_r_F),V_F))
	return(profit)}

## firm profit taking account of popularity index
## (this function takes in set of time series of gamma and alpha)
# input:
# lambda - proportion of consumer from G population 
# r_G0 - max. fractional increase in return probability of G population
# r_F0 - max. fractional increase in return probability of F population
# r_G - return probability of G population
# r_F - ... F population
# V_G - decay constant in exponential distribution of G population
# V_F - ... F population
# price_op - firm normal operation optimal price
# marg_cost - firm product or service marginal cost
# gamma - fraction of consumer taking up the discount vouchers
# alpha - 1 - voucher discount rate
# beta - fraction received by firm using the voucher service
# kappa - fraction of consumers recall the firm product or service
# output:
# profit - firm profit taking account of discount vouchers
profit_disc_vouc_kappa=function(lambda,r_G0,r_F0,r_G,r_F,V_G,V_F,price_op,marg_cost,
	gamma,alpha,beta,kappa){
		if (ncol(gamma)==1) gamma=repmat(gamma,1,ncol(alpha))
		if (ncol(alpha)==1) alpha=repmat(alpha,1,ncol(gamma))
		gamma[alpha==1]=0
		alpha[gamma==0]=1
		inc_r_G=r_G0*(1-r_G)*(1-exp(-(apply(gamma,2,sum)+
				apply(1-alpha,2,sum)/nrow(alpha))))
		inc_r_F=r_F0*(1-r_F)*(1-exp(-(apply(gamma,2,sum)+
				apply(1-alpha,2,sum)/nrow(alpha))))
	profit=0
	for (i in (1:nrow(gamma))){
		profit=profit+gamma[i,]*lambda*(r_G+inc_r_G)*(alpha[i,]*
			beta*price_op-marg_cost)*(1-exp_cdf(alpha[i,]*price_op/
			(r_G+inc_r_G),V_G))}
		profit=profit+kappa[1,]*(1-apply(gamma,2,sum))*lambda*
			(r_G+inc_r_G)*(price_op-marg_cost)*(1-exp_cdf(price_op/
			(r_G+inc_r_G),V_G))+kappa[1,]*(1-lambda)*(r_F+inc_r_F)*
			(price_op-marg_cost)*(1-exp_cdf(price_op/(r_F+inc_r_F),V_F))
	return(profit)}

## firm profit taking account of non-redeemed discount voucher and popularity
## index
## (this function takes in set of time series of gamma and alpha)
# input:
# lambda - proportion of consumer from G population 
# r_G0 - max. fractional increase in return probability of G population
# r_F0 - max. fractional increase in return probability of F population
# r_G - return probability of G population
# r_F - ... F population
# V_G - decay constant in exponential distribution of G population
# V_F - ... F population
# price_op - firm normal operation optimal price
# marg_cost - firm product or service marginal cost
# gamma - fraction of consumer taking up the discount vouchers
# alpha - 1 - voucher discount rate
# beta - fraction received by firm using the voucher service
# omega - fraction of non-redeemed discount vouchers
# kappa - fraction of consumers recall the firm product or service
# output:
# profit - firm profit taking account of discount vouchers
profit_disc_vouc_omega_kappa=function(lambda,r_G0,r_F0,r_G,r_F,V_G,V_F,price_op,marg_cost,
	gamma,alpha,beta,omega,kappa){
		if (ncol(gamma)==1) gamma=repmat(gamma,1,ncol(alpha))
		if (ncol(alpha)==1) alpha=repmat(alpha,1,ncol(gamma))
		gamma[alpha==1]=0
		alpha[gamma==0]=1
		inc_r_G=r_G0*(1-r_G)*(1-exp(-(apply(gamma,2,sum)+
				apply(1-alpha,2,sum)/nrow(alpha))))
		inc_r_F=r_F0*(1-r_F)*(1-exp(-(apply(gamma,2,sum)+
				apply(1-alpha,2,sum)/nrow(alpha))))
	profit=0
	for (i in (1:nrow(gamma))){
		profit=profit+gamma[i,]*lambda*(r_G+inc_r_G)*(alpha[i,]*beta*
			price_op-marg_cost)*(1-exp_cdf(alpha[i,]*price_op/(r_G+
			inc_r_G),V_G))+omega*gamma[i,]*lambda*(r_G+inc_r_G)*
			marg_cost*(1-exp_cdf(alpha[i,]*price_op/(r_G+
			inc_r_G),V_G))}
		profit=profit+kappa[1,]*(1-apply(gamma,2,sum))*lambda*
			(r_G+inc_r_G)*(price_op-marg_cost)*(1-exp_cdf(price_op/
			(r_G+inc_r_G),V_G))+kappa[1,]*(1-lambda)*(r_F+inc_r_F)*
			(price_op-marg_cost)*(1-exp_cdf(price_op/(r_F+inc_r_F),V_F))
	return(profit)}

## popularity index 
# input:
# q - consumer average ratings
# t - time
# T_k - time scale of the advertising and social effect
# output:
# kappa - popularity index, fraction of population recall the firm, product 
# 	    or service
popu_index=function(q,t,T_k){
	kappa=q*(1-exp(-q*t/T_k))
return(kappa)}

# price_re_op - firm price reoptimization with issuance of discount vouchers


## ------------------- firm optimal gamma and alpha -------------------------

## firm optimal alpha (fix a gamma)
## firm operation taking account of discount vouchers
# input:
# lambda - proportion of consumers from G population
# r_G0 - max. fractional increase in return probability of G population thanks
		# to issuance of discount vouchers
# r_F0 - ... F population ...
# r_G - return probability of G population
# r_F - ... of F population
# V_G - decay constant in exponential distribution of G population
# V_F - ... of F population
# price_op - firm product or service normal operation optimal price 
# marg_cost - firm product or service marginal cost
# gamma - fraction of consumers from G population taking up the discount 
		# vouchers 
# beta - fraction of revenue received by firm using the voucher service
# algo - 1: constant optimal alpha over time
#	   2: optimal alpha that varies with time
# output:
# alpha_op - optimal alpha (1 - discount rate) given the gamma 
firm_alpha_op=function(lambda,r_G0,r_F0,r_G,r_F,V_G,V_F,price_op,marg_cost,
	gamma,beta,algo){
	alpha_op=matrix(1,nrow=1,ncol=length(r_F))
	for (i in 1:length(r_F)){	
		r_F2=r_F[i]
		profit_deri_alpha=function(alpha,lambda1=lambda,r_G01=r_G0,
			r_F01=r_F0,r_G1=r_G,r_F1=r_F2,V_G1=V_G,V_F1=V_F,
			price_op1=price_op,marg_cost1=marg_cost,gamma1=gamma,
			beta1=beta){			
			gamma2=matrix(gamma1,ncol=length(alpha))
			gamma2[alpha==1]=0
			inc_r_G1=r_G01*(1-r_G1)*(1-exp(-(gamma2+1-alpha)))
			inc_r_F1=r_F01*(1-r_F1)*(1-exp(-(gamma2+1-alpha)))
			inc_r_G1_deri=-r_G01*(1-r_G1)*exp(-(gamma2+1-alpha))
			inc_r_F1_deri=-r_F01*(1-r_F1)*exp(-(gamma2+1-alpha))
			val_G1=price_op1/(r_G1+inc_r_G1)
			val_F1=price_op1/(r_F1+inc_r_F1)
			val_hat_G1=alpha*price_op1/(r_G1+inc_r_G1)
			val_G1_deri=-price_op1/(r_G1+inc_r_G1)^2*inc_r_G1_deri
			val_F1_deri=-price_op1/(r_F1+inc_r_F1)^2*inc_r_F1_deri
			val_hat_G1_deri=price_op1/(r_G1+inc_r_G1)-alpha*price_op1/
				(r_G1+inc_r_G1)^2*inc_r_G1_deri
			out=gamma2*lambda1*(r_G1+inc_r_G1)*beta1*price_op1*
				(1-exp_cdf(val_hat_G1,V_G1))+	gamma2*lambda1*
				inc_r_G1_deri*(alpha*beta1*price_op1-marg_cost1)*
				(1-exp_cdf(val_hat_G1,V_G1))-gamma2*lambda1*
				(r_G1+inc_r_G1)*(alpha*beta1*price_op1-marg_cost1)*
				exp_pdf(val_hat_G1,V_G1)*val_hat_G1_deri+
				(1-gamma2)*lambda1*inc_r_G1_deri*(price_op1-
				marg_cost1)*(1-exp_cdf(val_G1,V_G1))-(1-gamma2)*
				lambda1*(r_G1+inc_r_G1)*(price_op1-marg_cost1)*
				exp_pdf(val_G1,V_G1)*val_G1_deri+(1-lambda1)*
				inc_r_F1_deri*(price_op1-marg_cost1)*
				(1-exp_cdf(val_F1,V_F1))-(1-lambda1)*(r_F1+inc_r_F1)*
				(price_op1-marg_cost1)*exp_pdf(val_F1,V_F1)*
				val_F1_deri
			return(out)}
		alpha_op1<-tryCatch(uniroot.all(profit_deri_alpha,c(0,1),tol=0.001,maxiter=1000,n=100),error=function(e) 1)
		# print(alpha_op1)
		if ((alpha_op1[length(alpha_op1)]>=0) &&
			(alpha_op1[length(alpha_op1)]<=1) && length(alpha_op1)!=0) {
				alpha_op[i]<-alpha_op1[length(alpha_op1)]}}
	if (algo==1){
		profit_sum=0
		alpha_op1=0
		for (j in unique(alpha_op)){
			alpha=matrix(j,1,1)
			if (alpha==1) gamma1=matrix(0,1,1)
			if (alpha!=1) gamma1=gamma			
			profit_sum1=sum(profit_disc_vouc(lambda,r_G0,r_F0,r_G,r_F,
				V_G,V_F,price_op,marg_cost,gamma1,alpha,beta))
			if (profit_sum1>profit_sum){
				profit_sum=profit_sum1
				alpha_op1=alpha}}
		alpha_op=alpha_op1}
	return(alpha_op)}

## firm optimal alpha (given a gamma)
## firm operation taking account of discount vouchers
## (this algorithm uses the brute force method to search for the alpha that 
## maximizes the profit)
# input:
# lambda - proportion of consumers from G population
# r_G0 - max. fractional increase in return probability of G population thanks
		# to issuance of discount vouchers
# r_F0 - ... F population ...
# r_G - return probability of G population
# r_F - ... of F population
# V_G - decay constant in exponential distribution of G population
# V_F - ... of F population
# price_op - firm product or service normal operation optimal price 
# marg_cost - firm product or service marginal cost
# gamma - fraction of consumers from G population taking up the discount 
		# vouchers 
# beta - fraction of revenue received by firm using the voucher service
# algo - 1: constant optimal alpha over time
#	   2: optimal alpha that varies with time
# output:
# alpha_op - optimal alpha (1 - discount rate) given the gamma 
firm_alpha_op1=function(lambda,r_G0,r_F0,r_G,r_F,V_G,V_F,price_op,marg_cost,
	gamma,beta,algo){
	alpha_op=matrix(1,nrow=1,ncol=length(r_F))
	alpha=matrix(seq(0,1,0.01),nrow=1)
	gamma1=matrix(gamma,nrow=1,ncol=ncol(alpha))
	gamma1[alpha==1]=0
	for (i in (1:length(r_F))){
		profit1=profit_disc_vouc(lambda,r_G0,r_F0,r_G,r_F[i],
			V_G,V_F,price_op,marg_cost,gamma1,alpha,beta)
		alpha_op[i]=alpha[which.max(profit1)]}
	if (algo==1){
		profit_sum=0
		alpha_op1=0
		for (j in unique(alpha_op)){
			alpha=matrix(j,1,1)
			if (alpha==1) gamma1=matrix(0,1,1)
			if (alpha!=1) gamma1=gamma			
			profit_sum1=sum(profit_disc_vouc(lambda,r_G0,r_F0,r_G,r_F,
				V_G,V_F,price_op,marg_cost,gamma1,alpha,beta))
			if (profit_sum1>profit_sum){
				profit_sum=profit_sum1
				alpha_op1=alpha}}
		alpha_op=alpha_op1}
	return(alpha_op)}

## firm optimal gamma (given a gamma)
## firm operation taking account of discount vouchers
## (this algorithm uses the brute force method to search for the gamma that 
## maximizes the profit)
# input:
# lambda - proportion of consumers from G population
# r_G0 - max. fractional increase in return probability of G population thanks
		# to issuance of discount vouchers
# r_F0 - ... F population ...
# r_G - return probability of G population
# r_F - ... of F population
# V_G - decay constant in exponential distribution of G population
# V_F - ... of F population
# price_op - firm product or service normal operation optimal price 
# marg_cost - firm product or service marginal cost
# alpha - (1 - discount rate) of the firm product or service price 
# beta - fraction of revenue received by firm using the voucher service
# algo - 1: constant optimal gamma over time
#	   2: optimal gamma that varies with time
# output:
# gamma_op - optimal gamma (fraction of consumers from G population taking up 
		# the discount vouchers) given the alpha 
firm_gamma_op1=function(lambda,r_G0,r_F0,r_G,r_F,V_G,V_F,price_op,marg_cost,
	alpha,beta,algo){
	gamma_op=matrix(0,nrow=1,ncol=length(r_F))
	gamma=matrix(seq(0,1,0.01),nrow=1)
	alpha1=matrix(alpha,nrow=1,ncol=ncol(gamma))
	alpha1[gamma==0]=1
	for (i in (1:length(r_F))){
		profit1=profit_disc_vouc(lambda,r_G0,r_F0,r_G,r_F[i],
			V_G,V_F,price_op,marg_cost,gamma,alpha1,beta)
		gamma_op[i]=gamma[which.max(profit1)]}
	if (algo==1){
		profit_sum=0
		gamma_op1=0
		for (j in unique(gamma_op)){
			gamma=matrix(j,1,1)
			if (gamma==0) alpha1=matrix(1,1,1)
			if (gamma!=0) alpha1=alpha			
			profit_sum1=sum(profit_disc_vouc(lambda,r_G0,r_F0,r_G,r_F,
				V_G,V_F,price_op,marg_cost,gamma,alpha1,beta))
			if (profit_sum1>profit_sum){
				profit_sum=profit_sum1
				gamma_op1=gamma}}
		gamma_op=gamma_op1}
	return(gamma_op)}

## firm optimal gamma (fix an alpha)
## firm operation taking account of discount vouchers
# input:
# lambda - proportion of consumers from G population
# r_G0 - max. fractional increase in return probability of G population thanks
		# to issuance of discount vouchers
# r_F0 - ... F population ...
# r_G - return probability of G population
# r_F - ... of F population
# V_G - decay constant in exponential distribution of G population
# V_F - ... of F population
# price_op - firm product or service normal operation optimal price 
# marg_cost - firm product or service marginal cost
# alpha - (1 - discount rate)
# beta - fraction of revenue received by firm using the voucher service
# algo - 1: constant optimal gamma over time
#	   2: optimal gamma that varies with time
# output:
# gamma_op - optimal gamma (fraction of consumers from G population taking up 
		# the discount vouchers) given the alpha 
firm_gamma_op=function(lambda,r_G0,r_F0,r_G,r_F,V_G,V_F,price_op,marg_cost,
	alpha,beta,algo){
	gamma_op=matrix(0,nrow=1,ncol=length(r_F))
	for (i in 1:length(r_F)){	
		r_F2=r_F[i]
		profit_deri_gamma=function(gamma,lambda1=lambda,r_G01=r_G0,
			r_F01=r_F0,r_G1=r_G,r_F1=r_F2,V_G1=V_G,V_F1=V_F,
			price_op1=price_op,marg_cost1=marg_cost,alpha1=alpha,
			beta1=beta){			
			alpha2=matrix(alpha1,ncol=length(gamma))
			alpha2[gamma==0]=1
			inc_r_G1=r_G01*(1-r_G1)*(1-exp(-(gamma+1-alpha2)))
			inc_r_F1=r_F01*(1-r_F1)*(1-exp(-(gamma+1-alpha2)))
			inc_r_G1_deri=r_G01*(1-r_G1)*exp(-(gamma+1-alpha2))
			inc_r_F1_deri=r_F01*(1-r_F1)*exp(-(gamma+1-alpha2))
			val_G1=price_op1/(r_G1+inc_r_G1)
			val_F1=price_op1/(r_F1+inc_r_F1)
			val_hat_G1=alpha2*price_op1/(r_G1+inc_r_G1)
			val_G1_deri=-price_op1/(r_G1+inc_r_G1)^2*inc_r_G1_deri
			val_F1_deri=-price_op1/(r_F1+inc_r_F1)^2*inc_r_F1_deri
			val_hat_G1_deri=-alpha2*price_op1/
				(r_G1+inc_r_G1)^2*inc_r_G1_deri
			out=lambda1*(r_G1+inc_r_G1)*(alpha2*beta1*price_op1-
				marg_cost1)*(1-exp_cdf(val_hat_G1,V_G1))+
				gamma*lambda1*inc_r_G1_deri*(alpha2*beta1*price_op1-
				marg_cost1)*(1-exp_cdf(val_hat_G1,V_G1))-gamma*
				lambda1*(r_G1+inc_r_G1)*(alpha2*beta1*price_op1-
				marg_cost1)*exp_pdf(val_hat_G1,V_G1)*val_hat_G1_deri-
				lambda1*(r_G1+inc_r_G1)*(price_op1-
				marg_cost1)*(1-exp_cdf(val_G1,V_G1))+(1-gamma)*
				lambda1*inc_r_G1_deri*(price_op1-marg_cost1)*
				(1-exp_cdf(val_G1,V_G1))-(1-gamma)*lambda1*(r_G1+
				inc_r_G1)*(price_op1-marg_cost1)*exp_pdf(val_G1,V_G1)*
				val_G1_deri+(1-lambda1)*inc_r_F1_deri*(price_op1-
				marg_cost1)*(1-exp_cdf(val_F1,V_F1))-(1-lambda1)*
				(r_F1+inc_r_F1)*(price_op1-marg_cost1)*
				exp_pdf(val_F1,V_F1)*val_F1_deri
			return(out)}
		gamma_op1<-tryCatch(uniroot.all(profit_deri_gamma,c(0,1),tol=0.001,maxiter=1000,n=100),error=function(e) 0)
		# print(gamma_op1)
		if ((gamma_op1[length(gamma_op1)]>=0) && 
			(gamma_op1[length(gamma_op1)]<=1) && length(gamma_op1)!=0) {
				gamma_op[i]<-gamma_op1[length(gamma_op1)]}}
	if (algo==1){
		profit_sum=0
		gamma_op1=0
		for (j in unique(gamma_op)){
			gamma=matrix(j,1,1)
			if (gamma==0) alpha1=matrix(1,1,1)
			if (gamma!=0) alpha1=alpha			
			profit_sum1=sum(profit_disc_vouc(lambda,r_G0,r_F0,r_G,r_F,
				V_G,V_F,price_op,marg_cost,gamma,alpha1,beta))
			if (profit_sum1>profit_sum){
				profit_sum=profit_sum1
				gamma_op1=gamma}}
		gamma_op=gamma_op1}
	return(gamma_op)}

## firm optimal alpha (fix a gamma, taking account of non-redeemed discount vouchers)
# profit_deri_gamma_omega=function(){}
## firm optimal alpha (fix a gamma, taking account of popularity index)
# profit_deri_gamma_omega_kappa=function(){}

## ---------------- firm set of optimal gamma(s) and alpha(s) ----------------

## firm optimal alpha(s) (fix a set of gamma(s))
## firm operation taking account of discount vouchers
# input:
# lambda - proportion of consumers from G population
# r_G0 - max. fractional increase in return probability of G population thanks
		# to issuance of discount vouchers
# r_F0 - ... F population ...
# r_G - return probability of G population
# r_F - ... of F population
# V_G - decay constant in exponential distribution of G population
# V_F - ... of F population
# price_op - firm product or service normal operation optimal price 
# marg_cost - firm product or service marginal cost
# gammaS - set of fractions of consumers from G population taking up the 
		# discount vouchers 
# beta - fraction of revenue received by firm using the voucher service
# algo - 1: set of constant optimal alphas over time
#	   2: set of optimal alphas that vary with time
# output:
# alphaS_op - optimal alpha(s) (1 - discount rate) given the set of gamma(s) 
firm_alphaS_op=function(lambda,r_G0,r_F0,r_G,r_F,V_G,V_F,price_op,marg_cost,
	gammaS,beta,algo){
	alphaS_op=matrix(1,nrow=length(gammaS),ncol=length(r_F))
	for (i in 1:length(r_F)){	
		r_F2=r_F[i]
		profit_deri_alphaS=function(alphaS,lambda1=lambda,r_G01=r_G0,
			r_F01=r_F0,r_G1=r_G,r_F1=r_F2,V_G1=V_G,V_F1=V_F,
			price_op1=price_op,marg_cost1=marg_cost,gamma1=gammaS,
			beta1=beta){			
			gamma2=matrix(gamma1,ncol=length(alphaS))
			gamma2[alphaS==1]=0
			inc_r_G1=r_G01*(1-r_G1)*(1-exp(-(sum(gamma2)+mean(1-alphaS))))
			inc_r_F1=r_F01*(1-r_F1)*(1-exp(-(sum(gamma2)+mean(1-alphaS))))
			inc_r_G1_deri=-r_G01*(1-r_G1)*exp(-(sum(gamma2)+
				mean(1-alphaS)))/length(gamma2[gamma2!=0])
			inc_r_F1_deri=-r_F01*(1-r_F1)*exp(-(sum(gamma2)+
				mean(1-alphaS)))/length(gamma2[gamma2!=0])
			val_G1=price_op1/(r_G1+inc_r_G1)
			val_F1=price_op1/(r_F1+inc_r_F1)
			val_G1_deri=-price_op1/(r_G1+inc_r_G1)^2*inc_r_G1_deri
			val_F1_deri=-price_op1/(r_F1+inc_r_F1)^2*inc_r_F1_deri
			out=matrix(0,nrow=length(gamma2))
			for (j in (1:length(gamma2))){
				for (k in (1:length(gamma2))){ 
					out[j]=out[j]+pi_alphaS_deri(j,k,alphaS,lambda1,
					r_G01,r_F01,r_G1,r_F1,V_G1,V_F1,price_op1,
					marg_cost1,gamma2,beta1,inc_r_G1,inc_r_G1_deri)}
				out[j]=out[j]+(1-sum(gamma2))*lambda1*inc_r_G1_deri*
				(price_op1-marg_cost1)*(1-exp_cdf(val_G1,V_G1))-(1-
				sum(gamma2))*lambda1*(r_G1+inc_r_G1)*(price_op1-
				marg_cost1)*exp_pdf(val_G1,V_G1)*val_G1_deri+(1-
				lambda1)*inc_r_F1_deri*(price_op1-marg_cost1)*(1-
				exp_cdf(val_F1,V_F1))-(1-lambda1)*(r_F1+inc_r_F1)*
				(price_op1-marg_cost1)*exp_pdf(val_F1,V_F1)*
				val_F1_deri}
			return(out)}
		alphaS_op1<-multiroot(profit_deri_alphaS,
				start=matrix(0.5,nrow=length(gammaS)),
				atol=0.001,maxiter=1000)$root
		# print(alphaS_op1)
		alphaS_op1[(alphaS_op1<0) | (alphaS_op1>1)]=1 
		alphaS_op[,i]<-alphaS_op1}

	if (algo==1){
		profit_sum=0
		alphaS_op1=matrix(0,nrow=nrow(alphaS_op),ncol=1)
		for (j in (1:ncol(alphaS_op))){
			alpha=matrix(alphaS_op[,j],ncol=1)
			gamma1=matrix(gammaS,ncol=1)
			gamma1[alpha==1]=0
			profit_sum1=sum(profit_disc_vouc(lambda,r_G0,r_F0,r_G,r_F,
				V_G,V_F,price_op,marg_cost,gamma1,alpha,beta))
			if (profit_sum1>profit_sum){
				profit_sum=profit_sum1
				alphaS_op1=alpha}}
		alphaS_op=alphaS_op1}
	return(alphaS_op)}

	pi_alphaS_deri=function(j3,k3,alphaS3,lambda3,r_G03,r_F03,r_G3,
			r_F3,V_G3,V_F3,price_op3,marg_cost3,gammaS3,beta3,
			inc_r_G3,inc_r_G3_deri){
		val_hat_G3=alphaS3[j3]*price_op3/(r_G3+inc_r_G3)
		val_hat_G3_deri1=price_op3/(r_G3+inc_r_G3)-alphaS3[j3]*price_op3/
			(r_G3+inc_r_G3)^2*inc_r_G3_deri
		val_hat_G3_deri2=-alphaS3[j3]*price_op3/(r_G3+inc_r_G3)^2*
			inc_r_G3_deri
		if (j3==k3){
			out=gammaS3[j3]*lambda3*(r_G3+inc_r_G3)*beta3*price_op3*
			(1-exp_cdf(val_hat_G3,V_G3))+gammaS3[j3]*lambda3*
			inc_r_G3_deri*(alphaS3[j3]*beta3*price_op3-marg_cost3)*
			(1-exp_cdf(val_hat_G3,V_G3))-gammaS3[j3]*lambda3*
			(r_G3+inc_r_G3)*(alphaS3[j3]*beta3*price_op3-marg_cost3)*
			exp_pdf(val_hat_G3,V_G3)*val_hat_G3_deri1}
		if (j3!=k3){
			out=gammaS3[j3]*lambda3*inc_r_G3_deri*(alphaS3[j3]*beta3*
			price_op3-marg_cost3)*(1-exp_cdf(val_hat_G3,V_G3))-
			gammaS3[j3]*lambda3*(r_G3+inc_r_G3)*(alphaS3[j3]*beta3*
			price_op3-marg_cost3)*exp_pdf(val_hat_G3,V_G3)*
			val_hat_G3_deri2}
	return(out)}

## firm optimal gamma(s) (fix a set of alpha(s))
## firm operation taking account of discount vouchers
# input:
# lambda - proportion of consumers from G population
# r_G0 - max. fractional increase in return probability of G population thanks
		# to issuance of discount vouchers
# r_F0 - ... F population ...
# r_G - return probability of G population
# r_F - ... of F population
# V_G - decay constant in exponential distribution of G population
# V_F - ... of F population
# price_op - firm product or service normal operation optimal price 
# marg_cost - firm product or service marginal cost
# alphaS - set of (1 - discount rates)
# beta - fraction of revenue received by firm using the voucher service
# algo - 1: set of constant optimal gammas over time
#	   2: set of optimal gammas that vary with time
# output:
# gammaS_op - optimal gamma(s) (fraction of consumers from G population taking 
		# up the discount vouchers) given the set of alpha(s) 
firm_gammaS_op=function(lambda,r_G0,r_F0,r_G,r_F,V_G,V_F,price_op,marg_cost,
	alphaS,beta,algo){
	gammaS_op=matrix(0,nrow=length(alphaS),ncol=length(r_F))
	for (i in 1:length(r_F)){	
		r_F2=r_F[i]
		profit_deri_gammaS=function(gammaS,lambda1=lambda,r_G01=r_G0,
			r_F01=r_F0,r_G1=r_G,r_F1=r_F2,V_G1=V_G,V_F1=V_F,
			price_op1=price_op,marg_cost1=marg_cost,alpha1=alphaS,
			beta1=beta){			
			alpha2=matrix(alpha1,ncol=length(gammaS))
			alpha2[gammaS[1:length(alpha2)]==0]=1
			inc_r_G1=r_G01*(1-r_G1)*(1-exp(-(sum(gammaS[1:length(alpha2)])+
					mean(1-alpha2))))
			inc_r_F1=r_F01*(1-r_F1)*(1-exp(-(sum(gammaS[1:length(alpha2)])+
					mean(1-alpha2))))
			inc_r_G1_deri=r_G01*(1-r_G1)*exp(-(sum(gammaS[1:length(alpha2)])+
					mean(1-alpha2)))
			inc_r_F1_deri=r_F01*(1-r_F1)*exp(-(sum(gammaS[1:length(alpha2)])+
					mean(1-alpha2)))
			val_G1=price_op1/(r_G1+inc_r_G1)
			val_F1=price_op1/(r_F1+inc_r_F1)
			val_G1_deri=-price_op1/(r_G1+inc_r_G1)^2*inc_r_G1_deri
			val_F1_deri=-price_op1/(r_F1+inc_r_F1)^2*inc_r_F1_deri
			out=matrix(0,nrow=length(alpha2))
			gamma0=1-sum(gammaS[1:length(alpha2)])
			for (j in (1:length(alpha2))){
				for (k in (1:length(alpha2))){ 
					out[j]=out[j]+pi_gammaS_deri(j,k,alpha2,lambda1,
					r_G01,r_F01,r_G1,r_F1,V_G1,V_F1,price_op1,
					marg_cost1,gammaS,beta1,inc_r_G1,inc_r_G1_deri)}
				out[j]=out[j]-lambda1*(r_G1+inc_r_G1)*(price_op1-
				marg_cost1)*(1-exp_cdf(val_G1,V_G1))+gamma0*lambda1*
				inc_r_G1_deri*(price_op1-marg_cost1)*
				(1-exp_cdf(val_G1,V_G1))-gamma0*lambda1*(r_G1+
				inc_r_G1)*(price_op1-marg_cost1)*exp_pdf(val_G1,V_G1)*
				val_G1_deri+(1-lambda1)*inc_r_F1_deri*(price_op1-
				marg_cost1)*(1-exp_cdf(val_F1,V_F1))-(1-lambda1)*
				(r_F1+inc_r_F1)*(price_op1-marg_cost1)*
				exp_pdf(val_F1,V_F1)*val_F1_deri}
			return(out)}
		gammaS_op1<-multiroot(profit_deri_gammaS,
				start=matrix(0.5,nrow=length(alphaS)),
				atol=0.001,maxiter=1000)$root
		# print(gammaS_op1)
		gammaS_op1[(gammaS_op1<0) | (gammaS_op1>1)]=0
		gammaS_op[,i]<-gammaS_op1}

	if (algo==1){
		profit_sum=0
		gammaS_op1=matrix(0,nrow=nrow(gammaS_op),ncol=1)
		for (j in (1:ncol(gammaS_op))){
			gamma=matrix(gammaS_op[,j],ncol=1)
			alpha1=matrix(alphaS,ncol=1)
			alpha1[gamma==0]=1			
			profit_sum1=sum(profit_disc_vouc(lambda,r_G0,r_F0,r_G,r_F,
				V_G,V_F,price_op,marg_cost,gamma,alpha1,beta))
			if (profit_sum1>profit_sum){
				profit_sum=profit_sum1
				gammaS_op1=gamma}}
		gammaS_op=gammaS_op1}
	return(gammaS_op)}

	pi_gammaS_deri=function(j3,k3,alphaS3,lambda3,r_G03,r_F03,r_G3,
			r_F3,V_G3,V_F3,price_op3,marg_cost3,gammaS3,beta3,
			inc_r_G3,inc_r_G3_deri){
		val_hat_G3=alphaS3[j3]*price_op3/(r_G3+inc_r_G3)
		val_hat_G3_deri=-alphaS3[j3]*price_op3/(r_G3+inc_r_G3)^2*
			inc_r_G3_deri
		if (j3==k3){
			out=lambda3*(r_G3+inc_r_G3)*(alphaS3[j3]*beta3*
			price_op3-marg_cost3)*(1-exp_cdf(val_hat_G3,V_G3))+
			gammaS3[j3]*lambda3*inc_r_G3_deri*(alphaS3[j3]*beta3*
			price_op3-marg_cost3)*(1-exp_cdf(val_hat_G3,V_G3))-
			gammaS3[j3]*lambda3*(r_G3+inc_r_G3)*(alphaS3[j3]*beta3*
			price_op3-marg_cost3)*exp_pdf(val_hat_G3,V_G3)*
			val_hat_G3_deri}
		if (j3!=k3){
			out=gammaS3[j3]*lambda3*inc_r_G3_deri*(alphaS3[j3]*beta3*
			price_op3-marg_cost3)*(1-exp_cdf(val_hat_G3,V_G3))-
			gammaS3[j3]*lambda3*(r_G3+inc_r_G3)*(alphaS3[j3]*beta3*
			price_op3-marg_cost3)*exp_pdf(val_hat_G3,V_G3)*
			val_hat_G3_deri}
	return(out)}
		
## firm optimal gamma(s) (fix a set of alpha(s))
## ... (taking account of non-redeemed discount vouchers)
## ... (taking account of popularity index)

## firm optimal alpha(s) (fix a set gamma(s))
## ... (taking account of non-redeemed discount vouchers)
## ... (taking account of popularity index)

## ------------------------ miscellaneous functions --------------------------

## population exponential distribution - cumulative distribution function(cdf)
# input:
# val - valuation of consumer
# V_popu - decay constant in exponential distribution of valuation
# output:
# cdf - cumulative density function
exp_cdf=function(v,V_popu){
	cdf=v
	cdf[v<0]=0
	cdf[v>=0]=1-exp(-v[v>=0]/V_popu)
return(cdf)}

## population exponential distribution - probability density function (pdf) 
# input:
# val - valuation of consumer
# V_popu - decay constant in exponential distribution of valuation
# output:
# prob - probability density function
exp_pdf=function(v,V_popu){
	prob=v
	prob[v<0]=0
	prob[v>=0]=exp(-v[v>=0]/V_popu)/V_popu
return(prob)}