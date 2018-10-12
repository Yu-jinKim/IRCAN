library(ggplot2)
library(data.table)

setwd("~/Desktop/Yujin/Rstudio")

dat=read.table("captureRegions",sep="\t")

colnames(dat)=c("chr", "CS", "CE", "CloneName", "score", "strand", "# locs per clone", "# capReg alignments", "Type of region",
								 "Rep1: # seqs per clone","Rep1: # seqs at location",
								 "Rep2: # seqs per clone","Rep2: # seqs at location",
								 "Input: # seqs per clone","Input: # seqs at location",
								 "Lib: # seqs per clone","Lib: # seqs at location",
								 "CloseReg ratio")

Rep1=28716324
Rep2=22606756
Inp=39562439
Lib=21919204

# get values according to the values of another column
subset_Rep1=as.data.frame(sapply(1:150, function(x){subset(dat, (dat$`Rep1: # seqs per clone`==x) & (!is.na(dat$`Rep1: # seqs per clone`)))}))
subset_Rep2=as.data.frame(sapply(1:150, function(x){subset(dat, (dat$`Rep2: # seqs per clone`==x) & (!is.na(dat$`Rep2: # seqs per clone`)))}))
subset_Inp=as.data.frame(sapply(1:150, function(x){subset(dat, (dat$`Input: # seqs per clone`==x) & (!is.na(dat$`Input: # seqs per clone`)))}))
subset_Lib=as.data.frame(sapply(1:150, function(x){subset(dat, (dat$`Lib: # seqs per clone`==x) & (!is.na(dat$`Lib: # seqs per clone`)))}))

DF=subset_Rep1
column="Input....seqs.per.clone"
ratio=Rep1/Inp
DF_bin=subset_Rep1
column_bin="Rep1....seqs.per.clone"

DF=subset_Rep2
column="Lib....seqs.per.clone"
ratio=Rep2/Lib
DF_bin=subset_Rep2
column_bin="Rep2....seqs.per.clone"

d1=sapply(1:150, function(x){
	tmp=as.data.frame(DF[,x])
	tmp=tmp[!duplicated(tmp[,c("CloneName","strand")]),]
	tmp=tmp[[column]]*ratio
})

d1_bin=sapply(1:150, function(x){
	tmp=as.data.frame(DF_bin[,x])
	tmp=tmp[!duplicated(tmp[,c("CloneName","strand")]),]
	tmp=tmp[[column_bin]]
})

# Compute distribution of data and simulation

k=20
new_data=d1[[k]]
pois=rpois(1000, lambda = k)
sim_pois=rpois(1000, lambda = mean(new_data, na.rm=TRUE))
size=mean(new_data, na.rm=TRUE)^2/(var(new_data, na.rm=TRUE)-mean(new_data, na.rm=TRUE))
sim_nb=rnbinom(1000, size=size, mu=mean(new_data, na.rm=TRUE))

DF_d1=data.frame(values=new_data)
DF_pois=data.frame(values=pois)
DF_sim=data.frame(values=sim_pois)
DF_nb=data.frame(values=sim_nb)

DF_d1$dist="d1"
DF_pois$dist="DF_pois"
DF_sim$dist="DF_sim"
DF_nb$dist="DF_nb"

distLength=rbind(DF_pois,DF_d1,DF_sim,DF_nb)

ggplot(distLength, aes(values, fill = dist)) + geom_density(alpha = 0.5)

# Calculate mean and var of all input values

d_means=sapply(1:150, function(x){mean(d1[[x]], na.rm = TRUE)})
d_var=sapply(1:150, function(x){var(d1[[x]], na.rm = TRUE)})

varMeans=d_var/d_means

poissons=c()

for(i in 1:150) {
	loi_poisson=rpois(sort(d1[[i]]), lambda=mean(d1[[i]]))
	poissons=c(poissons,c(var(loi_poisson, na.rm=TRUE)/mean(loi_poisson, na.rm=TRUE)))
}

plot(1:150, varMeans, col="red", ylim=c(0,10))        # var/mean vs lib
points(1:150, poissons, col="blue")                   # poisson vs lib
plot(1:150, d_var, ylim = c(0,50))                    # var vs lib
plot(1:150, d_means)
plot(d_var, d_means)                                  # var vs mean
plot(varMeans, d_means)                               # var/mean vs mean
plot(range_var)
plot(varmean, ylim=c(0,10), ylab = "Var Rep2/Mean Lib")
title("Range of var and mean")
lines(range_mean_bin,col="blue")
lines(range_mean/range_mean_bin,col="red")

# group range of adjacents counts for better variance

div=1000
flat_d1=unlist(d1)
flat_d1_bin=unlist(d1_bin)
steps=length(flat_d1)/div
steps_bin=length(flat_d1_bin)/div

range_var=sapply(1:steps, function(x){
	floor=(x-1)*div+1
	ceil=x*div
	var(flat_d1[floor:ceil],na.rm=TRUE)
})

range_mean=sapply(1:steps, function(x){
	floor=(x-1)*div+1
	ceil=x*div
	mean(flat_d1[floor:ceil], na.rm=TRUE)
})

range_mean_bin=sapply(1:steps_bin, function(x){
	floor=(x-1)*div+1
	ceil=x*div
	mean(flat_d1_bin[floor:ceil])
})

varmean=range_var/range_mean_bin
mean_var=mean(varmean)
meangt1_var=mean(varmean[range_mean_bin>1])

# computing the chosen variance
cutoff_var=mean(d_var[1:20], na.rm = TRUE)
cutoff_mean=mean(d_means[1:20], na.rm = TRUE)

ranged_var=mean(range_var, na.rm = TRUE)

range_computed_var=size(cutoff_mean,ranged_var)
cutoff_computed_var=size(cutoff_mean,cutoff_var)

# functions
size=function(mean, var){
	mean^2/(var-mean)
}

vnbinom = function(x, rep, lib, scale_factor){
	if(!is.na(x[rep]) && !is.na(x[lib])){
		dnbinom(as.numeric(x[rep]),mu=as.numeric(x[lib]),size=size(as.numeric(x[lib]),as.numeric(x[lib])*scale_factor))
	}
	else{
		tmp=NA
	}
}

negBinom = function(rep, ctrl, activity){
	if(!is.na(rep) && !is.na(ctrl) && !is.na(activity)){
		value=as.numeric(rep)
		mu=as.numeric(ctrl)/activity
		var=mean_var*mu
		size=mu^2/(var-mu)
		dnbinom(value,mu=mu,size=size)
	}
}

# replace NA in Input by 0
dat[,14][is.na(dat[,14])] = 0
# calculate true abundance for Input == 0
trueAbundance_mean0 = mean(dat[,16][dat[,14] == 0], na.rm = TRUE)
# replace 0 by calculated true abundance
dat[,14][dat[,14] == 0] = trueAbundance_mean0

# range of true abundance
valuesInput=split(dat[,c(-4)], dat$`Input: # seqs per clone`)
valuesInput=lapply(1:length(valuesInput), function(x){droplevels(valuesInput[[x]])})
# calculate true abundance
trueAbundance_mean=sapply(1:length(valuesInput), function(x){
	mean(as.numeric(valuesInput[[x]][["Lib: # seqs per clone"]]), na.rm = TRUE)
})

# split the data
splitting=function(dat, column){
	tmp=split(dat[,c(-4)], column, drop=TRUE)
	tmp=lapply(1:length(tmp), function(x){droplevels(tmp[[x]])})
}

res=splitting(dat, dat$`Type of region`)
res_test=splitting(dat, dat$`Type of region`)

res_test=split(dat[,c(-4,-9)], dat$`Type of region`, drop=TRUE)

# label capture region with biggest number reads per clone
# replace by true abundance
# compute activity
# compute p-values
# likelihood/chi-squared likelihood ratio
# remove capture regions which don't have less than n fragments

# main function
process=function(l, cloneThreshold = 20, trueAbundance_means, use_trueAbundance = FALSE, removeReg = FALSE){
	if(removeReg == TRUE){
		tmp=list()	
	}
	
	for(i in 1:length(l)){
		
		if(use_trueAbundance == TRUE){
			l[[i]][["Input: # seqs per clone"]]=trueAbundance_means[l[[i]][["Input: # seqs per clone"]]+1]
		}

		l[[i]][["Max value"]]=sapply(l[i], function(x){max(c(max(x[["Rep1: # seqs per clone"]], na.rm = TRUE),max(x[["Rep2: # seqs per clone"]], na.rm = TRUE),max(x[["Input: # seqs per clone"]], na.rm = TRUE),max(x[["Lib: # seqs per clone"]], na.rm = TRUE)))})
		
		l[[i]][["capReg_activity"]]=(sum(l[[i]]["Lib: # seqs per clone"], na.rm = TRUE)+sum(l[[i]]["Input: # seqs per clone"], na.rm = TRUE))/(sum(l[[i]]["Rep2: # seqs per clone"], na.rm = TRUE)+sum(l[[i]]["Rep1: # seqs per clone"], na.rm = TRUE))
		
		for(j in 1:nrow(l[[i]])){
			if(!is.na(l[[i]][["Rep1: # seqs per clone"]][[j]]) && !is.na(l[[i]][["Input: # seqs per clone"]][[j]]) && !is.na(l[[i]][["capReg_activity"]][[j]])){
				l[[i]][["pval_obs1"]][[j]]=negBinom(l[[i]][["Rep1: # seqs per clone"]][[j]], l[[i]][["Input: # seqs per clone"]][[j]], l[[i]][["capReg_activity"]][[j]])
			}
			else{l[[i]][["pval_obs1"]][[j]]=NA
			}
			if(!is.na(l[[i]][["Rep2: # seqs per clone"]][[j]]) && !is.na(l[[i]][["Lib: # seqs per clone"]][[j]]) && !is.na(l[[i]][["capReg_activity"]][[j]])){
				l[[i]][["pval_obs2"]][[j]]=negBinom(l[[i]][["Rep2: # seqs per clone"]][[j]], l[[i]][["Lib: # seqs per clone"]][[j]], l[[i]][["capReg_activity"]][[j]])
			}
			else{l[[i]][["pval_obs2"]][[j]]=NA
			}
			if(!is.na(l[[i]][["Rep1: # seqs per clone"]][[j]]) && !is.na(l[[i]][["Input: # seqs per clone"]][[j]])){
				l[[i]][["pval_exp1"]][[j]]=negBinom(l[[i]][["Rep1: # seqs per clone"]][[j]], l[[i]][["Input: # seqs per clone"]][[j]], Inp/Rep1)
			}
			else{l[[i]][["pval_exp1"]][[j]]=NA
			}
			if(!is.na(l[[i]][["Rep2: # seqs per clone"]][[j]]) && !is.na(l[[i]][["Lib: # seqs per clone"]][[j]])){
				l[[i]][["pval_exp2"]][[j]]=negBinom(l[[i]][["Rep2: # seqs per clone"]][[j]], l[[i]][["Lib: # seqs per clone"]][[j]], Lib/Rep2)
			}
			else{l[[i]][["pval_exp2"]][[j]]=NA
			}
		}
		
		l[[i]][["likelihood_obs"]]=sapply(l[i], function(x){prod(x[["pval_obs1"]], na.rm=TRUE)*prod(x[["pval_obs2"]], na.rm=TRUE)})
		l[[i]][["likelihood_exp"]]=sapply(l[i], function(x){prod(x[["pval_exp1"]], na.rm=TRUE)*prod(x[["pval_exp2"]], na.rm=TRUE)})
		l[[i]][["chi2"]]=sapply(l[i], function(x){
			if(mean(x[["likelihood_exp"]])<mean(x[["likelihood_obs"]])){
				pchisq(-2*log(mean(x[["likelihood_exp"]])/mean(x[["likelihood_obs"]])),1,lower.tail = FALSE)
			}
			else{
				pchisq(-2*log(mean(x[["likelihood_exp"]])/mean(x[["likelihood_obs"]])),1,lower.tail = TRUE)
			}
		})
		
		if(removeReg == TRUE){
			if(nrow(l[[i]])>=cloneThreshold){
				tmp[[i]]=as.data.frame(l[[i]])
			}
		}
	}
	
	if(removeReg == TRUE){
		return(tmp)
	}
	else{
		return(l)
	}
}

res_processed=process(l = res, cloneThreshold = 20, trueAbundance_means = trueAbundance_mean, use_trueAbundance = FALSE, removeReg = FALSE)

res_test_processed=process(l = res_test, trueAbundance_means = trueAbundance_mean, use_trueAbundance = FALSE, removeReg = TRUE)

# plot in list of dataframes
chi2=sapply(res_processed, function(x){mean(x[["chi2"]], na.rm = TRUE)})
likelihood=sapply(res_processed, function(x){mean(x[["capReg_activity"]], na.rm = TRUE)})

tab=rbind(log(likelihood), -2*log(chi2))
plot(tab[1,], tab[2,], xlab = "Log(likelihood)", ylab = "-2*log(chi2)")
title("remove and true abundance")


chi2_test=sapply(res_test_processed, function(x){mean(x[["chi2"]], na.rm = TRUE)})
likelihood_test=sapply(res_test_processed, function(x){mean(x[["capReg_activity"]], na.rm = TRUE)})

tab_test=rbind(log(likelihood_test), -2*log(chi2_test))
plot(tab_test[1,], tab_test[2,], xlab = "Log(likelihood)", ylab = "-2*log(chi2)")
title("remove regions")

# ggplot
tab[tab==Inf]=NA
tab=as.data.frame(t(tab))
ggplot(data=tab, aes(x=tab[,1], y=tab[,2])) + geom_jitter()

# rbind into dataframe
res_df=rbindlist(res_processed)
res_df$`logLike`=log(res_df$capReg_activity)
res_df$`-logChi2`=-2*log(res_df$chi2)

res_test_df=rbindlist(res_test_processed)
res_test_df$`logLike`=log(res_test_df$capReg_activity)
res_test_df$`-logChi2`=-2*log(res_test_df$chi2)

# plot
plot(res_df$logLike, res_df$`-logChi2`, ylab="-2*log(chi2 pvalue)", xlab="log(Obs activity)")
plot(res_test_df$logLike, res_test_df$`-logChi2`, ylab="-2*log(chi2 pvalue)", xlab="log(Obs activity)")
title(main="True abundance for 1:100")

plot(remove_df$logLike, remove_df$`-logChi2`, ylab="-2*log(chi2 pvalue)", xlab="log(Obs activity)")
