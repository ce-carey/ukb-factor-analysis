#!/usr/bin/Rscript

install.packages("lavaan.tar.gz",repos=NULL,type="source")
install.packages("dplyr")
install.packages("truncnorm")
library(dplyr)
library(truncnorm)
library(lavaan)

core = read.csv("full_data.csv")

samples = core$userId

fitobj = readRDS("fit.Rds")

#varnames 
ov_names = unlist(fitobj@pta$vnames$ov)
ov_exo_names = unlist(fitobj@pta$vnames$ov.x)
ov_endog_names = unlist(fitobj@pta$vnames$ov.nox)
ov_numeric_names = unlist(fitobj@pta$vnames$ov.num)
ov_numeric_noexo_names = setdiff(ov_numeric_names,ov_exo_names)
ov_ordered_names = unlist(fitobj@pta$vnames$ov.ord)
ov_names_noexo = unlist(fitobj@pta$vnames$ov.nox)
factor_names = unlist(fitobj@pta$vnames$lv.regular)

core = core[,ov_names]

core[,ov_ordered_names] = lapply(core[,ov_ordered_names],ordered)

#standardize continuous vars (and covariates)
core[,ov_numeric_names] = lapply(core[,ov_numeric_names],as.numeric)
core[,ov_numeric_names] = scale(core[,ov_numeric_names,drop=FALSE])

core_resid = as.data.frame(matrix(nrow=length(rownames(core)),ncol=length(colnames(core))))
rownames(core_resid) = rownames(core)
colnames(core_resid) = colnames(core)
core_resid[,ov_exo_names] = core[,ov_exo_names]

#residualize continuous vars
#intercepts all set to 0? check with raymond...
gamma_mat = lavInspect(fitobj,what="est")$gamma
coefs = gamma_mat[ov_numeric_noexo_names,]
exo_vals = core[,ov_exo_names]
num_y = core[,ov_numeric_noexo_names]
num_yhat = as.matrix(exo_vals) %*% t(coefs)
core_resid[,ov_numeric_noexo_names] = num_y - num_yhat

ov_info = fitobj@Data@ov
ordered.idx = which(fitobj@Data@ov$type=="ordered")

th = lavInspect(fitobj, what="th")

#residualize ordered vars
for(i in ordered.idx) {
        name = ov_info$name[i]
	print(name)
        nth = ov_info$nlev[i]-1
        lnam = unlist(strsplit(ov_info$lnam[i],"\\|"))
        thresh_prefix = paste0(name,"|t")
        zeta = rep(NA,nth)
        for(level in 1:nth){
                zeta[level] = th[paste0(thresh_prefix,level)]
        }

        coefs = gamma_mat[name,]
        X = core[,ov_exo_names]
        fitted_vals = as.matrix(X) %*% coefs
	response = unname(as.integer(core[,name]))

	bounds = unname(c(-Inf,zeta-zeta[1L],Inf))
	mean_responses = fitted_vals - zeta[1L]
	
	resid = etruncnorm(a=bounds[response],b=bounds[response+1L],mean=mean_responses,sd=1)-mean_responses

	core_resid[, name] = resid        
}


loadings = lavInspect(fitobj, what="est")$beta[ov_names_noexo,factor_names]

residdiag = diag(lavInspect(fitobj, what="est")$psi)[ov_names_noexo]

varcov = lavInspect(fitobj,what="sampstat")$res.cov

sigma_hat = lavInspect(fitobj,what="implied")$res.cov

#correct resids for ordered variables 
for(i in ordered.idx) {
        name = ov_info$name[i]
        print(name)

	sd_obs = sd(core_resid[,name],na.rm=TRUE)	

	var_obs = var(core_resid[,name],na.rm=TRUE)

        nth = ov_info$nlev[i]
        lnam = unlist(strsplit(ov_info$lnam[i],"\\|"))
        thresh_prefix = paste0(name,"|t")
        zeta = rep(NA,nth-1)
        for(level in 1:(nth-1)){
                zeta[level] = th[paste0(thresh_prefix,level)]
        }
	
	response = unname(as.integer(core[,name]))
	
	residvar = residdiag[name]
	
	x_i = sort(unique(response[!is.na(response)]))
	p_i = c(pnorm(zeta),1) - c(0,pnorm(zeta))

	mu = sum(p_i*x_i)
	sigma2_y = sum(p_i * (x_i - mu)**2)
	
	adjresidvar = var_obs*(1-((1-residvar)*(sum(dnorm(zeta))^2/sigma2_y)))
	
	residdiag[name] = adjresidvar

	loadings_adj = sd_obs*(sum(dnorm(zeta))/sqrt(sigma2_y))

	loadings[name,] = loadings[name,]*loadings_adj

}

resid = diag(residdiag)

#go factor by factor
factor_mat = lavInspect(fitobj,what="free")$beta[-(1:35),1:35]>0

#construct output data frames
FS = as.data.frame(matrix(nrow=length(rownames(core)),ncol=length(colnames(factor_mat))))
regression_weights = as.data.frame(matrix(nrow=length(rownames(core)),ncol=length(colnames(factor_mat))))
miss_nomiss_cor = as.data.frame(matrix(nrow=length(rownames(core)),ncol=length(colnames(factor_mat))))

rownames(FS) = rownames(core)
rownames(regression_weights)=rownames(core)
rownames(miss_nomiss_cor) = rownames(core)
colnames(FS) = colnames(factor_mat)
colnames(regression_weights)=colnames(factor_mat)
colnames(miss_nomiss_cor) = colnames(factor_mat)

items = rownames(factor_mat)
mp = lavaan:::lav_data_missing_patterns(core_resid[items])
#mp = readRDS("mp.Rds")

n_factors = length(colnames(factor_mat))

coresamp_cov = cov(core_resid[,items],use="pairwise")

fullresid = diag(residdiag)

W_nomiss = solve(fullresid) %*% loadings %*% solve( diag(n_factors) + t(loadings) %*% solve(fullresid) %*% loadings )

for(i in 1:length(mp$case.idx)){
	print(paste0(i, " of ",length(mp$case.idx)))
	iids = as.numeric(unlist(mp$case.idx[i]))
        mp_data = core_resid[iids,items]
	n_samples = length(iids)
	misspat = mp$pat[i,]
	lambda = loadings[misspat,]
	W_nomiss_temp = W_nomiss
	mp_resid = diag(residdiag[misspat])
	emptyfact = which(colSums(lambda!=0)==0)

	if(length(emptyfact)>0){
		lambda = lambda[,-emptyfact]
		W_nomiss_temp = W_nomiss_temp[,-emptyfact]
	}

	W = solve(mp_resid) %*% lambda %*% solve( diag(length(colnames(lambda))) + t(lambda) %*% solve(mp_resid) %*% lambda )
	rownames(W) = rownames(lambda)
	omegas = t(W) %*% as.matrix(coresamp_cov[misspat,misspat]) %*% W
	omega_S_S = diag(omegas)
	f_S = as.matrix(mp_data[,misspat]) %*% as.matrix(W)
	inv_var_weight = 1/(omega_S_S)

        W_miss = matrix(0, length(items),length(colnames(lambda)))
        W_miss[misspat] = W
        miss_nomiss = diag( t(W_miss) %*% coresamp_cov %*% W_nomiss_temp) / (sqrt(diag( t(W_miss) %*% coresamp_cov %*% W_miss))*sqrt(diag( t(W_nomiss_temp) %*% coresamp_cov %*% W_nomiss_temp)))
	FS[iids,colnames(lambda)] = f_S
	regression_weights[iids,colnames(lambda)]=matrix(rep(as.vector(unlist(inv_var_weight)),n_samples),ncol=length(colnames(lambda)),byrow=TRUE)
	miss_nomiss_cor[iids,colnames(lambda)]=matrix(rep(as.vector(unlist(miss_nomiss)),n_samples),ncol=length(colnames(lambda)),byrow=TRUE)
}

rownames(FS) = samples
rownames(regression_weights) = rownames(FS)
rownames(miss_nomiss_cor) = rownames(FS)

colnames(FS) = colnames(factor_mat)
colnames(regression_weights) = colnames(factor_mat)
colnames(miss_nomiss_cor) = colnames(factor_mat)


saveRDS(FS, "FS_newcode_full_semifinal_ind.Rds")
write.csv(FS,"FS_newcode_full_semifinal_ind.csv",row.names=TRUE,quote=FALSE)
saveRDS(regression_weights,"WLSweight_newcode_full_semifinal_ind.Rds")
write.csv(regression_weights,"WLSweight_newcode_full_semifinal_ind.csv",row.names=TRUE,quote=FALSE)
saveRDS(miss_nomiss_cor,"MissCor_newcode_full_semifinal_ind.Rds")
write.csv(miss_nomiss_cor,"MissCor_newcode_full_semifinal_ind.csv",row.names=TRUE,quote=FALSE)

