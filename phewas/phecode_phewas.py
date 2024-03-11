import pandas as pd
import numpy as np
import statsmodels
import statsmodels.api as sm
import statsmodels.formula.api as smf
import sys

myfac = sys.argv[1]

fs = pd.read_csv('/stanley/robinson/ccarey/UKBB/factor_gwas/factor_analysis/cfa/FINAL_CFA/FACTOR_SCORES/core/newcode/semifinalfullind/FS_newcode_full_semifinal_ind.csv').set_index("s")

misscor = pd.read_csv('/stanley/robinson/ccarey/UKBB/factor_gwas/factor_analysis/cfa/FINAL_CFA/FACTOR_SCORES/core/newcode/semifinalfullind/MissCor_newcode_full_semifinal_ind.csv').set_index("s")

phenosum = pd.read_table('/stanley/robinson/ccarey/UKBB/factor_gwas/factor_analysis/cfa/FINAL_CFA/FACTOR_SCORES/PHEWAS/round3_phenos_summary.tsv')

phenodict_phecodes = {}
for pheno in phenosum[phenosum.trait_type=="phecode"].phenocode:
    phenodict_phecodes[pheno] = phenosum.loc[phenosum.phenocode==pheno,"description"].values[0]

misscor_ind = pd.read_csv('/stanley/robinson/ccarey/UKBB/factor_gwas/factor_analysis/cfa/FINAL_CFA/FACTOR_SCORES/core/wls_gwas/FINAL/MissCor_newcode_full_semifinal.csv').set_index("s")

thresh = np.sqrt(0.80)

fs_thresh = fs.copy()
fs_thresh[misscor<thresh]=np.nan

fs_thresh[misscor_ind<thresh]=np.nan

fs_thresh.f27[misscor_ind.f27<np.sqrt(0.90)]=np.nan

fs_thresh[misscor_ind.isna()] = np.nan

phecodes = pd.read_table('/stanley/robinson/ccarey/UKBB/factor_gwas/factor_analysis/cfa/FINAL_CFA/FACTOR_SCORES/PHEWAS/round3_phenos_phecodes.tsv').set_index("userId")
phecodes = phecodes.loc[fs.index]

covs = pd.read_table("/stanley/robinson/ccarey/UKBB/factor_gwas/factor_analysis/cfa/FINAL_CFA/FACTOR_SCORES/core/wls_gwas/to_gwas/full_covs.csv",sep=",",index_col="userId")
covs.isFemale = covs.isFemale.astype(int)

mycovs = covs.columns[0:21].tolist() + covs.columns[27:].tolist()

fs_thresh["f4"] = fs_thresh.f4*-1
fs_thresh["f10"] = fs_thresh.f10*-1
fs_thresh["f14"] = fs_thresh.f14*-1
fs_thresh["f19"] = fs_thresh.f19*-1
fs_thresh["f23"] = fs_thresh.f23*-1
fs_thresh["f24"] = fs_thresh.f24*-1
fs_thresh["f36"] = fs_thresh.f36*-1

input_data = fs_thresh.join(phecodes,how="left").join(covs,how="left")

phecodes_use = phecodes.columns[input_data[phecodes.columns].sum()>=250]

phecodes_use_female = phecodes_use[input_data.groupby("isFemale")[phecodes_use].sum().loc[0]>=25]
phecodes_use_use = phecodes_use_female[input_data.groupby("isFemale")[phecodes_use_female].sum().loc[1]>=25]

len(phecodes_use_use)

for factor in [fs.columns[int(myfac)-1]]:
    temp = pd.DataFrame(index=phecodes_use_use, columns=["name","n","n_case","coef","se","teststat","or","ci_lower","ci_upper","p","has_nan"])
    for i,item in enumerate(phecodes_use_use):
        y=input_data[item]
        x=input_data[[factor]+mycovs]
        x = x[x.columns[x.dropna().var()!=0]]
        x=sm.add_constant(x)
        model = sm.GLM(y, x, family=sm.families.Binomial(), missing="drop").fit(cov_type='HC0')
        temp.loc[item,"name"] = phenodict_phecodes[item]
        temp.loc[item,"n"] = model.nobs
        temp.loc[item,"n_case"] = int(y.dropna().sum())
        temp.loc[item,"coef"] = model.params[1]
        temp.loc[item,"se"] = model.bse[1]
        temp.loc[item,"teststat"] = model.tvalues[1]
        temp.loc[item,"or"] = np.exp(model.params[1])
        temp.loc[item,["ci_lower","ci_upper"]] = np.exp(model.conf_int().loc[factor]).values
        temp.loc[item,"p"] = model.pvalues[1]
        temp.loc[item,"has_nan"] = sum(np.isnan(model.bse))>0
        temp.to_csv('/stanley/robinson/ccarey/UKBB/factor_gwas/factor_analysis/cfa/FINAL_CFA/FACTOR_SCORES/PHEWAS/phecodes_gwassamp/'+factor+"_phecodes_phewas.csv")
        print(factor,i,item,phenodict_phecodes[item],model.nobs,int(y.dropna().sum()),round(model.params[1],3),round(model.tvalues[1],3),model.pvalues[1],sum(np.isnan(model.bse))>0)
