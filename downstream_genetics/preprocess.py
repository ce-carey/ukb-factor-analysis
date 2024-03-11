import pandas as pd
import numpy as np
import glob

fullsamp_wt = pd.read_csv('/stanley/robinson/ccarey/UKBB/factor_gwas/factor_analysis/cfa/FINAL_CFA/FACTOR_SCORES/core/wls_gwas/FINAL/WLSweight_newcode_full_semifinal.csv').set_index("s")

misscor = pd.read_csv('/stanley/robinson/ccarey/UKBB/factor_gwas/factor_analysis/cfa/FINAL_CFA/FACTOR_SCORES/core/wls_gwas/FINAL/MissCor_newcode_full_semifinal.csv').set_index("s")

threshold = np.sqrt(0.8)

fullsamp_wt[misscor<threshold] = np.nan

miss0wt = pd.read_table('/stanley/robinson/ccarey/UKBB/factor_gwas/factor_analysis/cfa/FINAL_CFA/FACTOR_SCORES/core/newcode/maxweights_full.txt',sep="\s+")

miss0wtdict = dict(zip(miss0wt.factor,miss0wt.max_weight))

neffdict_0missnotrunc = {}
for factor in  fullsamp_wt.columns:
    neffdict_0missnotrunc[factor] = round(sum(fullsamp_wt[factor].dropna()/miss0wtdict[factor]))

infiles = glob.glob("/stanley/robinson/ccarey/UKBB/factor_gwas/factor_analysis/cfa/FINAL_CFA/FACTOR_SCORES/core/wls_gwas/FINAL/outputs/*.tsv")

print("loading var_annot")
var_annot = pd.read_csv('/stanley/robinson/ccarey/UKBB/factor_gwas/factor_analysis/cfa/FINAL_CFA/FACTOR_SCORES/core/test_wls_gwas/variants.tsv',sep="\t")

for file in infiles:
    temp = pd.read_csv(file,sep="\t")
    temp_plusannot = pd.concat([temp,var_annot],axis=1)
    temp_plusannot_maf = temp_plusannot[temp_plusannot.minor_AF>0.001].dropna(subset=["locus","p_value"])
    temp_plusannot_maf["A1"] = temp_plusannot_maf["a2"]
    temp_plusannot_maf["A2"] = temp_plusannot_maf["a1"]
    temp_plusannot_maf["N"] = neffdict_0missnotrunc[file.split("/")[-1].split(".")[0]]
    temp_to_output = temp_plusannot_maf[["rsid","A1","A2","N","p_value","minor_AF","info","t_stat"]]
    temp_to_output.columns = ["SNP","A1","A2","N","P","MAF","INFO","Z"]
    temp_to_output.to_csv("/stanley/robinson/ccarey/UKBB/factor_gwas/factor_analysis/cfa/FINAL_CFA/FACTOR_SCORES/core/wls_gwas/FINAL/outputs/preprocessed/"+file.split("/")[-1],sep="\t",index=False)
