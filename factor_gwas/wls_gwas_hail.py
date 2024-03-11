import hail as hl

#GWAS were run in batches
#m below specifies which factors to be run
m=range(30,35)

bgen_files = 'gs://fc-7d5088b4-7673-45b5-95c2-17ae00a04183/imputed/ukb_imp_chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}_v3.bgen'

pheno_file = 'gs://ukbb_prs/factor_gwas/new/FINAL_full/FS_newcode_full_semifinal.csv'
weight_file = 'gs://ukbb_prs/factor_gwas/new/FINAL_full/WLSweight_newcode_full_semifinal.csv'
cor_file = 'gs://ukbb_prs/factor_gwas/new/FINAL_full/MissCor_newcode_full_semifinal.csv'

covs_file = 'gs://ukbb_prs/factor_gwas/new/full_covs.csv'

ht_variants = hl.read_table('gs://ukb31063/ukb31063.neale_gwas_variants.ht')

phenos = hl.import_table(pheno_file, delimiter=',',missing=["","NA"],impute=True,types={'s': hl.tstr}).key_by('s')
weights = hl.import_table(weight_file, delimiter=',',missing=["","NA"],impute=True,types={'s': hl.tstr}).key_by('s')
cors = hl.import_table(cor_file, delimiter=',',missing=["","NA"],impute=True,types={'s': hl.tstr}).key_by('s')
covs = hl.import_table(covs_file, delimiter=',',missing=["","NA"],impute=True,types={'userId': hl.tstr}).key_by('userId')

threshold=hl.eval(hl.sqrt(0.80))

df_phenos_and_cors = phenos.annotate(cors = cors[phenos.s])
to_select = {i : hl.if_else(df_phenos_and_cors.cors[i] < threshold, hl.missing(hl.tfloat64), df_phenos_and_cors[i]) for i in list(phenos.row_value)}
new_df_phenos = df_phenos_and_cors.select(**to_select)

pheno_names = [list(phenos.row_value)[x] for x in m]

mt = hl.import_bgen(bgen_files,entry_fields=["dosage"],sample_file="gs://ukb31063/ukb31063.autosomes.sample",variants=ht_variants)

mycovs = list(covs.row_value)[0:21] + list(covs.row_value)[27:]

mt = mt.annotate_cols(phenotypes=new_df_phenos[mt.s],covariates=covs[mt.s],weights=weights[mt.s])

mt = mt.filter_cols(hl.is_defined(mt.phenotypes))

ys = [mt['phenotypes'][x] for x in list(mt['phenotypes'].keys())]
wt = [mt['weights'][x] for x in list(mt['weights'].keys())]
covs = [mt['covariates'][x] for x in mycovs]

ys=[ys[x] for x in m]
wt=[wt[x] for x in m]

new_y = [[pheno] for pheno in ys]

linreg_results = hl.linear_regression_rows(y=new_y, x=mt.dosage, covariates = [1, *covs], weights=wt, block_size=4)

results = linreg_results.checkpoint('gs://ukbb_prs/factor_gwas/new/FINAL_full/results/results_gwas_FINAL_thresh_r2_80_asscovs7.mt')

results = results.annotate(a1 = results.alleles[0], a2 = results.alleles[1])

for i in range(len(m)):
  fields = {'n': results.n[i]}
  for field in ['beta', 'standard_error', 't_stat', 'p_value']:
    fields[field] = results[field][i][0]
  results.select(a1 = results.alleles[0], a2 = results.alleles[1], **fields).export(f'gs://ukbb_prs/factor_gwas/new/FINAL_full/results/{pheno_names[i]}.tsv')
