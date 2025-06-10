import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats
import math
from statsmodels.stats.multitest import multipletests
import sys


###########################################################################################################################################
OUT = "/PROJECTES/MALARIA_EPIGENETICS/Projects/ReCHIP_GDV1_Project_ene25/New_DiffHC/transgenicVSwt/M-T6vsNF54/M-T6vsNF54_results.csv"
fold_change_data = pd.read_csv('/PROJECTES/MALARIA_EPIGENETICS/Projects/ReCHIP_GDV1_Project_ene25/New_DiffHC/transgenicVSwt/M-T6vsNF54/genomewise.csv')
wt = 'NF54'
mut = 'mut'
sample_ids = fold_change_data.columns[1:]
strain = ['NF54', 'NF54','NF54','mut', 'mut']
treatment = ['Ctl', 'Ctl', 'Ctl','Ctl','Ctl']
experiment_number = [3, 1, 2, 1, 2]
#strain = ['NF54', 'NF54','mut', 'mut']
#treatment = ['Ctl', 'Ctl', 'Ctl','Ctl']
#experiment_number = [1, 2, 1, 2]


metadata = pd.DataFrame({
    'sample_id': sample_ids,
    'strain': strain,
    'treatment': treatment,
    'experiment_number': experiment_number  
})

##########################################################################################################################################


########## START ########################
## Reshape fold change data to long format

fold_change_long = fold_change_data.melt(id_vars='Gene', var_name='sample_id', value_name='fold_change')
print((fold_change_long))
# Merge with metadata
data = fold_change_long.merge(metadata, on='sample_id')
print(data)
# Ensure that fold_change column is numeric (convert non-numeric values to NaN)
data['fold_change'] = pd.to_numeric(data['fold_change'], errors='coerce')

# Drop rows with NaN fold_change values
data = data.dropna(subset=['fold_change'])

# Perform likelihood ratio tests for each gene
results = []
for gene in data['Gene'].unique():
    gene_data = data[data['Gene'] == gene]
    if all(gene_data['fold_change'] <= 0):
        continue
    else:
        # Null model: fold_change ~ experiment_number
        null_model = smf.ols('fold_change ~ 1', data=gene_data).fit()

        # Alternative model: fold_change ~ treatment + experiment_number
        alt_model = smf.ols('fold_change ~ strain', data=gene_data).fit()

        # LRT: Compare the models
        lr_stat = 2 * (alt_model.llf - null_model.llf)
        p_value = stats.chi2.sf(lr_stat, df=1)  # df=1 because 1 additional parameter (treatment)

        # Guardar resultados
        # Calcular la media de fold_change para este gen
        mean_wt = gene_data[gene_data['strain'] == wt]['fold_change'].mean()
        mean_mut = gene_data[gene_data['strain'] == mut]['fold_change'].mean()
        mean_lfc = mean_mut - mean_wt
        results.append({
            'Gene': gene,
            'log_fold_change': mean_lfc,
            'lr_stat': lr_stat,
            'p_value': p_value
        })

# Convert results to DataFrame and display significant results
lrt_results = pd.DataFrame(results)
# Ajustar los p-valores usando el método de FDR (Benjamini-Hochberg)
_, pvals_corrected, _, _ = multipletests(lrt_results['p_value'], method='fdr_bh')

# Añadir la columna de p-valor ajustado
lrt_results['adjusted_p_value'] = pvals_corrected
lrt_results = lrt_results.sort_values('adjusted_p_value')
# Guardar los resultados completos en un archivo CSV (sin filtrar)
lrt_results.to_csv(OUT, index=False)

# Mostrar los resultados sin filtrar
print(lrt_results)

