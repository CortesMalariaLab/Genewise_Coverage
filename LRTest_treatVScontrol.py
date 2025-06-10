import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats
import math
from statsmodels.stats.multitest import multipletests
import sys

###########################################################################################################################################
OUT = "/PROJECTES/MALARIA_EPIGENETICS/Projects/ReCHIP_GDV1_Project_ene25/New_DiffHC/treated/NF54_HS/NF54_HS_results.csv"
fold_change_data = pd.read_csv('/PROJECTES/MALARIA_EPIGENETICS/Projects/ReCHIP_GDV1_Project_ene25/New_DiffHC/treated/NF54_HS/genomewise.csv')
#treat = 'choDep'
#treat = 'DHA'
treat = 'HS'
sample_ids = fold_change_data.columns[1:]
strain = ['NF54', 'NF54','NF54','NF54']
#strain = ['NF54', 'NF54', 'NF54', 'NF54', 'NF54', 'NF54']
#treatment = ['choDep', 'choDep', 'Ctl','Ctl']
#treatment = ['Ctl', 'Ctl', 'DHA', 'Ctl', 'DHA', 'DHA']
#treatment = ['Ctl', 'choDep', 'Ctl', 'choDep', 'Ctl', 'choDep']
#treatment = ['choDep', 'Ctl', 'choDep', 'Ctl']
#treatment = ['Ctl', 'choDep', 'Ctl', 'choDep']
treatment = ['Ctl', 'Ctl','HS',  'HS']
experiment_number = [1, 2, 1, 2]
#experiment_number = [3, 1, 1, 2, 2, 3]

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
        #gene_data.loc[gene_data['fold_change'] <= 0, 'fold_change'] = 0
        lfc_values = list()
        for exp in gene_data['experiment_number'].unique():
            ctl_value = gene_data[(gene_data['experiment_number'] == exp) & (gene_data['treatment'] == 'Ctl')]['fold_change'].values
            choDep_value = gene_data[(gene_data['experiment_number'] == exp) & (gene_data['treatment'] == treat)]['fold_change'].values
            if len(ctl_value) > 0 and len(choDep_value) > 0:
                lfc = choDep_value[0] - ctl_value[0]  # Calcular LFC para cada experimento
                lfc_values.append(lfc)
                # Calcular la media de los LFCs
        if lfc_values:
            mean_lfc = sum(lfc_values) / len(lfc_values)
        else:
            mean_lfc = float('nan')  # Si no se pudieron calcular LFCs
        # Null model: fold_change ~ experiment_number
        null_model = smf.ols('fold_change ~ experiment_number', data=gene_data).fit()

        # Alternative model: fold_change ~ treatment + experiment_number
        alt_model = smf.ols('fold_change ~ experiment_number + treatment', data=gene_data).fit()

        # LRT: Compare the models
        lr_stat = 2 * (alt_model.llf - null_model.llf)
        p_value = stats.chi2.sf(lr_stat, df=1)  # df=1 because 1 additional parameter (treatment)

        # Guardar resultados
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

