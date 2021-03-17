import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

features_table_file = '../results/normal_samples_only/features_table.txt'
wilcoxon_table_file = '../results/normal_samples_only/wilcoxon_smoking_status_2_vs_others.txt'

#gene expression tables
counts_per_million_file = '../data/normal_counts_per_million.txt'
counts_per_million_diff_genes_file = '../data/normal_counts_per_million_diff_genes.txt'

#ate tables
ipw_file = '../results/normal_samples_only/ATE_by_IPW.txt'
ps_matching_ate_table_11_file = '../results/normal_samples_only/ATE_by_ps_match_11.txt'
ps_matching_ate_table_1n_file = '../results/normal_samples_only/ATE_by_ps_match_1n.txt'
cov_matching_ate_table_11_file = '../results/normal_samples_only/ATE_by_cov_match_11.txt'
cov_matching_ate_table_1n_file = '../results/normal_samples_only/ATE_by_cov_match_1n.txt'
genetic_ate_table_file = '../results/normal_samples_only/ate_genetic_sensitivity.txt'
ps_sensitivity_ate_table_file = '../results/normal_samples_only/ate_ps_sensitivity.txt'
cov_sensitivity_ate_table_file = '../results/normal_samples_only/ate_cov_sensitivity.txt'

# matching files
ps_matching_with_replacement_11_file = '../results/normal_samples_only/ps_matching_with_replacement_11.txt'
cov_matching_with_replacement_11_file = '../results/normal_samples_only/cov_matching_with_replacement_11.txt'


def pie(df, col, numerical=False, label_dict=None):
    val_lst = list(df[col])
    if numerical:
        try:
            temp = [int(x) for x in val_lst if x not in ["'--","Unknown",'Not Reported','not reported',np.nan]]
        except:
            temp = [float(x) for x in val_lst if x not in ["'--","Unknown",'Not Reported','not reported',np.nan]]
        val_lst = temp+(len(val_lst)-len(temp))*[np.nan]
        labels = sorted(set(val_lst))
    else:
        labels = sorted(set(val_lst))
    sizes = []
    for l in labels:
        sizes.append(val_lst.count(l))
    plt.figure(figsize=(10,10))
    if label_dict is not None:
        labels = [label_dict[l] for l in labels]
    plt.pie(x=sizes, labels=labels, autopct=lambda p : '{:.2f}%  ({:,.0f})'.format(p,p * sum(sizes)/100))

target_list_symbols_file_name = 'hsapiens_gene_names.csv'
target_list_symbols = pd.read_csv(target_list_symbols_file_name, index_col=0)

def add_gene_names(table):
    table['ENSG'] = [x[:x.find('.')] for x in table.index]
    table = pd.merge(table,target_list_symbols,how='left',left_on='ENSG',right_index=True)
    table = table[['external_gene_name']+list(table.columns[:-1])]
    table = table.drop('ENSG', axis=1)
    return table