#calculate_zRE.py for calculating z-normalized Release Effect

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import collections
from sklearn import linear_model
import time
import glob
import os

def ciber(li, an, rep):
    start = time.time()
    #select count files
    all_count_files = glob.glob('count*.csv')
    selected_count_files = [i for i in all_count_files if li in i and an in i]

    #calculate relative gRNA read
    df_all = [pd.read_csv(i, index_col = 0) for i in selected_count_files]
    df_all_composition = []
    for d in df_all:
        t = np.sum(d)
        e = len(d)
        df_all_composition.append(d/t*e)

    df_name = [i.lower() for i in selected_count_files]
    index = df_all[0].index

    #classify files 
    cell_cas = [i for i, j in zip(df_all_composition, df_name) if 'cell' in j and 'cas' in j]
    cell_non = [i for i, j in zip(df_all_composition, df_name) if 'cell' in j and 'non' in j]
    sEVs_cas = [i for i, j in zip(df_all_composition, df_name) if 'sevs' in j and 'cas' in j]
    sEVs_non = [i for i, j in zip(df_all_composition, df_name) if 'sevs' in j and 'non' in j]
    print('loaded:')
    for i, j in zip(selected_count_files, df_name):
        if 'cell' in j and 'cas' in j:
            print('\t' + i + '\tas cell_cas')
        if 'cell' in j and 'non' in j:
            print('\t' + i + '\tas cell_non')
        if 'sevs' in j and 'cas' in j:
            print('\t' + i + '\tas sEVs_cas')
        if 'sevs' in j and 'non' in j:
            print('\t' + i + '\tas sEVs_non')

    #calculate the average of Cas9- ('non') samples
    cell_non_mean = pd.DataFrame(np.mean([i.values for i in cell_non], axis = 0), index = index)
    sEVs_non_mean = pd.DataFrame(np.mean([i.values for i in sEVs_non], axis = 0), index = index)

    #remove gRNAs with compositions lower than 0.005
    cell_cas = [i[i>=0.005].dropna() for i in cell_cas]
    cell_non_mean = cell_non_mean[cell_non_mean >= 0.005].dropna()
    sEVs_cas = [i[i>=0.005].dropna() for i in sEVs_cas]
    sEVs_non_mean = sEVs_non_mean[sEVs_non_mean >= 0.005].dropna()

    #calculate the Viability and Release for each replicate
    def calculate_ratio(cas, non, rep):
        merged = pd.concat([cas[rep - 1], non], axis = 1, join = 'inner')
        merged.columns = ['cas', 'non']
        ratio = pd.DataFrame(np.log2(merged['cas'])-np.log2(merged['non']))      
        return ratio

    ratio_cell = calculate_ratio(cell_cas, cell_non_mean, rep = rep)
    ratio_sEVs = calculate_ratio(sEVs_cas, sEVs_non_mean, rep = rep)
    gRNA_ratio = pd.concat([ratio_cell, ratio_sEVs], axis = 1, join = 'inner', sort = True)
    gRNA_ratio.columns = ['FCcells', 'FCsEVs']
    
    #remove genes with gRNAs less than 3
    row_genes = [i.split('_')[0] for i in gRNA_ratio.index]
    count_genes = collections.Counter(row_genes)
    genes_with_3gRNAs = [i for i, j in count_genes.items() if j > 2]
    index_valid = [i for i in gRNA_ratio.index if i.split('_')[0] in genes_with_3gRNAs]
    gRNA_ratio = gRNA_ratio.loc[index_valid, :]        

    #calculate Scores
    #plot Viability and Release
    def linear_regression(viability, release, title):
        fig, ax = plt.subplots(figsize = (4, 4))
        ax.grid(axis = 'both', color = '0.8', linewidth = 1, zorder = 0)
        plt.rcParams['axes.axisbelow'] = True
        ax.set_xlabel(r'FC$_{cells}$', size = 14)
        ax.set_ylabel(r'FC$_{sEVs}$', size = 14)
        ax.set_title(title, size = 16, verticalalignment = 'bottom')
        ax.scatter(viability, release, color = '0.2', alpha = 0.3, s = 15, linewidth = 0, label = None)
        safe_id = [guide for guide in viability.index if 'safe' in guide]
        ax.scatter(viability.loc[safe_id], release.loc[safe_id], color = 'b', s = 15, linewidth = 0, label = 'safe')
        none_id = [guide for guide in viability.index if 'none' in guide]
        ax.scatter(viability.loc[none_id], release.loc[none_id], color = 'r', s = 15, linewidth = 0, label = 'none')
        maxi = max(max(viability), max(release))
        mini = min(min(viability), min(release))
        ax.set_xlim(mini-0.3, maxi+0.3)
        ax.set_ylim(mini-0.3, maxi+0.3)
        ticks = np.arange(round(mini-0.3, 0), maxi+0.3, 2)
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)

        #linear regression
        clf = linear_model.LinearRegression()
        x_cor = np.matrix(viability).T 
        y_cor = np.matrix(release).T
        clf.fit(x_cor, y_cor)
        alpha = clf.intercept_[0]
        beta = clf.coef_[0, 0]
        r2 = clf.score(x_cor, y_cor)
        r = np.sqrt(r2) 
        ax.plot(x_cor, clf.predict(x_cor), c="lightgreen", linestyle="solid", linewidth = 2.5, 
        label = f"y = {round(alpha, 2)} + {round(beta, 2)}x\nR = {round(r, 2)}")   
        ax.legend(bbox_to_anchor=(1, 1), loc = 'upper left', fontsize = 12, facecolor = '1')
        plt.show()
        return alpha, beta

    print('\ncalculating Release Effect at gRNA level...')
    v = gRNA_ratio['FCcells']
    s = gRNA_ratio['FCsEVs']
    alpha, beta = linear_regression(v, s, 'Linear regression_not trimmed')
    gRNA_ratio['RE at gRNA level_not trimmed'] = s-alpha-beta*v
    
    print('removing gRNAs with top/bottom REs for each gene...')
    genes = list(dict.fromkeys([i.split('_')[0] for i in gRNA_ratio.index]))
    for g in genes:
        gRNA_for_gene = [guide for guide in gRNA_ratio.index if guide.split('_')[0] == g]
        ranked_id = gRNA_ratio.loc[gRNA_for_gene, :].sort_values(by='RE at gRNA level_not trimmed').index
        gRNA_ratio.drop(ranked_id[0], inplace = True)
        gRNA_ratio.drop(ranked_id[-1], inplace = True)

    v = gRNA_ratio['FCcells']
    s = gRNA_ratio['FCsEVs']
    alpha_trimmed, beta_trimmed = linear_regression(v, s, 'Linear regression_trimmed')
    gRNA_ratio['RE at gRNA level_trimmed'] = s-alpha_trimmed-beta_trimmed*v
    gRNA_ratio.to_excel(f'RE at gRNA level_{li}_{an}_rep{rep}.xlsx')

    print('\ncalculating RE at gene level...')
    gRNA_ratio.index = [v.split("_")[0] for v in gRNA_ratio.index]
    d = gRNA_ratio.groupby(axis = 0, level = 0).median()
    m = np.mean(d['RE at gRNA level_trimmed'])
    sd = np.std(d['RE at gRNA level_trimmed'])
    d['z-RE_trimmed'] = (d['RE at gRNA level_trimmed']-m)/sd
    d['z-RE_trimmed'].to_excel(f'RE at gene level_{li}_{an}_rep{rep}.xlsx')

    print('completed')
    end = time.time()
    h, mod = divmod(end - start, 3600)
    m, s = divmod(mod, 60)
    print(f"{int(h)} hr {int(m)} min {round(s, 4)} sec")
    return

path = r'' #dir of count files
os.chdir(path)

library = ['ACOC', 'DTKP', 'PROT', 'TMMO']
anchor = ['CD63', 'CD9']
for rep in [1,2]:
    ciber(library[1], anchor[0], rep)

