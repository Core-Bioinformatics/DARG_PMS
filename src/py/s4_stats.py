import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import fdrcorrection

pd.options.display.float_format = '{:.10f}'.format

absinta = pd.read_csv('absinta.csv')
schirmer = pd.read_csv('schirmer.csv')

cols = []

for col in absinta.columns[:-1]:
    mat = [
        [absinta[col].iloc[2], absinta[col].iloc[0] - absinta[col].iloc[2]],
        [absinta[col].iloc[3], absinta[col].iloc[1] - absinta[col].iloc[3]],
    ]
    print([col, 'senescence', fisher_exact(mat).pvalue])
    cols.append([col, 'senescence', fisher_exact(mat).pvalue])

    mat = [
        [absinta[col].iloc[4], absinta[col].iloc[0] - absinta[col].iloc[4]],
        [absinta[col].iloc[5], absinta[col].iloc[1] - absinta[col].iloc[5]],
    ]
    print([col, 'interferon', fisher_exact(mat).pvalue])
    cols.append([col, 'interferon', fisher_exact(mat).pvalue])

for col in schirmer.columns[:-1]:
    mat = [
        [schirmer[col].iloc[2], schirmer[col].iloc[0] - schirmer[col].iloc[2]],
        [schirmer[col].iloc[3], schirmer[col].iloc[1] - schirmer[col].iloc[3]],
    ]
    print([col, 'senescence', fisher_exact(mat).pvalue])
    cols.append([col, 'senescence', fisher_exact(mat).pvalue])

    mat = [
        [schirmer[col].iloc[4], schirmer[col].iloc[0] - schirmer[col].iloc[4]],
        [schirmer[col].iloc[5], schirmer[col].iloc[1] - schirmer[col].iloc[5]],
    ]
    print([col, 'interferon', fisher_exact(mat).pvalue])
    cols.append([col, 'interferon', fisher_exact(mat).pvalue])

df = pd.DataFrame.from_records(cols, columns=['condition', 'set', 'pval'])
df['adj_pval'] = fdrcorrection(df['pval'])[1]
df.to_csv('stats.csv', index=False)

cols = []

for df_name, df in zip(['absinta', 'schirmer'], [absinta, schirmer]):
    for col1 in df.columns[:-1]:
        for col2 in df.columns[:-1]:
            if col1 == col2:
                continue
            mat = [
                [df[col1].iloc[2], df[col1].iloc[0] - df[col1].iloc[2]],
                [df[col2].iloc[2], df[col2].iloc[0] - df[col2].iloc[2]],
            ]
            entry = [df_name, col1, col2, 'senescence-rglike', fisher_exact(mat).pvalue]
            print(entry)
            cols.append(entry)

            mat = [
                [df[col1].iloc[3], df[col1].iloc[1] - df[col1].iloc[3]],
                [df[col2].iloc[3], df[col2].iloc[1] - df[col2].iloc[3]],
            ]
            entry = [df_name, col1, col2, 'senescence-nonrglike', fisher_exact(mat).pvalue]
            print(entry)
            cols.append(entry)
    
            mat = [
                [df[col1].iloc[4], df[col1].iloc[0] - df[col1].iloc[4]],
                [df[col2].iloc[4], df[col2].iloc[0] - df[col2].iloc[4]],
            ]
            entry = [df_name, col1, col2, 'interferon-rglike', fisher_exact(mat).pvalue]
            print(entry)
            cols.append(entry)
    
            mat = [
                [df[col1].iloc[5], df[col1].iloc[1] - df[col1].iloc[5]],
                [df[col2].iloc[5], df[col2].iloc[1] - df[col2].iloc[5]],
            ]
            entry = [df_name, col1, col2, 'interferon-nonrglike', fisher_exact(mat).pvalue]
            print(entry)
            cols.append(entry)
    
df = pd.DataFrame.from_records(cols, columns=['experiment', 'col1', 'col2', 'set', 'pval'])
df = df.drop_duplicates(subset=['experiment', 'set', 'pval'])
df['adj_pval'] = fdrcorrection(df['pval'])[1]
df.to_csv('stats2.csv', index=False)
