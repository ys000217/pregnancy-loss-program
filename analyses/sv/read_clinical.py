import pandas as pd
df = pd.read_excel('/mnt/e/genotype_data/临床信息表0.2.xlsx')
print('SHAPE', df.shape)
print('COLS:', list(df.columns))
print(df.head(8).to_string())
print('---- dtypes ----')
print(df.dtypes.to_string())
