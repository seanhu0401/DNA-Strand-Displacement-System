"""
The data input should have the response column to be the last column.
"""
import itertools as it
import os
from functools import reduce

import pandas as pd

os.chdir("./")
FNAME = "dna_sdr/stats/HigginsABC.art.csv"  # 3-way sample data w/ result
df_full = pd.read_csv(FNAME)
df = df_full.loc[:, ["A", "B", "C", "Y"]].copy()

"""
2-way sample data
row = [1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2]
col = [1, 1, 1, 2, 2, 2, 3, 3, 3, 1, 1, 1, 2, 2, 2, 3, 3, 3]
resp = [5, 4, 3, 1, 1, 2, 4, 3, 2, 4, 3, 4, 3, 2, 4, 1, 3, 3]
df = pd.DataFrame({"row": row, "col": col, "resp": resp})
"""

resp = df.iloc[:, -1]
mean = resp.mean()

factors = list(df.columns)[:-1]
n_ways = len(factors)
levels_lst = [list(set(df[factor])) for factor in factors]
level_combo = list(it.product(*levels_lst))

level_full_lst = []
for index, combo_set in enumerate(level_combo):
    combo_lst = [combo_set]
    for j in reversed(range(1, len(combo_set))):
        combo_lst.extend(list(it.combinations(combo_set, r=j)))
    level_full_lst.append(combo_lst)

factor_full_lst = [tuple(factors)]
for i in reversed(range(1, len(factors))):
    factor_full_lst.extend(list(it.combinations(factors, r=i)))

factor_combo_lst = []
for factor in factor_full_lst:
    factor_combo = []
    for length in reversed(range(1, len(factor) + 1)):
        factor_combo.extend(list(it.combinations(list(factor), r=length)))
    factor_combo_lst.append(factor_combo)

invert_levels = list(map(list, zip(*level_full_lst)))
invert_levels = [list(set(levels)) for levels in invert_levels]
factor_level_dict = dict(zip(factor_full_lst, invert_levels))

df_lst = []
for key, value in factor_level_dict.items():
    sublist = []
    for condition in value:
        cond_df = df[df[list(key)].apply(tuple, 1).isin([condition])].copy()
        cond_df[f"{key}"] = cond_df.iloc[:, -1].mean()
        sublist.append(cond_df)
    df_lst.append(pd.concat(sublist).sort_index())

condition_mean_df = reduce(
    lambda left, right: pd.merge(left, right, on=list(df.columns)), df_lst
)
cell_mean = condition_mean_df[str(factor_full_lst[0])]
residual = resp - cell_mean

aligned_df = df.copy()
aligned_column_name = []
for factor_combo in factor_combo_lst:
    n_way_comp = len(factor_combo[0])
    for factor in factor_combo:
        h = n_way_comp - len(factor)
        if h == 0:
            VALUE = condition_mean_df[str(factor)]
        elif h % 2 == 0:
            VALUE += condition_mean_df[str(factor)]
        else:
            VALUE -= condition_mean_df[str(factor)]
    if n_way_comp % 2 == 0:
        VALUE += mean
    else:
        VALUE -= mean
    aligned_column_name.append(f"Aligned {factor_combo[0]}")
    aligned_df[f"Aligned {factor_combo[0]}"] = VALUE + residual

rank_column_name = [name.replace("Aligned", "Rank") for name in aligned_column_name]

rank_df = aligned_df[aligned_column_name].rank()
rank_df.rename(dict(zip(aligned_column_name, rank_column_name)), axis=1, inplace=True)
final_df = pd.concat([aligned_df, rank_df], axis=1)
print(final_df)
