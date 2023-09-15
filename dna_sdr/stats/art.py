import numpy as np
import pandas as pd
import scipy.stats as stats
import itertools as it
import math
import os

"""
The data input should have the response column to be the last column.
"""
os.chdir("./")
fname = "./Ignore/code/HigginsABC.art.csv"  # 3-way sample data w/ result
df_full = pd.read_csv(fname)
df = df_full.loc[:, ["A", "B", "C", "Y"]].copy()
# print(df)

"""
2-way sample data
row = [1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2]
col = [1, 1, 1, 2, 2, 2, 3, 3, 3, 1, 1, 1, 2, 2, 2, 3, 3, 3]
resp = [5, 4, 3, 1, 1, 2, 4, 3, 2, 4, 3, 4, 3, 2, 4, 1, 3, 3]
df = pd.DataFrame({"row": row, "col": col, "resp": resp})
"""

resp = df.iloc[:, -1]
mean = resp.mean()

# Step 1: compute residual
factors = list(df.columns)[:-1]
n_ways = len(factors)
levels_lst = [list(set(df[factor])) for factor in factors]
level_combo = list(it.product(*levels_lst))

term_lst = list()
for i in range(n_ways):
    term = int(
        math.factorial(n_ways) / (math.factorial(i) * math.factorial(n_ways - i))
    )
    if i == 0:
        term_loc = term
    else:
        term_loc = term_lst[-1] + term
    term_lst.append(term_loc)

level_full_lst = list()
for i in range(len(level_combo)):
    combo_set = level_combo[i]
    combo_lst = [combo_set]
    for j in reversed(range(1, len(combo_set))):
        combo_lst.extend(list(it.combinations(combo_set, r=j)))
    level_full_lst.append(combo_lst)

factor_full_lst = [tuple(factors)]
for i in reversed(range(1, len(factors))):
    factor_full_lst.extend(list(it.combinations(factors, r=i)))

df_lst = list()
for cond in level_combo:
    cond_df = df[df[factors].apply(tuple, 1).isin([cond])].copy()
    mean_response = cond_df.iloc[:, -1].mean()
    cond_df["residual"] = cond_df.iloc[:, -1] - mean_response
    df_lst.append(cond_df)

df = pd.concat(df_lst).sort_index()

matrix = np.zeros((term_lst[-1], term_lst[-1] + 1))
matrix[:, -1] = mean

full_est_effect_mean_lst = list()

# Mean calculation for each level of effect combination
for factor_level in level_full_lst:
    est_effect_mean_lst = list()
    for location in range(len(factor_level)):
        factor_combo = factor_full_lst[location]
        level = factor_level[location]
        df_copy = df.loc[:, factor_combo].copy()
        factor_copy = df_copy[df_copy.loc[:, factor_combo] == level].dropna()
        factor_indexes = factor_copy.index
        est_effect_mean_lst.append(resp[factor_indexes].mean())
    full_est_effect_mean_lst.append(est_effect_mean_lst)

test_set = full_est_effect_mean_lst[0]
param_array = np.zeros(len(test_set))

for i in range(len(term_lst)):
    if i - 1 < 0:
        param_array[0] = test_set[0]
    else:
        param_array[term_lst[i - 1]] = sum(test_set[term_lst[i - 1] : term_lst[i]])

for i in range(len(test_set)):
    matrix[0:i, i] = param_array[i]
    matrix[i, i] = test_set[i]

print(matrix)
print(term_lst)
# print(matrix[term_lst[0] : term_lst[1]])

for i in range(len(term_lst)):
    print()

"""
If the leading number of the set is zero, N-way is decrease by 1
h increase each interval as the location group is updated and resets to
0 after each loop. If leading number of the set is zero, h does not increase.
"""

"""
number_of_ways = len(factors)
counts = 0
for r in range(number_of_ways):
    num = math.factorial(number_of_ways)
    den = math.factorial(r) * math.factorial(number_of_ways - r)
    counts += num / den

coef_matrix = [np.zeros((int(counts), int(counts)))]
matrice_lst = coef_matrix * len(combo)
"""

# Step 2: compute estimated effects for main + interaction effects
# df_lst = list()

"""
for cond_num in range(len(combo)):
    matrix = matrice_lst[cond_num]
    i_index = 0
    cond = combo[cond_num]
    print(cond)
    cond_df = df[df[factors].apply(tuple, 1).isin([cond])].copy()
    factor_mean = cond_df["resp"].mean()
    for i in range(len(cond)):
        i_index += 1
        print(cond[i])
        cond_i_df = df[df[factors[i]] == cond[i]]
        cond_i_mean = cond_i_df["resp"].mean()
"""
