"""
The data input should have the response column to be the last column.
"""
import math
import os
import itertools as it
import numpy as np
import pandas as pd

# import scipy.stats as stats


os.chdir("./")
FNAME = "dna_sdr/stats/HigginsABC.art.csv"  # 3-way sample data w/ result
df_full = pd.read_csv(FNAME)
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

term_lst = []
for i in range(n_ways):
    term = int(
        math.factorial(n_ways) / (math.factorial(i) * math.factorial(n_ways - i))
    )
    if i == 0:
        term_loc = term
    else:
        term_loc = term_lst[-1] + term
    term_lst.append(term_loc)

level_full_lst = []
for index, combo_set in enumerate(level_combo):
    combo_lst = [combo_set]
    for j in reversed(range(1, len(combo_set))):
        combo_lst.extend(list(it.combinations(combo_set, r=j)))
    level_full_lst.append(combo_lst)


factor_full_lst = [tuple(factors)]
for i in reversed(range(1, len(factors))):
    factor_full_lst.extend(list(it.combinations(factors, r=i)))

df_lst = []
for cond in level_combo:
    cond_df = df[df[factors].apply(tuple, 1).isin([cond])].copy()
    mean_response = cond_df.iloc[:, -1].mean()
    cond_df["residual"] = cond_df.iloc[:, -1] - mean_response
    df_lst.append(cond_df)

df = pd.concat(df_lst).sort_index()

matrix = np.zeros((term_lst[-1], term_lst[-1] + 1))
matrix[:, -1] = mean

full_est_effect_mean_lst = []

# Mean calculation for each level of effect combination
for factor_level in level_full_lst:
    est_effect_mean_lst = []
    for index, level in enumerate(factor_level):
        factor_combo = factor_full_lst[index]
        df_copy = df.loc[:, factor_combo].copy()
        factor_copy = df_copy[df_copy.loc[:, factor_combo] == level].dropna()
        factor_indexes = factor_copy.index
        est_effect_mean_lst.append(resp[factor_indexes].mean())
    full_est_effect_mean_lst.append(est_effect_mean_lst)

# Step 2: compute estimated effects for main + interaction effects
# TODO: generalize it - for loop over each level combo
test_set = full_est_effect_mean_lst[0]
param_array = np.zeros(len(test_set))

for index, _ in enumerate(term_lst):
    if index - 1 < 0:
        param_array[0] = test_set[0]
    else:
        param_array[term_lst[index - 1]] = sum(
            test_set[term_lst[index - 1] : term_lst[index]]
        )

for index, result in enumerate(test_set):
    matrix[0:index, index] = param_array[index]
    matrix[index, index] = result

# TODO: Try to reduce the nested loop
value_lst = []
for index, _ in enumerate(term_lst):
    ways = n_ways - index
    if index == 0:
        sub_matrix = [matrix[index]]
    else:
        sub_matrix = matrix[term_lst[index - 1] : term_lst[index]]

    for array in sub_matrix:
        VALUE = 0
        COUNTER = 1
        for h in range(index, n_ways):
            if VALUE == 0 and array[0] != 0:
                VALUE = array[0]
            elif VALUE == 0 and array[0] == 0:
                VALUE = sum(array[: term_lst[h]])

            if COUNTER % 2 != 0 and h < (n_ways - 1):
                VALUE -= sum(array[term_lst[h] : term_lst[h + 1]])
            elif COUNTER % 2 == 0 and h < (n_ways - 1):
                VALUE += sum(array[term_lst[h] : term_lst[h + 1]])
            else:
                pass
            COUNTER += 1

        if ways % 2 == 0:
            VALUE += array[-1]
        else:
            VALUE -= array[-1]

        value_lst.append(VALUE)

# value_test_lst = [-1.25, -130.25, -124.75, -118.25, -10.25, -4.5, -5.0]
# print(value_lst)
# print(value_test_lst == value_lst)

# df_lst = list()

# for cond_num in range(len(combo)):
#     matrix = matrice_lst[cond_num]
#     i_index = 0
#     cond = combo[cond_num]
#     print(cond)
#     cond_df = df[df[factors].apply(tuple, 1).isin([cond])].copy()
#     factor_mean = cond_df["resp"].mean()
#     for i in range(len(cond)):
#         i_index += 1
#         print(cond[i])
#         cond_i_df = df[df[factors[i]] == cond[i]]
#         cond_i_mean = cond_i_df["resp"].mean()
