import pandas as pd


if __name__ == "__main__":
    trig_pickle = "./dna_sdr/pickles/trig_param.pkl"  # Fitting result
    conc_pickle = "./dna_sdr/pickles/concentration.pkl"  # Nupack analysis
    trig_info_pickle = "./dna_sdr/pickles/trig_info.pkl"  # Trigger design information

    trig_cf_df = pd.read_pickle(trig_pickle)
    trig_info_df = pd.read_pickle(trig_info_pickle)
    conc_df = pd.read_pickle(conc_pickle)

    # Reordering the columns for the dataframe
    cols = trig_cf_df.columns.tolist()
    cols.insert(0, "plate_loc")
    cols.pop(6)
    trig_cf_df = trig_cf_df[cols]

    trig_info_cf_df = pd.merge(trig_info_df, trig_cf_df, on="plate_loc")

    # Recover the trigger name from the nupack analysis
    name_lst = list()
    for index in conc_df.index:
        name = index.split("+")[0][1:]
        name_lst.append(name)

    conc_df["name"] = name_lst
    conc_df.reset_index(inplace=True)
    merged_df = pd.merge(trig_info_cf_df, conc_df, on="name")
    merged_df.drop("index", axis=1, inplace=True)
    merged_df.rename(columns={"Conc (nM)": "Nupack"}, inplace=True)
    merged_df["percent diff"] = (
        (merged_df["plateau"] - merged_df["Nupack"]) / merged_df["Nupack"] * 100
    )

    print(merged_df)
    # merged_df.to_csv("merged_df.csv")
