import pandas as pd
from nupack import *


""" 
1X PBS Solution Concentrations:
    - NaCl: 137 mM
    - KCl: 2.7 mM
    - Na2HPO4: 10 mM
    - KH2PO4: 1.8 mM
1X Folding Buffer Concentration:
    - Tris-base: 40 mM
    - Acetic acid: 20 mM
    - EDTA: 2 mM
    - MgCl2: 12 mM
"""


if __name__ == "__main__":
    folding_buffer_mg = 0.012  # [M] 1X folding buffer
    pbs_na = (137 + 10 * 2) / 1000  # [M] 1X PBS
    final_mg_concentration = folding_buffer_mg / 10  # 10x dilution in qPCR

    # Salt is in [M]
    model1 = Model(
        material="dna",
        ensemble="stacking",
        celsius=40,
        sodium=pbs_na,
        magnesium=final_mg_concentration,
    )

    # trig_lst = list()
    # name_lst = list()
    # fname = "./dna_sdr/DNA_Strands.xlsx"
    # xls = pd.ExcelFile(fname)
    # for sname in xls.sheet_names:
    #     dna_seq = pd.read_excel(fname, sheet_name=sname)
    #     for seq in dna_seq["Seqences (5' to 3')"]:
    #         Trig_list.append(seq)
    #     for name in dna_seq["Name"]:
    #         Name_list.append(name)

    fname = "./dna_sdr/pickles/trig_info.pkl"
    trig_info_df = pd.read_pickle(fname)
    trig_lst = trig_info_df["seq"].tolist()
    name_lst = trig_info_df["name"].tolist()

    trig_dict = dict(zip(name_lst, trig_lst))

    base = Strand(
        "TG GGTGG TG AGA TG GATTG TG AGA TG TG AGA CAT ACA GCG CCG ACC GTA".replace(
            " ", ""
        ),
        name="base",
    )
    incumb = Strand(("CA CAATC CA TCT CA CCACC CA".replace(" ", "")), name="incumb")

    complex_list = list()
    free_energy_list = list()
    seq_list = list()
    conc_list = list()

    for i in range(len(trig_lst)):
        trig = Strand((trig_lst[i]), name=name_lst[i])
        trig_base_complex = "({}+base)".format(name_lst[i])
        t1 = Tube(
            strands={base: 5e-7, incumb: 5e-7, trig: 5e-7},
            complexes=SetSpec(max_size=2),
            name="t1",
        )
        complex_result = complex_analysis(complexes=t1, model=model1, compute=["pfunc"])
        complex_list.append(trig_base_complex)
        free_energy_list.append(complex_result[trig_base_complex].free_energy)
        seq_list.append((trig_lst[i].replace(" ", "")))
        complex_conc = complex_concentrations(tube=t1, data=complex_result)
        for complex, conc in complex_conc["t1"].complex_concentrations.items():
            if complex.name == trig_base_complex:
                conc_list.append(conc * 10**9)

    conc_dict = dict(zip(complex_list, conc_list))
    conc_df = pd.DataFrame.from_dict(conc_dict, orient="index", columns=["Conc (nM)"])
    conc_df.to_pickle("./dna_sdr/pickles/concentration.pkl")

    complex_list.append("(incumb+base)")
    seq_list.append(("CA CAATC CA TCT CA CCACC CA".replace(" ", "")))

    free_energy_list.append(complex_result["(incumb+base)"].free_energy)
    free_energy_dict = dict(zip(complex_list, free_energy_list))
    free_energy_df = pd.DataFrame.from_dict(
        free_energy_dict, orient="index", columns=["dG (kcal/mol)"]
    )
    # free_energy_df["Seqences (5' to 3')"] = seq_list
    free_energy_df.to_pickle("./dna_sdr/pickles/free_energy.pkl")
