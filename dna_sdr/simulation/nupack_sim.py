"""
DNA analysis using NUPACK

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

import nupack
import pandas as pd

if __name__ == "__main__":
    # The sum of the concentrations of (monovalent) sodium, potassium, and ammonium ions,
    # [Na+]+[K+]+[NH+4], is specified in units of molar using the keyword sodium.

    MONOVALENT = (137 + 10 * 2 + 1.8 + 2.7) / 1000  # [M] 1X PBS
    DIVALENT = 0.0012  # 0.1x folding buffer in qPCR

    # Salt is in [M]
    model1 = nupack.Model(
        material="dna",
        celsius=40,
        sodium=MONOVALENT,
        magnesium=DIVALENT,
    )

    FNAME = "./dna_sdr/pickles/trig_info.pkl"
    trig_info_df: pd.DataFrame = pd.read_pickle(FNAME)
    trig_lst: list[str] = trig_info_df["seq"].tolist()
    name_lst: list[str] = trig_info_df["name"].tolist()

    trig_dict = dict(zip(name_lst, trig_lst))

    T1_trig_dict = {k: v for k, v in trig_dict.items() if "T3" not in k}
    # T3_trig_dict = {k: v for k, v in trig_dict.items() if "T3" in k}

    base = nupack.Strand(
        "TG GGTGG TG AGA TG GATTG TG AGA TG TG AGA CAT ACA GCG CCG ACC GTA".replace(
            " ", ""
        ),
        name="base",
    )

    # base_T3 = nupack.Strand(
    #     "TG TTGGA TG AGA TG AGTTT TG AGA TG GG AGA CTC TGA AGG GAT GGC GCA".replace(
    #         " ", ""
    #     ),
    #     name="base_t3",
    # )

    incumb = nupack.Strand(
        ("CA CAATC CA TCT CA CCACC CA".replace(" ", "")), name="incumb"
    )

    # incumb_T3 = nupack.Strand(
    #     ("CA AAACT CA TCT CA TCCAA CA".replace(" ", "")), name="incumb_t3"
    # )

    complex_list: list[str] = []
    free_energy_list: list[float] = []
    seq_list: list[str] = []
    conc_list: list[float] = []

    for tname, trigger in T1_trig_dict.items():
        trig = nupack.Strand(trigger, name=tname)
        trig_base_complex = f"({tname}+base)"
        t1 = nupack.Tube(
            strands={base: 5e-7, incumb: 5e-7, trig: 5e-7},
            complexes=nupack.SetSpec(max_size=3),
            name="t1",
        )
        complex_result = nupack.complex_analysis(
            complexes=t1, model=model1, compute=["pfunc"]
        )
        complex_list.append(trig_base_complex)
        free_energy_list.append(complex_result[trig_base_complex].free_energy)
        seq_list.append((trigger.replace(" ", "")))
        complex_conc = nupack.complex_concentrations(tube=t1, data=complex_result)
        for complex_sturc, conc in complex_conc["t1"].complex_concentrations.items():
            if complex_sturc.name == trig_base_complex:
                conc_list.append(conc * 10**9)

    conc_dict = dict(zip(complex_list, conc_list))
    conc_df = pd.DataFrame.from_dict(conc_dict, orient="index", columns=["Conc (nM)"])

    # conc_df.to_pickle("./dna_sdr/pickles/concentration_T1.pkl")
    # conc_df.to_csv("concentration_T1.csv")

    complex_list.append("(incumb+base)")
    seq_list.append(("CA CAATC CA TCT CA CCACC CA".replace(" ", "")))

    free_energy_list.append(complex_result["(incumb+base)"].free_energy)
    free_energy_dict = dict(zip(complex_list, free_energy_list))
    free_energy_df = pd.DataFrame.from_dict(
        free_energy_dict, orient="index", columns=["dG (kcal/mol)"]
    )
    free_energy_df["Seqences (5' to 3')"] = seq_list
    # free_energy_df.to_pickle("./dna_sdr/pickles/free_energy_T1.pkl")
    # free_energy_df.to_csv("free_energy_T1.csv")
