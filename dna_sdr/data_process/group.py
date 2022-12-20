import regex as re
import pandas as pd


def group_query(fname):
    re_query = "[a-zA-Z]\d\d"
    match = re.findall(re_query, fname)
    match = [re.sub("0+(?!$)", "", loc) for loc in match]
    group_lst = ["T1"]

    # If the samples are contiunes in the file,i.e., A01-A03 or A09-B03
    if len(match) <= 2:
        row = (match[0][0], match[1][0])
        col = (int(match[0][1:]), int(match[1][1:]))

        if row[0] == row[1]:
            group_quant = col[1] - col[0] + 2
            for i in range(group_quant - 1):
                group_label = row[0] + str(col[0] + i)
                group_lst.append(group_label)

        elif row[0] != row[1]:
            group_quant = 12 - col[0] + col[1] + 2
            for i in range(group_quant - 1):
                if col[0] + i > 12:
                    group_label = row[1] + str(col[0] + i - 12)
                elif col[0] + i <= 12:
                    group_label = row[0] + str(col[0] + i)
                group_lst.append(group_label)

        else:
            raise Exception("Error processing file: {}".format(fname))

    # If the samples are not continues in the file, i.e., A01, A04, A09, B03
    elif len(match) > 2:
        group_quant = len(match) + 1
        group_lst.extend(match)

    return group_quant, group_lst


def group_generation(df, con_tube_number):
    sorted_non_duplicate_tube_number = sorted(
        [int(i) for i in list(set(list(df.columns)))]
    )
    presented_tube_number = [str(i) for i in sorted_non_duplicate_tube_number]
    # Create subgroups of individual conditions (4 tubes per condition continuously)
    counter1 = 1
    counter2 = 1
    counter3 = 0
    number_of_samples = list(range(1, con_tube_number))
    sample_group = []
    while counter1 <= (con_tube_number - 1) / 4:
        subgroup = []
        while counter3 != 4:
            subgroup.append(str(number_of_samples[counter2 - counter1]))
            counter2 += 1
            counter3 += 1
        sample_group.append(subgroup)
        counter1 += 1
        counter2 += 1
        counter3 = 0

    presented_groups = []
    for i in sample_group:
        presented_subgroups = []
        for j in i:
            if j in presented_tube_number:
                presented_subgroups.append(j)
        presented_groups.append(presented_subgroups)
    return sample_group, presented_groups


def screen_grouping(fname):
    df = pd.read_pickle(fname)
    fname_split = fname.split(".")[0].split("_")
    plate_num = fname_split[3]

    group_quant, group_lst = group_query(fname)

    control_tube = group_quant * 4 + 1
    _, presented_group = group_generation(df, control_tube)
    presented_group_dict = dict(zip(group_lst, presented_group))

    return plate_num, presented_group_dict
