import pandas as pd


def total_interaction_of_residue(interaction_info):
    total_score = 0

    for _, _, score in interaction_info:
        total_score += score

    return total_score


def total_interaction(interaction_info):
    interaction = []
    for segid, resid, resn, pairs in interaction_info, :
        for atom1, atom2, dist in pairs:
            interaction.append((segid, resid, resn, atom1, atom2, dist))

    return interaction


def pose_feature(all_interactions):
    feature = []
    for interaction in all_interactions:
        feature += total_interaction(interaction)
        #score = total_interaction_of_residue(interaction[3])
        #feature.append((interaction[0], interaction[1], interaction[2], score))

    return feature


def create_feature_table(**kwargs):

    #column = ["Segid", "Resid", "Resn", "PairI", "PairJ"]
    column = ["Segid", "Resid", "Resn", "PairI", "PairJ","Surface_d"]
    table = pd.DataFrame(columns=column)

    for key, value in kwargs.items():
        column.append(key)
        #df1 = pd.DataFrame(value, columns=["Segid", "Resid", "Resn", "PairI", "PairJ", key, "Surface_d"])
        # Add surface distance
        df1 = pd.DataFrame(value, columns=["Segid", "Resid", "Resn", "PairI", "PairJ", key, "Surface_d"])
        table = pd.merge(table, df1, how='outer', sort=False)

    table = table[column]
    table = table.sort_values(["Segid", "Resid"])
    table = table.fillna(0)

    return table
