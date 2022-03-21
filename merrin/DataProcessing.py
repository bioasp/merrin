#! /usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Kerian Thuillier"
__email__ = "kerian.thuillier@ens-rennes.fr"

#IMPORT#########################################################################

import pandas as pd
from .Parameters import DATA_TYPE_COLUMN, THRESHOLD

#CODE###########################################################################


def drop_repeats(df: pd.DataFrame) -> pd.DataFrame:
    """
    Function:
    Params:
    Return:
    """
    return df.loc[(df.shift() != df).any(1)]


def drop_repeats_transition(df: pd.DataFrame) -> pd.DataFrame:
    """
    Function:
    Params:
    Return:
    """

    # inputs = set(['Carbon1', 'Carbon2', 'Hext', 'Fext', 'Oxygen', 'Growth'])
    # regulators = set(['RPcl', 'RPb', 'RPh', 'RPO2'])

    # cols = sorted(inputs.union(regulators).intersection(df.columns))

    # Filter duplicated transitions
    rows = list(df.iterrows())
    transitions = dict()
    last_row = (None, None)
    for i in range(len(rows) - 1):
        r_current = tuple(rows[i][1].to_list())
        r_next = tuple(rows[i+1][1].to_list())
        t = (r_current, r_next)
        if t not in transitions.keys():
            if last_row[1] == r_current:
                transitions[t] = (last_row[0], i+1)
            else:
                transitions[t] = (i, i+1)
            last_row = (i+1, r_next)
            # transitions[t] = (i, i+1)
    transitions_set = list(set([e for t in transitions.values() for e in t]))
    transitions_set.sort()
    df = df.iloc[transitions_set]
    return df


def bin_normalized_threshold(df: pd.DataFrame) -> pd.DataFrame:
    """
    Function:
    Params:
    Return:
    """
    df_ = df.loc[:, df.columns != DATA_TYPE_COLUMN] \
        / df.loc[:, df.columns != DATA_TYPE_COLUMN].max()    # normalization
    db = (df_ > THRESHOLD).astype(int)                       # binarize
    return drop_repeats_transition(db)


def preprocess_simulations(simulations: dict) -> pd.DataFrame:
    """
    Function:
    Params:
    Return:
    """
    sims = pd.concat({k+1: df for k, df in simulations.items()})

    df = pd.concat({k+1: bin_normalized_threshold(df)
                    for k, df in simulations.items()})

    data = {}
    for et, _ in df.iterrows():
        obs = {k: float(v) if k != DATA_TYPE_COLUMN else v
               for (k, v) in sims.loc[et].items()}
        data[(et[0], et[1])] = obs

    df = pd.DataFrame(data).T

    df = df.where(df.notnull(), None)

    return df


#CODE###########################################################################


def main():
    """
    Function: main function executing the project
    Params: None
    Return: None
    """


#RUN############################################################################
if __name__ == '__main__':
    main()
