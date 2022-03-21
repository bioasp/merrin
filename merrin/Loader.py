#! /usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Kerian Thuillier"
__email__ = "kerian.thuillier@ens-rennes.fr"

#IMPORT#########################################################################

import bonesis
import pandas as pd
import networkx as nx
from libsbml import SBMLReader
from .Parameters import DATA_TYPE_COLUMN
from .MetabolicNetworks import MetabolicNetwork

#CODE###########################################################################


def load_sbml(sbml_file: str):
    """
    Function: 
    Params: 
    Return: 
    """
    sbmld = SBMLReader().readSBML(sbml_file)
    sbmlm = sbmld.getModel()

    metabolic_network = MetabolicNetwork()
    metabolic_network.init_sbml(sbmlm)

    return metabolic_network


def load_pkn(pkn_file: str, inputs):
    """
    Function: 
    Params: 
    Return: 
    """

    pkn = nx.DiGraph()
    with open(pkn_file) as fp:
        for line in fp:
            line = line.strip()
            if not line:
                continue
            a, b = line.split()
            pkn.add_edge(a, b, sign=0)

    pkn = nx.DiGraph()
    with open(pkn_file) as fp:
        for line in fp:
            line = line.strip()
            if not line:
                continue
            a, b = line.split()
            if b not in inputs:
                pkn.add_edge(a, b, sign=0)
            if a not in inputs:
                pkn.add_edge(b, a, sign=0)

    dom = bonesis.InfluenceGraph(pkn, allow_skipping_nodes=True)

    return dom


def read_simulation(csv_file: str) -> pd.DataFrame:
    """
    Function: load the file_path simulations into a panda dataframe
    Params: file_path, str
    Return: panda dataframe
    """

    df = pd.read_csv(csv_file, sep='\t')
    # df['Time'] = [math.floor(i * 100) if i * 100 - math.floor(i * 100)
    #               < 0.5 else math.ceil(i * 100) for i in df['Time']]
    df.set_index('Time', inplace=True)
    return df


def load_simulations(files_path: str, k: int = 0, tag: str = None) -> pd.DataFrame:
    """
    Function: load the file_path simulations into a panda dataframe
    Params: file_path, str
            k, int, initial value of the counter
    Return: panda dataframe
    """

    simulations = {k + i: read_simulation(f)
                   for i, f in enumerate(sorted(files_path))}

    if tag is not None:
        for sim in simulations.values():
            sim[DATA_TYPE_COLUMN] = tag

    return simulations


def main():
    """
    Function: main function executing the project
    Params: None
    Return: None
    """


#RUN############################################################################
if __name__ == '__main__':
    main()
