#! /usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Kerian Thuillier"
__email__ = "kerian.thuillier@ens-rennes.fr"

#IMPORT#########################################################################
import os
import math
import biolqm
import random
import pandas as pd
from libsbml import SBMLReader
from scipy import rand

#CONSTANT#######################################################################

BIOMASS = 'biomass'
OBJECTIVE = 'Growth'

THRESHOLD = 10**-5

#CODE###########################################################################


#LOADER#########################################################################


def read_simulation(csv_file: str) -> pd.DataFrame:
    df = pd.read_csv(csv_file, sep='\t')
    df['Time'] = [math.floor(i * 100) if i * 100 - math.floor(i * 100)
                  < 0.5 else math.ceil(i * 100) for i in df['Time']]
    df.set_index('Time', inplace=True)
    return df


def shift_external_metabolites(simulations, input_metabolites):
    to_shift = list(input_metabolites) + ['biomass']
    for _, df in simulations.items():
        df[to_shift] = df[to_shift].shift(1)
        df.loc[0, to_shift] = df.loc[1, to_shift]


def get_file_name(file_path: str) -> str:
    return os.path.split(os.path.splitext(file_path)[0])[1]


def load_simulations(files_path: str) -> pd.DataFrame:
    """
    Function: load the file_path simulations into a panda dataframe
    Params: file_path, str
    Return: panda dataframe
    """

    simulations = {get_file_name(f): read_simulation(f)
                   for i, f in enumerate(sorted(files_path))}

    return simulations


def load_regulations(reg_file: str):
    """
    """

    lqm = biolqm.load(f'{reg_file}')
    bn = biolqm.to_minibn(lqm)

    return bn


def load_sbml(sbml_file: str) -> tuple:
    """
    """

    sbmld = SBMLReader().readSBML(f'{sbml_file}')
    sbmlm = sbmld.getModel()

    reactants = set()
    products = set()
    reactions = set()
    reaction_metabo_relation = set()
    for reaction in sbmlm.getListOfReactions():
        name = reaction.getId()
        reactions.add(name)
        reactants.update([a.getSpecies()
                         for a in reaction.getListOfReactants()])
        reaction_metabo_relation.update(
            (name, a.getSpecies()) for a in reaction.getListOfReactants())
        products.update([a.getSpecies() for a in reaction.getListOfProducts()])
        reaction_metabo_relation.update(
            (name, a.getSpecies()) for a in reaction.getListOfProducts())
        assert not reaction.getListOfModifiers(), 'Not implemented'

    inputs = reactants.difference(products)
    outputs = products.difference(reactants)
    metabolites = reactants.union(products).difference(inputs.union(outputs))

    input_reactions = set(r
                          for (r, m) in reaction_metabo_relation
                          if m in inputs or m in outputs)

    return reactions, inputs, outputs, metabolites, input_reactions


#GENERATOR######################################################################


def binarize_cols(simulations: dict, cols: list, threshold: float = 10**-5):
    """
    """
    bin_simulations = {}
    for i, sim in simulations.items():
        sim[cols] = (sim[cols] > threshold).astype(int)
        bin_simulations[i] = sim
    return bin_simulations


def remove_cols(simulations, cols):
    """
    """
    bin_simulations = {}
    for i, sim in simulations.items():
        sim.drop(cols, inplace=True, axis=1)
        bin_simulations[i] = sim
    return bin_simulations


def generate_data(sim_files, sbml_file, reg_file,
                  fluxomic=True, transcriptomic=True, cinetic=True,
                  hide_reactions=False, obj=OBJECTIVE) -> dict:
    """
    """
    simulations = load_simulations(sim_files)
    reactions, inputs, outputs, metabolites, input_r = load_sbml(sbml_file)
    bn = load_regulations(reg_file)

    shift_external_metabolites(simulations, inputs)

    gen_simulations = simulations

    #--------------------------------------------------------------------------#
    # Inputs metabolites
    #--------------------------------------------------------------------------#
    if not fluxomic and not cinetic:
        inputs_l = list(inputs)
        outputs_l = list(outputs)
        gen_simulations = binarize_cols(gen_simulations, inputs_l)
        gen_simulations = binarize_cols(gen_simulations, outputs_l)
        gen_simulations = remove_cols(gen_simulations, [BIOMASS])

    #--------------------------------------------------------------------------#
    # Reactions
    #--------------------------------------------------------------------------#
    reactions_l = list(r for r in reactions if r != obj)
    internal_reactions_l = [r for r in reactions_l
                            if r not in input_r and r != obj]
    if (not fluxomic and not transcriptomic):
        gen_simulations = remove_cols(gen_simulations, reactions_l)
    elif hide_reactions:
        gen_simulations = remove_cols(gen_simulations, internal_reactions_l)
    elif not fluxomic:
        gen_simulations = binarize_cols(gen_simulations, reactions_l)

    #--------------------------------------------------------------------------#
    # Objective reaction
    #--------------------------------------------------------------------------#
    if not fluxomic and not cinetic:
        gen_simulations = remove_cols(gen_simulations, [obj])

    #--------------------------------------------------------------------------#
    # Regulatory proteins
    #--------------------------------------------------------------------------#
    if not transcriptomic:
        mn_elements = reactions.union(inputs).union(outputs).union(metabolites)
        reg_proteins = list(set(bn).difference(mn_elements))
        gen_simulations = remove_cols(gen_simulations, reg_proteins)

    return gen_simulations


#NOISE##########################################################################
def add_qualitative_noise(simulations: list, cols: list,
                          error_rate: float = 0.0, seed: int = 0) -> dict:
    """
    """
    random.seed(seed)
    for i in range(len(list(simulations.iterrows()))):
        for c in cols:
            noise = random.random()
            if noise < error_rate:
                simulations[c].iloc[i] = None
    return simulations


def add_quantitative_noise(simulations: pd.DataFrame, cols: list,
                           error_range: float = 0.0, seed: int = 0) -> dict:
    """
    """
    random.seed(seed)
    for i in range(len(list(simulations.iterrows()))):
        for c in cols:
            noise = 2 * error_range * random.random() - error_range
            simulations[c].iloc[i] *= 1 + noise
    return simulations


def sample_random(simulations: dict, n: int, seed: int = 0, first: bool = True, last: bool = True) -> dict:
    """
    """
    random.seed(seed)
    for k, sim in simulations.items():
        index = sorted(sim.index)
        index_to_sample = index.copy()
        if first:
            index_to_sample = index_to_sample[1:]
        if last:
            index_to_sample = index_to_sample[:-1]
        sample_index = random.sample(index_to_sample, n)
        if first:
            sample_index = [index[0]] + sample_index
        if last:
            sample_index = sample_index + [index[-1]]
        sample_index = sorted(sample_index)
        simulations[k] = sim.loc[sample_index]
    return simulations


def sample_freq(simulations: dict, freq: int, noise_range: int = 0, seed: int = 0) -> dict:
    """
    """
    random.seed(seed)
    for k, sim in simulations.items():
        index = sorted(sim.index)
        sample = [index[0]]
        counter = freq + random.randint(-noise_range, noise_range)
        for i in index:
            counter -= 1
            if counter == 0:
                sample.append(index[i])
                counter = freq - (counter) + \
                    random.randint(-noise_range, noise_range)
        if index[-1] not in sample:
            sample.append(index[-1])
        simulations[k] = sim.loc[sorted(sample)]
    return simulations

#AUXILIARY######################################################################


def write_simulations(output_dir: str, simulations: dict) -> None:
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if output_dir[-1] != '/':
        output_dir = output_dir + '/'
    for i, sim in simulations.items():
        sim.to_csv(f'{output_dir}{i}.csv', sep='\t', header=True, na_rep='')


#OTHERS##########################################################################

def main():
    """
    Function: main function executing the project
    Params: None
    Return: None
    """


#RUN############################################################################
if __name__ == '__main__':
    main()
