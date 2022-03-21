#! /usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Kerian Thuillier"
__email__ = "kerian.thuillier@ens-rennes.fr"

#IMPORT#########################################################################
import argparse
import pandas as pd
from pandas.core.series import Series

from .generators import (
    add_qualitative_noise, add_quantitative_noise, generate_data, write_simulations)

#CODE###########################################################################


def parsing():
    """
    Function: parse python arguments
    Params: None
    Return: Argparse object
    """
    parser = \
        argparse.ArgumentParser(
            description='Python3 Script generating realistics data',
            epilog='',
            prog='')

    #****************************************************************************#
    # Input simulations
    #****************************************************************************#
    parser.add_argument('-sim', '--simulations',
                        help='FlexFlux simulations.',
                        required=True,
                        dest='sim',
                        nargs='+',
                        type=str)
    parser.add_argument('-sbml', '--SBML',
                        help='Metabolic network in SBML file.',
                        required=True,
                        dest='sbml',
                        type=str)
    parser.add_argument('-reg', '--reg-SBML',
                        help='Regulatory network in SBML file.',
                        required=True,
                        dest='reg',
                        type=str)

    #****************************************************************************#
    # Output simulations
    #****************************************************************************#
    parser.add_argument('-out', '--output',
                        help='Output generated file.',
                        required=True,
                        dest='out',
                        type=str)

    #****************************************************************************#
    # Data type
    #****************************************************************************#
    parser.add_argument('--fluxomic',
                        help='Generate fluxomic data.',
                        required=False,
                        dest='fluxomic',
                        action='store_true',
                        default=False)

    parser.add_argument('--cinetic',
                        help='Generate cinetic data.',
                        required=False,
                        dest='cinetic',
                        action='store_true',
                        default=False)

    parser.add_argument('--transcriptomic',
                        help='Generate transcriptomic data.',
                        required=False,
                        dest='transcriptomic',
                        action='store_true',
                        default=False)

    #****************************************************************************#
    # Noise
    #****************************************************************************#
    parser.add_argument('--noise-cinetic-range',
                        help='Hide internal reactions.',
                        required=False,
                        dest='noise_cinetic_range',
                        type=float,
                        default=0.0)

    parser.add_argument('--noise-flux-range',
                        help='Add noise to flux.',
                        required=False,
                        dest='noise_flux_range',
                        type=float,
                        default=0.0)

    parser.add_argument('--noise-observation-freq',
                        help='Remove some observations at a given frequence.',
                        required=False,
                        dest='noise_observation_freq',
                        type=float,
                        default=0.0)

    parser.add_argument('--sample-observations-rate',
                        help='Select a subset of the total observations.',
                        required=False,
                        dest='noise_observation_freq',
                        type=float,
                        default=0.0)

    #****************************************************************************#
    # Other options
    #****************************************************************************#
    parser.add_argument('--hide-reactions',
                        help='Hide internal reactions.',
                        required=False,
                        dest='hide_reactions',
                        action='store_true',
                        default=False)

    args = parser.parse_args()

    return args


def main():
    """
    Function: main function executing the project
    Params: None
    Return: None
    """
    args = parsing()
    gen_sims = generate_data(args.sim, args.sbml, args.reg,
                             fluxomic=args.fluxomic,
                             transcriptomic=args.transcriptomic,
                             cinetic=args.cinetic,
                             hide_reactions=args.hide_reactions)

    write_simulations(args.out, gen_sims)

    # gen_sims2 = add_qualitative_noise(gen_sims,  ['RPcl', 'RPb', 'RPh', 'RPO2'], error_rate=0.75, seed=0)

    # gen_sims2 = add_quantitative_noise(gen_sims,  ['RPcl', 'RPb', 'RPh', 'RPO2'], error_range=0.05, seed=0)
    # gen_sims2 = sample(gen_sims, 10, first=True, last=True, seed=0)
    # print(gen_sims2)


#RUN############################################################################
if __name__ == '__main__':
    main()
