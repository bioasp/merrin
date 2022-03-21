#! /usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Kerian Thuillier"
__email__ = "kerian.thuillier@ens-rennes.fr"

#IMPORT#########################################################################

import argparse
import pandas as pd

import merrin.Loader as loader
import merrin.ASPApplication as caspoflux
import merrin.DataProcessing as DataProcessing
from .Parameters import FLUXOMIC_DATA, PROTEOMIC_DATA
from .ASPApplication import ProjectionMode, ResolutionMode

#CODE###########################################################################


def parsing():
    """
    Function: parse python arguments
    Params: None
    Return: Argparse object
    """
    parser = \
        argparse.ArgumentParser(
            description='Inferring regulatory rules from metabolic traces.',
            epilog='',
            prog='')

    #**************************************************************************#
    # REQUIRED
    #**************************************************************************#
    parser.add_argument('-fluxomic', '--fluxomic',
                        help='Input fluxomic observations.',
                        required=False,
                        dest='fluxomic',
                        nargs='+',
                        default=[],
                        type=str)

    parser.add_argument('-proteomic', '--proteomic',
                        help='Input proteomic observations.',
                        required=False,
                        dest='proteomic',
                        nargs='+',
                        default=[],
                        type=str)

    parser.add_argument('-sbml', '--SBML',
                        help='Metabolic network in SBML file.',
                        required=True,
                        dest='sbml',
                        type=str)
    parser.add_argument('-pkn', '--PKN',
                        help='Prior Knowledge Network.',
                        required=True,
                        dest='pkn',
                        type=str)
    parser.add_argument('-obj', '--objective-reaction',
                        help='Objective reaction.',
                        required=True,
                        dest='obj',
                        type=str)

    #**************************************************************************#
    # Output simulations
    #**************************************************************************#
    parser.add_argument('--output',
                        help='Output generated file.',
                        required=False,
                        dest='out',
                        default='./out',
                        type=str)

    #**************************************************************************#
    # OPTIONAL
    #**************************************************************************#
    parser.add_argument('--projection',
                        help='Projection mode: Network, Component or Node.',
                        required=False,
                        dest='projection',
                        default='Network',
                        type=str)

    args = parser.parse_args()

    # TODO: gérer l'erreur si le mode de projection demandé n'existe pas.

    if args.projection == 'Network':
        args.projection = ProjectionMode.NETWORK
    elif args.projection == 'Component':
        args.projection = ProjectionMode.COMPONENT
    elif args.projection == 'Node':
        args.projection = ProjectionMode.NODE
    else:
        assert(False)  # FIXME

    assert(len(args.proteomic) != 0 or len(args.fluxomic) != 0)  # FIXME

    return args


def main():
    """
    Function: main function executing the project
    Params: None
    Return: None
    """
    args = parsing()

    #--------------------------------------------------------------------------#
    # DATA PROCESSING
    #--------------------------------------------------------------------------#
    mn = loader.load_sbml(args.sbml)
    inputs = sorted(mn.get_inputs())
    pkn = loader.load_pkn(args.pkn, inputs)
    obj = args.obj
    flux_obs = loader.load_simulations(
        args.fluxomic,
        tag=FLUXOMIC_DATA)
    if len(args.fluxomic) != 0:
        prot_obs = loader.load_simulations(
            args.proteomic,
            k=1+max(flux_obs.keys()),
            tag=PROTEOMIC_DATA)
    else:
        prot_obs = loader.load_simulations(
            args.proteomic,
            tag=PROTEOMIC_DATA)

    observations = flux_obs | prot_obs

    processed_observations = DataProcessing.preprocess_simulations(
        observations)

    processed_observations.to_csv('temp.csv')

    print('Data are processed...')

    #--------------------------------------------------------------------------#
    # MODEL BUILDING
    #--------------------------------------------------------------------------#

    model = caspoflux.Model(mn, pkn, obj, n=0, thread=1, statistics=True,
                            resolution_mode=ResolutionMode.SUBSET_MINIMAL)
    model.build(processed_observations)

    with open('temp.lp', 'w') as file:
        file.write(model.export_model())

    print('Model is built...')

    #--------------------------------------------------------------------------#
    # MODEL SOLVING
    #--------------------------------------------------------------------------#

    if args.projection == ProjectionMode.NETWORK:
        status = model.solve()
        results = model.get_results()
        print(f'Model is solved... (status: {status})')
        if len(results) == 0:
            print('UNSATISFIABLE')
        else:
            print('SATISFIABLE')
            df = pd.DataFrame(results)
            df = df.fillna('')
            print('-'*10)
            print(df)
            print('-'*10)

    elif args.projection == ProjectionMode.COMPONENT:
        results = model.solve_projection(ProjectionMode.COMPONENT)
        if len(results) == 0:
            print('UNSATISFIABLE')
        else:
            print('SATISFIABLE')
            print('-'*5, end='')
            for cc in results:
                print('-'*5)
                df = pd.DataFrame(results[cc])
                df = df.fillna('')
                print(df)
            print('-'*10)

    elif args.projection == ProjectionMode.NODE:
        results = model.solve_projection(ProjectionMode.NODE)
        if len(results) == 0:
            print('UNSATISFIABLE')
        else:
            print('SATISFIABLE')
            print('-'*10)
            df = pd.DataFrame(dict((k, pd.Series(v, dtype='str'))
                                   for k, v in results.items()))
            df = df.fillna('')
            print(df)
            print('-'*10)


#RUN############################################################################
if __name__ == '__main__':
    main()
