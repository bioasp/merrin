# ==============================================================================
# Imports
# ==============================================================================
from __future__ import annotations
from os import path, makedirs
from argparse import Namespace, ArgumentParser
from datetime import datetime
from pandas import DataFrame

from merrin.learn import MerrinLearner
from merrin.datastructure import Observation


# ==============================================================================
# Functions
# ==============================================================================
def parse() -> Namespace:
    parser: ArgumentParser = ArgumentParser(
        description='Inferring regulatory rules from metabolic traces.',
        epilog='',
        prog=''
    )
    # --------------------------------------------------------------------------
    # Required arguments
    # --------------------------------------------------------------------------
    parser.add_argument(
        '-sbml', '--SBML',
        help='Metabolic network in SBML file format.',
        required=True,
        dest='sbml',
        type=str
    )
    parser.add_argument(
        '-pkn', '--PKN',
        help='Prior Knowledge Network (see format in README).',
        required=True,
        dest='pkn',
        type=str
    )
    parser.add_argument(
        '-obj', '--objective-reaction',
        help='Objective reaction.',
        required=True,
        dest='obj',
        type=str
    )
    parser.add_argument(
        '-obs', '--observations',
        help='JSON file describing the input timeseries' +
        ' (see format in README).',
        required=True,
        dest='obs',
        type=str
    )
    # --------------------------------------------------------------------------
    # Optional arguments
    # --------------------------------------------------------------------------
    parser.add_argument(
        '-out', '--output-file',
        help='Output CSV file',
        required=False,
        default=None,
        dest='output',
        type=str
    )
    parser.add_argument(
        '--lpsolver',
        help='Linear solver to use (default: glpk)',
        required=False,
        default='glpk',
        dest='lpsolver',
        type=str,
        choices=['glpk', 'gurobi']
    )
    parser.add_argument(
        '--timelimit',
        help='Timelimit for each resolution, -1 if none (default: -1)',
        required=False,
        default=-1,
        dest='timelimit',
        type=int
    )
    parser.add_argument(
        '--optimisation',
        help='Select optimisation mode:' +
             ' all networks or subset minimal ones only (default: subsetmin)',
        required=False,
        default='subsetmin',
        dest='optimisation',
        type=str,
        choices=('all', 'subsetmin')
    )
    parser.add_argument(
        '--projection',
        help='Select project mode (default: network):' +
             '\nnode: enumerate the candidate rules for each node'
             '\nnetwork: enumerate all the rules of each network',
        required=False,
        default='network',
        dest='projection',
        type=str,
        choices=('network', 'node')
    )
    return parser.parse_args()


def create_folder(folder):
    if folder == '' or path.exists(folder):
        return
    makedirs(folder)


def main() -> None:
    args: Namespace = parse()
    # --------------------------------------------------------------------------
    # Load inputs
    # --------------------------------------------------------------------------
    pkn: list[tuple[str, int, str]] = []
    with open(args.pkn, 'r', encoding='utf-8') as file:
        for line in file.readlines():
            line = line.strip()
            u, s, v = line.split('\t')
            pkn.append((u, int(s), v))
    observations: list[Observation] = Observation.load_json(args.obs)
    # --------------------------------------------------------------------------
    # Initialise Merrin Learn object
    # --------------------------------------------------------------------------
    learner: MerrinLearner = MerrinLearner(
        args.sbml,
        args.obj,
        pkn
    )
    # ~ Set optimisation mode
    if args.optimisation == 'all':
        learner.set_optimisation(MerrinLearner.Optimisation.ALL)
    elif args.optimisation == 'subsetmin':
        learner.set_optimisation(MerrinLearner.Optimisation.SUBSETMIN)
    # ~ Set projection mode
    if args.projection == 'network':
        learner.set_projection(MerrinLearner.Projection.NETWORK)
    elif args.projection == 'node':
        learner.set_projection(MerrinLearner.Projection.NODE)
    results: DataFrame = learner.learn(
        observations,
        timelimit=args.timelimit,
        lpsolver=args.lpsolver
    )
    # --------------------------------------------------------------------------
    # Export Output CSV file
    # --------------------------------------------------------------------------
    output_csv: str | None = args.output
    if output_csv is None:
        # ~ Timestamp
        now: datetime = datetime.now()
        timestamp: str = now.strftime("%Y%m%d-%H%M%S")
        output_csv = \
            f'merrin-{args.optimisation}-{args.projection}-{timestamp}.csv'
    create_folder(path.dirname(output_csv))
    results.to_csv(
        output_csv,
        sep=',',
        index=False
    )


# ==============================================================================
#
# ==============================================================================
if __name__ == '__main__':
    main()
