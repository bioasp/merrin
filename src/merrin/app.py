# ==============================================================================
# Imports
# ==============================================================================
from __future__ import annotations
from os import path, makedirs
from argparse import Namespace, ArgumentParser

from merrin.learn import MerrinLearner
from merrin.datastructure import Observation, MetabolicNetwork


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
    # ~ PKN
    pkn: list[tuple[str, int, str]] = []
    with open(args.pkn, 'r', encoding='utf-8') as file:
        for line in file.readlines():
            line = line.strip()
            u, s, v = line.split('\t')
            pkn.append((u, int(s), v))
    # ~ Observations
    observations: list[Observation] = Observation.load_json(args.obs)
    # ~ MN
    mn: MetabolicNetwork = MetabolicNetwork.read_sbml(args.sbml)
    # --------------------------------------------------------------------------
    # Initialise Merrin Learn object
    # --------------------------------------------------------------------------
    learner: MerrinLearner = MerrinLearner()
    learner.load_instance(mn, args.obj, pkn, observations)
    if args.projection == 'network':
        learner.learn(
            nbsol=0, display=True, lp_solver=args.lpsolver,
            max_error=0.3, max_gap=10, timelimit=args.timelimit,
            subsetmin=args.optimisation == 'subsetmin'
        )
    elif args.projection == 'node':
        learner.learn_per_node(
            nbsol=0, display=True, lp_solver=args.lpsolver,
            max_error=0.3, max_gap=10, timelimit=args.timelimit,
            subsetmin=args.optimisation == 'subsetmin'
        )


# ==============================================================================
# Main
# ==============================================================================
if __name__ == '__main__':
    main()
