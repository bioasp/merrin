# ==============================================================================
# Imports
# ==============================================================================
from __future__ import annotations

from libsbml import (  # type: ignore
    Model,
    SBMLDocument,
    SBMLReader,
    Reaction,
    FbcModelPlugin,
    FbcReactionPlugin,
    XMLNode,
    GeneProductAssociation
)


# ==============================================================================
# Metabolic Network
# ==============================================================================
# > Data structure use to model Metabolic Network
# > Associated with all the basic functions needed for Flux Balance Analysis
# ==============================================================================
class MetabolicNetwork:

    # ==========================================================================
    # Initialization
    # ==========================================================================
    def __init__(self: MetabolicNetwork):
        self.irreversible: bool = True

        # ----------------------------------------------------------------------
        # Metabolism
        # ----------------------------------------------------------------------
        self.__metabolites: set[str] = set()
        self.__reactions: set[str] = set()
        self.__exchanges: dict[str, str | None] = {}
        self.__coeff: dict[tuple[str, str], float] = {}
        self.__bounds: dict[str, tuple[float, float]] = {}
        self.__reversible_reactions: set[str] = set()

        # ----------------------------------------------------------------------
        # Genes
        # ----------------------------------------------------------------------
        self.__genes: set[str] = set()
        self.__genes_association: dict[str, str] = {}

        # ----------------------------------------------------------------------
        # Reversible reaction cache
        # ----------------------------------------------------------------------
        self.__reversible_reactions_mapping: dict[str, tuple[str, str]] = {}

    # ==========================================================================
    # Getters
    # ==========================================================================
    def reactions(self: MetabolicNetwork, only_exchange: bool = False) \
            -> set[str]:
        if only_exchange:
            return set(self.__exchanges.values())  # type: ignore
        return self.__reactions

    def previously_reversible_reactions(self: MetabolicNetwork) \
            -> dict[str, tuple[str, str]]:
        return self.__reversible_reactions_mapping

    def metabolites(self: MetabolicNetwork, external: bool = True,
                    internal: bool = True) -> set[str]:
        if external and internal:
            return self.__metabolites | set(self.__exchanges.keys())
        if external:
            return set(self.__exchanges.keys())
        if internal:
            return self.__metabolites
        return set()

    def exchange(self: MetabolicNetwork, metabolite: str) -> str | None:
        assert metabolite in self.__exchanges
        return self.__exchanges[metabolite]

    def bound(self: MetabolicNetwork, r: str) -> tuple[float, float]:
        assert r in self.__bounds
        return self.__bounds[r]

    def bounds(self: MetabolicNetwork) -> dict[str, tuple[float, float]]:
        return self.__bounds

    def stoichiometry(self: MetabolicNetwork) -> dict[tuple[str, str], float]:
        return self.__coeff

    def is_reversible(self: MetabolicNetwork, r: str) -> bool:
        return r in self.__reversible_reactions

    def genes(self: MetabolicNetwork) -> set[str]:
        return self.__genes

    # ==========================================================================
    # Setters
    # ==========================================================================
    def set_bound(self: MetabolicNetwork, r: str,
                  lb: float, ub: float) -> None:
        for r_ in self.__reversible_reactions_mapping.get(r, [r]):
            assert r_ in self.__bounds
            self.__bounds[r_] = (lb, ub)

    # ==========================================================================
    # Auxiliary functions
    # ==========================================================================
    def to_irreversible(self: MetabolicNetwork) -> MetabolicNetwork:
        def __irreversible_renaming(r: str) -> tuple[str, str]:
            return f'{r}_forward', f'{r}_reverse'

        if self.irreversible:
            return self
        # ----------------------------------------------------------------------
        # Update Exchange reactions
        # ----------------------------------------------------------------------
        for m, r in self.__exchanges.items():
            assert r is not None
            if not self.is_reversible(r):
                continue
            rf, rr = __irreversible_renaming(r)
            self.__reversible_reactions_mapping[r] = (rf, rr)
            if self.__coeff[(m, r)] < 0:  # case: m is a reactant
                self.__exchanges[m] = rf
            else:  # case: m is a product
                self.__exchanges[m] = rr
        # ----------------------------------------------------------------------
        # Update Reactions
        # ----------------------------------------------------------------------
        for r in self.__reactions:
            if not self.is_reversible(r):
                continue
            rf, rr = __irreversible_renaming(r)
            self.__reversible_reactions_mapping[r] = (rf, rr)
            # ------------------------------------------------------------------
            # Update reaction set
            # ------------------------------------------------------------------
            self.__reactions.remove(r)
            self.__reactions.add(rf)
            self.__reactions.add(rr)
            # ------------------------------------------------------------------
            # Update bounds
            # ------------------------------------------------------------------
            del self.__bounds[r]
            self.__bounds[rf] = (0, self.__bounds[r][1])
            self.__bounds[rr] = (0, -self.__bounds[r][0])
            # ------------------------------------------------------------------
            # Update stoichiometry
            # ------------------------------------------------------------------
            for m, r_ in self.__coeff:
                if r != r_:
                    continue
                del self.__coeff[(m, r)]
                self.__coeff[(m, rf)] = self.__coeff[(m, r)]
                self.__coeff[(m, rr)] = -self.__coeff[(m, r)]
            # ------------------------------------------------------------------
            # Update gene association
            # ------------------------------------------------------------------
            del self.__genes_association[r]
            self.__genes_association[rf] = self.__genes_association[r]
            self.__genes_association[rr] = self.__genes_association[r]
        # ----------------------------------------------------------------------
        # Return
        # ----------------------------------------------------------------------
        self.irreversible = True
        return self

    # ==========================================================================
    # Parsing
    # ==========================================================================
    @classmethod
    def __parse_gene_formula_rec(cls, xml_node: XMLNode) -> str:
        name: str = xml_node.getName()
        if name == 'geneProductRef':
            gene_product_name: str | None = None
            for i in range(xml_node.getAttributesLength()):
                if xml_node.getAttrName(i) == 'geneProduct':
                    gene_product_name = xml_node.getAttrValue(i)
                    break
            if gene_product_name is None:
                print('No gene product found...')
                assert False
            return f'{gene_product_name}'
        if xml_node.getNumChildren() == 0:
            print('Structural error...')
            assert False
        children: list[str] = [
            cls.__parse_gene_formula_rec(xml_node.getChild(i))
            for i in range(xml_node.getNumChildren())
        ]
        if len(children) == 1:
            return children[0]
        if name == 'and':
            return '(' + ' & '.join(children) + ')'
        if name == 'or':
            return '(' + ' | '.join(children) + ')'
        print('Unknown fbc operator')
        assert False

    @classmethod
    def __parse_gene_formula(cls, r: Reaction) -> str:
        r_fbc: FbcReactionPlugin | None = r.getPlugin('fbc')
        if r_fbc is None:
            return ''
        gpa: GeneProductAssociation = r_fbc.getGeneProductAssociation()
        if gpa is None:
            return ''
        xml_node: XMLNode = gpa.toXMLNode()
        if xml_node.getNumChildren() != 1:
            print('I do not know what happen here...', f'{r.getName()}')
            assert False
        return cls.__parse_gene_formula_rec(xml_node.getChild(0))

    @classmethod
    def read_sbml(cls, sbml_file: str) -> MetabolicNetwork:
        name: str

        # Object initialization
        mn: MetabolicNetwork = cls()

        # Open SBML file
        sbmld: SBMLDocument = SBMLReader().readSBML(sbml_file)  # type: ignore
        sbmlm: Model = sbmld.getModel()
        sbml_fbc: FbcModelPlugin = sbmlm.getPlugin('fbc')

        # Param
        param: dict[str, float] = {}
        for parameter in sbmlm.getListOfParameters():
            param[parameter.getId()] = parameter.getValue()

        # Genes
        if sbml_fbc is not None:
            for gene_product in sbml_fbc.getListOfGeneProducts():
                mn.__genes.add(gene_product.getId())

        # Species
        for species in sbmlm.getListOfSpecies():
            name = species.getId()
            if not species.getBoundaryCondition():
                mn.__metabolites.add(name)
            else:
                mn.__exchanges[name] = None

        # Reactions
        for reaction in sbmlm.getListOfReactions():
            name = reaction.getId()
            mn.__reactions.add(name)
            stoechiometry: float
            # Reactants
            for reactant in reaction.getListOfReactants():
                stoechiometry = float(reactant.getStoichiometry())
                reactant_name: str = reactant.getSpecies()
                mn.__coeff[(reactant_name, name)] = -stoechiometry
                if reactant_name in mn.__exchanges:
                    mn.__exchanges[reactant_name] = name

            # Products
            for product in reaction.getListOfProducts():
                stoechiometry = float(product.getStoichiometry())
                product_name: str = product.getSpecies()
                mn.__coeff[(product_name, name)] = stoechiometry
                if product_name in mn.__exchanges:
                    mn.__exchanges[product_name] = name

            # Bounds
            reaction_fbc: FbcReactionPlugin | None = reaction.getPlugin('fbc')
            low_bound: float
            up_bound: float
            if reaction_fbc is not None:
                low_bound = param[reaction_fbc.getLowerFluxBound()]
                up_bound = param[reaction_fbc.getUpperFluxBound()]
            else:
                low_bound = 0
                up_bound = float('inf')
                if reaction.getReversible():
                    low_bound = float('-inf')
            mn.__bounds[name] = (low_bound, up_bound)

            # Is Reversible
            if reaction.getReversible():
                mn.__reversible_reactions.add(name)
                mn.irreversible = False

            # Genes formula
            mn.__genes_association[name] = \
                MetabolicNetwork.__parse_gene_formula(reaction)

        return mn
