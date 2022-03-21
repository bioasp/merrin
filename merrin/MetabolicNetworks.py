#! /usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Kerian Thuillier"
__email__ = "kerian.thuillier@ens-rennes.fr"

#IMPORT#########################################################################

import libsbml
from .Parameters import SIMULATION_TIMESTEP


#CLASS##########################################################################


class MetabolicNetwork:

    def __init__(self):
        """
        Function: 
        Params: 
        Return: 
        """
        self.__reactions = set()
        self.__bounds = {}
        self.__inputs = set()
        self.__outputs = set()
        self.__metabolites = set()
        self.__stoichiometry = {}

    #INITIALISATION#############################################################

    def init_sbml(self, sbml):
        """
        Function: 
        Params: 
        Return: 
        """

        flux_bounds_parameters = {}
        for parameter in sbml.getListOfParameters():
            flux_bounds_parameters[parameter.id] = parameter.value

        reactants = set()
        products = set()
        self.__reactions = set()
        self.__stoichiometry = {}
        for reaction in sbml.getListOfReactions():
            name = reaction.getId()
            self.__reactions.add(name)
            for a in reaction.getListOfReactants():
                coeff = a.getStoichiometry()
                species = a.getSpecies()
                reactants.add(species)
                self.__stoichiometry[(name, species)] = -coeff
            for a in reaction.getListOfProducts():
                coeff = a.getStoichiometry()
                species = a.getSpecies()
                products.add(species)
                self.__stoichiometry[(name, species)] = coeff

            xml_attributes = libsbml.XMLNode.getAttributes(
                reaction.toXMLNode())
            lower_bound = flux_bounds_parameters[xml_attributes.getValue(
                'lowerFluxBound')]
            upper_bound = flux_bounds_parameters[xml_attributes.getValue(
                'upperFluxBound')]
            self.__bounds[name] = (lower_bound, upper_bound)

        self.__inputs = reactants.difference(products)
        self.__outputs = products.difference(reactants)
        output_metabolites = self.__inputs.union(self.__outputs)
        self.__metabolites = products.union(
            reactants).difference(output_metabolites)

    #GETTER#####################################################################

    def get_reactions(self):
        return self.__reactions

    def get_metabolites(self):
        return self.__metabolites

    def get_inputs(self):
        return self.__inputs

    def get_outputs(self):
        return self.__outputs

    def get_stoichiometry(self, reaction, metabolite):
        assert(reaction in self.__reactions)
        assert(metabolite in self.__metabolites
               or metabolite in self.__inputs
               or metabolite in self.__outputs)
        if (reaction, metabolite) in self.__stoichiometry:
            return self.__stoichiometry[(reaction, metabolite)]
        return 0

    def get_stoichiometric_coefficients(self):
        return self.__stoichiometry

    def get_reaction_bounds(self, reaction):
        assert(reaction in self.__reactions)
        if reaction in self.__bounds:
            return self.__bounds[reaction]
        return (None, None)

    def get_input_reactions(self):
        return set((r, m) for r, m in self.__stoichiometry if m in self.__inputs)

    #OTHER######################################################################

    def get_transport_reaction_bounds(self, reaction, metabolite, C, biomass):
        assert(reaction in self.__reactions)

        upper_bound = C / (biomass * SIMULATION_TIMESTEP)

        upper_bound = min(upper_bound, self.__bounds[reaction][1])

        return (self.__bounds[reaction][0], upper_bound)
