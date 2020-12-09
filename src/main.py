import libsbml
import sys


def obj2str(obj):
    args = ",".join([f"{k}={v}" for k, v in obj.__dict__.items()])
    return f"{obj.__class__.__name__}({args})"


class Specie:
    def __init__(self, nome, compartment, ivalue, constant):
        self.nome = nome
        self.compartment = compartment
        self.ivalue = ivalue
        self.constant = constant
        self.involved_as_reactant = []
        self.involved_as_product = []
        self.involved_as_modifier = []

    def __str__(self):
        return obj2str(self)

    def add_reaction_as_reactant(self, reaction, stoichiometry_value):
        self.involved_as_reactant.append((reaction, stoichiometry_value))
    
    def add_reaction_as_product(self, reaction, stoichiometry_value):
        self.involved_as_product.append((reaction, stoichiometry_value))
    
    def add_reaction_as_modifier(self, reaction):
        self.involved_as_modifier.append(reaction)
    

class Compartment:
    def __init__(self, nome, size, volume, unit):
        self.nome = nome
        self.size = size
        self.volume = volume
        self.unit = unit

    def __str__(self):
        return obj2str(self)


class Parameter:
    def __init__(self, nome, value, constant):
        self.nome = nome
        self.value = value
        self.constant = constant

    def __str__(self):
        return obj2str(self)


class Reaction:
    def __init__(self, nome, second_name, reactant, product, modifier, local_parameters, math_formula):
        self.nome = nome
        self.second_name = second_name
        self.reactants = reactant
        self.products = product
        self.modifiers = modifier
        self.local_parameters = local_parameters
        self.math_formula = math_formula

    def __str__(self):
        return obj2str(self)


class Rule:
    def __init__(self, lhs, rhs):
        self.lhs = lhs
        self.rhs = rhs
    
    def __str__(self):
        return obj2str(self)


class AssignmentRule(Rule):
    def __init__(self, lhs, rhs):
        super().__init__(lhs, rhs)

    
class RateRule(Rule):
    def __init__(self, lhs, rhs):
        super().__init__(f"der({lhs})", rhs)


class SBMLModel:
    def __init__(self, name, compartmnents, species, parameters, assignment_rules, reactions):
        self.name = name
        self.compartments = compartmnents
        self.species = species
        self.parameters = parameters
        self.assignment_rules = assignment_rules
        self.reactions = reactions
        self.create_rate_rule()

    def create_sum_from_reactant(self, specie_name, specie_obj):
        formula_list = []
        for reaction_id, stoichiometry_value in specie_obj.involved_as_reactant:
            formula_list.append(f"({stoichiometry_value} * {self.reactions[reaction_id].math_formula})")
        return " - ".join(formula_list)

    def create_sum_from_products(self, specie_name, specie_obj):
        formula_list = []
        for reaction_id, stoichiometry_value in specie_obj.involved_as_product:
            formula_list.append(f"({stoichiometry_value} * {self.reactions[reaction_id].math_formula})")
        return " - ".join(formula_list)

    def create_rate_rule(self):
        self.rate_rules_dict = dict()
        for specie_id, specie_obj in self.species.items():
            rate_rule = RateRule(specie_id, "0.0")
            if not specie_obj.constant:
                reactant_formula = self.create_sum_from_reactant(specie_id, specie_obj)
                product_formula = self.create_sum_from_products(specie_id, specie_obj)
                rate_rule = RateRule(specie_id, f"{product_formula} - {reactant_formula}")
            self.rate_rules_dict[specie_id] = rate_rule

    def __str__(self):
        printable = ""
        for v in self.__dict__.values():
            if isinstance(v, str):
                printable += f"Model Name: {v}\n"
            if isinstance(v, dict):
                printable += "\n".join([x.__str__() for x in v.values()])
            printable += "\n\n"
        return printable


class SBMLTranslator:
    def __init__(self, *args, **kargs):
        pass


class SBMLExtrapolator:
    def __init__(self, sbmlfile):
        self.model = libsbml.readSBMLFromFile(sbmlfile).getModel()
        self.extraploate()

    def extraploate(self):
        self.getmodelname()
        self.getcompartments()
        self.getspecies()
        self.getparameters()
        self.getrules()
        self.getreactions()
    
    def getmodelname(self):
        self.nome = self.model.getName()

    def getcompartments(self):
        """ Prendiamo tutti i compartments del modello """
        self.comp_dict = dict()
        for comp in self.model.getListOfCompartments():
            self.comp_dict[comp.getId()] = Compartment(
                comp.getId(), # Nome
                comp.getSize(), # Taglia
                comp.getVolume(), # Volume
                comp.getUnits() # Unità di misura
            )

    def getspecies(self):
        """ Prendiamo tutte le Specie con i relativi compartment """
        self.species_dict = dict() 
        for sp in self.model.getListOfSpecies():
            self.species_dict[sp.getId()] = Specie(
                sp.getId(), # Nome
                self.comp_dict[sp.getCompartment()], # Oggetto Compartment
                sp.getInitialConcentration(), # Concentrazione iniziale
                sp.getConstant() # Se è costante oppure no
            )

    def getparameters(self):
        """ Prendiamo tutti i parametri del modello in listOfParameters """
        self.parameter_dict = dict()
        for param in self.model.getListOfParameters():
            self.parameter_dict[param.getId()] = Parameter(
                param.getId(), # Nome
                param.getValue(), # Valore
                param.getConstant() # Se è costante oppure no
            )
           
    def getrules(self):
        self._assignment_rule()

    def _assignment_rule(self):
        self.assignment_dict = dict()
        for rule in self.model.getListOfRules():
            if isinstance(rule, libsbml.AssignmentRule):
                self.assignment_dict[rule.getVariable()] = AssignmentRule(
                    rule.getVariable(), rule.getFormula())

    def getreactions(self):
        self.reaction_dict = dict()
        for reaction in self.model.getListOfReactions():
            reaction_name = reaction.getId()
            second_reaction_name = reaction.getName()
            kinetic_law = reaction.getKineticLaw().getFormula()
            reactants = self.get_and_set(reaction_name, reaction.getListOfReactants(), "r")
            products = self.get_and_set(reaction_name, reaction.getListOfProducts(), "p")
            modifiers = self.get_and_set_modifier(reaction, reaction_name)
            parameters = []
            list_of_parameters = list(filter(lambda x: isinstance(x, libsbml.Parameter), reaction.getListOfAllElements()))
            # Devo vedere se ci sono parametri locali alla reazione
            for param in list_of_parameters:
                self.parameter_dict[param.getId()] = Parameter(
                    param.getId(), # Nome
                    param.getValue(), # Valore
                    param.getConstant() # Se è costante oppure no
                )
                parameters.append(param.getId())
            self.reaction_dict[reaction_name] = Reaction(
                reaction_name, second_reaction_name, reactants, 
                products, modifiers, parameters, kinetic_law
            )
        
    def get_and_set(self, reaction_name, lista, category):
        # Prendo tutti i reagenti con i loro valori stechiometrici
        elements = []
        for element in lista:
            # Ottengo un oggetto SpeciesReference
            # Devo ottenere il valore stechiometrico
            stoichiometry_value = element.getStoichiometry()
            if stoichiometry_value == 1.0 and element.isSetStoichiometryMath():
                tmp_kinetic = libsbml.KineticLaw(2, 4)
                tmp_kinetic.setMath(element.getStoichiometryMath().getMath())
                stoichiometry_value = tmp_kinetic.getFormula()
            if category == "r":
                self.species_dict[element.getSpecies()].add_reaction_as_reactant(
                    reaction_name, stoichiometry_value
                )
            elif category == "p":
                self.species_dict[element.getSpecies()].add_reaction_as_product(
                    reaction_name, stoichiometry_value
                )
            elements.append(element.getSpecies())
        return elements
    
    def get_and_set_modifier(self, libsbml_reaction_obj, reaction_name):
        modifiers = []
        for modif in libsbml_reaction_obj.getListOfModifiers():
            modifiers.append(modif.getSpecies())
            self.species_dict[modif.getSpecies()].add_reaction_as_modifier(reaction_name)
        return modifiers


if __name__ == "__main__":
    try:
        modelname = sys.argv[1]
        sbmlext = SBMLExtrapolator(modelname)
        sbmlmodel = SBMLModel(sbmlext.nome, 
                              sbmlext.comp_dict, 
                              sbmlext.species_dict,
                              sbmlext.parameter_dict,
                              sbmlext.assignment_dict,
                              sbmlext.reaction_dict
                             )
        print(sbmlmodel)
    except Exception as e:
        print(e)
        print("Devi inserire come argomenti il path del modello")