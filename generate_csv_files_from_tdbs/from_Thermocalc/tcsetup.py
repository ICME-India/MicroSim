import os
from tc_python import *


class setup:
    def __init__(self, start, Thermodynamic_database, Kinetic_database, solutes, solvent, composition, comp_choice, Tmax, Tmin, material):
        self.Thermodata = Thermodynamic_database
        self.Kineticdata = Kinetic_database
        self.solutes = solutes
        self.solvent = solvent
        self.composition = composition
        self.Tmax = Tmax
        self.Tmin = Tmin
        self.Fname = material
        self.start = start
        self.comp_choice = comp_choice

    def _set_multiple_conditions(self):
        c = "W(" if self.comp_choice == 'mass' else "X("
        temp=["s-c"]
        for i in range(len(self.solutes)):
            temp.append(c+self.solutes[i]+")="+str(self.composition[i]))
        return (' '.join(temp))

    def Database_selection(self):
        calculation_setup= (
        self.start
            .set_cache_folder(os.path.basename(self.Fname[:-4]) + "_cache")
            .select_thermodynamic_and_kinetic_databases_with_elements(self.Thermodata, self.Kineticdata, self.solutes+[self.solvent])
            )
        return calculation_setup
    
    def Eq_composition(self, calculation_setup, num_steps):
        condition = setup._set_multiple_conditions(self)
        Tinit = self.Tmin + int(0.5*(self.Tmax - self.Tmin))
        propertyDia=(calculation_setup
                .get_system()
                .with_property_diagram_calculation()
                .run_poly_command(condition)          
                .with_axis(CalculationAxis(ThermodynamicQuantity.temperature()).
                              set_min(self.Tmin).
                              set_max(self.Tmax).
                              with_axis_type(Linear().
                                             set_min_nr_of_steps(num_steps)))
                .set_condition(ThermodynamicQuantity.temperature(), Tinit)
                .enable_global_minimization()
                .calculate())
        return propertyDia
    
    def Hessian_setup(self, calculation_setup, phase):
        Hessian=(calculation_setup
                .without_default_phases()
                .select_phase(phase)
                .get_system()
                .with_single_equilibrium_calculation())
        return Hessian

    def Hessian_calculation(self, Hessian_setup, temperature):
        condition = setup._set_multiple_conditions(self)
        l=(Hessian_setup
            .set_condition("T", temperature)
            .run_poly_command(condition)
            .disable_global_minimization()
            .calculate())
        return l
    
    def order(self, phase, some_value):
        order = ['Temperature'] if some_value == 1 else []
        for i in range(len(self.solutes)):
            order.append(self.solutes[i] + phase)
        return order
