###############
# User Inputs #
###############

# Material System Name
Material 	= "Hynes282"

# Which ThermoCalc Database to use?
Thermo_Database = 'TCNI11'

# Solvent Element Name
solvent		= 'NI'

# Solute Element Names List
solutes		= ['CR','CO','MO','TI','AL']

# Concentration of Solute Elements
concentration	= [0.0196, 0.103,0.084,0.02, 0.0154]

# Min Temperature for calculation
T_min 		= 300

# Min No. of time steop
Num_time_steps 	= 20

################
# Code Section #
################

import matplotlib.pyplot as plt
import pandas as pd
from tc_python import *
import shutil

def set_multiple_conditions(solutes,conc):
    temp=["s-c"]
    for i in range(len(solutes)):
        temp.append('X('+solutes[i]+")="+str(conc[i]))
    return (' '.join(temp))

with TCPython() as session:
    system =(session.set_cache_folder(os.path.basename(Material) + "_cache")
                    .select_database_and_elements(Thermo_Database, [solvent] + solutes))
    
    singleP_eq = (system.get_system()
                    .with_single_equilibrium_calculation()
                    .run_poly_command(set_multiple_conditions(solutes,concentration))
                    .disable_global_minimization())
    

    #Calculating the liquidus temperature and stable phases
    Liquidus = singleP_eq.set_phase_to_fixed("LIQUID",1).calculate()
    Tl = Liquidus.get_value_of(ThermodynamicQuantity.temperature())
    
    T_init = (T_min + int(Tl+10))*0.5

    PhaseDiagram = (system.get_system()
                    .with_property_diagram_calculation()
                    .run_poly_command(set_multiple_conditions(solutes,concentration))
                    .set_condition("P", 101325)
                    .set_condition("n", 1)
                    .with_axis(CalculationAxis(ThermodynamicQuantity.temperature())
                        .set_min(T_min)
                        .set_max(int(Tl+10))
                        .with_axis_type(Linear().set_min_nr_of_steps(Num_time_steps)))
                    .set_condition(ThermodynamicQuantity.temperature(), T_init)
                    .enable_global_minimization()
                    .calculate())

    groups = PhaseDiagram.get_values_grouped_by_quantity_of(ThermodynamicQuantity.temperature(), ThermodynamicQuantity.mass_fraction_of_a_phase(ALL_PHASES))

PhaseMap = pd.DataFrame(columns=["Phase", "Min T", "Max T"])

for group in groups.values():
    header = group.label
    plt.plot(group.x, group.y, label= header)
    new_Phase = {
        "Phase": group.label,
        "Min T": min(group.x),
        "Max T": max(group.x)
    }
    PhaseMap = pd.concat([PhaseMap, pd.DataFrame([new_Phase])], ignore_index=True)

#print the liquidus temperature
print(f"Liquidus Temperature : {Tl} K")

# Print the updated DataFrame for available phases and temperature range
print(PhaseMap)

# Amound of phases vs temperature plot
plt.xlabel("Temperature (K)")
plt.axis([T_min, Tl+10, 0, 1.1])
plt.ylabel("Mass fraction of phases")
plt.legend(loc="upper left", fontsize=8)
plt.savefig("Phase_diagram"+Material+".png")
plt.show()

shutil.rmtree(os.path.basename(Material) + "_cache")
