from tc_python import *
import pandas as pd
import json
import csv
import shutil
import logging
from tcsetup import *

logging.basicConfig(filename='tdbs_microsim.log', level=logging.INFO)
#, format='%(asctime)s %(message)s', datefmt='%m/%d/%y %I:%M:%S %p')
logging.info('The Calculation starts for composition and hessian files')


#open and read json file and store the variables
with open(sys.argv[1], 'r') as myfile:
    data=myfile.read()

# parse file
obj = json.loads(data)

Thermodatabase	= obj['Thermodynamic_Database']
Kineticdatabase	= obj['Kintic_Database']
solvent		= obj['Base_element']
solutes		= obj['Solutes']
conc		= obj['alloy_composition']
comp_choice	= obj["Composition_unit"]
T_min		= obj['Tmin']
T_max		= obj['Tmax']
matrix_pahse	= obj['Matrix_phase']
Active_phases	= obj['Solid_Phases']
material	= obj['Material']

logging.info('The variables are imported from JSON script')

num_steps = int(T_max-T_min)

current_dir = os.getcwd()

# Create a new directory because it does not exist
isExist = os.path.exists('tdbs_encrypted_'+material)
if not isExist:
   os.mkdir('tdbs_encrypted'+material)
os.chdir('tdbs_encrypted'+material)

logging.info('Directory is created : tdbs_encrypted_'+material)

#Strart the thermocalc Calculation
with TCPython() as start:
    # create equilibrium calculation object
    S1 = setup(start, Thermodatabase, Kineticdatabase, solutes, solvent, conc, comp_choice, T_max, T_min, material)
    calculation_setup = S1.Database_selection()
    logging.info('System is generated with Thermodynamic database and Kinetic database')

    #Set condition for equilobrium thermodynamic calculation
    comp_system = S1.Eq_composition(calculation_setup, num_steps)
    logging.info('The instace for equilibrium composition is created from '+str(T_min)+' to '+str(T_max)+' with_property_diagram_calculation')

    for i in range(len(Active_phases)):
        temp = comp_system.get_values_grouped_by_quantity_of(ThermodynamicQuantity.temperature(), "X("+Active_phases[i]+",*)")
        logging.info('The equilibrium composition calculation is performed for '+Active_phases[i])
        if i == 0:
            df1= pd.DataFrame({'Temperature': (temp[solutes[0]+' in '+Active_phases[i]]).x})
        
        dft = pd.DataFrame()
        for group in temp.values():
            colHeader = group.label[0: 2:].replace(" ", "")+ "@" + Active_phases[i]
            dft[colHeader] = group.y
    
        dft = dft.drop(solvent+'@'+Active_phases[i], axis=1)
        dft = dft[S1.order('@'+Active_phases[i], 0)]
        df1 = pd.concat([df1, dft], axis=1)

    temp = comp_system.get_values_grouped_by_quantity_of(ThermodynamicQuantity.temperature(), "X("+matrix_pahse+",*)")
    dft = pd.DataFrame()
    for group1 in temp.values():
        colHeader = group1.label[0: 2:].replace(" ", "") +"@"+matrix_pahse
        dft[colHeader] = group1.y
    
    dft = dft.drop(solvent+'@'+matrix_pahse, axis=1)
    dft = dft[S1.order('@'+matrix_pahse,0)]
    df1 = pd.concat([df1, dft], axis=1)


    df1.to_csv('Composition_'+material+'.csv', index=False)
    logging.info('The composition data is saved as Composition_'+material+'.csv')

    #Calculating the Hessian values from the Hessian_setup
    logging.info('Instance for hessian calculation is created with_single_equilibrium_calculation')
    phase_list=[matrix_pahse]+Active_phases
    for phase in phase_list:
        dfphase=pd.DataFrame()
        #get system for Hessian calculation with phase
        Hessian_setup= S1.Hessian_setup(calculation_setup, phase)
        logging.info('The instance for hessian calculation for the '+phase+' is created')


        for i in range(len(df1.index)):
            calc_temp = df1.at[i,"Temperature"]
            dfphase.at[i,"Temperature"]= calc_temp
            
            #set condition for Hessian calculation
            l = S1.Hessian_calculation(Hessian_setup, calc_temp)

            for j1, s1 in zip(range(len(solutes)),solutes):
                for j2, s2 in zip(range(len(solutes)),solutes):
                    if j1==j2:
                        ddgdxdx=l.get_value_of("mu("+s1+").x("+s2+")")-l.get_value_of("mu("+solvent+").x("+s2+")")
                        dfphase.at[i,'HSN('+s1+','+s2+')@'+phase]=ddgdxdx  

            for j1, s1 in zip(range(len(solutes)),solutes):
                for j2, s2 in zip(range(len(solutes)),solutes):
                    if j1<j2:
                        ddgdxdx=l.get_value_of("mu("+s1+").x("+s2+")")-l.get_value_of("mu("+solvent+").x("+s2+")")
                        dfphase.at[i,'HSN('+s1+','+s2+')@'+phase]=ddgdxdx
                    
            dfphase.to_csv("HSN_"+phase+".csv",index=False)
            logging.info('Hessian data for '+phase+' is saved to .csv file')

#Delete the cache folders
logging.info('The calculation is completed!')
shutil.rmtree(os.path.basename(material[:-4]) + "_cache")
os.chdir(current_dir)
shutil.rmtree("__pycache__")
logging.info('The cache files are delelted')
