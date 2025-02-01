# Thermodynamic data for Microsim

## Purpose of the script
The TC-Python script generate thermodynamic data for microsim simulation.


### Input
- Relevant Thermodynamic and kinetic database
- Phases
- Elements and composition
- Temperature range

### Ouput
- Equilibrium composition of alloying elements for the phases in the given temperature range
- Hessian files for each phases

## How to use

There are three files in the script 
- input_tdbs.json
- tdbs_microsim.py
- tcsetup.py

The files are splited into three based on the functions

- The system can be defined by modifing the _input_tdbs.json_ file with the _Inputs_ **(Modification is required for every different system)**
- _tbds_microsim.py_ is the main python script for the simulation **(Modify if you need additional calculations or plotting to done)**
- _tcsetup.py_ contain function calls to the Thermocalc interface **(Don't modify)**



```sh
cd tdbs_microsim
python3 tdbs_microsim.py input_tdbs.json
```

## Modules

Using anacoda is the prefered way to setup the TC-Python environment.
The detailed discription on installing Anaconda is available [here](https://docs.anaconda.com/anaconda/install/linux/).

| Module | Purpose |
| ------ | ------ |
| tc_python | interface to the thermocalc |
| pandas | Manege and rearrange the data in a specific order |
| json | Read the input .json file |
| csv | save data in as a .csv file |
| os | create and manege the directory |
| shutill | force remove the cache files |
| loggings | create a log file for debugging |

## Aknowledgement

Dr. Dasari Mohan and Dr. Hariharan for their contribution in the work 

If you find any bugs or need any help in running the script, please reach me throgh the email

**_Arun C_**
**_mm23d002@smail.iitm.ac.in_**
