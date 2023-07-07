# Importing required libraries

import sys
import vtk
import pandas as pd
import numpy as np
import math
from vtkmodules.util.numpy_support import vtk_to_numpy
from vtkmodules.numpy_interface import dataset_adapter as dsa
from scipy.interpolate import interp1d

# Input VTK 
input_vtk = sys.argv[1]

# Output VTK
output_vtk = sys.argv[2]

# Number of cells along each direction
numcell = int(sys.argv[3])
numpt = (numcell+1)**2

# Temperature for evaluating composition during solidification
temp_n1 = float(sys.argv[4])

# Temperature for evaluating chemical potential during coarsening
temp_n = float(sys.argv[5])

# Read the binary VTK file using vtk.vtkUnstructuredGridReader
reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName(input_vtk)
reader.Update()

# Get the data from the reader using GetOutput()
output = reader.GetOutput()

# Convert the VTK data to a numpy array using vtk_to_numpy()
data = vtk_to_numpy(output.GetPoints().GetData())


reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName(input_vtk)
reader.ReadAllScalarsOn()
reader.ReadAllVectorsOn()
reader.Update()
usg = dsa.WrapDataObject(reader.GetOutput() )
# Getting the data arrays
array_phi = usg.PointData['phi'] # Assuming you know the name of the array
# array_phi is a child class type of numpy.ndarray type
array_mu1 = usg.PointData['mu_1']
array_mu2 = usg.PointData['mu_2']
#array_theta = usg.PointData['theta']
#array_psi = usg.PointData['psi']

# Dividing grid into quarters. 1-> upper left 2-> upper right 3-> bottom left 4-> bottom right
res = int(input("Please enter required resolution (integer power of 2 for bisecting multiple times): "))
location = []

print(" 1-> upper left 2-> upper right 3-> bottom left 4-> bottom right ")

for i in range(int(math.log(res, 2))):
    choice = int(input("Please choose the required quadrant within previous chosen quadrant: "))
    location.append(choice)

# initialising start and end point

start_point = 0
end_point = 0

# Defining function to return data having increased resolution by a magnitude of 2

#def increase_resolution(points_array, phi_array, mu1_array, mu2_array, psi_array, theta_array,quad):
def increase_resolution(points_array, phi_array, mu1_array, mu2_array, quad):
    # takes initial data and required quadrant as input and returns final data with increased resolution.
    points_total = len(points_array)
    total_x = len(set(points_array[:, 0]))
    total_y = len(set(points_array[:, 1]))
    side_x = int(math.ceil(total_x/2))
    side_y = int(math.ceil(total_y/2))

    if quad == 1:
        # upper left quadrant
        start_point = int(0.5*total_x*(total_x-1))
        end_point = int((total_x+0.5)*(total_x-1))

    elif quad == 2:
        # upper right quadrant
        start_point = int(0.5*(total_x*total_x - 1))
        end_point = total_x*total_x - 1

    elif quad == 3:
        # lower left quadrant
        start_point = 0
        end_point = int(0.5*(total_x*total_x - 1))

    elif quad == 4:
        # lower right quadrant
        start_point = int(0.5*(total_x - 1))
        end_point = int(0.5*total_x*(total_x-1)) - 1

    final_array = np.zeros((int(math.pow((side_x * 2 - 1), 2)), 3))

    tempx = np.linspace(points_array[start_point][0], points_array[start_point + side_x - 1][0], (2*side_x - 1))
    tempy = np.linspace(points_array[start_point][1], points_array[end_point][1], (2*side_x - 1))
    flag = 0
    for i in range(2*side_x - 1):
        for j in range(2*side_x - 1):
            final_array[flag][0] = tempx[j]
            final_array[flag][1] = tempy[i]
            flag += 1

    # Creating the final data array files
    phi_original = np.zeros(side_x * side_y)
    mu1_original = np.zeros(side_x * side_y)
    mu2_original = np.zeros(side_x * side_y)
    #psi_original = np.zeros(side_x * side_y)
    #theta_original = np.zeros(side_x * side_y)

    flag1 = 0
    for i in range(side_x):
        for j in range(side_y):
            phi_original[flag1] = phi_array[flag1 + start_point + (i * (side_x - 1))]
            mu1_original[flag1] = mu1_array[flag1 + start_point + (i * (side_x - 1))]
            mu2_original[flag1] = mu2_array[flag1 + start_point + (i * (side_x - 1))]
            #psi_original[flag1] = psi_array[flag1 + start_point + (i * (side_x - 1))]
            #theta_original[flag1] = theta_array[flag1 + start_point + (i * (side_x - 1))]
            flag1 += 1

    phix_array = np.zeros((2 * side_x - 1) * side_y)
    mu1x_array = np.zeros((2 * side_x - 1) * side_y)
    mu2x_array = np.zeros((2 * side_x - 1) * side_y)
    #psix_array = np.zeros((2 * side_x - 1) * side_y)
    #thetax_array = np.zeros((2 * side_x - 1) * side_y)

    flag_or = 0
    flag2 = 0
    for i in range(side_y):
        for j in range(2 * side_x - 1):
            if j % 2 == 0:
                phix_array[flag2] = phi_original[flag_or]
                mu1x_array[flag2] = mu1_original[flag_or]
                mu2x_array[flag2] = mu2_original[flag_or]
                #psix_array[flag2] = psi_original[flag_or]
                #thetax_array[flag2] = theta_original[flag_or]
                flag_or += 1
            else:
                phix_array[flag2] = 0.5 * (phi_original[flag_or - 1] + phi_original[flag_or])
                mu1x_array[flag2] = 0.5 * (mu1_original[flag_or - 1] + mu1_original[flag_or])
                mu2x_array[flag2] = 0.5 * (mu2_original[flag_or - 1] + mu2_original[flag_or])
                #psix_array[flag2] = 0.5 * (psi_original[flag_or - 1] + psi_original[flag_or])
                #thetax_array[flag2] = 0.5 * (theta_original[flag_or - 1] + theta_original[flag_or])
            flag2 += 1


    final_phi_array = np.zeros((2 * side_x - 1) * (2 * side_y - 1))
    final_mu1_array = np.zeros((2 * side_x - 1) * (2 * side_y - 1))
    final_mu2_array = np.zeros((2 * side_x - 1) * (2 * side_y - 1))
    #final_psi_array = np.zeros((2 * side_x - 1) * (2 * side_y - 1))
    #final_theta_array = np.zeros((2 * side_x - 1) * (2 * side_y - 1))

    flag_new = 0
    flag_old = 0
    for i in range(2 * side_x - 1):
        for j in range(2 * side_y - 1):
            if i % 2 == 0:
                final_phi_array[flag_new] = phix_array[flag_old]
                final_mu1_array[flag_new] = mu1x_array[flag_old]
                final_mu2_array[flag_new] = mu2x_array[flag_old]
                #final_psi_array[flag_new] = psix_array[flag_old]
                #final_theta_array[flag_new] = thetax_array[flag_old]
                flag_old += 1
            else:
                final_phi_array[flag_new] = 0.5 * (
                            phix_array[flag_old + j] + phix_array[flag_old - (2 * side_x - 1) + j])
                final_mu1_array[flag_new] = 0.5 * (
                        mu1x_array[flag_old + j] + mu1x_array[flag_old - (2 * side_x - 1) + j])
                final_mu2_array[flag_new] = 0.5 * (
                        mu2x_array[flag_old + j] + mu2x_array[flag_old - (2 * side_x - 1) + j])
                #final_psi_array[flag_new] = 0.5 * (
                        #psix_array[flag_old + j] + psix_array[flag_old - (2 * side_x - 1) + j])
                #final_theta_array[flag_new] = 0.5 * (
                        #thetax_array[flag_old + j] + thetax_array[flag_old - (2 * side_x - 1) + j])
            flag_new += 1

    #return final_array, final_phi_array, final_mu1_array, final_mu2_array, final_psi_array, final_theta_array
    return final_array, final_phi_array, final_mu1_array, final_mu2_array


# Calling function as many times as required based on user input resolution with input data being output data of previous iteration


for i in range(int(math.log(res, 2))):
    if i == 0:
        final_array, final_phi_array, final_mu1_array, final_mu2_array = increase_resolution(data[0:numpt], array_phi, array_mu1, array_mu2, location[i])
    else:
        final_array, final_phi_array, final_mu1_array, final_mu2_array = increase_resolution(final_array, final_phi_array, final_mu1_array, final_mu2_array, location[i])


# Creating new vtk file with required data
shiftx = np.min(final_array[:,0])
shifty = np.min(final_array[:,1])
#print(shiftx)
#print(shifty)

# Create the points
points = vtk.vtkPoints()
n_points = len(final_array)
n_points_per_side = int(np.sqrt(n_points))
spacing = 1.0 / (n_points_per_side - 1)
for i in range(n_points):
    final_array[i][0] = final_array[i][0] - shiftx
    final_array[i][1] = final_array[i][1] - shifty
    x = final_array[i][0]
    y = final_array[i][1]
    z = final_array[i][2]
    points.InsertNextPoint(x, y, z)

# Create the unstructured grid
grid = vtk.vtkUnstructuredGrid()
grid.SetPoints(points)

# Create the cells
n_cells = (n_points_per_side - 1) * (n_points_per_side - 1)
cell_type = vtk.VTK_QUAD
grid.Allocate(n_cells, n_cells)
for i in range(n_points_per_side - 1):
    for j in range(n_points_per_side - 1):
        cell = vtk.vtkQuad()
        p1 = j + i * n_points_per_side
        p2 = p1 + 1
        p3 = p2 + n_points_per_side
        p4 = p1 + n_points_per_side
        cell.GetPointIds().SetId(0, p1)
        cell.GetPointIds().SetId(1, p2)
        cell.GetPointIds().SetId(2, p3)
        cell.GetPointIds().SetId(3, p4)
        grid.InsertNextCell(cell_type, cell.GetPointIds())

# Create the data arrays

data_array = vtk.vtkFloatArray()
data_array.SetName("phi")
data_array.SetNumberOfComponents(1)
data_array.SetNumberOfTuples(n_points)
for i in range(n_points):
    phi = final_phi_array[i]
    data_array.SetValue(i, phi)

# Add the data array to the grid
grid.GetPointData().AddArray(data_array)

data_array = vtk.vtkFloatArray()
data_array.SetName("mu_1")
data_array.SetNumberOfComponents(1)
data_array.SetNumberOfTuples(n_points)
for i in range(n_points):
    mu1 = final_mu1_array[i]
    data_array.SetValue(i, mu1)

# Add the data array to the grid
grid.GetPointData().AddArray(data_array)

data_array = vtk.vtkFloatArray()
data_array.SetName("mu_2")
data_array.SetNumberOfComponents(1)
data_array.SetNumberOfTuples(n_points)
for i in range(n_points):
    mu2 = final_mu2_array[i]
    data_array.SetValue(i, mu2)

# Add the data array to the grid
grid.GetPointData().AddArray(data_array)

#data_array = vtk.vtkFloatArray()
#data_array.SetName("psi")
#data_array.SetNumberOfComponents(1)
#data_array.SetNumberOfTuples(n_points)
#for i in range(n_points):
#    psi = final_psi_array[i]
#    data_array.SetValue(i, psi)

# Add the data array to the grid
#grid.GetPointData().AddArray(data_array)

#data_array = vtk.vtkFloatArray()
#data_array.SetName("theta")
#data_array.SetNumberOfComponents(1)
#data_array.SetNumberOfTuples(n_points)
#for i in range(n_points):
#    theta = final_theta_array[i]
#    data_array.SetValue(i, theta)

# Add the data array to the grid
#grid.GetPointData().AddArray(data_array)

# Write the grid to a binary VTK file
#writer = vtk.vtkUnstructuredGridWriter()
#writer.SetFileName(output_vtk)
#writer.SetInputData(grid)
#writer.SetFileTypeToBinary()
#writer.Write()
    

tdb_c = pd.read_csv('Composition_FCC_A1.csv');
T_tdb_c = np.array(tdb_c.loc[:,'Temp']).squeeze()
cs_1 = np.array(tdb_c.loc[:,'Al@fcc']).squeeze()
cs_2 = np.array(tdb_c.loc[:,'Mo@fcc']).squeeze()
cl_1 = np.array(tdb_c.loc[:,'Al@liq']).squeeze()
cl_2 = np.array(tdb_c.loc[:,'Mo@liq']).squeeze()

cs_1_func = interp1d(T_tdb_c, cs_1, kind='cubic')
cs_2_func = interp1d(T_tdb_c, cs_2, kind='cubic')
cl_1_func = interp1d(T_tdb_c, cl_1, kind='cubic')
cl_2_func = interp1d(T_tdb_c, cl_2, kind='cubic')

cs_1_final = cs_1_func(temp_n1).item();
cs_2_final = cs_2_func(temp_n1).item();
cl_1_final = cl_1_func(temp_n1).item();
cl_2_final = cl_2_func(temp_n1).item();
# use temp to find compo

#compo_tdb_func = interp1d(c_L_tdb, c_S_tdb, kind='cubic')

tdb_hliq = pd.read_csv('HSN_LIQUID.csv');
T_tdb_hl = np.array(tdb_hliq.loc[:,'Temp']).squeeze()
hl_11 = np.array(tdb_hliq.loc[:,'HSN(Al,Al)@liq']).squeeze()
hl_22 = np.array(tdb_hliq.loc[:,'HSN(Mo,Mo)@liq']).squeeze()
hl_12 = np.array(tdb_hliq.loc[:,'HSN(Al,Mo)@liq']).squeeze()

hl_11_func = interp1d(T_tdb_hl, hl_11, kind='cubic')
hl_22_func = interp1d(T_tdb_hl, hl_22, kind='cubic')
hl_12_func = interp1d(T_tdb_hl, hl_12, kind='cubic')

hl_11_final = hl_11_func(temp_n1).item();
hl_22_final = hl_22_func(temp_n1).item();
hl_12_final = hl_12_func(temp_n1).item();
# use temp to find hessian

tdb_hsol = pd.read_csv('HSN_FCC_A1.csv');
T_tdb_hs = np.array(tdb_hsol.loc[:,'Temp']).squeeze()
hs_11 = np.array(tdb_hsol.loc[:,'HSN(Al,Al)@fcc']).squeeze()
hs_22 = np.array(tdb_hsol.loc[:,'HSN(Mo,Mo)@fcc']).squeeze()
hs_12 = np.array(tdb_hsol.loc[:,'HSN(Al,Mo)@fcc']).squeeze()

hs_11_func = interp1d(T_tdb_hs, hs_11, kind='cubic')
hs_22_func = interp1d(T_tdb_hs, hs_22, kind='cubic')
hs_12_func = interp1d(T_tdb_hs, hs_12, kind='cubic')

hs_11_final = hs_11_func(temp_n1).item();
hs_22_final = hs_22_func(temp_n1).item();
hs_12_final = hs_12_func(temp_n1).item();
# use temp to find hessian

B_s1 = (hl_11_final*cl_1_final - hs_11_final*cs_1_final) + (hl_12_final*cl_2_final - hs_12_final*cs_2_final)
B_s2 = (hl_22_final*cl_2_final - hs_22_final*cs_2_final) + (hl_12_final*cl_1_final - hs_12_final*cs_1_final)

detdmudc_s = hs_11_final*hs_22_final - hs_12_final*hs_12_final

dcdmu_s11 = hs_22_final/detdmudc_s
dcdmu_s22 = hs_11_final/detdmudc_s
dcdmu_s12 = -hs_12_final/detdmudc_s

cs1_array = np.zeros(n_points)
cs2_array = np.zeros(n_points)

ppt_mu1_array = np.zeros(n_points)
ppt_mu2_array = np.zeros(n_points)

tdb_hs1 = pd.read_csv('HSN_L12_FCC_A1_1.csv');
T_tdb = np.array(tdb_hs1.loc[:,'Temp']).squeeze()
hs1_11 = np.array(tdb_hs1.loc[:,'HSN(Al,Al)@fcc1']).squeeze()
hs1_22 = np.array(tdb_hs1.loc[:,'HSN(Mo,Mo)fcc1']).squeeze()
hs1_12 = np.array(tdb_hs1.loc[:,'HSN(Al,Mo)fcc1']).squeeze()

hs1_11_func = interp1d(T_tdb, hs1_11, kind='cubic')
hs1_22_func = interp1d(T_tdb, hs1_22, kind='cubic')
hs1_12_func = interp1d(T_tdb, hs1_12, kind='cubic')

hs1_11_final = hs1_11_func(temp_n).item();
hs1_22_final = hs1_22_func(temp_n).item();
hs1_12_final = hs1_12_func(temp_n).item();

#hs1_11_final = hs1_11[-1]
#hs1_22_final = hs1_22[-1]
#hs1_12_final = hs1_12[-1]
# use temp to find hessian

cs1_sum = 0.
cs2_sum = 0.

for i in range(n_points):
    cs1_array[i] = dcdmu_s11*(final_mu1_array[i] - B_s1) + dcdmu_s12*(final_mu2_array[i] - B_s2)
    cs2_array[i] = dcdmu_s12*(final_mu1_array[i] - B_s1) + dcdmu_s22*(final_mu2_array[i] - B_s2)
    
    cs1_sum = cs1_sum + cs1_array[i]
    cs2_sum = cs2_sum + cs2_array[i]
    
    #cs1_array[i] = 0.0684597
    #cs2_array[i] = 0.00114417
    
    #cs1_array[i] = 0.0987016
    #cs2_array[i] = 0.0741531

    #cs1_array[i] = 0.12326377737601378
    #cs2_array[i] = 0.005887372806647023
    
    ppt_mu1_array[i] = hs1_11_final*cs1_array[i] + hs1_12_final*cs2_array[i]
    ppt_mu2_array[i] = hs1_12_final*cs1_array[i] + hs1_22_final*cs2_array[i]


cs1_avg = cs1_sum/n_points
cs2_avg = cs2_sum/n_points

print(" average c_1 "+str(cs1_avg))
print(" average c_2 "+str(cs2_avg))


with open('mu_1initial', 'w') as f0:
    f0.write("internalField   nonuniform List<scalar>\n"+str(numcell**2)+"\n(")
    for i in range(n_points_per_side - 1):
        for j in range(n_points_per_side - 1):
            p1 = j + i * n_points_per_side
            p2 = p1 + 1
            p3 = p2 + n_points_per_side
            p4 = p1 + n_points_per_side
            mu1 = 0.25*(ppt_mu1_array[p1] + ppt_mu1_array[p2] + ppt_mu1_array[p3] + ppt_mu1_array[p4])
            #mu1 = 0.25*(cs1_array[p1] + cs1_array[p2] + cs1_array[p3] + cs1_array[p4])
            f0.write(str(mu1)+"\n")
    f0.write(")\n;")
            
with open('mu_2initial', 'w') as f1:
    f1.write("internalField   nonuniform List<scalar>\n"+str(numcell**2)+"\n(")
    for i in range(n_points_per_side - 1):
        for j in range(n_points_per_side - 1):
            p1 = j + i * n_points_per_side
            p2 = p1 + 1
            p3 = p2 + n_points_per_side
            p4 = p1 + n_points_per_side
            mu2 = 0.25*(ppt_mu2_array[p1] + ppt_mu2_array[p2] + ppt_mu2_array[p3] + ppt_mu2_array[p4])
            #mu2 = 0.25*(cs2_array[p1] + cs2_array[p2] + cs2_array[p3] + cs2_array[p4])
            f1.write(str(mu2)+"\n")
    f1.write(")\n;")


data_array = vtk.vtkFloatArray()
data_array.SetName("c1")
data_array.SetNumberOfComponents(1)
data_array.SetNumberOfTuples(n_points)
for i in range(n_points):
    c1 = cs1_array[i]
    data_array.SetValue(i, c1)

# Add the data array to the grid
grid.GetPointData().AddArray(data_array)

data_array = vtk.vtkFloatArray()
data_array.SetName("c2")
data_array.SetNumberOfComponents(1)
data_array.SetNumberOfTuples(n_points)
for i in range(n_points):
    c2 = cs2_array[i]
    data_array.SetValue(i, c2)

# Add the data array to the grid
grid.GetPointData().AddArray(data_array)

# Write the grid to a binary VTK file
writer = vtk.vtkUnstructuredGridWriter()
writer.SetFileName(output_vtk)
writer.SetInputData(grid)
writer.SetFileTypeToBinary()
writer.Write()
