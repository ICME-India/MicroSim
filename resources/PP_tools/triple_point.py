import numpy as np
from vtk import vtkDataSetReader
import vtk


def triple_point(vtkData,timeItretion,Scalar_name):
    overall_tp = np.empty(len(timeItretion) ,  dtype = object)
    
    for t in range(len(timeItretion)):
        
        pf2 = vtkData[t].GetPointData().GetArray(Scalar_name)
        grid_shape = vtkData[t].GetDimensions()
        if grid_shape[0] == 1:
            grid_reshape = ( grid_shape[1], grid_shape[2]  )

        elif grid_shape[1] == 1:
            grid_reshape = ( grid_shape[0], grid_shape[2]  )

        elif grid_shape[2] == 1:
            grid_reshape = ( grid_shape[1], grid_shape[0]  )
            
        pf2 = np.copy(np.reshape(pf2, grid_reshape))
            
        
        thisset = []
        triple_point = []
        x_list = []
        y_list = []
        for i in range(1,grid_reshape[0]):
            for j in range(1,grid_reshape[1]):
                thisset = []
                thisset.append(pf2[i-1][j-1])
                thisset.append(pf2[i-1][j])
                thisset.append(pf2[i][j-1])
                thisset.append(pf2[i][j])
                thisset = set(thisset)

                if len(thisset) >= 3:
                    triple_point.append((i,j))
                    x_list.append(i)
                    y_list.append(j)

                    # print("total triple points ", len(triple_point))
                    #print(triple_point)
                    overall_tp[t] = triple_point
        
                    
    return overall_tp
