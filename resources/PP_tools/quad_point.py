import numpy as np
from vtk import vtkDataSetReader


def quad_point(vtkData,dataset,infileDimension,timeItretion,Scalar_name):
    overall_tp = np.empty(len(timeItretion) ,  dtype = object)
    
    for t in range(len(timeItretion)):

        if(dataset == "UNSTRUCTURED_GRID"):
            grid_shape = infileDimension
            pf2 = vtkData[t].GetCellData().GetArray(Scalar_name)
        else:
            grid_shape = vtkData[t].GetDimensions()
            pf2 = vtkData[t].GetPointData().GetArray(Scalar_name)

        
        if grid_shape[0] == 1:
            grid_reshape = ( grid_shape[1], grid_shape[2]  )

        elif grid_shape[1] == 1:
            grid_reshape = ( grid_shape[0], grid_shape[2]  )

        elif grid_shape[2] == 1:
            grid_reshape = ( grid_shape[1], grid_shape[0]  )
            
        pf2 = np.copy(np.reshape(pf2, grid_reshape))
            
        
        thisset = []
        triple_point = []
        quad_point = []
        x_list = []
        y_list = []
        z_list = []
        for i in range(1,grid_reshape[0]):
            for j in range(1,grid_reshape[1]):
                for k in range(1,grid_reshape[2]):
                    thisset = []
                    thisset.append(pf2[i-1][j-1][k-1])
                    thisset.append(pf2[i-1][j-1][k])
                    thisset.append(pf2[i-1][j][k-1])
                    thisset.append(pf2[i-1][j][k])
                    thisset.append(pf2[i][j-1][k-1])
                    thisset.append(pf2[i][j-1][k])
                    thisset.append(pf2[i][j][k-1])
                    thisset.append(pf2[i][j][k])
                    thisset = set(thisset)

                if len(thisset) >= 3:
                    triple_point.append((i,j))
                    x_list.append(i)
                    y_list.append(j)
                    z_list.append(j)
                    overall_tp[t] = triple_point
                if len(thisset) >= 4:
                    quad_point.append((i,j))
                    x_list.append(i)
                    y_list.append(j)
                    z_list.append(j)
                    overall_tp[t] = quad_point

                    # print("total triple points ", len(triple_point))
                    #print(triple_point)
        
                    
    return overall_tp
