import numpy as np
from vtk import vtkDataSetReader



def getShiftData( vtkData,dataset,infileDimension,timeItretion,scalerValue,savet,arr ,Is3d,depth_plot):
      

    if(dataset == "UNSTRUCTURED_GRID"):
        grid_shape = infileDimension
        vtkPointData = vtkData[0].GetCellData().GetArray(scalerValue)
    else:
        grid_shape = vtkData[0].GetDimensions()
        vtkPointData = vtkData[0].GetPointData().GetArray(scalerValue)
        
    if grid_shape[0] == 1:
        grid_reshape = ( grid_shape[1], grid_shape[2]  )
        
    elif grid_shape[1] == 1:
        grid_reshape = ( grid_shape[0], grid_shape[2]  )
        
    elif grid_shape[2] == 1:
        grid_reshape = ( grid_shape[0], grid_shape[1]  )


    if Is3d == 0:
        pf = np.copy(np.reshape(vtkPointData, grid_reshape))
        
    
    elif Is3d == 1 :
        pf = np.copy(np.reshape(vtkPointData,  (grid_shape[2], grid_shape[1], grid_shape[0])))
        pf = pf[depth_plot,:,:]

    
    elif Is3d == 2 :
        pf = np.copy(np.reshape(vtkPointData,  (grid_shape[2], grid_shape[1], grid_shape[0])))
        pf = pf[:,depth_plot,:]
        
    
    elif Is3d == 3 :
        pf = np.copy(np.reshape(vtkPointData,  (grid_shape[2], grid_shape[1], grid_shape[0])))
        pf = pf[:,:,depth_plot]
    
    shiftData = np.transpose(pf).tolist()
    #print(shiftData[0])

    for t in np.arange(1,len(timeItretion)-1,1):
        

        if(dataset == "UNSTRUCTURED_GRID"):
            vtkPointData = vtkData[t].GetCellData().GetArray(scalerValue)
        else:
            vtkPointData = vtkData[t].GetPointData().GetArray(scalerValue)
    
        if Is3d == 0:
            pf = np.copy(np.reshape(vtkPointData, grid_reshape))
            
      
        elif Is3d == 1 :
            pf = np.copy(np.reshape(vtkPointData,  (grid_shape[2], grid_shape[1], grid_shape[0])))
            pf = pf[depth_plot,:,:]

        
        elif Is3d == 2 :
            pf = np.copy(np.reshape(vtkPointData,  (grid_shape[2], grid_shape[1], grid_shape[0])))
            pf = pf[:,depth_plot,:]
            
        
        elif Is3d == 3 :
            pf = np.copy(np.reshape(vtkPointData,  (grid_shape[2], grid_shape[1], grid_shape[0])))
            pf = pf[:,:,depth_plot]
            
        newshiftData = np.transpose(pf).tolist()

        if(int(savet*(t+1)) == int(arr[t][0])):
            
            shiftData[int(arr[t][1]) : ] = newshiftData[:-int(arr[t][1])]
            shiftData.extend(newshiftData[-int(arr[t][1]) : ])
        
        
    return np.transpose(np.array(shiftData))
