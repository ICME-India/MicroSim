import numpy as np
from vtk import vtkDataSetReader
from skimage import measure


def getContourData( vtkData,dataset,infileDimension,timeItretion,scalerValue,contour_value ,Is3d,depth_plot):
    

    pointData = [ [] for _ in range(timeItretion)]
  
    for t in range(timeItretion):

        if(dataset == "UNSTRUCTURED_GRID"):
            grid_shape = infileDimension
            vtkPointData = vtkData[t].GetCellData().GetArray(scalerValue)
        else:
            grid_shape = vtkData[t].GetDimensions()
            vtkPointData = vtkData[t].GetPointData().GetArray(scalerValue)
        
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

        
        contours = measure.find_contours(pf, contour_value)
        
        pointData[t] = np.copy(np.array(contours, dtype=object)) 
        
        
    return pointData
