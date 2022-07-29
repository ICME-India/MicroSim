import numpy as np
from skimage import measure
from vtk import vtkDataSetReader
import vtk

def front_Velocity(timeItretion, vtkData ,scalerValue,dt,saveT,dx, Is3d,depth_plot):
    ppt_major_axis =  [ [] for _ in range(timeItretion) ]
    ppt_count = np.empty(timeItretion)
    fvelocity = [ [] for _ in range(timeItretion) ]


    for t in range(timeItretion):
        grid_shape = vtkData[t].GetDimensions()
        
        if grid_shape[0] == 1:
            grid_reshape = ( grid_shape[2], grid_shape[1]  )
            
        elif grid_shape[1] == 1:
            grid_reshape = ( grid_shape[0], grid_shape[2]  )
            
        elif grid_shape[2] == 1:
            grid_reshape = ( grid_shape[1], grid_shape[0]  )
        
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
        
        
        thresh_bin = 0.5
        pf[pf < thresh_bin] = 0.0
        pf[pf >= thresh_bin] = 1.0

        ppt = pf
        ppt_labels, count_labels = measure.label(ppt, background=0, return_num=True)
        properties = measure.regionprops(ppt_labels)

        major_axis_list = ppt_major_axis[t]
        
        for props in properties:
            major_axis_list.append(props.major_axis_length )
        
        ppt_count[t] = len(properties)
        
    for t in range(1,timeItretion):
        velocitylist = fvelocity[t]
        for v in range(int(ppt_count[t])):
            velocitylist.append( float(dx)*  ((ppt_major_axis[t][v] - ppt_major_axis[t-1][v])/(2*float(dt)*float(saveT)))  )
    
    return fvelocity;
        
