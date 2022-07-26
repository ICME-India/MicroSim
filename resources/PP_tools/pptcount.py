import numpy as np
from skimage import io
from skimage import measure
import matplotlib.pyplot as plt
from vtk import vtkDataSetReader
from vtkmodules.vtkIOLegacy import vtkUnstructuredGridReader

def load_ppt_property(timeItretion, vtkData ,scalerValue ,Is3d,depth_plot):
    ppt_count = np.empty(len(timeItretion) ,  dtype = int)
    ppt_area = [ [] for _ in range(len(timeItretion)) ]
    ppt_radius = [ [] for _ in range(len(timeItretion)) ]
    ppt_coords = [ [] for _ in range(len(timeItretion)) ]
    ppt_major_axis =  [ [] for _ in range(len(timeItretion)) ]
    ppt_minor_axis = [ [] for _ in range(len(timeItretion)) ]
    ppt_perimeter = [ [] for _ in range(len(timeItretion)) ]

    for t in range(len(timeItretion)):
        
        grid_shape = vtkData[t].GetDimensions()
        
        if grid_shape[0] == 1:
            grid_reshape = ( grid_shape[2], grid_shape[1]  )

        elif grid_shape[1] == 1:
            grid_reshape = ( grid_shape[0], grid_shape[2]  )

        elif grid_shape[2] == 1:
            grid_reshape = ( grid_shape[1], grid_shape[0]  )
        

        #pf = np.copy(np.reshape(pf, grid_reshape))
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

        ppt_area_list = ppt_area[t]
        ppt_radius_list = ppt_radius[t]
        major_axis_list = ppt_major_axis[t]
        ppt_coords_list = ppt_coords[t]
        minor_axis_list = ppt_minor_axis[t]
        ppt_perimeter_list = ppt_perimeter[t]
        
        
        for props in properties:
            ppt_area_list.append(props.area)
            ppt_radius_list.append(np.sqrt(props.area/np.pi))
            ppt_coords_list.append(props.centroid)
            major_axis_list.append(round(props.major_axis_length , 3))
            minor_axis_list.append(round(props.minor_axis_length , 3))
            ppt_perimeter_list.append(round( props.perimeter ,3))
        
        ppt_count[t] = int(len(properties))
   
    
    return [ppt_area, ppt_radius, ppt_coords, ppt_major_axis , ppt_minor_axis , ppt_count, ppt_perimeter];



def volFrac_SA_Vol(timeItretion, vtkData ,scalerValue ,Is3d,depth_plot):
   
    ppt_area = [ [] for _ in range(len(timeItretion)) ]
    ppt_perimeter = [ [] for _ in range(len(timeItretion)) ]
    
    total_volume = np.empty(len(timeItretion) )
    total_SA = np.empty(len(timeItretion) )
    volume_fraction = np.empty(len(timeItretion) )

    for t in range(len(timeItretion)):
        
        grid_shape = vtkData[t].GetDimensions()
        
        if grid_shape[0] == 1:
            grid_reshape = ( grid_shape[2], grid_shape[1]  )

        elif grid_shape[1] == 1:
            grid_reshape = ( grid_shape[0], grid_shape[2]  )

        elif grid_shape[2] == 1:
            grid_reshape = ( grid_shape[1], grid_shape[0]  )
            
        
        vtkPointData = vtkData[t].GetPointData().GetArray(scalerValue)
            
        if Is3d == 0:
            pf_blob = np.copy(np.reshape(vtkPointData, grid_reshape))
            
        elif Is3d == 1 :
            pf_blob = np.copy(np.reshape(vtkPointData,  (grid_shape[2], grid_shape[1], grid_shape[0])))
            pf_blob = pf_blob[depth_plot,:,:]
            
        
        elif Is3d == 2 :
            pf_blob = np.copy(np.reshape(vtkPointData,  (grid_shape[2], grid_shape[1], grid_shape[0])))
            pf_blob = pf_blob[:,depth_plot,:]
            
        
        elif Is3d == 3 :
            pf_blob = np.copy(np.reshape(vtkPointData,  (grid_shape[2], grid_shape[1], grid_shape[0])))
            pf_blob = pf_blob[:,:,depth_plot]
        
        
        thresh_bin = 0.5
        pf_blob[pf_blob < thresh_bin] = 0.0
        pf_blob[pf_blob >= thresh_bin] = 1.0

        ppt = pf_blob
        ppt_labels, count_labels = measure.label(ppt, background=0, return_num=True)
       
        properties = measure.regionprops(ppt_labels)

        ppt_area_list = ppt_area[t]
        ppt_perimeter_list = ppt_perimeter[t]
        
        
        for props in properties:
            ppt_area_list.append(props.area)
            ppt_perimeter_list.append(round( props.perimeter ,3))
            
        total_volume[t] = sum(ppt_area_list)
        total_SA[t] = sum(ppt_perimeter_list)
        volume_fraction[t]  = total_volume[t] /  (grid_reshape[0]*grid_reshape[1])

    
    return [total_volume, volume_fraction , total_SA];
