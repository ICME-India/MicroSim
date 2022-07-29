import numpy as np
from skimage import measure
from vtk import vtkDataSetReader
import math

def front_undercooling_cal(vtkData, timeItretion ,Scalar_name,Scalar_name_all, Is3d,depth_plot):
    front_undercooling = np.empty(len(timeItretion) ,  dtype = float)
    
    grid_shape = vtkData[0].GetDimensions()

    if grid_shape[0] == 1:
        grid_reshape = ( grid_shape[2], grid_shape[1]  )

    elif grid_shape[1] == 1:
        grid_reshape = ( grid_shape[0], grid_shape[2]  )

    elif grid_shape[2] == 1:
        grid_reshape = ( grid_shape[1], grid_shape[0]  )
    
        
    for t in timeItretion:
 
        pf1 = vtkData[t].GetPointData().GetArray(Scalar_name)
        
        if('liquid' in Scalar_name_all):
            pf3 = vtkData[t].GetPointData().GetArray('liquid')
        elif('LIQUID' in Scalar_name_all):
            pf3 = vtkData[t].GetPointData().GetArray('LIQUID')
        elif('liq' in Scalar_name_all):
            pf3 = vtkData[t].GetPointData().GetArray('liq')
        elif('LIQ' in Scalar_name_all):
            pf3 = vtkData[t].GetPointData().GetArray('LIQ')
        else:
            print("liquid , liq, LIQ or LIQUID phase not present")
            return False
        
        
        if Is3d == 0:
            pf1 = np.copy(np.reshape(pf1, grid_reshape))
            pf3 = np.copy(np.reshape(pf3, grid_reshape))
            
      
        elif Is3d == 1 :
            pf1 = np.copy(np.reshape(pf1,  (grid_shape[2], grid_shape[1], grid_shape[0])))
            pf1 = pf1[depth_plot,:,:]
            
            pf3 = np.copy(np.reshape(pf3,  (grid_shape[2], grid_shape[1], grid_shape[0])))
            pf3 = pf3[depth_plot,:,:]
            
        elif Is3d == 2 :
            pf1 = np.copy(np.reshape(pf1,  (grid_shape[2], grid_shape[1], grid_shape[0])))
            pf1 = pf1[:,depth_plot,:]
            
            pf3 = np.copy(np.reshape(pf3,  (grid_shape[2], grid_shape[1], grid_shape[0])))
            pf3 = pf3[:,depth_plot,:]
            
        
        elif Is3d == 3 :
            pf1 = np.copy(np.reshape(pf1,  (grid_shape[2], grid_shape[1], grid_shape[0])))
            pf1 = pf1[:,:,depth_plot]
            
            pf3 = np.copy(np.reshape(pf3,  (grid_shape[2], grid_shape[1], grid_shape[0])))
            pf3 = pf3[:,:,depth_plot]
        

        Alpha_liq = pf1 - pf3

        contours = measure.find_contours(Alpha_liq,0 )

        flag = 0
        for contour in contours:
            if flag == 0:
                contours_data =  np.array(contour)
                flag = 1

            if flag == 1:
                contours_data = np.vstack(( contours_data,np.array(contour) ))


        pf5 = vtkData[t].GetPointData().GetArray('T')
        
        
        if Is3d == 0:
            pf5 = np.copy(np.reshape(pf5, grid_reshape))
            
      
        elif Is3d == 1 :
            pf5 = np.copy(np.reshape(pf5,  (grid_shape[2], grid_shape[1], grid_shape[0])))
            pf5 = pf5[depth_plot,:,:]
           
            
        elif Is3d == 2 :
            pf5 = np.copy(np.reshape(pf5,  (grid_shape[2], grid_shape[1], grid_shape[0])))
            pf5 = pf5[:,depth_plot,:]
           
        
        elif Is3d == 3 :
            pf5 = np.copy(np.reshape(pf5,  (grid_shape[2], grid_shape[1], grid_shape[0])))
            pf5 = pf5[:,:,depth_plot]
          

        sumT2 = 0
        for i in contours_data:
            x = int(i[0])
            y = int(i[1])

            distance = np.sqrt(  (i[0] - x )**2 + ( i[1] - y)**2   )
           
            #sumT2 = sumT2 + pf5[x][y]*(distance/np.sqrt(2)) + pf5[x+1][y+1]*(  1- (distance/np.sqrt(2))   ) 
            sumT2 = sumT2 + pf5[x][y] + distance*(pf5[x][y+1] - pf5[x][y] )

        approxT2  = sumT2/len(contours_data)

        front_undercooling[t] = approxT2
    
    return front_undercooling
