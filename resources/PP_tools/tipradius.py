import numpy as np
from skimage import measure
from vtk import vtkDataSetReader
import math

def tip_radius_calculate(vtkData, timeItretion ,Scalar_name,dx, Is3d, depth_plot):
    tip_radius = np.empty(len(timeItretion) ,  dtype = float)
    
    
    for t in range(len(timeItretion)):
        grid_shape = vtkData[t].GetDimensions()
        
        if grid_shape[0] == 1:
            grid_reshape = ( grid_shape[2], grid_shape[1]  )

        elif grid_shape[1] == 1:
            grid_reshape = ( grid_shape[0], grid_shape[2]  )

        elif grid_shape[2] == 1:
            grid_reshape = ( grid_shape[1], grid_shape[0]  )

        vtkPointData = vtkData[t].GetPointData().GetArray(Scalar_name)
        
        
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
        
        
        contours = measure.find_contours(pf, 0.5)
        
        data = contours[0]
        distance = np.add(np.multiply(data[:,0],data[:,0]),np.multiply(data[:,1],data[:,1]))
        d = math.sqrt(np.amax(distance))
        index = np.argmax(distance)
        
        
        #check for dendrite location
        
        if index == 0:
            tip_radius[t] = d
            continue
        
        
        elif data[index][0] < 20 or data[index][1] < 20:
            data = data[index:]
            data = data[:-index]
            data = [x - data[0][0] for x in data]
            
            data2 = np.copy(data)
            data2[:,1] = data2[:,1]*(-1)
            data2 = data2[::-1]
            half_dendrite = np.vstack((data,data2))
            half_X=half_dendrite[:,0]
            half_Y=half_dendrite[:,1]
            
            half2_dendrite = np.copy(half_dendrite)
            half2_dendrite[:,0] = half2_dendrite[:,0]*(-1)
            half2_dendrite = half2_dendrite[::-1]
            full_dendrite = np.vstack((half_dendrite,half2_dendrite))
            full_X=full_dendrite[:,0]
            full_Y=full_dendrite[:,1]
            
            X = np.array(full_dendrite[:,0])
            Y = np.array(full_dendrite[:,1])
            
            distance = np.add(np.multiply(X,X),np.multiply(Y,Y))
            d = math.sqrt(np.amax(distance))
            index = np.argmax(distance)
            
            if index <= 20:
                X_ = np.hstack((X[len(X) - 20 + index :] , X[ :index + 20] ))
                Y_ = np.hstack((Y[len(X) - 20 + index :] , Y[ :index + 20] ))

            elif index >= lex(X)-20:
                X_ = np.hstack((X[len(X) - 20:] , X[ : 20 - index] ))
                Y_ = np.hstack((Y[len(X) - 20:] , Y[ : 20 - index] ))
        else:
            X=np.array(data[:,0])
            Y=np.array(data[:,1])
            X_ = X[index-20:index+20]
            Y_ = Y[index-20:index+20]
        
        
        p_tip = np.polyfit(X_,Y_,2)
        d2ydx2 = np.polyder(p_tip,2)
        dydx   = np.polyder(p_tip,1)

        p      =  np.poly1d(d2ydx2)

        r_tip = p(X[index])

        r_tip = -1.0/r_tip

        tip_radius[t] = float(dx)*round(r_tip,1)
    
    return tip_radius
