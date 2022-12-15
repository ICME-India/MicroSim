import numpy as np
import matplotlib.pyplot as plt
from pymks.stats import autocorrelate
from pymks.stats import correlate
from pymks import PrimitiveBasis
from sklearn.decomposition import PCA
import sys
import os
import glob
from vtk import vtkDataSetReader


def two_point_correlation( vtkData,dataset,infileDimension,timeItretion ,scalerValue,Is3d,depth_plot):
      
    pointData = [ [] for _ in range(timeItretion) ]
    domain = [0, 1]
    n_states = 2
    p_basis = PrimitiveBasis(n_states=n_states, domain=domain)
    i, j = 0, 0
    all_correlations = [(i, j)]
    n_components = 3
    
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

        thresh_bin = 0.5
        pf[pf < thresh_bin] = 0.0
        pf[pf >= thresh_bin] = 1.0

        pointData[t] = np.copy(np.flipud(np.array(pf)))
    
    X_stats = correlate(np.array(pointData), p_basis, correlations=all_correlations, periodic_axes=(0,1))
    X_reshaped = X_stats.reshape((X_stats.shape[0], -1))
    pca = PCA(n_components)
    n_samples = timeItretion
    pca_out = pca.fit_transform(X_reshaped)

    return X_stats, pca_out;
