
from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()
output = XDMFReader(FileNames=source_files)
output.PointArrayStatus = [' liquid', 'Composition_Al', 'Mu_Al', 'alpha']
animationScene1 = GetAnimationScene()
animationScene1.UpdateAnimationUsingDataTimeSteps()
SaveData('/home/sagar/Desktop/testing_h5tovtk/DATA/bcbcbcbcbbbcvb.vtk', proxy=output, Writetimestepsasfileseries=1,Firsttimestep=0, Lasttimestep=-1,Timestepstride=1)
Delete(output)
del output
