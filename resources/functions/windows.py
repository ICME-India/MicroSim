import sys, os, glob
from xml.dom import minidom

def H5toVTK_Func(self):
    if self.h5tovtk_h5radio.isChecked():
            
        h5_outhead, h5_outtail = os.path.split(self.h5tovtk_outputloc.text())
        
        split_text = "_" + self.h5tovtk_sTime.text()+".h5"
        
        h5_outputfilename = h5_outtail.split(split_text)


        h5tovtk_infileLoc_head, h5tovtk_infileLoc_tail  = os.path.split(self.h5tovtk_infileLoc.text())
        
        #reading vtk output file dir
        
        if self.h5tovtk_vtkoutput.text() == "":
            vtk_output_fname = h5_outhead + "/" + h5_outputfilename[0]+ ".vtk"
        
        elif os.path.isdir(self.h5tovtk_vtkoutput.text()):
            vtk_output_fname = self.h5tovtk_vtkoutput.text() + "/" + h5_outputfilename[0]+ ".vtk"
            
        else:
            return False    
        

        h5_outheadLinux = "/mnt/" +  str. lower(h5_outhead[0])+h5_outhead[2:]
        h5toxmf_cmd = "cd " +h5_outheadLinux + "; cd ..;" + "cp /mnt/c/Users/%username%/Documents/MicroSim/Grand_potential_Finite_difference_2D_MPI/write_xdmf write_xdmf ; ./write_xdmf " + h5tovtk_infileLoc_tail + " "+ h5_outputfilename[0] + " " +  self.h5tovtk_sTime.text() + " "+ self.h5tovtk_eTime.text() 
        
        os.system(f'wsl ~ -e sh -c "{h5toxmf_cmd }" ')
        
        xmlfile = minidom.parse( h5_outhead +"/" + h5_outputfilename[0] +"_" +self.h5tovtk_sTime.text() + ".xmf")
        models = xmlfile.getElementsByTagName('Attribute')
        scalar_names_xml = []
        
        for elem in models:
            scalar_names_xml.append(elem.attributes['Name'].value)
            
            
        source_files = glob.glob(h5_outhead +"/"+ h5_outputfilename[0] +'*.xmf')
        xmf_names = []
        for files in source_files:
            xmf_names.append(str(files))
            
        xmf_names = sorted(xmf_names,key=lambda x: int(x.split('_')[-1].split('.')[0]))
            
            
        python_script_for_vtk = "from paraview.simple import *\nparaview.simple._DisableFirstRenderCameraReset()\noutput = XDMFReader(FileNames=" +str(xmf_names) +")\noutput.PointArrayStatus = " +str(scalar_names_xml) + "\nanimationScene1 = GetAnimationScene()\nanimationScene1.UpdateAnimationUsingDataTimeSteps()\nprint('Converting to vtk. Please wait...')\nSaveData('" + vtk_output_fname + "', proxy=output, Writetimestepsasfileseries=1,Firsttimestep=0, Lasttimestep=-1,Timestepstride="+  str(int(self.h5tovtk_eskip.text()) + 1) +")\nprint('Done')"
        
        fileName = h5_outhead + "/h5tovtk.py"
        f = open(fileName, "w")
        f.write( python_script_for_vtk )
        f.close()
        
        xmftovtk_cmd = "cd " +h5_outhead + "& pvpython h5tovtk.py & del h5tovtk.py & cd .. & del write_xdmf"
        print(xmftovtk_cmd)
        os.system("start \"\" cmd /k \""+ xmftovtk_cmd +"\"")
        
    elif self.h5tovtk_xmfradio.isChecked():
        
        h5_outhead, h5_outtail = os.path.split(self.h5tovtk_outputloc.text())
        
        split_text = "_"+self.h5tovtk_sTime.text()+".xmf"
        
        h5_outputfilename = h5_outtail.split(split_text)
        
        #reading vtk output file dir
        
        if self.h5tovtk_vtkoutput.text() == "":
            vtk_output_fname = h5_outhead + "/" + h5_outputfilename[0]+ "_.vtk"
        
        elif os.path.isdir(self.h5tovtk_vtkoutput.text()):
            vtk_output_fname = self.h5tovtk_vtkoutput.text() + "/" + h5_outputfilename[0]+ "_.vtk"
            
        else:
            return False
        
        xmlfile = minidom.parse( h5_outhead +"/" + h5_outputfilename[0] + "_" +self.h5tovtk_sTime.text() + ".xmf")
        models = xmlfile.getElementsByTagName('Attribute')
        scalar_names_xml = []
        
        for elem in models:
            scalar_names_xml.append(elem.attributes['Name'].value)
        
        source_files = glob.glob(h5_outhead +"/"+ h5_outputfilename[0] +'*.xmf')
        xmf_names = []
        for files in source_files:
            xmf_names.append(str(files))
            
        xmf_names = sorted(xmf_names,key=lambda x: int(x.split('_')[-1].split('.')[0]))
        
        python_script_for_vtk = "from paraview.simple import *\nparaview.simple._DisableFirstRenderCameraReset()\noutput = XDMFReader(FileNames=" +str(xmf_names) +")\noutput.PointArrayStatus = " +str(scalar_names_xml) + "\nanimationScene1 = GetAnimationScene()\nanimationScene1.UpdateAnimationUsingDataTimeSteps()\nprint('Converting to vtk. Please wait...')\nSaveData('" + vtk_output_fname + "', proxy=output, Writetimestepsasfileseries=1,Firsttimestep=0, Lasttimestep=-1,Timestepstride="+  str(int(self.h5tovtk_eskip.text()) + 1) +")\nprint('Done')"
        
        fileName = h5_outhead + "/h5tovtk.py"
        f = open(fileName, "w")
        f.write( python_script_for_vtk )
        f.close()
        
        xmftovtk_cmd = "cd " +h5_outhead + ";pvpython h5tovtk.py;rm h5tovtk.py;cd ..;rm write_xdmf"
        os.system("start \"\" cmd /k \"cd "+ xmftovtk_cmd +"\"")
        
def paraviewFunc(self):
    list_of_files = glob.glob(self.runDir + "/DATA/"+ self.output.text() +"*.*")
    if len(list_of_files) == 0:
        self.finish_error.setText("Sorry, output file not found.")
        return
    
    latest_file = max(list_of_files, key=os.path.getctime)

    paraviewcmd = "cd "+self.runDir  +"/DATA/ & paraview " +latest_file 
    os.system("start \"\" cmd /k \""+ paraviewcmd +"\"")
    
    
        
def SolverExecute(self):
    if self.radio_GP.isChecked():
            
        runDirGE = "/mnt/c" + self.runDir[2:]
        commandLine ="cd /mnt/c/Users/%username%/Documents/MicroSim/Grand_potential_Finite_difference_2D_MPI/; python3 GEdata_writer.py " +runDirGE +"/"+self.infile.text() + " ;make clean;make; cp microsim_gp " + runDirGE +  ";cd " + runDirGE + ";mpirun.mpich -np 4 ./microsim_gp "  +self.infile.text()+" "+self.filling.text()+" "+self.output.text() + " 2 2"
        
        os.system(f'cmd /c start cmd /c wsl ~ -e sh -c "{commandLine }" ')
    
    elif self.radio_KKR.isChecked():

        runDirGE = "/mnt/c" + self.runDir[2:]
        commandLine ="cd /mnt/c/Users/%username%/Documents/MicroSim/KKS_CuFFT/; python3 GEdata_writer.py " +runDirGE +"/"+self.infile.text() + " ;make clean;make; cp microsim_kks_cufft " + runDirGE +  ";cd " + runDirGE + ";mpirun -np 4 ./microsim_kks_cufft "  +self.infile.text()+" "+self.filling.text()+" "+self.output.text() 
        
        os.system(f'cmd /c start cmd /c wsl ~ -e sh -c "{commandLine }" ')


    elif self.radio_KKS2.isChecked():

        runDirGE = "/mnt/c" + self.runDir[2:]
        commandLine ="cd /mnt/c/Users/%username%/Documents/MicroSim/KKS_OpenCl/; python3 GEdata_writer.py " +runDirGE +"/"+self.infile.text() + " ;make clean;make; cp microsim_kks_opencl " + runDirGE +  ";cd " + runDirGE + ";mpirun -np 4 ./microsim_kks_opencl "  +self.infile.text()+" "+self.filling.text()+" "+self.output.text() 
        
        os.system(f'cmd /c start cmd /c wsl ~ -e sh -c "{commandLine }" ')


    elif self.radio_CH.isChecked():

        runDirGE = "/mnt/c" + self.runDir[2:]
        commandLine ="cd /mnt/c/Users/%username%/Documents/MicroSim/Cahn_Hilliard_FFT_2D/; python3 GEdata_writer.py " +runDirGE +"/"+self.infile.text() + " ;make clean;make; cp microsim_ch_fft " + runDirGE +  ";cd " + runDirGE + ";mpirun -np 4 ./microsim_ch_fft "  +self.infile.text()+" "+self.filling.text()+" "+self.output.text() 
        
        os.system(f'cmd /c start cmd /c wsl ~ -e sh -c "{commandLine }" ')


def SolverExecuteHelp(self):
    runHelpmsg = QMessageBox()
    runHelpmsg.setWindowTitle("Run Help")
    runDirGE = "/mnt/c" + self.runDir[2:]
    
    if self.radio_GP.isChecked():

        Model_Folder= "Grand_potential_Finite_difference_2D_MPI"
        Model_code ="microsim_gp"
    
    elif self.radio_KKR.isChecked():

        Model_Folder= "KKS_CuFFT"
        Model_code ="microsim_kks_cufft"
    

    elif self.radio_KKS2.isChecked():

        Model_Folder= "KKS_OpenCl"
        Model_code ="microsim_kks_opencl"

    elif self.radio_CH.isChecked():
        Model_Folder = "Cahn_Hilliard_FFT_2D"
        Model_code = "microsim_ch_fft"

    runcmdhelp = "This action will execute the following command.\n\n\n     1) Solver compilation\n\n       cd ~/MicroSim/"+Model_Folder+"/  \n \n      python3 GEdata_writer.py " +self.runDir +"/"+self.infile.text() + "\n\n      make clean\n\n      make\n\n\n     2) Solver execution\n\n      cp "+ Model_code+" ~/MicroSim/bin/\n\n      cd " + self.runDir + "\n\n      mpirun.mpich -np 4 ~/MicroSim/bin/"+ Model_code+" "  +self.infile.text()+" "+self.filling.text()+" "+self.output.text() + " 2 2"

    runHelpmsg.setText(runcmdhelp )
    runHelpmsg.setStyleSheet("QLabel{min-width: 650px;}");
    runHelpmsg.exec_()


def generateJobscript(self):
    os.system("copy " + self.runDir + "/" + self.infile.text() +" "+ self.runDir + "/JOB_FILE/" + self.infile.text())
    os.system("copy " + self.runDir + "/" + self.filling.text() +" "+ self.runDir + "/JOB_FILE/" + self.filling.text())


    if self.radio_GP.isChecked():
            
        commandLine ="cd /mnt/c/Users/%username%/Documents/MicroSim/Grand_potential_Finite_difference_2D_MPI/; python3 GEdata_writer.py " +self.runDir +"/"+self.infile.text() + " ;make clean;make; cp microsim_gp ~/MicroSim/bin/;cp microsim_gp "+ self.runDir +"/JOB_FILE/;cp reconstruct "+ self.runDir +"/JOB_FILE/;cp write_xdmf "+ self.runDir +"/JOB_FILE/"
        
        os.system("gnome-terminal -e 'bash -c  \""+commandLine+";bash\"'")
    
    elif self.radio_KKR.isChecked():
        commandLine ="cd /mnt/c/Users/%username%/Documents/MicroSim/KKS_CuFFT/; python3 GEdata_writer.py " +self.runDir +"/"+self.infile.text() + " ;make clean;make; cp microsim_kks_cufft ~/MicroSim/bin/;cp microsim_kks_cufft "+ self.runDir +"/JOB_FILE/;cp reconstruct "+ self.runDir +"/JOB_FILE/;cp write_xdmf "+ self.runDir +"/JOB_FILE/"
        
        os.system("gnome-terminal -e 'bash -c \""+commandLine+";bash\"'")

    elif self.radio_KKS2.isChecked():
        commandLine ="cd /mnt/c/Users/%username%/Documents/MicroSim/KKS_OpenCl/; python3 GEdata_writer.py " +self.runDir +"/"+self.infile.text() + " ;make clean;make; cp microsim_kks_opencl ~/MicroSim/bin/;cp microsim_kks_opencl "+ self.runDir +"/JOB_FILE/;cp reconstruct "+ self.runDir +"/JOB_FILE/;cp write_xdmf "+ self.runDir +"/JOB_FILE/"
        
        os.system("gnome-terminal -e 'bash -c  \""+commandLine+";bash\"'")

    elif self.radio_CH.isChecked():
        commandLine ="cd /mnt/c/Users/%username%/Documents/MicroSim/Cahn_Hilliard_FFT_2D/; python3 GEdata_writer.py " +self.runDir +"/"+self.infile.text() + " ;make clean;make; cp cp microsim_ch_fft ~/MicroSim/bin/;cp microsim_ch_fft "+ self.runDir +"/JOB_FILE/;cp reconstruct "+ self.runDir +"/JOB_FILE/;cp write_xdmf "+ self.runDir +"/JOB_FILE/"
        
        os.system("gnome-terminal -e 'bash -c \""+commandLine+";bash\"'")