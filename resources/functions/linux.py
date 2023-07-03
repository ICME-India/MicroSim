import sys, os, glob
from xml.dom import minidom
import numpy as np
import yt
from PyQt5.QtWidgets import QMessageBox

def H5toVTK_Func(self):
    if self.h5tovtk_h5radio.isChecked():
            
        h5_outhead, h5_outtail = os.path.split(self.h5tovtk_outputloc.text())
        
        split_text = "_" + self.h5tovtk_sTime.text()+".h5"
        
        h5_outputfilename = h5_outtail.split(split_text)
        
        #reading vtk output file dir
        
        if self.h5tovtk_vtkoutput.text() == "":
            vtk_output_fname = h5_outhead + "/" + h5_outputfilename[0]+ ".vtk"
        
        elif os.path.isdir(self.h5tovtk_vtkoutput.text()):
            vtk_output_fname = self.h5tovtk_vtkoutput.text() + "/" + h5_outputfilename[0]+ ".vtk"
            
        else:
            return False
            
        
        h5toxmf_cmd = "cd " +h5_outhead + "; cd ..;" + "cp ../Grand_potential_Finite_difference_2D_MPI/write_xdmf write_xdmf ; ./write_xdmf " + self.h5tovtk_infileLoc.text() + " "+ h5_outputfilename[0] + " " +  self.h5tovtk_sTime.text() + " "+ self.h5tovtk_eTime.text() 
        

        os.system("gnome-terminal -e 'bash -c \""+h5toxmf_cmd+";\"'")
        
        
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
        
        xmftovtk_cmd = "cd " +h5_outhead + ";pvpython h5tovtk.py;rm h5tovtk.py;cd ..;rm write_xdmf"
        os.system("gnome-terminal -e 'bash -c \""+xmftovtk_cmd+";\"'")
        
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
        os.system("gnome-terminal -e 'bash -c \""+xmftovtk_cmd+";\"'")

    elif self.h5tovtk_pltradio.isChecked():
        
        Zeroth_file = (self.h5tovtk_outputloc.text())
            

        header_dir = Zeroth_file.split('plt')[0] + 'plt'

        list_of_files =list(filter( os.path.isdir, glob.glob(header_dir +  '*') ) )

        list_of_files2 =np.copy(list_of_files)

        Int_list = [int(i.split('plt')[1]) for i in list_of_files2]

        Shorted_dir =[x for _, x in sorted(zip(Int_list, list_of_files))]

        for i in range(int(self.h5tovtk_sTime.text()), len(Shorted_dir)  , int(self.h5tovtk_eskip.text())):
           
            filename =Shorted_dir[i]
            ds = yt.load(filename)
            Scaler_fields = []
            Grid_data =  ds.index.select_grids(ds.index.max_level)
            Grid_dimension = int(np.sqrt(Grid_data.shape[0]))

            
            
            if self.h5tovtk_vtkoutput.text() == "":
                f = open( header_dir +'_'+str(i) +'.vtk','w') # change your vtk file name

            else:
                f = open(  self.h5tovtk_vtkoutput.text() + '/plt_'+str(i) +'.vtk','w') # change your vtk file name

            self.h5tovtk_error.setText("Saving : plt_ " + str(i) + ".vtk" )
            self.h5tovtk_error.repaint()

            f.write('# vtk DataFile Version 2.0\n')
            f.write(filename+'\n')
            f.write('ASCII\n')
            f.write('DATASET STRUCTURED_POINTS\n')
            f.write('DIMENSIONS '+ str(ds.domain_dimensions[0])+ ' '+ str(ds.domain_dimensions[1]) +' '+  str(ds.domain_dimensions[2]) +'\n')
            f.write('SPACING 0.100000 0.100000 0.100000\n') # task - find dx
            f.write('ORIGIN 0 0 0\n')
            f.write('POINT_DATA ' + str(ds.domain_dimensions[0]*ds.domain_dimensions[1]*ds.domain_dimensions[2]) +'\n')

            for i in range( len(ds.field_list)):
                scalar_field = ds.field_list[i][1]
                f.write('SCALARS '+ scalar_field +' float 1\n')
                f.write('LOOKUP_TABLE default\n')

                #merging all blocks
                row_data = []
                for i in range(Grid_dimension):
                    vs_data = Grid_data[int(Grid_dimension*i)]
                    data_vs = vs_data[scalar_field][:,:,0]

                    for j in range(1,Grid_dimension):
                        gs_data = Grid_data[j+ int(Grid_dimension*i) ]
                        data2 = gs_data[scalar_field][:,:,0]
                        data_vs = np.concatenate((data_vs,data2),axis=0)

                    row_data.append(data_vs)
                full_data = np.concatenate(row_data,axis=1)
                np.savetxt(f, full_data, fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ', encoding=None)
            f.close()
        
        self.h5tovtk_infileLoc.setText("")
        self.h5tovtk_infileLoc.setText("")
        self.h5tovtk_infileLoc.setText("")
        self.h5tovtk_infileLoc.setText("")

        self.h5tovtk_error.setText("Done" )
        self.h5tovtk_error.repaint()


        
def paraviewFunc(self):
    if os.path.isfile("/opt/paraviewopenfoam56/bin/paraview"):  ##checking paraview installation status
        list_of_files = glob.glob(self.runDir + "/DATA/"+ self.output.text() +"*.*")
        if len(list_of_files) == 0:
            self.finish_error.setText("Sorry, output file not found.")
            return
        latest_file = max(list_of_files, key=os.path.getctime)

        paraviewcmd = "gnome-terminal -e 'bash -c \"/opt/paraviewopenfoam56/bin/paraview " +latest_file +"; bash\" '"
        os.system(paraviewcmd)
    elif os.path.isfile(os.path.expanduser("./.Paraview")): ## Checking for paraview saved path
        readPath = open(os.path.expanduser('./.Paraview'), "r")
        readPathParaview = readPath.read().replace("\n", "")
        list_of_files = glob.glob(self.runDir + "/DATA/"+ self.output.text() +"*.*")
        if len(list_of_files) == 0:
            self.finish_error.setText("Sorry, output file not found.")
            return
        latest_file = max(list_of_files, key=os.path.getctime)

        paraviewcmd = "gnome-terminal -e 'bash -c \"" + readPathParaview +" " +latest_file +"; bash\" '"
        os.system(paraviewcmd)

    else:
        self.paraviewError.show()
        
def SolverExecute(self):
    commandLine ="mkdir -p ~/MicroSim/bin"
        
    os.system("gnome-terminal -e 'bash -c  \""+commandLine+";bash\"'")

    if self.radio_GP.isChecked():
            
        commandLine ="cd Grand_potential_Finite_difference_2D_MPI/; python3 GEdata_writer.py " +self.runDir +"/"+self.infile.text() + " ;make clean;make; cp microsim_gp ~/MicroSim/bin/;cd " + self.runDir + ";mpirun.mpich -np 4 ~/MicroSim/bin/microsim_gp "  +self.infile.text()+" "+self.filling.text()+" "+self.output.text() + " 2 2"
        
        os.system("gnome-terminal -e 'bash -c  \""+commandLine+";bash\"'")

    elif self.radio_of.isChecked():
            
        commandLine ="cd " + self.runDir + "; cd ../../solver ; wclean; wmake; cd " + self.runDir + " ; ./Allclean; ./Allrun "
        
        os.system("gnome-terminal -e 'bash -c  \""+commandLine+";bash\"'")

    elif self.radio_amrex.isChecked():
            
        commandLine ="cd Grand_potential_AMReX/Exec; make clean; make;  g++ -o Replace Replace.cpp; ./Replace "  +self.infile.text()+" "+self.filling.text()+"; cp main2d.gnu.MPI.ex ~/MicroSim/bin/;cd " + self.runDir + ";mpirun -np 2  ~/MicroSim/bin/main2d.gnu.MPI.ex input2.in"
        
        os.system("gnome-terminal -e 'bash -c  \""+commandLine+";bash\"'")
    
    elif self.radio_KKR.isChecked():
        commandLine ="cd KKS_FD_CUDA_MPI/; python3 GEdata_writer.py " +self.runDir +"/"+self.infile.text() + " ;make clean;make; cp microsim_kks_fd_cuda_mpi ~/MicroSim/bin/;cd " + self.runDir + "; make run KKS_OUT=~/MicroSim/bin/microsim_kks_fd_cuda_mpi INPUT="+self.infile.text()+" FILLING="+self.filling.text()+" OUTPUT="+self.output.text()+" NPROCS=1"
        
        os.system("gnome-terminal -e 'bash -c \""+commandLine+";bash\"'")

    elif self.radio_KKS2.isChecked():
        commandLine ="cd KKS_OpenCl/; python3 GEdata_writer.py " +self.runDir +"/"+self.infile.text() + " ;make clean;make; cp microsim_kks_opencl ~/MicroSim/bin/;cd " + self.runDir + ";~/MicroSim/bin/microsim_kks_opencl "  +self.infile.text()+" "+self.filling.text()+" "+self.output.text()
        
        os.system("gnome-terminal -e 'bash -c  \""+commandLine+";bash\"'")

    elif self.radio_CH.isChecked():
        commandLine ="cd Cahn_Hilliard_FFT_2D/; python3 GEdata_writer.py " +self.runDir +"/"+self.infile.text() + " ;make clean;make; cp microsim_ch_fft ~/MicroSim/bin/;cd " + self.runDir + ";~/MicroSim/bin/microsim_ch_fft "  +self.infile.text()+" "+self.filling.text()+" "+self.output.text()
        
        os.system("gnome-terminal -e 'bash -c \""+commandLine+";bash\"'")

def SolverExecuteHelp(self):
    runHelpmsg = QMessageBox()
    runHelpmsg.setWindowTitle("Run Help")
    
    if self.radio_GP.isChecked():

        Model_Folder= "Grand_potential_Finite_difference_2D_MPI"
        Model_code ="microsim_gp"
        runcmdhelp = "This action will execute the following command:\n\n\n     1) Solver compilation\n\n       cd MicroSim/"+Model_Folder+"/  \n \n      python3 GEdata_writer.py " +self.runDir +"/"+self.infile.text() + "\n\n      make clean\n\n      make\n\n\n     2) Solver execution\n\n      cp "+ Model_code+" ~/MicroSim/bin/\n\n      cd " + self.runDir + "\n\n      mpirun.mpich -np 4 ~/MicroSim/bin/"+ Model_code+" "  +self.infile.text()+" "+self.filling.text()+" "+self.output.text() + " 2 2"
        
    elif self.radio_of.isChecked():

        Model_Folder= "Grand_potential_OpenFOAM"
        Model_code ="Allrun"
        runcmdhelp = "This action will execute the following command:\n\n\n     1) Solver compilation\n\n       cd MicroSim/"+Model_Folder+"/solver\n\n      wclean\n\n      wmake\n\n\n     2) Solver execution\n\n      ./Allclean\n\n      ./Allrun"
        
    elif self.radio_amrex.isChecked():

        Model_Folder= "Grand_potential_AMReX"
        Model_code ="main2d.gnu.MPI.ex"
        runcmdhelp = "This action will execute the following command:\n\n\n     1) Solver compilation\n\n       cd MicroSim/"+Model_Folder+"/  \n \n            make clean\n\n      make\n      g++ -o Replace Replace.cpp\n      ./Replace "  +self.infile.text()+" "+self.filling.text()+" \n\n\n     2) Solver execution\n\n      cp "+ Model_code+" ~/MicroSim/bin/\n\n      cd " + self.runDir + "\n\n      mpirun -np 4 ~/MicroSim/bin/main2d.gnu.MPI.ex input2.in"
    
    elif self.radio_KKR.isChecked():

        Model_Folder= "KKS_FD_CUDA_MPI"
        Model_code ="microsim_kks_fd_cuda_mpi"
        runcmdhelp = "This action will execute the following command:\n\n\n     1) Solver compilation\n\n       cd MicroSim/"+Model_Folder+"/  \n \n      python3 GEdata_writer.py " +self.runDir +"/"+self.infile.text() + "\n\n      make clean\n\n      make\n\n\n     2) Solver execution\n\n      cp "+ Model_code+" ~/MicroSim/bin/\n\n      cd " + self.runDir + "\n\n      mpirun -np 4 ~/MicroSim/bin/"+ Model_code+" "  +self.infile.text()+" "+self.filling.text()+" "+self.output.text()
    

    elif self.radio_KKS2.isChecked():

        Model_Folder= "KKS_OpenCl"
        Model_code ="microsim_kks_opencl"
        runcmdhelp = "This action will execute the following command:\n\n\n     1) Solver compilation\n\n       cd MicroSim/"+Model_Folder+"/  \n \n      python3 GEdata_writer.py " +self.runDir +"/"+self.infile.text() + "\n\n      make clean\n\n      make\n\n\n     2) Solver execution\n\n      cp "+ Model_code+" ~/MicroSim/bin/\n\n      cd " + self.runDir + "\n\n      ~/MicroSim/bin/"+ Model_code+" "  +self.infile.text()+" "+self.filling.text()+" "+self.output.text()

    elif self.radio_CH.isChecked():
        Model_Folder = "Cahn_Hilliard_FFT_2D"
        Model_code = "microsim_ch_fft"
        runcmdhelp = "This action will execute the following command:\n\n\n     1) Solver compilation\n\n       cd MicroSim/"+Model_Folder+"/  \n \n      python3 GEdata_writer.py " +self.runDir +"/"+self.infile.text() + "\n\n      make clean\n\n      make\n\n\n     2) Solver execution\n\n      cp "+ Model_code+" ~/MicroSim/bin/\n\n      cd " + self.runDir + "\n\n      ~/MicroSim/bin/"+ Model_code+" "  +self.infile.text()+" "+self.filling.text()+" "+self.output.text()

    #runcmdhelp = "This action will execute the following command.\n\n\n     1) Solver compilation\n\n       cd ~/MicroSim/"+Model_Folder+"/  \n \n      python3 GEdata_writer.py " +self.runDir +"/"+self.infile.text() + "\n\n      make clean\n\n      make\n\n\n     2) Solver execution\n\n      cp "+ Model_code+" ~/MicroSim/bin/\n\n      cd " + self.runDir + "\n\n      mpirun.mpich -np 4 ~/MicroSim/bin/"+ Model_code+" "  +self.infile.text()+" "+self.filling.text()+" "+self.output.text() + " 2 2"

    runHelpmsg.setText(runcmdhelp )
    runHelpmsg.setStyleSheet("QLabel{min-width: 650px;}");
    runHelpmsg.exec_()


def generateJobscript(self):
    os.system("cp " + self.runDir + "/" + self.infile.text() +" "+ self.runDir + "/JOB_FILE/" + self.infile.text())
    os.system("cp " + self.runDir + "/" + self.filling.text() +" "+ self.runDir + "/JOB_FILE/" + self.filling.text())


    if self.radio_GP.isChecked():
            
        commandLine ="cd Grand_potential_Finite_difference_2D_MPI/; python3 GEdata_writer.py " +self.runDir +"/"+self.infile.text() + " ;make clean;make; cp microsim_gp ~/MicroSim/bin/;cp microsim_gp "+ self.runDir +"/JOB_FILE/;cp reconstruct "+ self.runDir +"/JOB_FILE/;cp write_xdmf "+ self.runDir +"/JOB_FILE/"
        
        os.system("gnome-terminal -e 'bash -c  \""+commandLine+";bash\"'")
    
    elif self.radio_of.isChecked():
            
        commandLine ="cd " + self.runDir + "; cd ../../solver ; wclean; wmake"
        
        os.system("gnome-terminal -e 'bash -c  \""+commandLine+";bash\"'")
        
    elif self.radio_amrex.isChecked():
            
        commandLine ="cd Grand_potential_AMReX/Exec; make clean; make;  g++ -o Replace Replace.cpp; ./Replace "  +self.infile.text()+" "+self.filling.text()+"; cp main2d.gnu.MPI.ex ~/MicroSim/bin/;cp main2d.gnu.MPI.ex "+ self.runDir +"/JOB_FILE/;cp input2.in "+ self.runDir +"/JOB_FILE/"
        
        os.system("gnome-terminal -e 'bash -c  \""+commandLine+";bash\"'")
    
    elif self.radio_KKR.isChecked():
        commandLine ="cd KKS_FD_CUDA_MPI/; python3 GEdata_writer.py " +self.runDir +"/"+self.infile.text() + " ;make clean;make; cp microsim_kks_fd_cuda_mpi ~/MicroSim/bin/;cp microsim_kks_fd_cuda_mpi "+ self.runDir +"/JOB_FILE/;cp reconstruct "+ self.runDir +"/JOB_FILE/;cp write_xdmf "+ self.runDir +"/JOB_FILE/"
        
        os.system("gnome-terminal -e 'bash -c \""+commandLine+";bash\"'")

    elif self.radio_KKS2.isChecked():
        commandLine ="cd KKS_OpenCl/; python3 GEdata_writer.py " +self.runDir +"/"+self.infile.text() + " ;make clean;make; cp microsim_kks_opencl ~/MicroSim/bin/;cp microsim_kks_opencl "+ self.runDir +"/JOB_FILE/;cp reconstruct "+ self.runDir +"/JOB_FILE/;cp write_xdmf "+ self.runDir +"/JOB_FILE/"
        
        os.system("gnome-terminal -e 'bash -c  \""+commandLine+";bash\"'")

    elif self.radio_CH.isChecked():
        commandLine ="cd Cahn_Hilliard_FFT_2D/; python3 GEdata_writer.py " +self.runDir +"/"+self.infile.text() + " ;make clean;make; cp microsim_ch_fft ~/MicroSim/bin/;cp microsim_ch_fft "+ self.runDir +"/JOB_FILE/;cp reconstruct "+ self.runDir +"/JOB_FILE/;cp write_xdmf "+ self.runDir +"/JOB_FILE/"
        
        os.system("gnome-terminal -e 'bash -c \""+commandLine+";bash\"'")
