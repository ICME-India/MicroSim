#ifndef CHKPNT_H_
#define CHKPNT_H_

#include <AMReX_VisMF.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_VisMF.H>
#include "Variables.H"

using namespace std;
using namespace amrex;


void SetBA(int lev, const BoxArray& ba_in)
{
    grids.resize(lev+1);
    if (grids[lev] != ba_in) grids[lev] = ba_in;
}


void GotoNextLine (std::istream& is)
{
    constexpr std::streamsize bl_ignore_max { 100000 };
    is.ignore(bl_ignore_max, '\n');
}



void Writechkfile (amrex::MultiFab& phi_new, amrex::MultiFab& mu_new, int n, std::string chk_file)
{
	BL_PROFILE("Writechkfile()");       

 const std::string& checkpointname = amrex::Concatenate(chk_file,n,1);

        //Print()<<"n is: "<<n<<"\n";

        amrex::Print() << "Writing checkpoint " << checkpointname << "\n";

        const int nlevels = 1;

        amrex::PreBuildDirectorHierarchy(checkpointname, "Level_", nlevels, true);

        if (ParallelDescriptor::IOProcessor()) {

       std::string HeaderFileName(checkpointname + "/Header");
       VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
       std::ofstream HeaderFile;
       HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
       HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                                               std::ofstream::trunc |
                                               std::ofstream::binary);
       if( ! HeaderFile.good()) {
           amrex::FileOpenFailed(HeaderFileName);
       }

       
       HeaderFile.precision(17);

       HeaderFile << "Checkpoint file for heat eq 1\n";

       HeaderFile <<  n << "\n";

       HeaderFile << timee << "\n";
       
       HeaderFile << ParallelDescriptor::second()-strt_time << "\n";
       
       (phi_new.boxArray()).writeOn(HeaderFile);
       HeaderFile << "\n";


        }

        VisMF::Write(phi_new,
                    amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "phi"));
        
        VisMF::Write(mu_new,
                    amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "mu"));
}



void Readchkfile (MultiFab& phi_new, MultiFab& mu_new)
{   
    amrex::Print() << "Restart from checkpoint " << restart_chkfile << "\n";

    std::string File(restart_chkfile + "/Header");

    VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);

    std::string line, word;

    GotoNextLine(is);

    std::getline(is, line);
    {
        std::istringstream lis(line);
        //int i = 0;
        while (lis >> word) {
            stepnum = std::stoi(word);
        }
    }

    std::getline(is, line);
    {
        std::istringstream lis(line);
        //int i = 0;
        while (lis >> word) {
            timee = std::stod(word);
        }
    }
    
    std::getline(is, line);
    {
        std::istringstream lis(line);
        //int i = 0;
        while (lis >> word) {
            rst_time = std::stod(word);
        }
    }

    BoxArray ba;
        ba.readFrom(is);
        GotoNextLine(is);
        int lev = 0;

        SetBA(lev, ba);
        
        VisMF::Read(phi_new,
                    amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "phi"));
        
        VisMF::Read(mu_new,
                    amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "mu"));
    
}



#endif
