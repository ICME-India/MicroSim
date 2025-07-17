#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_Utility.H>
#include <AMReX_Print.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include <AmrCoreGP.H>
#include "CalculateGrad.H"
#include "Boundary_conditions.h"


using namespace amrex;
using namespace std;


//####################################################################################################################

AmrCoreGP::~AmrCoreGP() 
{
    // Destructor implementation
}

//####################################################################################################################

AmrCoreGP::AmrCoreGP()
{   
    ReadInput();
    int nlevs_max = max_level + 1;

    istep.resize(nlevs_max, 0);
    dtlev.resize(nlevs_max,0);
    nsubsteps.resize(nlevs_max, 1);

    phi_new.resize(nlevs_max);
    phi_old.resize(nlevs_max);

    mu_new.resize(nlevs_max);
    mu_old.resize(nlevs_max);

    comp_old.resize(nlevs_max);
    comp_new.resize(nlevs_max);

    term1.resize(nlevs_max);
    term2.resize(nlevs_max);
    term3.resize(nlevs_max);
    term4.resize(nlevs_max);

    psi.resize(nlevs_max);
    lambad.resize(nlevs_max);

    Grad.resize(AMREX_SPACEDIM);

    plot.resize(nlevs_max);

    for (int lev = 1; lev <= max_level; ++lev) 
    {
        nsubsteps[lev] = MaxRefRatio(lev-1);
    }

    bcs_phi.resize(nump); 
    bcs_mu.resize(numcom-1); 
    bcs_comp.resize(numcom-1); 
	
    Vector<int> bc_lo_phi(AMREX_SPACEDIM,0);
    Vector<int> bc_lo_mu(AMREX_SPACEDIM,0);
    Vector<int> bc_lo_comp(AMREX_SPACEDIM,0);

    Vector<int> bc_hi_phi(AMREX_SPACEDIM,0);
    Vector<int> bc_hi_mu(AMREX_SPACEDIM,0);
    Vector<int> bc_hi_comp(AMREX_SPACEDIM,0);

    for(int w=0; w<bound.size(); w++)
    {    
        if(bound[w][0]=="phi")
        {   
            bc_hi_phi[0]=stod(bound[w][1]);
            bc_lo_phi[0]=stod(bound[w][2]);
            bc_hi_phi[1]=stod(bound[w][3]);
            bc_lo_phi[1]=stod(bound[w][4]);
            #if(AMREX_SPACEDIM>2)
			bc_hi_phi[2]=stod(bound[w][5]);
			bc_lo_phi[2]=stod(bound[w][6]);
			#endif
        }
        else if(bound[w][0]=="mu")
        {   
            bc_hi_mu[0]=stod(bound[w][1]);
            bc_lo_mu[0]=stod(bound[w][2]);
            bc_hi_mu[1]=stod(bound[w][3]);
            bc_lo_mu[1]=stod(bound[w][4]);
            #if(AMREX_SPACEDIM>2)
			bc_hi_mu[2] = stod(bound[w][5]);
			bc_lo_mu[2] = stod(bound[w][6]);
			#endif
        }
        else if(bound[w][0]=="c")
        {   
            bc_hi_comp[0]=stod(bound[w][1]);
            bc_lo_comp[0]=stod(bound[w][2]);
            bc_hi_comp[1]=stod(bound[w][3]);
            bc_lo_comp[1]=stod(bound[w][4]);
            #if(AMREX_SPACEDIM>2)
			bc_hi_comp[2] = stod(bound[w][5]);
			bc_lo_comp[2] = stod(bound[w][6]);
			#endif
        }
    }
    bound_cond(bcs_phi,  bc_hi_phi,  bc_lo_phi,nump); 
	bound_cond(bcs_mu,   bc_hi_mu,   bc_lo_mu,numcom-1); 
	bound_cond(bcs_comp, bc_hi_comp, bc_lo_comp,numcom-1);
}

//####################################################################################################################

amrex::Real 
AmrCoreGP::FunctionTau(int i, int j, int k, int numphase, Real Tauu, amrex::Array4<Real const> const& phi)
{
    double sum=0.0, sum1=0.0;
    long a, b;
    for (a=0; a<numphase; a++) 
    {
        for (b=0; b<numphase; b++) 
        {
            if (a<b)
            {
                sum  += Tauu*phi(i,j,k,a)*phi(i,j,k,b);
                sum1 += phi(i,j,k,a)*phi(i,j,k,b);
            }
         }
    }
    if (sum1) 
    {
        return sum/sum1;
    } 
    else 
    {
        return Tauu;
    }
}

//####################################################################################################################

void 
AmrCoreGP::Calculate_Tau()
{ 
    Real min_tau{0.0};
    Vector<Real> deltac(numcom-1,0.0);
    Vector<Real> deltamu(numcom-1,0.0);
    Vector<Vector<Real>>prod(numcom-1,Vector<Real>(numcom-1,0.0));
    Vector<Vector<Real>>inv_dcdmu(numcom-1,Vector<Real>(numcom-1,0.0));
    Vector<Vector<Real>> tau_ab(nump,Vector<Real>(nump,0.0));
    for(int a=0; a<nump-1; a++)
    {
        for(int k=0; k<numcom-1; k++)
        {
            deltac[k] = conceq[nump-1][k]-conceq[a][k];
        }
        if(numcom ==2)
        {
            prod[0][0] = dcdmu[nump-1][0][0]*diff[nump-1][0][0];
            inv_dcdmu[0][0] = 1/prod[0][0];
            deltamu[0] = inv_dcdmu[0][0]*deltac[0];
        }

        if(numcom ==3)
        {
            prod[0][0] = dcdmu[nump-1][0][0]*diff[nump-1][0][0] + dcdmu[nump-1][0][1]*diff[nump-1][1][0];
            prod[0][1] = dcdmu[nump-1][0][1]*diff[nump-1][0][0] + dcdmu[nump-1][1][1]*diff[nump-1][0][1];
            prod[1][0] = dcdmu[nump-1][0][0]*diff[nump-1][1][0] + dcdmu[nump-1][1][0]*diff[nump-1][1][1];
            prod[1][1] = dcdmu[nump-1][1][0]*diff[nump-1][0][1] + dcdmu[nump-1][1][1]*diff[nump-1][1][1];
            Real det = prod[0][0]*prod[1][1] - prod[0][1]*prod[1][0];
            inv_dcdmu[0][0] = prod[1][1]/det;
            inv_dcdmu[0][1] = -prod[0][1]/det;
            inv_dcdmu[1][0] = -prod[1][0]/det;
            inv_dcdmu[1][1] = prod[0][0]/det; 
            deltamu[0] = inv_dcdmu[0][0]*deltac[0] + inv_dcdmu[0][1]*deltac[1];
            deltamu[1] = inv_dcdmu[1][0]*deltac[0] + inv_dcdmu[1][1]*deltac[1]; 
        }

        if(numcom>3)
        {
            int sz = numcom-1;
            mat_mul2D(dcdmu,diff,prod,sz);
            Array2D<Real,0,compcount-2,0,compcount-2, Order::C> prd{};
            Array2D<Real,0,compcount-2,0,compcount-2, Order::C> inv_prd{};
            for(int i=0; i<prod.size();i++)
            {
                for(int j=0; j<prod[0].size();j++)
                {
                    prd(i,j) = prod[i][j];
                }
            }
            mat_inv(prd,inv_prd,numcom);
            for(int i=0; i<prod.size();i++)
            {
                for(int j=0; j<prod[0].size();j++)
                {
                    inv_dcdmu[i][j] = inv_prd(i,j);
                }
            }
            mat_mul1D(inv_dcdmu,deltac,deltamu,sz);
        }
        double sum=0.0;
        for (int k=0; k<numcom-1; k++) 
        {
            sum += deltamu[k]*deltac[k];
        }
        tau_ab[a][nump-1] = sum*eps*(0.2222)/Vm;
        tau_ab[nump-1][a] = tau_ab[a][nump-1];
        if (a==0) 
        {
            min_tau = tau_ab[a][nump-1];
        }
        if (tau_ab[a][nump-1] < min_tau) 
        {
            min_tau = tau_ab[a][nump-1];
        }
        deltac.clear();
        deltamu.clear();
        prod.clear();
        inv_dcdmu.clear();
        deltac = Vector<Real>(numcom-1,0.0);
        deltamu = Vector<Real>(numcom-1,0.0);
        prod = Vector<Vector<Real>>(numcom-1,Vector<Real>(numcom-1,0.0));
        inv_dcdmu = Vector<Vector<Real>>(numcom-1,Vector<Real>(numcom-1,0.0));
    }
    for (int a=0; a<nump; a++) 
    {
        for (int b=0; b<nump; b++) 
        {
            tau_ab[a][b] = min_tau;
        }
    }
    tau = min_tau;  
}

//####################################################################################################################

void 
AmrCoreGP::MakeNewLevelFromScratch (int lev, Real time, const BoxArray& ba,const DistributionMapping& dm)
{   const int ncomp_phi = nump;
    const int ncomp_mu = numcom-1;
    const int ncomp_comp = numcom-1;
    const int nghost = 0;

    
    phi_new[lev].define(ba,dm,ncomp_phi,nghost);
    phi_old[lev].define(ba,dm,ncomp_phi,nghost);
    mu_new[lev].define(ba,dm,ncomp_mu,nghost);
    mu_old[lev].define(ba,dm,ncomp_mu,nghost);
    comp_new[lev].define(ba,dm,ncomp_comp,nghost);
    comp_old[lev].define(ba,dm,ncomp_comp,nghost);
    term1[lev].define(ba,dm,ncomp_phi,nghost);
    term2[lev].define(ba,dm,ncomp_phi,nghost);
    term3[lev].define(ba,dm,ncomp_phi,nghost);
    term4[lev].define(ba,dm,ncomp_phi,nghost);
    psi[lev].define(ba,dm,ncomp_phi,nghost);
    lambad[lev].define(ba,dm,ncomp_phi,nghost);
  
    MultiFab& phi_state = phi_new[lev];
    MultiFab& mu_state = mu_new[lev];
    MultiFab& comp_state = comp_new[lev];
    phi_state.setVal(0);
    mu_state.setVal(0);
    comp_state.setVal(0);

    if (lev==0)
    {
        if(cylinder.size()>0)
        {
            Init_phi_cyl(phi_state);
        }
        if(sphere.size()>0)
        {
            Init_phi_sph(phi_state);
        }
        
        Init_liq(phi_state);
        if(funcf == 4)
            {
                Init_mu(mu_state,phi_state);
            }
        Init_comp(phi_state,comp_state);
    }
   
    if (lev>0)
    {   
        FillCoarsePatch ("P", lev, time, phi_state, 0, ncomp_phi);
        FillCoarsePatch ("M", lev, time, mu_state, 0, ncomp_mu);
        FillCoarsePatch ("C", lev, time, comp_state, 0, ncomp_comp);
    };
}

//####################################################################################################################

void
AmrCoreGP::MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba,const DistributionMapping& dm)
{   
    // Not required for this example.
}

//####################################################################################################################

void
AmrCoreGP::ClearLevel (int lev)
{   
    // Not required for this example.
}

//####################################################################################################################

void
AmrCoreGP::RemakeLevel (int lev, Real time, const BoxArray& ba,const DistributionMapping& dm)
{  
    const int ncomp_phi  = phi_new[lev-1].nComp();
    const int ncomp_mu  = mu_new[lev-1].nComp();
    const int ncomp_comp  = comp_new[lev-1].nComp();
    const int nghost  = phi_new[lev-1].nGrow();
  
    MultiFab new_phi(ba,dm,ncomp_phi,nghost);
    MultiFab old_phi(ba,dm,ncomp_phi,nghost);
    MultiFab new_mu(ba,dm,ncomp_mu,nghost);
    MultiFab old_mu(ba,dm,ncomp_mu,nghost);
    MultiFab new_comp(ba,dm,ncomp_comp,nghost);
    MultiFab old_comp(ba,dm,ncomp_comp,nghost);
    MultiFab state_term1(ba,dm,ncomp_phi,nghost);
    MultiFab state_term2(ba,dm,ncomp_phi,nghost);
    MultiFab state_term3(ba,dm,ncomp_phi,nghost);
    MultiFab state_term4(ba,dm,ncomp_phi,nghost);
    MultiFab state_lambad(ba,dm,ncomp_phi,nghost);
    MultiFab state_psi(ba,dm,ncomp_phi,nghost);
   
    FillPatch("P", lev, time, new_phi, 0, ncomp_phi);
    FillPatch("M", lev, time, new_mu, 0, ncomp_mu);
    FillPatch("C", lev, time, new_comp, 0, ncomp_comp);

    std::swap(new_phi, phi_new[lev]);
    std::swap(old_phi, phi_old[lev]);
    std::swap(new_mu, mu_new[lev]);
    std::swap(old_mu, mu_old[lev]);
    std::swap(new_comp, comp_new[lev]);
    std::swap(old_comp, comp_old[lev]);
    std::swap(state_term1, term1[lev]);
    std::swap(state_term2, term2[lev]);
    std::swap(state_term3, term3[lev]);
    std::swap(state_term4, term4[lev]);
    std::swap(state_psi, psi[lev]);
    std::swap(state_lambad, lambad[lev]);
}

//####################################################################################################################

void
AmrCoreGP::CalculateGrad(amrex::Real t, int lev)
{
    MultiFab state1(grids[lev], dmap[lev], phi_new[lev].nComp(),1);
    FillPatch("P",lev, t, state1 , 0, state1.nComp());
    for (int i=0;i<AMREX_SPACEDIM;i++)
        {
            Grad[i].define(grids[lev],dmap[lev],1,0);
        }    
    
    for (MFIter mfi(phi_new[lev]); mfi.isValid(); ++mfi)
        {   const Box& bx      = mfi.tilebox();
            const Array4<Real>& statefab = state1.array(mfi);
            const Array4<Real>& grad_x   = Grad[0].array(mfi);
            const Array4<Real>& grad_y   = Grad[1].array(mfi);
            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {   
                Gradient_x(i, j, k, statefab, grad_x);
            });

            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {   
                Gradient_y(i, j, k, statefab, grad_y);
            });

#if(AMREX_SPACEDIM > 2)
            const Array4<Real>& grad_z   = Grad[2].array(mfi);
            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {   
                Gradient_z(i, j, k, statefab, grad_z);
            });
#endif
        }
}

//####################################################################################################################

void
AmrCoreGP::ErrorEst (int lev, TagBoxArray& tags, Real /*time*/, int /*ngrow*/)
{   
    const int tagval  = TagBox::SET;
    CalculateGrad(time, lev);
    for (MFIter mfi(phi_new[lev]); mfi.isValid(); ++mfi)
    {   
        const Box& bx       = mfi.validbox();
        const auto gradx    = Grad[0].array(mfi);
        const auto grady    = Grad[1].array(mfi);
        const auto tagfab   = tags.array(mfi);

        #if(AMREX_SPACEDIM == 2)
        amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {   
            if ((gradx(i,j,k)>0.1)||(grady(i,j,k))>0.1)
            {
                tagfab(i,j,k) = tagval;
            }
        });
        #else
        const auto gradz    = Grad[2].array(mfi);
        amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {   
            if ((gradx(i,j,k)>0.1)||(grady(i,j,k))>0.1 || (gradz(i,j,k)>0.1))
            {
                tagfab(i,j,k) = tagval;
            }
        });
        #endif
    }
}

//####################################################################################################################

void
AmrCoreGP::FillPatch (std::string ref,int lev, Real time, MultiFab& mf, int icomp, int ncomp)
{   
    if(ref=="M")
    {   
        bcs.resize(numcom);
        bcs = bcs_mu;
    }
    else if(ref=="C")
    {   
        bcs.resize(numcom);
        bcs = bcs_comp;
    }
    else
    {   
        bcs.resize(nump);
        bcs = bcs_phi;
    }
    if (lev == 0)
    {   Vector<MultiFab*> smf;
        Vector<Real> stime;
        GetData(ref,0, time, smf, stime);
        if(Gpu::inLaunchRegion())
        {   GpuBndryFuncFab<AmrCoreFill> gpu_bndry_func(AmrCoreFill{});
            PhysBCFunct<GpuBndryFuncFab<AmrCoreFill> > physbc(geom[lev],bcs,gpu_bndry_func);
            amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, icomp, ncomp,geom[lev], physbc, 0);
        }
        else
        {   CpuBndryFuncFab bndry_func(nullptr);  // Without EXT_DIR, we can pass a nullptr.
            PhysBCFunct<CpuBndryFuncFab> physbc(geom[lev],bcs,bndry_func);
            amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, icomp, ncomp,geom[lev], physbc, 0);
        }
    }
    else
    {   
        Vector<MultiFab*> cmf, fmf;
        Vector<Real> ctime, ftime;
        GetData(ref,lev-1, time, cmf, ctime);
        GetData(ref,lev  , time, fmf, ftime);
        Interpolater* mapper = &lincc_interp;
        if(Gpu::inLaunchRegion())
        {   GpuBndryFuncFab<AmrCoreFill> gpu_bndry_func(AmrCoreFill{});
            PhysBCFunct<GpuBndryFuncFab<AmrCoreFill> > cphysbc(geom[lev-1],bcs,gpu_bndry_func);
            PhysBCFunct<GpuBndryFuncFab<AmrCoreFill> > fphysbc(geom[lev],bcs,gpu_bndry_func);
            amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,0, icomp, ncomp, geom[lev-1], geom[lev],cphysbc, 0, fphysbc, 0, refRatio(lev-1),mapper, bcs, 0);
        }
        else
        {   CpuBndryFuncFab bndry_func(nullptr);  // Without EXT_DIR, we can pass a nullptr.
            PhysBCFunct<CpuBndryFuncFab> cphysbc(geom[lev-1],bcs,bndry_func);
            PhysBCFunct<CpuBndryFuncFab> fphysbc(geom[lev],bcs,bndry_func);
            amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,0, icomp, ncomp, geom[lev-1], geom[lev],cphysbc, 0, fphysbc, 0, refRatio(lev-1),mapper, bcs, 0);
        }
    }
}

//####################################################################################################################

void
AmrCoreGP::FillCoarsePatch (std::string ref, int lev, Real time, MultiFab& mf, int icomp, int ncomp)
{   
    BL_ASSERT(lev > 0);
    Vector<MultiFab*> cmf;
    Vector<Real> ctime;
    GetData(ref,lev-1, time, cmf, ctime);
    if(ref=="M")
    {   
        bcs.resize(numcom);
        bcs = bcs_mu;
    }
    else if(ref=="C")
    {   
        bcs.resize(numcom);
        bcs = bcs_comp;
    }
    else
    {   
        bcs.resize(nump);
        bcs = bcs_phi;
    }
    Interpolater* mapper = &lincc_interp;
    if (cmf.size() != 1) 
    {amrex::Abort("FillCoarsePatch: how did this happen?");}
    if(Gpu::inLaunchRegion()) 
    {   GpuBndryFuncFab<AmrCoreFill> gpu_bndry_func(AmrCoreFill{});
        PhysBCFunct<GpuBndryFuncFab<AmrCoreFill> > cphysbc(geom[lev-1],bcs,gpu_bndry_func);
        PhysBCFunct<GpuBndryFuncFab<AmrCoreFill> > fphysbc(geom[lev],bcs,gpu_bndry_func);
        amrex::InterpFromCoarseLevel(mf, time, *cmf[0], 0, icomp, ncomp, geom[lev-1], geom[lev],cphysbc, 0, fphysbc, 0, refRatio(lev-1),mapper, bcs, 0);
    }
    else
    {   CpuBndryFuncFab bndry_func(nullptr);   // Without EXT_DIR, we can pass a nullptr.
        PhysBCFunct<CpuBndryFuncFab> cphysbc(geom[lev-1],bcs,bndry_func);
        PhysBCFunct<CpuBndryFuncFab> fphysbc(geom[lev],bcs,bndry_func);
        amrex::InterpFromCoarseLevel(mf, time, *cmf[0], 0, icomp, ncomp, geom[lev-1], geom[lev],cphysbc, 0, fphysbc, 0, refRatio(lev-1),mapper, bcs, 0);
    }
}


//####################################################################################################################

void
AmrCoreGP::GetData (std::string ref,int lev, Real time, Vector<MultiFab*>& data, Vector<Real>& datatime)
{   data.clear();
    datatime.clear();
    if(ref=="P") {data.push_back(&phi_new[lev]);};
    if(ref=="M") {data.push_back(&mu_new[lev]);};
    if(ref=="C") {data.push_back(&comp_new[lev]);};
    if(ref=="T1") {data.push_back(&term1[lev]);};
    if(ref=="T2") {data.push_back(&term2[lev]);};
    if(ref=="T3") {data.push_back(&term3[lev]);};
    if(ref=="T4") {data.push_back(&term4[lev]);};
    if(ref=="L") {data.push_back(&lambad[lev]);};
    if(ref=="PS") {data.push_back(&psi[lev]);};
    datatime.push_back(dtlev[0]);
}

//####################################################################################################################

void
AmrCoreGP::AverageDown ()
{   for (int lev = finest_level-1; lev >= 0; --lev)
    {   amrex::average_down(phi_new[lev+1], phi_new[lev],geom[lev+1], geom[lev],0, phi_new[lev].nComp(), refRatio(lev));
        amrex::average_down(mu_new[lev+1], mu_new[lev],geom[lev+1], geom[lev],0, mu_new[lev].nComp(), refRatio(lev));
        amrex::average_down(comp_new[lev+1], comp_new[lev],geom[lev+1], geom[lev],0, mu_new[lev].nComp(), refRatio(lev));
  }
}

//####################################################################################################################
std::string
AmrCoreGP::PlotFileName (int lev) const
{
    return amrex::Concatenate(plot_file, lev, 5);
}

//####################################################################################################################
Vector<const MultiFab*>
AmrCoreGP::PlotFileMF () const
{   Vector<const MultiFab*> r;
    for (int i = 0; i <= finest_level; ++i) 
    {r.push_back(&plot[i]);}
    return r;
}

//####################################################################################################################
Vector<std::string>
AmrCoreGP::PlotFileVarNames () const
{   
    return {"phi","liquid","mu","comp"};
}

//####################################################################################################################
void
AmrCoreGP::WritePlotFile () 
{
    const std::string& plotfilename = PlotFileName(istep[0]);
    for(int m=0; m < numcom-1; m++)
    {
		phase.push_back("mu_"+comp[m]);
	}
	
	for(int m=0; m < numcom-1; m++)
    {
		phase.push_back("comp_"+comp[m]);
	}
    const auto& mf = PlotFileMF();
    //const auto& varnames = PlotFileVarNames();
    amrex::WriteMultiLevelPlotfile(plotfilename, finest_level+1, mf, phase,Geom(), dtlev[0], istep, refRatio());
    //amrex::WriteMultiLevelPlotfile(plotfilename, finest_level+1, mf, varnames,Geom(), dtlev[0], istep, refRatio());
}

//####################################################################################################################

void
AmrCoreGP::InitData()
{
    const Real time = 0.0;
    if(funcf == 4)
    {
        function_F_04_function_A();
        function_F_04_function_B();
        function_F_04_function_C();
        function_F_04_dc_dmu();
    }
    Calculate_Tau();
    InitFromScratch(time);
    AverageDown();
    for (int lev = 0; lev <= finest_level; lev++)
    {
        //plot[lev].define(grids[lev],dmap[lev],nump+2*(numcom-1),0);
        plot[lev].define(grids[lev],dmap[lev],nump+2*numcom+2,0);
        MultiFab::Copy(plot[lev], phi_new[lev], 0, 0, nump, 0);
	    MultiFab::Copy(plot[lev], mu_new[lev],  0, nump, numcom-1, 0);
        MultiFab::Copy(plot[lev], comp_new[lev], 0, nump+numcom-1, numcom-1, 0);
    }; 
    if (plot_int > 0) 
    {
        WritePlotFile();
    }
}

//####################################################################################################################
void 
AmrCoreGP::timeStepWithSubcycling(amrex::Real time, amrex::Real dt, int lev)
{   if (regrid_int > 0)  
    { 
        static Vector<int> last_regrid_step(max_level+1, 0);
        if (lev < max_level && istep[lev] > last_regrid_step[lev])
        { 
            if (istep[lev] % regrid_int == 0)
            { 
                int old_finest = finest_level;
                regrid(lev, time);
                for (int k = lev; k <= finest_level; ++k) 
                {last_regrid_step[k] = istep[k];}
            }
        }
    }
    AdvancePhiAtLevel(time, dt, lev);
    istep[lev] = istep[lev]+1;
    if (lev<finest_level)
    {
        amrex::Real time_step = dt/nsubsteps[lev+1];
        for(int i=1;i<=nsubsteps[lev+1];i++)
        {
            timeStepWithSubcycling(time,time_step,lev+1);
            time = time+time_step;
        }
        AverageDown ();
    }
}

//####################################################################################################################

void 
AmrCoreGP::timeStepNoSubcycling(amrex::Real time, amrex::Real dt)
{
    if (max_level > 0 && regrid_int > 0)  
    {   
        if ((istep[0] % regrid_int) == 0)
        {   
            regrid(0, time);
        }
    }
    AdvancePhiAllLevels(time, dt);
    AverageDown ();
    for (int lev = 0; lev <= finest_level; lev++)
    {   
        istep[lev] = istep[lev] + 1;
    };    
}

//####################################################################################################################

void
AmrCoreGP::Evolve()
{   
    dtlev.resize(max_level+1,dt);
    int last_plot_file_step = 0;  
    for (int step = istep[0]; step < max_step; ++step)
    {   
        if (do_subcycle==1)
        {
            int lev = 0;
            timeStepWithSubcycling(time,dt,lev);
        }
        else
        {
            timeStepNoSubcycling(time, dt);
        }
        time = time + dt;
        dtlev[0]=dtlev[0]+dt;
    
        if (plot_int > 0 && istep[0] % plot_int == 0) 
        {   
            last_plot_file_step = istep[0] + 1;
            for (int lev = 0; lev <= finest_level; lev++)
            {   
                plot[lev].define(grids[lev],dmap[lev],nump+2*(numcom-1),0);
                MultiFab::Copy(plot[lev], phi_new[lev], 0, 0, nump, 0);
	            MultiFab::Copy(plot[lev], mu_new[lev],  0, nump, numcom-1, 0);
                MultiFab::Copy(plot[lev], comp_new[lev], 0, nump+numcom-1, numcom-1, 0);
            };
            WritePlotFile();
        }
    }
    if (plot_int > 0 && istep[0] > last_plot_file_step) 
    {
        WritePlotFile();
    }
}

//####################################################################################################################

