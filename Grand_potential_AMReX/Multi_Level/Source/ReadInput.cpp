#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <iostream>
#include <fstream>
#include <string>
#include "AmrCoreGP.H"


using namespace std;
using namespace amrex; 

void
AmrCoreGP::ReadInput()
{
     ParmParse pp;
        pp.get("regrid_int",regrid_int);
        pp.get("do_subcycle",do_subcycle);
        pp.get("max_step",max_step);
        pp.get("DELTA_t",dt);
        pp.get("DIM",dim);
        pp.get("NUMPHASES",nump);
        pp.get("NUMCOMPONENTS",numcom);
        pp.getarr("COMPONENTS",comp);
        pp.getarr("PHASES",phase);
        pp.getarr("GAMMA",gammaa);

        int ss = pp.countname("DIFFUSIVITY");
    	for(int i=0;i<ss;i++)
        {   
            Vector<Real> n1;
    	    pp.getktharr ("DIFFUSIVITY", i, n1, 0); 
    	    diffu.push_back(n1);
    	}
        ss=0;

        pp.query("R",R);
        pp.get("V",Vm);
        pp.query("ELASTICITY",ELASTICITY);
        pp.query("rho",ro);
        pp.query("damping_factor",dampfac);
        pp.query("max_iterations",MAX_ITERATIONS);

        ss = pp.countname("EIGEN_STRAIN");
    	for(int i=0;i<ss;i++)
        {
            Vector<Real> n1;
    	    pp.queryktharr ("EIGEN_STRAIN", i, n1, 0); 
    	    egstr.push_back(n1);
    	}
        ss=0;

        ss = pp.countname("VOIGT_ISOTROPIC");
    	for(int i=0;i<ss;i++)
        {
            Vector<Real> n1;
    	    pp.queryktharr ("VOIGT_ISOTROPIC", i, n1, 0); 
    	    voigiso.push_back(n1);
    	}
        ss=0;

        ss = pp.countname("BOUNDARY");
    	for(int i=0;i<ss;i++)
        {
            Vector<std::string> n2;
    	    pp.getktharr ("BOUNDARY", i, n2, 0); 
    	    bound.push_back(n2);
    	}
        ss=0;

        ss = pp.countname("BOUNDARY_VALUE");
    	for(int i=0;i<ss;i++)
        {
            Vector<std::string> n2;
    	    pp.getktharr ("BOUNDARY_VALUE", i, n2, 0); 
    	    boundval.push_back(n2);
    	}
        ss=0;

        pp.get("ISOTHERMAL",isothermal);
        pp.query("BINARY",binary);
        pp.query("TERNARY",ternary);
        pp.query("DILUTE",dilute);
        pp.get("T",T);
        pp.query("WRITEFORMAT",writeformat);
        pp.get("WRITEHDF5",writehdf5);
        pp.get("TRACK_PROGRESS",trackprog);
        pp.get("epsilon",eps);
        pp.get("tau",tau);
        pp.get("Tau",Tau);
        pp.get("Function_anisotropy",funcANI);
        pp.get("Anisotropy_type",ANItype);
        pp.getarr("dab",dab);

        ss = pp.countname("Rotation_matrix");
    	for(int i=0;i<ss;i++)
        {
            Vector<Real> n1;
    	    pp.getktharr ("Rotation_matrix", i, n1, 0); 
    	    rotmat.push_back(n1);
    	}
        ss=0;

        pp.get("Function_W",funcW);
        pp.query("Gamma_abc",gamma_abc);
         pp.query("Shift",shiftdom);
        pp.query("Shiftj",shiftj);
        pp.get("Equilibrium_temperature",Teq);
        pp.get("Filling_temperature",Tfill);
        pp.queryarr("Tempgrady",tempgrady);
        pp.get("Function_F",funcf);
        pp.get("num_thermo_phases",ntp);
        pp.get("tdbfname",tdbname);
        pp.getarr("tdb_phases",tdb_phase);
        pp.getarr("phase_map",phasemap);

        ss = pp.countname("FILLCUBE");
    	for(int i=0;i<ss;i++)
        {
            Vector<Real> n1;
    	    pp.getktharr ("FILLCUBE", i, n1, 0); 
    	    cube.push_back(n1);
    	}
        ss=0;

        ss = pp.countname("FILLCYLINDER");
    	for(int i=0;i<ss;i++)
        {
            Vector<Real> n1;
    	    pp.getktharr ("FILLCYLINDER", i, n1, 0); 
    	    cylinder.push_back(n1);
    	}
        ss=0;

        ss = pp.countname("FILLSPHERE");
    	for(int i=0;i<ss;i++)
        {
            Vector<Real> n1;
    	    pp.getktharr ("FILLSPHERE", i, n1, 0); 
    	    sphere.push_back(n1);
    	}
        ss=0;

        ss = pp.countname("FILLELLIPSE");
    	for(int i=0;i<ss;i++)
        {
            Vector<Real> n1;
    	    pp.getktharr ("FILLELLIPSE", i, n1, 0); 
    	    ellipse.push_back(n1);
    	}
        ss=0;

        ss = pp.countname("FILLCYLINDERRANDOM");
    	for(int i=0;i<ss;i++)
        {
            Vector<Real> n1;
    	    pp.getktharr ("FILLCYLINDERRANDOM", i, n1, 0); 
    	    cylrand.push_back(n1);
    	}
        ss=0;

        ss = pp.countname("FILLSPHERERANDOM");
    	for(int i=0;i<ss;i++)
        {
            Vector<Real> n1;
    	    pp.getktharr ("FILLSPHERERANDOM", i, n1, 0); 
    	    sphrand.push_back(n1);
    	}
        ss=0;

        pp.queryarr("FILLCUBERANDOM",cuberand);
        pp.queryarr("FILLCUBEPATTERN",cubepat);
        
        ss = pp.countname("ceq");
    	for(int i=0;i<ss;i++)
        {
            Vector<Real> n1;
    	    pp.queryktharr ("ceq", i, n1, 0); 
    	    ceq.push_back(n1);
    	}
        ss=0;

        ss = pp.countname("cfill");
    	for(int i=0;i<ss;i++)
        {
            Vector<Real> n1;
    	    pp.queryktharr ("cfill", i, n1, 0); 
    	    cfill.push_back(n1);
    	}
        ss=0;
    {
        ParmParse pp("amr"); // Traditionally, these have prefix, amr.
        pp.query("regrid_int", regrid_int);
        pp.query("plot_file", plot_file);
        pp.query("plot_int", plot_int);  
    }

    std::stringstream ss1;
    unsigned int count{0};
    gam = Vector<Vector<Real>>(nump,Vector<Real>(nump,0.0));
    int m=0;
    for (int a=0; a<nump; a++)
    {
        for(int b=0; b<nump; b++)
        {
            if(a!=b && a<b)
            {
                gam[a][b] = gammaa[m];
                m++;
            }
        }
    }

    diff = Vector<Vector<Vector<Real>>>(nump,Vector<Vector<Real>>(numcom-1,Vector<Real>(numcom-1,0.0)));
    int cnt=0;
    for (int a=0; a<nump; a++)
    {
        for(int k=0; k<numcom-1; k++)
        {
            for(int l=0; l<numcom-1; l++)
            {
                if(k==l)
                {
                    diff[a][k][l] = diffu[a][2+cnt];
                    cnt++;
                } 
            }
        }
        cnt=0;
    }
    
    count=0;
    dabb = Vector<Vector<Real>>(nump,Vector<Real>(nump,0.0));
    m=0;
    for (int a=0; a<nump; a++)
    {
        for(int b=0; b<nump; b++)
        {
            if(a!=b && a<b)
            {
                dabb[a][b] = dab[m];
                m++;
            }
        }
    }

    if(nump>3)
    {
        // for(char c:chars){
        //         gamma_abc.erase(remove(gamma_abc.begin(),gamma_abc.end(),c),gamma_abc.end());
        //     }
    
        ss1.str("");
        ss1.clear();
        ss1.str(gamma_abc);
        while(ss1.good())
        {
            std::string substr;
            getline(ss1,substr,',');
            val.push_back(substr);
        }

    ss1.str("");
    ss1.clear();
    for (int i = 0; i < val.size(); i++)
    {
        gammaa_abc.push_back(std::stod(val[i]));
    }
    val.clear();
    }

    gam_abc = Vector<Vector<Vector<Real>>>(nump,Vector<Vector<Real>>(nump,Vector<Real>(nump,0)));
    m=0;
    for(int a=0; a<nump; a++)
    {
        for(int b=0; b<nump; b++)
        {
            for(int d=0; d<nump; d++)
            {
                if(a<b && b<d)
                {
                    gam_abc[a][b][d] = gammaa_abc[m];
                    m++;
                }
            }
        }
    }  


    c_eq = Vector<Vector<Vector<Real>>>(nump,Vector<Vector<Real>>(nump,Vector<Real>(numcom,0)));
    for(auto k=2; k<ceq[0].size(); k++)
    {
    	for(auto i=0; i<(int)ceq.size(); i++)
        {
            auto v = (long unsigned int)ceq[i][0];
            auto w = (long unsigned int)ceq[i][1]; 
    		c_eq[v][w][count] = ceq[i][k];
        }
    count++;
    }	  
    count=0;
    
    c_fill = Vector <Vector<Vector<Real>>>(nump,Vector<Vector<Real>>(nump,Vector<Real>(numcom,0)));
    for(auto k=2; k<cfill[0].size(); k++)
    {
    	for(auto i=0; i<(int)cfill.size(); i++)
        {
            auto v = (long unsigned int)cfill[i][0];
            auto w = (long unsigned int)cfill[i][1]; 
    		c_fill[v][w][count] = cfill[i][k];
        }
    count++;
    }	  
    count=0;

}



































                        