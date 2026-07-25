#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
using namespace std;

int main()
{
    int buffer = 1024;
    ifstream inpfph("constant/phases");
    if (!inpfph.is_open())
    {
	cout << "No phase name found" << endl;
    }

    string line;
    getline(inpfph, line);
    int n_compo = stoi(line);
    if ((n_compo != 2)&&(n_compo != 3))
    {
        cout << "This solver is only for binary and ternary alloys" << endl;
	return 0;
    }

    getline(inpfph, line);
    int n_phase = stoi(line);
    if (n_phase != 2)
    {
        cout << "This solver is only for two phase systems" << endl;
	return 0;
    }

    //! Reading phase names
    string phase[n_phase];
    int i_phase = 0;

    while (getline(inpfph, line))
    {
	phase[i_phase] = line;
    ++i_phase;
    }

    inpfph.close();

    //const fileName pathToFile1 = "constant/HSN_FCC_A1.csv";
    /*ifstream inpfl("constant/HSN_"+phase[0]+".csv");
    if (!inpfs.is_open())
    {
        cout << "No solid data found" << endl;
    }
    string line;
    getline(inpfs, line);
    double T1[buffer], ASol[buffer];
    int np = 0;
    
    while (getline(inpfs, line))
    {
        string line_value;
        istringstream ss(line);
        getline(ss, line_value, ',');
        T1[np] = stod(line_value);
        getline(ss, line_value);
        ASol[np] = 0.5*stod(line_value);
    
    //Info << T1[np] << " " << ASol[np] << endl;
    cout << T1[np] << " " << ASol[np] << endl;
            
    ++np;
    }
    
    inpfs.close();
    */
    
    //const fileName pathToFile2 = "constant/HSN_LIQUID.csv";
    ifstream inpfl("constant/HSN_"+phase[1]+".csv");
    if (!inpfl.is_open())
    {
        cout << "No liquid data found" << endl;
    }
    //string line;
    getline(inpfl, line);
    double T1[buffer]; 
    double ALiq[buffer], HLiq11[buffer], HLiq12[buffer], HLiq22[buffer];
    int np = 0;

    while (getline(inpfl, line))
    {
        string line_value;
        istringstream ss(line);
        getline(ss, line_value, ',');
        T1[np] = stod(line_value);
        if (n_compo == 2)
        {
        getline(ss, line_value);
        ALiq[np] = 0.5*stod(line_value);
        }
        else if (n_compo == 3)
        {
        getline(ss, line_value, ',');
        HLiq11[np] = stod(line_value);
        getline(ss, line_value, ',');
        HLiq22[np] = stod(line_value);
        getline(ss, line_value);
        HLiq12[np] = stod(line_value);
        }
    
    //Info << T1[np] << " " << ALiq[np] << endl;
    //cout << T1[np] << " " << ALiq[np] << endl;
            
    ++np;
    }
    
    inpfl.close();
    
    //const fileName pathToFile3 = "constant/Composition_FCC_A1.csv";
    ifstream inpfc("constant/Composition_"+phase[0]+".csv");
    if (!inpfc.is_open())
    {
        cout << "No composition data found" << endl;
    }
    getline(inpfc, line);
    double cSol[buffer], cSol1[buffer], cSol2[buffer];
    double cLiq[buffer], cLiq1[buffer], cLiq2[buffer];
    np = 0;
    
    while (getline(inpfc, line))
    {
        string line_value;
        istringstream ss(line);
        getline(ss, line_value, ',');
        T1[np] = stod(line_value);
        getline(ss, line_value, ',');
        if (n_compo == 2)
        {
        cSol[np] = stod(line_value);
        getline(ss, line_value);
        cLiq[np] = stod(line_value);
        }
        else if (n_compo == 3)
        {
        cSol1[np] = stod(line_value);
        getline(ss, line_value, ',');
        cSol2[np] = stod(line_value);
        getline(ss, line_value, ',');
        cLiq1[np] = stod(line_value);
        getline(ss, line_value);
        cLiq2[np] = stod(line_value);
        }
    
    //Info << T1[np] << " " << cSol[np] << " " << cLiq[np] << endl;
    //cout << T1[np] << " " << cLiq[np] << endl;
            
    ++np;
    }

    inpfc.close();
    
    //cout << np << endl;
    
    ifstream inpft("constant/temperature");
    if (!inpft.is_open())
    {
        cout << "No temperature data found" << endl;
    }
    
    getline(inpft, line);
    string line_value;
    istringstream ss1(line);
    getline(ss1, line_value, ' ');
    getline(ss1, line_value);
    double Tinit = stod(line_value);
    
    cout << "Fill in temperature = " << Tinit << " K" << endl;
    
    getline(inpft, line);
    istringstream ss2(line);
    getline(ss2, line_value, ' ');
    getline(ss2, line_value);
    double T0 = stod(line_value);
    
    cout << "Equilibrium temperature = " << T0 << " K" << endl;

    inpft.close();
    
    ofstream outpf("0/muInitial");
    if (!outpf.is_open())
    {
		cout << "Initial mu was not printed" << endl;
	}
    
    if (n_compo == 2)
    {
    gsl_spline *spline1 = gsl_spline_alloc (gsl_interp_cspline, np);
    gsl_interp_accel *acc1 = gsl_interp_accel_alloc ();
    
    gsl_spline *spline2 = gsl_spline_alloc (gsl_interp_cspline, np);
    gsl_interp_accel *acc2 = gsl_interp_accel_alloc ();
   
    gsl_spline *spline3 = gsl_spline_alloc (gsl_interp_cspline, np);
    gsl_interp_accel *acc3 = gsl_interp_accel_alloc ();

    gsl_spline_init (spline1, T1, ALiq, np);
    double ALiqin = gsl_spline_eval (spline1, Tinit, acc1);
    //cout << ALiqin << endl;
    
    gsl_spline_init (spline2, T1, cLiq, np);
    double cLiq0 = gsl_spline_eval (spline2, T0, acc2);
    //cout << cLiq0 << endl;
    
    gsl_spline_init (spline3, T1, cSol, np);
    //double cSol0 = gsl_spline_eval (spline3, T0, acc3);
    //double m0 = gsl_spline_eval_deriv (spline3, T0, acc3);
    //cout << cSol0 << endl;
    //cout << "slope solidus, m = " << 1/m0 << endl;

    gsl_spline_free (spline1);
    gsl_spline_free (spline2);
    gsl_spline_free (spline3);
    gsl_interp_accel_free (acc1);
    gsl_interp_accel_free (acc2);
    gsl_interp_accel_free (acc3);
	
	double mu = 2.0*ALiqin*cLiq0;
    //ALiq at undercooling and cLiq at equilibrium

    outpf.precision(15);
    outpf << "muInit " << mu << ";" << endl;
    //cout << "muInit " << mu << ";" << endl;
    }
    
    else if (n_compo == 3)
    {
    gsl_spline *spline4 = gsl_spline_alloc (gsl_interp_cspline, np);
    gsl_interp_accel *acc4 = gsl_interp_accel_alloc ();
    gsl_spline_init (spline4, T1, HLiq11, np);
    
    gsl_spline *spline5 = gsl_spline_alloc (gsl_interp_cspline, np);
    gsl_interp_accel *acc5 = gsl_interp_accel_alloc ();
    gsl_spline_init (spline5, T1, HLiq12, np);
    
    gsl_spline *spline6 = gsl_spline_alloc (gsl_interp_cspline, np);
    gsl_interp_accel *acc6 = gsl_interp_accel_alloc ();
    gsl_spline_init (spline6, T1, HLiq22, np);
    
    gsl_spline *spline9 = gsl_spline_alloc (gsl_interp_cspline, np);
    gsl_interp_accel *acc9 = gsl_interp_accel_alloc ();
    gsl_spline_init (spline9, T1, cLiq1, np);

    gsl_spline *spline10 = gsl_spline_alloc (gsl_interp_cspline, np);
    gsl_interp_accel *acc10 = gsl_interp_accel_alloc ();
    gsl_spline_init (spline10, T1, cLiq2, np);

    //gsl_spline_init (spline1, T1, ALiq, np);
    double HLiq11in = gsl_spline_eval (spline4, Tinit, acc4);
    
    double HLiq12in = gsl_spline_eval (spline5, Tinit, acc5);
    
    double HLiq22in = gsl_spline_eval (spline6, Tinit, acc6);
    //cout << HLiq11in << " " << HLiq22in << " " << HLiq12in << endl;
    
    //gsl_spline_init (spline2, T1, cLiq, np);
    double cLiq10 = gsl_spline_eval (spline9, T0, acc9);
    
    double cLiq20 = gsl_spline_eval (spline10, T0, acc10);
    //cout << cLiq10 << " " << cLiq20 << endl;
    
    //gsl_spline_init (spline3, T1, cSol, np);
    //double cSol0 = gsl_spline_eval (spline3, T0, acc3);
    //double m0 = gsl_spline_eval_deriv (spline3, T0, acc3);
    //cout << cSol0 << endl;
    //cout << "slope solidus, m = " << 1/m0 << endl;

    gsl_spline_free (spline4);
    gsl_spline_free (spline5);
    gsl_spline_free (spline6);
    gsl_spline_free (spline9);
    gsl_spline_free (spline10);
    gsl_interp_accel_free (acc4);
    gsl_interp_accel_free (acc5);
    gsl_interp_accel_free (acc6);
    gsl_interp_accel_free (acc9);
    gsl_interp_accel_free (acc10);
	
	double mu1 = HLiq11in*cLiq10 + HLiq12in*cLiq20;
    double mu2 = HLiq12in*cLiq10 + HLiq22in*cLiq20;
    //ALiq at undercooling and cLiq at equilibrium

    outpf.precision(15);
    outpf << "mu1Init " << mu1 << ";" << endl;
    outpf << "mu2Init " << mu2 << ";" << endl;
    //cout << "muInit " << mu << ";" << endl;
    }

	outpf.close();
    
    return 0;
}
