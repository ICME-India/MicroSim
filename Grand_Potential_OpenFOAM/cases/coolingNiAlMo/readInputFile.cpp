//This code is used for reading MicroSim input files and converting those to OpenFOAM format
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

int main(int argc, char *argv[])
{
    
    int DIMENSION, MESH_X, MESH_Y, MESH_Z, NUMPHASES, NUMCOMPONENTS, i_phase, Function_anisotropy;
    int ELASTICITY = 0, GRAIN_GROWTH = 0, swch; 
    double DELTA_X, DELTA_t, NTIMESTEPS, SAVET, STARTTIME, DIFFUSIVITY00, DIFFUSIVITY01;
    double GAMMA1, GAMMA2, GAMMA3, GAMMA4, V, DIFFUSIVITY[2], EIGEN_STRAIN[4][6];
    double VOIGT[4][6], T, epsilon;
    double DIFFUSIVITYl0, DIFFUSIVITYl1, dab, Amp_Noise_Phase, Equilibrium_temperature, Filling_temperature;
    double theta_x, theta_y, theta_z, center_x[3], center_y[3], center_z[3], seed_radius[3], volume_fraction[3], shield_dist[3], spread[3];
    double DIFFUSIVITY10 = 0, DIFFUSIVITY11 = 0, DIFFUSIVITY20 = 0, DIFFUSIVITY21 = 0;
    
    //ifstream inpf("Input_tdb_new.in");
    ifstream inpf(argv[1]);
    if (!inpf.is_open())
    {
        cout << "No input data found" << endl;
    }
    
    //ifstream fill("Filling.in");
    ifstream fill(argv[2]);
    if (!fill.is_open())
    {
        cout << "No filling data found" << endl;
    }
    
    ofstream outpf("generatedInput");
    if (!outpf.is_open())
    {
	cout << "Input was not generated" << endl;
	}
	
    ofstream outpft("constant/temperature");
    if (!outpft.is_open())
    {
	cout << "temperature file was not generated" << endl;
	}

    ofstream outpfp("constant/phases");
    if (!outpfp.is_open())
    {
	cout << "phases file was not generated" << endl;
	}
    
    string line;
    while (getline(inpf, line))
    {
    if (line[0] != '#')
    {

    string line_value1, line_value2, line_value;
    istringstream ss1(line);
    getline(ss1, line_value1, ' ');
    getline(ss1, line_value, ' ');
    getline(ss1, line_value2, ';');
    
    //reading input variables
    if (line_value1 == "DIMENSION")
    {
    DIMENSION = stoi(line_value2);
    }
    
    else if (line_value1 == "MESH_X")
    {
    MESH_X = stoi(line_value2);
    }
    
    else if (line_value1 == "MESH_Y")
    {
    MESH_Y = stoi(line_value2);
    }
    
    else if (line_value1 == "MESH_Z")
    {
    MESH_Z = stoi(line_value2);
    }
    
    else if (line_value1 == "DELTA_X")
    {
    DELTA_X = stod(line_value2);
    }
    
    else if (line_value1 == "DELTA_t")
    {
    DELTA_t = stod(line_value2);
    }
    
    else if (line_value1 == "NUMPHASES")
    {
    NUMPHASES = stoi(line_value2);
    }
    
    else if (line_value1 == "NUMCOMPONENTS")
    {
    NUMCOMPONENTS = stoi(line_value2);
    outpfp << NUMCOMPONENTS << endl;
    outpfp << NUMPHASES << endl;
    }
    
    else if (line_value1 == "NTIMESTEPS")
    {
    NTIMESTEPS = stod(line_value2);
    }
    
    else if (line_value1 == "SAVET")
    {
    SAVET = stod(line_value2);
    }
    
    else if (line_value1 == "STARTTIME")
    {
    STARTTIME = stod(line_value2);
    }
    
    else if (line_value1 == "phase_map")
    {
    istringstream ss2(line_value2);
    getline(ss2, line_value2, '{');
    
    int i_phase = 0;
    while (i_phase < (NUMPHASES-1)){
    getline(ss2, line_value2, ',');
    if (line_value2[0] == ' '){
    for (int i = 1; i < line_value2.size(); i++){
    outpfp << line_value2[i];
    }
    }
    else{
    outpfp << line_value2;
    }
    outpfp << endl;
    i_phase = i_phase + 1;
    }
    
    getline(ss2, line_value2, '}');
    if (line_value2[0] == ' '){
    for (int i = 1; i < line_value2.size(); i++){
    outpfp << line_value2[i];
    }
    }
    else{
    outpfp << line_value2;
    }
    outpfp << endl;
    }
    
    else if (line_value1 == "GAMMA")
    {
    istringstream ss2(line_value2);
    getline(ss2, line_value2, '{');
    getline(ss2, line_value2, ',');
    
    GAMMA1 = stod(line_value2);

    getline(ss2, line_value2, ',');
    
    GAMMA2 = stod(line_value2);

    getline(ss2, line_value2, ',');
    
    GAMMA3 = stod(line_value2);

    getline(ss2, line_value2, '}');
    
    GAMMA4 = stod(line_value2);
    }
    
    else if (line_value1 == "DIFFUSIVITY")
    {
    istringstream ss2(line_value2);
    getline(ss2, line_value2, '{');
    getline(ss2, line_value2, ',');
    
    getline(ss2, line_value2, ',');
    i_phase = stod(line_value2);
    
    getline(ss2, line_value2, ',');
    DIFFUSIVITY[0] = stod(line_value2);
    
    getline(ss2, line_value2, '}');
    DIFFUSIVITY[1] = stod(line_value2);
    
    if (i_phase == 0)
    {
    DIFFUSIVITY00 = DIFFUSIVITY[0];
    DIFFUSIVITY01 = DIFFUSIVITY[1];
    }
    
    else if (i_phase == (NUMPHASES-1))
    {
    DIFFUSIVITYl0 = DIFFUSIVITY[0];
    DIFFUSIVITYl1 = DIFFUSIVITY[1];
    }

    if (NUMPHASES > 2)
    {
    if (i_phase == 1)
    {
    DIFFUSIVITY10 = DIFFUSIVITY[0];
    DIFFUSIVITY11 = DIFFUSIVITY[1];
    }
    }

    if (NUMPHASES > 3)
    {
    if (i_phase == 2)
    {
    DIFFUSIVITY20 = DIFFUSIVITY[0];
    DIFFUSIVITY21 = DIFFUSIVITY[1];
    }
    }

    }
    
    else if (line_value1 == "V")
    {
    V = stod(line_value2);
    }

    else if (line_value1 == "ELASTICITY")
    {
    ELASTICITY = stoi(line_value2);
    }

    else if (line_value1 == "GRAIN_GROWTH")
    {
    GRAIN_GROWTH = stoi(line_value2);
    }
    
    else if (line_value1 == "EIGEN_STRAIN")
    {
    istringstream ss2(line_value2);
    getline(ss2, line_value2, '{');
    getline(ss2, line_value2, ',');
    i_phase = stod(line_value2);
    
    getline(ss2, line_value2, ',');
    EIGEN_STRAIN[i_phase][0] = stod(line_value2);
    
    getline(ss2, line_value2, ',');
    EIGEN_STRAIN[i_phase][1] = stod(line_value2);
    
    getline(ss2, line_value2, ',');
    EIGEN_STRAIN[i_phase][2] = stod(line_value2);
    
    getline(ss2, line_value2, ',');
    EIGEN_STRAIN[i_phase][3] = stod(line_value2);
    
    getline(ss2, line_value2, ',');
    EIGEN_STRAIN[i_phase][4] = stod(line_value2);
    
    getline(ss2, line_value2, '}');
    EIGEN_STRAIN[i_phase][5] = stod(line_value2);
    }
    
    else if (line_value1 == "VOIGT_ISOTROPIC")
    {
    istringstream ss2(line_value2);
    getline(ss2, line_value2, '{');
    getline(ss2, line_value2, ',');
    i_phase = stod(line_value2);
    
    getline(ss2, line_value2, ',');
    VOIGT[i_phase][0] = stod(line_value2);
    
    getline(ss2, line_value2, ',');
    VOIGT[i_phase][1] = stod(line_value2);
    
    getline(ss2, line_value2, '}');
    VOIGT[i_phase][2] = stod(line_value2);
    }
    
    else if (line_value1 == "T")
    {
    T = stod(line_value2);
    }
    
    else if (line_value1 == "epsilon")
    {
    epsilon = stod(line_value2);
    }
    
    else if (line_value1 == "Function_anisotropy")
    {
    Function_anisotropy = stoi(line_value2);
    }
    
    else if (line_value1 == "dab")
    {
    istringstream ss2(line_value2);
    getline(ss2, line_value2, '{');
    getline(ss2, line_value2, '}');
    
    dab = stod(line_value2);
    
    if (Function_anisotropy == 0)
    {dab = 0;}
    }

    else if (line_value1 == "Rotation_matrix")
    {
    istringstream ss2(line_value2);
    getline(ss2, line_value2, '{');
    getline(ss2, line_value2, ',');
    if (stod(line_value2) == 0)
    {

    getline(ss2, line_value2, ',');
    if (stod(line_value2) == 1)
    {
    getline(ss2, line_value2, ',');
    theta_x = stod(line_value2);

    getline(ss2, line_value2, ',');
    theta_y = stod(line_value2);

    getline(ss2, line_value2, '}');
    theta_z = stod(line_value2);
    }
    }
    }
    
    else if (line_value1 == "Amp_Noise_Phase")
    {
    Amp_Noise_Phase = stod(line_value2);
    }
    
    else if (line_value1 == "Equilibrium_temperature")
    {
    Equilibrium_temperature = stod(line_value2);
    }
    
    else if (line_value1 == "Filling_temperature")
    {
    Filling_temperature = stod(line_value2);
    }
    
    
    }
    }

    inpf.close();

    while (getline(fill, line))
    {
    if (line[0] != '#')
    {

    string line_value1, line_value2, line_value;
    istringstream ss1(line);
    getline(ss1, line_value1, ' ');
    getline(ss1, line_value, ' ');
    getline(ss1, line_value2, ';');

    //reading fill in variables
    if (line_value1 == "FILLCYLINDER")
    {
    istringstream ss2(line_value2);
    getline(ss2, line_value2, '{');
    getline(ss2, line_value2, ',');
    i_phase = stod(line_value2);
    
    getline(ss2, line_value2, ',');
    center_x[i_phase] = stod(line_value2);

    getline(ss2, line_value2, ',');
    center_y[i_phase] = stod(line_value2);

    getline(ss2, line_value2, ',');
    center_z[i_phase] = stod(line_value2);

    getline(ss2, line_value2, ',');

    getline(ss2, line_value2, '}');
    seed_radius[i_phase] = stod(line_value2);
    volume_fraction[i_phase] = 0;
    shield_dist[i_phase] = 0;
    spread[i_phase] = 0;
    }

    else if (line_value1 == "FILLCYLINDERRANDOM")
    {
    istringstream ss2(line_value2);
    getline(ss2, line_value2, '{');
    getline(ss2, line_value2, ',');
    i_phase = stod(line_value2);
    
    getline(ss2, line_value2, ',');
    seed_radius[i_phase] = stod(line_value2);

    getline(ss2, line_value2, ',');
    volume_fraction[i_phase] = stod(line_value2);

    center_x[i_phase] = 0;
    center_y[i_phase] = 0;
    center_z[i_phase] = 0;

    getline(ss2, line_value2, ',');
    shield_dist[i_phase] = stod(line_value2);

    getline(ss2, line_value2, '}');
    spread[i_phase] = stod(line_value2);
    }

    else if (line_value1 == "FILLSPHERE")
    {
    istringstream ss2(line_value2);
    getline(ss2, line_value2, '{');
    getline(ss2, line_value2, ',');
    i_phase = stod(line_value2);

    getline(ss2, line_value2, ',');
    center_x[i_phase] = stod(line_value2);

    getline(ss2, line_value2, ',');
    center_y[i_phase] = stod(line_value2);

    getline(ss2, line_value2, ',');
    center_z[i_phase] = stod(line_value2);

    getline(ss2, line_value2, '}');
    seed_radius[i_phase] = stod(line_value2);
    volume_fraction[i_phase] = 0;
    shield_dist[i_phase] = 0;
    spread[i_phase] = 0;
    }

    else if (line_value1 == "FILLSPHERERANDOM")
    {
    istringstream ss2(line_value2);
    getline(ss2, line_value2, '{');
    getline(ss2, line_value2, ',');
    i_phase = stod(line_value2);

    getline(ss2, line_value2, ',');
    seed_radius[i_phase] = stod(line_value2);

    getline(ss2, line_value2, ',');
    volume_fraction[i_phase] = stod(line_value2);

    center_x[i_phase] = 0;
    center_y[i_phase] = 0;
    center_z[i_phase] = 0;

    getline(ss2, line_value2, ',');
    shield_dist[i_phase] = stod(line_value2);

    getline(ss2, line_value2, '}');
    spread[i_phase] = stod(line_value2);
    }

    }
    }

    fill.close();

    //cout << MESH_X << endl;

	NTIMESTEPS = NTIMESTEPS*DELTA_t;
    SAVET = SAVET*DELTA_t;
    theta_x = theta_x*3.14159/180;
    theta_y = theta_y*3.14159/180;
    theta_z = theta_z*3.14159/180;

    if (ELASTICITY == 0)
    {
        if (GRAIN_GROWTH == 0)
        {
            swch = 0;
        }
        else if (GRAIN_GROWTH == 1)
        {
            swch = 1;
        }
    }
    else if (ELASTICITY == 1)
    {
        swch = 2;
    }

	//writing to include in openfoam dictionaries
    outpf.precision(15);
    //outpfp << "NUMPHASES " << NUMPHASES << ";" << endl;
    //outpfp << "NUMCOMPONENTS " << NUMCOMPONENTS << ";" << endl;
    
    outpf << "DIMENSION " << DIMENSION << ";" << endl;
    outpf << "NUMPHASES " << NUMPHASES << ";" << endl;
    outpf << "NUMCOMPONENTS " << NUMCOMPONENTS << ";" << endl;
    outpf << "MESH_X " << MESH_X << ";" << endl;
    outpf << "MESH_Y " << MESH_Y << ";" << endl;
    outpf << "MESH_Z " << MESH_Z << ";" << endl;
    outpf << "DELTA_X " << DELTA_X << ";" << endl;
    
    outpf << "DELTA_t " << DELTA_t << ";" << endl;    
    outpf << "NTIMESTEPS " << NTIMESTEPS << ";" << endl;
    outpf << "SAVET " << SAVET << ";" << endl;
    outpf << "STARTTIME " << STARTTIME << ";" << endl;
    
    outpf << "GAMMA1 " << GAMMA1 << ";" << endl;
    outpf << "GAMMA2 " << GAMMA2 << ";" << endl;
    outpf << "GAMMA3 " << GAMMA3 << ";" << endl;
    outpf << "GAMMA4 " << GAMMA4 << ";" << endl;
    outpf << "epsilon " << epsilon << ";" << endl;
    outpf << "dab " << dab << ";" << endl;
    outpf << "theta_x " << theta_x << ";" << endl;
    outpf << "theta_y " << theta_y << ";" << endl;
    outpf << "theta_z " << theta_z << ";" << endl;
    outpf << "Amp_Noise_Phase " << Amp_Noise_Phase << ";" << endl;
    outpf << "Vm " << V << ";" << endl;
    
    outpf << "DIFFUSIVITY00 " << DIFFUSIVITY00 << ";" << endl;
    outpf << "DIFFUSIVITY01 " << DIFFUSIVITY01 << ";" << endl;
    outpf << "DIFFUSIVITY10 " << DIFFUSIVITY10 << ";" << endl;
    outpf << "DIFFUSIVITY11 " << DIFFUSIVITY11 << ";" << endl;
    outpf << "DIFFUSIVITY20 " << DIFFUSIVITY20 << ";" << endl;
    outpf << "DIFFUSIVITY21 " << DIFFUSIVITY21 << ";" << endl;
    outpf << "DIFFUSIVITYl0 " << DIFFUSIVITYl0 << ";" << endl;
    outpf << "DIFFUSIVITYl1 " << DIFFUSIVITYl1 << ";" << endl;

    outpft << "initial " << Filling_temperature << ";" << endl;    
    outpft << "T0 " << Equilibrium_temperature << ";" << endl;

    outpf << "swch " << swch << ";" << endl;
    
    outpf << "EIGEN_STRAIN1 (" << EIGEN_STRAIN[0][0] << " " << EIGEN_STRAIN[0][5] << " " << EIGEN_STRAIN[0][4] <<
    " " << EIGEN_STRAIN[0][1] << " " << EIGEN_STRAIN[0][3] << " " << EIGEN_STRAIN[0][2] << ");" << endl;

    outpf << "EIGEN_STRAIN2 (" << EIGEN_STRAIN[1][0] << " " << EIGEN_STRAIN[1][5] << " " << EIGEN_STRAIN[1][4] <<
    " " << EIGEN_STRAIN[1][1] << " " << EIGEN_STRAIN[1][3] << " " << EIGEN_STRAIN[1][2] << ");" << endl;

    outpf << "EIGEN_STRAIN3 (" << EIGEN_STRAIN[2][0] << " " << EIGEN_STRAIN[2][5] << " " << EIGEN_STRAIN[2][4] <<
    " " << EIGEN_STRAIN[2][1] << " " << EIGEN_STRAIN[2][3] << " " << EIGEN_STRAIN[2][2] << ");" << endl;
    
    outpf << "C11_1 " << VOIGT[0][0] << ";" << endl;
    outpf << "C12_1 " << VOIGT[0][1] << ";" << endl;
    outpf << "C44_1 " << VOIGT[0][2] << ";" << endl;
    outpf << "C11_2 " << VOIGT[1][0] << ";" << endl;
    outpf << "C12_2 " << VOIGT[1][1] << ";" << endl;
    outpf << "C44_2 " << VOIGT[1][2] << ";" << endl;
    
    double mu1_elast = VOIGT[0][2];
    double lambda1 = VOIGT[0][1];
    double mu1_elast_ = (VOIGT[0][0] - VOIGT[0][1]) - 2*VOIGT[0][2];
    
    double mu2_elast = VOIGT[1][2];
    double lambda2 = VOIGT[1][1];
    double mu2_elast_ = (VOIGT[1][0] - VOIGT[1][1]) - 2*VOIGT[1][2];
    
    outpf << "mu1_elast " << mu1_elast << ";" << endl;
    outpf << "lambda1 " << lambda1 << ";" << endl;
    outpf << "mu1_elast_ " << mu1_elast_ << ";" << endl;
    
    outpf << "mu2_elast " << mu2_elast << ";" << endl;
    outpf << "lambda2 " << lambda2 << ";" << endl;
    outpf << "mu2_elast_ " << mu2_elast_ << ";" << endl;
    
    for (i_phase = 0; i_phase < (NUMPHASES-1); i_phase++) {
    outpf << "centerX[" << i_phase << "] " << DELTA_X*center_x[i_phase] << ";" << endl;
    outpf << "centerY[" << i_phase << "] " << DELTA_X*center_y[i_phase] << ";" << endl;
    outpf << "centerZ[" << i_phase << "] " << DELTA_X*center_z[i_phase] << ";" << endl;
    outpf << "seedRadius[" << i_phase << "] " << DELTA_X*seed_radius[i_phase] << ";" << endl;
    outpf << "volumeFraction[" << i_phase << "] " << volume_fraction[i_phase] << ";" << endl;
    outpf << "shieldDist[" << i_phase << "] " << DELTA_X*shield_dist[i_phase] << ";" << endl;
    outpf << "spread[" << i_phase << "] " << spread[i_phase] << ";" << endl;
    }

	outpf.close();
    
    return 0;
}
