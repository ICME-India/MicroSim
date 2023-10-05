    //! Reading thermodynamic database
    
    int buffer = 2048;
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
    if ((n_phase != 2)&&(n_phase != 3)&&(n_phase != 4))
    {
        cout << "This solver is only for two, three and four phase systems" << endl;
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

    int np;
    double T1[buffer], A_arr[buffer][n_phase], H11_arr[buffer][n_phase], H12_arr[buffer][n_phase], H22_arr[buffer][n_phase];

    for (i_phase = 0; i_phase < n_phase; i_phase++){
    const fileName pathToFile1 = "constant/HSN_"+phase[i_phase]+".csv";
    ifstream inpfs(pathToFile1);
    if (!inpfs.is_open())
    {
        Info << "No hessian data found" << endl;
    }
    //string line;
    getline(inpfs, line);
    //double T1[buffer], ASol[buffer], HSol11[buffer], HSol12[buffer], HSol22[buffer];
    np = 0;

    while (getline(inpfs, line))
    {
        string line_value;
        std::istringstream ss(line);
        getline(ss, line_value, ',');
        T1[np] = stod(line_value);
        if (n_compo == 2)
        {
        getline(ss, line_value);
        A_arr[np][i_phase] = 0.5*stod(line_value);
        }
        else if (n_compo == 3)
        {
        getline(ss, line_value, ',');
        H11_arr[np][i_phase] = stod(line_value);
        getline(ss, line_value, ',');
        H22_arr[np][i_phase] = stod(line_value);
        getline(ss, line_value);
        H12_arr[np][i_phase] = stod(line_value);
        }
    
    //Info << T1[np] << " " << H11_arr[np][i_phase] << endl;
    //cout << T1[np] << " " << ASol[np] << endl;
            
    ++np;
    }
    
    inpfs.close();
    }
    
    /*
    double ALiq[buffer], HLiq11[buffer], HLiq12[buffer], HLiq22[buffer];
    
    const fileName pathToFile2 = "constant/HSN_"+phase[n_phase-1]+".csv";
    ifstream inpfl(pathToFile2);
    if (!inpfl.is_open())
    {
        Info << "No liquid data found" << endl;
    }
    getline(inpfl, line);
    np = 0;

    while (getline(inpfl, line))
    {
        string line_value;
        std::istringstream ss(line);
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

    inpfl.close();*/
    
    double c_arr[buffer][n_phase], c1_arr[buffer][n_phase], c2_arr[buffer][n_phase];
    //double cLiq[buffer], cLiq1[buffer], cLiq2[buffer];
    
    for (i_phase = 0; i_phase < (n_phase-1); i_phase++){
    const fileName pathToFile3 = "constant/Composition_"+phase[i_phase]+".csv";
    ifstream inpfc(pathToFile3);
    if (!inpfc.is_open())
    {
        Info << "No composition data found" << endl;
    }
    getline(inpfc, line);
    np = 0;
    
    while (getline(inpfc, line))
    {
        string line_value;
        std::istringstream ss(line);
        getline(ss, line_value, ',');
        T1[np] = stod(line_value);
        getline(ss, line_value, ',');
        if (n_compo == 2)
        {
        c_arr[np][i_phase] = stod(line_value);
        getline(ss, line_value);
        c_arr[np][int(phases)-1] = stod(line_value);
        }
        else if (n_compo == 3)
        {
        c1_arr[np][i_phase] = stod(line_value);
        getline(ss, line_value, ',');
        c2_arr[np][i_phase] = stod(line_value);
        getline(ss, line_value, ',');
        c1_arr[np][int(phases)-1] = stod(line_value);
        getline(ss, line_value);
        c2_arr[np][int(phases)-1] = stod(line_value);
        }
    
    //Info << T1[np] << " " << c1_arr[np][i_phase] << " " << c2_arr[np][i_phase] << endl;
    //cout << T1[np] << " " << cSol[np] << " " << cLiq[np] << endl;
            
    ++np;
    }
    
    inpfc.close();
    }
    
    double T2[buffer], Lf_arr[buffer], Cp_arr[buffer];
    int np1 = 0;
    
    const fileName pathToFile4 = "constant/Thermal.csv";
    ifstream inpft(pathToFile4);
    if (!inpft.is_open())
    {
        Info << "No thermal data found" << endl;
    }
    getline(inpft, line);
    getline(inpft, line);
    
    while (getline(inpft, line))
    {
        string line_value;
        std::istringstream ss(line);
        getline(ss, line_value, ',');
        T2[np1] = stod(line_value);
        getline(ss, line_value, ',');
        Cp_arr[np1] = stod(line_value);
        getline(ss, line_value);
        Lf_arr[np1] = stod(line_value);
    
    //Info << T2[np1] << " " << Lf_arr[np1] << " " << Cp_arr[np1] << endl;
    //cout << T2[np1] << " " << Lf_arr[np1] << " " << Cp_arr[np1] << endl;
            
    ++np1;
    }
    
    inpft.close();
    
