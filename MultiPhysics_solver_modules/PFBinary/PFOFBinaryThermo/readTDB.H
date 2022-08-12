    //! Reading thermodynamic database
    
    int buffer = 1024;
    ifstream inpfph("constant/phases");
    if (!inpfph.is_open())
    {
	cout << "No phase name found" << endl;
    }

    string line;
    getline(inpfph, line);
    int n_compo = stoi(line);
    if (n_compo != 2)
    {
        cout << "This solver is only for binary alloys" << endl;
	return 0;
    }

    getline(inpfph, line);
    int n_phase = stoi(line);
    if (n_phase != 2)
    {
        cout << "This solver is only for two phase systems" << endl;
	return 0;
    }

    string phase[n_phase];
    int i_phase = 0;

    while (getline(inpfph, line))
    {
	phase[i_phase] = line;
    ++i_phase;
    }

    inpfph.close();

    const fileName pathToFile1 = "constant/HSN_"+phase[0]+".csv";
    ifstream inpfs(pathToFile1);
    if (!inpfs.is_open())
    {
        Info << "No solid data found" << endl;
    }
    //string line;
    getline(inpfs, line);
    double T1[buffer], ASol[buffer];
    int np = 0;

    while (getline(inpfs, line))
    {
        string line_value;
        std::istringstream ss(line);
        getline(ss, line_value, ',');
        T1[np] = stod(line_value);
        getline(ss, line_value);
        ASol[np] = 0.5*stod(line_value);
    
    //Info << T1[np] << " " << ASol[np] << endl;
    //cout << T1[np] << " " << ASol[np] << endl;
            
    ++np;
    }
    
    inpfs.close();
    
    const fileName pathToFile2 = "constant/HSN_"+phase[1]+".csv";
    ifstream inpfl(pathToFile2);
    if (!inpfl.is_open())
    {
        Info << "No liquid data found" << endl;
    }
    getline(inpfl, line);
    double ALiq[buffer];
    np = 0;

    while (getline(inpfl, line))
    {
        string line_value;
        std::istringstream ss(line);
        getline(ss, line_value, ',');
        T1[np] = stod(line_value);
        getline(ss, line_value);
        ALiq[np] = 0.5*stod(line_value);
    
    //Info << T1[np] << " " << ALiq[np] << endl;
    //cout << T1[np] << " " << ALiq[np] << endl;
            
    ++np;
    }

    inpfl.close();
    
    const fileName pathToFile3 = "constant/Composition_"+phase[0]+".csv";
    ifstream inpfc(pathToFile3);
    if (!inpfc.is_open())
    {
        Info << "No composition data found" << endl;
    }
    getline(inpfc, line);
    double cSol[buffer], cLiq[buffer];
    np = 0;
    
    while (getline(inpfc, line))
    {
        string line_value;
        std::istringstream ss(line);
        getline(ss, line_value, ',');
        T1[np] = stod(line_value);
        getline(ss, line_value, ',');
        cSol[np] = stod(line_value);
        getline(ss, line_value);
        cLiq[np] = stod(line_value);
    
    //Info << T1[np] << " " << cSol[np] << " " << cLiq[np] << endl;
    //cout << T1[np] << " " << cSol[np] << " " << cLiq[np] << endl;
            
    ++np;
    }
    
    inpfc.close();
    
    const fileName pathToFile4 = "constant/Thermal.csv";
    ifstream inpft(pathToFile4);
    if (!inpft.is_open())
    {
        Info << "No thermal data found" << endl;
    }
    getline(inpft, line);
    getline(inpft, line);
    double T2[buffer], Lf_arr[buffer], Cp_arr[buffer];
    int np1 = 0;
    
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
