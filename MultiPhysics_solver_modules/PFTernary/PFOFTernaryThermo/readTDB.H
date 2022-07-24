    //! Reading thermodynamic database
    
    int buffer = 1024;
    const fileName pathToFile1 = "constant/HSN_FCC_A1.csv";
    ifstream inpfs(pathToFile1);
    if (!inpfs.is_open())
    {
        Info << "No solid data found" << endl;
    }
    string line;
    getline(inpfs, line);
    double T1[buffer], HSol11[buffer], HSol12[buffer], HSol22[buffer];
    int np = 0;

    while (getline(inpfs, line))
    {
        string line_value;
        std::istringstream ss(line);
        getline(ss, line_value, ',');
        T1[np] = stod(line_value);
        getline(ss, line_value, ',');
        HSol11[np] = stod(line_value);
        getline(ss, line_value, ',');
        HSol22[np] = stod(line_value);
        getline(ss, line_value);
        HSol12[np] = stod(line_value);
    
    //Info << T1[np] << " " << HSol11[np] << " " << HSol22[np] << " " << HSol12[np] << endl;
    //cout << T1[np] << " " << ASol[np] << endl;
            
    ++np;
    }
    
    inpfs.close();
    
    const fileName pathToFile2 = "constant/HSN_LIQUID.csv";
    ifstream inpfl(pathToFile2);
    if (!inpfl.is_open())
    {
        Info << "No liquid data found" << endl;
    }
    getline(inpfl, line);
    double HLiq11[buffer], HLiq12[buffer], HLiq22[buffer];
    np = 0;

    while (getline(inpfl, line))
    {
        string line_value;
        std::istringstream ss(line);
        getline(ss, line_value, ',');
        T1[np] = stod(line_value);
        getline(ss, line_value, ',');
        HLiq11[np] = stod(line_value);
        getline(ss, line_value, ',');
        HLiq22[np] = stod(line_value);
        getline(ss, line_value);
        HLiq12[np] = stod(line_value);
    
    //Info << T1[np] << " " << HLiq11[np] << " " << HLiq22[np] << " " << HLiq12[np] << endl;
    //cout << T1[np] << " " << ALiq[np] << endl;
            
    ++np;
    }

    inpfl.close();
    
    const fileName pathToFile3 = "constant/Composition_FCC_A1.csv";
    ifstream inpfc(pathToFile3);
    if (!inpfc.is_open())
    {
        Info << "No composition data found" << endl;
    }
    getline(inpfc, line);
    double cSol1[buffer], cSol2[buffer], cLiq1[buffer], cLiq2[buffer];
    np = 0;
    
    while (getline(inpfc, line))
    {
        string line_value;
        std::istringstream ss(line);
        getline(ss, line_value, ',');
        T1[np] = stod(line_value);
        getline(ss, line_value, ',');
        cSol1[np] = stod(line_value);
        getline(ss, line_value, ',');
        cSol2[np] = stod(line_value);
        getline(ss, line_value, ',');
        cLiq1[np] = stod(line_value);
        getline(ss, line_value);
        cLiq2[np] = stod(line_value);
    
    //Info << T1[np] << " " << cSol1[np] << " " << cSol2[np] << " " << cLiq1[np] << " " << cLiq2[np] << endl;
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
    
