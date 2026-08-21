#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cstdio>

using namespace std;

string getfile(ifstream& is)
{
    string contents;

    for (char ch; is.get(ch); contents.push_back(ch)) {}

    return contents;
}

void find_and_replace(string& file_contents,
                      const string& morn,
                      const string& night)
{
    size_t pos = file_contents.find(morn);

    while (pos != string::npos)
    {
        file_contents.replace(pos, morn.length(), night);

        pos = file_contents.find(morn, pos);
    }
}


// ---------------------------------------------------------
// Swap FILLCUBE:
// (phase, x_lo, y_lo, z_lo, x_hi, y_hi, z_hi)
//
// becomes:
//
// (phase, y_lo, x_lo, z_lo, y_hi, x_hi, z_hi)
// ---------------------------------------------------------
string swapFillCubeXY(string line)
{
    if (line.find("FILLCUBE") == string::npos)
        return line;

    size_t start = line.find('{');
    size_t end = line.find('}', start);

    if (start == string::npos || end == string::npos)
        return line;

    string inside = line.substr(start + 1, end - start - 1);

    vector<string> values;
    string value;

    stringstream ss(inside);

    while (getline(ss, value, ','))
    {
        // Remove leading whitespace
        value.erase(0, value.find_first_not_of(" \t"));

        // Remove trailing whitespace
        value.erase(value.find_last_not_of(" \t") + 1);

        values.push_back(value);
    }

    // phase, x_lo, y_lo, z_lo, x_hi, y_hi, z_hi
    if (values.size() == 7)
    {
        swap(values[1], values[2]);
        swap(values[4], values[5]);

        string newInside;

        for (size_t i = 0; i < values.size(); ++i)
        {
            if (i > 0)
                newInside += ", ";

            newInside += values[i];
        }

        line.replace(start + 1,
                     end - start - 1,
                     newInside);
    }

    return line;
}


string swapBoundaryXY(string line)
{
    if (line.find("BOUNDARY") == string::npos)
        return line;

    size_t start = line.find('{');
    size_t end = line.find('}', start);

    if (start == string::npos || end == string::npos)
        return line;

    string inside = line.substr(start + 1, end - start - 1);

    vector<string> values;
    string value;

    stringstream ss(inside);

    while (getline(ss, value, ','))
    {
        // Remove leading whitespace
        value.erase(0, value.find_first_not_of(" \t"));

        // Remove trailing whitespace
        value.erase(value.find_last_not_of(" \t") + 1);

        values.push_back(value);
    }

    // BOUNDARY:
    //
    // 0 = phase
    // 1 = x_lo
    // 2 = x_hi
    // 3 = y_lo
    // 4 = y_hi
    // 5 = z_lo
    // 6 = z_hi

    if (values.size() == 7)
    {
        swap(values[1], values[3]);  // x_lo <-> y_lo
        swap(values[2], values[4]);  // x_hi <-> y_hi

        string newInside;

        for (size_t i = 0; i < values.size(); ++i)
        {
            if (i > 0)
                newInside += ", ";

            newInside += values[i];
        }

        line.replace(start + 1,
                     end - start - 1,
                     newInside);
    }

    return line;
}


// ---------------------------------------------------------
// Main
// ---------------------------------------------------------
int main(int argc, char *argv[])
{
    string s1, s2;

    int num{0};

    ifstream ifile1(argv[1]);
    ifstream ifile2(argv[2]);

    // -----------------------------------------------------
    // First read the filling file completely.
    // Determine whether FILLCUBE exists.
    // -----------------------------------------------------

    vector<string> fillingLines;

    bool hasFillCube = false;

    while (getline(ifile2, s2))
    {
        fillingLines.push_back(s2);

        if (s2.find("FILLCUBE") != string::npos)
        {
            hasFillCube = true;
        }
    }

    // -----------------------------------------------------
    // Create input1.in
    // -----------------------------------------------------

    ofstream ofile("input1.in");

    string meshXLine;
    string meshYLine;

    // -----------------------------------------------------
    // Read input file
    // -----------------------------------------------------

    while(getline(ifile1, s1))
{
    if(s1 == "NUMPHASES = 2;")
    {
        num = 2;
    }

    if((s1 == "Gamma_abc = {};") && (num == 2))
    {
        continue;
    }

    // ---------------------------------------------
    // If FILLCUBE exists in filling file,
    // swap MESH_X and MESH_Y.
    // Keep them in their original positions.
    // ---------------------------------------------
    if(hasFillCube)
    {
        if(s1.find("MESH_X") != string::npos)
        {
            size_t pos = s1.find("MESH_X");

            if(pos != string::npos)
            {
                s1.replace(pos, 6, "TEMP_MESH_X");
            }
        }
        else if(s1.find("MESH_Y") != string::npos)
        {
            size_t pos = s1.find("MESH_Y");

            if(pos != string::npos)
            {
                s1.replace(pos, 6, "MESH_X");
            }
        }

        // Change temporary MESH_X into MESH_Y
        size_t pos = s1.find("TEMP_MESH_X");

        if(pos != string::npos)
        {
            s1.replace(pos, 11, "MESH_Y");
        }
    }

    // ---------------------------------------------
    // Swap BOUNDARY x/y coordinates
    // ---------------------------------------------
    s1 = swapBoundaryXY(s1);

    ofile << s1 << endl;
}

    // -----------------------------------------------------
    // If FILLCUBE exists:
    //
    // MESH_X line becomes MESH_Y
    // MESH_Y line becomes MESH_X
    // -----------------------------------------------------

    if (hasFillCube)
    {
        if (!meshXLine.empty())
        {
            // Change MESH_X to MESH_Y
            size_t pos = meshXLine.find("MESH_X");

            meshXLine.replace(pos, 6, "MESH_Y");

            ofile << meshXLine << endl;
        }

        if (!meshYLine.empty())
        {
            // Change MESH_Y to MESH_X
            size_t pos = meshYLine.find("MESH_Y");

            meshYLine.replace(pos, 6, "MESH_X");

            ofile << meshYLine << endl;
        }
    }

    // -----------------------------------------------------
    // Add filling file
    // -----------------------------------------------------

    for (string line : fillingLines)
    {
        if (hasFillCube)
        {
            line = swapFillCubeXY(line);
        }

        ofile << line << endl;
    }

    ifile1.close();
    ifile2.close();
    ofile.close();


    // -----------------------------------------------------
    // Create input2.in
    // -----------------------------------------------------

    ifstream filein("input1.in");
    ofstream fileout("input2.in");

    string contents = getfile(filein);

    find_and_replace(contents, "{", " ");
    find_and_replace(contents, ",", " ");
    find_and_replace(contents, "}", " ");
    find_and_replace(contents, ";", " ");

    fileout << contents;

    filein.close();
    fileout.close();

    remove("input1.in");

    return 0;
}