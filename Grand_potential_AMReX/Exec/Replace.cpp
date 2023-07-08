#include <iostream>
#include <fstream>
#include <string>
#include<bits/stdc++.h>

using namespace std;

std::string getfile(std::ifstream& is) {
  std::string contents;
  // Here is one way to read the whole file
  for (char ch; is.get(ch); contents.push_back(ch)) {}
  return contents;
}

void find_and_replace(std::string& file_contents, 
    const std::string& morn, const std::string& night) {
  // This searches the file for the first occurence of the morn string.
  size_t pos = file_contents.find(morn);
  while (pos != std::string::npos) {
    file_contents.replace(pos, morn.length(), night);
    // Continue searching from here.
    pos = file_contents.find(morn, pos);
  }
}

int main(int argc, char *argv[])
{   
    string s1,s2;
    ifstream ifile1(argv[1]);
    ifstream ifile2(argv[2]);
    ofstream ofile("input1.in");
    
    while(getline(ifile1,s1)){
    	
    	ofile<<s1;
    	ofile<<endl;
    }
    
    while(getline(ifile2,s2)){
    	
    	ofile<<s2;
    	ofile<<endl;
    }
    
    ifile1.close();
    ifile2.close();
    ofile.close();
   
    
    ifstream filein("input1.in");
    ofstream fileout("input2.in");
    std::string contents = getfile(filein);
    find_and_replace(contents, "{", "(");
    find_and_replace(contents, "}", ")");
    fileout << contents;

    filein.close();
    fileout.close();

    
    
   
}
