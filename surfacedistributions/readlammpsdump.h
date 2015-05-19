#ifndef READLAMMPSDUMP_H
#define READLAMMPSDUMP_H

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <time.h>
#include "atom.h"

using namespace std;

class ReadLAMMPSdump
{
public:
    ReadLAMMPSdump(string &path);
    void read_initialfile(string &filename, vector <atom> &Atoms);
    void read_newfile(string &filename, vector <atom> &Atoms);
    vector <string> split_on_whitespace(const string &str);
    int get_Natoms();
    string get_path();
    vector <double> get_systemsize();
    double get_timeusage();
    int get_timestep();

private:
    string PATH;
    string ENDCHAR;
    vector <int> TIME;
    vector <int> NATOMS;
    double Xlo,Xhi,Ylo,Yhi,Zlo,Zhi;
    clock_t CLOCK;
    double timeusage;
};

#endif // READLAMMPSDUMP_H
