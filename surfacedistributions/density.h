#ifndef DENSITY_H
#define DENSITY_H

#include <iostream>
#include <time.h>
#include <math.h>
#include "atom.h"
#include "readlammpsdump.h"

class density
{
public:
    density(vector<atom> Atoms, vector <int> &types, vector<double> &molarmass, int &Natoms,
            vector<double> systemsize, string &dir, ReadLAMMPSdump &PATH);

private:
    bool isintypes(int &type);
    void writetofile(string &filename, string &path);
    int DIR;                             // direction in which to divide system
    int Ntypes;                          // Number of types to follow
    int Nboxes;                          // Number of boxes to create along DIR
    double L;                            // system length (defined from SystemStartValue and SystemEndValue) along DIR
    double dL = 0.5;                     // spacing of boxed along DIR
    double center;                       // center pos of cleavage (SystemStartValue + L/2)
    double Area;                         // area of box along surface normal
    double systemLo, systemHi;           // DIR origin and end of system in direction DIR
    bool truthvalue;                     // true or false. Test type in TYPES
    vector <int> TYPES;                  // Types to follow
    vector <vector <int> > cont1;        // number of atoms in box for every type
    vector <vector <double> > cont2;     // -//- divided by volume of box
    vector <double> MOLARMASS;           // molar masses of types


};

#endif // DENSITY_H
