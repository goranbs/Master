#ifndef DIFFUSION_H
#define DIFFUSION_H

#include <iostream>
#include <vector>
#include <stdexcept>
#include <math.h>
#include <sstream>
#include <time.h>
#include "atom.h"
#include "readlammpsdump.h"

class diffusion
{
public:
    diffusion(vector<atom> &atoms, int &type, int &Natoms, double &SystemStartValue, double &SystemEndValue, string &dir, double &DL);
    int counted_particles();
    vector <int> particle_distribution();
    void calculate_diffusion(vector <atom> &atoms, ReadLAMMPSdump &PATH, string &frontname, double &Lx, double &Ly, double &Lz, int &dt, int &Tstart, int &Tend);

private:
    void MinimumImageConvention(double &delta, double &systemsize);
    void write_to_file(string &filename, string &path, vector <vector <double> > &MSD);

    int DIR;                  // surface normal cleavage
    int Counted;              // number of atoms of Type counted
    int NATOMS;               // total number of atoms
    int Type;                 // atom type to follow
    int Nboxes;               // number of boxes
    int Ntimesteps;           // number of time steps
    int TimeStart;            // timestep after timestep at origin!!! origin = TimeStart - Dt
    int TimeEnd;              // last timestep
    int Dt;                   // timestep
    double L;                 // L/2 = center
    double dL;          // spacing bins
    double center;            // center of cleavage
    double Start;             // System start value
    double End;               // System end value
    vector <int> Container;   // number of particles in every container at the beginning
    vector <int> Container_e; // number of particles in every container at the end
};

#endif // DIFFUSION_H
