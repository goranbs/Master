#include <iostream>
#include <vector>
#include "atom.h"
#include "readlammpsdump.h"
#include "diffusion.h"
#include "density.h"

using namespace std;

int main(int argc, char * argv[]){

    int Tstart, Tend, dt, typetofollow, Natoms;
    double Lx,Ly,Lz,Lower,Upper,dL;
    string path, filename, frontname, dir, calculate_diffusion, calculate_density;

    if (argc < 14){
        cout << "ERROR! Number of command line arguments is < 10 " << endl;
        for (int arg=0; arg<argc; arg++){ cout << arg << " " << argv[arg] << endl; }
        cout << "SYSTEM EXIT..." << endl;
        return 0;
    }

    // path_executable = argv[0];
    Tstart = atoi(argv[1]);
    Tend = atoi(argv[2]);
    dt = atoi(argv[3]);
    typetofollow = atoi(argv[4]);
    dir = argv[5];
    path = argv[6];
    frontname = argv[7];
    Lower = atof(argv[8]);
    Upper = atof(argv[9]);
    dL = atof(argv[10]);
    calculate_diffusion = argv[11];
    calculate_density = argv[12];
    vector < double > molarmass;
    vector < int > types;
    int atype = 1;
    for (int arg=13; arg<argc; arg++){
        types.push_back(atype);
        atype = atype + 1;
        molarmass.push_back(atof(argv[arg]));
    }

    cout << "Tstart                = " << Tstart << endl;
    cout << "Tend                  = " << Tend << endl;
    cout << "dt                    = " << dt << endl;
    cout << "typetofollow          = " << typetofollow << endl;
    cout << "dir                   = " << dir << endl;
    cout << "path                  = " << path << endl;
    cout << "frontname             = " << frontname << endl;
    cout << "Lower surface         = " << Lower << endl;
    cout << "Upper surface         = " << Upper << endl;
    cout << "bin size dL           = " << dL << endl;
    cout << "calculate_diffusion   = " << calculate_diffusion << endl;
    cout << "calculate_density     = " << calculate_density << endl;
    for (int i=0; i<molarmass.size();i++){
        cout << "molarmass " << i << "           = " <<  molarmass[i] << endl;
    }

    // ------------------------- Read the initial file ------------------------------

    ReadLAMMPSdump PATH(path);                            // set path
    vector <atom> Atoms;                                  // holds atom information

    filename = frontname + "." + to_string(Tstart) + ".txt";   // should find a better way!!!
    cout << "initial file filename = " << filename << endl;

    PATH.read_initialfile(filename, Atoms);               // filename = file-origo for one estimate!
    Natoms = PATH.get_Natoms();                           // number of atoms in system

    vector <double> systemsize = PATH.get_systemsize();   // get and set the system size
    Lx = systemsize[1] - systemsize[0];
    Ly = systemsize[3] - systemsize[2];
    Lz = systemsize[5] - systemsize[4];

    cout << "# ----------------------------------------------------------------- #" << endl;
    cout << " Number of atoms: " << Natoms << endl;
    cout << " System size: " << Lx << " " << Ly << " " << Lz << endl;         // print the system size
    cout << "# ----------------------------------------------------------------- #" << endl;
     // ----------------------- Run the diffusion calculation if calculate diffusion == yes -------

    if (calculate_diffusion == "yes"){
        cout << "Calculating the diffusion..." << endl;
        if (Lower == Upper){
            // then set upper and lower surface pos to system boundaries
            if (dir == "xdir"){
                Lower = systemsize[0];  // xlo
                Upper = systemsize[1];  // xhi
            }
            if (dir == "ydir"){
                Lower = systemsize[2];  // ylo
                Upper = systemsize[3];  // yhi
            }
            if (dir == "zdir"){
                Lower = systemsize[4]; // zlo
                Upper = systemsize[5]; // zhi
            }
        }
        diffusion diffuse(Atoms,typetofollow,Natoms,Lower,Upper,dir,dL);                   // initialize diffusion calc from this origin
        diffuse.calculate_diffusion(Atoms,PATH,frontname,Lx,Ly,Lz,dt,Tstart,Tend);         // calculate diffusion
        cout << "# ----------------------------------------------------------------- #" << endl <<  "Done calculating the diffusion!" << endl;
    }

    if (calculate_density == "yes"){
        cout << "Calculating the density of the system at timestep T =  " << PATH.get_timestep() << endl;
        density systemdensity(Atoms,types,molarmass,Natoms,systemsize,dir,PATH);
    }
    cout << "# --------------------------------------------------------------------------------------------------- #" << endl;
    return 0;
}

