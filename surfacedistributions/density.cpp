#include "density.h"

density::density(
        vector<atom> Atoms, vector <int> &types, vector <double> &molarmass, int &Natoms,
        vector <double> systemsize,string &dir, ReadLAMMPSdump &PATH){

    if (dir == "xdir"){
        DIR = 0;
        Area = (systemsize[3] - systemsize[2])*(systemsize[5] - systemsize[4]);
        systemLo = systemsize[0];  // xlo
        systemHi = systemsize[1];  // xhi
    }
    if (dir == "ydir"){
        DIR = 1;
        Area = (systemsize[1] - systemsize[0])*(systemsize[5] - systemsize[4]);
        systemLo = systemsize[2];  // ylo
        systemHi = systemsize[3];  // yhi
    }
    if (dir == "zdir"){
        DIR = 2;
        Area = (systemsize[1] - systemsize[0])*(systemsize[3] - systemsize[2]);
        systemLo = systemsize[4]; // zlo
        systemHi = systemsize[5]; // zhi
    }

    TYPES = types;                  // types to follow
    MOLARMASS = molarmass;          // molar masses of types
    Ntypes = TYPES.size();          // number of atom types to follow

    int Type;
    int boxnr = 0;
    double Pos = 0;
    vector <double> pos (3,0.0);

    L = (systemHi - systemLo);      // whole system size
    Nboxes = int(L/dL);             // number of boxes to divide the system into

    vector <vector <int> > container1 (Nboxes,vector <int> (Ntypes,0.0));       // Number of atoms of every type in box
    vector <vector <double> > container2 (Nboxes,vector <double> (Ntypes,0.0)); // atoms/boxvolume -//-

    for (int i=0;i<Natoms;i++){
        Type = Atoms[i].get_atom_type_number();
        if (isintypes(Type) == true){
            pos = Atoms[i].get_position();
            Pos = pos[DIR] - systemLo;        // adjust system so that we start at the origin!

            if (Pos < 0){
                //cout << "Neg position! Pos= " << Pos << " systemLo= " << systemLo <<  endl;
                Pos = 0;
            }

            boxnr = floor(Pos/dL);

            if (boxnr >= Nboxes){
                //cout << "boxnr= (" << boxnr << "/" << (Nboxes -1) << ") pos= " << Pos << " systemHi= " << (systemHi - systemLo) <<  " dL= " << dL << " pos/dL= " << (Pos/dL) << endl;
                boxnr = Nboxes -1;
            }

            container1[boxnr][Type-1] = container1[boxnr][Type-1] + 1;


        }
    }

    double boxvolume = dL*Area;
    for (int i=0;i<Nboxes;i++){
        for (int j=0; j<Ntypes;j++){
            container2[i][j] = container1[i][j]/boxvolume;
        }
    }

    cont1 = container1;
    cont2 = container2;

    string path = PATH.get_path();              // get path of database
    string filename;                            //
    int tid = PATH.get_timestep();              // return timestep of density evaluation
    stringstream ss;                            //
    ss << "density_results." << tid << ".txt";  // create filename
    ss >> filename;                             // stream from stringstream into filename
    writetofile(filename,path);                 // write to file

}

void density::writetofile(string &filename, string &path){

    ofstream thisfile;          // write to file
    thisfile.open(filename);
    thisfile << "# Densities in system: types: ";                                                       // line 1
    for (int i=0; i<Ntypes; i++){ thisfile << TYPES[i] << " " ;}
    thisfile << endl;
    thisfile << "# Database: " << path << endl;                                                         // line 2
    thisfile << "# time:  " << endl;                                                                    // line 3
    thisfile << "# cleavage: L dL Nboxes " << L << " " << dL << " " << Nboxes << endl;                  // line 4
    thisfile << "# Systemsize: systemLo systemHi " << systemLo << " " << systemHi << endl;              // line 5
    thisfile << "# boxID Type1 Nparticles1 Nparticles1/volume rho1 Type2 Nparticles2 Nparticles2/volume rho2..." << endl;                   // line 6

    for (int i=0; i<Nboxes; i++ ){
        thisfile << i << " ";
        for (int j=0;j<Ntypes;j++){
            // assuming that molarmass is in g/mol and that the dimentions of the system is in Angstrom
            thisfile << TYPES[j] << " " << cont1[i][j] << " " << cont2[i][j] << " " << (cont2[i][j]*MOLARMASS[j]/0.6022) << " ";
        }
        thisfile << endl;
    }
    thisfile.close();
}


bool density::isintypes(int &type){
    // if type is in types, return true!
    truthvalue = false;
    for (int i=0;i<Ntypes;i++){
        if (type == TYPES[i]){
            truthvalue = true;
        }
    }
    return truthvalue;
}
