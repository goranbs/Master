#include "moleculesystem.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <atom.h>

using namespace std;

MoleculeSystem::MoleculeSystem(){
}

void MoleculeSystem::Cout_Info(){
    cout << "Useage: " << endl;
    cout << "MoleculeSystem <system>;" << endl;
    cout << "Setup_<choose_system>(unit_cells_x_dir, unit_cells_y_dir, unit_cells_z_dir);" << endl;
    cout << "how to get the <system> file?" << endl;
}

void MoleculeSystem::Setup_Calcium(int n_unit_cells_x, int n_unit_cells_y, int n_unit_cells_z){

    double unit_cell_length_x, unit_cell_length_y, unit_cell_length_z;
    unit_cell_length_x = 2.3;
    unit_cell_length_y = 2.3;
    unit_cell_length_z = 2.3;

    vector <double> r = {0,0,0};
    vector <double> v = r;
    vector <double> f = v;
    double u = 0;

    vector <Atom> Calcium;
    Atom Ca(r, v, f, u);
    Ca.set_type("Ca");
    Calcium.push_back(Ca);

    Initialize(Calcium, n_unit_cells_x, n_unit_cells_y, n_unit_cells_z, unit_cell_length_x, unit_cell_length_y, unit_cell_length_z);
    Write_Initial_State(Calcium, "calcium.txt");
}

void MoleculeSystem::Setup_Kaolinite(int n_unit_cells_x, int n_unit_cells_y, int n_unit_cells_z){

    double unit_cell_length_x, unit_cell_length_y, unit_cell_length_z, unit_cell_angle_alpha, unit_cell_angle_beta, unit_cell_angle_gamma;

    unit_cell_length_x = 5.154;
    unit_cell_length_y = 8.942;
    unit_cell_length_z = 7.401;
    unit_cell_angle_alpha = 91.69;
    unit_cell_angle_beta = 104.65;
    unit_cell_angle_gamma = 89.82;

    vector <Atom> kaolinite; // holds the atoms in the unit cell of kaolinite

    vector <double> Al1 = {0.6650, 4.3300, 3.4040};
    vector <double> Al2 = {3.2140, 2.8550, 3.3960};
    vector <double> Si1 = {4.9750, 3.0050, 0.6610};
    vector <double> Si2 = {2.4620, 1.4720, 0.6710};
    vector <double> O1 =  {-0.3210, 3.0970, 2.2630};
    vector <double> O2 =  {0.0550, 5.8590, 2.2660};
    vector <double> O3 =  {0.0140, 4.4710, 0.0000};
    vector <double> O4 =  {1.0450, 2.0680, 0.1750};
    vector <double> O5 =  {1.0710, 6.8310, 0.0020};
    vector <double> O6 =  {-0.3200, 8.5930, 2.3290};
    vector <double> O7 =  {3.8190, 1.3530, 4.3260};
    vector <double> O8 =  {-0.9220, 4.1030, 4.3240};
    vector <double> O9 =  {-0.9230, 7.5290, 4.3520};
    vector <double> H1 = {0.0860, 0.2420, 2.4870};
    vector <double> H2 = {-1.0150, 1.4610, 5.0180};
    vector <double> H3 = {-1.1230, 4.1950, 5.068};
    vector <double> H4 = {-1.1110, 6.9610, 4.9970};
    vector <double> Al3 = {3.2560, 8.8010, 3.4040};
    vector <double> Al4 = {0.6510, 7.3260, 3.3960};
    vector <double> Si3 = {2.4120, 7.4760, 0.6610};
    vector <double> Si4 = {-0.1010, 5.9430, 0.6710};
    vector <double> O10 = {2.2700, 7.5680, 2.2630};
    vector <double> O11 = {2.6180, 1.3880, 2.2660};
    vector <double> O12 = {2.5770, 0.0000, 0.0000};
    vector <double> O13 = {3.6360, 6.5390, 0.1750};
    vector <double> O14 = {3.6340, 2.3600, 0.0020};
    vector <double> O15 = {2.2430, 4.1220, 2.3290};
    vector <double> O16 = {1.2560, 5.8240, 4.3260};
    vector <double> O17 = {1.6690, 8.5740, 4.3240};
    vector <double> O18 = {1.6400, 3.0580, 4.3520};
    vector <double> H5 = {2.6770, 4.7130, 2.4870};
    vector <double> H6 = {1.5760, 5.9320, 5.0180};
    vector <double> H7 = {1.4690, 8.6660, 5.0680};
    vector <double> H8 = {1.4520, 2.4900, 4.9970};

    // If you want to write the above directly to a file:
/*
    // test writing to file:
    ofstream myfile;
    myfile.open("kaolinite.txt");
    myfile << 34 << endl;
    myfile << 0 << endl;
    myfile << "Al" << " " << Al1[0] << " " << Al1[1] << " " << Al1[2] << endl;
    myfile << "Al" << " " << Al2[0] << " " << Al2[1] << " " << Al2[2] << endl;
    myfile << "Si" << " " << Si1[0] << " " << Si1[1] << " " << Si1[2] << endl;
    myfile << "Si" << " " << Si2[0] << " " << Si2[1] << " " << Si2[2] << endl;
    myfile << "O" << " " << O1[0] << " " << O1[1] << " " << O1[2] << endl;
    myfile << "O" << " " << O2[0] << " " << O2[1] << " " << O2[2] << endl;
    myfile << "O" << " " << O3[0] << " " << O3[1] << " " << O3[2] << endl;
    myfile << "O" << " " << O4[0] << " " << O4[1] << " " << O4[2] << endl;
    myfile << "O" << " " << O5[0] << " " << O5[1] << " " << O5[2] << endl;
    myfile << "O" << " " << O6[0] << " " << O6[1] << " " << O6[2] << endl;
    myfile << "O" << " " << O7[0] << " " << O7[1] << " " << O7[2] << endl;
    myfile << "O" << " " << O8[0] << " " << O8[1] << " " << O8[2] << endl;
    myfile << "O" << " " << O9[0] << " " << O9[1] << " " << O9[2] << endl;
    myfile << "H" << " " << H1[0] << " " << H1[1] << " " << H1[2] << endl;
    myfile << "H" << " " << H2[0] << " " << H2[1] << " " << H2[2] << endl;
    myfile << "H" << " " << H3[0] << " " << H3[1] << " " << H3[2] << endl;
    myfile << "H" << " " << H4[0] << " " << H4[1] << " " << H4[2] << endl;
    myfile << "Al" << " " << Al3[0] << " " << Al3[1] << " " << Al3[2] << endl;
    myfile << "Al" << " " << Al4[0] << " " << Al4[1] << " " << Al4[2] << endl;
    myfile << "Si" << " " << Si3[0] << " " << Si3[1] << " " << Si3[2] << endl;
    myfile << "Si" << " " << Si4[0] << " " << Si4[1] << " " << Si4[2] << endl;
    myfile << "O" << " " << O10[0] << " " << O10[1] << " " << O10[2] << endl;
    myfile << "O" << " " << O11[0] << " " << O11[1] << " " << O11[2] << endl;
    myfile << "O" << " " << O12[0] << " " << O12[1] << " " << O12[2] << endl;
    myfile << "O" << " " << O13[0] << " " << O13[1] << " " << O13[2] << endl;
    myfile << "O" << " " << O14[0] << " " << O14[1] << " " << O14[2] << endl;
    myfile << "O" << " " << O15[0] << " " << O15[1] << " " << O15[2] << endl;
    myfile << "O" << " " << O16[0] << " " << O16[1] << " " << O16[2] << endl;
    myfile << "O" << " " << O17[0] << " " << O17[1] << " " << O17[2] << endl;
    myfile << "O" << " " << O18[0] << " " << O18[1] << " " << O18[2] << endl;
    myfile << "H" << " " << H5[0] << " " << H5[1] << " " << H5[2] << endl;
    myfile << "H" << " " << H6[0] << " " << H6[1] << " " << H6[2] << endl;
    myfile << "H" << " " << H7[0] << " " << H7[1] << " " << H7[2] << endl;
    myfile << "H" << " " << H8[0] << " " << H8[1] << " " << H8[2] << endl;

    myfile.close();
    */

    vector <double> v (3,0.0);
    double p = 0.0;

    Atom new_atom(Al1, v, v, p);  // #1
    new_atom.set_type("Al");
    kaolinite.push_back(new_atom);

    new_atom.set_position(Al2);  // #2
    new_atom.set_type("Al");
    kaolinite.push_back(new_atom);

    new_atom.set_position((Si1));  // #3

    new_atom.set_type("Si");
    kaolinite.push_back(new_atom);

    new_atom.set_position(Si2);    // #4
    new_atom.set_type("Si");
    kaolinite.push_back(new_atom);

    new_atom.set_position(O1);      // #5
    new_atom.set_type("O");
    kaolinite.push_back(new_atom);

    new_atom.set_position(O2);       // #6
    new_atom.set_type("O");
    kaolinite.push_back(new_atom);

    new_atom.set_position(O3);       // #7
    new_atom.set_type("O");
    kaolinite.push_back(new_atom);

    new_atom.set_position(O4);       // #8
    new_atom.set_type("O");
    kaolinite.push_back(new_atom);

    new_atom.set_position(O5);       // #9
    new_atom.set_type("O");
    kaolinite.push_back(new_atom);

    new_atom.set_position(O6);      // #10
    new_atom.set_type("O");
    kaolinite.push_back(new_atom);

    new_atom.set_position(O7);       // #11
    new_atom.set_type("O");
    kaolinite.push_back(new_atom);

    new_atom.set_position(O8);      // #12
    new_atom.set_type("O");
    kaolinite.push_back(new_atom);

    new_atom.set_position(O9);      // #13
    new_atom.set_type("O");
    kaolinite.push_back(new_atom);

    new_atom.set_position(H1);        // #14
    new_atom.set_type("H");
    kaolinite.push_back(new_atom);

    new_atom.set_position(H2);       // #15
    new_atom.set_type("H");
    kaolinite.push_back(new_atom);

    new_atom.set_position(H3);        // #16
    new_atom.set_type("H");
    kaolinite.push_back(new_atom);

    new_atom.set_position(H4);       // #17
    new_atom.set_type("H");
    kaolinite.push_back(new_atom);

    new_atom.set_position(Al3);       // #18
    new_atom.set_type("Al");
    kaolinite.push_back(new_atom);

    new_atom.set_position(Al4);       // #19
    new_atom.set_type("Al");
    kaolinite.push_back(new_atom);

    new_atom.set_position(Si3);       // #20
    new_atom.set_type("Si");
    kaolinite.push_back(new_atom);

    new_atom.set_position(Si4);      // #21
    new_atom.set_type("Si");
    kaolinite.push_back(new_atom);

    new_atom.set_position(O10);      // #22
    new_atom.set_type("O");
    kaolinite.push_back(new_atom);

    new_atom.set_position(O11);      // #23
    new_atom.set_type("O");
    kaolinite.push_back(new_atom);

    new_atom.set_position(O12);      // #24
    new_atom.set_type("O");
    kaolinite.push_back(new_atom);

    new_atom.set_position(O13);      // #25
    new_atom.set_type("O");
    kaolinite.push_back(new_atom);

    new_atom.set_position(O14);      // #26
    new_atom.set_type("O");
    kaolinite.push_back(new_atom);

    new_atom.set_position(O15);      // #27
    new_atom.set_type("O");
    kaolinite.push_back(new_atom);

    new_atom.set_position(O16);      // #28
    new_atom.set_type("O");
    kaolinite.push_back(new_atom);

    new_atom.set_position(O17);      // #29
    new_atom.set_type("O");
    kaolinite.push_back(new_atom);

    new_atom.set_position(O18);      // #30
    new_atom.set_type("O");
    kaolinite.push_back(new_atom);

    new_atom.set_position(H5);        // #31
    new_atom.set_type("H");
    kaolinite.push_back(new_atom);

    new_atom.set_position(H6);        // #32
    new_atom.set_type("H");
    kaolinite.push_back(new_atom);

    new_atom.set_position(H7);        // #33
    new_atom.set_type("H");
    kaolinite.push_back(new_atom);

    new_atom.set_position(H8);        // #34
    new_atom.set_type("H");
    kaolinite.push_back(new_atom);

    Initialize(kaolinite, n_unit_cells_x, n_unit_cells_y, n_unit_cells_z, unit_cell_length_x, unit_cell_length_y, unit_cell_length_z);
    Write_Initial_State(kaolinite, "kaolinite.txt");

}

void MoleculeSystem::Initialize(vector <Atom> &atoms, int &lx, int &ly, int &lz, double &xi, double &yi, double &zi){

    vector <double> ri (3,0.0);
    vector <double> r (3,0.0);
    string atom_type;
    int N_unit_cell_atoms = atoms.size();
    Atom new_atom = atoms[0];

    vector <Atom> AtomContainer;

    // unit cell loop:
    for (int i=0; i<lx; i++){          // unit cell x
        for (int j=0; j<ly; j++){      // unit cell y
            for (int k=0; k<lz; k++){  // unit cell z

                ri[0] = i*xi; ri[1] = j*yi; ri[2] = k*zi;

                for (int atom=0; atom<N_unit_cell_atoms ;atom++){
                    r = atoms[atom].get_position();
                    atom_type = atoms[atom].get_type();

                    for (int cor=0; cor<3; cor++){
                        r[cor] = r[cor] + ri[cor];
                    }
                    atom_type = atoms[atom].get_type();
                    new_atom.set_position(r);
                    new_atom.set_type(atom_type);
                    AtomContainer.push_back(new_atom);
                }
            }
        }
    }
    //cout << AtomContainer.size() << endl;
    atoms = AtomContainer;
}

void MoleculeSystem::Write_Initial_State(vector <Atom> &atoms, string filename){

    vector <double> r (3,0.0);

    ofstream myfile;
    myfile.open(filename);
    myfile << atoms.size() << endl;          // number of atoms
    myfile << 0 << endl;                      // timestep 0

    for (int i=0; i<atoms.size(); i++){
        r = atoms[i].get_position();
        myfile << atoms[i].get_type() << " " << r[0] << " " << r[1] << " " << r[2] << endl;
    }
    myfile.close();

}
