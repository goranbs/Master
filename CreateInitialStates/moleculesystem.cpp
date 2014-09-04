#include "moleculesystem.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <atom.h>
#include <cmath>

double PI = 3.14159265;

using namespace std;

MoleculeSystem::MoleculeSystem(){
}

void MoleculeSystem::Cout_Info(){
    cout << "Useage: " << endl;
    cout << "MoleculeSystem <system>;" << endl;
    cout << "Setup_<choose_system>(unit_cells_x_dir, unit_cells_y_dir, unit_cells_z_dir);" << endl;
    cout << "how to get the <system> file?" << endl;
}

void MoleculeSystem::Setup_calciumhydroxide(int n_unit_cells_x, int n_unit_cells_y, int n_unit_cells_z, string output){
    /* Portlandite is a oxide mineral mainly consisting of Calcium Hydroxide ( Ca(OH)_2 )
     *  A unit cell has according to Cygan et al. in Molecular Models of Hydroxide, Oxyhydroxide and Clay
     *  Phases and the Developent of a General Force Field, a,b,c specified as below. The angles between
     *  the axises is given as 120deg between a and b, and that the c-axis is located normal to the ab-plane.
     */

    double unit_cell_length_x = 3.589;      // a
    double unit_cell_length_y = 3.589;      // b
    double unit_cell_length_z = 4.911;      // c

    double OH_dist, Ca_OH_dist;
    OH_dist = 0.942;
    Ca_OH_dist = 2.369;

    vector <double> Ca = {0,0,0};
    vector <double> O1 = {1.1845,-2.0516,0};
    vector <double> O2 = {1.1845, 2.0516,0};
    vector <double> H1 = {1.1845,-2.0516,OH_dist};
    vector <double> H2 = {1.1845,2.0516,OH_dist};

    xlo = 0.0;
    ylo = -2.0;
    zlo = 0.0;
    xhi = 1.18;
    yhi = 2.05;
    zhi = 0.0;


    vector <double> v = {0,0,0}; // initial velocity
    vector <double> f = {0,0,0}; // initial force on atom
    double u = 0.0;              // initial potential
    int n_atom_types = 3;        // number of atom types

    vector <Atom> Portlandite;

    Atom new_atom(Ca,v,f,u);            // #1 Ca
    new_atom.set_type("Ca");
    new_atom.set_type_number(1);
    Portlandite.push_back(new_atom);

    new_atom.set_position(O1);          // #2 O
    new_atom.set_type("O");
    new_atom.set_type_number(2);
    Portlandite.push_back(new_atom);

    new_atom.set_position(O2);          // #3 O
    new_atom.set_type("O");
    new_atom.set_type_number(2);
    Portlandite.push_back(new_atom);

    new_atom.set_position(H1);          // #4 H
    new_atom.set_type("H");
    new_atom.set_type_number(3);
    Portlandite.push_back(new_atom);

    new_atom.set_position(H2);          // #5 H
    new_atom.set_type("H");
    new_atom.set_type_number(3);
    Portlandite.push_back(new_atom);

    // Ca = 1; O = 2; H = 3;
    masses.push_back(mass_calcium);
    masses.push_back(mass_oxygen);
    masses.push_back(mass_hydrogen);

    vector <Bond> bonds;
    Bond bond;
    bond.set_new_bond(1,1,1);
    bonds.push_back(bond);

    Initialize(Portlandite,n_unit_cells_x,n_unit_cells_y, n_unit_cells_z, unit_cell_length_x,unit_cell_length_y, unit_cell_length_z);
    if (output == "Ovito"){
        Write_Initial_State_Ovito(Portlandite, "portlandite_ovito.txt", n_atom_types);
    }
    else {
        Write_Initial_State_LAMMPS(Portlandite, "calciumhydroxide_LAMMPS.dat", n_atom_types);
    }

}

void MoleculeSystem::Setup_Calcium(int n_unit_cells_x, int n_unit_cells_y, int n_unit_cells_z, string output){

    double unit_cell_length_x, unit_cell_length_y, unit_cell_length_z;
    unit_cell_length_x = 2.3;
    unit_cell_length_y = 2.3;
    unit_cell_length_z = 2.3;

    vector <double> r = {0,0,0};
    vector <double> v = r;
    vector <double> f = v;
    double u = 0;
    int n_atom_types = 1;       // number of atom types

    vector <Atom> Calcium;
    Atom Ca(r, v, f, u);
    Ca.set_type("Ca");
    Ca.set_type_number(1);
    Calcium.push_back(Ca);

    // Ca = 1;
    masses.push_back(mass_calcium);

    xlo = 0;
    ylo = 0;
    zlo = 0;
    xhi = 0;
    yhi = 0;
    zhi = 0;

    vector <Bond> bonds;

    Initialize(Calcium, n_unit_cells_x, n_unit_cells_y, n_unit_cells_z, unit_cell_length_x, unit_cell_length_y, unit_cell_length_z);
    if (output == "Ovito"){
        Write_Initial_State_Ovito(Calcium, "calcium_ovito.txt", n_atom_types);
    }
    else {
        Write_Initial_State_LAMMPS(Calcium, bonds, "calcium_LAMMPS.dat", n_atom_types);
    }

}

void MoleculeSystem::Setup_Kaolinite(int n_unit_cells_x, int n_unit_cells_y, int n_unit_cells_z, string output){

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

    xlo = -0.321;
    xhi = 4.9;
    ylo = 1.4;
    yhi = 8.5;
    zlo = 0.0;
    zhi = 5.06;

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
    int n_atom_types = 4; // number of atom types

    Atom new_atom(Al1, v, v, p);  // #1
    new_atom.set_type("Al");
    new_atom.set_type_number(1);
    kaolinite.push_back(new_atom);

    new_atom.set_position(Al2);  // #2
    new_atom.set_type("Al");
    new_atom.set_type_number(1);
    kaolinite.push_back(new_atom);

    new_atom.set_position((Si1));  // #3

    new_atom.set_type("Si");
    new_atom.set_type_number(2);
    kaolinite.push_back(new_atom);

    new_atom.set_position(Si2);    // #4
    new_atom.set_type("Si");
    new_atom.set_type_number(2);
    kaolinite.push_back(new_atom);

    new_atom.set_position(O1);      // #5
    new_atom.set_type("O");
    new_atom.set_type_number(3);
    kaolinite.push_back(new_atom);

    new_atom.set_position(O2);       // #6
    new_atom.set_type("O");
    new_atom.set_type_number(3);
    kaolinite.push_back(new_atom);

    new_atom.set_position(O3);       // #7
    new_atom.set_type("O");
    new_atom.set_type_number(3);
    kaolinite.push_back(new_atom);

    new_atom.set_position(O4);       // #8
    new_atom.set_type("O");
    new_atom.set_type_number(3);
    kaolinite.push_back(new_atom);

    new_atom.set_position(O5);       // #9
    new_atom.set_type("O");
    new_atom.set_type_number(3);
    kaolinite.push_back(new_atom);

    new_atom.set_position(O6);      // #10
    new_atom.set_type("O");
    new_atom.set_type_number(3);
    kaolinite.push_back(new_atom);

    new_atom.set_position(O7);       // #11
    new_atom.set_type("O");
    new_atom.set_type_number(3);
    kaolinite.push_back(new_atom);

    new_atom.set_position(O8);      // #12
    new_atom.set_type("O");
    new_atom.set_type_number(3);
    kaolinite.push_back(new_atom);

    new_atom.set_position(O9);      // #13
    new_atom.set_type("O");
    new_atom.set_type_number(3);
    kaolinite.push_back(new_atom);

    new_atom.set_position(H1);        // #14
    new_atom.set_type("H");
    new_atom.set_type_number(4);
    kaolinite.push_back(new_atom);

    new_atom.set_position(H2);       // #15
    new_atom.set_type("H");
    new_atom.set_type_number(4);
    kaolinite.push_back(new_atom);

    new_atom.set_position(H3);        // #16
    new_atom.set_type("H");
    new_atom.set_type_number(4);
    kaolinite.push_back(new_atom);

    new_atom.set_position(H4);       // #17
    new_atom.set_type("H");
    new_atom.set_type_number(4);
    kaolinite.push_back(new_atom);

    new_atom.set_position(Al3);       // #18
    new_atom.set_type("Al");
    new_atom.set_type_number(1);
    kaolinite.push_back(new_atom);

    new_atom.set_position(Al4);       // #19
    new_atom.set_type("Al");
    new_atom.set_type_number(1);
    kaolinite.push_back(new_atom);

    new_atom.set_position(Si3);       // #20
    new_atom.set_type("Si");
    new_atom.set_type_number(2);
    kaolinite.push_back(new_atom);

    new_atom.set_position(Si4);      // #21
    new_atom.set_type("Si");
    new_atom.set_type_number(2);
    kaolinite.push_back(new_atom);

    new_atom.set_position(O10);      // #22
    new_atom.set_type("O");
    new_atom.set_type_number(3);
    kaolinite.push_back(new_atom);

    new_atom.set_position(O11);      // #23
    new_atom.set_type("O");
    new_atom.set_type_number(3);
    kaolinite.push_back(new_atom);

    new_atom.set_position(O12);      // #24
    new_atom.set_type("O");
    new_atom.set_type_number(3);
    kaolinite.push_back(new_atom);

    new_atom.set_position(O13);      // #25
    new_atom.set_type("O");
    new_atom.set_type_number(3);
    kaolinite.push_back(new_atom);

    new_atom.set_position(O14);      // #26
    new_atom.set_type("O");
    new_atom.set_type_number(3);
    kaolinite.push_back(new_atom);

    new_atom.set_position(O15);      // #27
    new_atom.set_type("O");
    new_atom.set_type_number(3);
    kaolinite.push_back(new_atom);

    new_atom.set_position(O16);      // #28
    new_atom.set_type("O");
    new_atom.set_type_number(3);
    kaolinite.push_back(new_atom);

    new_atom.set_position(O17);      // #29
    new_atom.set_type("O");
    new_atom.set_type_number(3);
    kaolinite.push_back(new_atom);

    new_atom.set_position(O18);      // #30
    new_atom.set_type("O");
    new_atom.set_type_number(3);
    kaolinite.push_back(new_atom);

    new_atom.set_position(H5);        // #31
    new_atom.set_type("H");
    new_atom.set_type_number(4);
    kaolinite.push_back(new_atom);

    new_atom.set_position(H6);        // #32
    new_atom.set_type("H");
    new_atom.set_type_number(4);
    kaolinite.push_back(new_atom);

    new_atom.set_position(H7);        // #33
    new_atom.set_type("H");
    new_atom.set_type_number(4);
    kaolinite.push_back(new_atom);

    new_atom.set_position(H8);        // #34
    new_atom.set_type("H");
    new_atom.set_type_number(4);
    kaolinite.push_back(new_atom);


    // Al = 1; Si = 2, O = 3; H = 4.
    masses.push_back(mass_aluminium);
    masses.push_back(mass_silicon);
    masses.push_back(mass_oxygen);
    masses.push_back(mass_hydrogen);

    Initialize(kaolinite, n_unit_cells_x, n_unit_cells_y, n_unit_cells_z, unit_cell_length_x, unit_cell_length_y, unit_cell_length_z);
    vector <Bond> bonds;

    if (output == "Ovito"){
        Write_Initial_State_Ovito(kaolinite, "kaolinite_ovito.txt", n_atom_types);
    }
    else {
        Write_Initial_State_LAMMPS(kaolinite, bonds, "kaolinite_LAMMPS.dat", n_atom_types);
    }
}

void MoleculeSystem::Setup_Portlandite(int n_unit_cells_x, int n_unit_cells_y, int n_unit_cells_z, string output){
    /* Portlandite is a oxide mineral mainly consisting of Calcium Hydroxide ( Ca(OH)_2 )
     *  A unit cell has according to Cygan et al. in Molecular Models of Hydroxide, Oxyhydroxide and Clay
     *  Phases and the Developent of a General Force Field, a,b,c specified as below. The angles between
     *  the axises is given as 120deg between a and b, and that the c-axis is located normal to the ab-plane.
     *
     *  4 molecules per unit cell.
     */

    // unit cell sizes:
    double a = 3.589;      // a
    double b = 3.589;      // b  // 3.585
    double c = 4.911;      // c  // 4.871 ? webmineral.com

    double x,y,z,alpha,beta,gamma, theta, factor, ca_factor;
    double CaOH,OH,CaCa;

    // plane angles (deg):
    alpha = 90;
    beta = 90;
    gamma = 120;
    theta = 30;

    ca_factor = sin(6*PI/18);
    factor = tan(theta*PI/180);

    // bond distances:
    OH = 0.942;            // observed bond distance O-H
    CaOH = 2.369;          // observed bond distance Ca-OH
    //z = 0.93;
    //z = sqrt(pow(CaOH,2.0) - pow(a*tan(3*PI/18),2.0)); // z = 1.14826
    z = 1.148268;
    x = a/2.0;
    y = (b/2.0)*tan(3*PI/18.0);  // lenght scale

    vector <double> Ca1 = {0,0,0};                                  // # 1   v
    vector <double> Ca2 = {a,0,0};                                  // # 2   v
    vector <double> Ca3 = {a*cos(6*PI/18),a*sin(6*PI/18),0};        // # 3   v
    vector <double> Ca4 = {a*(1+cos(6*PI/18)),a*sin(6*PI/18),0};    // # 4   v

    vector <double> O1 = {-a/2.0,-(b/2.0)*tan(3*PI/18.0),z};        // # 5   v
    vector <double> O2 = {0,b*tan(3*PI/18.0),z};                    // # 6   v
    vector <double> O3 = {a/2.0,(b/2.0)*tan(3*PI/18.0),-z};         // # 7   v
    vector <double> O4 = {a/2.0,-y,z};                              // # 8   v
    vector <double> O5 = {a,2*b*tan(3*PI/18.0),-z};                 // # 9   v
    vector <double> O6 = {a,b*tan(3*PI/18.0),z};                    // # 10  v
    vector <double> O7 = {3*a/2.0, (b/2.0)*tan(3*PI/18.0),-z};      // # 11  v
    vector <double> O8 = {2*a,2*b*tan(3*PI/18.0),-z};               // # 12  v

    vector <double> H1 = O1;             // # 13
    H1[2] = O1[2] + OH;
    vector <double> H2 = O2;             // # 14
    H2[2] = O2[2] + OH;
    vector <double> H3 = O3;             // # 15
    H3[2] = O3[2] - OH;
    vector <double> H4 = O4;             // # 16
    H4[2] = O4[2] + OH;
    vector <double> H5 = O5;             // # 17
    H5[2] = O5[2] - OH;
    vector <double> H6 = O6;             // # 18
    H6[2] = O6[2] + OH;
    vector <double> H7 = O7;             // # 19
    H7[2] = O7[2] - OH;
    vector <double> H8 = O8;             // # 20
    H8[2] = O8[2] - OH;

    xlo = -0.5*a;
    ylo = 0.0;
    zlo = 0.0;
    xhi = 2*a;
    yhi = a;
    zhi = 0.0;


    vector <double> v = {0,0,0}; // initial velocity
    vector <double> f = {0,0,0}; // initial force on atom
    double u = 0.0;              // initial potential
    int n_atom_types = 3;        // number of atom types

    vector <Atom> Portlandite;
    vector <Bond> bonds;

    // molecule number 1
    Atom new_atom(Ca1,v,f,u);           // #1 Ca
    new_atom.set_type("Ca");
    new_atom.set_type_number(1);
    new_atom.set_part_of_molecule(1);
    Portlandite.push_back(new_atom);

    new_atom.set_position(O1);          // #5 O
    new_atom.set_type("O");
    new_atom.set_type_number(2);
    new_atom.set_part_of_molecule(1);
    Portlandite.push_back(new_atom);

    new_atom.set_position(O4);          // #8 O
    new_atom.set_type("O");
    new_atom.set_type_number(2);
    new_atom.set_part_of_molecule(1);
    Portlandite.push_back(new_atom);

    new_atom.set_position(H1);          // #13 H
    new_atom.set_type("H");
    new_atom.set_type_number(3);
    new_atom.set_part_of_molecule(1);
    Portlandite.push_back(new_atom);

    new_atom.set_position(H4);          // #16 H
    new_atom.set_type("H");
    new_atom.set_type_number(3);
    new_atom.set_part_of_molecule(1);
    Portlandite.push_back(new_atom);

    Bond bond;
    bond.set_new_bond(1,2,4);
    bonds.push_back(bond);
    bond.set_new_bond(1,3,5);
    bonds.push_back(bond);
    //--------end molecule #1 -------------------
    // molecule number 2
    new_atom.set_position(Ca2);         // #2 Ca
    new_atom.set_type("Ca");
    new_atom.set_type_number(1);
    new_atom.set_part_of_molecule(2);
    Portlandite.push_back(new_atom);

    new_atom.set_position(O2);          // #6 O
    new_atom.set_type("O");
    new_atom.set_type_number(2);
    new_atom.set_part_of_molecule(2);
    Portlandite.push_back(new_atom);

    new_atom.set_position(O6);          // #10 O
    new_atom.set_type("O");
    new_atom.set_type_number(2);
    new_atom.set_part_of_molecule(2);
    Portlandite.push_back(new_atom);

    new_atom.set_position(H2);          // #14 H
    new_atom.set_type("H");
    new_atom.set_type_number(3);
    new_atom.set_part_of_molecule(2);
    Portlandite.push_back(new_atom);

    new_atom.set_position(H6);          // #18 H
    new_atom.set_type("H");
    new_atom.set_type_number(3);
    new_atom.set_part_of_molecule(2);
    Portlandite.push_back(new_atom);

    bond.set_new_bond(1,7,9);
    bonds.push_back(bond);
    bond.set_new_bond(1,8,10);
    bonds.push_back(bond);
    //--------end molecule #2 -------------------
    // molecule number 3
    new_atom.set_position(Ca3);         // #3 Ca
    new_atom.set_type("Ca");
    new_atom.set_type_number(1);
    new_atom.set_part_of_molecule(3);
    Portlandite.push_back(new_atom);

    new_atom.set_position(O3);          // #7 O
    new_atom.set_type("O");
    new_atom.set_type_number(2);
    new_atom.set_part_of_molecule(3);
    Portlandite.push_back(new_atom);

    new_atom.set_position(O7);          // #11 O
    new_atom.set_type("O");
    new_atom.set_type_number(2);
    new_atom.set_part_of_molecule(3);
    Portlandite.push_back(new_atom);

    new_atom.set_position(H3);          // #15 H
    new_atom.set_type("H");
    new_atom.set_type_number(3);
    new_atom.set_part_of_molecule(3);
    Portlandite.push_back(new_atom);

    new_atom.set_position(H7);          // #19 H
    new_atom.set_type("H");
    new_atom.set_type_number(3);
    new_atom.set_part_of_molecule(3);
    Portlandite.push_back(new_atom);

    bond.set_new_bond(1,12,14);
    bonds.push_back(bond);
    bond.set_new_bond(1,13,15);
    bonds.push_back(bond);
    //--------end molecule #4 -------------------
    // molecule number 4
    new_atom.set_position(Ca4);         // #4 Ca
    new_atom.set_type("Ca");
    new_atom.set_type_number(1);
    new_atom.set_part_of_molecule(4);
    Portlandite.push_back(new_atom);

    new_atom.set_position(O5);          // #9 O
    new_atom.set_type("O");
    new_atom.set_type_number(2);
    new_atom.set_part_of_molecule(4);
    Portlandite.push_back(new_atom);

    new_atom.set_position(O8);          // #12 O
    new_atom.set_type("O");
    new_atom.set_type_number(2);
    new_atom.set_part_of_molecule(4);
    Portlandite.push_back(new_atom);

    new_atom.set_position(H5);          // #17 H
    new_atom.set_type("H");
    new_atom.set_type_number(3);
    new_atom.set_part_of_molecule(4);
    Portlandite.push_back(new_atom);

    new_atom.set_position(H8);          // #20 H
    new_atom.set_type("H");
    new_atom.set_type_number(3);
    new_atom.set_part_of_molecule(4);
    Portlandite.push_back(new_atom);

    bond.set_new_bond(1,17,19);
    bonds.push_back(bond);
    bond.set_new_bond(1,18,20);
    bonds.push_back(bond);
    //------ end molecule #4 ---------------------

    // Ca = 1; O = 2; H = 3;
    masses.push_back(mass_calcium);   // 1
    masses.push_back(mass_oxygen);    // 2
    masses.push_back(mass_hydrogen);  // 3

    int molecules_per_unit_cell = 4;
    int bonds_per_unit_cell = 8;
    double unit_cell_length_x, unit_cell_length_y, unit_cell_length_z;
    unit_cell_length_x = 2*a;
    unit_cell_length_y = 2*b*sin(6*PI/18);
    unit_cell_length_z = c;

    Initialize_2pointO(Portlandite,bonds,n_unit_cells_x,n_unit_cells_y, n_unit_cells_z, unit_cell_length_x,unit_cell_length_y, unit_cell_length_z,molecules_per_unit_cell,bonds_per_unit_cell);

    if (output == "Ovito"){
        Write_Initial_State_Ovito(Portlandite, "portlandite_ovito.txt", n_atom_types);
    }
    else {
        Write_Initial_State_LAMMPS(Portlandite, bonds, "portlandite_LAMMPS.dat", n_atom_types);

    }

}

// Initialize 2.0 is working now, and updates molecule number as it should. Just set the molecule number to the atoms in the
// function where the unit cell i created!
void MoleculeSystem::Initialize_2pointO(vector <Atom> &atoms, vector <Bond> &bonds, int &lx, int &ly, int &lz, double &xi, double &yi, double &zi, int molecules_per_unit_cell, int bonds_per_unit_cell){

    vector <double> ri (3,0.0);
    vector <double> r (3,0.0);
    vector <int> ANB (2,0);
    string atom_type;
    int atom_type_nr;
    int molecule_nr;
    int unit_cell_nr = 0;
    int N_unit_cell_atoms = atoms.size();  // number of atoms in one unit cell.
    int N_bonds = bonds.size();
    int bond_nr = 0;
    int bond_type = 0;
    Atom new_atom = atoms[0];
    Bond new_bond;
    vector <Atom> AtomContainer;
    vector <Bond> BondContainer;

    // unit cell loop:
    for (int i=0; i<lx; i++){          // unit cell x
        for (int j=0; j<ly; j++){      // unit cell y
            for (int k=0; k<lz; k++){  // unit cell z

                ri[0] = i*xi; ri[1] = j*yi; ri[2] = k*zi;

                if (ri[0] < xlo){ xlo = ri[0];}
                if (ri[0] > xhi){ xhi = ri[0];}
                if (ri[1] < ylo){ ylo = ri[1];}
                if (ri[1] > yhi){ yhi = ri[1];}
                if (ri[2] < zlo){ zlo = ri[2];}
                if (ri[2] > zhi){ zhi = ri[2];}

                for (int atom=0; atom<N_unit_cell_atoms ;atom++){
                    r = atoms[atom].get_position();
                    atom_type = atoms[atom].get_type();
                    atom_type_nr = atoms[atom].get_type_number();

                    molecule_nr =  molecules_per_unit_cell*unit_cell_nr + atoms[atom].get_molecule_number();

                    /*
                    cout << "------------------------------------------" << endl;
                    cout << "molecule nr     : " << molecule_nr <<  " atom nr: " << 1 + atom + unit_cell_nr*N_unit_cell_atoms << endl;
                    cout << "unit_cell_nr    : " << unit_cell_nr << " i=" << i <<" j=" << j << " k=" << k << endl;
                    cout << "atom molecule nr: " << atoms[atom].get_molecule_number() << endl;
                    */

                    for (int cor=0; cor<3; cor++){ r[cor] = r[cor] + ri[cor];}

                    atom_type = atoms[atom].get_type();
                    new_atom.set_position(r);
                    new_atom.set_type(atom_type);
                    new_atom.set_type_number(atom_type_nr);
                    new_atom.set_part_of_molecule(molecule_nr);
                    AtomContainer.push_back(new_atom);

                }
                for (int bond=0; bond<N_bonds;bond++){
                    ANB = bonds[bond].get_atoms();
                    bond_nr = bonds_per_unit_cell*unit_cell_nr; // ANB[] + bond_nr = atom number
                    bond_type = bonds[bond].get_bond_type();
                    new_bond.set_new_bond(bond_type,ANB[0]+bond_nr,ANB[1]+bond_nr);
                    BondContainer.push_back(new_bond);

                }
            unit_cell_nr = unit_cell_nr + 1;
            }
        }
    }
    //cout << AtomContainer.size() << endl;
    atoms = AtomContainer;
    bonds = BondContainer;
}

void MoleculeSystem::Initialize(vector <Atom> &atoms, int &lx, int &ly, int &lz, double &xi, double &yi, double &zi){

    vector <double> ri (3,0.0);
    vector <double> r (3,0.0);
    string atom_type;
    int atom_type_nr;
    int molecule_nr;                                  // we assume that one unit cell contains one molecule. This is not necessarly correct!!!!
    int N_unit_cell_atoms = atoms.size();
    Atom new_atom = atoms[0];

    vector <Atom> AtomContainer;

    // unit cell loop:
    for (int i=0; i<lx; i++){          // unit cell x
        for (int j=0; j<ly; j++){      // unit cell y
            for (int k=0; k<lz; k++){  // unit cell z
                molecule_nr = i + j + k + 1;

                ri[0] = i*xi; ri[1] = j*yi; ri[2] = k*zi;

                if (ri[0] < xlo){ xlo = ri[0];}
                if (ri[0] > xhi){ xhi = ri[0];}
                if (ri[1] < ylo){ ylo = ri[1];}
                if (ri[1] > yhi){ yhi = ri[1];}
                if (ri[2] < zlo){ zlo = ri[2];}
                if (ri[2] > zhi){ zhi = ri[2];}

                for (int atom=0; atom<N_unit_cell_atoms ;atom++){
                    r = atoms[atom].get_position();
                    atom_type = atoms[atom].get_type();
                    atom_type_nr = atoms[atom].get_type_number();

                    for (int cor=0; cor<3; cor++){
                        r[cor] = r[cor] + ri[cor];
                    }
                    atom_type = atoms[atom].get_type();
                    new_atom.set_position(r);
                    new_atom.set_type(atom_type);
                    new_atom.set_type_number(atom_type_nr);
                    new_atom.set_part_of_molecule(molecule_nr);
                    AtomContainer.push_back(new_atom);
                }
            }
        }
    }
    //cout << AtomContainer.size() << endl;
    atoms = AtomContainer;
}

void MoleculeSystem::Write_Initial_State_Ovito(vector <Atom> &atoms, string filename, int &n_atom_types){

    vector <double> r (3,0.0);

    ofstream myfile;
    myfile.open(filename);
    myfile << atoms.size() << endl;          // number of atoms
    myfile << 0 << endl;                      // timestep 0

    for (uint i=0; i<atoms.size(); i++){
        r = atoms[i].get_position();
        myfile << atoms[i].get_type() << " " << r[0] << " " << r[1] << " " << r[2] << endl;
    }
    myfile.close();

    cout << "--------------------------------------------------------" << endl;
    cout << "****   Initial state file for Ovito is generated    ****" << endl;
    cout << "Filename            : " << filename << endl;
    cout << "Number of atoms     : " << atoms.size() << endl;
    cout << "Number of atom types: " << n_atom_types << endl;

}

void MoleculeSystem::Write_Initial_State_LAMMPS(vector<Atom> &atoms, string &filename, int &n_atom_types){

    vector <double> r (3,0.0);

    int bond_types = 1; // for now...

    ofstream myfile;
    myfile.open(filename);
    myfile << "# LAMMPS data set file. Genereated by Goeran Brekke Svaland" << endl;
    myfile << "           " << atoms.size() << " atoms" << endl;                    // number of atoms
    myfile << "           " << n_atom_types << " atom types" << endl;        // number of atom types
    myfile << "           " << bonds.size() << " bonds" << endl;
    myfile << "           " << bond_types   << " bond types" << endl;

    /*
    myfile << number of bonds << endl;
    myfile << angles << endl;
    myfile << dihedrals << endl;
    myfile << impropers << endl;
    myfile << atom types << endl;
    myfile << bond types << endl;
    myfile << angle types << endl;
    myfile << dihedral types << endl;
    myfile << improper types << endl;

    */

    myfile << "    " <<  xlo << "    " << xhi << "    xlo xhi" << endl;
    myfile << "    " <<  ylo << "    " << yhi << "    ylo yhi" << endl;
    myfile << "    " <<  zlo << "    " << zhi << "    zlo zhi" << endl;

    myfile << endl << "Masses" << endl;
    //myfile << "# atom-type  mass(au)" << endl << endl;

    for (int j=0; j<n_atom_types; j++){
        myfile << "       " << j+1 << "        " << masses[j] << endl;
    }

    myfile << endl <<  "Atoms" << endl;
    //myfile << "# atom-ID  molecule-ID  atom-type   x      y      z" << endl << endl;

    for (uint i=0; i<atoms.size(); i++){
        r = atoms[i].get_position();
        myfile << "       "  << i+1 << "       " << atoms[i].get_molecule_number() << "       " << atoms[i].get_type_number() <<  "       " << r[0] << "       " << r[1] << "       " << r[2] << endl;
    }


    // bonds :-)
    myfile << "Bonds" << endl << endl;
    myfile << "# Bond-nr   bond-type   atom1   atom2" << endl;
    vector <int> ANB (2,0);

    for (int k=0;k<bonds.size();k++){
        ANB = bonds[k].get_atoms();
        myfile << "       "  << k+1 << "       " << bonds[k].get_bond_type() << "       " << ANB[0] << "       " << ANB[1] << endl;
    }


    myfile.close();

    cout << "--------------------------------------------------------" << endl;
    cout << "****   Initial state file for Ovito is generated    ****" << endl;
    cout << "Filename            : " << filename << endl;
    cout << "Number of atoms     : " << atoms.size() << endl;
    cout << "Number of atom types: " << n_atom_types << endl;
    cout << xlo << " " << xhi << "     xlo xhi" << endl;
    cout << ylo << " " << yhi << "     ylo yhi" << endl;
    cout << zlo << " " << zhi << "     zlo zhi" << endl;
    cout << "--------------------------------------------------------" << endl;

} // end;


