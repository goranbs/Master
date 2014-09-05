#ifndef MOLECULESYSTEM_H
#define MOLECULESYSTEM_H

#include <vector>
#include "atom.h"
#include <string>
#include "bond.h"

using namespace std;


class MoleculeSystem
{
public:
    MoleculeSystem();
    void Cout_Info();
    void Setup_Argon();
    void Setup_CalciumCarbonite();
    void Setup_CarbonDioxide();
    void Setup_Water();
    void Setup_Kaolinite(int n_unit_cells_x, int n_unit_cells_y, int n_unit_cells_z, string output);
    void Setup_Calcium(int n_unit_cells_x, int n_unit_cells_y, int n_unit_cells_z, string output);
    void Setup_Portlandite(int n_unit_cells_x, int n_unit_cells_y, int n_unit_cells_z, string output);
    void Setup_calciumhydroxide(int n_unit_cells_x, int n_unit_cells_y, int n_unit_cells_z, string output);



private:
    void Initialize(vector <Atom> &atoms, int &lx, int &ly, int &lz, double &xi, double &yi, double &zi);
    void Initialize_2pointO(vector <Atom> &atoms, int &lx, int &ly, int &lz, double &xi, double &yi, double &zi);
    void Initialize_3pointO(vector <Atom> &atoms, vector <Bond> &bonds, int &lx, int &ly, int &lz, double &xi, double &yi, double &zi);
    void Write_Initial_State_Ovito(vector <Atom> &atoms, string filename, int &n_atom_types);
    void Write_Initial_State_LAMMPS(vector <Atom> &atoms, string filename, int &n_atom_types);
    void Write_Initial_State_LAMMPS_2pointO(vector<Atom> &atoms, vector <Bond> &bonds, string filename, int &n_atom_types);
    double xlo, xhi, ylo, yhi, zlo, zhi;  // system size.
    vector <double> masses;

    // masses (au):
    double mass_hydrogen = 1.00794;    // H
    double mass_oxygen = 15.9994;      // O
    double mass_aluminium = 26.9815;   // Al
    double mass_silicon = 28.065;      // Si
    double mass_calcium = 40.078;      // Ca



};

#endif // MOLECULESYSTEM_H
