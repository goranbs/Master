#ifndef MOLECULESYSTEM_H
#define MOLECULESYSTEM_H

#include <vector>
#include <atom.h>
#include <string>

class MoleculeSystem
{
public:
    MoleculeSystem();
    void Cout_Info();
    void Setup_Argon();
    void Setup_CalciumCarbonite();
    void Setup_CarbonDioxide();
    void Setup_Water();
    void Setup_Kaolinite(int n_unit_cells_x, int n_unit_cells_y, int n_unit_cells_z);
    void Setup_Calcium(int n_unit_cells_x, int n_unit_cells_y, int n_unit_cells_z);


private:
    void Initialize(vector <Atom> &atoms, int &lx, int &ly, int &lz, double &xi, double &yi, double &zi);
    void Write_Initial_State(vector <Atom> &atoms, string filename);

};

#endif // MOLECULESYSTEM_H
