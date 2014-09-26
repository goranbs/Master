#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include "atom.h"
#include "bond.h"
#include "moleculesystem.h"

using namespace std;

int main(){

    // should be able to choose between generating a LAMMPS file or a ovito file!

    //string output = "Ovito";
    //string output = "LAMMPS";
    string output = "Charge";

    /*
    // setup Kaolinite:
    MoleculeSystem Kaolinite;
    Kaolinite.Setup_Kaolinite(1,1,1, output);

    // setup a calcium surface:
    MoleculeSystem Calcium;
    Calcium.Setup_Calcium(1,10,10,output);
*/

    // setup a Portlandite surface:
    MoleculeSystem Portlandite;
    Portlandite.Setup_Portlandite(1,1,2, output);



    return 0;
}

