#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <atom.h>
#include <moleculesystem.h>

using namespace std;

int main(){

    // setup Kaolinite:
    MoleculeSystem Kaolinite;
    Kaolinite.Setup_Kaolinite(1,2,4);

    // setup a calcium surface:
    MoleculeSystem Calcium;
    Calcium.Setup_Calcium(3,10,10);

    return 0;
}

