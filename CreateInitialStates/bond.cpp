#include "bond.h"

Bond::Bond(){
    /*
     * Make sure to set an integer value for different bondtypes.
     * For example set bondtype = 1 for bond between O-H and
     * bondtype = 2 for O-O bonds.
    */
}

void Bond::set_new_bond(int bondtype, int atom1, int atom2){
    Atom1 = atom1;
    Atom2 = atom2;
    BOND_TYPE = bondtype;
}

void Bond::set_atom_type_number(int type_atom1, int type_atom2){
    TYPE1 = type_atom1;
    TYPE2 = type_atom2;
}

//---------------------------------------------------------------
int Bond::get_atom_type1(){
    return TYPE1;
}

int Bond::get_atom_type2(){
    return TYPE2;
}

int Bond::get_bond_type(){
    return BOND_TYPE;
}

vector <int> Bond::get_atoms(){
    vector <int> atoms = {Atom1, Atom2};
    return atoms;
}
