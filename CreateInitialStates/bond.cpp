#include "bond.h"

Bond::Bond(){

}

void Bond::set_bond(int &bondtype, Atom &atom1, Atom &atom2){
    Atom1 = atom1;
    Atom2 = atom2;
    BOND_TYPE = bondtype;
    TYPE1 = Atom1.get_type_number();
    TYPE2 = Atom2.get_type_number();
}

void Bond::set_new_atom2(Atom &atom2){
    Atom2 = atom2;
    TYPE2 = Atom2.get_type_number();
}

int Bond::get_bondtype(){
    return BOND_TYPE;
}

int Bond::get_type_atom1(){
    return TYPE1;
}

int Bond::get_type_atom2(){
    return TYPE2;
}
