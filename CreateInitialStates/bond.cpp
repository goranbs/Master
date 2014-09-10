#include "bond.h"

Bond::Bond(){

}

void Bond::set_bond(int &bondtype, Atom &atom1, Atom &atom2){
    Atom1 = atom1;
    Atom2 = atom2;
    BOND_TYPE = bondtype;
    TYPE1 = Atom1.get_type_number();
    TYPE2 = Atom2.get_type_number();
    ATOM_TYPE1 = Atom1.get_type();
    ATOM_TYPE2 = Atom2.get_type();
}

void Bond::set_new_atom2(Atom &atom2){
    Atom2 = atom2;
    TYPE2 = Atom2.get_type_number();
}

void Bond::set_atom_index_number1(int index_number){
    INDEX1 = index_number;
}

void Bond::set_atom_index_number2(int index_number){
    INDEX2 = index_number;
}

int Bond::get_bondtype(){
    return BOND_TYPE;
}

int Bond::get_atom1_typenumber(){
    return TYPE1;
}

int Bond::get_atom2_typenumber(){
    return TYPE2;
}

int Bond::get_atom_index_number1(){
    return INDEX1;
}

int Bond::get_atom_index_number2(){
    return INDEX2;
}

string Bond::get_atomtype1(){
    return ATOM_TYPE1;
}

string Bond::get_atomtype2(){
    return ATOM_TYPE2;
}
