#ifndef BOND_H
#define BOND_H

#include <atom.h>

using namespace std;

class Bond{

public:

    Bond();
    void set_new_bond(int bondtype, int atom1, int atom2);
    void set_atom_type_number(int type_atom1, int type_atom2);
    int get_bond_type();
    int get_atom_type1();
    int get_atom_type2();
    vector <int> get_atoms();

private:
    int Atom1;
    int Atom2;
    int TYPE1;
    int TYPE2;
    int BOND_TYPE;
};

#endif // BOND_H
