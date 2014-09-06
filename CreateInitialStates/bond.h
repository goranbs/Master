#ifndef BOND_H
#define BOND_H

#include "atom.h"

using namespace std;

class Bond
{
public:
    Bond();
    void set_bond(int &bondtype, Atom &atom1, Atom &atom2);
    void set_new_atom2(Atom &atom2);
    void set_atom_index_number1(int &index_number);
    void set_atom_index_number2(int &index_number);

    int get_bondtype();
    int get_type_atom1();
    int get_type_atom2();
    int get_atom_index_number1();
    int get_atom_index_number2();

private:
    Atom Atom1;
    Atom Atom2;
    int BOND_TYPE;
    int TYPE1;
    int TYPE2;
    int INDEX1;
    int INDEX2;
};

#endif // BOND_H
