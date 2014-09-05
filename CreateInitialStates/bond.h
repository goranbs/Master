#ifndef BOND_H
#define BOND_H

#include <atom.h>
#include <moleculesystem.h>

using namespace std;

class Bond
{
public:
    Bond();
    void set_bond(int &bondtype, Atom &atom1, Atom &atom2);
    void set_new_atom2(Atom &atom2);

    int get_bondtype();
    int get_type_atom1();
    int get_type_atom2();

private:
    Atom Atom1;
    Atom Atom2;
    int BOND_TYPE;
    int TYPE1;
    int TYPE2;
};

#endif // BOND_H
