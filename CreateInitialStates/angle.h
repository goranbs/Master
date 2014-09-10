#ifndef ANGLE_H
#define ANGLE_H

#include "atom.h"
#include <string>

using namespace std;

class Angle
{
public:
    Angle();

    void set_angle(int &angletype, Atom &atom1, Atom &atom2, Atom &atom3);
    void reset_atom1(Atom &atom1);
    void reset_atom2(Atom &atom2);
    void reset_atom3(Atom &atom3);
    void set_angletype(int &angletype);
    void set_atom_index_numbers(int &index1, int &index2, int &index3);

    string get_atomtype1();
    string get_atomtype2();
    string get_atomtype3();

    int get_atomtype_number1();
    int get_atomtype_number2();
    int get_atomtype_number3();

    int get_atom1_index_number();
    int get_atom2_index_number();
    int get_atom3_index_number();

    int get_angletype();

private:

    string ATOM_TYPE1;
    string ATOM_TYPE2;
    string ATOM_TYPE3;
    int ATOMTYPE_NUMBER1;
    int ATOMTYPE_NUMBER2;
    int ATOMTYPE_NUMBER3;
    int ANGLE_TYPE;
    int INDEX1;
    int INDEX2;
    int INDEX3;

};

#endif // ANGLE_H
