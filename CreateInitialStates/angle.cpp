#include "angle.h"

Angle::Angle()
{
}

void Angle::set_angle(int &angletype, Atom &atom1, Atom &atom2, Atom &atom3){
    ANGLE_TYPE = angletype;
    ATOM_TYPE1 = atom1.get_type();
    ATOM_TYPE2 = atom2.get_type();
    ATOM_TYPE3 = atom3.get_type();
    ATOMTYPE_NUMBER1 = atom1.get_type_number();
    ATOMTYPE_NUMBER2 = atom2.get_type_number();
    ATOMTYPE_NUMBER3 = atom3.get_type_number();
    INDEX1 = atom1.get_index_number();
    INDEX2 = atom2.get_index_number();
    INDEX3 = atom3.get_index_number();
}

void Angle::set_atom_index_numbers(int &index1, int &index2, int &index3){
    INDEX1 = index1;
    INDEX2 = index2;
    INDEX3 = index3;

}

void Angle::set_angletype(int &angletype){
    ANGLE_TYPE = angletype;
}

void Angle::reset_atom1(Atom &atom1){
    ATOM_TYPE1 = atom1.get_type();
    ATOMTYPE_NUMBER1 = atom1.get_type_number();
    INDEX1 = atom1.get_index_number();
}

void Angle::reset_atom2(Atom &atom2){
    ATOM_TYPE2 = atom2.get_type();
    ATOMTYPE_NUMBER2 = atom2.get_type_number();
    INDEX2 = atom2.get_index_number();
}

void Angle::reset_atom3(Atom &atom3){
    ATOM_TYPE3 = atom3.get_type();
    ATOMTYPE_NUMBER3 = atom3.get_type_number();
    INDEX3 = atom3.get_index_number();
}

string Angle::get_atomtype1(){       // atom type (string)
    return ATOM_TYPE1;
}

string Angle::get_atomtype2(){
    return ATOM_TYPE2;
}

string Angle::get_atomtype3(){
    return ATOM_TYPE3;
}

int Angle::get_atomtype_number1(){   // atom type (number)
    return ATOMTYPE_NUMBER1;
}

int Angle::get_atomtype_number2(){
    return ATOMTYPE_NUMBER2;
}

int Angle::get_atomtype_number3(){
    return ATOMTYPE_NUMBER3;
}

int Angle::get_atom1_index_number(){  // atom index number
    return INDEX1;
}

int Angle::get_atom2_index_number(){
    return INDEX2;
}

int Angle::get_atom3_index_number(){
    return INDEX3;
}

int Angle::get_angletype(){           // angle type
    return ANGLE_TYPE;
}
