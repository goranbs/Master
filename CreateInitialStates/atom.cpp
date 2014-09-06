#include "atom.h"

/*****************************************************
 *                   Class Atom
 * contains information about atom i such as
 * - Atomnumber/index    i          Atom.index(i) = i
 * - Position            r(x,y,x)   Atom.position(i) = r(x,y,z) for atom i
 * - Velocity            v(x,y,z)   Atom.velocity(i) = v(x,y,z) for atom i
 * - Initial position    ri(x,y,z)  Atom.initial(i)  = ri(x,y,z) for atom i
 * - Travel              r^2(t)     Atom.travel(i)   = r(t) - r_initial
 */


Atom::Atom(){

}

void Atom::set_initial_state(vector <double> &r_, vector <double> &v_, vector <double> &f_, double &u_){
    // construct the Atom object that holds the r, v, f, u, n_crossings and initial position r0.
    //dist = vector <double> (3);
    n_crossings = vector <double> (3);
    //N = 0;
    u = u_;
    r = r_;
    v = v_;
    f = f_;
    r0 = r_;

    setIs_matrix(false);  // default
}

void Atom::set_all(const vector<double> &r_, const vector<double> &v_, const vector<double> &f_, const double &u_, int &atom_type_nr, double &q, int &index_number){
    u = u_;
    r = r_;
    v = v_;
    f = f_;
    r0 = r_;
    setIs_matrix(false);  // default
    //ATOM_TYPE = atom_type;
    ATOM_TYPE_NUMBER = atom_type_nr;
    ATOM_CHARGE = q;
    INDEX_NUMBER = index_number;

}

void Atom::set_position(vector <double> &r_){
    for (int cor=0; cor <3; cor++){
        r[cor] = r_[cor];
    }
}

void Atom::set_type(string atom_type){

    ATOM_TYPE = atom_type;
/*
    string val = "NO";

    for (int i=0; i<atomlist.size(); i++){
        if (atomlist[i] == atom_type){
            val == "OK";
        }
    }

    if (val == "OK"){
        ATOM_TYPE = atom_type;
    }
    else {
        ATOM_TYPE = "NONE";
        cout << "ATOM TYPE is not in vocabulary" << endl;
    }
    */
}

void Atom::set_index_number(int &index_number){
    INDEX_NUMBER = index_number;
}

void Atom::set_type_number(int atom_type_nr){
    ATOM_TYPE_NUMBER = atom_type_nr;
}

void Atom::set_charge(double q){
    ATOM_CHARGE = q;
}

void Atom::set_part_of_molecule(int molecule_nr){
    MOLECULE_NUMBER = molecule_nr;
}

void Atom::cross_boundary(int i, int j, int k){
    /* i = 1, if the atom has crossed the boundary in positive x-direction, i = -1 for negative direction
     * and i is zero if it has not crossed.
     * Same applies for j for y-direction, and k for z. */
    n_crossings[0] += i;
    n_crossings[1] += j;
    n_crossings[2] += k;
}

bool Atom::getIs_matrix() const
{
    return is_matrix;
}

void Atom::setIs_matrix(bool value){
    if (value == true) {
    update_velocity({0,0,0});
    }
    is_matrix = value;
}

vector <double> Atom::get_position(){
    return r;
}

string Atom::get_type(){
    return ATOM_TYPE;
}

int Atom::get_type_number(){
    return ATOM_TYPE_NUMBER;
}

double Atom::get_charge(){
    return ATOM_CHARGE;
}

int Atom::get_molecule_number(){
    return MOLECULE_NUMBER;
}

int Atom::get_index_number(){
    return INDEX_NUMBER;
}

void Atom::update_position(const vector <double> &r_){          // R. reset postion
    for (int i = 0; i < 3; ++i) {
        r[i] = r_[i];
    }
}
void Atom::update_velocity(const vector <double> &v_){          // V. reset velocity
    for (int i = 0; i < 3; ++i) {
        v[i] = v_[i];
    }
}
void Atom::update_force(const vector <double> &f_){             // F. reset force
    for (int i = 0; i < 3; ++i) {
        f[i] = f_[i];
    }
}
void Atom::subtract_force(const vector<double> &f_){            // F. subtract force from existing
    for (int i = 0; i < 3; ++i) {
        f[i] -= f_[i];
    }
}

void Atom::add_force(const vector <double> &f_){                // F. add force to existing
    for (int i = 0; i < 3; ++i) {
        f[i] += f_[i];
    }
}
void Atom::clear_force(){                                       // F. clear force vector. Force set to zero
    for (int i = 0; i < 3; ++i) {
        f[i]= 0;
    }
}
void Atom::update_potential(double &u_){                        // P. reset potential
    u = u_;
}
void Atom::add_potential(double &u_){                           // P. add potential to existing
    u += u_;
}

void Atom::clear_potential(){                                   // P. clear potential. Potential set to zero.
    u = 0;
}
