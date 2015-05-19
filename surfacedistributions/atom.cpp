#include "atom.h"

atom::atom(){
}

/*****************************************************
 *                   Class atom
 */

void atom::set_atom_initial_state(int &atomID, int &atomTYPE, vector <double> &position, vector <double> &velocity, vector <double> &forces){
    // construct the Atom object that holds the r, v, f, n_crossings and initial position r0.

    n_crossings = vector <double> (3);  // Crossings out of the system x,y,z
    INDEX = atomID;                            // atomID
    TYPENR = atomTYPE;                         // atom type
    POS = position;                            // position of atom
    VEL = velocity;                            // velocity of atom
    FOR = forces;                              // forces on atom
    IPOS = position;                           // initial position
    vector <double> sqd (3,0.0);
    vector <double> szize (6,0.0);
    SQUARED_TRAVELED_DISTANCE = sqd;
    SystemSize = szize;
    PastPos = sqd;

}

void atom::set_past_position(vector <double> &past_position){
    for (int i=0;i<3;i++){ PastPos[i] = past_position[i];}
}

void atom::set_box_number(int &boxnumber){
    BOXnr = boxnumber;
}

void atom::set_index_number(int &atomID){
    INDEX = atomID;
}

void atom::set_type_number(int &atomTYPE){     // integer value atom type
    TYPENR = atomTYPE;
}

void atom::set_type_element(char &atomTYPE){   // character atom type
    ATOM_TYPE = atomTYPE;
}

void atom::set_charge(double &q){
    ATOM_CHARGE = q;
}

void atom::set_part_of_molecule(int &moleculeNR){
    MOLECULE_NUMBER = moleculeNR;
}

void atom::set_systemsize(double &xlo, double &xhi, double &ylo, double &yhi, double &zlo, double &zhi){
    SystemSize[0] = xlo;
    SystemSize[1] = xhi;
    SystemSize[2] = ylo;
    SystemSize[3] = yhi;
    SystemSize[4] = zlo;
    SystemSize[5] = zhi;
}

//-----------------------------------------------------------------------------------------//
// get values:

vector <double> atom::get_past_position(){
    return PastPos;
}

vector <double> atom::get_systemsize(){
    return SystemSize;
}

vector <double> atom::get_squared_traveled_distance(){
    return SQUARED_TRAVELED_DISTANCE;
}

vector <double> atom::get_crossings(){
    return n_crossings;
}

vector <double> atom::get_position(){
    return POS;
}

vector <double> atom::get_initial_position(){
    return IPOS;
}

vector <double> atom::get_velocity(){
    return VEL;
}

vector <double> atom::get_force(){
    return FOR;
}

int atom::get_atom_type_number(){
    return TYPENR;
}

char atom::get_atom_type(){
    return ATOM_TYPE;
}

double atom::get_charge(){
    return ATOM_CHARGE;
}

int atom::get_molecule_number(){
    return MOLECULE_NUMBER;
}

int atom::get_atomID(){
    return INDEX;
}

int atom::get_box_number(){
    return BOXnr;
}

//-----------------------------------------------------------------------------------------//
//  update values :

void atom::cross_periodic_boundary(int ni, int nj, int nk){
    /* ni = 1, if the atom has crossed the boundary in positive x-direction, ni = -1 for negative direction
     * and ni is zero if it has not crossed.
     * Same applies for nj for y-direction, and nk for z. */
    n_crossings[0] += ni;
    n_crossings[1] += nj;
    n_crossings[2] += nk;
}

void atom::update_position(vector <double> &r_){          // R. reset postion
    for (int i = 0; i < 3; ++i) {
        POS[i] = r_[i];
        SQUARED_TRAVELED_DISTANCE[i] = (r_[i] - IPOS[i])*(r_[i] - IPOS[i]);
    }
}
void atom::update_velocity(vector <double> &v_){          // V. reset velocity
    for (int i = 0; i < 3; ++i) {
        VEL[i] = v_[i];
    }
}
void atom::update_force(vector <double> &f_){             // F. reset force
    for (int i = 0; i < 3; ++i) {
        FOR[i] = f_[i];
    }
}
void atom::subtract_force(vector<double> &f_){            // F. subtract force from existing
    for (int i = 0; i < 3; ++i) {
        FOR[i] -= f_[i];
    }
}
void atom::add_force(vector <double> &f_){                // F. add force to existing
    for (int i = 0; i < 3; ++i) {
        FOR[i] += f_[i];
    }
}
void atom::clear_force(){                                       // F. clear force vector. Force set to zero
    for (int i = 0; i < 3; ++i) {
        FOR[i]= 0;
    }
}

void atom::update_systemsize(double &xlo, double &xhi, double &ylo, double &yhi, double &zlo, double &zhi){
    SystemSize[0] = xlo;
    SystemSize[1] = xhi;
    SystemSize[2] = ylo;
    SystemSize[3] = yhi;
    SystemSize[4] = zlo;
    SystemSize[5] = zhi;
}

