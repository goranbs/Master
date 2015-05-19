#ifndef ATOM_H
#define ATOM_H

#include <vector>
#include <string>
#include <iostream>

using namespace std;

class atom
{
public:
    atom();
    void set_atom_initial_state(int &atomID, int &atomTYPE, vector <double> &position, vector <double> &velocity, vector <double> &forces);
    void set_index_number(int &atomID);
    void set_type_number(int &atomTYPE);
    void set_type_element(char &atomTYPE);
    void set_charge(double &q);
    void set_part_of_molecule(int &moleculeNR);
    void set_systemsize(double &xlo, double &xhi, double &ylo, double &yhi, double &zlo, double &zhi);
    void set_box_number(int &boxnumber);
    void set_past_position(vector <double> &past_position);
    vector <double> get_past_position();
    vector <double> get_squared_traveled_distance();
    vector <double> get_crossings();
    vector <double> get_position();
    vector <double> get_initial_position();
    vector <double> get_velocity();
    vector <double> get_force();
    vector <double> get_systemsize();
    int get_box_number();
    int get_atom_type_number();
    char get_atom_type();
    double get_charge();
    int get_molecule_number();
    int get_atomID();
    void update_systemsize(double &xlo, double &xhi, double &ylo, double &yhi, double &zlo, double &zhi);
    void cross_periodic_boundary(int ni, int nj, int nk);
    void update_position(vector <double> &r_);
    void update_velocity(vector <double> &v_);
    void update_force(vector <double> &f_);
    void subtract_force(vector<double> &f_);
    void add_force(vector <double> &f_);
    void clear_force();

private:
    vector <double> n_crossings;
    int INDEX;
    int TYPENR;
    vector <double> POS;
    vector <double> VEL;
    vector <double> FOR;
    vector <double> IPOS;
    vector <double> SQUARED_TRAVELED_DISTANCE;
    vector <double> SystemSize;
    vector <double> PastPos;
    char ATOM_TYPE;
    int ATOM_CHARGE;
    int MOLECULE_NUMBER;
    int BOXnr;

};

#endif // ATOM_H
