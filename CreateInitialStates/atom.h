#ifndef ATOM_H
#define ATOM_H

#include <vector>
#include <string>
#include <iostream>

using namespace std;

class Atom{

public:

    Atom();
    void set_all(const vector <double> &r_, const vector <double> &v_ , const vector <double> &f_, const double &u_, int &atom_type_nr, double &q, int &index_number);
    void set_initial_state(vector <double> &r_, vector <double> &v_, vector <double> &f_, double &u_);
    void update_position(const vector <double> &r_);
    void update_velocity(const vector <double> &v_);
    void update_force(const vector <double> &f_);
    void subtract_force(const vector <double> &f_);
    void add_force(const vector <double> &f_);
    void update_potential(double &u_);
    void add_potential(double &u_);
    void clear_force();
    void clear_potential();
    inline vector <double> &position();
    inline vector <double> &velocity();
    inline vector <double> &force();
    inline double potential();
//    inline const double potential() const;

    inline const vector<double> &return_initial_position();
    void cross_boundary(int i, int j , int k);
    inline const vector<double> &return_n_crossings() const;

    void set_index_number(int &index_number);
    void set_part_of_molecule(int molecule_nr);
    void set_type(string atom_type);
    void set_type_number(int atom_type_nr);
    void set_charge(double q);
    void setIs_matrix(bool value);
    void set_position(vector <double> &r_);

    int get_index_number();
    bool getIs_matrix() const;
    string get_type();
    int get_type_number();
    int get_molecule_number();
    double get_charge();
    vector <double> get_position();

private:
    vector <double> r;            // position of atom
    vector <double> v;            // velocity
    vector <double> f;            // total force felt from other particles
    vector <double> r0;           // initial position
    vector <double> n_crossings ; // number of crossings out of the system (+1 pos dir, -1 neg dir.)
    //vector <double> dist;         // distance traveled.
    double u;                     // total potential
    //int N;
    bool is_matrix;
    string ATOM_TYPE;
    int ATOM_TYPE_NUMBER;
    double ATOM_CHARGE;
    int MOLECULE_NUMBER;
    int INDEX_NUMBER;
    vector <string> atomlist = {"H", "He", "Li", "Be", "C", "O", "Al", "Si", "Ca"};

};

inline vector < double > &Atom::position(){               // R. return position
    return r;
}
inline vector < double > &Atom::velocity(){               // V. return velocity
    return v;
}
inline vector < double > &Atom::force(){                  // F. return force
    return f;
}
inline double Atom::potential() {                        // P. return potential
    return u;
}
inline const vector<double> &Atom::return_initial_position(){   // R0. return initial position
    return r0;
}
inline const vector < double > &Atom::return_n_crossings() const{           // number of crossings out of system oundaries
    return n_crossings;
}


#endif // ATOM_H
