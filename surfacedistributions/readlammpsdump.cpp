#include "readlammpsdump.h"

ReadLAMMPSdump::ReadLAMMPSdump(string &path){
    PATH = path;
    ENDCHAR = PATH.back();
}

void ReadLAMMPSdump::read_initialfile(string &filename, vector <atom> &Atoms){
    string fullpath;
    ENDCHAR = PATH.back();
    if (ENDCHAR == "/"){fullpath = PATH + filename;}
    else {
        fullpath = PATH + "/" + filename;
        PATH = PATH + "/";
    }

    ifstream myfile;
    myfile.open(fullpath.c_str());
    if (myfile.is_open()) {

        vector <double> r (3);
        vector <double> v (3);
        vector <double> f (3);      // new value
        vector <double> f_ (3);     // old value
        int atomID;
        int molID;
        int type;
        int timestep;
        int Natoms;
        string line;
        string Time;
        string Natomz;
        string xloxhi;
        string yloyhi;
        string zlozhi;
        atom some_atom;

        getline(myfile,line);
        getline(myfile,Time);
        getline(myfile,line);
        getline(myfile,Natomz);
        getline(myfile,line);
        getline(myfile,xloxhi);
        getline(myfile,yloyhi);
        getline(myfile,zlozhi);
        getline(myfile,line);

        timestep = stoi(Time);
        Natoms = stoi(Natomz);
        vector <string> XloXhi = split_on_whitespace(xloxhi);
        vector <string> YloYhi = split_on_whitespace(yloyhi);
        vector <string> ZloZhi = split_on_whitespace(zlozhi);

        Xlo = stof(XloXhi[0]);
        Xhi = stof(XloXhi[1]);
        Ylo = stof(YloYhi[0]);
        Yhi = stof(YloYhi[1]);
        Zlo = stof(ZloZhi[0]);
        Zhi = stof(ZloZhi[1]);

        while(!myfile.eof()){
            myfile >> atomID;
            myfile >> molID;
            myfile >> type;

            for (int i = 0; i < 3; i++){ myfile >> r[i];}
            for (int i = 0; i < 3; i++){ myfile >> v[i];}
            for (int i = 0; i < 3; i++){ myfile >> f[i];}

            if (f[0] == f_[0]){
                if (f[1] == f_[1]){
                    if (f[2] == f_[2]){
                        break;
                    }
                }
            }
            for (int i = 0; i < 3; ++i) {
                f_[i] = f[i];
            }

            some_atom.set_atom_initial_state(atomID,type,r,v,f);
            some_atom.set_part_of_molecule(molID);
            //some_atom.set_systemsize(xlo,xhi,ylo,yhi,zlo,zhi);
            Atoms.push_back(some_atom);
        }
        myfile.close();
        TIME.push_back(timestep);
        NATOMS.push_back(Natoms);
    }
}

void ReadLAMMPSdump::read_newfile(string &filename, vector <atom> &Atoms){
    string fullpath;
    ENDCHAR = PATH.back();
    if (ENDCHAR == "/"){fullpath = PATH + filename;}
    else {fullpath = PATH + "/" + filename;}

    int Natoms;
    vector <atom> newAtoms;
    vector <double> r (3);
    vector <double> v (3);
    vector <double> f (3);          // new value

    ifstream myfile;
    myfile.open(fullpath.c_str());
    CLOCK = clock();
    if (myfile.is_open()) {

        vector <double> f_ (3);     // old value
        int atomID;
        int molID;
        int type;
        int timestep;
        string line;
        string Time;
        string Natomz;
        string xloxhi;
        string yloyhi;
        string zlozhi;
        atom some_atom;

        getline(myfile,line);
        getline(myfile,Time);
        getline(myfile,line);
        getline(myfile,Natomz);
        getline(myfile,line);
        getline(myfile,xloxhi);
        getline(myfile,yloyhi);
        getline(myfile,zlozhi);
        getline(myfile,line);

        timestep = stoi(Time);
        Natoms = stoi(Natomz);

        while(!myfile.eof()){
            myfile >> atomID;
            myfile >> molID;
            myfile >> type;

            for (int i = 0; i < 3; i++){ myfile >> r[i];}
            for (int i = 0; i < 3; i++){ myfile >> v[i];}
            for (int i = 0; i < 3; i++){ myfile >> f[i];}

            if (f[0] == f_[0]){
                if (f[1] == f_[1]){
                    if (f[2] == f_[2]){
                        break;
                    }
                }
            }
            for (int i = 0; i < 3; ++i) {
                f_[i] = f[i];
            }

            some_atom.set_atom_initial_state(atomID,type,r,v,f);
            some_atom.set_part_of_molecule(molID);
            newAtoms.push_back(some_atom);
        }
        myfile.close();
        TIME.push_back(timestep);
        NATOMS.push_back(Natoms);
    }

    int index1; // index of atom in this file
    int index2; // index of atom in initial file
    for (int i=0;i< Natoms;i++){
        index1 = newAtoms[i].get_atomID();
        for (int j=0; j< Natoms; j++){
            index2 = Atoms[j].get_atomID();
            if (index2 == index1){
                r = newAtoms[i].get_position();
                v = newAtoms[i].get_velocity();
                f = newAtoms[i].get_force();
                Atoms[j].update_position(r);
                Atoms[j].update_velocity(v);
                Atoms[j].update_force(f);

            } // end if
        }  // end for j
    }  // end for i

    CLOCK = clock() - CLOCK;
    timeusage = float(CLOCK)/CLOCKS_PER_SEC;

} // end read_newfile

vector <double> ReadLAMMPSdump::get_systemsize(){
    vector <double> systemsize (6,0.0);
    systemsize[0] = Xlo;
    systemsize[1] = Xhi;
    systemsize[2] = Ylo;
    systemsize[3] = Yhi;
    systemsize[4] = Zlo;
    systemsize[5] = Zhi;
    return systemsize;
}

int ReadLAMMPSdump::get_Natoms(){
    return NATOMS[0];
}


string ReadLAMMPSdump::get_path(){
    return PATH;
}

vector <string> ReadLAMMPSdump::split_on_whitespace(const string &str){
    stringstream ss(str); // Insert the string into a stream
    vector<string> tokens; // Create vector to hold our words
    string buf;
    while (ss >> buf){
        tokens.push_back(buf);
    }
    return tokens;
}

double ReadLAMMPSdump::get_timeusage(){
    return timeusage;
}

int ReadLAMMPSdump::get_timestep(){
    return TIME.back();  // return last element in TIME vector
}

