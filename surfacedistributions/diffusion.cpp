#include "diffusion.h"

diffusion::diffusion(
        vector <atom> &atoms, int &atomtype,int &Natoms,
        double &SystemStartValue, double &SystemEndValue, string &dir, double &DL){
    /* dir == xdir,ydir,zdir = the direction to divide the system into bins
       track all atoms within a box and calculate the diffusion of the particles of type type
       SystemStartValue is the location of a lower surface (if there is any) in the dir direction. It defines the inner region
       SystemEndValue is the location of a upper surface (if there is any) in the dir direction. It defines the inner region
    */

    Start = SystemStartValue;
    End = SystemEndValue;
    Type = atomtype;
    int type;
    int boxID;
    double ldist,udist,dist;
    double Pos;
    vector <double> pos (3,0.0);

    if ( Start > End) { throw std::invalid_argument( "recieved SytstemStartValue > SystemEndValue" );}
    if (dir == "xdir"){DIR = 0;}
    if (dir == "ydir"){DIR = 1;}
    if (dir == "zdir"){DIR = 2;}

    L = (End - Start);                               // cleavage widht
    if (DL < 0.001){ dL = int(L/2.0);}               // dL is set to half system size <=> one box
    else{ dL = DL;}                                  // else, dL is set to whatever DL is given as

    Nboxes = ceil(L/(2.0*dL));                       // create more boxes than dL allows, to ensure that all particles are accouted for
    //Nboxes = int((L/dL)/2.0);                        // number of boxes to use // only need boxes up to the center of system
    dL = (N/(2.0*L));                                // the actual dL becomes
    center = SystemStartValue + L/2.0;               // center of Start-End system
    vector <int> counter (Nboxes,0);                 // holds number of particles in every box

    for (int i=0; i<Natoms; i++){
        type = atoms[i].get_atom_type_number();
        pos = atoms[i].get_position();
        if (type == Type){
            atoms[i].set_past_position(pos);
            Pos = pos[DIR];
            ldist = Pos - Start;                      // distance from lower surface
            udist = End - Pos;                        // distance from upper surface
            if (ldist < udist && ldist > 0){          // the particle is nearer the lower surface
                dist = ldist;
            }
            else{
                dist = udist;                         // the particle is closer to the upper surface
            }

            boxID = floor(dist/dL);                   // put particle in its right boxIDainer

            if (boxID >= Nboxes){                     // if for some reason particle is out of region, put it in last container
                //cout << "atom outside system! Pos= " << Pos << " dist=" << dist << " boxID=" << boxID << " Nboxes= " << Nboxes <<  endl;
                boxID = Nboxes-1;
            }
            if (boxID < 0) { boxID = 0;}

            atoms[i].set_box_number(boxID);           // set the index of what box the atom belongs to!
            counter[boxID] = counter[boxID] + 1;      // add one particle to container!

        } // end if
        atoms[i].set_past_position(pos);
    } // end for


    for (int i=0; i<Nboxes;i++){
        Counted = Counted + counter[i];
    }
    Container = counter;

    NATOMS = Natoms;               // number of atoms in total (total in system)

}

void diffusion::calculate_diffusion(
        vector <atom> &atoms, ReadLAMMPSdump &PATH, string &frontname,
        double &Lx, double &Ly, double &Lz, int &dt, int &Tstart, int &Tend){

    clock_t TIME,TIME2;
    TIME = clock();                  // start time
    string filename;                 // filename
    int boxID, NEWboxID;             // initial index and index now
    int type;                        // atom type
    int t = 0;                       // timestep
    Dt = dt;                         // timestep used
    TimeStart = Tstart;              // starting timestep after the timestep at the origin!
    TimeEnd = Tend;                  // timestep at end
    Ntimesteps = (Tend-Tstart)/dt;   // number of timesteps
    double Pos;                      // position of interesting direction dir
    double ldist,udist,dist;         // distance to surface in direction of dir. ldist, udist = dist to lower and upper surface
    double dx,dy,dz;                 // displacement
    double t_readfile;               // time usage reading file
    double accumulate_time_readfile; // for averaging time usage read file
    double accumulate_time_calculate;// for averaging time usage diffusion calculation
    vector <double> ipos (3,0.0);    // position at timestep t0
    vector <double> pos (3,0.0);     // position at timestep t
    vector <double> pastpos (3,0.0); // position at timesteo t-1
    vector <vector <double> > MSD (Nboxes,vector <double> (Ntimesteps,0.0));  // msd in boxes from surface
    double avglostatoms = 0;
    int lost = 0;

    vector <int> counter (Nboxes,0);      // number of particles within a box
    vector <double> msd (Nboxes,0.0);     // mean square displacement in a box

    for (int time=Tstart; time<Tend; time = time + dt){
        stringstream ss;
        ss << frontname << "." << time << ".txt" ;
        ss >> filename;
        //cout << "  Reading file: " << filename << endl;             //" path= " << endl; //path << " time=" << time << endl;
        PATH.read_newfile(filename,atoms);                        // update particle information!

        lost = 0;                             // number of lost atoms set to zero
        TIME2 = clock();
        for (int i=0; i<NATOMS; i++){
            type = atoms[i].get_atom_type_number();            // follow correct atomtype
            pos = atoms[i].get_position();
            if (type == Type){
                boxID = atoms[i].get_box_number();             // what box does the atom belong to?
                pastpos = atoms[i].get_past_position();
                ipos = atoms[i].get_initial_position();
                Pos = pos[DIR];
                ldist = Pos - Start;
                udist = End - Pos;

                if (ldist < udist && ldist > 0){
                    dist = ldist;
                }
                else{
                    dist = udist;
                }

                NEWboxID = floor(dist/dL);

                if (NEWboxID == boxID){
                    // we will only go further with the atoms that are still in the same boxIDainer!
                    dx = pos[0] - ipos[0];
                    dy = pos[1] - ipos[1];
                    dz = pos[2] - ipos[2];
                    MinimumImageConvention(dx,Lx);
                    MinimumImageConvention(dy,Ly);
                    MinimumImageConvention(dz,Lz);

                    msd[boxID] = msd[boxID] + (dx*dx + dy*dy + dz*dz);
                    counter[boxID] = counter[boxID] + 1; // count the number of particles within every box

                }
                else {
                    lost = lost + 1;                     // we lost an atom!
                }
            } // follow correct types
            atoms[i].set_past_position(pos);
        } // particle-loop
        Container_e = counter;                           // number of particles in every container at end
        for (int c=0; c<Nboxes; c++){
            if (counter[c] > 0){ MSD[c][t] = MSD[c][t] + msd[c]/(6*counter[c]);} // add square displacement
        }

        fill(counter.begin(), counter.end(), 0);  // empty vector
        fill(msd.begin(), msd.end(), 0);          // empty vector

        t = t + 1;
        avglostatoms = avglostatoms + lost;
        TIME2 = clock() - TIME2;
        t_readfile = PATH.get_timeusage();
        accumulate_time_readfile = accumulate_time_readfile + t_readfile;
        accumulate_time_calculate = accumulate_time_calculate + float(TIME2)/CLOCKS_PER_SEC;
        //cout << "    Time usage timeloop: " << float(TIME2)/CLOCKS_PER_SEC << " s. Time useage readfile: " << t_readfile << " s." << endl;
        //cout << "# ------------------------------------------------------------------------------------------------ #" << endl;
    }  // time-loop

    // time-loop is finished! Now we have to store the data! write to file :-)

    avglostatoms = avglostatoms/t;
    stringstream ss;
    ss << "msd_results_" << frontname << "." << TimeStart << "-" << TimeEnd << ".txt";
    ss >> filename;
    string path = PATH.get_path();
    write_to_file(filename, path, MSD);
    TIME = clock() - TIME;

    cout << "# _-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_ #" << endl;
    cout << "Total timeusage: " << float(TIME)/CLOCKS_PER_SEC << " seconds.  Avg nr of lost atoms per timestep = " << avglostatoms << endl;
    cout << "Avg time usage reading files: " << accumulate_time_readfile/t << " s. Avg time usage calculating the diffusion: " << accumulate_time_calculate/t << " s." << endl;

}

void diffusion::write_to_file(string &filename, string &path, vector <vector <double> > &MSD){

    ofstream thisfile;          // write to file
    thisfile.open(filename);
    thisfile << "# Diffusion of atoms: type " << Type << endl;
    thisfile << "# Database: " << path << endl;
    thisfile << "# time: Tstart Tend dt Ntimesteps " << TimeStart << " " << TimeEnd << " " << Dt << " " << Ntimesteps << endl;
    thisfile << "# cleavage: L dL Nboxes " << L << " " << dL << " " << Nboxes << endl;
    double Poresize = End - Start;
    thisfile << "# Poresize systemstart systemend " << Poresize << " " << Start << " " << End << endl;
    thisfile << "# boxID Nstart Nend MSD(t) " << endl;

    for (int i=0; i<Nboxes; i++ ){
        thisfile << i << " " << Container[i] << " " << Container_e[i] << " " ;
        for (int t=0; t<Ntimesteps; t++ ){
            thisfile << MSD[i][t] << " " ;
        }
        thisfile << endl;
    }
    thisfile.close();
}

int diffusion::counted_particles(){
    return Counted;
}

vector <int> diffusion::particle_distribution(){
    return Container;
}

void diffusion::MinimumImageConvention(double &delta, double &systemsize){
    if (delta < -(systemsize/2.0) ){ delta = delta + systemsize;}
    if (delta > (systemsize/2.0) ){ delta = delta - systemsize;}
}


