/*
 * simulation.h
 *
 *  Created on: 14/01/2014
 *      Author: Ran
 */

#ifndef SIMULATION_H_
#define SIMULATION_H_
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
//#include <vector>
using namespace std;
/*
 * #include <iostream>
 * #include <string>
 * #include <fstream>
 * #include <iostream>
 * //#include <gsl/gsl_rng.h>
 * //#include <gsl/gsl_randist.h>
 */
struct neuron_t {
    int* conn_id_out_deg;    // list of neurons which are connected
    double* conn_J_out_deg;  // list of strenght of connections
    int nb_conn_out_deg;  // NUMBER OF CONNECTIONS
    
    double* conn_J_in_degE;  // list of strenght of connections
    double* conn_rho_in_degE;  // list of rhos
    double* conn_alpha_in_degE;  // list of rhos
    int nb_conn_in_degE;  // NUMBER OF CONNECTIONS
    int* conn_id_in_degE;    // list of neurons which are connected
    
    double* conn_J_in_degI;  // list of strenght of connections
    double* conn_rho_in_degI;  // list of rhos
        double* conn_alpha_in_degI;  // list of rhos
    int nb_conn_in_degI;  // NUMBER OF CONNECTIONS
    int* conn_id_in_degI;    // list of neurons which are connected
    
    
    // std::string type; // Type: E/I optional: which network
    double h_EE;
    double h_IE;
    double h_EI;
    double h_II;
    double V;
};


class Simulation {
    
    
private:
    // Variables:
    
    double dt;
    int steps_tot;
    double time_sim;
    
    neuron_t * C;
    FILE * fid_connectivity;
    FILE * fid_spikee;
    FILE * fid_spikei;
    FILE * fid_voltagee;
    FILE * fid_voltagei;
    FILE * fid_he;
    FILE * fid_hi;
    FILE * fid_currente;
    FILE * fid_currenti;
    
    int save_current,save_voltage,save_tot_current;
    
    //  std::string strE ("E");
    // std::string strE ("I");
    int n, ni, ne;
    int TstopStim;
    int TstopStim2;
    int TstartStim2;
    double RLTrial;
    int ki, ke, kff;
    double* V; double* h_E; double* h_I;
    double* hE_factor1; double* hI_factor1;double* hE_factor2;double* hI_factor2;double* v_factor1;double* v_factor2;
    double* thresh;
    double* V_old;
    double* tau_E; double* tau_I; double* tau_mem;
    double  h_E_i, h_I_i, v_i;
    double* spk_t;
    double* Iext;
    double* Iext2;
    double  Ie, Ii, ie, ii;
    double  I2e, I2i, i2e, i2i;
    double  sd_Ie,sd_Ii;
    double  thresh_par;
    double  tme, tmi, tee,tei, tie, tii;
    double  gee,gei, gie, gii;
    double  VresetE,VresetI, Vspike;
    double* C_m; double* V_L;double* Vreset; double* g_L; double* g_Na; double* V_Na; double* g_K; double* V_K;
    double* g_Adap; double* V_Adap;
    double C_mE, C_mI, V_LE, g_LE, V_NaE, g_NaE, V_KE, g_KE, V_AdapE, g_AdapE;
    double V_LI, g_LI, V_NaI, g_NaI, V_KI, g_KI, V_AdapI, g_AdapI;
    double tau_Adap;
    double  V_E, V_I;
    double* h_Na;
    double* n_K;
    double* z_Aadp;
    double* I_inp;
    
    double sdV_LI, sdV_LE;
    double sd_thresh_I, sd_thresh_E;
    
    double rho;
    
    double* u0;
    double* u1;
    double t0;
    int 	 P ; // u=(V,m,h,n,z) HERE????
    // random number
    gsl_rng * rand_eng;
    
    //save variables
//    double* Ve_save;
    vector<double> inp;
    vector< vector<double> > Ve_save;
    vector< vector<double> > Vi_save;
    vector< vector<double> > he_save;
    vector< vector<double> > hi_save;
    vector< vector<double> > Cure_save;
    vector< vector<double> > Curi_save;
    
    vector< vector<double> > gsyn;
    
    
    int cells_to_save;
    
public:
    Simulation();
    //virtual ~Simulation();
    void load();
    void init();
        void init_ic(int iter);
    void run();
    void close();
};

#endif /* SIMULATION_H_ */
