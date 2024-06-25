/*
 * simulation.cpp
 *
 *  Created on: 05/17/2018
 *  main changes: Added distribution of VL, added input as Input.txt, Added O(1) distribution of Iexternal
 *  Main changes 24/1/2019: Conductance and Current based IF. Added uniform distributon for thresholds
 *      Author: Ran
 */
#include <iostream>
#include <cstring>
#include <cmath>
#include <ctime>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include "Simulation.h"
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <yaml-cpp/yaml.h>

// add the rk4 functions
# include "my_rk4.h"
double *rk4_fun ( double t, int n, double u[], double C_m, double g_L, double V_L, double g_Na, double V_Na,
        double g_K, double V_K, double g_Adap, double V_Adap, double tau_Adap, double I_inp  );

using namespace std;


Simulation::Simulation() {
}

void Simulation::load() {
    
//Choose the average rates of E and I neurons
//E/I>gei/gii>gee/gie Balance conditions
//Firing rate of the E and I neurons (For large K) should be:
//v_I=(gie*E-gee*I)/(gie*gei-gii*gee)
//v_E=(gii*E-gei*I)/(gie*gei-gii*gee)
    
    string yaml_file_name = "par.yaml";
    
    YAML::Node doc = YAML::LoadFile(yaml_file_name.c_str());
    
    time_sim  = doc["time_sim"].as<int>();
    dt   	  = doc["dt"].as<double>();
    ni 	  = doc["ni"].as<int>();
    ne 	  = doc["ne"].as<int>();
    n=ne+ni;
    ke         = doc["ke"].as<int>();
    ki         = doc["ki"].as<int>();
    kff         = doc["kff"].as<int>();
    thresh_par= doc["thresh_par"].as<double>();
    gee 	  = doc["gee"].as<double>();
    gei 	  = doc["gei"].as<double>();
    gie 	  = doc["gie"].as<double>();
    gii       = doc["gii"].as<double>();
    tme       = doc["tme"].as<double>();
    tmi 	  = doc["tmi"].as<double>();
    tee 	  = doc["tee"].as<double>();
    tei 	  = doc["tei"].as<double>();
    tie 	  = doc["tie"].as<double>();
    tii       = doc["tii"].as<double>();
    ie        = doc["ie"].as<double>();
    ii 	  = doc["ii"].as<double>();
    i2e        = doc["i2e"].as<double>();
    i2i 	  = doc["i2i"].as<double>();
    sd_Ie 	  = doc["sd_Ie"].as<double>();
    sd_Ii 	  = doc["sd_Ii"].as<double>();
    
    
    VresetE 	  = doc["VresetE"].as<double>();
    VresetI 	  = doc["VresetI"].as<double>();
    Vspike    = doc["Vspike"].as<double>();
    
    C_mE	  = doc["C_mE"].as<double>();
    C_mI	  = doc["C_mI"].as<double>();
    V_LE	  = doc["V_LE"].as<double>();
    g_LE	  = doc["g_LE"].as<double>();
    V_NaE	  = doc["V_NaE"].as<double>();
    g_NaE	  = doc["g_NaE"].as<double>();
    V_KE	  = doc["V_KE"].as<double>();
    g_KE	  = doc["g_KE"].as<double>();
    V_AdapE	  = doc["V_AdapE"].as<double>();
    g_AdapE	  = doc["g_AdapE"].as<double>();
    
    V_LI	  = doc["V_LI"].as<double>();
    g_LI	  = doc["g_LI"].as<double>();
    V_NaI	  = doc["V_NaI"].as<double>();
    g_NaI	  = doc["g_NaI"].as<double>();
    V_KI	  = doc["V_KI"].as<double>();
    g_KI	  = doc["g_KI"].as<double>();
    V_AdapI	  = doc["V_AdapI"].as<double>();
    g_AdapI	  = doc["g_AdapI"].as<double>();
    
    V_E	  = doc["V_E"].as<double>();
    V_I	  = doc["V_I"].as<double>();
    //  rho	  = doc["rho"].as<double>();
    cells_to_save = doc["n_save"].as<int>();
    
    sdV_LI	  = doc["sdV_LI"].as<double>();
    sdV_LE	  = doc["sdV_LE"].as<double>();
    
    sd_thresh_I	  = doc["sd_thresh_I"].as<double>();
    sd_thresh_E	  = doc["sd_thresh_E"].as<double>();
    TstopStim	  = doc["TstopStim"].as<int>();
    TstopStim2= doc["TstopStim2"].as<int>();
    TstartStim2= doc["TstartStim2"].as<int>();
    
    RLTrial	  = doc["RLTrial"].as<double>();
    
    
    
    tau_Adap  = doc["tau_Adap"].as<double>(); // tau_Adap is te same fo E and I!!!
    P=4; ///????
    
    save_current = doc["save_current"].as<int>();
    save_voltage = doc["save_voltage"].as<int>();
    save_tot_current = doc["save_tot_current"].as<int>();
    Ii =ii*sqrt(kff);
    Ie =ie*sqrt(kff);
    I2i =i2i*sqrt(kff);
    I2e =i2e*sqrt(kff);
    
    gee=-gee*(V_LE-V_E);
    gie=-gie*(V_LE-V_E);
    gii=-gii*(V_LI-V_I);
    gei=-gei*(V_LI-V_I);
    
    double tmp_v_I=(gie*ie*(V_E-V_LE)-gee*ii*(V_E-V_LE))/(gie*abs(gei)-abs(gii)*gee);
    double tmp_v_E=(abs(gii)*ie*(V_E-V_LE)-abs(gei)*ii*(V_E-V_LE))/(gie*abs(gei)-abs(gii)*gee);
    
    //  double tmp_v_I=(gie*ie-gee*ii)/(gie*abs(gei)-abs(gii)*gee);
    //  double tmp_v_E=(abs(gii)*ie-abs(gei)*ii)/(gie*abs(gei)-abs(gii)*gee);
    
    cout<<"Average firing rates:\n"<<"E="<<tmp_v_E<<"\t"<<"I="<<tmp_v_I<<"\n";        
} // end load

void Simulation::init() {
    
// open files to read\write
//string file_name = "./" + folder + "/data/" + name_id + "_spike.txt";




    
// Allocate vectors
    
    tau_mem    = new double[n];
    V_old      = new double[n];
    h_E        = new double[n];
    h_I        = new double[n];
    hE_factor1 = new double[n];
    hE_factor2 = new double[n];
    hI_factor1 = new double[n];
    hI_factor2 = new double[n];
    v_factor1  = new double[n];
    v_factor2  = new double[n];
    thresh     = new double[n];
    spk_t      = new double[n];
    Iext 	   = new double[n];
    Iext2 	   = new double[n];
    Vreset       = new double[n];
    
//Leak
    C_m        = new double[n];
    V_L        = new double[n];
    g_L        = new double[n];
//Na
    g_Na       = new double[n];
    V_Na 	   = new double[n];
//K
    g_K 	   = new double[n];
    V_K 	   = new double[n];
//Adapt
    g_Adap     = new double[n];
    V_Adap     = new double[n];
    
    V          = new double[n];
    h_Na       = new double[n];
    n_K        = new double[n];
    z_Aadp     = new double[n];
    
    
    I_inp      = new double[n];
    u0  	   = new double[P]; //Pointer to u1????????
    u0[0]      = 0.0;
    u0[1]      = 0.0;
    u0[2]      = 0.0;
    u0[3]      = 0.0;
    
    
    
    
    // Load input Input
    
    string line;
    ifstream myfile ("Input.txt");
    if (myfile.is_open())
    {
        while ( getline (myfile,line) )
        {
            stringstream ss(line);
            double ti;
            while(ss >> ti)
                inp.push_back(ti);
        }
        myfile.close();
    }
    else cout << "Unable to open file";
    for (int ii=299;ii<310;ii++)
        cout<<inp[ii]<<'\n';
    
    
    
    
//start the random generator
    const gsl_rng_type * T;
    gsl_rng_env_setup();
    T = gsl_rng_mt19937;
    rand_eng = new gsl_rng;
    rand_eng = gsl_rng_alloc(T);
//gsl_rng_set(rand_eng, (unsigned long) time(NULL) * getpid());
    double seed =1;
    gsl_rng_set(rand_eng, seed);
       
    steps_tot = int(time_sim / dt);
    cout<<steps_tot<<"\n";
    
    
    
    // CHANGE IT BECAUSE CURE AND I ARE NOT THAT BIG!
    if (save_voltage){
        //construct a matrix to save the files
        //cells_to_save = 30;
        Ve_save.resize(steps_tot+1);
        Vi_save.resize(steps_tot+1);
        Cure_save.resize(steps_tot+1);
        Curi_save.resize(steps_tot+1);
        
        for(int i = 0 ; i < steps_tot+1 ; ++i){
            Ve_save[i].resize(cells_to_save+1);
            Vi_save[i].resize(cells_to_save+1);
            Cure_save[i].resize( cells_to_save);
            Curi_save[i].resize( cells_to_save);
        }
    }
    
    
    if (save_tot_current){
        //construct a matrix to save the files
        //cells_to_save = 30;
        he_save.resize(steps_tot+1);
        hi_save.resize(steps_tot+1);
        
        for(int i = 0 ; i < steps_tot+1 ; ++i){
            he_save[i].resize(cells_to_save+1);
            hi_save[i].resize(cells_to_save+1);
        }
    }
    
    
    
    
    // Construct g_E, g_I matrix
    gsyn.resize(n+1);
    for(int i = 0 ; i < n ; ++i){
        gsyn[i].resize(n);
    }
    
    
    
    
    
    
// init synaptic and membrane time const and inputs
    tau_E = new double[n];
    tau_I = new double[n];
//ex neurons
    for (int i = 0; i < ne; i++) {
        Iext[i]      = Ie+gsl_ran_gaussian_ziggurat(rand_eng,sd_Ie);
        Iext2[i]      = I2e;
        tau_E[i]   = tee;
        tau_I[i]   = tei;
        tau_mem[i] = tme;
    }
//inh neurons
    for (int i = ne; i < n; i++) {
        Iext[i]      = Ii+gsl_ran_gaussian_ziggurat(rand_eng,sd_Ii);
        Iext2[i]      = I2i;
        tau_E[i]   = tie;
        tau_I[i]   = tii;
        tau_mem[i] = tmi;
    }
    
    
// inialize neuron stuct
    C = new neuron_t[n];
// init E neurons
    for (int i = 0; i < ne; i++) {
        C[i].nb_conn_out_deg =0;
        C[i].nb_conn_in_degE =0;
        C[i].nb_conn_in_degI =0;
        thresh[i]    =thresh_par-sqrt(3.0)*sd_thresh_E+2*sqrt(3.0)*sd_thresh_E*gsl_rng_uniform(rand_eng);;
        
        //Leak: ???  Check: ge an gi in HvV
        C_m[i]        = C_mE;
        V_L[i]        = V_LE+gsl_ran_gaussian_ziggurat(rand_eng,sdV_LE);
        Vreset[i]     = VresetE;
        g_L[i]        = g_LE;
        
        //Na
        g_Na[i]       = g_NaE;
        V_Na[i]       = V_NaE;
        
        //K
        g_K[i]        = g_KE;
        V_K[i]        = V_KE;
        
        g_Adap[i]    = g_AdapE;
        V_Adap[i]    = V_AdapE;
    }
    
// init I neurons %% BUG it was ne+1
    for (int i = ne; i < n; i++) {
        C[i].nb_conn_out_deg =0;
        C[i].nb_conn_in_degE =0;
        C[i].nb_conn_in_degI =0;
        thresh[i]    =thresh_par-sqrt(3.0)*sd_thresh_I+2*sqrt(3.0)*sd_thresh_I*gsl_rng_uniform(rand_eng);;
        
        //Leak: ???  Check: ge an gi in HvV
        C_m[i]        = C_mI;
        V_L[i]        = V_LI+gsl_ran_gaussian_ziggurat(rand_eng,sdV_LI);
        Vreset[i]     = VresetI;
        g_L[i]        = g_LI;
        
        //Na
        g_Na[i]       = g_NaI;
        V_Na[i]       = V_NaI;
        
        //K
        g_K[i]        = g_KI;
        V_K[i]        = V_KI;
        
        g_Adap[i]    = g_AdapI;
        V_Adap[i]    = V_AdapI;
    }
    
    // Read OutDeg matrix and store
    int tmp_from = 0;
    int tmp_to[n];
    double tmp_j[n];
    int tmp_nb_conn = 0;
    
    
    ifstream file("connectivityOutDeg.txt",ios::in);
    
    if (file.good())
    {
        string str;
        while(getline(file, str))
        {
            istringstream ss(str);
            int buf_from;
            int buf_to;
            double buf_J;
            //string buf_type;
            // read a line from file
            ss>>buf_from>>buf_to>>buf_J;//>>buf_type;
            // assign the type of the neuron (not efficient...)
            //C[buf_from].type=buf_type;
            //
            
            //check if from_neuron is new - if so, store the old neuron
            if (buf_from!=tmp_from) {
                //store the old info
                
                C[tmp_from].nb_conn_out_deg = tmp_nb_conn;
                C[tmp_from].conn_id_out_deg = new int[tmp_nb_conn];
                C[tmp_from].conn_J_out_deg = new double[tmp_nb_conn];
                for (int j = 0; j < tmp_nb_conn; j++) {
                    C[tmp_from].conn_id_out_deg[j] = tmp_to[j];
                    C[tmp_from].conn_J_out_deg[j] = tmp_j[j];
                }
                //init new
                tmp_from = buf_from;
                tmp_nb_conn = 0;
            }
            //store new connection
            tmp_to[tmp_nb_conn] = buf_to;
            tmp_j[tmp_nb_conn] = buf_J;
            tmp_nb_conn++;
        }
        
        //save the last connection
        C[tmp_from].nb_conn_out_deg = tmp_nb_conn;
        C[tmp_from].conn_id_out_deg = new int[tmp_nb_conn];
        C[tmp_from].conn_J_out_deg = new double[tmp_nb_conn];
        for (int j = 0; j < tmp_nb_conn; j++) {
            C[tmp_from].conn_id_out_deg[j] = tmp_to[j];
            C[tmp_from].conn_J_out_deg[j] = tmp_j[j];
        }
    } // end read connectivity
    
    
    // Read OutDeg matrix and store
    tmp_from = 0;
    int tmp_toE[ne];
    double tmp_jE[ne];
    double tmp_rhoE[ne];
    double tmp_alphaE[ne];
    
    int tmp_toI[ni];
    double tmp_jI[ni];
    double tmp_rhoI[ni];
    double tmp_alphaI[ni];
    int tmp_nb_connE = 0;
    int tmp_nb_connI = 0;
    
    
    ifstream fileB("connectivityInDeg.txt",ios::in);
    
    if (fileB.good()){
        string str;
        while(getline(fileB, str)){
            
            istringstream ss(str);
            int buf_from;
            int buf_to;
            double buf_J;
            double buf_rho;
            double buf_alpha;
            // read a line from file
            ss>>buf_from>>buf_to>>buf_J>>buf_rho>>buf_alpha;
            
            //check if from_neuron is new - if so, store the old neuron
            
            if (buf_from==tmp_from){
                if(buf_to<ne){
                    tmp_toE[tmp_nb_connE] = buf_to;
                    tmp_jE[tmp_nb_connE] = buf_J;
                    tmp_rhoE[tmp_nb_connE] = buf_rho;
                    tmp_alphaE[tmp_nb_connE] = buf_alpha;
                    tmp_nb_connE++;
                } else{
                    tmp_toI[tmp_nb_connI] = buf_to;
                    tmp_jI[tmp_nb_connI] = buf_J;
                    tmp_rhoI[tmp_nb_connI] = buf_rho;
                    tmp_alphaI[tmp_nb_connI] = buf_alpha;
                    tmp_nb_connI++;
                }
            }else{
                // Store the previous neuron
                C[tmp_from].nb_conn_in_degE = tmp_nb_connE;
                C[tmp_from].conn_id_in_degE = new int[tmp_nb_connE];
                C[tmp_from].conn_J_in_degE = new double[tmp_nb_connE];
                C[tmp_from].conn_rho_in_degE = new double[tmp_nb_connE];
                C[tmp_from].conn_alpha_in_degE = new double[tmp_nb_connE];
                C[tmp_from].nb_conn_in_degI = tmp_nb_connI;
                C[tmp_from].conn_id_in_degI = new int[tmp_nb_connI];
                C[tmp_from].conn_J_in_degI = new double[tmp_nb_connI];
                C[tmp_from].conn_rho_in_degI = new double[tmp_nb_connI];
                C[tmp_from].conn_alpha_in_degI = new double[tmp_nb_connI];
                
                for (int j = 0; j < tmp_nb_connE; j++) {
                    C[tmp_from].conn_id_in_degE[j] = tmp_toE[j];
                    C[tmp_from].conn_J_in_degE[j] = tmp_jE[j];
                    C[tmp_from].conn_rho_in_degE[j] = tmp_rhoE[j];
                    C[tmp_from].conn_alpha_in_degE[j] = tmp_alphaE[j];
                }
                for (int j = 0; j < tmp_nb_connI; j++) {
                    C[tmp_from].conn_id_in_degI[j] = tmp_toI[j];
                    C[tmp_from].conn_J_in_degI[j] = tmp_jI[j];
                    C[tmp_from].conn_rho_in_degI[j] = tmp_rhoI[j];
                    C[tmp_from].conn_alpha_in_degI[j] = tmp_alphaI[j];
                }
                tmp_nb_connE=0;
                tmp_nb_connI=0;
 // Store the new neuron (should be only E)               
                if(buf_to<ne){
                    tmp_toE[tmp_nb_connE] = buf_to;
                    tmp_jE[tmp_nb_connE] = buf_J;
                    tmp_rhoE[tmp_nb_connE] = buf_rho;
                    tmp_alphaE[tmp_nb_connE] = buf_alpha;
                    tmp_nb_connE++;
                } else{
                    tmp_toI[tmp_nb_connI] = buf_to;
                    tmp_jI[tmp_nb_connI] = buf_J;
                    tmp_rhoI[tmp_nb_connI] = buf_rho;
                    tmp_alphaI[tmp_nb_connI] = buf_alpha;
                    tmp_nb_connI++;
                }
 
            }
            tmp_from=buf_from;
        }
    } // end read connectivity
} // end init

void Simulation::init_ic(int iter) {
 
  // HERE START THE LOOP AND CHANGE THE SEED
  

  //  int iter=2;
    char numstr[21];
    sprintf(numstr,"%d",iter);
    
  // CHANGE THE NAMES!
      string file_name_tmp  = "./Data/spikee";
      string ending =".txt";
      string file_name_e;
      file_name_e= file_name_tmp+numstr+ending;  
            
      cout<<file_name_e<<'\n';
    fid_spikee          = fopen(file_name_e.c_str(), "w");
    string file_name_tmpi  = "./Data/spikei";
    string file_name_i;
     file_name_i= file_name_tmpi+numstr+ending;
     
    fid_spikei  	    = fopen(file_name_i.c_str(), "w");
   if (save_current) {
        string file_name_tmp = "./Data/cure";
        string file_name;
        file_name= file_name_tmp+numstr+ending;
        fid_currente     = fopen(file_name.c_str(), "w");
        
        file_name_tmp        = "./Data/curi";
        file_name= file_name_tmp+numstr+ending;
        fid_currenti     = fopen(file_name.c_str(), "w");
    }
    if (save_voltage){
        string file_name_tmp = "./Data/ve";
        string file_name;
        file_name= file_name_tmp+numstr+ending;        
        fid_voltagee     = fopen(file_name.c_str(), "w");
        file_name_tmp        = "./Data/vi";
        file_name= file_name_tmp+numstr+ending;        
        fid_voltagei	   = fopen(file_name.c_str(), "w");
    }
    
    if (save_tot_current){
        string file_name_tmp = "./Data/he";
        string file_name;
        file_name= file_name_tmp+numstr+ending;        
        fid_he     = fopen(file_name.c_str(), "w");
        file_name_tmp        = "./Data/hi";
        file_name= file_name_tmp+numstr+ending;        
        fid_hi	   = fopen(file_name.c_str(), "w");
    }
 
// initial conditions for currents and voltage: Changed in Conductance Vs LIF: BE CARFULL!!!
    for (int i = 0; i < n; i++) {
        V[i]      = 20.0*gsl_rng_uniform(rand_eng)+Vreset[i]-10.0;
//                V[i]      = 0.1*gsl_rng_uniform(rand_eng)+Vreset[i]+0.5;
        h_Na[i]   = gsl_rng_uniform(rand_eng)*0.05+0.4;
        n_K[i]    = gsl_rng_uniform(rand_eng)*0.05+0.4;
        z_Aadp[i] = gsl_rng_uniform(rand_eng)*0.05+0.4;
//     		h_Na[i]   = 0.9;
        //   		n_K[i]    = 0.1;
        //  		z_Aadp[i] = 0.1;
        h_I[i]    = gsl_rng_uniform(rand_eng); ///AAA
        h_E[i]    = gsl_rng_uniform(rand_eng); ///AAA
        
        for (int j = 0; j < n; j++) {
            gsyn[i][j]    = gsl_rng_uniform(rand_eng)/sqrt(kff); ///AAA
        }                
        I_inp[i]  = gsl_rng_uniform(rand_eng);
        hE_factor1[i] = -(dt )*(1.0 / tau_E[i]); // 0:ne-1 tau_EE; ne:n tau_EI
        hI_factor1[i] = -(dt )*(1.0 / tau_I[i]); // 0:ne-1 tau_IE; ne:n tau_II
        v_factor1[i]  = (dt )*(1.0 / tau_mem[i]); // 0:ne-1 tau_m_E; ne:n tau_m_I
    }
    std::cout << " Neurons loaded" << endl;
} // end init_ic




void Simulation::run() {
    
    std::cout << "Simulation started" << endl;
    clock_t ticks = clock();    
    for (int step = 1; step < (steps_tot + 1); step++) {
        if ((step % (steps_tot / 100)) == 0) {
            std::cout << "Time elapsed: " << ((clock() - ticks) / CLOCKS_PER_SEC) << endl;
        }
       // cout<<step<<endl;        
// Euler step for E neurons with external input that depends on time
        for (int i = 0; i < ne; i++) {
            // Calculate the recurrent input current for neuron i,E
            // Exc current
            h_I[i]=0.0;h_E[i]=0.0;            
            for (int j=0; j<C[i].nb_conn_in_degE;j++) {
                int tmp_conn_id = C[i].conn_id_in_degE[j];
                gsyn[i][tmp_conn_id] += hE_factor1[i] * (gsyn[i][tmp_conn_id]);
                h_E[i] +=-gsyn[i][tmp_conn_id]*C[i].conn_alpha_in_degE[j]*(   C[i].conn_rho_in_degE[j]*(V[i]-V_E)+(1.0-C[i].conn_rho_in_degE[j])*(V_L[i]-V_E)    );
            }            
            // Inh current
            for (int j=0; j<C[i].nb_conn_in_degI;j++) {
                int tmp_conn_id = C[i].conn_id_in_degI[j];
                gsyn[i][tmp_conn_id] += hE_factor1[i] * (gsyn[i][tmp_conn_id]);
                h_I[i] +=-gsyn[i][tmp_conn_id]*C[i].conn_alpha_in_degI[j]*(   C[i].conn_rho_in_degI[j]*(V[i]-V_I)+(1.0-C[i].conn_rho_in_degI[j])*(V_L[i]-V_I)    );
            }
            // Add the feedforward
            I_inp[i]  = -(Iext[i]+inp[step]*(Iext2[i]+( 1.0*RLTrial*(step<(TstopStim/dt)  ) - 1.0*RLTrial*(  (step<(TstopStim2/dt))&&(step>(TstartStim2/dt))  )    )*((i > (ne+ni/2)) - (i < (ne+ni/2))))  )*( (V_L[i]-V_E)    );
            //V_old[i]	=V[i];
            V[i]   += v_factor1[i]  * (-(V[i]-V_L[i]) +I_inp[i]+h_E[i]+h_I[i]); //voltage
        }        
// Euler step for I neurons
        for (int i = ne; i < n; i++) {
            // Exc current
            h_I[i]=0.0;h_E[i]=0.0;
            for (int j=0; j<C[i].nb_conn_in_degE;j++) {
                int tmp_conn_id = C[i].conn_id_in_degE[j];
                gsyn[i][tmp_conn_id] += hE_factor1[i] * (gsyn[i][tmp_conn_id]);
                h_E[i] +=-gsyn[i][tmp_conn_id]*C[i].conn_alpha_in_degE[j]*(   C[i].conn_rho_in_degE[j]*(V[i]-V_E)+(1.0-C[i].conn_rho_in_degE[j])*(V_L[i]-V_E)    );
            }
            // Inh current            
            for (int j=0; j<C[i].nb_conn_in_degI;j++) {
                int tmp_conn_id = C[i].conn_id_in_degI[j];
                gsyn[i][tmp_conn_id] += hE_factor1[i] * (gsyn[i][tmp_conn_id]);
                h_I[i] +=-gsyn[i][tmp_conn_id]*C[i].conn_alpha_in_degI[j]*(   C[i].conn_rho_in_degI[j]*(V[i]-V_I)+(1.0-C[i].conn_rho_in_degI[j])*(V_L[i]-V_I)    );
            }            
            // External current
            I_inp[i] =-(Iext[i]+inp[step]*(Iext2[i]+( 1.0*RLTrial*(step<(TstopStim/dt)  ) - 1.0*RLTrial*(  (step<(TstopStim2/dt))&&(step>(TstartStim2/dt))  )    )*((i > (ne+ni/2)) - (i < (ne+ni/2))))  )*( (V_L[i]-V_E)    );
            //V_old[i]	=V[i];
            V[i]   += v_factor1[i]  * (-(V[i]-V_L[i]) +I_inp[i]+h_E[i]+h_I[i]); //voltage
        }
// end of Euler step
        
// RK4 step P=4 (V,h_Na,n_K,z_Adap)
        /*
         * for (int i = 0; i < n; i++) {
         *
         * //generate u0
         * u0[0]   = V[i];
         * u0[1]   = h_Na[i];
         * u0[2]   = n_K[i];
         * u0[3]   = z_Aadp[i];
         * V_old[i]	= V[i];
         *
         *
         * // MAYA: How rk4vec_fun knows who are its variables?
         * t0=dt*step;
         * u1= my_rk4vec ( t0, P, u0, dt,C_m[i], g_L[i], V_L[i], g_Na[i], V_Na[i], g_K[i], V_K[i], g_Adap[i], V_Adap[i], tau_Adap, I_inp[i], rk4_fun);
         *
         *
         *
         *
         * //store the RK4 step
         * V[i]       = u1[0];
         * h_Na[i]    = u1[1];
         * n_K[i]     = u1[2];
         * z_Aadp[i]  = u1[3];
         *
         * delete [] u1;
         * }
         */
        
        //Update currents due to spiking neurons
        for (int i = 0; i < n; i++) {
            //is the neuron above threshold (Different with the H&H)
            if ( (V[i] > thresh[i])) {                
                // correct voltage at t+dt due to the spike- for the LIF!!!
               //        double correction    = (  V[i]-thresh[i]  )*(  1+(dt/tau_mem[i])*( (V_old[i]-Vreset[i] )/( V[i]-V_old[i] ) )   ) + Vreset[i];
               //        V[i]                 = correction;
                // Do  not use corrections as for the conductance it is irrelevant
                V[i]                 = Vreset[i];
                //Update efferents
                if (i<ne){ // E neuron fired spike //HOW DID IT LOOK BEFORE??? ISNT IT A BUG WITH h_E and h_I????
                    //save spike time
                    spk_t[i] = dt * (step - 1);
                    fprintf(fid_spikee, "%.2f %u\n", spk_t[i], i);
                    
                    for (int j=0; j<C[i].nb_conn_out_deg;j++) {
                        int tmp_conn_id = C[i].conn_id_out_deg[j];
                        gsyn[tmp_conn_id][i] += (C[i].conn_J_out_deg[j]);  //////// AAAAAAA
                    }
                }
                else{     //  I neuron fired spike
                    //save spike time
                    spk_t[i] = dt * (step - 1);
                    fprintf(fid_spikei, "%.2f %u\n", spk_t[i], i);
                    
                    for (int j=0; j<C[i].nb_conn_out_deg;j++) {
                        int tmp_conn_id = C[i].conn_id_out_deg[j];
                        gsyn[tmp_conn_id][i] += (C[i].conn_J_out_deg[j]);//////// AAAAAAA
                    }
                }
            }// end of if
        } // end of update efferents
        
        //write to memory
        if (save_voltage) {
            Ve_save[step][0]=dt*(step-1);
            Vi_save[step][0]=dt*(step-1);
            for (int l=1; l<=cells_to_save; l++){
                if ( (spk_t[l+(ne-cells_to_save)/2]==dt * (step - 1)) ){
                    Ve_save[step][l]=Vspike;
                }
                else{
                    Ve_save[step][l]=V[l+(ne-cells_to_save)/2];
                }
                if ( (spk_t[ne+l+(ni-cells_to_save)/2]==dt * (step - 1)) ){
                    Vi_save[step][l]=Vspike;
                }
                else{
                    Vi_save[step][l]=V[ne+l+(ni-cells_to_save)/2];
                }
            }
        }
        
        if (save_tot_current) {
            he_save[step][0]=dt*(step-1);
            hi_save[step][0]=dt*(step-1);
            for (int l=1; l<=cells_to_save; l++){
                he_save[step][l]=h_E[l]+Iext[l];
                hi_save[step][l]=h_I[l];
            }
        }

        if (save_current) {
            Cure_save[step][0]=dt*(step-1);
            Curi_save[step][0]=dt*(step-1);
            Cure_save[step][1]=h_E[1];
            Curi_save[step][1]=h_E[ne+1];
            Cure_save[step][2]=h_I[1];
            Curi_save[step][2]=h_I[ne+1];
            Cure_save[step][3]=Iext[1];
            Curi_save[step][3]=Iext[ne+1];
            Cure_save[step][4]=V[1];
            Curi_save[step][4]=V[ne+1];
        }
    }// end of steps loop
        
    if (save_voltage) {
        for (int step = 1; step < (steps_tot + 1); step++) {
            fprintf(fid_voltagee, "%.2f \t", Ve_save[step][0]);
            fprintf(fid_voltagei, "%.2f \t", Vi_save[step][0]);
            for(int i=1; i<  cells_to_save+1; i++){
                fprintf(fid_voltagee, "%.2f \t", Ve_save[step][i]);
                fprintf(fid_voltagei, "%.2f \t", Vi_save[step][i]);
            }
            fprintf(fid_voltagee, "\n");
            fprintf(fid_voltagei, "\n");
        }
    }
    

    if (save_tot_current) {
        for (int step = 1; step < (steps_tot + 1); step++) {
            fprintf(fid_he, "%.2f \t", he_save[step][0]);
            fprintf(fid_hi, "%.2f \t", hi_save[step][0]);
            for(int i=1; i<  cells_to_save+1; i++){
                fprintf(fid_he, "%.2f \t", he_save[step][i]);
                fprintf(fid_hi, "%.2f \t", hi_save[step][i]);
            }
            fprintf(fid_he, "\n");
            fprintf(fid_hi, "\n");
        }
    }
    
    
    if (save_current) {
        for (int step = 1; step < (steps_tot + 1); step++) {
            fprintf(fid_currente, "%.2f %.2f %.2f %.2f %.2f \n", Cure_save[step][0], Cure_save[step][1], Cure_save[step][2], Cure_save[step][3], Cure_save[step][4]);
            fprintf(fid_currenti, "%.2f %.2f %.2f %.2f %.2f \n", Curi_save[step][0], Curi_save[step][1], Curi_save[step][2], Curi_save[step][3], Curi_save[step][4]);
        }
    }
    
    
} //end of run

void Simulation::close() {
    fclose(fid_spikee);
    fclose(fid_spikei);
    if (save_voltage) {
        fclose(fid_voltagei);
        fclose(fid_voltagee);
    }
    if (save_tot_current) {
        fclose(fid_hi);
        fclose(fid_he);
    }
    if (save_current) {
        fclose(fid_currente);
        fclose(fid_currenti);
    }
}

//****************************************************************************

double *rk4_fun ( double t, int n, double u[], double C_m, double g_L, double V_L, double g_Na, double V_Na,
        double g_K, double V_K, double g_Adap, double V_Adap, double tau_Adap, double I_inp  )
        
//****************************************************************************
//
//  Purpose:
//
//    RK4VEC_FUN evaluates the right hand side of a vector ODE.
//
//  Parameters:
//
//    Input, double T, the current time.
//
//    Input, int N, the dimension of the system.
//
//    Input, double U[N], the current solution value.
//
//    Input, after U: parameters for the H&H equations
//
//    Output, double RK4VEC_FUN[N], the value of the derivative, dU/dT.
//
// U[0]=Voltage; U[1]=h; U[2]=n; U[3]=z;
{
    double *uprime;
    
    uprime = new double[n];
    double m_inf;
    double z_inf;
    double alpha_m;
    double beta_m;
    double alpha_n;
    double beta_n;
    double alpha_h;
    double beta_h;
    double I_L;
    double I_Na;
    double I_K;
    double I_Adap;
    double V, h_Na, n_K, z_Adap;
    
    
// Eq.1 for the voltage:
// CdV 		= -I_L-I_Na-I_K-I_Adap+I_inp
// I_L  	=  g_L*(V-V_L)
// I_Na 	=  g_Na*m_Na_inf^3*h_Na(V-V_Na)
// I_K		=  g_K*n_K^4*(V-V_K)
// I_Adap 	=  g_Adap*z_Adap*(V-V_Adap)
    
// Eq.2-3 for the gating param (m,n,h):
// dx=alpha_x*(1-x)-beta_x*x
    
// Eq.4 for the adapt param z:
// dz=(z_inf-z)/tau_adap
    
    V      = u[0];
    h_Na   = u[1];
    n_K    = u[2];
    z_Adap = u[3];
    
    
    alpha_m= (0.1*(V+30.0))/(1.0-exp(-0.1*(V+30.0)));
    beta_m= 4.0*exp(-(V+55.0)/18.0);
    
    m_inf=alpha_m/(alpha_m+beta_m);
    
    //alpha_h= 0.7*exp(-(V+58.0)/20.0);       David Type 2 old model
    //beta_h= 10.0/(exp(-0.1*(V+28.0))+1.0);
    
    alpha_h= 0.7*exp(-(V+44.0)/20.0);
    beta_h= 10.0/(exp(-0.1*(V+14.0))+1.0);
    
    alpha_n= (0.1*(V+34.0))/(1.0-exp(-0.1*(V+34.0)));
    beta_n= 1.25*exp(-(V+44.0)/80.0);
    
    z_inf=1.0/(1.0+exp(-0.7*(V+30.0)));
    
    I_L  		= g_L*(V-V_L);
    I_Na 		=  g_Na*pow(m_inf,3)*h_Na*(V-V_Na);
    I_K		=  g_K*pow(n_K,4)*(V-V_K);
    I_Adap 	=  g_Adap*z_Adap*(V-V_Adap);
    
    uprime[0] = (-I_L-I_Na-I_K-I_Adap+I_inp)/C_m;
    uprime[1] = alpha_h*(1-u[1])-beta_h*u[1];
    uprime[2] = alpha_n*(1-u[2])-beta_n*u[2];
    uprime[3] = (z_inf-u[3])/tau_Adap;
    
    return uprime;
}


