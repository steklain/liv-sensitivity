#include <iostream>
#include <cmath>
#include <string.h>
#include <float.h>
#include <complex.h>
#include <vector>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <globes/globes.h>
#include <fstream>
#include <string>
#include <algorithm>
#include <sys/stat.h>
#include <filesystem>
#define GLB_SIGMA_E 6

const int GLB_AEE = 51;
const int GLB_ABS_AEMU = 52;
const int GLB_ARG_AEMU = 53;
const int GLB_ABS_AETAU = 54;
const int GLB_ARG_AETAU = 55;
const int GLB_AMUMU = 56;
const int GLB_ABS_AMUTAU = 57;
const int GLB_ARG_AMUTAU = 58;
const int GLB_ATAUTAU = 59;

const int GLB_CEE = 60;
const int GLB_ABS_CEMU = 61;
const int GLB_ARG_CEMU = 62;
const int GLB_ABS_CETAU = 63;
const int GLB_ARG_CETAU = 64;
const int GLB_CMUMU = 65;
const int GLB_ABS_CMUTAU = 66;
const int GLB_ARG_CMUTAU = 67;
const int GLB_CTAUTAU = 68;

extern "C"
{
    #include "bsm.h"
}
using namespace std;

char AEDLFILE[] = "DUNE_GLoBES.glb";
int main(int argc, char * argv[])
{
    glbInit(argv[0]);
    glbInitExperiment(AEDLFILE, &glb_experiment_list[0], &glb_num_of_exps);
    
    ofstream outBackup;
    ofstream outBSM;
    ofstream outBSMn;
    ofstream outBSM_co;
    ofstream outchi;

    /*Arquivo de entrada chi2*/
    outchi.open("chi2_co.dat");
    outBSM.open("betachi2.dat");
    outBSM_co.open("parametros.dat");
    outBSMn.open("BSMprobability_a_mue.dat");
    
    
    /*nome base para criação de todos arquivos e pastas
    o restante nome_novo etc. para os especificos*/
    ifstream in;

    double dm21 = 7.42e-5;
    double dm31 = 2.515e-3;
    double theta12 = asin(sqrt(0.304));
    double theta23 = asin(sqrt(0.573));
    double theta13 = asin(sqrt(0.0222));
    double deltacp = -0.68 * M_PI; 

    /*ERRORS*/

    double theta12error = 0.14*theta12;
    double theta13error = 0.09*theta13;
	double theta23error = 0.27*theta23;
	double deltacperror = 1.0*deltacp;
    double dm21error    = 0.16*dm21;
	double dm31error    = 0.07*dm31;   

    bsm_init_probability_engine_3();
    glbRegisterProbabilityEngine(8 * 9 - 3,
                               &bsm_probability_matrix,
                               &bsm_set_oscillation_parameters,
                               &bsm_get_oscillation_parameters,
                               NULL);

    /* Define "true" oscillation parameter vector */
    glb_params true_values = glbAllocParams();
    glb_params test_values = glbAllocParams();
    glb_params input_errors = glbAllocParams();
    glb_projection theta13_projection = glbAllocProjection();
    glb_projection myprojection = glbAllocProjection();
    
    for(unsigned int i=0; i < 69; i++)
    {
    glbSetOscParams(true_values, 0.0, i);
    glbSetOscParams(test_values, 0.0, i);
    glbSetOscParams(input_errors, 0.0, i);
    }

    for(unsigned int i=0; i < 69; i++)
	{
	glbSetProjectionFlag(myprojection, GLB_FIXED, i);
	}

    glbSetProjectionFlag(myprojection, GLB_FREE, GLB_THETA_13);
    //glbSetProjectionFlag(myprojection, GLB_FREE, GLB_THETA_23);
    glbSetProjectionFlag(myprojection, GLB_FREE, GLB_DELTA_CP);
    glbSetProjectionFlag(myprojection, GLB_FREE, GLB_ARG_AEMU);

    glbSetDensityProjectionFlag(myprojection,GLB_FIXED,GLB_ALL);
    glbSetProjection(myprojection);


    glbDefineParams(true_values,theta12,theta13,theta23,deltacp,dm21,dm31);
    glbDefineParams(test_values,theta12,theta13,theta23,deltacp,dm21,dm31);
   

    /*define test_values for chi²*/
    glbSetRates();

    /*Termos para Matriz a*/
    double abs_a_ee = 0;
    double abs_a_mue = -2.0e-23 / 1.0e-9;
    double arg_a_mue =  0;
    double abs_a_etau = 0;
    double arg_a_etau = 0;
    double abs_a_mumu = 0;
    double abs_a_mutau = 0;
    double arg_a_mutau = 0;
    double abs_a_tautau = 0;

    /*Termos para Matriz c*/
    double abs_c_ee = 0;
    double abs_c_mue = 0;//-2.0e-32 / 1.0e-9;
    double arg_c_mue = 0;
    double abs_c_etau = 0;
    double arg_c_etau = 0;
    double abs_c_mumu = 0;
    double abs_c_mutau = 0;
    double arg_c_mutau = 0;
    double abs_c_tautau = 0;


    //############ LIV Parameter #################################//

    /*A*/
    glbSetOscParams(test_values, abs_a_ee, GLB_AEE);  // a_ee
    glbSetOscParams(test_values, abs_a_mue, GLB_ABS_AEMU);  // a_mue magnitude
    glbSetOscParams(test_values, arg_a_mue, GLB_ARG_AEMU);  // a_mue phase
    glbSetOscParams(test_values, abs_a_etau, GLB_ABS_AETAU);  // a_etau
    glbSetOscParams(test_values, arg_a_etau, GLB_ARG_AETAU);  // a_etau phase
    glbSetOscParams(test_values, abs_a_mumu, GLB_AMUMU);  // a_mumu
    glbSetOscParams(test_values, abs_a_mutau, GLB_ABS_AMUTAU);  // a_mutau
    glbSetOscParams(test_values, arg_a_mutau, GLB_ARG_AMUTAU);  // a_mutau phase
    glbSetOscParams(test_values, abs_a_tautau, GLB_ATAUTAU);  //a_tautau

    /*C*/
    glbSetOscParams(test_values, abs_c_ee, GLB_CEE);  // c_ee
    glbSetOscParams(test_values, abs_c_mue, GLB_ABS_CEMU);  // c_mue magnitude
    glbSetOscParams(test_values, arg_c_mue, GLB_ARG_CEMU);  // c_mue phase
    glbSetOscParams(test_values, abs_c_etau, GLB_ABS_CETAU);  // c_etau
    glbSetOscParams(test_values, arg_c_etau, GLB_ARG_CETAU);  // c_etau phase
    glbSetOscParams(test_values, abs_c_mumu, GLB_CMUMU);  // c_mumu
    glbSetOscParams(test_values, abs_c_mutau, GLB_ABS_CMUTAU);  // c_mutau
    glbSetOscParams(test_values, arg_c_mutau, GLB_ARG_CMUTAU);  // c_mutau phase
    glbSetOscParams(test_values, abs_c_tautau, GLB_CTAUTAU);// c_tautest

    glbSetDensityParams(true_values,1.0,GLB_ALL);
    glbSetDensityParams(test_values,1.0,GLB_ALL);
    glbSetDensityParams(input_errors, 0.05, GLB_ALL);

    glbSetInputErrors(input_errors);
    
    glbSetOscillationParameters(true_values);
    glbSetOscillationParameters(test_values);
    
    glbCopyParams(true_values,test_values);

    /*Correlações*/
	glbDefineParams(input_errors,theta12error,theta13error,theta23error,deltacperror,dm21error,dm31error);
	glbSetDensityParams(input_errors,0.05,GLB_ALL);
	glbSetCentralValues(true_values);
	glbSetInputErrors(input_errors);

    /* Define my own two-parameter projection for glbChiNP: Only deltacp is free! */
	glbDefineProjection(theta13_projection,GLB_FIXED,GLB_FIXED,GLB_FIXED,
    GLB_FREE,GLB_FIXED,GLB_FIXED);

      for(unsigned int i=6; i < 69; i++)
	{
	glbSetProjectionFlag(theta13_projection, GLB_FIXED, i);
	}

    glbSetDensityProjectionFlag(theta13_projection,GLB_FIXED,GLB_ALL);
	glbSetProjection(theta13_projection);

    for(unsigned int i=0; i < 69; i++)
	{
	glbSetProjectionFlag(myprojection, GLB_FIXED, i);
	}

    glbSetProjectionFlag(myprojection, GLB_FREE, GLB_THETA_13);
    //glbSetProjectionFlag(myprojection, GLB_FREE, GLB_THETA_23);
    glbSetProjectionFlag(myprojection, GLB_FREE, GLB_DELTA_CP);
    glbSetProjectionFlag(myprojection, GLB_FREE, GLB_ARG_AEMU);	

    glbSetRates();
    
    double energy, prob_BSM;
    double emin= 0.25 ; //GeV
    double emax=8 ; //GeV
    double step= 1000;
    //double L = 1300; // km
    for (energy=emin;energy<=emax;energy+=(emax-emin)/step)
    {
    glbSetOscillationParameters(true_values);
    //prob_BSM=glbVacuumProbability(2,1,+1,energy,L);
    prob_BSM=glbProfileProbability(0,2,1,+1,energy);
    outBSMn<<energy<<"  "<<prob_BSM<<endl;
    }
    
    /* Compute chi^2 sem as correlações */
    double a, x, y, z, chi;
    int i,j;
    double ai = 0;
    double af = 5.0;
    int points = 200;

    for(i=1; i<2; i++)        /* th13 loop */
    for(j=1; j<2; j++)  /*loop */
    for(a=ai;a<af;a=a+(af-ai)/points)
    {
        /*É preciso definir o parametro de espaço, utilizando 100 steps
        o alcance do theta13 é 0.122-0.175 e do deltacp é 0-2pi*/

        x=0.13+i*(0.17-0.13)/100;
        y=0+j*(2*M_PI)/100;
		z= a*abs_a_mue;

	    /* Seta o vetor dos valores que vão filtrar os test_values */
	    glbSetOscParams(test_values,x,GLB_THETA_13);
	    glbSetOscParams(test_values,y,GLB_DELTA_CP);
		glbSetOscParams(test_values,z,GLB_ABS_AEMU);
	
	    /* Compute Chi^2 com os sistematicos apenas, para todos os experimentos e regras*/
	    chi=glbChiNP(test_values,NULL,GLB_ALL);
		//chi2=glbChiTheta13(test_values,NULL,GLB_ALL);
        outchi<<a<<" "<<chi<<endl;
	    outBSM_co<<x<<" "<<y<<" "<<z<<endl;
    }
    outchi.close();
    outBSM_co.close();
    outBSMn.close();
    outBackup.close();
    outBSM.close();
    glbFreeParams(true_values);
    glbFreeParams(test_values);
    glbFreeParams(input_errors); 
    glbFreeProjection(theta13_projection);
    glbFreeProjection(myprojection);
    return 0;
}