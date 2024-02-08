#include <iostream>
#include <cmath>
#include <string>
#include<float.h>
#include<complex.h>
#include <vector>
#include<gsl/gsl_complex.h>
#include<gsl/gsl_complex_math.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include <globes/globes.h>
#include<fstream>

#include <algorithm>

const int GLB_AEE = 51;
const int GLB_ABS_AEMU = 52;
const int GLB_ARG_AEMU = 53;

extern "C"
{
	#include "bsm.h"
}

using namespace std;

char AEDLFILE[] = "DUNE_GLoBES.glb";
string OUTFILEe = "DUMMYeventrate_e_BSM.dat";
string OUTFILEmu = "DUMMYeventrate_mu_BSM.dat";
FILE * outstde = NULL;
FILE * outstdmu = NULL;

int main(int argc, char * argv[])
{

	glbInit(argv[0]);
	glbInitExperiment(AEDLFILE, &glb_experiment_list[0], &glb_num_of_exps);

	outstde = fopen(OUTFILEe.c_str(), "w");
	if (outstde == NULL) 
    {
        printf("Error opening output file.\n");
		return -1;
	}

	outstdmu = fopen(OUTFILEmu.c_str(), "w");
	if (outstdmu == NULL) 
    {
        printf("Error opening output file.\n");
		return -1;
	}

    ofstream outChi2stat;
	outChi2stat.open("BSMchi2_statisticsonly.dat");

    ofstream outChi2sys;
	outChi2sys.open("BSMchi2_systematics.dat");

    ofstream outChi2np;
	outChi2np.open("BSMchi2_NP.dat");

	double dm21 = 7.42e-5;//7.55e-5;
	double dm31 = 2.51e-3;//2.50e-3;
	double theta12 = 33.45*M_PI/180;// asin(sqrt(0.320));
	double theta23 = 42.1*M_PI/180;//asin(sqrt(0.547));
	double theta13 = 8.62*M_PI/180;//asin(sqrt(0.02160));
	double deltacp = 230*M_PI/180;//-0.68 * M_PI;

    bsm_init_probability_engine_3();

	glbRegisterProbabilityEngine(8 * 9 - 3,
                               &bsm_probability_matrix,
							   &bsm_set_oscillation_parameters,
  							   &bsm_get_oscillation_parameters,
  							   NULL);

	/* Initialize patameter and projection vectors */
	glb_params true_values = glbAllocParams();
    glb_params test_values = glbAllocParams();
    glb_params input_errors = glbAllocParams();
    glb_projection myprojection = glbAllocProjection();


    for(unsigned int i=0; i < 69; i++)
	{
	glbSetOscParams(true_values, 0.0, i);
    glbSetOscParams(test_values, 0.0, i);
    glbSetOscParams(input_errors, 0.0, i);
	}

    glbDefineParams(true_values,theta12,theta13,theta23,deltacp,dm21,dm31);
    glbDefineParams(test_values,theta12,theta13,theta23,deltacp,dm21,dm31);
	
    double abs_a_ee = 0;
	double abs_a_mue = 2.0e-23* 1.0e9;
	double arg_a_mue = 0;
	double abs_a_etau = 0;
	double arg_a_etau = 0;
	double abs_a_mumu = 0;
	double abs_a_mutau = 0;
	double arg_a_mutau = 0;
	double abs_a_tautau = 0;

	double abs_c_ee = 0;
	double abs_c_mue = 0;//-2.0e-25 * 1.0e9;
	double arg_c_mue = 0;
	double abs_c_etau = 0;
	double arg_c_etau = 0;
	double abs_c_mumu = 0;//1.0e-22;
	double abs_c_mutau = 0;
	double arg_c_mutau = 0;
	double abs_c_tautau = 0;
    
    //############ LIV Parameters #################################//
	glbSetOscParams(test_values, abs_a_ee, GLB_AEE);  // a_ee 
	glbSetOscParams(test_values, abs_a_mue, GLB_ABS_AEMU);  // a_mue magnitude
    glbSetOscParams(test_values, arg_a_mue, GLB_ARG_AEMU);  // a_mue phase
    glbSetOscParams(test_values, abs_a_etau, 54);  // a_etau 
    glbSetOscParams(test_values, arg_a_etau, 55);  // a_etau phase
    glbSetOscParams(test_values, abs_a_mumu, 56);  // a_mumu
    glbSetOscParams(test_values, abs_a_mutau, 57);  // a_mutau
    glbSetOscParams(test_values, arg_a_mutau, 58);  // a_mutau phase
    glbSetOscParams(test_values, abs_a_tautau, 59);  // a_tautau

	glbSetOscParams(test_values, abs_c_ee, 60);  // c_ee 
	glbSetOscParams(test_values, abs_c_mue, 61);  // c_mue magnitude
    glbSetOscParams(test_values, arg_c_mue, 62);  // c_mue phase
    glbSetOscParams(test_values, abs_c_etau, 63);  // c_etau 
    glbSetOscParams(test_values, arg_c_etau, 64);  // c_etau phase
    glbSetOscParams(test_values, abs_c_mumu, 65);  // c_mumu
    glbSetOscParams(test_values, abs_c_mutau, 66);  // c_mutau
    glbSetOscParams(test_values, arg_c_mutau, 67);  // c_mutau phase
    glbSetOscParams(test_values, abs_c_tautau, 68);  // c_tautau
    
    glbSetDensityParams(true_values, 1.0, GLB_ALL);
    glbSetDensityParams(test_values, 1.0, GLB_ALL);

	glbSetOscillationParameters(true_values);
	glbSetRates();

	int glbNameToValue(int exp, const char* context, const char *name);
    int glbShowRuleRates(FILE* stream, int exp, int rule, int pos, int effi, int bgi, int coeffi, int signal);

    int ruleosce = glbNameToValue(0, "rule", "#nue_app");
    //printf("%d\n",ruleosce);
    glbShowRuleRates(outstde, 0, ruleosce, GLB_ALL, GLB_W_EFF, GLB_W_BG, GLB_W_COEFF, GLB_SIG);

	int ruleoscmu = glbNameToValue(0, "rule", "#numu_dis");
    //printf("%d\n",ruleoscmu);
    glbShowRuleRates(outstdmu, 0, ruleoscmu, GLB_ALL, GLB_W_EFF, GLB_W_BG, GLB_W_COEFF, GLB_SIG);


    glbDefineParams(input_errors,theta12*0.1,0,0,0,dm21*0.1,0);
    //glbSetOscParams(input_errors, M_PI, GLB_ARG_AEMU);
    glbSetDensityParams(input_errors,0.05,GLB_ALL);
    glbSetCentralValues(true_values);
    glbSetInputErrors(input_errors);


    double a,res;
    double ai = 0;
    double af = 5.0;

    //glbDefineProjection(myprojection,GLB_FIXED, GLB_FIXED, GLB_FIXED,GLB_FIXED, GLB_FIXED, GLB_FIXED);
    
    for(unsigned int i=0; i < 69; i++)
	{
	glbSetProjectionFlag(myprojection, GLB_FIXED, i);
	}

    glbSetProjectionFlag(myprojection, GLB_FREE, GLB_THETA_23);
    glbSetProjectionFlag(myprojection, GLB_FREE, GLB_DELTA_CP);
    glbSetProjectionFlag(myprojection, GLB_FREE, GLB_ARG_AEMU);
    
    glbSetDensityProjectionFlag(myprojection,GLB_FREE,GLB_ALL);
    glbSetProjection(myprojection);
    //res = glbChiNP(test_values,minimum,GLB_ALL);
    //printf("chi2 with correlation only with a_mue:%g \n\n",res);
    
    for(a=ai;a<af;a=a+(af-ai)/100)
      {
       /* Set vector of test values */
       glbSetOscParams(test_values, a*1.0e-23*1.0e9, 52);
       glbSetCentralValues(test_values);
  
       /* Compute Chi^2 for all loaded experiments and all rules */
       res=glbChiNP(test_values,NULL,GLB_ALL);
       outChi2np<<a<<"  "<<res<<endl;
      }
        
    glbFreeProjection(myprojection);


    for(a=ai;a<af;a=a+(af-ai)/100)
      {
       /* Set vector of test values */
       glbSetOscParams(test_values, a*1.0e-23*1.0e9, 52);
  
       /* Compute Chi^2 for all loaded experiments and all rules */
       res=glbChiSys(test_values,GLB_ALL,GLB_ALL);
       outChi2sys<<a<<"  "<<res<<endl;
      }
     
    glbSwitchSystematics(GLB_ALL,GLB_ALL,GLB_OFF);

    for(a=ai;a<af;a=a+(af-ai)/100)
      {
       /* Set vector of test values */
       glbSetOscParams(test_values, a*1.0e-23*1.0e9, 52);
  
       /* Compute Chi^2 for all loaded experiments and all rules */
       res=glbChiSys(test_values,GLB_ALL,GLB_ALL);
       outChi2stat<<a<<"  "<<res<<endl;
      }
    
    glbSwitchSystematics(GLB_ALL,GLB_ALL,GLB_ON);




    //const char *glbValueToName(int exp,const char* context, int value); 
    //const char* name = glbValueToName(0, "rule", 2);
    //printf("%s\n",name);

	fclose(outstde);
	fclose(outstdmu);

    outChi2stat.close();
    outChi2sys.close();
    outChi2np.close();
	
	glbFreeParams(true_values);
    glbFreeParams(test_values); 
    glbFreeProjection(myprojection);
 	return 0;

}