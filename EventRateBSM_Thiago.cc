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
#include <iomanip>

#include <algorithm>

extern "C"
{
	#include "bsm.h"
}

using namespace std;

char AEDLFILE[] 	= "./lib/DUNE_GLoBES.glb";

string OUTFILESTDe	= "eventrate_e_BSM.dat";
string OUTFILESTDmu = "eventrate_mu_BSM.dat";
FILE * outstde 	= NULL;
FILE * outstdmu = NULL;
/*
string OUTFILEBSMe	= "eventrate_e_BSM.dat";
string OUTFILEBSMmu = "eventrate_mu_BSM.dat";
FILE * outbsme 	= NULL;
FILE * outbsmmu = NULL;
*/

//string OUTBSM_e = "eventvalues_e.dat";
//string OUTBSM_mu = "eventvalues_mu.dat";
string OUTFILE = "eventfilter_e.dat";
string OUTFILMU = "eventfilter_mu.dat";

int main(int argc, char * argv[])
{

	glbInit(argv[0]);
	glbInitExperiment(AEDLFILE, &glb_experiment_list[0], &glb_num_of_exps);

	outstde = fopen(OUTFILESTDe.c_str(), "w");
	if (outstde == NULL) 
    {
        printf("Error opening output file.\n");
		return -1;
	}

	outstdmu = fopen(OUTFILESTDmu.c_str(), "w");
	if (outstdmu == NULL) 
    {
        printf("Error opening output file.\n");
		return -1;
	}

	
	ofstream outfilter_e;
	ofstream outfilter_mu;

	outfilter_e.open(OUTFILE);
	outfilter_mu.open(OUTFILMU);

	double dm21 = 7.41e-5;//double theta12  = 33.45*M_PI/180.0;// asin(sqrt(0.320));
  	double dm31 = 2.498e-3;//double theta13  = 8.62*M_PI/180.0;//asin(sqrt(0.02160));
	double theta12 = 33.45*M_PI/180;//double theta23  = 42.1*M_PI/180.0;//asin(sqrt(0.547));
	double theta23 = 42.1*M_PI/180;//double deltacp  = 230*M_PI/180.0;//-0.68 * M_PI;
  	double theta13 = 8.62;//double dm21     = 7.42e-5;//7.55e-5;
	double deltacp = 230*M_PI/180;//double dm31     = 2.51e-3;//2.50e-3;

	int n_bins = glbGetNumberOfBins(0);

    bsm_init_probability_engine_3();

	glbRegisterProbabilityEngine(8 * 9 - 3,
                               &bsm_probability_matrix,
							   &bsm_set_oscillation_parameters,
  							   &bsm_get_oscillation_parameters,
  							   NULL);

	/* Define "true" oscillation parameter vector */
	glb_params true_values = glbAllocParams();

    for(unsigned int i=0; i < 69; i++)
	{
	glbSetOscParams(true_values, 0.0, i);
	}
	

    glbDefineParams(true_values,theta12,theta13,theta23,deltacp,dm21,dm31);
    glbSetDensityParams(true_values, 1.0, GLB_ALL);
	glbSetOscillationParameters(true_values);
	glbSetRates();

    int chanosce = glbNameToValue(0, "channel", "#FHC_app_osc_nue");
	int chanoscmu = glbNameToValue(0, "channel", "#FHC_dis_sig_numu");
	
	double *myratesstde = glbGetChannelRatePtr(0, chanosce, GLB_PRE);
    glbShowChannelRates(outstde, 0, chanosce, GLB_PRE, GLB_W_EFF, GLB_W_BG);
	double *myratesstdmu = glbGetChannelRatePtr(0, chanoscmu, GLB_PRE);
    glbShowChannelRates(outstdmu, 0, chanoscmu, GLB_PRE, GLB_W_EFF, GLB_W_BG);

	double stdratese[n_bins];
	double stdratesmu[n_bins];

	for(int i=0;i<n_bins;i++){
		stdratese[i]	=	myratesstde[i];
		stdratesmu[i]	=	myratesstdmu[i];
	}

	printf("%g\n", stdratese[n_bins-1]);
	printf("%g\n", stdratesmu[n_bins-1]);


	double a,b;
	double bi = -25;
	double bf = -23;
	double ai = 0;
	double af = 10;
	double points = 500;
	double y;
	double x;
	
	
	glb_params test_values = glbAllocParams();
    for(unsigned int i=0; i < 69; i++)
	{
	glbSetOscParams(test_values, 0.0, i);
	}
	
    glbDefineParams(test_values,theta12,theta13,theta23,deltacp,dm21,dm31);
    glbSetDensityParams(test_values, 1.0, GLB_ALL);

	for(b=bi; b<bf; b++)
		for(a=ai;a<af;a=a+(af-ai)/points){


			glbSetOscParams(test_values,a*pow(10,b)*1.0e9, 61);
			glbSetOscillationParameters(test_values);
			glbSetRates();

			double *myratesbsme =  glbGetChannelRatePtr(0, chanosce, GLB_PRE);
			double *myratesbsmmu = glbGetChannelRatePtr(0, chanoscmu, GLB_PRE);

			for(int i=0; i<n_bins;i++){
				y = ((myratesbsme[i]*100) /stdratese[i])-100;
				x = ((myratesbsmmu[i]*100)/stdratesmu[i])-100;

				if(y > 10){
					outfilter_e<< a <<"  "<< b <<"  "<< y << endl;
				}

				if(x < -10){
					outfilter_mu<< a <<"  "<< b <<"  "<< x << endl;
				}
			}

}

	fclose(outstde);
	fclose(outstdmu);
	//fclose(outbsme);
	//fclose(outbsmmu);

	outfilter_e.close();
	outfilter_mu.close();
	glbFreeParams(true_values);
 	return 0;

}


