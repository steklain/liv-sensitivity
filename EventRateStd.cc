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

using namespace std;

char AEDLFILE[] = "DUNE_GLoBES.glb";
//char AEDLFILE[] = "T2HK.glb";
string OUTFILEe = "eventrate_e.dat";
string OUTFILEmu = "eventrate_mu.dat";
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

	double theta12  = 33.45*M_PI/180.0;// asin(sqrt(0.320));
  	double theta13  = 8.62*M_PI/180.0;//asin(sqrt(0.02160));
	double theta23  = 42.1*M_PI/180.0;//asin(sqrt(0.547));
	double deltacp  = 230*M_PI/180.0;//-0.68 * M_PI;
  	double dm21     = 7.42e-5;//7.55e-5;
	double dm31     = 2.51e-3;//2.50e-3;

	/* Define "true" oscillation parameter vector */
	glb_params true_values = glbAllocParams();
    glbDefineParams(true_values,theta12,theta13,theta23,deltacp,dm21,dm31);
	glbSetDensityParams(true_values, 1.0, GLB_ALL);
	glbSetOscillationParameters(true_values);
	glbSetRates();

	int glbNameToValue(int exp, const char* context, const char *name);
    int glbShowRuleRates(FILE* stream, int exp, int rule, int pos, int effi, int bgi, int coeffi, int signal);

    //int ruleosc = glbNameToValue(0, "rule", "#Nu_Mu_Appearance");
    int ruleosce = glbNameToValue(0, "rule", "#nue_app");
    printf("%d\n",ruleosce);
    glbShowRuleRates(outstde, 0, ruleosce, GLB_ALL, GLB_W_EFF, GLB_W_BG, GLB_W_COEFF, GLB_SIG);

	int ruleoscmu = glbNameToValue(0, "rule", "#numu_dis");
    printf("%d\n",ruleoscmu);
    glbShowRuleRates(outstdmu, 0, ruleoscmu, GLB_ALL, GLB_W_EFF, GLB_W_BG, GLB_W_COEFF, GLB_SIG);


    //const char *glbValueToName(int exp,const char* context, int value); 
    //const char* name = glbValueToName(0, "rule", 2);
    //printf("%s\n",name);

	fclose(outstde);
	fclose(outstdmu);
	
	glbFreeParams(true_values);
 	return 0;

}