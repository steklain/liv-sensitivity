/* GLoBES -- General LOng Baseline Experiment Simulator
 * (C) 2002 - 2007,  The GLoBES Team
 *
 * GLoBES is mainly intended for academic purposes. Proper
 * credit must be given if you use GLoBES or parts of it. Please
 * read the section 'Credit' in the README file.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/* 
 * Example: Projection of two- and n-dimensional manifold onto stheta-axis
 * Compile with ``make example2''
 */ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <globes/globes.h>   /* GLoBES library */
//#include "myio.h"             /* my input-output routines */

#include <iostream>
#include <cmath>
#include<float.h>
#include<complex.h>
#include <vector>
#include<gsl/gsl_complex.h>
#include<gsl/gsl_complex_math.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<fstream>
#include <algorithm>

#define EXP_FAR  0
#define EXP_NEAR 1

#define MAX_SYS 200
double sys_errors[MAX_SYS];       /* Uncertainties of systematics params */
double sys_startval[MAX_SYS];     /* Starting values for systematics minimizer */
double sigma_binbin = 0.0;        /* Bin-to-bin error */

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

using namespace std;

array<double,4> chiscan(glb_params tv, double a)
{
  double res,th23,dcp,phi,min;
  double th23i = 39.7;
  double th23f = 50.9;
  double dcpi = 144;
  double dcpf = 350;

  array<double,4> params;

  int points = 100;

  //outparams<<"a"<<"  "<<"theta23"<<"  "<<"delta_CP"<<"  "<<"phi"<<endl;
  min=1000;
  for(th23=th23i;th23<th23f;th23=th23+(th23f-th23i)/points)
    for(dcp=dcpi;dcp<dcpf;dcp=dcp+(dcpf-dcpi)/points)
      for(phi=0;phi<2*M_PI;phi=phi+2*M_PI/points)
      {
        glbSetOscParams(tv,th23*M_PI/180.0,GLB_THETA_23);
        glbSetOscParams(tv,dcp*M_PI/180.0,GLB_DELTA_CP);
        glbSetOscParams(tv,phi,GLB_ARG_AEMU);

        /* Compute Chi^2 for all loaded experiments and all rules */
        res=glbChiSys(tv,GLB_ALL,GLB_ALL);
        if(res<min)
        {
            min=res;
            params = {th23,dcp,phi,min};
            //cout<<"chi = "<<min<<endl;
        } 
      }
  return params;

}

extern "C"
{
	#include "bsm.h"
}

/* If filename given, write to file; for empty filename write to screen */
char AEDLFILE[] = "DUNE_GLoBES.glb";
//char AEDLFILE[] = "T2HK.glb";
string OUTFILEchi = "BSMchi2_NP.dat";
string OUTFILEparams = "BSMchi2_params.dat";
string mylog = "BSMlog.dat";

int main(int argc, char *argv[])
{ 
  /* Initialize libglobes */
  glbInit(argv[0]); 

  /* Initialize experiment AEDFILE */
  glbInitExperiment(AEDLFILE,&glb_experiment_list[0],&glb_num_of_exps); 

  /* Intitialize output */
  ofstream outChi2np;
	outChi2np.open(OUTFILEchi);

  ofstream outparams;
	outparams.open(OUTFILEparams);

  ofstream outlog;
	outlog.open(mylog);

  /* Define standard oscillation parameters */
	double theta12  = 33.45*M_PI/180.0;// asin(sqrt(0.320));
  double theta13  = 8.62*M_PI/180.0;//asin(sqrt(0.02160));
	double theta23  = 42.1*M_PI/180.0;//asin(sqrt(0.547));
	double deltacp  = 230*M_PI/180.0;//-0.68 * M_PI;
  double dm21     = 7.42e-5;//7.55e-5;
	double dm31     = 2.51e-3;//2.50e-3;

  [[maybe_unused]] double theta12error = 0.14*theta12;
  [[maybe_unused]] double theta13error = 0.09*theta13;
	[[maybe_unused]] double theta23error = 0.27*theta23;
	[[maybe_unused]] double deltacperror = 1.0*deltacp;
  [[maybe_unused]] double dm21error    = 0.16*dm21;
	[[maybe_unused]] double dm31error    = 0.07*dm31;


  bsm_init_probability_engine_3();

	glbRegisterProbabilityEngine(8 * 9 - 3,
                             &bsm_probability_matrix,
                             &bsm_set_oscillation_parameters,
                             &bsm_get_oscillation_parameters,
                             NULL);
  

  /* Initialize parameter and projection vector(s) */
  glb_params true_values = glbAllocParams();
  glb_params test_values = glbAllocParams();
  glb_params input_errors = glbAllocParams();
  glb_params minimum = glbAllocParams();
  glb_projection theta13_projection = glbAllocProjection();
  glb_projection myprojection = glbAllocProjection();

  for(unsigned int i=0; i < 69; i++)
	{
	  glbSetOscParams(true_values, 0.0, i);
    glbSetOscParams(test_values, 0.0, i);
    glbSetOscParams(input_errors, 0.0, i);
	}

  glbDefineParams(true_values,theta12,theta13,theta23,deltacp,dm21,dm31);
  glbSetDensityParams(true_values,1.0,GLB_ALL);
  glbDefineParams(test_values,theta12,theta13,theta23,deltacp,dm21,dm31);  
  glbSetDensityParams(test_values,1.0,GLB_ALL);

  double abs_a_ee = 0;
	double abs_a_mue = 0;//2.0e-23* 1.0e9;
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
  glbSetOscParams(test_values, abs_a_etau, GLB_ABS_AETAU);  // a_etau 
  glbSetOscParams(test_values, arg_a_etau, GLB_ARG_AETAU);  // a_etau phase
  glbSetOscParams(test_values, abs_a_mumu, GLB_AMUMU);  // a_mumu
  glbSetOscParams(test_values, abs_a_mutau, GLB_ABS_AMUTAU);  // a_mutau
  glbSetOscParams(test_values, arg_a_mutau, GLB_ARG_AMUTAU);  // a_mutau phase
  glbSetOscParams(test_values, abs_a_tautau, GLB_ATAUTAU);  // a_tautau

	glbSetOscParams(test_values, abs_c_ee, GLB_CEE);  // c_ee 
	glbSetOscParams(test_values, abs_c_mue, GLB_ABS_CEMU);  // c_mue magnitude
  glbSetOscParams(test_values, arg_c_mue, GLB_ARG_CEMU);  // c_mue phase
  glbSetOscParams(test_values, abs_c_etau, GLB_ABS_CETAU);  // c_etau 
  glbSetOscParams(test_values, arg_c_etau, GLB_ARG_CETAU);  // c_etau phase
  glbSetOscParams(test_values, abs_c_mumu, GLB_CMUMU);  // c_mumu
  glbSetOscParams(test_values, abs_c_mutau, GLB_ABS_CMUTAU);  // c_mutau
  glbSetOscParams(test_values, arg_c_mutau, GLB_ARG_CMUTAU);  // c_mutau phase
  glbSetOscParams(test_values, abs_c_tautau, GLB_CTAUTAU);  // c_tautau
    

  /* The simulated data are computed */
  glbSetOscillationParameters(true_values);
  glbSetRates();

  /* Set starting values and input errors for all projections */  
  //glbDefineParams(input_errors,theta12error,theta13error,theta23error,deltacperror,dm21error,dm31error);

  //glbDefineParams(input_errors,0,0,theta23error,deltacperror,0,0);
  glbSetOscParams(input_errors,theta23error, GLB_THETA_23);
  glbSetOscParams(input_errors,deltacperror, GLB_DELTA_CP);
  glbSetOscParams(input_errors,M_PI, GLB_ARG_AEMU);

  glbSetDensityParams(input_errors,0.05,GLB_ALL);
  glbSetCentralValues(true_values);
  glbSetInputErrors(input_errors);

  /* Set two-parameter projection onto s22th13-axis: only deltacp free! 
  glbDefineProjection(theta13_projection,GLB_FIXED,GLB_FIXED,GLB_FIXED,
    GLB_FREE,GLB_FIXED,GLB_FIXED);
  for(unsigned int i=6; i < 69; i++)
	{
	glbSetProjectionFlag(theta13_projection, GLB_FIXED, i);
	}
  glbSetDensityProjectionFlag(theta13_projection, GLB_FIXED, GLB_ALL);
  glbSetProjection(theta13_projection); 

  */

  for(unsigned int i=0; i < 69; i++)
	{
	glbSetProjectionFlag(myprojection, GLB_FIXED, i);
	}

  //glbSetProjectionFlag(myprojection, GLB_FREE, GLB_THETA_13);
  glbSetProjectionFlag(myprojection, GLB_FREE, GLB_THETA_23);
  glbSetProjectionFlag(myprojection, GLB_FREE, GLB_DELTA_CP);
  glbSetProjectionFlag(myprojection, GLB_FREE, GLB_ARG_AEMU);

  glbSetDensityProjectionFlag(myprojection,GLB_FIXED,GLB_ALL);
  glbSetProjection(myprojection);

  /* Iteration over all values to be computed */
  double a,res,minth23,mindcp,minargamue;
  [[maybe_unused]] double mymin,mytest,mytrue,restest,resmin;
  double ai = 0.7;
  double af = 5.0;
  int points = 200;
  array<double,4> myparams;
  int counter = 1;

  

  outparams<<"a"<<"  "<<"theta23"<<"  "<<"delta_CP"<<"  "<<"phi"<<endl;
  outlog<<"This is the log"<<endl;

  for(a=ai;a<=af;a=a+(af-ai)/points)
  {
      /* Set vector of test=fit values */

      glbSetOscParams(test_values,a*1.0e-23*1.0e9,GLB_ABS_AEMU);

      if(counter == points/20)
      {
        outlog<<"Starting scan... a = "<<a<<endl;
        myparams=chiscan(test_values,a);
        outlog<<"Minimum found: "<<myparams[3]<<endl;
        outlog<<"Parameters: theta23 = "<<myparams[0]<<" deltacp = "<<myparams[1]<<" phi = "<<myparams[2]<<endl;
        glbSetOscParams(test_values,myparams[0]*M_PI/180,GLB_THETA_23);
        glbSetOscParams(test_values,myparams[1]*M_PI/180,GLB_DELTA_CP);
        glbSetOscParams(test_values,myparams[2],GLB_ARG_AEMU);
        counter = 0;
      }



      /* Guess fit value for deltacp in order to safely find minimum */
       
      res=glbChiNP(test_values,minimum,GLB_ALL);

      /* Compute Chi^2 for all loaded experiments and all rules */
      //res=glbChiNP(test_values,minimum,GLB_ALL);
      //restest=glbChiSys(test_values,GLB_ALL,GLB_ALL);
      //resmin=glbChiSys(minimum,GLB_ALL,GLB_ALL);
      //res=glbChiAll(test_values,minimum,GLB_ALL);

      minth23 = glbGetOscParams(minimum,GLB_THETA_23);
      mindcp =  glbGetOscParams(minimum,GLB_DELTA_CP);
      minargamue = glbGetOscParams(minimum,GLB_ARG_AEMU);

      //mymin = glbGetOscParams(minimum,GLB_ABS_AEMU)/(1.0e-23*1.0e9);
      //mytest = glbGetOscParams(test_values,GLB_ABS_AEMU)/(1.0e-23*1.0e9);
      //mytrue = glbGetOscParams(true_values,GLB_ABS_AEMU)/(1.0e-23*1.0e9);

      outChi2np<<a<<"  "<<res<<endl;
      outparams<<a<<"  "<<minth23*180/M_PI<<"  "<<mindcp*180/M_PI<<"  "<<minargamue*180/M_PI<<endl;

      counter++;
  }
  
  /* Destroy parameter and projection vector(s) */
  glbFreeParams(true_values);
  glbFreeParams(test_values); 
  glbFreeParams(input_errors); 
  glbFreeProjection(theta13_projection);
  glbFreeProjection(myprojection);
    
  exit(0);
}
