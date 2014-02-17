#include "mex.h"
#include<math.h>

/* APF 2010-03-10 */

double logfactorial(double i)
{
    double logfact = 0;
    for(i;i>1;i--)
        logfact += log(i);
    return logfact;
}

/* ARGUMENTS:
 * double photonrate0dt
 * double photonrate1dt
 * double k01dt
 * double k10dt
 * double *photonbins (N elements)
 * 
 * RETURNS:
 * double loglikelihood
 * mxLogical *states (N elements)
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double photonrate0dt = mxGetScalar(prhs[0]);
    double photonrate1dt = mxGetScalar(prhs[1]);
    double k01dt = mxGetScalar(prhs[2]);
    double k10dt = mxGetScalar(prhs[3]);
    double *photonbins = mxGetPr(prhs[4]);
    double log0 = log(photonrate0dt);
    double log1 = log(photonrate1dt);
    double log00 = log(1-k01dt);
    double log01 = log(k01dt);
    double log10 = log(k10dt);
    double log11 = log(1-k10dt);
    const mwSize numbins = mxGetM(prhs[4]);
    double logfact = logfactorial(photonbins[0]);
    double delta0 = photonbins[0]*log0-photonrate0dt-logfact;
    double delta1 = photonbins[0]*log1-photonrate1dt-logfact;
    double delta00,delta01,delta10,delta11;
    mxLogical *states;
    mxLogical *backtracker = mxCalloc(2*numbins,sizeof(mxLogical));
    mwIndex i;
    
    plhs[1] = mxCreateLogicalMatrix(numbins,1); /*states matrix*/
    
    states = mxGetLogicals(plhs[1]);
    
    if(states == NULL)
        mexErrMsgTxt("Not enough memory to create states vector!");

    backtracker[0] = 0;
    backtracker[1] = 0;
    for(i = 1; i < numbins; i++)
    {
        delta00 = delta0+log00;
        delta01 = delta0+log01;
        delta10 = delta1+log10;
        delta11 = delta1+log11;
        logfact = logfactorial(photonbins[i]);
        if(delta00 >= delta10)
        {
            delta0 = delta00+photonbins[i]*log0-photonrate0dt-logfact;
            backtracker[2*i] = 0;
        }
        else
        {
            delta0 = delta10+photonbins[i]*log0-photonrate0dt-logfact;
            backtracker[2*i] = 1;
        }
        if(delta11 >= delta01)
        {
            delta1 = delta11+photonbins[i]*log1-photonrate1dt-logfact;
            backtracker[2*i+1] = 1;
        }
        else
        {
            delta1 = delta01+photonbins[i]*log1-photonrate1dt-logfact;
            backtracker[2*i+1] = 0;
        }
    }
    if(delta0 > delta1)
    {
        plhs[0] = mxCreateDoubleScalar(delta0);
        states[numbins-1] = 0;
    }
    else
    {
        plhs[0] = mxCreateDoubleScalar(delta1);
        states[numbins-1] = 1;
    }
    for(i = numbins-1;i >= 1; i--)
        states[i-1] = backtracker[2*i+states[i]];
    mxFree(backtracker);
}
