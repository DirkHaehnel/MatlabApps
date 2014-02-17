#include "mex.h"
#include<math.h>

/* APF 2010-03-10 */

/* ARGUMENTS:
 * double photonrate0dt
 * double photonrate1dt
 * double k01dt
 * double k10dt
 * double *photonbins (N elements)
 * 
 * RETURNS:
 * double newphotonrate0dt
 * double newphotonrate1dt
 * double newk01dt
 * double newk10dt
 * double *odds (N elements)
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double photonrate0dt = mxGetScalar(prhs[0]);
    double photonrate1dt = mxGetScalar(prhs[1]);
    double k01dt = mxGetScalar(prhs[2]);
    double k10dt = mxGetScalar(prhs[3]);
    double *photonbins = mxGetPr(prhs[4]);
    const mwSize numbins = mxGetM(prhs[4]);
    double photonodds = photonrate0dt/photonrate1dt;
    double expphotonratediff = exp(photonrate1dt-photonrate0dt);
    double *odds;
    double backwards;
    double currtranscounts[4] = {0,0,0,0};
    double sumcurrtranscounts;
    double transcounts[4] = {0,0,0,0};
    double photoncounts[2] = {0,0};
    mwIndex i;

    plhs[4] = mxCreateDoubleMatrix(numbins,1,mxREAL); /*odds matrix*/
    
    odds = mxGetPr(plhs[4]);
    
    if(odds == NULL)
        mexErrMsgTxt("Not enough memory to create odds vector!");

    odds[0] = expphotonratediff*pow(photonodds,photonbins[0]);
    for(i = 1; i < numbins; i++)
    {
        odds[i] = expphotonratediff*pow(photonodds,photonbins[i])*(odds[i-1]*(1-k01dt)+k10dt)/(odds[i-1]*k01dt+1-k10dt);
    }
    
    backwards = expphotonratediff*pow(photonodds,photonbins[numbins-1]);
    for(i = numbins-2; i >= 0; i--)
    {
        currtranscounts[0] = (1-k01dt)*odds[i]/(1+odds[i])*backwards/(1+backwards);
        currtranscounts[1] = k01dt*odds[i]/(1+odds[i])/(1+backwards);
        currtranscounts[2] = k10dt/(1+odds[i])*backwards/(1+backwards);
        currtranscounts[3] = (1-k10dt)/(1+odds[i])/(1+backwards);
        sumcurrtranscounts = currtranscounts[0]+currtranscounts[1]+currtranscounts[2]+currtranscounts[3];
        transcounts[0] += currtranscounts[0]/sumcurrtranscounts;
        transcounts[1] += currtranscounts[1]/sumcurrtranscounts;
        transcounts[2] += currtranscounts[2]/sumcurrtranscounts;
        transcounts[3] += currtranscounts[3]/sumcurrtranscounts;
        odds[i] = odds[i]*(backwards*(1-k01dt)+k01dt)/(backwards*k10dt+1-k10dt);
        photoncounts[0] += photonbins[i]*odds[i]/(1+odds[i]);
        photoncounts[1] += photonbins[i]*1/(1+odds[i]);
        backwards = expphotonratediff*pow(photonodds,photonbins[i])*(backwards*(1-k01dt)+k01dt)/(backwards*k10dt+1-k10dt);
    }
    plhs[0] = mxCreateDoubleScalar(photoncounts[0]/(transcounts[0]+transcounts[1]));
    plhs[1] = mxCreateDoubleScalar(photoncounts[1]/(transcounts[2]+transcounts[3]));
    plhs[2] = mxCreateDoubleScalar(transcounts[1]/(transcounts[0]+transcounts[1]));
    plhs[3] = mxCreateDoubleScalar(transcounts[2]/(transcounts[2]+transcounts[3]));
}
