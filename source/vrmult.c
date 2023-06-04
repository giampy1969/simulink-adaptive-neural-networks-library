/*  S-Function for Simulink Real Matrices Multiplication **************************************/
/*  Giampiero Campa 27-August-00 **************************************************************/

#define S_FUNCTION_NAME vrmult
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"

/* mdlCheckParameters, check parameters, this routine is called later from mdlInitializeSizes */
#define MDL_CHECK_PARAMETERS
static void mdlCheckParameters(SimStruct *S)
{
    /* Basic check : All parameters must be real positive vectors                             */
    real_T *pr;                            

    int_T  i, el, nEls;
    for (i = 0; i < 1; i++) {
        if (mxIsEmpty(    ssGetSFcnParam(S,i)) || mxIsSparse(   ssGetSFcnParam(S,i)) ||
            mxIsComplex(  ssGetSFcnParam(S,i)) || !mxIsNumeric( ssGetSFcnParam(S,i))  )
                  { ssSetErrorStatus(S,"Parameters must be real finite vectors"); return; } 
        pr   = mxGetPr(ssGetSFcnParam(S,i));
        nEls = mxGetNumberOfElements(ssGetSFcnParam(S,i));
        for (el = 0; el < nEls; el++) {
            if (!mxIsFinite(pr[el])) 
                  { ssSetErrorStatus(S,"Parameters must be real finite vectors"); return; }
        }
    }

    /* Check number of elements in parameter: [no nm ni]                                      */
    if ( mxGetNumberOfElements(ssGetSFcnParam(S,0)) != 3 )
    { ssSetErrorStatus(S,"The parameter must be a 3 elements vector"); return; }

    /* get the basic parameters and check them                                                */
    pr=mxGetPr(ssGetSFcnParam(S,0));
    if ( (pr[0] < 1) | (pr[1] < 1) | (pr[2] < 1) )
    { ssSetErrorStatus(S,"Dimensions must be greater than zero"); return; }

}

/* mdlInitializeSizes - initialize the sizes array ********************************************/
static void mdlInitializeSizes(SimStruct *S)
{
    real_T *n;                            

    ssSetNumSFcnParams(S,1);                          /* number of expected parameters        */

    /* Check the number of parameters and then calls mdlCheckParameters to see if they are ok */
    if (ssGetNumSFcnParams(S) == ssGetSFcnParamsCount(S))
    { mdlCheckParameters(S); if (ssGetErrorStatus(S) != NULL) return; } else return;
    n=mxGetPr(ssGetSFcnParam(S,0));

    ssSetNumContStates(S,0);                          /* number of continuous states          */
    ssSetNumDiscStates(S,0);                          /* number of discrete states            */

    if (!ssSetNumInputPorts(S,2)) return;             /* number of input ports                */
    ssSetInputPortWidth(S,0,(int_T)(n[0]*n[1]));      /* first input port width               */
    ssSetInputPortWidth(S,1,(int_T)(n[1]*n[2]));      /* second input port width              */
    ssSetInputPortDirectFeedThrough(S,0,1);           /* first port direct feedthrough flag   */
    ssSetInputPortDirectFeedThrough(S,1,1);           /* second port direct feedthrough flag  */

    if (!ssSetNumOutputPorts(S,1)) return;            /* number of output ports               */
    ssSetOutputPortWidth(S,0,(int_T)(n[0]*n[2]));     /* first output port width              */
   
    ssSetNumSampleTimes(S,0);                         /* number of sample times               */

    ssSetNumRWork(S,0);                               /* number real work vector elements     */
    ssSetNumIWork(S,3);                               /* number int_T work vector elements    */
    ssSetNumPWork(S,0);                               /* number ptr work vector elements      */
    ssSetNumModes(S,0);                               /* number mode work vector elements     */
    ssSetNumNonsampledZCs(S,0);                       /* number of nonsampled zero crossing   */
}

/* mdlInitializeSampleTimes - initialize the sample times array *******************************/
static void mdlInitializeSampleTimes(SimStruct *S)
{
    /* Set things up to run with inherited sample time                                        */
    ssSetSampleTime(S, 0, INHERITED_SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0);
}

/* mdlStart - initialize work vectors *********************************************************/
#define MDL_START
static void mdlStart(SimStruct *S)
{
real_T     *pr = mxGetPr(ssGetSFcnParam(S,0));
int_T   i, *iv = ssGetIWork(S);

for (i=0;i<3;i++) iv[i]=(int_T)(pr[i]);
}

/* mdlOutputs - compute the outputs ***********************************************************/
static void mdlOutputs(SimStruct *S, int_T tid)
{
int_T     i, j, k, *n  = ssGetIWork(S);
real_T             *y  = ssGetOutputPortRealSignal(S,0);
InputRealPtrsType up1  = ssGetInputPortRealSignalPtrs(S,0);
InputRealPtrsType up2  = ssGetInputPortRealSignalPtrs(S,1);

for(i = 0; i < n[0]; i++)
   for(j = 0; j < n[2]; j++)
	   for( y[i+j*n[0]] = 0, k = 0; k < n[1]; k++)
            y[i+j*n[0]] += (*up1[i+k*n[0]])*(*up2[k+j*n[1]]);
}

/* mdlTerminate - called when the simulation is terminated ***********************************/
static void mdlTerminate(SimStruct *S) {}

/* Trailer information to set everything up for simulink usage *******************************/
#ifdef  MATLAB_MEX_FILE                      /* Is this file being compiled as a MEX-file?   */
#include "simulink.c"                        /* MEX-file interface mechanism                 */
#else
#include "cg_sfun.h"                         /* Code generation registration function        */
#endif
