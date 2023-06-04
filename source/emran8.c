/*  S-Function for Extended Minimum Allocation Resource Network   *****************************/
/*  Giampiero Campa, June 2007   **************************************************************/

#define S_FUNCTION_NAME emran8
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"
#include <math.h>

/* mdlCheckParameters, check parameters, this routine is called later from mdlInitializeSizes */
#define MDL_CHECK_PARAMETERS
static void mdlCheckParameters(SimStruct *S) {
	
	/* declaration of the variables used in mdlCheckParameters                                */
	real_T   *pr, ovl, radius, limW, limS, limC, T; 
	real_T E1thr, E2min, E2max, E2gam, E3thr, E3gam, etaW, etaS, etaC, gamW, gamS, gamC;
	int_T Ni, No, Nmax, i, el, nEls, Nstates, prune;
	
	/* Basic check : All parameters must be real finite vectors, this is a standard check     */
	for (i = 0; i < ssGetNumSFcnParams(S); i++) {
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
	
	
	/* Check parameter # 1 : Number of Inputs *************************************************/
	if ( mxGetNumberOfElements(ssGetSFcnParam(S,0)) != 2 )
	{ ssSetErrorStatus(S,"The Dimension parameter must be a 2 elements vector"); return; }
	
	pr=mxGetPr(ssGetSFcnParam(S,0));      Ni=(int_T)pr[0];No=(int_T)pr[1];
	
	if ((Ni < 1))
	{ ssSetErrorStatus(S,"The Number of inputs must be greater than 1"); return; }
	
	if ((No < 1))
	{ ssSetErrorStatus(S,"The Number of inputs must be greater than 1"); return; }
	
	
	/* Check parameter # 2: [Nmax ovl radius prune] *******************************************/
	if ( mxGetNumberOfElements(ssGetSFcnParam(S,1)) != 4 )
	{ ssSetErrorStatus(S,"The overlapping parameter must be scalar"); return; }
	
	pr=mxGetPr(ssGetSFcnParam(S,1));
	Nmax=(int_T)pr[0];ovl=pr[1];radius=pr[2];prune=(int_T)pr[3];
	
	if ((ovl <= 0))
	{ ssSetErrorStatus(S,"The overlapping parameter must be > 0"); return; }
	
	if ((Nmax < 2))
	{ ssSetErrorStatus(S,"Nmax must be greater than 1"); return; }
	
	if ((radius < 0))
	{ ssSetErrorStatus(S,"Radius must not be less than 0"); return; }
	
	if ((prune != 0) && (prune != 1))
	{ ssSetErrorStatus(S,"Prune must be either 0 or 1"); return; }
	
	
	
	/* Check parameter # 3: etaW etaS etaC ****************************************************/
	if ( mxGetNumberOfElements(ssGetSFcnParam(S,2)) != 3 )
	{ ssSetErrorStatus(S," eta parameter must be a vector"); return; }
	
	pr=mxGetPr(ssGetSFcnParam(S,2));      etaW=pr[0];  etaS=pr[1]; etaC=pr[2];
	
	if ((etaW < 0) | (etaS < 0)  | (etaC < 0) )
	{ ssSetErrorStatus(S,"eta must not be less than 0"); return; }
	
	
	
	/* Check parameter # 4: gamW gamS gamC ****************************************************/
	if ( mxGetNumberOfElements(ssGetSFcnParam(S,3)) != 3 )
	{ ssSetErrorStatus(S," gam parameter must be a vector"); return; }

	pr=mxGetPr(ssGetSFcnParam(S,3));      gamW=pr[0];  gamS=pr[1]; gamC=pr[2];

	if ((gamW < 0) | (gamS < 0)  | (gamS < 0) )
	{ ssSetErrorStatus(S,"eta must not be less than 0"); return; }
	

	/* Check parameter # 5: limW limS limC ****************************************************/
	if ( mxGetNumberOfElements(ssGetSFcnParam(S,4)) != 3 )
	{ ssSetErrorStatus(S," lim parameter must be a vector"); return; }

	pr=mxGetPr(ssGetSFcnParam(S,4));      limW=pr[0];  limS=pr[1]; limC=pr[2];

	if ((limW < 0) | (limS < 0)  | (limS < 0) )
	{ ssSetErrorStatus(S,"lim must not be less than 0"); return; }
	
	
	/* Check parameter # 6: E1 threshold ******************************************************/
	if ( mxGetNumberOfElements(ssGetSFcnParam(S,5)) != 1 )
	{ ssSetErrorStatus(S,"The E1thr parameter must be a scalar"); return; }
	
	pr=mxGetPr(ssGetSFcnParam(S,5));      E1thr=pr[0];
	
	if ((E1thr < 0))
	{ ssSetErrorStatus(S,"E1thr must not be less than 0"); return; }
	
	
		
	/* Check parameter # 7: E2 max (initial) threshold, min (final) threshold, decay rate *****/
	if ( mxGetNumberOfElements(ssGetSFcnParam(S,6)) != 3 )
	{ ssSetErrorStatus(S,"The E2 parameter must be a vector [E2max E2min E2gam]"); return; }
	
	pr=mxGetPr(ssGetSFcnParam(S,6));      E2max=pr[0],E2min=pr[1],E2gam=pr[2];
	
	if ((E2max < 0) | (E2min < 0)  | (E2gam < 0) )
	{ ssSetErrorStatus(S,"The elements of E2 must not be less than 0"); return; }
	
	
	
	/* Check parameter # 8: E3 threshold, E3 filter decay rate ********************************/
	if ( mxGetNumberOfElements(ssGetSFcnParam(S,7)) != 2 )
	{ ssSetErrorStatus(S,"The E3 parameter must be a vector [E3Thr E3Pole]"); return; }
	
	pr=mxGetPr(ssGetSFcnParam(S,7));      E3thr=pr[0];E3gam=pr[1];
	
	if ((E3thr < 0) | (E3gam < 0))
	{ ssSetErrorStatus(S,"The elements of E3 must not be less than 0"); return; }		
		
	
	
	/* Check number of elements in parameter x0                                               */
	Nstates=(Ni*Nmax+2*Nmax+2)*No;
	if (mxGetNumberOfElements(ssGetSFcnParam(S,8)) != Nstates )
	{ ssSetErrorStatus(S,"x0 must be a [(Nx+2)*Nmax+2]*No size vector"); return; }
	
	
	
	/* Check parameter # 9: sampling time *****************************************************/
	if ( mxGetNumberOfElements(ssGetSFcnParam(S,9)) != 1 )
	{ ssSetErrorStatus(S,"The T parameter must be a scalar"); return; }
	
	pr=mxGetPr(ssGetSFcnParam(S,9));      T=pr[0];
	
	if ((T <= 0))
	{ ssSetErrorStatus(S,"T must be greater than 0"); return; }
	
}



/* mdlInitializeSizes - initialize the sizes array ********************************************/
static void mdlInitializeSizes(SimStruct *S) {
	real_T *pr ;
	int_T  Nstates, Ni, No, Nmax ;
	
    ssSetNumSFcnParams(S,10);                           /* number of expected parameters       */
	
    /* Check the number of parameters and then calls mdlCheckParameters to see if they are ok */
    if (ssGetNumSFcnParams(S) == ssGetSFcnParamsCount(S))
	{ mdlCheckParameters(S); if (ssGetErrorStatus(S) != NULL) return; } else return;
	
	pr=mxGetPr(ssGetSFcnParam(S,0));     Ni=(int_T)pr[0];No=(int_T)pr[1];
	pr=mxGetPr(ssGetSFcnParam(S,1));   Nmax=(int_T)pr[0];
	
	Nstates=(Ni*Nmax+2*Nmax+1+1)*No;
	
/*  For EACH OUTPUT there is a vector of Nmax*(Ni+2)+2 states organized like this: 
/*  X=[1.............................Nmax*(Ni),           % neurons coordinates
	   Nmax*(Ni)+1...................Nmax*(Ni+1),         % neurons sigma
	   Nmax*(Ni+1)+1.................Nmax*(Ni+2),         % neurons weigths
	   Nmax*(Ni+2)+1                                      % number of current neurons
	   Nmax*(Ni+2)+2                             ]        % state of error E3 filter          */
	
	ssSetNumContStates(S,0);                          /* number of continuous states          */
	
    /* number of discrete states                                                              */
    ssSetNumDiscStates(S,Nstates);
	
    if (!ssSetNumInputPorts(S,3)) return;             /* number of input ports                */
	
    ssSetInputPortWidth(S,0,Ni);                      /* first input port width (# NN inputs) */
	ssSetInputPortDirectFeedThrough(S,0,1);           /* first port direct feedthrough flag   */
	
    ssSetInputPortWidth(S,1,No);                      /* second input port width (e)          */
	ssSetInputPortDirectFeedThrough(S,1,0);           /* second port direct feedthrough flag  */
	
    ssSetInputPortWidth(S,2,1);                       /* third input port width  (LE)         */
    ssSetInputPortDirectFeedThrough(S,2,0);           /* third port direct feedthrough flag   */
	
	if (!ssSetNumOutputPorts(S,2)) return;	          /* number of output ports               */
	
	/* first and second output ports width                                                    */
	ssSetOutputPortWidth(S,0,No);
	ssSetOutputPortWidth(S,1,Nstates);
	
	ssSetNumSampleTimes(S,1);                          /* number of sample times              */
	ssSetNumRWork(S,2*Nmax);                           /* number real_T work vector elements  */
	ssSetNumIWork(S,0);                                /* number int_T work vector elements   */
	ssSetNumPWork(S,0);                                /* number ptr work vector elements     */
	ssSetNumModes(S,0);                                /* number mode work vector elements    */
	ssSetNumNonsampledZCs(S,0);                        /* number of nonsampled zero crossing  */
}

/* mdlInitializeSampleTimes - initialize the sample times array *******************************/
static void mdlInitializeSampleTimes(SimStruct *S) {
	real_T  *T=mxGetPr(ssGetSFcnParam(S,9));		              /* pointer to T             */
	
	/* Set things up to run with inherited sample time                                        */
	ssSetSampleTime(S, 0, T[0]);
	ssSetOffsetTime(S, 0, 0);
}

/* mdlStart - initialize work vectors *********************************************************/
#define MDL_START
static void mdlStart(SimStruct *S) {
	real_T *pr, *WORK;
	int_T     i,Nmax;
	
	pr=mxGetPr(ssGetSFcnParam(S,1));      Nmax=(int_T)pr[0];
	WORK = ssGetRWork(S);
	
	for (i=0;i<2*Nmax;i++)             /* [d2(0)....d2(Nmax-1),activ(0).....activ(Nmax-1)]    */
		WORK[i]=0;
}

/* mdlInitializeConditions - initialize the states ********************************************/
#define MDL_INITIALIZE_CONDITIONS
static void mdlInitializeConditions(SimStruct *S) {
	
	real_T  *x0  = ssGetRealDiscStates(S);
	real_T  *pr, *px0, tot;
	int_T Ni, No, Nmax, Nstates, i, k, ofY, ofC , ofS, ofW, ofN;
	
	pr=mxGetPr(ssGetSFcnParam(S,0));      Ni=(int_T)pr[0];No=(int_T)pr[1];
	pr=mxGetPr(ssGetSFcnParam(S,1));      Nmax=(int_T)pr[0];
	px0=mxGetPr(ssGetSFcnParam(S,8));
	
	Nstates=((Ni+2)*Nmax+2)*No;
	
	tot=0;for(i=0;i<Nstates;i++) tot+=fabs(px0[i]);
	
	if(tot<1e-10) {
		ofC=0;
		ofS=Nmax*Ni;
		ofW=Nmax*(Ni+1);
		ofN=Nmax*(Ni+2);

		for (k=0;k<No;k++) {
			ofY=k*((Ni+2)*Nmax+2);
			for (i=0;i<(Ni+2)*Nmax+2;i++) x0[ofY+i]=0;

		    /* initialization of the first neuron                                             */
			x0[ofY+ofS]=0.0001;
			x0[ofY+ofN]=1;
		}
	}
	else {
		for(i=0;i<Nstates;i++)  x0[i]=px0[i];
	}
	
}

/* mdlOutputs - compute the outputs ***********************************************************/
static void mdlOutputs(SimStruct *S, int_T tid) {
	
	real_T   *x   = ssGetRealDiscStates(S);
	real_T   *y   = ssGetOutputPortRealSignal(S,0);
	real_T   *Y   = ssGetOutputPortRealSignal(S,1);
	
	InputRealPtrsType u = ssGetInputPortRealSignalPtrs(S,0);
	
	real_T  *pr, ys, d2, *WORK;
	int_T Ni, No, Nmax, Nstates, i, j, k, ofY, ofC, ofS, ofW, ofN;
	
	pr=mxGetPr(ssGetSFcnParam(S,0));      Ni=(int_T)pr[0];No=(int_T)pr[1];
	pr=mxGetPr(ssGetSFcnParam(S,1));      Nmax=(int_T)pr[0];
	WORK=ssGetRWork(S);
	
    Nstates=((Ni+2)*Nmax+2)*No;
	
	ofC=0;
	ofS=Nmax*Ni;
	ofW=Nmax*(Ni+1);
	ofN=Nmax*(Ni+2);
	
	for (k=0;k<No;k++) {
        
		ys=0;
		ofY=k*((Ni+2)*Nmax+2);

		for (i=0;i<(int_T)x[ofY+ofN];i++) {
			d2=0;
			for (j=0;j<Ni;j++)
				d2=d2+(*u[j]-x[ofY+j*Nmax+i])*(*u[j]-x[ofY+j*Nmax+i]);
			
			ys=ys+x[ofY+ofW+i]*exp(-d2/(2*x[ofY+ofS+i]*x[ofY+ofS+i]));
		}
		y[k]=ys;
	}
	
	for(i=0;i<Nstates;i++) Y[i]=x[i];
}

/* mdlUpdate - perform action at major integration time step **********************************/
#define MDL_UPDATE
static void mdlUpdate(SimStruct *S, int_T tid)
{
	real_T  *pr, *WORK, E1thr, E2thr, E2min, E2max, E2gam, E3thr, E3gam, e1, e2, e3; 
	real_T  etaW, etaS, etaC, gamW, gamS, gamC, limW, limS, limC;
	real_T  ovl, radius, d2, d_min, xf, act_min;
	int_T   ofY, ofC, ofS, ofW, ofN, Ni, No, Nmax, i, j, k, option, N, I_min, prune;
	
	real_T  *x= ssGetRealDiscStates(S), t=ssGetT(S);
	
	InputRealPtrsType u  = ssGetInputPortRealSignalPtrs(S,0);
	InputRealPtrsType e  = ssGetInputPortRealSignalPtrs(S,1);
	InputRealPtrsType LE = ssGetInputPortRealSignalPtrs(S,2);
	
	if (*LE[0]!=0) {

		WORK=ssGetRWork(S);
		
		pr=mxGetPr(ssGetSFcnParam(S,0));      Ni=(int_T)pr[0];No=(int_T)pr[1];
		pr=mxGetPr(ssGetSFcnParam(S,1));      Nmax=(int_T)pr[0];ovl=pr[1];radius=pr[2];prune=(int_T)pr[3];
		
		pr=mxGetPr(ssGetSFcnParam(S,2));      etaW=pr[0];  etaS=pr[1]; etaC=pr[2];
		pr=mxGetPr(ssGetSFcnParam(S,3));      gamW=pr[0];  gamS=pr[1]; gamC=pr[2];
		pr=mxGetPr(ssGetSFcnParam(S,4));      limW=pr[0];  limS=pr[1]; limC=pr[2];
		pr=mxGetPr(ssGetSFcnParam(S,5));      E1thr=pr[0];
		pr=mxGetPr(ssGetSFcnParam(S,6));      E2max=pr[0],E2min=pr[1],E2gam=pr[2];
		pr=mxGetPr(ssGetSFcnParam(S,7));      E3thr=pr[0];E3gam=pr[1];
		
		ofC=0;
		ofS=Nmax*Ni;
		ofW=Nmax*(Ni+1);
		ofN=Nmax*(Ni+2);
		
		for (k=0;k<No;k++) {
			ofY=k*((Ni+2)*Nmax+2);
			N=(int_T)x[ofY+ofN];
		
			d_min=1e100;
			act_min=1e100;
			I_min=0;
			
			/* filter for e  (third criteria E3) */
			x[ofY+ofN+1]=E3gam*x[ofY+ofN+1]+(1-E3gam)*(*e[k]);  

			/* calculate activation function */
			
			for (i=0;i<N;i++) {
				d2=0;
				for (j=0;j<Ni;j++)
					d2+=(*u[j]-x[ofY+j*Nmax+i])*(*u[j]-x[ofY+j*Nmax+i]);
				
				WORK[i]=sqrt(d2);
				if(d_min>WORK[i]) d_min=WORK[i];
				
				WORK[Nmax+i]=exp(-d2/(2*x[ofY+ofS+i]*x[ofY+ofS+i]));    
				
				/* pruning index */
				if( (WORK[Nmax+i]<act_min) && (N==Nmax) && (prune==1) ) {I_min=i; act_min=WORK[Nmax+i];}
			}
			

			/* parameters for the error criteria */
			
			e1=fabs(*e[k]);
			e2=d_min;
			e3=fabs(x[ofY+ofN+1]);
			
			E2thr=pow(E2gam,t)*E2max;
			if(E2thr<E2min)  E2thr=E2min;


			/* test for the addition of a new neuron */
			
			if( (e1>E1thr) & (e2>E2thr) &  (e3>E3thr) & (N<Nmax)   )  
				option=1;  /* must add one  */
			else 
				option=0;  /* update weights only  */
			
			if( (e1>E1thr) & (e2>E2thr) &  (e3>E3thr) & (N==Nmax) & (prune==1) ) 
				option=2;  /* must prune before adding */
						   
			
			/* actually update the state */

			switch (option) {
			case 0:  /* only update the neurons weights (INSIDE the radius) */
				
				if (radius<d_min) radius=d_min; /* if radius=0 == EMRAN, if radius -> inf == MRAN */
				
				for (i=0;i<N;i++) {
					if(WORK[i]<=radius) {

						/* sigmas */
						x[ofY+ofS+i] += etaS*x[ofY+ofW+i]*WORK[i]*WORK[i]*WORK[Nmax+i]/(x[ofY+ofS+i]*x[ofY+ofS+i]*x[ofY+ofS+i])*(*e[k])-etaS*gamS*x[ofY+ofS+i];
						if(x[ofY+ofS+i] < 1e-4) x[ofY+ofS+i]=1e-4;
						if(x[ofY+ofS+i] > limS) x[ofY+ofS+i]=limS;

						/* weights */
						x[ofY+ofW+i] += etaW*WORK[Nmax+i]*(*e[k])-etaW*gamW*x[ofY+ofW+i];
						if(x[ofY+ofW+i] < -limW) x[ofY+ofW+i]=-limW;
						if(x[ofY+ofW+i] >  limW) x[ofY+ofW+i]= limW;

						/* centers */
							if(etaC>0)
								for (j=0;j<Ni;j++) {
									x[ofY+i+j*Nmax] += etaC*x[ofY+ofW+i]*(*u[j]-x[ofY+i+j*Nmax])*WORK[Nmax+i]/(x[ofY+ofS+i]*x[ofY+ofS+i])*(*e[k])-etaC*gamC*x[ofY+i+j*Nmax];
									if(x[ofY+i+j*Nmax]< -limC) x[ofY+i+j*Nmax]= -limC;
									if(x[ofY+i+j*Nmax]>  limC) x[ofY+i+j*Nmax]=  limC;
								}
					}
				}
				break;
				
			case 1:  /* just add a new neuron */

				/* center */
				for (j=0;j<Ni;j++) {
					x[ofY+ofC+j*Nmax+N]=(*u[j]);
					if(x[ofY+ofC+j*Nmax+N]< -limC) x[ofY+ofC+j*Nmax+N]= -limC;
					if(x[ofY+ofC+j*Nmax+N]>  limC) x[ofY+ofC+j*Nmax+N]=  limC;
				}
				
				/* sigma */
				x[ofY+ofS+N]=ovl*d_min;
				if(x[ofY+ofS+N] < 1e-4) x[ofY+ofS+N]=1e-4;
				if(x[ofY+ofS+N] > limS) x[ofY+ofS+N]=limS;

				/* weight */
				x[ofY+ofW+N]=(*e[k]);
				if(x[ofY+ofW+N] < -limW) x[ofY+ofW+N]=-limW;
				if(x[ofY+ofW+N] >  limW) x[ofY+ofW+N]= limW;

				/* number of active neurons */
				x[ofY+ofN]=N+1;

				break;
				
			case 2:  /* prune the worst neuron and add the new neuron  */
				
				/* center */
				for (j=0;j<Ni;j++) {
					x[ofY+ofC+j*Nmax+I_min]=(*u[j]);
					if(x[ofY+ofC+j*Nmax+I_min]< -limC) x[ofY+ofC+j*Nmax+I_min]= -limC;
					if(x[ofY+ofC+j*Nmax+I_min]>  limC) x[ofY+ofC+j*Nmax+I_min]=  limC;
				}

				/* sigma */
				x[ofY+ofS+I_min]=ovl*d_min;
				if(x[ofY+ofS+I_min] < 1e-4) x[ofY+ofS+I_min]=1e-4;
				if(x[ofY+ofS+I_min] > limS) x[ofY+ofS+I_min]=limS;

				/* weight */
				x[ofY+ofW+I_min]=(*e[k]);
				if(x[ofY+ofW+I_min] < -limW) x[ofY+ofW+I_min]=-limW;
				if(x[ofY+ofW+I_min] >  limW) x[ofY+ofW+I_min]= limW;

				break;
				
			} /* close switch ( state updated at this point) */


		} /* close cycle trough different outputs */
	} /* close if LE!=0 */
} /* close mdlupdate */


/* mdlTerminate - called when the simulation is terminated ************************************/
static void mdlTerminate(SimStruct *S) {}


/* Trailer information to set everything up for simulink usage ********************************/
#ifdef  MATLAB_MEX_FILE                       /* Is this file being compiled as a MEX-file?   */
#include "simulink.c"                         /* MEX-file interface mechanism                 */
#else
#include "cg_sfun.h"                          /* Code generation registration function        */
#endif

