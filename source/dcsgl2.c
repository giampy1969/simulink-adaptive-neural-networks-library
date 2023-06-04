/*  S-Function for Dynamic cell structure NN    ***********************************************/
/*  Giampiero Campa, June 2007 ****************************************************************/

#define S_FUNCTION_NAME dcsgl2
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"
#include <math.h>

/* mdlCheckParameters, check parameters, this routine is called later from mdlInitializeSizes */
#define MDL_CHECK_PARAMETERS
static void mdlCheckParameters(SimStruct *S) {
    /* Basic check : All parameters must be integer positive vectors                          */
	real_T  *pr, T;
	int_T   i, el, nEls, Ni, No, Nmax;

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


	/* Check number parameter # 1 : [Ni No] = [Nx Ny]                                        */
	if ( mxGetNumberOfElements(ssGetSFcnParam(S,0)) != 2 )
	{ ssSetErrorStatus(S,"[Nx Ny] must be a 2 dimensional vector"); return; }
	
	pr=mxGetPr(ssGetSFcnParam(S,0));     
	
	for (i=0;i<2;i++)
	{ if ((int_T)pr[i] < 1) { ssSetErrorStatus(S,"Nx and Ny must be greater than 0"); return; }}
	Ni= (int_T) pr[0];No= (int_T) pr[1];


	/* Check number parameter # 2 : norlat = [Nmax Overlap RsThr Lambda Alpha Theta]         */
	if ( mxGetNumberOfElements(ssGetSFcnParam(S,1)) != 6 )
	{ ssSetErrorStatus(S,"The second parameter must be a 6 dimensional vector"); return; }

	pr=mxGetPr(ssGetSFcnParam(S,1));     
	for (i=0;i<6;i++)
	{ if (pr[i] < 0) { ssSetErrorStatus(S,"The second parameter must not be less than 0"); return; }}
	Nmax= (int_T) pr[0];


	/* Check number parameter # 3 : eps = [Epsb Epsn etaW etaS etaC]                         */
	if ( mxGetNumberOfElements(ssGetSFcnParam(S,2)) != 5 )
	{ ssSetErrorStatus(S,"The third parameter must be a 5 dimensional vector"); return; }

	pr=mxGetPr(ssGetSFcnParam(S,2));     
	for (i=0;i<5;i++)
	{ if (pr[i] < 0) { ssSetErrorStatus(S,"The third parameter must not be less than 0"); return; }}


	/* Check number parameter # 4 : lim = [limW limS limC]                         */
	if ( mxGetNumberOfElements(ssGetSFcnParam(S,3)) != 3 )
	{ ssSetErrorStatus(S,"The fourth parameter must be a 3 dimensional vector"); return; }

	pr=mxGetPr(ssGetSFcnParam(S,3));     
	for (i=0;i<3;i++)
	{ if (pr[i] < 0) { ssSetErrorStatus(S,"The fourth parameter must not be less than 0"); return; }}

	/* Check number parameter # 5 : act                         */
	if ( mxGetNumberOfElements(ssGetSFcnParam(S,4)) != 1 )
	{ ssSetErrorStatus(S,"The fifth parameter must be a scalar"); return; }

	pr=mxGetPr(ssGetSFcnParam(S,4));     
	if (pr[0] < 0 || pr[0]>1) { ssSetErrorStatus(S,"The fifth parameter must be either 1 or 0"); return; }


	/* Check number parameter # 6 : X0							      						 */
	if ( mxGetNumberOfElements(ssGetSFcnParam(S,5)) !=  (Nmax*(Nmax+Ni+5)+2)*No )
	{ ssSetErrorStatus(S,"The initial state must have (Nmax*(Nmax+Nx+5)+2)*Ny elements"); return; }


	/* Check number parameter # 7 : T														 */
    pr=mxGetPr(ssGetSFcnParam(S,6));      T=pr[0];
    if (T < 0) { ssSetErrorStatus(S," Sample Time (T) must be greater than 0"); return; }

}


/* mdlInitializeSizes - initialize the sizes array ********************************************/
static void mdlInitializeSizes(SimStruct *S) {
	 real_T *pr ;
	 int_T  Ni, No, Nmax ;

     ssSetNumSFcnParams(S,7);                          /* number of expected parameters       */
  
    /* Check the number of parameters and then calls mdlCheckParameters to see if they are ok */
    if (ssGetNumSFcnParams(S) == ssGetSFcnParamsCount(S))
	{ mdlCheckParameters(S); if (ssGetErrorStatus(S) != NULL) return; } else return;

	/* Get Ni No Nmax   																	  */
	pr=mxGetPr(ssGetSFcnParam(S,0));     Ni=(int_T)pr[0];No=(int_T)pr[1];
	pr=mxGetPr(ssGetSFcnParam(S,1));     Nmax=(int_T)pr[0];

	ssSetNumContStates(S,0);                          /* number of continuous states          */

    /* number of discrete states                                                              */
    ssSetNumDiscStates(S,No*(2+Nmax*(Nmax+Ni+5)));

    if (!ssSetNumInputPorts(S,3)) return;             /* number of input ports                */

    ssSetInputPortWidth(S,0,Ni);                      /* first input port width (# NN inputs) */
	ssSetInputPortDirectFeedThrough(S,0,1);           /* first port direct feedthrough flag   */

    ssSetInputPortWidth(S,1,No);                      /* second input port width (e)          */
	ssSetInputPortDirectFeedThrough(S,1,0);           /* second port direct feedthrough flag  */

    ssSetInputPortWidth(S,2,1);                       /* third input port width  (LE)         */
    ssSetInputPortDirectFeedThrough(S,2,0);           /* third port direct feedthrough flag   */


	if (!ssSetNumOutputPorts(S,2)) return;			  /* number of output ports               */

	/* first and second output ports width                                                    */
	ssSetOutputPortWidth(S,0,No);
	ssSetOutputPortWidth(S,1,No*(2+Nmax*(Nmax+Ni+5)));

	ssSetNumSampleTimes(S,1);                         /* number of sample times               */
	ssSetNumRWork(S,Nmax);                            /* number real_T work vector elements   */
	ssSetNumIWork(S,0);                               /* number int_T work vector elements    */
	ssSetNumPWork(S,0);                               /* number ptr work vector elements      */
	ssSetNumModes(S,0);                               /* number mode work vector elements     */
	ssSetNumNonsampledZCs(S,0);                       /* number of nonsampled zero crossing   */
}

/* mdlInitializeSampleTimes - initialize the sample times array *******************************/
static void mdlInitializeSampleTimes(SimStruct *S) {
	 real_T  *T=mxGetPr(ssGetSFcnParam(S,6));		   /* pointer to T */
	 /* Set things up to run with inherited sample time                                       */
	 ssSetSampleTime(S, 0, T[0]);
	 ssSetOffsetTime(S, 0, 0);
}


/* mdlStart - initialize work vectors *********************************************************/
#define MDL_START
static void mdlStart(SimStruct *S) {
	 real_T *pr, *D;
	 int_T     i,Nmax;

	 pr=mxGetPr(ssGetSFcnParam(S,1));      Nmax=(int_T)pr[0];
	 D = ssGetRWork(S);

	 for (i=0;i<Nmax;i++)   
		 D[i]=0;                 

}

/* mdlInitializeConditions - initialize the states ********************************************/
#define MDL_INITIALIZE_CONDITIONS
static void mdlInitializeConditions(SimStruct *S) {

	real_T  *x0  = ssGetRealDiscStates(S);
	real_T  *pr, *px0, tot;

	int_T   Ni, No, Nmax, n_states, i, j, k;
    int_T   ofY, ofct, ofCn, ofC, ofW, ofS, ofR, ofRm, ofT, ofN;
 
	/* Get Ni No Nmax   																	  */
	pr=mxGetPr(ssGetSFcnParam(S,0));     Ni=(int_T)pr[0];No=(int_T)pr[1];
	pr=mxGetPr(ssGetSFcnParam(S,1));     Nmax=(int_T)pr[0];

	n_states=(2+Nmax*(Nmax+Ni+5))*No;

 	px0=mxGetPr(ssGetSFcnParam(S,5));
	
	tot=0;for(i=0;i<n_states;i++) tot=tot+fabs(px0[i]);

	if(tot<0.0001)	{

		/* printf("standard initial conditions 2 nodes and %d states\n",n_states); */
	
		ofct=0;
		ofCn=ofct+1;
		ofC=ofCn+Nmax*Nmax;
 		ofW=ofC+Nmax*Ni;
		ofS=ofW+Nmax;
		ofR=ofS+Nmax;
		ofRm=ofR+Nmax;
		ofT=ofRm+Nmax;
		ofN=ofT+Nmax;

		for (k=0;k<No;k++) {

			ofY=k*(2+Nmax*(Nmax+Ni+5));

			for (j=1;j<=Ni;j++) {
				x0[ofY+ofC+(j-1)*Nmax]= 0.5;
				x0[ofY+ofC+(2-1)+(j-1)*Nmax]=-0.5;
			}      
    
			x0[ofY+ofW]=1;
			x0[ofY+ofW+1]=1;
			x0[ofY+ofS]=1;
			x0[ofY+ofS+1]=1;
			x0[ofY+ofCn+(1-1)+(2-1)*Nmax]=1;
			x0[ofY+ofCn+(2-1)]=1;
			x0[ofY+ofN]=2;
		}
     
    }
	else {
		for(i=0;i<n_states;i++)  
			x0[i]=px0[i];
	}
}


/* mdlOutputs - compute the outputs ***********************************************************/
static void mdlOutputs(SimStruct *S, int_T tid) {

	real_T   *x   = ssGetRealDiscStates(S);
	real_T   *y   = ssGetOutputPortRealSignal(S,0);
	real_T   *Y   = ssGetOutputPortRealSignal(S,1);

	InputRealPtrsType u = ssGetInputPortRealSignalPtrs(S,0);

    int_T   ofY, ofct, ofCn, ofC, ofW, ofS, ofR, ofRm, ofT, ofN;
	real_T  *pr, *D, ys, sum2, activ, d, dmn;
	int_T   Ni, No, Nmax, n_states, i, j, k, N; 
	int_T   bmu,sec,lga,r2;
	
	D=ssGetRWork(S);
	
	/* Get Ni No Nmax   																	  */
	pr=mxGetPr(ssGetSFcnParam(S,0));     Ni=(int_T)pr[0];No=(int_T)pr[1];
	pr=mxGetPr(ssGetSFcnParam(S,1));     Nmax=(int_T)pr[0];

	pr=mxGetPr(ssGetSFcnParam(S,4));      
	lga   = (int_T) pr[0];     // linear or gaussian activation 

	n_states=(2+Nmax*(Nmax+Ni+5))*No;

	ofct=0;
	ofCn=ofct+1;
	ofC=ofCn+Nmax*Nmax;
 	ofW=ofC+Nmax*Ni;
	ofS=ofW+Nmax;
	ofR=ofS+Nmax;
	ofRm=ofR+Nmax;
	ofT=ofRm+Nmax;
	ofN=ofT+Nmax;
	
	for (k=0;k<No;k++) {

		ofY=k*(2+Nmax*(Nmax+Ni+5));
  

		/* evaluates the bmu and the sbu to the current input */
		sum2=0;
		for(j=1;j<=Ni;j++)
			sum2+=(*u[j-1]-x[ofY+ofC+(1-1)+(j-1)*Nmax])*(*u[j-1]-x[ofY+ofC+(1-1)+(j-1)*Nmax]);
		D[0]=sqrt(sum2);

		sum2=0;
		for(j=1;j<=Ni;j++) 
			sum2+=(*u[j-1]-x[ofY+ofC+(2-1)+(j-1)*Nmax])*(*u[j-1]-x[ofY+ofC+(2-1)+(j-1)*Nmax]);
		D[1]=sqrt(sum2);


		if (D[1-1]<D[2-1]) { bmu=1;sec=2;} else { bmu=2;sec=1;}

		N=(int_T)x[ofY+ofN];

		for(i=3;i<=N;i++) {
		   sum2=0;
		   for (j=1;j<=Ni;j++)
				sum2+=(*u[j-1]-x[ofY+ofC+(i-1)+(j-1)*Nmax])*(*u[j-1]-x[ofY+ofC+(i-1)+(j-1)*Nmax]);
			D[i-1]=sqrt(sum2);  
			
			if (D[i-1]<D[bmu-1]) { sec=bmu;bmu=i;} else if (D[i-1]<D[sec-1]) sec=i;
		}
			

		/* calculate output */

		ys=0;
		if (lga!=0) {
			/* linear case */

			r2=1;
			dmn=1e37;
			for (i=1;i<=N;i++)
				if ( x[ofY+ofCn+(bmu-1)+(i-1)*Nmax]>0 )
					if ( D[i-1]<dmn ) {
						dmn=D[i-1];r2=i;
					}		

			d=0;
			for (j=1;j<=Ni;j++)
				d += (x[ofY+ofC+(j-1)*Nmax+(bmu-1)]-x[ofY+ofC+(j-1)*Nmax+(r2-1)]) * (x[ofY+ofC+(j-1)*Nmax+(bmu-1)]-x[ofY+ofC+(j-1)*Nmax+(r2-1)]);
			d=sqrt(d);

			if (dmn<d) {
				activ=D[bmu-1]/(D[bmu-1]+dmn);
				ys=activ*x[ofY+ofW+r2-1]+(1-activ)*x[ofY+ofW+bmu-1];
			} 
			else {
				ys=x[ofY+ofW+bmu-1];
			}

		} /* close linear case */
		else {
			/* gaussian case */

			for (j=1;j<=N;j++) {
				activ=exp( -(D[j-1]/(2*x[ofY+ofS+j-1]))*(D[j-1]/(2*x[ofY+ofS+j-1]))  );
				ys+=x[ofY+ofW+j-1]*activ;
			}

		} /* close gaussian case */
		y[k]=ys;

	} /* close cycle trough outputs */

	/* output current states as a second output*/
	for(i=0;i<n_states;i++) Y[i]=x[i];

} /* close mdloutput */


/* mdlUpdate - perform action at major integration time step *************************/
#define MDL_UPDATE
static void mdlUpdate(SimStruct *S, int_T tid) {

    int_T   ofY, ofct, ofCn, ofC, ofW, ofS, ofR, ofRm, ofT, ofN;
    real_T  alpha, gamma, theta, eps_b, eps_n, etaW, etaS, etaC;
	real_T  limW, limS, limC, trR, overl, sum2, lamda;
	int_T   Ni, No, Nmax, N, n_states, i, j, k, h, r1, r2, lga; 
	int_T   bmu, sec, conta;
	real_T  mx1, mx2, dmn, d, activ, rm1 ;

	real_T  *pr, w, *D, ys;
	real_T   *x= ssGetRealDiscStates(S), t=ssGetT(S);

	InputRealPtrsType u  = ssGetInputPortRealSignalPtrs(S,0);
	InputRealPtrsType e  = ssGetInputPortRealSignalPtrs(S,1);
	InputRealPtrsType LE = ssGetInputPortRealSignalPtrs(S,2);

	if (*LE[0]!=0) {

		D=ssGetRWork(S);

		pr=mxGetPr(ssGetSFcnParam(S,0));      
		Ni=(int_T)pr[0];		   // number of inputs
		No=(int_T)pr[1];		   // number of outputs

		pr=mxGetPr(ssGetSFcnParam(S,1));      
		Nmax  =(int_T)pr[0];       // maximum number of neurons 
		overl =pr[1];			   // overlapping factor for a new neuron  
		trR   =pr[2];              // threshlod for neuron insertion ( resource based)
		lamda =(int_T)pr[3];       // number of time between insertions  
		alpha =pr[4];              // threschold for connection deletion
		theta =pr[5];              // decay constant

		pr=mxGetPr(ssGetSFcnParam(S,2));      
		eps_b =pr[0];              // BMU  updata coefficient
		eps_n =pr[1];              // NEIG updata coefficient
		etaW  =pr[2];              // learning rate 
		etaS  =pr[3];              // learning rate
		etaC  =pr[4];              // learning rate

		pr=mxGetPr(ssGetSFcnParam(S,3));      
		limW  =pr[0];              // norm limiter 
		limS  =pr[1];              // norm limiter
		limC  =pr[2];              // norm limiter

		pr=mxGetPr(ssGetSFcnParam(S,4));      
		lga   = (int_T) pr[0];     // linear or gaussian activation 

		n_states=No*(2+Nmax*(Nmax+Ni+5));

		ofct=0;
		ofCn=ofct+1;
		ofC=ofCn+Nmax*Nmax;
 		ofW=ofC+Nmax*Ni;
		ofS=ofW+Nmax;
		ofR=ofS+Nmax;
		ofRm=ofR+Nmax;
		ofT=ofRm+Nmax;
		ofN=ofT+Nmax;
	
		for (k=0;k<No;k++) {

			ofY=(int_T) k*(2+Nmax*(Nmax+Ni+5));

			conta=(int_T)x[ofY+ofct];
			N=(int_T)x[ofY+ofN];

			/* evaluates the bmu and the sbu to the current input */
			sum2=0;
			for(j=1;j<=Ni;j++)  
				sum2+=(*u[j-1]-x[ofY+ofC+(1-1)+(j-1)*Nmax])*(*u[j-1]-x[ofY+ofC+(1-1)+(j-1)*Nmax]);
			D[0]=sqrt(sum2);

			sum2=0;
			for(j=1;j<=Ni;j++)   
				sum2+=(*u[j-1]-x[ofY+ofC+(2-1)+(j-1)*Nmax])*(*u[j-1]-x[ofY+ofC+(2-1)+(j-1)*Nmax]);
			D[1]=sqrt(sum2);

			if (D[1-1]<D[2-1]) { bmu=1;sec=2;} else { bmu=2;sec=1;}

			for(i=3;i<=N;i++) {
				sum2=0;
				for (j=1;j<=Ni;j++)
				sum2+=(*u[j-1]-x[ofY+ofC+(i-1)+(j-1)*Nmax])*(*u[j-1]-x[ofY+ofC+(i-1)+(j-1)*Nmax]);

				D[i-1]=sqrt(sum2);  

				if (D[i-1]<D[bmu-1]) {sec=bmu; bmu=i;} else if (D[i-1]<D[sec-1]) sec=i;
			}


			/* HEBBIAN LEARNING FOR THE CONNECTION MATRIX */
			j=bmu;   //updating only NEIG
			for(i=1;i<=N;i++) {
				if(x[ofY+ofCn+(j-1)+(i-1)*Nmax]<theta) {  
					x[ofY+ofCn+(j-1)+(i-1)*Nmax]=0;
					x[ofY+ofCn+(i-1)+(j-1)*Nmax]=0;
				}

			  else 
				if ( x[ofY+ofCn+(bmu-1)+(i-1)*Nmax]>0 ) {
					x[ofY+ofCn+(j-1)+(i-1)*Nmax]=alpha*x[ofY+ofCn+(j-1)+(i-1)*Nmax];
					x[ofY+ofCn+(i-1)+(j-1)*Nmax]=x[ofY+ofCn+(j-1)+(i-1)*Nmax];   
				}

			}

			x[ofY+ofCn+(bmu-1)+(sec-1)*Nmax]=1;
			x[ofY+ofCn+(sec-1)+(bmu-1)*Nmax]=1;


			/* KOHONEN LEARNING FOR CENTERS ADAPTATION */

			for(j=1;j<=Ni;j++) {
				x[ofY+ofC+(bmu-1)+(j-1)*Nmax] += eps_b*(*u[j-1]-x[ofY+ofC+(bmu-1)+(j-1)*Nmax]);
				if(x[ofY+ofC+(bmu-1)+(j-1)*Nmax]< -limC) x[ofY+ofC+(bmu-1)+(j-1)*Nmax]= -limC;
				if(x[ofY+ofC+(bmu-1)+(j-1)*Nmax]>  limC) x[ofY+ofC+(bmu-1)+(j-1)*Nmax]=  limC;
			}

			for(j=1;j<=N;j++)
			if( x[ofY+ofCn+(bmu-1)+(j-1)*Nmax]>0 )
				for (h=1;h<=Ni;h++) {
					x[ofY+ofC+(j-1)+(h-1)*Nmax] += eps_n*(*u[h-1]-x[ofY+ofC+(j-1)+(h-1)*Nmax]);
					if(x[ofY+ofC+(j-1)+(h-1)*Nmax]< -limC) x[ofY+ofC+(j-1)+(h-1)*Nmax]= -limC;
					if(x[ofY+ofC+(j-1)+(h-1)*Nmax]>  limC) x[ofY+ofC+(j-1)+(h-1)*Nmax]=  limC;
				}

   
			/* BMU RESORCE UPDATING */
			x[ofY+ofR+bmu-1]+=fabs(*e[k]);
			x[ofY+ofT+bmu-1]+=1;
			x[ofY+ofRm+bmu-1]=x[ofY+ofR+bmu-1]/x[ofY+ofT+bmu-1];
					

			sum2=0;
			for(i=1;i<=N;i++)
			  sum2+=x[ofY+ofRm+(i-1)];
				 
			rm1=sum2/N;
			conta++;
			if(conta>lamda) conta=0;


			/* TEST FOR A NEW NEURON */
			if( rm1>trR && N<Nmax && conta==lamda) {
				conta=0;
				N=N+1;   

				r1=1;   
				mx1=x[ofY+ofRm+(1-1)];

				for(i=2;i<=N;i++)
					if( x[ofY+ofRm+(i-1)]>mx1 ) {
						mx1=x[ofY+ofRm+(i-1)];
						r1=i;
				 }

				r2=1;
				mx2=0;
				for(j=1;j<=N;j++)
					if( x[ofY+ofCn+(r1-1)+(j-1)*Nmax]>0 )
						if ( x[ofY+ofRm+(j-1)]>mx2 ) {
							mx2=x[ofY+ofRm+(j-1)];
							r2=j;
						}

				gamma=x[ofY+ofRm+r1-1]/( x[ofY+ofRm+r1-1]+x[ofY+ofRm+r2-1] );

				/* center */
				for(j=1;j<=Ni;j++) {
				 x[ofY+ofC+(N-1)+(j-1)*Nmax]=x[ofY+ofC+(r1-1)+(j-1)*Nmax]+
					gamma*( x[ofY+ofC+(r2-1)+(j-1)*Nmax]-x[ofY+ofC+(r1-1)+(j-1)*Nmax]); 
					if(x[ofY+ofC+(N-1)+(j-1)*Nmax]< -limC) x[ofY+ofC+(N-1)+(j-1)*Nmax]= -limC;
					if(x[ofY+ofC+(N-1)+(j-1)*Nmax]>  limC) x[ofY+ofC+(N-1)+(j-1)*Nmax]=  limC;
				}

				/* weight */
				x[ofY+ofW+(N-1)]=x[ofY+ofW+(r1-1)]+gamma*(x[ofY+ofW+(r2-1)]-x[ofY+ofW+(r1-1)]);
				if(x[ofY+ofW+(N-1)] < -limW) x[ofY+ofW+(N-1)]=-limW;
				if(x[ofY+ofW+(N-1)] >  limW) x[ofY+ofW+(N-1)]= limW;

				/* distance for sigma */
				sum2=0;
				for(j=1;j<=Ni;j++)
				  sum2+=(x[ofY+ofC+(r1-1)+(j-1)*Nmax]-x[ofY+ofC+(r2-1)+(j-1)*Nmax])*
				  (x[ofY+ofC+(r1-1)+(j-1)*Nmax]-x[ofY+ofC+(r2-1)+(j-1)*Nmax]);
				d=sqrt(sum2);

				/* sigma */
				x[ofY+ofS+(N-1)]=overl*d;
				if(x[ofY+ofS+(N-1)] < 1e-3) x[ofY+ofS+(N-1)]=1e-3;
				if(x[ofY+ofS+(N-1)] > limS) x[ofY+ofS+(N-1)]=limS;

				x[ofY+ofR+(N-1)]=x[ofY+ofRm+(r1-1)]+gamma*(x[ofY+ofRm+(r2-1)]-x[ofY+ofRm+(r1-1)]);
				x[ofY+ofR+(r1-1)]=x[ofY+ofRm+(r1-1)];
				x[ofY+ofR+(r2-1)]=x[ofY+ofRm+(r2-1)];
				x[ofY+ofT+(N-1)]=1;
				x[ofY+ofT+(r1-1)]=1;
				x[ofY+ofT+(r2-1)]=1;

				x[ofY+ofCn+(r1-1)+(r2-1)*Nmax]=0;
				x[ofY+ofCn+(r2-1)+(r1-1)*Nmax]=0;
				x[ofY+ofCn+(N-1)+(r1-1)*Nmax]=1;
				x[ofY+ofCn+(r1-1)+(N-1)*Nmax]=1;
				x[ofY+ofCn+(N-1)+(r2-1)*Nmax]=1;
				x[ofY+ofCn+(r2-1)+(N-1)*Nmax]=1;
			}

			/* WEIGHTS ADAPTATION */

			if (lga!=0) {
				/* linear case */

				r2=1;
				dmn=1e37;
				for (i=1;i<=N;i++)
					if ( x[ofY+ofCn+(bmu-1)+(i-1)*Nmax]>0 )
						if ( D[i-1]<dmn ) {
							dmn=D[i-1];r2=i;
						}		

				d=0;
				for (j=1;j<=Ni;j++)
					d += (x[ofY+ofC+(j-1)*Nmax+(bmu-1)]-x[ofY+ofC+(j-1)*Nmax+(r2-1)]) * (x[ofY+ofC+(j-1)*Nmax+(bmu-1)]-x[ofY+ofC+(j-1)*Nmax+(r2-1)]);
				d=sqrt(d);

				if (dmn<d) {
					activ=D[bmu-1]/(D[bmu-1]+dmn);

					x[ofY+ofW+bmu-1] += etaW*(1-activ)*(*e[k]);
					if(x[ofY+ofW+bmu-1] < -limW) x[ofY+ofW+bmu-1]=-limW;
					if(x[ofY+ofW+bmu-1] >  limW) x[ofY+ofW+bmu-1]= limW;

					x[ofY+ofW+r2-1] += etaW*activ*(*e[k]);
					if(x[ofY+ofW+r2-1] < -limW) x[ofY+ofW+r2-1]=-limW;
					if(x[ofY+ofW+r2-1] >  limW) x[ofY+ofW+r2-1]= limW;
				}


			}
			else {
				/* gaussian case */

				for(j=1;j<=N;j++) {
					activ=exp( -(D[j-1]/(2*x[ofY+ofS+j-1]))*(D[j-1]/(2*x[ofY+ofS+j-1]))  );

					x[ofY+ofW+j-1] += etaW*activ*(*e[k]);
					if(x[ofY+ofW+j-1] < -limW) x[ofY+ofW+j-1]=-limW;
					if(x[ofY+ofW+j-1] >  limW) x[ofY+ofW+j-1]= limW;

					x[ofY+ofS+j-1] += etaS*(*e[k])*(x[ofY+ofW+j-1]*activ*D[j-1]*D[j-1])/(x[ofY+ofS+j-1]*x[ofY+ofS+j-1]*x[ofY+ofS+j-1]);
					if(x[ofY+ofS+j-1] < 1e-3) x[ofY+ofS+j-1]=1e-3;
					if(x[ofY+ofS+j-1] > limS) x[ofY+ofS+j-1]=limS;
				}
			}


			x[ofY+ofct]=conta;
			x[ofY+ofN]=N;

		} /* close cycle trough No */
	} /* close if LE!=0 */
} /* close mdloutput */

/* mdlTerminate - called when the simulation is terminated ***********************************/
static void mdlTerminate(SimStruct *S) {}

/* Trailer information to set everything up for simulink usage *******************************/
#ifdef  MATLAB_MEX_FILE                      /* Is this file being compiled as a MEX-file?   */
#include "simulink.c"                        /* MEX-file interface mechanism                 */
#else
#include "cg_sfun.h"                         /* Code generation registration function        */
#endif

#undef DIM
