#define S_FUNCTION_LEVEL 2 
#define S_FUNCTION_NAME FHGO
#include "simstruc.h" 
#include <math.h> 

#define U(element) (*uPtrs[element]) /*Pointer to Input Port0*/ 

static void mdlInitializeSizes(SimStruct *S){ 
	ssSetNumContStates(S, 4); // 4 CONTINUOUS STATE
	if (!ssSetNumInputPorts(S, 1)) return; 
	ssSetInputPortWidth(S, 0, 2); // 2 INPUT
	ssSetInputPortDirectFeedThrough(S, 0, 1); 
	ssSetInputPortOverWritable(S, 0, 1); 
	if (!ssSetNumOutputPorts(S, 1)) return; 
	ssSetOutputPortWidth(S, 0, 2); // 2 OUTPUT
	ssSetNumSampleTimes(S, 1); 
	ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE); 
} 

static void mdlInitializeSampleTimes(SimStruct *S) {
	ssSetSampleTime(S, 0, CONTINUOUS_SAMPLE_TIME); 
	ssSetOffsetTime(S, 0, 0.0);
} 

#define MDL_INITIALIZE_CONDITIONS 
static void mdlInitializeConditions(SimStruct *S) { 

	real_T *X0 = ssGetContStates(S); 
	int_T nStates = ssGetNumContStates(S); 
	int_T i; 

	/* initialize the states to 0.0 */ 
	for (i=0; i < nStates; i++) {
		X0[i] = 0.0;
	}
} 

static void mdlOutputs(SimStruct *S, int_T tid) { 
	real_T *Y = ssGetOutputPortRealSignal(S,0); 
	real_T *X = ssGetContStates(S); 
	InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S,0); 
	
    real_T omega_hat, d_hat;
    real_T zeta1, zeta2;
    
     // STATE VARIABLE
    omega_hat = X[0];
    d_hat = X[1];
    zeta1 = X[2];
    zeta2 = X[3];
   
    Y[0] = d_hat;
    Y[1] = omega_hat;
} 

#define MDL_DERIVATIVES 
static void mdlDerivatives(SimStruct *S) { 
	
	real_T *dX = ssGetdX(S); 
	real_T *X = ssGetContStates(S); 
	InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S,0); 
	
	// PMSM MODEL'S PARAMETER
    real_T N    = 4;
    real_T psi  = 0.121;
    real_T Lsd  = 16.61e-3;
    real_T Lsq  = 16.22e-3;
    real_T Rs   = 0.55;
    real_T J    = 0.01;
    real_T B    = 0.08;
    
    
    real_T omega_hat_dot, d_hat_dot;
    real_T omega_hat, d_hat;
    real_T zeta1_dot, zeta2_dot;
    real_T zeta1, zeta2;
    real_T omega_dot, omega, d;
    
    
    // FHGO parameter
    real_T theta, k01, k02, Dn1, Dn2;
    theta = 20;
    k01 = 5;
    k02 = 200;
    Dn1 = -70;
    Dn2 = 60;
    
    
    // INPUT
    real_T omega_act, Isq;
	omega_act = U(0);
    Isq = U(1);
    
    // STATE VARIABLE
    omega_hat = X[0];
    d_hat = X[1];
    zeta1 = X[2];
    zeta2 = X[3];
    
    real_T mt, nt;
    nt = 1;
    mt = omega_act - nt;
    
    
    omega_hat_dot = (-B / J) * omega_hat + d_hat + (3 * N * psi / (2 * J)) * Isq - theta * k01 * zeta1;
    d_hat_dot = -theta * k02 * zeta2;
    
    zeta1_dot = -theta * Dn1 * zeta1 - (pow(theta, 2) * B * zeta1) / J + theta * (omega_hat - mt);
    zeta2_dot = -theta * Dn2 * zeta2 + pow(theta, 2) * zeta1;
    
    // State Derivatives
    dX[0] = omega_hat_dot;
    dX[1] = d_hat_dot;
    dX[2] = zeta1_dot;
    dX[3] = zeta2_dot;
} 

static void mdlTerminate(SimStruct *S) 
{} /*Keep this function empty since no memory is allocated*/ 

#ifdef MATLAB_MEX_FILE 
/* Is this file being compiled as a MEX-file? */ 
#include "simulink.c" /* MEX-file interface mechanism */ 
#else 
#include "cg_sfun.h" /*Code generation registration function*/ 
#endif 
