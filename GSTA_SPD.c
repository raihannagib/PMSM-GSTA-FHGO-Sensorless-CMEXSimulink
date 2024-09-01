#define S_FUNCTION_LEVEL 2 
#define S_FUNCTION_NAME GSTA_SPD
#include "simstruc.h" 
#include <math.h> 

#define U(element) (*uPtrs[element]) /*Pointer to Input Port0*/ 

static void mdlInitializeSizes(SimStruct *S){ 
	ssSetNumDiscStates(S, 5); // 5 DISCRETE STATE USED
	if (!ssSetNumInputPorts(S, 1)) return; 
	ssSetInputPortWidth(S, 0, 6); // 6 INPUT
	ssSetInputPortDirectFeedThrough(S, 0, 1); 
	ssSetInputPortOverWritable(S, 0, 1); 
	if (!ssSetNumOutputPorts(S, 1)) return; 
	ssSetOutputPortWidth(S, 0, 2); // 2 OUTPUT
	ssSetNumSampleTimes(S, 1); 

	ssSetOptions(S, (SS_OPTION_EXCEPTION_FREE_CODE 
	| SS_OPTION_DISCRETE_VALUED_OUTPUT));
} 

static void mdlInitializeSampleTimes(SimStruct *S){ 
	ssSetSampleTime(S, 0, 1e-3); 
	ssSetOffsetTime(S, 0, 0.0);
} 

#define MDL_INITIALIZE_CONDITIONS 
static void mdlInitializeConditions(SimStruct *S){ 
	real_T *X0 = ssGetRealDiscStates(S); 
	int_T nXStates = ssGetNumDiscStates(S); 
	InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S,0); 
	int_T i; 

	/* initialize the states to 0.0 */ 
	for (i=0; i < nXStates; i++) { 
	X0[i] = 0.0; 
	}
} 

static void mdlOutputs(SimStruct *S, int_T tid){ 
	real_T *Y = ssGetOutputPortRealSignal(S,0); 
	real_T *X = ssGetRealDiscStates(S); 
    
    real_T control_signal, sliding_surface;
    
    control_signal = X[0];
    sliding_surface = X[4];
    
    Y[0] = control_signal;
    
    Y[1] = sliding_surface;
} 

#define MDL_UPDATE 
static void mdlUpdate(SimStruct *S, int_T tid) {
	real_T *X = ssGetRealDiscStates(S); 
	InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S,0); 

	real_T dt = 1e-3;

    
    // PMSM'S MODEL PARAMETER
    // PMSM MODEL'S PARAMETER
    real_T N    = 4;
    real_T psi  = 0.121;
    real_T Lsd  = 16.61e-3;
    real_T Lsq  = 16.22e-3;
    real_T Rs   = 0.55;
    real_T J    = 0.01;
    real_T B    = 0.08;
    
    
    // PARAMETER OF SPEED TRACKING CONTROLLER
    real_T nw1, nw2, kwm, uw, eps_w;
    nw1 = 1000;
    nw2 = 100000;
    kwm = 100;
    uw = 10;
    eps_w = 100;
    
    // INPUT
    real_T alpha1, alpha2, m, speed_ref, theta_m_act, d_hat, speed_act;
    alpha1 = U(0);
	alpha2 = U(1);
    m = U(2);
    speed_ref = U(3);
    speed_act = U(4);
    d_hat = U(5);
    
    
    real_T control_signal;
    real_T vw, vw_dot; 
    real_T psi_1, psi_2;
    real_T kw1, kw1_dot;
    real_T kw2;
    real_T s_w;
    real_T sign_sw;
    real_T sign_sw_uw;
    
    
    // STATE
    real_T kw1_old, vw_old, theta_m_old, speed_ref_old, error_integral_old;
    kw1_old = X[1];
    vw_old = X[2];
    speed_ref_old = X[3];
    
    
    // SPEED TRACKING ALGORITHM
	s_w = speed_act - speed_ref; 
   
    
    if ((fabs(s_w) - uw) > 0) {
        sign_sw_uw = 1;
    }
    
    else if ((fabs(s_w) - uw) < 0) {
        sign_sw_uw = -1;
    }
    
    else {
        sign_sw_uw = 0;
    }
    
    if (fabs(kw1_old) > kwm) {
        kw1_dot = nw1 * sign_sw_uw;
    }
    
    else {
        kw1_dot = nw2;
    }
    
    kw1 = kw1_old + kw1_dot * dt;
    
    kw2 = 2 * eps_w * kw1;
    
    if (s_w > 0) {
        sign_sw = 1;
    }
    
    else if (s_w < 0) {
        sign_sw = -1;
    }
    
    else {
        sign_sw = 0;
    }
    
    psi_1 = alpha1 * pow(fabs(s_w), (m - 1) / m) * sign_sw + alpha2 * s_w;
    psi_2 = ((m - 1) / m) * pow(alpha1, 2) * pow(fabs(s_w), (m - 2) / m) * sign_sw + 
            ((2 * m - 1) / m) * alpha1 * alpha2 * pow(fabs(s_w), (m - 1) / m) * sign_sw + pow(alpha2, 2) * s_w;

    vw_dot = -kw2 * psi_2;
    
    vw = vw_old + vw_dot * dt;
    
    control_signal = ( (2 * J) / (3 * N * psi) ) * ( (B * speed_act) / J + ((speed_ref - speed_ref_old) / dt) - 0 - kw1 * psi_1 + vw );
    
    // Update State
    X[0] = control_signal;
    X[1] = kw1;
    X[2] = vw;
    X[3] = speed_ref;
    X[4] = s_w;
}

static void mdlTerminate(SimStruct *S) 
{ } /*Keep this function empty since no memory is allocated*/ 

#ifdef MATLAB_MEX_FILE 
/* Is this file being compiled as a MEX-file? */ 
#include "simulink.c" /*MEX-file interface mechanism*/ 
#else 
#include "cg_sfun.h" /*Code generation registration function*/ 
#endif
