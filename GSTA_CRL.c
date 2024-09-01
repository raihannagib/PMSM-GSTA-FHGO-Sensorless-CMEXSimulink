#define S_FUNCTION_LEVEL 2 
#define S_FUNCTION_NAME GSTA_CRL
#include "simstruc.h" 
#include <math.h> 

#define U(element) (*uPtrs[element]) /*Pointer to Input Port0*/ 

static void mdlInitializeSizes(SimStruct *S){
	ssSetNumDiscStates(S, 8); // 8 DISCRETE STATE USED
	if (!ssSetNumInputPorts(S, 1)) return; 
	ssSetInputPortWidth(S, 0, 9); // 9 INPUT
	ssSetInputPortDirectFeedThrough(S, 0, 1); 
	ssSetInputPortOverWritable(S, 0, 1); 
	if (!ssSetNumOutputPorts(S, 1)) return; 
	ssSetOutputPortWidth(S, 0, 3); // 3 OUTPUT
	ssSetNumSampleTimes(S, 1); 

	ssSetOptions(S, (SS_OPTION_EXCEPTION_FREE_CODE 
	| SS_OPTION_DISCRETE_VALUED_OUTPUT));
	
} 

static void mdlInitializeSampleTimes(SimStruct *S){ 
	ssSetSampleTime(S, 0, 1e-4); 
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
	
    // OUTPUT
    real_T Va, Vb, Vc;
    
    Va = X[5];
    Vb = X[6];
    Vc = X[7];
    
    Y[0] = Va;
    Y[1] = Vb;
    Y[2] = Vc;

} 

#define MDL_UPDATE 
static void mdlUpdate(SimStruct *S, int_T tid) {
 
	real_T *X = ssGetRealDiscStates(S); 
	InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S,0); 

	real_T dt = 1e-4;
	
	// PMSM MODEL'S PARAMETER
    real_T N    = 4;
    real_T psi  = 0.121;
    real_T Lsd  = 0.01661;
    real_T Lsq  = 0.01622;
    real_T Rs   = 0.55;
    real_T J    = 0.01;
    real_T B    = 0.08;
    
    /*
    real_T N    = 4;
    real_T psi  = 0.432;
    real_T Lsd  = 0.00932;
    real_T Lsq  = 0.001414;
    real_T Rs   = 0.602;
    real_T J    = 0.007;
    real_T B    = 0.08;
    */
    // PARAMETER OF CURRENT REGULATION CONTROLLER
    real_T ni1, ni2, kim, ui, eps_i;
    ni1 = 10000;
    ni2 = 1000000;
    kim = 100;
    ui = 0.5;
    eps_i = 20;
    
    real_T Ia, Ib, Ic;
    real_T I_alpha, I_beta;
    real_T Isd_act, Isq_act;
    
    real_T Va, Vb, Vc;
    real_T V_alpha, V_beta;
    real_T Vsd, Vsq;
    
    real_T pi = 3.141592654;

    real_T K = 0.8164965809; // akar(2/3)
    real_T L = 0.8660254038; // akar(3/2)
	
    
    // INPUT
    real_T Isd_ref, Isq_ref, alpha1, alpha2, m, theta_e_act;
	// Current feedback
	Ia = U(0);
    Ib = U(1);
    Ic = U(2);
    
    // Current reference
    Isd_ref = U(3);
    Isq_ref = U(4);
    
    // Controller Parameter
    alpha1 = U(5);
    alpha2 = U(6);
    m = U(7);
    
    // Position Feedback
    theta_e_act = U(8);
    
	// TRANSFORM ABC TO ALPHA-BETA
	I_alpha = K * (Ia - 0.5 * Ib - 0.5 * Ic);
    I_beta = K * L * (Ib - Ic);
	
	// TRANSFORM ALPHA-BETA TO DQ
	Isd_act = cos(theta_e_act) * I_alpha + sin(theta_e_act) * I_beta;
    Isq_act = -sin(theta_e_act) * I_alpha + cos(theta_e_act) * I_beta;
    
    
    real_T error_Isd, error_Isq;
    real_T theta_e_old;
    real_T omega_e;
    real_T k_Isd1_old, k_Isd1_dot, k_Isd1;
    real_T k_Isd2;
    real_T k_Isq1_old, k_Isq1_dot, k_Isq1;
    real_T k_Isq2;
    real_T psi_1_Isd;
    real_T psi_1_Isq;
    real_T psi_2_Isd_old, psi_2_Isd, psi_2_Isd_new;
    real_T psi_2_Isq_old, psi_2_Isq, psi_2_Isq_new;
    real_T ro_sd;
    real_T ro_sq;
    real_T sign_Isd, sign_Isd_ui;
    real_T sign_Isq, sign_Isq_ui;
    
    
    // STATE
    theta_e_old = X[0];
    k_Isd1_old = X[1];
    k_Isq1_old = X[2];
    psi_2_Isd_old = X[3];
    psi_2_Isq_old = X[4];
    
    
    // CURRENT REGULATION ALGORITHM  
    error_Isd = Isd_act - Isd_ref;
    error_Isq = Isq_act - Isq_ref;
    
    omega_e = (theta_e_act - theta_e_old) / ( dt * N );
    
    if ((fabs(error_Isd) - ui) > 0) {
        sign_Isd_ui = 1;
    }
    
    else if ((fabs(error_Isd) - ui) < 0) {
        sign_Isd_ui = -1;
    }
    
    else {
        sign_Isd_ui = 0;
    }
    
    if (k_Isd1_old > kim) {
        k_Isd1_dot = ni1 * sign_Isd_ui;
    }
    
    else {
        k_Isd1_dot = ni2;
    }
    
    k_Isd1 = k_Isd1_old + k_Isd1_dot * dt;
    
    k_Isd2 = 2 * eps_i * k_Isd1;
    
    if ((fabs(error_Isq) - ui) > 0) {
        sign_Isq_ui = 1;
    }
    
    else if ((fabs(error_Isq) - ui) < 0) {
        sign_Isq_ui = -1;
    }
    
    else {
        sign_Isq_ui = 0;
    } 
    
    if (k_Isq1_old > kim) {
        k_Isq1_dot = ni1 * sign_Isq_ui;
    }
    
    else {
        k_Isq1_dot = ni2;
    }
    
    k_Isq1 = k_Isq1_old + k_Isq1_dot * dt;
    
    k_Isq2 = 2 * eps_i * k_Isq1;
    
    if (error_Isd > 0) {
        sign_Isd = 1;
    }
    
    else if (error_Isd < 0) {
        sign_Isd = -1;
    }
    
    else {
        sign_Isd = 0;
    }
    
    if (error_Isq > 0) {
        sign_Isq = 1;
    }
    
    else if (error_Isq < 0) {
        sign_Isq = -1;
    }
    
    else {
        sign_Isq = 0;
    }
    
    psi_1_Isd =  alpha1 * pow(fabs(error_Isd), (m - 1) / m) * sign_Isd + alpha2 * error_Isd;
    psi_1_Isq =  alpha1 * pow(fabs(error_Isq), (m - 1) / m) * sign_Isq + alpha2 * error_Isq;
    
    psi_2_Isd = ((m - 1) / m) * pow(alpha1, 2) * pow(fabs(error_Isd), ((m - 2) / m)) * sign_Isd + 
            ((2 * m - 1) / m) * alpha1 * alpha2 * pow(fabs(error_Isd), ((m - 1) / m)) * sign_Isd + pow(alpha2, 2) * error_Isd;
    
    psi_2_Isq = ((m - 1) / m) * pow(alpha1, 2) * pow(fabs(error_Isq), ((m - 2) / m)) * sign_Isq + 
            ((2 * m - 1) / m) * alpha1 * alpha2 * pow(fabs(error_Isq), ((m - 1) / m)) * sign_Isq + pow(alpha2, 2) * error_Isq;
    
    ro_sd = k_Isd1 * psi_1_Isd + k_Isd2 * (psi_2_Isd_old + psi_2_Isd * dt);
    ro_sq = k_Isq1 * psi_1_Isq + k_Isq2 * (psi_2_Isq_old + psi_2_Isq * dt);
    
    psi_2_Isd_new = psi_2_Isd_old + psi_2_Isd * dt;
    psi_2_Isq_new = psi_2_Isq_old + psi_2_Isq * dt;
    
    Vsd = Lsd * ((-ro_sd) + (Rs / Lsd) * Isd_act + (omega_e * Lsq * Isq_act) / Lsd);
    Vsq = Lsq * ((-ro_sq) + (Rs / Lsq) * Isq_act + ((omega_e * (Lsd * Isd_act + psi)) / Lsq));
    
    // TRANSFORM DQ TO ALFA-BETA
    V_alpha = Vsd * cos(theta_e_act) - Vsq * sin(theta_e_act);
    V_beta = Vsd * sin(theta_e_act) + Vsq * cos(theta_e_act);
    
    // TRANSFORM ALFA-BETA TO ABC  
    Va = K * V_alpha;
    Vb = K * (-0.5 * V_alpha + L * V_beta);
    Vc = K * (-0.5 * V_alpha - L * V_beta);
    
    
    // UPDATE STATE
    X[0] = theta_e_act;
    X[1] = k_Isd1;
    X[2] = k_Isq1;
    X[3] = psi_2_Isd_new;
    X[4] = psi_2_Isq_new;
    X[5] = Va;
    X[6] = Vb;
    X[7] = Vc;
}

static void mdlTerminate(SimStruct *S) 
{ } /*Keep this function empty since no memory is allocated*/ 

#ifdef MATLAB_MEX_FILE 
/* Is this file being compiled as a MEX-file? */ 
#include "simulink.c" /*MEX-file interface mechanism*/ 
#else 
#include "cg_sfun.h" /*Code generation registration function*/ 
#endif
