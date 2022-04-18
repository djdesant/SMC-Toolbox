#define S_FUNCTION_NAME  differentiator_sFunction
#define S_FUNCTION_LEVEL 2
#include <mex.h>
#include "matrix.h"
#include "math.h"
#include "blas.h"
#include "simstruc.h"

#define n(S)  (int)*mxGetPr(ssGetSFcnParam(S, 0))
#define nf(S)  (int)*mxGetPr(ssGetSFcnParam(S, 1))
#define d(S)  (double)*mxGetPr(ssGetSFcnParam(S, 2))
#define r(S)  (double)*mxGetPr(ssGetSFcnParam(S, 3))
#define m(S)  (int)*mxGetPr(ssGetSFcnParam(S, 4))
#define mu(S)  (double)*mxGetPr(ssGetSFcnParam(S, 5))
#define Ts(S)  (double)*mxGetPr(ssGetSFcnParam(S, 6))
#define tau(S)  *(double*)ssGetDWork(S, 2)

static void CheckParameters(SimStruct *S) {
    if(mxIsNaN(n(S)) || n(S) <= 0)   mexErrMsgIdAndTxt("Differentiator:parameterError", "n must be an integer value greater than 0.");
    if(mxIsNaN(nf(S)) || nf(S) < 0)  mexErrMsgIdAndTxt("Differentiator:parameterError", "nf must be an integer value greater than or equal to 0.");
    if(mxIsNaN(d(S)) || d(S) > 0 || d(S) < -1) mexErrMsgIdAndTxt("Parameter d", "d must be a value between -1 and 0.");
    if(mxIsNaN(r(S)) || r(S) < 0)    mexErrMsgIdAndTxt("Differentiator:parameterError", "r must be a positive value.");
    if(mxIsNaN(m(S)) || (m(S) != 0 && m(S) != 1 && m(S) != 2)) mexErrMsgIdAndTxt("Differentiator:parameterError", "m must be out of the set {0,1,2}.");
    if(mxIsNaN(mu(S)) || mu(S) < 0)  mexErrMsgIdAndTxt("Differentiator:parameterError", "mu must be a positive value.");
    if(mxIsNaN(Ts(S)) || (Ts(S) <= 0 && Ts(S) != -1)) mexErrMsgIdAndTxt("Differentiator:parameterError", "Ts must be a positive value or -1 for inherited step size.");
    
    if(n(S) + nf(S) > 10) mexErrMsgIdAndTxt("Differentiator:parameterError", "The sum n+nf must not exceed 10.");
    if(m(S) == 2 && d(S) != -1) mexErrMsgIdAndTxt("Differentiator:parameterError", "Uniform RED (m = 2) is only allowed with RED (d = -1).");
}
static void CheckTau(SimStruct *S) {
    if(mxIsNaN(tau(S)) || tau(S) <= 0) mexErrMsgIdAndTxt("Differentiator:parameterError", "The simulation solver type must be FixedStep, if step size is inherited (Ts = -1).");
}

static void multMat(double* C, const double* A, const double* B, ptrdiff_t rowsA, ptrdiff_t colsA, ptrdiff_t colsB) {
    double one = 1.0, zero = 0.0;
    char *chn = (char*)"N";
    dgemm(chn, chn, &rowsA, &colsB, &colsA, &one, A, &rowsA, B, &colsA, &zero, C, &rowsA);
}

static real_T err(const int nf, const real_T u, const real_T first_state) {
    real_T x0;
    if(nf > 0) {
        x0 = -first_state;
    } else {
        x0 = u - first_state;
    }
    return x0;
}

static real_T cont_eigenvalues(const real_T x0, const int sys_n, const real_T d, const real_T r) {
    return -r * pow(fabs(x0), d / (1-d*(sys_n-1)));
}

static real_T disc_eigenvalues(const int m, const real_T s, const real_T Ts) {
    real_T z = 0;
    if(m == 0) {
        if(mxIsNaN(s)) {
            z = 0;
        } else {
            z = 1 + Ts*s;
        }
    } else if(m == 1) {
        z = exp(Ts*s);
    }
    return z;
}

static real_T disc_eigenvalues_URED(const real_T x0, const int sys_n, const real_T r, const real_T mu, const real_T Ts) {
    const real_T c = pow(fabs(x0), (1.0/sys_n));
    return c / (Ts * r * mu * fabs(x0) + c + Ts * r);
}

static ackerman_precomputed(const int sys_n, const double Ts, const double z, double *lambda) {
    const double z0 = z;
    double t2, t3, t4, t5, t6, t7, t8;
    
    switch (sys_n-1) {
        case 1:
            lambda[0] = -z-z0+2.0;
            lambda[1] = ((z-1.0)*(z0-1.0))/Ts;
            break;
        case 2:
            t2 = z-1.0;
            lambda[0] = z*-2.0-z0+3.0;
            lambda[1] = (t2*(z+z0*3.0+z*z0-5.0))/(Ts*2.0);
            lambda[2] = -1.0/(Ts*Ts)*(t2*t2)*(z0-1.0);
            break;
        case 3:
            t2 = z*z;
            t3 = z-1.0;
            t4 = t3*t3;
            lambda[0] = z*-3.0-z0+4.0;
            lambda[1] = (t3*(t2+z*7.0+z0*1.1E1+t2*z0*2.0+z*z0*5.0-2.6E1))/(Ts*6.0);
            lambda[2] = -1.0/(Ts*Ts)*t4*(z0*2.0+z*z0-3.0);
            lambda[3] = 1.0/(Ts*Ts*Ts)*t3*t4*(z0-1.0);
            break;
        case 4:
            t2 = z*z;
            t3 = z-1.0;
            t4 = t3*t3;
            lambda[0] = z*-4.0-z0+5.0;
            lambda[1] = (t3*(t2*5.0+z*2.3E1+z0*2.5E1+t2*z+t2*z0*7.0+z*z0*1.3E1+t2*z*z0*3.0-7.7E1))/(Ts*1.2E1);
            lambda[2] = 1.0/(Ts*Ts)*t4*(t2-z*2.0+z0*3.5E1+t2*z0*1.1E1+z*z0*2.6E1-7.1E1)*(-1.0/1.2E1);
            lambda[3] = 1.0/(Ts*Ts*Ts)*t3*t4*(z-z0*5.0-z*z0*3.0+7.0)*(-1.0/2.0);
            lambda[4] = -1.0/(Ts*Ts*Ts*Ts)*(t4*t4)*(z0-1.0);
            break;
        case 5:
            t2 = z*z;
            t3 = t2*t2;
            t4 = z-1.0;
            t5 = t4*t4;
            t6 = t5*t5;
            lambda[0] = z*-5.0-z0+6.0;
            lambda[1] = (t4*(t2*4.3E1+t3*3.0+z*1.63E2+z0*1.37E2+t2*z*1.3E1+t2*z0*4.7E1+t3*z0*1.2E1+z*z0*7.7E1+t2*z*z0*2.7E1-5.22E2))/(Ts*6.0E1);
            lambda[2] = 1.0/(Ts*Ts)*t5*(t2*2.0-z*7.0+z0*4.5E1+t2*z+t2*z0*2.5E1+z*z0*4.0E1+t2*z*z0*1.0E1-1.16E2)*(-1.0/1.2E1);
            lambda[3] = 1.0/(Ts*Ts*Ts)*t4*t5*(t2+z*8.0-z0*1.7E1-t2*z0*7.0-z*z0*1.6E1+3.1E1)*(-1.0/4.0);
            lambda[4] = 1.0/(Ts*Ts*Ts*Ts)*t6*(z-z0*3.0-z*z0*2.0+4.0);
            lambda[5] = 1.0/(Ts*Ts*Ts*Ts*Ts)*t4*t6*(z0-1.0);
            break;
        case 6:
            t2 = z*z;
            t3 = t2*t2;
            t4 = z-1.0;
            t5 = t4*t4;
            t6 = t5*t5;
            lambda[0] = z*-6.0-z0+7.0;
            lambda[1] = (t4*(t2*6.3E1+t3*8.0+z*2.13E2+z0*1.47E2+t2*z*2.3E1+t2*z0*5.7E1+t3*z*2.0+t3*z0*2.2E1+z*z0*8.7E1+t2*z*z0*3.7E1+t3*z*z0*1.0E1-6.69E2))/(Ts*6.0E1);
            lambda[2] = 1.0/(Ts*Ts)*t5*(t2*3.3E1+t3*1.3E1-z*2.32E2+z0*8.12E2+t2*z*3.8E1+t2*z0*5.97E2+t3*z0*1.37E2+z*z0*8.02E2+t2*z*z0*3.52E2-2.552E3)*(-1.0/1.8E2);
            lambda[3] = 1.0/(Ts*Ts*Ts)*t4*t5*(t2*9.0+z*3.9E1-z0*4.9E1+t2*z-t2*z0*3.9E1-z*z0*5.7E1-t2*z*z0*1.5E1+1.11E2)*(-1.0/8.0);
            lambda[4] = (1.0/(Ts*Ts*Ts*Ts)*t6*(t2*5.0+z*2.6E1-z0*3.5E1-t2*z0*1.7E1-z*z0*3.8E1+5.9E1))/6.0;
            lambda[5] = 1.0/(Ts*Ts*Ts*Ts*Ts)*t4*t6*(z*3.0-z0*7.0-z*z0*5.0+9.0)*(-1.0/2.0);
            lambda[6] = -1.0/(Ts*Ts*Ts*Ts*Ts*Ts)*t5*t6*(z0-1.0);
            break;
        case 7:
            t2 = z*z;
            t3 = t2*t2;
            t4 = z-1.0;
            t5 = t4*t4;
            t6 = t5*t5;
            lambda[0] = z*-7.0-z0+8.0;
            lambda[1] = (t4*(t2*5.91E2+t3*1.01E2+z*1.851E3+z0*1.089E3+t2*t3*1.0E1+t2*z*2.41E2+t2*z0*4.59E2+t3*z*3.8E1+t3*z0*2.14E2+z*z0*6.69E2+t2*t3*z0*6.0E1+t2*z*z0*3.19E2+t3*z*z0*1.3E2-5.772E3))/(Ts*4.2E2);
            lambda[2] = 1.0/(Ts*Ts)*t5*(t2*1.6E1+t3*3.6E1-z*4.14E2+z0*9.38E2+t2*z*6.1E1+t2*z0*8.19E2+t3*z*1.1E1+t3*z0*3.29E2+z*z0*9.94E2+t2*z*z0*5.74E2+t3*z*z0*1.26E2-3.49E3)*(-1.0/1.8E2);
            lambda[3] = 1.0/(Ts*Ts*Ts)*t4*t5*(t2*3.57E2+t3*7.0+z*1.127E3-z0*9.67E2+t2*z*7.7E1-t2*z0*1.077E3-t3*z0*2.32E2-z*z0*1.277E3-t2*z*z0*6.47E2+2.632E3)*(-1.0/1.2E2);
            lambda[4] = (1.0/(Ts*Ts*Ts*Ts)*t6*(t2*2.3E1+z*6.8E1-z0*5.6E1+t2*z*4.0-t2*z0*5.6E1-z*z0*7.7E1-t2*z*z0*2.1E1+1.15E2))/6.0;
            lambda[5] = 1.0/(Ts*Ts*Ts*Ts*Ts)*t4*t6*(t2*1.0E1+z*4.3E1-z0*4.6E1-t2*z0*2.5E1-z*z0*5.5E1+7.3E1)*(-1.0/6.0);
            lambda[6] = 1.0/(Ts*Ts*Ts*Ts*Ts*Ts)*t5*t6*(z*2.0-z0*4.0-z*z0*3.0+5.0);
            lambda[7] = 1.0/(Ts*Ts*Ts*Ts*Ts*Ts*Ts)*t4*t5*t6*(z0-1.0);
            break;
        case 8:
            t2 = z*z;
            t3 = t2*t2;
            t4 = z-1.0;
            t5 = t4*t4;
            t6 = t5*t5;
            lambda[0] = z*-8.0-z0+9.0;
            lambda[1] = (t4*(t2*1.497E3+t3*3.07E2+z*4.437E3+z0*2.283E3+t2*t3*5.5E1+t2*z*6.57E2+t2*z0*1.023E3+t3*z*1.39E2+t3*z0*5.33E2+z*z0*1.443E3+t2*t3*z*1.5E1+t2*t3*z0*2.25E2+t2*z*z0*7.43E2+t3*z*z0*3.65E2+t2*t3*z*z0*1.05E2-1.3827E4))/(Ts*8.4E2);
            lambda[2] = 1.0/(Ts*Ts)*t5*(t2*-7.33E2+t3*1.787E3-z*1.8254E4+z0*2.9531E4+t2*t3*2.61E2+t2*z*2.172E3+t2*z0*2.9013E4+t3*z*8.98E2+t3*z0*1.5293E4+z*z0*3.2926E4+t2*t3*z0*3.267E3+t2*z*z0*2.2468E4+t3*z*z0*8.622E3-1.27251E5)*(-1.984126984126984E-4);
            lambda[3] = 1.0/(Ts*Ts*Ts)*t4*t5*(t2*1.462E3+t3*8.7E1+z*3.777E3-z0*2.403E3+t2*z*4.42E2-t2*z0*3.302E3+t3*z*5.0-t3*z0*1.367E3-z*z0*3.457E3-t2*z*z0*2.442E3-t3*z*z0*4.69E2+7.667E3)*(-1.0/2.4E2);
            lambda[4] = (1.0/(Ts*Ts*Ts*Ts)*t6*(t2*2.522E3+t3*1.27E2+z*5.572E3-z0*3.207E3+t2*z*7.72E2-t2*z0*4.682E3-t3*z0*9.67E2-z*z0*5.092E3-t2*z*z0*2.852E3+7.807E3))/2.4E2;
            lambda[5] = 1.0/(Ts*Ts*Ts*Ts*Ts)*t4*t6*(t2*5.0E1+z*1.22E2-z0*8.1E1+t2*z*1.0E1-t2*z0*9.5E1-z*z0*1.25E2-t2*z*z0*3.5E1+1.54E2)*(-1.0/6.0);
            lambda[6] = (1.0/(Ts*Ts*Ts*Ts*Ts*Ts)*t5*t6*(t2*1.1E1+z*4.2E1-z0*3.9E1-t2*z0*2.3E1-z*z0*5.0E1+5.9E1))/4.0;
            lambda[7] = 1.0/(Ts*Ts*Ts*Ts*Ts*Ts*Ts)*t4*t5*t6*(z*5.0-z0*9.0-z*z0*7.0+1.1E1)*(-1.0/2.0);
            lambda[8] = -1.0/(Ts*Ts*Ts*Ts*Ts*Ts*Ts*Ts)*(t6*t6)*(z0-1.0);
            break;
        case 9:
            t2 = z*z;
            t3 = t2*t2;
            t4 = t3*t3;
            t5 = z-1.0;
            t6 = t5*t5;
            t7 = t6*t6;
            t8 = t7*t7;
            lambda[0] = z*-9.0-z0+1.0E1;
            lambda[1] = (t5*(t2*5.471E3+t3*1.271E3+t4*3.5E1+z*1.5551E4+z0*7.129E3+t2*t3*3.05E2+t2*z*2.531E3+t2*z0*3.349E3+t3*z*6.41E2+t3*z0*1.879E3+t4*z0*2.8E2+z*z0*4.609E3+t2*t3*z*1.25E2+t2*t3*z0*9.55E2+t2*z*z0*2.509E3+t3*z*z0*1.375E3+t2*t3*z*z0*5.95E2-4.861E4))/(Ts*2.52E3);
            lambda[2] = 1.0/(Ts*Ts)*t6*(t2*-2.712E3+t3*2.566E3-z*2.6477E4+z0*3.2575E4+t2*t3*7.88E2+t2*z*2.321E3+t2*z0*3.4905E4+t3*z*1.677E3+t3*z0*2.1689E4+z*z0*3.7754E4+t2*t3*z*2.23E2+t2*t3*z0*8.095E3+t2*z*z0*2.8864E4+t3*z*z0*1.4514E4+t2*t3*z*z0*3.044E3-1.59826E5)*(-1.984126984126984E-4);
            lambda[3] = (1.0/(Ts*Ts*Ts)*t5*t6*(t2*-1.61922E5-t3*1.7337E4-z*3.63543E5+z0*1.8092E5+t2*t3*1.6E1-t2*z*6.0422E4+t2*z0*2.87607E5-t3*z*2.931E3+t3*z0*1.65702E5+z*z0*2.76981E5+t2*t3*z0*2.9531E4+t2*z*z0*2.40602E5+t3*z*z0*8.8737E4-6.63941E5))/1.512E4;
            lambda[4] = (1.0/(Ts*Ts*Ts*Ts)*t7*(t2*5.368E3+t3*6.38E2+z*9.853E3-z0*4.275E3+t2*z*2.198E3-t2*z0*7.938E3+t3*z*1.01E2-t3*z0*3.363E3-z*z0*7.488E3-t2*z*z0*6.108E3-t3*z*z0*1.068E3+1.2082E4))/2.4E2;
            lambda[5] = 1.0/(Ts*Ts*Ts*Ts*Ts)*t5*t7*(t2*3.534E3+t3*2.29E2+z*6.428E3-z0*3.013E3+t2*z*1.244E3-t2*z0*5.334E3-t3*z0*1.069E3-z*z0*5.444E3-t2*z*z0*3.284E3+6.709E3)*(-1.0/1.44E2);
            lambda[6] = (1.0/(Ts*Ts*Ts*Ts*Ts*Ts)*t6*t7*(t2*6.0E1+z*1.29E2-z0*7.5E1+t2*z*1.3E1-t2*z0*9.9E1-z*z0*1.26E2-t2*z*z0*3.6E1+1.34E2))/4.0;
            lambda[7] = 1.0/(Ts*Ts*Ts*Ts*Ts*Ts*Ts)*t5*t6*t7*(t2*4.9E1+z*1.72E2-z0*1.45E2-t2*z0*9.1E1-z*z0*1.96E2+2.11E2)*(-1.0/1.2E1);
            lambda[8] = 1.0/(Ts*Ts*Ts*Ts*Ts*Ts*Ts*Ts)*t8*(z*3.0-z0*5.0-z*z0*4.0+6.0);
            lambda[9] = 1.0/(Ts*Ts*Ts*Ts*Ts*Ts*Ts*Ts*Ts)*t5*t8*(z0-1.0);
            break;
        case 10:
            t2 = z*z;
            t3 = t2*t2;
            t4 = t3*t3;
            t5 = z-1.0;
            t6 = t5*t5;
            t7 = t6*t6;
            t8 = t7*t7;
            lambda[0] = z*-1.0E1-z0+1.1E1;
            lambda[1] = (t5*(t2*6.479E3+t3*1.649E3+t4*9.8E1+z*1.7819E4+z0*7.381E3+t2*t3*4.73E2+t2*z*3.119E3+t2*z0*3.601E3+t3*z*8.93E2+t3*z0*2.131E3+t4*z*2.8E1+t4*z0*5.32E2+z*z0*4.861E3+t2*t3*z*2.33E2+t2*t3*z0*1.207E3+t2*z*z0*2.761E3+t3*z*z0*1.627E3+t4*z*z0*2.52E2+t2*t3*z*z0*8.47E2-5.5991E4))/(Ts*2.52E3);
            lambda[2] = 1.0/(Ts*Ts)*t6*(t2*-2.7739E4+t3*1.6165E4+t4*9.62E2-z*1.81196E5+z0*1.77133E5+t2*t3*7.611E3+t2*z*1.0278E4+t2*z0*2.02949E5+t3*z*1.2728E4+t3*z0*1.40985E5+t4*z0*1.4258E4+z*z0*2.11686E5+t2*t3*z*3.454E3+t2*t3*z0*6.8899E4+t2*z*z0*1.75852E5+t3*z*z0*1.04102E5+t2*t3*z*z0*3.8136E4-9.76263E5)*(-3.968253968253968E-5);
            lambda[3] = (1.0/(Ts*Ts*Ts)*t5*t6*(t2*-5.14467E5-t3*8.0075E4-z*1.040321E6+z0*4.20475E5-t2*t3*2.669E3-t2*z*2.22035E5+t2*z0*7.44585E5-t3*z*2.1303E4+t3*z0*5.26605E5+z*z0*6.75075E5+t2*t3*z*4.27E2+t2*t3*z0*1.80175E5+t2*z*z0*6.76405E5+t3*z*z0*3.46845E5+t2*t3*z*z0*5.8635E4-1.748357E6))/3.024E4;
            lambda[4] = (1.0/(Ts*Ts*Ts*Ts)*t7*(t2*6.1743E5+t3*1.18155E5+z*9.94506E5-z0*3.41693E5+t2*t3*5.084E3+t2*z*3.0404E5-t2*z0*7.5099E5+t3*z*3.3126E4-t3*z0*4.62765E5-z*z0*6.43092E5-t2*t3*z0*7.2368E4-t2*z*z0*6.6566E5-t3*z*z0*2.38632E5+1.102859E6))/1.512E4;
            lambda[5] = 1.0/(Ts*Ts*Ts*Ts*Ts)*t5*t7*(t2*1.601E4+t3*2.445E3+z*2.4135E4-z0*8.591E3+t2*z*7.55E3-t2*z0*1.983E4+t3*z*4.27E2-t3*z0*8.555E3-z*z0*1.7305E4-t2*z*z0*1.573E4-t3*z*z0*2.565E3+2.2009E4)*(-1.0/2.88E2);
            lambda[6] = (1.0/(Ts*Ts*Ts*Ts*Ts*Ts)*t6*t7*(t2*1.1418E4+t3*8.53E2+z*1.8188E4-z0*7.513E3+t2*z*4.388E3-t2*z0*1.5378E4-t3*z0*3.013E3-z*z0*1.4948E4-t2*z*z0*9.548E3+1.5553E4))/2.4E2;
            lambda[7] = 1.0/(Ts*Ts*Ts*Ts*Ts*Ts*Ts)*t5*t6*t7*(t2*5.81E2+z*1.139E3-z0*6.05E2+t2*z*1.33E2-t2*z0*8.75E2-z*z0*1.085E3-t2*z*z0*3.15E2+1.027E3)*(-1.0/2.4E1);
            lambda[8] = (1.0/(Ts*Ts*Ts*Ts*Ts*Ts*Ts*Ts)*t8*(t2*1.7E1+z*5.6E1-z0*4.4E1-t2*z0*2.9E1-z*z0*6.2E1+6.2E1))/3.0;
            lambda[9] = 1.0/(Ts*Ts*Ts*Ts*Ts*Ts*Ts*Ts*Ts)*t5*t8*(z*7.0-z0*1.1E1-z*z0*9.0+1.3E1)*(-1.0/2.0);
            lambda[10] = -1.0/pow(Ts,1.0E1)*t6*t8*(z0-1.0);
            break;
        default:
            int i;
            for (i=0; i<sys_n; i++) {
                lambda[i] = 0;
            } break;
    }
}

static void step(real_T* z, const double x0, const double u, const int sys_n, const double *Phi, const double *lambda, const double *bD) {
    int i;
    
    double z_new[sys_n];
    multMat(z_new,  Phi, z, sys_n, sys_n, 1);
    for (i=0; i<sys_n; i++) {
        z[i] = z_new[i] + bD[i] * u + lambda[i] * x0;
    }
}

static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumSFcnParams(S, 7);
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
        return;
    }
    CheckParameters(S);
    
    const int n = n(S);
    const int nf = nf(S);
    const int sys_n = n+nf+1;
    
    ssSetNumContStates(S, 0);
    ssSetNumDiscStates(S, sys_n);

    if (!ssSetNumInputPorts(S, 1)) return;
    ssSetInputPortWidth(S, 0, 1);
    ssSetInputPortDirectFeedThrough(S, 0, 1);

    if (!ssSetNumOutputPorts(S, 3)) return;
    ssSetOutputPortWidth(S, 0, 1);
    ssSetOutputPortWidth(S, 1, 1);
    ssSetOutputPortWidth(S, 2, n);

    ssSetNumSampleTimes(S, 1);
    
    ssSetNumDWork(S, 3);
    ssSetDWorkDataType(S, 0, SS_DOUBLE);
    ssSetDWorkDataType(S, 1, SS_DOUBLE);
    ssSetDWorkDataType(S, 2, SS_DOUBLE);
    ssSetDWorkUsageType(S, 0, SS_DWORK_USED_AS_DWORK);
    ssSetDWorkUsageType(S, 1, SS_DWORK_USED_AS_DWORK);
    ssSetDWorkUsageType(S, 2, SS_DWORK_USED_AS_DWORK);
    ssSetDWorkWidth(S, 0, sys_n*sys_n);
    ssSetDWorkWidth(S, 1, sys_n);
    ssSetDWorkWidth(S, 2, 1);

    ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE);
}

static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, Ts(S));
    ssSetOffsetTime(S, 0, 0.0);
    ssSetModelReferenceSampleTimeDefaultInheritance(S); 
}

#define MDL_INITIALIZE_CONDITIONS
static void mdlInitializeConditions(SimStruct *S)
{
    const int n = n(S);
    const int nf = nf(S);
    const int sys_n = n+nf+1;
    if(Ts(S) == -1) {
        ((double*)ssGetDWork(S, 2))[0] = ssGetFixedStepSize(S);
        CheckTau(S);
    } else {
        ((double*)ssGetDWork(S, 2))[0] = Ts(S);
    }
    const double Ts = tau(S);
    
    int_T i, j;
    
    //Phi
    double *Phi = (double*)ssGetDWork(S, 0);
    int diff = 0;
    double prev = 0;
    double val = 0;
    for (i=0; i<sys_n; i++) {
        for (j=0; j<sys_n; j++) {
            diff = j - i;
            if(diff > 0) {
                val = prev * Ts / (double)diff;
            } else if(diff == 0) {
                val = 1;
            } else {
                val = 0;
            }
            
            Phi[i + j*(sys_n)] = val;
            prev = val;
        }
    }
    
    //bD
    prev = 0;
    val = 0;
    double *bD = (double*)ssGetDWork(S, 1);
    for (i=sys_n-1; i>=0; i--) {
        diff = nf - i;
        if(diff > 0) {
            val = prev * Ts / (double)diff;
            prev = val;
        } else if(diff == 0) {
            prev = -1;
        }

        bD[i] = val;
    }
    
    //init states
    real_T *z0 = ssGetRealDiscStates(S);
    int_T  lp;
    for (lp=0;lp<sys_n;lp++) { 
        *z0++=0.0;
    }
}

#define MDL_OUTPUTS
static void mdlOutputs(SimStruct *S, int_T tid)
{
    const int n = n(S);
    const int nf = nf(S);
    
    InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S,0);
    const real_T u = *(uPtrs[0]);
    real_T *outX0 = ssGetOutputPortRealSignal(S,0);
    real_T *outZ0 = ssGetOutputPortRealSignal(S,1);
    real_T *outZ = ssGetOutputPortRealSignal(S,2);
    real_T *z = ssGetRealDiscStates(S);

    *outX0 = err(nf, u, z[0]);
    *outZ0 = z[nf];
    int_T i;
    for (i=0; i<n; i++) {
        *outZ++ = z[nf+i+1];
    }
}

static void mdlTerminate(SimStruct *S)
{
    
}

#define MDL_UPDATE
static void mdlUpdate(SimStruct *S, int_T tid)
{
    const int n = n(S);
    const int nf = nf(S);
    const int sys_n = n+nf+1;
    const real_T d = d(S);
    const real_T r = r(S);
    const int m = m(S);
    const real_T mu = mu(S);
    const double Ts = tau(S);
    
    const double *Phi = (double*)ssGetDWork(S, 0);
    const double *bD = (double*)ssGetDWork(S, 1);
    
    InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S,0);
    const real_T u = *(uPtrs[0]);
    real_T *z = ssGetRealDiscStates(S);
    
    const double x0 = err(nf, u, z[0]);
    double z_eig = 0;
    if(m == 2) {
        z_eig = disc_eigenvalues_URED(x0, sys_n, r, mu, Ts);
    } else {
        const double s = cont_eigenvalues(x0, sys_n, d, r);
        z_eig = disc_eigenvalues(m, s, Ts);
    }
    
    double lambda[sys_n];
    ackerman_precomputed(sys_n, Ts, z_eig, lambda);
    step(z, x0, u, sys_n, Phi, lambda, bD);
}

#ifdef MATLAB_MEX_FILE
#include "simulink.c"
#else
#include "cg_sfun.h"
#endif