#include <math.h>
#include "mex.h"


/* computes obj = sum(sum(X.*WH-expWH)) */

void mexFunction( int nlhs, mxArray *plhs[], 
                  int nrhs, const mxArray *prhs[] )
{ 

    unsigned int indn,indx,next_indx, indm,indnM;
    int *ir_x,*jc_x;
    double *X,*WH,*expWH,*objective;
    unsigned int M,N;


    if (!(nrhs==3)||!(mxIsSparse(prhs[0]))||mxIsSparse(prhs[1])){
	printf("usage : obj = epca_objective(X,WH,expWH)\n");
	printf("where : X is sparse, WH,expWH full but of same size as X\n");
	printf("result: obj = sum(sum(X.*WH-expWH));\n");
	return;
    }

    if ((mxGetN(prhs[0])!=mxGetN(prhs[1]))||((mxGetM(prhs[0])!=mxGetM(prhs[1])))){
	printf("dimension mismatch of input matrices\n");
	return;
    }


    jc_x = mxGetJc(prhs[0]);
    ir_x = mxGetIr(prhs[0]);
    X = mxGetPr(prhs[0]);
    WH = mxGetPr(prhs[1]);
    expWH = mxGetPr(prhs[2]);

    N = mxGetN(prhs[0]);
    M = mxGetM(prhs[0]);

    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    objective = mxGetPr(plhs[0]);

    indx = 0;

    for (indn = 0; indn < N; indn++){
	next_indx = jc_x[indn+1];/* next column index */
	indnM = indn*M;

	for (indm = 0; indm < M; indm++){

	    if (indx < next_indx){ /* while nonzero entries in the matrix row */
		objective[0] += X[indx] * WH[indnM+ir_x[indx]] - expWH[indnM+indm];
		indx++;
	    }
	    else{
		objective[0] -= expWH[indnM+indm];
	    }
	}
    }
    return;
}

