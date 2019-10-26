#include <math.h>
#include "mex.h"


void mexFunction( int nlhs, mxArray *plhs[], 
                  int nrhs, const mxArray *prhs[] )
{ 

    unsigned int indN,indx,next_indx,sparse_count;
    int *ir_x,*jc_x;
    double *X,*W,*varB, *H, av, *dW, *dH, xwh, expav;
    unsigned int M,N,D,indD, indM;
    char* distr;
    unsigned int getdH, getdW, ind1,ind2;
    

    if (((nrhs<3)||(nrhs>6))||!(mxIsSparse(prhs[0]))||mxIsSparse(prhs[1])||(mxIsSparse(prhs[2]))){
	printf("THIS IS INEFFICIENT COMPARED TO THE CURRENT IMLPEMENTATION OF\n");
	printf("IT IS ALSO NOT CALLED AND LEFT ONLY AS A SAMPLE\n");
	printf("usage : [obj,dW,dH] = epca_obj(X,W,H,distr,dw,dh)\n");
	printf("where : X is sparse [MxN], W,H full [MxD],[DxN] resp\n");
	printf("        distr is 'poisson', 'normal', 'bernoulli'\n");
	printf("        dw 1 if gradient for W should be calculated, dh same\n");
	printf("Examples:\n");
	printf("    [obj,dW] = epca_obj(X,W,H,distr,1,0);\n");
	printf("    [obj,dH] = epca_obj(X,W,H,distr,0,1);\n");
	printf("    [obj,dW,dH] = epca_obj(X,W,H,distr,1,1);\n");
	return;
    }

    if (mxGetN(prhs[0])!=mxGetM(prhs[1])){
	printf("dimension mismatch of 1st and 2nd argument\n");
	return;
    }

    if ((mxGetN(prhs[0])!=mxGetM(prhs[1]))||((mxGetM(prhs[0])!=mxGetN(prhs[2])))||(mxGetN(prhs[1])!=mxGetM(prhs[2]))){
	printf("dimension mismatch of input matrices\n");
	return;
    }

    getdW = getdH = 0; 
    if (nrhs>4) getdW = (int)(mxGetPr(prhs[4])[0]>0);
    if (nrhs>5) getdH = (int)(mxGetPr(prhs[5])[0]>0);

    /*    printf("getdW = %d\n",getdW);
	  printf("getdH = %d\n",getdH);*/

    if (getdW+getdH+1 != nlhs){
	printf("wrong number of output arguments\n");
	return;
    }


    distr = mxArrayToString(prhs[3]);

    if (!( (strcmp(distr,"poisson")==0)||(strcmp(distr,"normal")==0)||(strcmp(distr,"bernoulli")==0))){
	printf("error: unknown distribution\n");
	printf("distr is 'poisson', 'normal', 'bernoulli'\n");
	return;
    }

    if (getdW&getdH){
	plhs[1] = mxCreateDoubleMatrix(mxGetM(prhs[1]),mxGetN(prhs[1]),mxREAL);
	dW = mxGetPr(plhs[1]);
	plhs[2] = mxCreateDoubleMatrix(mxGetM(prhs[2]),mxGetN(prhs[2]),mxREAL);
	dH = mxGetPr(plhs[2]);
    }
    if (getdW&!getdH){
	plhs[1] = mxCreateDoubleMatrix(mxGetM(prhs[1]),mxGetN(prhs[1]),mxREAL);
	dW = mxGetPr(plhs[1]);
    }
    if (getdH&!getdW){
	plhs[1] = mxCreateDoubleMatrix(mxGetM(prhs[2]),mxGetN(prhs[2]),mxREAL);
	dH = mxGetPr(plhs[1]);
    }

    jc_x = mxGetJc(prhs[0]);
    ir_x = mxGetIr(prhs[0]);

    X = mxGetPr(prhs[0]);
    W = mxGetPr(prhs[1]);
    H = mxGetPr(prhs[2]);

	
    N = mxGetN(prhs[0]);
    M = mxGetM(prhs[0]);
    D = mxGetN(prhs[1]);

    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    varB = mxGetPr(plhs[0]);

    indx = 0;
    sparse_count = 0;
    next_indx = jc_x[0];


    printf("mxGetN(X) = %d\nmxGetM(X) = %d\n",mxGetN(prhs[0]),mxGetM(prhs[0]));
    printf("mxGetN(W) = %d\nmxGetM(W) = %d\n",mxGetN(prhs[1]),mxGetM(prhs[1]));
    printf("mxGetN(H) = %d\nmxGetM(H) = %d\n",mxGetN(prhs[2]),mxGetM(prhs[2]));

    if  (getdW)
	printf("mxGetN(dW) = %d\nmxGetM(dW) = %d\n",mxGetN(plhs[1]),mxGetM(plhs[1]));
    if (getdW & getdH){
	printf("am here\n");
	printf("mxGetN(dH) = %d\nmxGetM(dH) = %d\n",mxGetN(plhs[2]),mxGetM(plhs[2]));
    }
    if (getdH & !getdW)
	printf("mxGetN(dH) = %d\nmxGetM(dH) = %d\n",mxGetN(plhs[1]),mxGetM(plhs[1]));


    if (!getdH & !getdW){
	if (strcmp(distr,"poisson") == 0){
	    for (indN = 0; indN < N; indN++){
		next_indx = jc_x[indN+1]; /* next column index */

		for (indM = 0; indM<M; indM++){
		    av = 0;
		    
		    ind1 = indN;
		    for (ind2=indM*D;ind2<indM*D+D;ind2++){
			av += W[ind1] * H[ind2];
			ind1+=N;
		    }
		    /*for (indD =0; indD<D; indD++){
			av += W[indN + indD*N] * H[indD+indM*D];
			}*/
		    expav = exp(av);
		    if ((indx<next_indx)&&(ir_x[indx]==indM)){
			varB[0] += X[sparse_count] * av - expav;
			sparse_count++;
			indx++;
		    }
		    else
			varB[0] -= expav;
		}
	    }
	}else if (strcmp(distr,"normal") == 0 ){
	    for (indN = 0; indN < N; indN++){
		next_indx = jc_x[indN+1]; /* next column index */
		
		for (indM = 0; indM<M; indM++){
		    av = 0;
		    for (indD =0; indD<D; indD++){
			av += W[indN + indD*N] * H[indD+indM*D];
		    }
		    /* X(i,j) not empty */
		    if ((indx<next_indx)&&(ir_x[indx]==indM)){
			varB[0] += X[sparse_count] * av - .5*av*av;
			sparse_count++;
			indx++;
		    }
		    /* only exp(W*H) */
		    else
			varB[0] -= .5*av*av;
		}
	    } /*indN < N */
	} /* normal */
    } /* no gradients */
    else if (getdW & !getdH){
	if (strcmp(distr,"poisson") == 0){
	    for (indN = 0; indN < N; indN++){
		next_indx = jc_x[indN+1]; /* next column index */

		for (indM = 0; indM<M; indM++){
		    av = 0;
		    ind1 = indN;
		    for (ind2=indM*D;ind2<indM*D+D;ind2++){
			av += W[ind1] * H[ind2];
			ind1+=N;
		    }
		    expav = exp(av);
		    /* X(i,j) not empty */
		    if ((indx<next_indx)&&(ir_x[indx]==indM)){
			varB[0] += X[sparse_count] * av - expav;
			xwh = X[sparse_count] - expav;
			sparse_count++;
			indx++;
		    }
		    /* only exp(W*H) */
		    else{
			varB[0] -= expav;
			xwh = - expav;
		    }
		    
		    ind2 = indN;
		    for (ind1=indM*D;ind1<indM*D+D;ind1++){
			dW[ind2] += xwh * H[ind1];
			ind2+=N;
		    }
		}
	    

	    } /* indN < N */
	} /* poisson */
    } /* only dW computed */
    else if (getdH & !getdW){
	if (strcmp(distr,"poisson") == 0){

	    for (indN = 0; indN < N; indN++){
		next_indx = jc_x[indN+1]; /* next column index */

		for (indM = 0; indM<M; indM++){
		    av = 0;
		    ind1 = indN;
		    for (ind2=indM*D;ind2<indM*D+D;ind2++){
			  av += W[ind1] * H[ind2];
			ind1+=N;
		    }
		    expav = exp(av);


		    /* X(i,j) not empty */
		    if ((indx<next_indx)&&(ir_x[indx]==indM)){
			varB[0] += X[sparse_count] * av - expav;
			xwh = X[sparse_count] - expav;
			sparse_count++;
			indx++;
		    }
		    /* only exp(W*H) */
		    else{
			varB[0] -= expav;
			xwh = - expav;
		    }
		    ind2 = indN;
		    for (ind1=indM*D;ind1<D+indM*D;ind1++){
			dH[ind1] += xwh * W[ind2];
			ind2+=N;
		    }
		}
	    

	    } /* indN < N */
	} /* poisson */

    } /* only dH */
    else if (getdH & getdW){
	if (strcmp(distr,"poisson") == 0){

	    for (indN = 0; indN < N; indN++){
		next_indx = jc_x[indN+1]; /* next column index */

		for (indM = 0; indM<M; indM++){
		    av = 0;
		    ind1 = indN;
		    for (ind2=indM*D;ind2<indM*D+D;ind2++){
			  av += W[ind1] * H[ind2];
			ind1+=N;
		    }
		    expav = exp(av);

		    /* X(i,j) not empty */
		    if ((indx<next_indx)&&(ir_x[indx]==indM)){
			varB[0] += X[sparse_count] * av - expav;
			xwh = X[sparse_count] - expav;
			sparse_count++;
			indx++;
		    }
		    /* only exp(W*H) */
		    else{
			varB[0] -= expav;
			xwh = -expav;
		    }
		    ind2 = indN;
		    for (ind1=indM*D;ind1<D+indM*D;ind1++){
			dW[ind2] += xwh * H[ind1];
			dH[ind1] += xwh * W[ind2];
			ind2+=N;
		    }
			}
	    

	    } /* indN < N */
	} /* poisson */

    }

    return;
}


