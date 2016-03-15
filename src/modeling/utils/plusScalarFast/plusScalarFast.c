/**
 * [ids, coefs] = plusScalarFast(id1, id2, coefs1, coefs2)
 *
 * Compile with
 * mex plusScalarFast.c
 */

#include "mex.h"
#include <string.h>

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    /* Check nbr of inputs and outputs */
    if (nrhs != 4 || nlhs != 2)
        mexErrMsgTxt("plusScalarFast needs 4 inputs and 2 outputs") ;
    
    const mxArray* pid1 = prhs[0] ;
    const mxArray* pid2 = prhs[1] ;
    const mxArray* pcoefs1 = prhs[2] ;
    const mxArray* pcoefs2 = prhs[3] ;
    
    if (NULL == pid1 || NULL == pid2 || NULL == pcoefs1 || NULL == pcoefs2)
        mexErrMsgTxt("At least one of the 4 inputs arguments is NULL.") ;
    
    /* Inputs should be column vectors */
    if (mxGetN(pid1) != 1 || mxGetN(pid2) != 1 || mxGetN(pcoefs1) != 1 || mxGetN(pcoefs2) != 1)
        mexErrMsgTxt("inputs arguments should be column vectors") ;
    
    mwSize dims1 = mxGetM(pid1) ;
    mwSize dims2 = mxGetM(pid2) ;
    
    /* They should have the same size */
    if (mxGetM(pcoefs1) != dims1 || mxGetM(pcoefs2) != dims2)
        mexErrMsgTxt("corresponding inputs arguments should have the same length") ;
    
    /* Inputs should be double */
    if ( !mxIsDouble(pcoefs1) || !mxIsDouble(pcoefs2) || !mxIsDouble(pid1) || !mxIsDouble(pid2) )
        mexErrMsgTxt("Inputs should be double") ;
    
    /* Retreive data */    
    double* id1    = mxGetPr(pid1) ;
    double* id2    = mxGetPr(pid2) ;
    double* coefs1 = mxGetPr(pcoefs1) ;
    double* coefs2 = mxGetPr(pcoefs2) ;
    
    mwSize sizeCommon ;
    mwIndex i, j, nCommon, nOnly ;
    int match ;
    
    /* 0/1 array */
    int* id1bis = malloc(dims1 * sizeof(int)) ;
    int* id2bis = malloc(dims2 * sizeof(int)) ;
    if (NULL == id1bis || NULL == id2bis)
        mexErrMsgTxt("Error allocating memory.") ;
    for(i = 0 ; i < dims1 ; i++)
        id1bis[i] = 0 ;
    for(i = 0 ; i < dims2 ; i++)
        id2bis[i] = 0 ;
    mwIndex* whereToLook = malloc(dims2 * sizeof(mwIndex)) ;
    if (NULL == whereToLook)
        mexErrMsgTxt("Error allocating memory.") ;
    
    /* Prepare output */
    double* coefs = malloc((dims1+dims2) * sizeof(double)) ; /* At most dims1 + dims2 elements */
    double* ids   = malloc((dims1+dims2) * sizeof(double)) ;
    if (NULL == coefs || NULL == ids)
        mexErrMsgTxt("Error allocating memory.") ;
    
    /* Find intersection */
    sizeCommon = 0 ;
    for(i = 0 ; i < dims1 ; i++) {
        match = 0 ;
        for(j = 0 ; j < dims2 ; j++) {
            if (id1[i] == id2[j]) { /* Match */
                id1bis[i] = 1 ;
                id2bis[j] = 1 ;
                whereToLook[j] = sizeCommon ;
                sizeCommon ++ ;
                break ;
            }
        }
    }
    
    /* Now, fill vectors */
    nCommon = 0 ; nOnly = 0 ;
    /* The first one */
    for(i = 0 ; i < dims1 ; i++) {
        if (id1bis[i] == 1) { /* A matching value */
            coefs[nCommon] = coefs1[i] ;
            ids[nCommon] = id1[i] ;
            nCommon ++ ;
        } else {
            coefs[sizeCommon + nOnly] = coefs1[i] ;
            ids[sizeCommon + nOnly] = id1[i] ;
            nOnly ++ ;
        }
    }
    /* The second one */
    for(i = 0 ; i < dims2 ; i++) {
        if (id2bis[i] == 1) { /* A matching value */
            coefs[whereToLook[i]] += coefs2[i] ;
        } else {
            coefs[sizeCommon + nOnly] = coefs2[i] ;
            ids[sizeCommon + nOnly] = id2[i] ;
            nOnly ++ ;
        }
    }
    
    plhs[0] = mxCreateDoubleMatrix(nCommon + nOnly, 1, mxREAL) ;
    plhs[1] = mxCreateDoubleMatrix(nCommon + nOnly, 1, mxREAL) ;            
    if (plhs[0] == NULL || plhs[1] == NULL)
        mexErrMsgTxt("Could not create mxArray.\n");
    
    memcpy(mxGetData(plhs[0]), ids,   (nCommon + nOnly) * sizeof(double)) ;
    memcpy(mxGetData(plhs[1]), coefs, (nCommon + nOnly) * sizeof(double)) ;
        
    free(id1bis) ;
    free(id2bis) ;
    free(whereToLook) ;           
}