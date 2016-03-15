/**
 * idx = forwardMapping(backwardMapping, ids)
 * Returns the indices of the elements of ids in backwardMapping
 *
 * Compile with
 * mex forwardMapping.c
 */

#include "mex.h"
#include <string.h>

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    /* Check nbr of inputs and outputs */
    if (nrhs != 2 || nlhs != 1)
        mexErrMsgTxt("forwardMapping needs 1 outputs and 2 inputs") ;
    
    const mxArray* pBackwardMapping = prhs[0] ;
    const mxArray* pIds = prhs[1] ;
    
    if (NULL == pBackwardMapping || NULL == pIds)
        mexErrMsgTxt("One or both of the input arguments are NULL.") ;
    
    /* Inputs should be column vectors */
    if (mxGetN(pBackwardMapping) != 1 || mxGetN(pIds) != 1)
        mexErrMsgTxt("inputs arguments should be column vectors") ;
    
    mwSize nBack = mxGetM(pBackwardMapping) ;
    mwSize nIds = mxGetM(pIds) ;
    
    double* back = mxGetPr(pBackwardMapping) ;
    double* ids = mxGetPr(pIds) ;
    
    mwIndex* idx = malloc(nIds * sizeof(mwIndex)) ;
    mwIndex i,j,k ;
    
    for(i = 0 ; i < nIds ; i++) {
        k = -1 ;
        for(j = 0 ; j < nBack ; j++) {
            if (back[j] == ids[i]) {
                k = j ;
                break ;
            }
        }
        idx[i] = k+1 ;
    }
    
    /* Return results */
    if (sizeof(mwIndex) == 64) {
        plhs[0] = mxCreateNumericMatrix(nIds, 1, mxUINT64_CLASS, mxREAL);
    } else {
        plhs[0] = mxCreateNumericMatrix(nIds, 1, mxUINT32_CLASS, mxREAL);
    }
    if (plhs[0] == NULL)
        mexErrMsgTxt("Could not create mxArray.\n");
    
    memcpy(mxGetData(plhs[0]), idx, nIds * sizeof(mwIndex)) ;
    free(idx) ;
}