/**
 * out = mergeCells(in)
 * Merge the n-D double arrays in the n-D cell in
 *
 * Compile with
 * mex mergeCells.c
 */

#include "mex.h"
#include <string.h>

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    /* Check nbr of inputs and outputs */
    if (nrhs != 1 || nlhs != 1)
        mexErrMsgTxt("mergeCells needs 1 output and 1 input.") ;
    const mxArray* cells = prhs[0] ;    
    if (NULL == cells)
        mexErrMsgTxt("First argument is NULL.") ;
    if (! mxIsCell(cells))
        mexErrMsgTxt("The unique argument should be a cell of n-D double arrays") ;            
    
    /* Step 1 : get the total size */
    mwSize nCells = mxGetNumberOfElements(cells) ;
    mxArray* cellContent ;
    mwSize nTotal = 0 ;
    mwIndex i, j, k ;
    for(i = 0 ; i < nCells ; i++) {
        cellContent = mxGetCell(cells, i) ;
        if (NULL == cellContent)
            continue ;        
        /* Get the size of cellContent, which should be a numeric array */
        if (! mxIsDouble(cellContent))
            mexErrMsgTxt("All the elements of input should be n-D double arrays") ;
        nTotal += mxGetNumberOfElements(cellContent) ;
    }
    
    /* Step 2 : prepare output array */
    plhs[0] = mxCreateDoubleMatrix(nTotal, 1, mxREAL);
    if (NULL == plhs[0])
        mexErrMsgTxt("Error allocating array.") ;
    
    /* Step 3 : fill it */
    double* data = mxGetPr(plhs[0]) ;
    double* data2 ;
    k = 0 ;
    for(i = 0 ; i < nCells ; i++) {
        cellContent = mxGetCell(cells, i) ;
        if (NULL == cellContent)
            continue ;        
        data2 = mxGetPr(cellContent) ;
        for(j = 0 ; j < mxGetNumberOfElements(cellContent) ; j++) {
            data[k++] = data2[j] ;
        }
    }
}