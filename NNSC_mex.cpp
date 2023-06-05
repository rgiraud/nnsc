#include<mex.h>
#include<matrix.h>
#include"./include/NNSC.h"


void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    unsigned char* R,* G,* B;
    int SuperpixelNum=0;
    int nRows=0;int nCols=0;
    int* label;
    int* output;
    unsigned char* img;
    
    //NNSC parameters
    float ratio = 0.075;  //Can be given in input
    int patch_w = 3;
    int Kpm = 4;  
    int iterationNum = 8;
    
    if(nrhs>=3)
    {
        if(mxGetNumberOfDimensions(prhs[0])!=3)
            mexErrMsgTxt("The input image must be in CIERGB form");
        if(mxGetClassID(prhs[0])!=mxUINT8_CLASS)
            mexErrMsgTxt("The input image must be in CIERGB form");
        nRows=mxGetM(prhs[0]);
        nCols=mxGetN(prhs[0])/3;
        
        img = (unsigned char*)mxGetPr(prhs[0]);
        SuperpixelNum = (int) mxGetScalar(prhs[1]);
        ratio = (float) mxGetScalar(prhs[2]);     
        
    }
    
    int pixel=nRows*nCols;
    R=new unsigned char[pixel];
    G=new unsigned char[pixel];
    B=new unsigned char[pixel];
    label=new int[pixel];
    
    for(int i=0;i<pixel;i++)
    {
        R[i]=img[i];
        G[i]=img[i+pixel];
        B[i]=img[i+pixel+pixel];
    }
    
    
    //NNSC
    NNSC(R, G, B, nCols, nRows, SuperpixelNum, ratio, label, Kpm, patch_w, iterationNum);
    
    
    //Output label map
    plhs[0]=mxCreateNumericMatrix(nRows,nCols,mxINT32_CLASS,mxREAL);
    output=(int*)mxGetPr(plhs[0]);
    for(int i=0;i<pixel;i++)
        output[i]=label[i];
    
    delete [] R;
    delete [] G;
    delete [] B;
    delete [] label;
}
