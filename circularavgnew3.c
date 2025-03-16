#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

# include "mex.h"
# include "math.h"
# define MAX_Array_SIZE 5000

int validatedata(double mask, double data, double maxlimit, double minlimit)
{
int retval = 0;
if (mask > 0) if ((data > minlimit) & (data < maxlimit)) retval = 1;
return retval;
}

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{ 
    /* circularavgnew(img, mask, qcMap(solid angle correction map),  qrMap(q-index map), qArray(q values), offset, limits(lowest, highest valid counts));
     * Xiaobing Zuo, 3/2011
    */
    
    
    double *img, *qrMap, *qcMap, *limits, *mask, *qArray, *data;
    double offset;
    int row, col;
    double pi = 3.141592;
    int isMask = 1;

    const int DEF_MAXL = 20000000;
    const int DEF_MINL = 0;
    double maxlimit = 10000000;
    double minlimit = 0;
    double maskdata = 0.0;
    double imgvEff = 0.0;
    double imgvEff2 = 0.0;
    double qCorrVal = 0.0;
    
    int qRMapVal=1;
    int qN, qNum;
    long indimg;
    
    int r2;
    int i, j, m;
    int isvalid, isvalid1;
    double imgv=0.0;
    double npxN = 1.0;
    //double avguc[MAX_Array_SIZE]; /* solid angle uncorrected */
    double avg[MAX_Array_SIZE];
    double npx[MAX_Array_SIZE];
    //double avg2[MAX_Array_SIZE];
    double avg3[MAX_Array_SIZE];
    
    
	if (nrhs != 7) {
        mexPrintf("%i\n", nrhs);
        mexErrMsgTxt("7 inputs: image, qRMap, qCorrMap, mask, qArray, offset, [highL, lowL] are required.");
    }
    

	if (nlhs != 1)
	    mexErrMsgTxt("1 output is required.");
    
	img   = mxGetPr(prhs[0]);
    mask  = mxGetPr(prhs[1]);   
    qcMap = mxGetPr(prhs[2]);
    qrMap = mxGetPr(prhs[3]);
    

	qArray = mxGetPr(prhs[4]);
	offset = mxGetScalar(prhs[5]);
    
	limits = mxGetPr(prhs[6]);
    
	row = (int) mxGetM(prhs[0]);
	col = (int) mxGetN(prhs[0]);	
    qNum = (int) mxGetNumberOfElements(prhs[4]);
    
    if (mxGetNumberOfElements(prhs[6]) >1) {
        if (limits[0]>limits[1]){
            maxlimit = limits[0];    
            minlimit = limits[1];
		}else{
            minlimit = limits[0];    
            maxlimit = limits[1];
        }
	}else {
        maxlimit = DEF_MAXL;
        minlimit = DEF_MINL;
    }
    
    plhs[0] = mxCreateDoubleMatrix(qNum, 6, mxREAL);
    data = mxGetPr(plhs[0]);


    for (i=0;i<qNum;i++){
        //avguc[i]=0.;
        avg[i]=0.;
        //avg2[i]=0.;
        avg3[i]=0.;
        npx[i]=0;
        data[i]=0.0;
        data[i+qNum]=0.0;
        data[i+2*qNum]=0.0;
        data[i+3*qNum]=0.0;
        data[i+4*qNum]=0.0;
    }
    

    
    for (i = 0; i < row; i++) {
        for (j = 0; j < col; j++) {
            indimg = i + j*row;
            imgv = *(img + indimg);
            isvalid = 1;
            if (isMask) {
                maskdata = *(mask + indimg);
            } 
            else 
                maskdata = 1.0;

            isvalid1 = validatedata(maskdata, imgv, maxlimit, minlimit);
            isvalid = isvalid*isvalid1;
            if (isvalid) {
                qCorrVal = *(qcMap + indimg);
                qRMapVal = (int)(*(qrMap + indimg));
                
                imgvEff2 = imgv - offset;
                /* imgvEff = imgvEff2 * qCorrVal * 0.000001; */
                /* imgvEff = imgvEff2 * qCorrVal;  */
                //imgvEff = imgvEff2 ; 
                /* mexPrintf("here3 qRMapVal %f\n.",imgvEff); */
                //avguc[qRMapVal] += imgv;
                avg[qRMapVal] += imgvEff2;
                /* avg2[qRMapVal] += imgvEff * imgvEff; */
                //avg2[qRMapVal] += imgvEff2 * imgvEff2;
                avg3[qRMapVal] += qCorrVal;
                npx[qRMapVal] += 1.0;
             
            }
        }
    }

    for (i=0;i<qNum;i++)
	{
        npxN = npx[i];
    	data[i] = qArray[i];                /* data(1) = q_i  */
        if (npx[i]==0) { 
            data[i+qNum] = 0.0;
            data[i+2*qNum] = 0.0;
            data[i+5*qNum] = 0.0;
        }
        else {
            data[i+qNum] = avg[i]/npxN ;          /* data(2) =average of iq_i  */
            
            data[i+2*qNum] = sqrt(avg[i]) /npxN;  /* data(3) = sqrt(Iq/N) */  
            data[i+5*qNum] = avg3[i]/npxN;  /* data(6) =average of corr_i */
        }
        //data[i+3*qNum] = avg2[i] /npxN;;     /* data(4) = average of iq_i ^2 */
        data[i+4*qNum] = npxN;               /* data(5) = pixel#_q_i */
        /* mexPrintf("%f   %f   %f   %f   %f\n",data[i],data[i+qNum],data[i+2*qNum],data[i+3*qNum],data[i+4*qNum] ); */
	}   
 
}