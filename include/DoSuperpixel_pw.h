#ifndef DOSUPERPIXEL
#define DOSUPERPIXEL

#include<vector>
#include"preEnforceConnectivity.h"
#include<algorithm>
#include"EnforceConnectivity.h"
#include"point.h"
#include<string.h>

#ifndef MAX
#define MAX(a, b) ((a)>(b)?(a):(b))
#endif

#ifndef MIN
#define MIN(a, b) ((a)<(b)?(a):(b))
#endif

#define EPSILON 0.0000001
#define omega 3

#ifndef enc_type
#define enc_type unsigned char
#endif



using namespace std;


//Perform weighted kmeans iteratively in the ten dimensional feature space.

void lookup_table2(float * g_s, float * g_r, int max_hw, float sigma_s, float sigma_c) {
    
    float  s2  = sigma_s*sigma_s;
    for(int i=0; i< max_hw; i++){
        float v = exp(-0.5*(i*i)/(s2));
        if(v<0.1){
            //hs = i-1;
            break;
        }
        g_s[i]=v;
    }
    //pre-compute range gaussian
    
    float r2 = sigma_c*sigma_c;
    for(int i=0; i<256; i++){
        g_r[i] = exp(-0.5*(i*i)/r2);
    }
    
}


// Measure distance between entire 2 patches,
// terminating early if we exceed a cutoff distance
float dist_2D_p_p(float* L1_1, float* L2_1, float* a1_1, float* a2_1, float* b1_1, float* b2_1,
        int ax, int ay, int bx, int by,
        int ha, int wa, int hb, int pw) {
    
    float ans = 0;
    int dx, dy;
    float ld = 0, ad = 0, bd = 0, ld2 = 0, ad2 = 0, bd2 = 0;
    
    float count = 0;
    for (dx = -pw; dx <= pw; dx++) {
        for (dy = -pw; dy <= pw; dy++) {
            
            if ( ((ax+dx)>=0) && ((ax+dx)<wa) && ((ay+dy)>=0) && ((ay+dy)<ha) ) {
                
                ld  = L1_1[(ax+dx)*ha+ay+dy] - L1_1[(bx+dx)*hb+by+dy];
                ad  = a1_1[(ax+dx)*ha+ay+dy] - a1_1[(bx+dx)*hb+by+dy];
                bd  = b1_1[(ax+dx)*ha+ay+dy] - b1_1[(bx+dx)*hb+by+dy];
                ld2 = L2_1[(ax+dx)*ha+ay+dy] - L2_1[(bx+dx)*hb+by+dy];
                ad2 = a2_1[(ax+dx)*ha+ay+dy] - a2_1[(bx+dx)*hb+by+dy];
                bd2 = b2_1[(ax+dx)*ha+ay+dy] - b2_1[(bx+dx)*hb+by+dy];
                ans += ld*ld + ld2*ld2    + ad*ad + bd*bd + ad2*ad2 + bd2*bd2;
                
                count += 1;
            }
        }
    }
    ans /= count;
    
    return ans;
}



/* PROPAGATION STEP: Measure distance between 2 patches.
 * Only computes the two limit lines that result from the shifting.
 */
// Measure distance between entire 2 patches,
// terminating early if we exceed a cutoff distance

float dist_2D_p_p_prop(float* L1_1, float* L2_1, float* a1_1, float* a2_1, float* b1_1, float* b2_1,
        int ax, int ay, int bx, int by, int ha, int hb, int patch_w, int x_offset, int y_offset) {
    
    float band_add_d=0, band_cut_d=0, ans = 0;
    int off=0;
    int dx, dy;
    float band_add_L1 = 0, band_add_L2 = 0, band_add_a1 = 0, band_add_a2 = 0, band_add_b1 = 0, band_add_b2 = 0;
    float band_cut_L1 = 0, band_cut_L2 = 0, band_cut_a1 = 0, band_cut_a2 = 0, band_cut_b1 = 0, band_cut_b2 = 0;
    
    if (x_offset !=0){
        
        if (x_offset == 1)
            off = 1;
        else
            off = 0;
        
        
        //Bands to withdraw and to add
        //Computation on only one dimension
        for (dy = -patch_w; dy <= patch_w; dy++) {
            
            band_cut_L1 = L1_1[(ax-patch_w-off)*ha+ay+dy] -  L1_1[(bx-patch_w-off)*hb+by+dy];
            band_cut_L2 = L2_1[(ax-patch_w-off)*ha+ay+dy] -  L2_1[(bx-patch_w-off)*hb+by+dy];
            band_cut_a1 = a1_1[(ax-patch_w-off)*ha+ay+dy] -  a1_1[(bx-patch_w-off)*hb+by+dy];
            band_cut_a2 = a2_1[(ax-patch_w-off)*ha+ay+dy] -  a2_1[(bx-patch_w-off)*hb+by+dy];
            band_cut_b1 = b1_1[(ax-patch_w-off)*ha+ay+dy] -  b1_1[(bx-patch_w-off)*hb+by+dy];
            band_cut_b2 = b2_1[(ax-patch_w-off)*ha+ay+dy] -  b2_1[(bx-patch_w-off)*hb+by+dy];
            
            band_cut_d += band_cut_L1*band_cut_L1 + band_cut_L2*band_cut_L2 +
                    band_cut_a1*band_cut_a1 + band_cut_a2*band_cut_a2 +
                    band_cut_b1*band_cut_b1 + band_cut_b2*band_cut_b2;
            
            
            band_add_L1 = L1_1[(ax+patch_w+1-off)*ha+ay+dy] -  L1_1[(bx+patch_w+1-off)*hb+by+dy];
            band_add_L2 = L2_1[(ax+patch_w+1-off)*ha+ay+dy] -  L2_1[(bx+patch_w+1-off)*hb+by+dy];
            band_add_a1 = a1_1[(ax+patch_w+1-off)*ha+ay+dy] -  a1_1[(bx+patch_w+1-off)*hb+by+dy];
            band_add_a2 = a2_1[(ax+patch_w+1-off)*ha+ay+dy] -  a2_1[(bx+patch_w+1-off)*hb+by+dy];
            band_add_b1 = b1_1[(ax+patch_w+1-off)*ha+ay+dy] -  b1_1[(bx+patch_w+1-off)*hb+by+dy];
            band_add_b2 = b2_1[(ax+patch_w+1-off)*ha+ay+dy] -  b2_1[(bx+patch_w+1-off)*hb+by+dy];
            
            band_add_d += band_add_L1*band_add_L1 + band_add_L2*band_add_L2 +
                    band_add_a1*band_add_a1 + band_add_a2*band_add_a2 +
                    band_add_b1*band_add_b1 + band_add_b2*band_add_b2;
            
        }
        //The bands to withdraw and to add are not the same according to the offset direction
        ans = x_offset*(band_add_d-band_cut_d);
    }
    
    else  {
        
        if (y_offset == 1)
            off = 1;
        else
            off = 0;
        
        //Bands to withdraw and to add
        //Computation on only one dimension
        for (dx = -patch_w; dx <= patch_w; dx++) {
            
            band_cut_L1 = L1_1[(ax+dx)*ha+ay-patch_w-off] -  L1_1[(bx+dx)*hb+by-patch_w-off];
            band_cut_L2 = L2_1[(ax+dx)*ha+ay-patch_w-off] -  L2_1[(bx+dx)*hb+by-patch_w-off];
            band_cut_a1 = a1_1[(ax+dx)*ha+ay-patch_w-off] -  a1_1[(bx+dx)*hb+by-patch_w-off];
            band_cut_a2 = a2_1[(ax+dx)*ha+ay-patch_w-off] -  a2_1[(bx+dx)*hb+by-patch_w-off];
            band_cut_b1 = b1_1[(ax+dx)*ha+ay-patch_w-off] -  b1_1[(bx+dx)*hb+by-patch_w-off];
            band_cut_b2 = b2_1[(ax+dx)*ha+ay-patch_w-off] -  b2_1[(bx+dx)*hb+by-patch_w-off];
            
            band_cut_d += band_cut_L1*band_cut_L1 + band_cut_L2*band_cut_L2 +
                    band_cut_a1*band_cut_a1 + band_cut_a2*band_cut_a2 +
                    band_cut_b1*band_cut_b1 + band_cut_b2*band_cut_b2;
            
            
            band_add_L1 = L1_1[(ax+dx)*ha+ay+patch_w+1-off] -  L1_1[(bx+dx)*hb+by+patch_w-off+1];
            band_add_L2 = L2_1[(ax+dx)*ha+ay+patch_w+1-off] -  L2_1[(bx+dx)*hb+by+patch_w-off+1];
            band_add_a1 = a1_1[(ax+dx)*ha+ay+patch_w+1-off] -  a1_1[(bx+dx)*hb+by+patch_w-off+1];
            band_add_a2 = a2_1[(ax+dx)*ha+ay+patch_w+1-off] -  a2_1[(bx+dx)*hb+by+patch_w-off+1];
            band_add_b1 = b1_1[(ax+dx)*ha+ay+patch_w+1-off] -  b1_1[(bx+dx)*hb+by+patch_w-off+1];
            band_add_b2 = b2_1[(ax+dx)*ha+ay+patch_w+1-off] -  b2_1[(bx+dx)*hb+by+patch_w-off+1];
            
            band_add_d += band_add_L1*band_add_L1 + band_add_L2*band_add_L2 +
                    band_add_a1*band_add_a1 + band_add_a2*band_add_a2 +
                    band_add_b1*band_add_b1 + band_add_b2*band_add_b2;
            
            
        }
        //The bands to withdraw and to add are not the same according to the offset direction
        ans = y_offset*(band_add_d-band_cut_d);
        
    }
    return ans;
}



// Measure distance between entire 2 patches,
// terminating early if we exceed a cutoff distance
float dist_2D_p_S(float* L1_1, float* L2_1, float* a1_1, float* a2_1, float* b1_1, float* b2_1, float Lab_2, int ax, int ay, int ha,
        float centerL1, float centerL2, float centera1, float centera2, float centerb1, float centerb2, float Lab_cc) {
    
    float ans = 0;
    int pos = ay+ax*ha;
    
    ans = Lab_cc + Lab_2 - 2*(centerL1*L1_1[pos] + centerL2*L2_1[pos] + centera1*a1_1[pos] + centera2*a2_1[pos] + centerb1*b1_1[pos] + centerb2*b2_1[pos]);
    
    return ans;
}



// Compare two distance and change the optimal ptr if the new distance is
// lower than the previous best
float improve_guess_2D_prop(int ax, int ay, int bx, int by, int ha, int wa, int hb, int wb, int pw, int *xbest_ptr, int *ybest_ptr, float *distbest_ptr,
        float* L1_1, float* L2_1, float* a1_1, float* a2_1, float* b1_1, float* b2_1, float* x1_1, float* x2_1, float* y1_1, float* y2_1, float Lab_2, float xy_2,
        float centerL1, float centerL2, float centera1, float centera2, float centerb1, float centerb2, float centerx1, float centerx2, float centery1, float centery2, float Lab_cc, float xy_cc,
        float cx, float cy, float step, float beta, float sp_var, int type,
        float * nnfd_pw, int x_offset, int y_offset)  {
    
    int pos = ay+ax*ha;
    float d = 0, dist_pm = 0;
    
    if ( (type > 1) && ((ax-pw-1)>=0) && ((ax+pw+1)<wa) && ((ay-pw-1)>=0) && ((ay+pw+1)<ha) &&
            ((bx-pw-1)>=0) && ((bx+pw+1)<wb) && ((by-pw-1)>=0) && ((by+pw+1)<hb) ) {
        
        //optim dist prop
        d = nnfd_pw[(ax-x_offset)*ha + ay-y_offset];
        dist_pm = d + dist_2D_p_p_prop(L1_1, L2_1, a1_1, a2_1, b1_1, b2_1, ax, ay, bx, by, ha, hb, pw, x_offset, y_offset)/((2*pw+1)*(2*pw+1));
    }
    else {
        dist_pm = dist_2D_p_p(L1_1, L2_1, a1_1, a2_1, b1_1, b2_1, ax, ay, bx, by, ha, wa, hb, pw);
    }
    
    float dist_c  = dist_2D_p_S(L1_1, L2_1, a1_1, a2_1, b1_1, b2_1, Lab_2, ax, ay, ha, centerL1, centerL2, centera1, centera2, centerb1, centerb2, Lab_cc);
    float dist_s  = xy_cc + xy_2 - 2*(centerx1*x1_1[pos] + centerx2*x2_1[pos] + centery1*y1_1[pos] + centery2*y2_1[pos]);
    
    float sigma = step*step;
    float dist_pm_s = (1 - exp(- ((cx-bx)*(cx-bx) + (cy-by)*(cy-by))/(sigma)));
    dist_pm_s *= 2*beta;
    
    
    float dist = dist_c*0.5 + dist_s*sp_var*sp_var + dist_pm + dist_pm_s*sp_var*sp_var;
    
    //Comparison
    if (dist < *distbest_ptr) {
        *distbest_ptr = dist;
        *xbest_ptr = bx;
        *ybest_ptr = by;
        nnfd_pw[pos] = dist_pm;
        
    }
    
    return dist;
    
    
}



// Compare two distance and change the optimal ptr if the new distance is
// lower than the previous best
float improve_guess_2D(int ax, int ay, int bx, int by, int ha, int wa, int hb, int pw, int *xbest_ptr, int *ybest_ptr, float *distbest_ptr,
        float* L1_1, float* L2_1, float* a1_1, float* a2_1, float* b1_1, float* b2_1, float* x1_1, float* x2_1, float* y1_1, float* y2_1, float Lab_2, float xy_2,
        float centerL1, float centerL2, float centera1, float centera2, float centerb1, float centerb2, float centerx1, float centerx2, float centery1, float centery2, float Lab_cc, float xy_cc,
        float cx, float cy, float step, float beta, float sp_var)  {
    
    int pos = ay+ax*ha;
    
    //Distance computation
    float dist_pm = dist_2D_p_p(L1_1, L2_1, a1_1, a2_1, b1_1, b2_1, ax, ay, bx, by, ha, wa, hb, pw);
    float dist_c  = dist_2D_p_S(L1_1, L2_1, a1_1, a2_1, b1_1, b2_1, Lab_2, ax, ay, ha, centerL1, centerL2, centera1, centera2, centerb1, centerb2, Lab_cc);
    float dist_s  = xy_cc + xy_2 - 2*(centerx1*x1_1[pos] + centerx2*x2_1[pos] + centery1*y1_1[pos] + centery2*y2_1[pos]);    //((cy-ay)*(cy-ay)+(cx-ax)*(cx-ax));
    
    float sigma = step*step;
    float dist_pm_s = (1 - exp(- ((cx-bx)*(cx-bx) + (cy-by)*(cy-by))/(sigma)));
    dist_pm_s *= 2*beta;
    
    
    float dist = dist_c*0.5 + dist_s*sp_var*sp_var + dist_pm + dist_pm_s*sp_var*sp_var;
    
    
    //Comparison
    if (dist < *distbest_ptr) {
        *distbest_ptr = dist;
        *xbest_ptr = bx;
        *ybest_ptr = by;
        
    }
    
    return dist;
    
    
}




int
        pm_rand(unsigned long* next)
{
    (*next) = (*next) * 1103515245 + 12345;
    return (unsigned int)((*next)/65536) % 32768;
}




typedef struct{
    
    unsigned long next;
    
    float** L1;
    float** L2;
    float** a1;
    float** a2;
    float** b1;
    float** b2;
    float** x1;
    float** x2;
    float** y1;
    float** y2;
    float *L1_1;
    float* L2_1;
    float* a1_1;
    float* a2_1;
    float* b1_1;
    float* b2_1;
    float* x1_1;
    float* x2_1;
    float* y1_1;
    float* y2_1;
    float* xy_2;
    float* Lab_2;
    float** W;
    
    unsigned char* R;
    unsigned char* G;
    unsigned char* B;
    
    float* centerL1_in;
    float* centerL2_in;
    float* centera1_in;
    float* centera2_in;
    float* centerb1_in;
    float* centerb2_in;
    float* centerx1_in;
    float* centerx2_in;
    float* centery1_in;
    float* centery2_in;
    point * seedArray_in;
    float* Lab_cc_in;
    float* xy_cc_in;
    float* ss_in;
    int* label_in;
    
    int nRows;
    int nCols;
    int seedNum;
    int Kpm_nbr;
    int * label_out;
    float * dist_out;
    int pw;
    int StepX;
    int StepY;
    int iterationNum;
    
}pm_struct;


void*
        pm_core(void *arg)
{
    pm_struct inputs = *(pm_struct*) arg;
    
    float** L1 = inputs.L1;
    float** L2 = inputs.L2;
    float** a1 = inputs.a1;
    float** a2 = inputs.a2;
    float** b1 = inputs.b1;
    float** b2 = inputs.b2;
    float** x1 = inputs.x1;
    float** x2 = inputs.x2;
    float** y1 = inputs.y1;
    float** y2 = inputs.y2;
    
    float *L1_1 = inputs.L1_1;
    float* L2_1 = inputs.L2_1;
    float* a1_1 = inputs.a1_1;
    float* a2_1 = inputs.a2_1;
    float* b1_1 = inputs.b1_1;
    float* b2_1 = inputs.b2_1;
    float* Lab_2 = inputs.Lab_2;
    
    unsigned char* R = inputs.R;
    unsigned char* G = inputs.G;
    unsigned char* B = inputs.B;
    
    float* x1_1 = inputs.x1_1;
    float* x2_1 = inputs.x2_1;
    float* y1_1 = inputs.y1_1;
    float* y2_1 = inputs.y2_1;
    float* xy_2 = inputs.xy_2;
    float** W = inputs.W;
    
    float* centerL1_in = inputs.centerL1_in;
    float* centerL2_in = inputs.centerL2_in;
    float* centera1_in = inputs.centera1_in;
    float* centera2_in = inputs.centera2_in;
    float* centerb1_in = inputs.centerb1_in;
    float* centerb2_in = inputs.centerb2_in;
    float* centerx1_in = inputs.centerx1_in;
    float* centerx2_in = inputs.centerx2_in;
    float* centery1_in = inputs.centery1_in;
    float* centery2_in = inputs.centery2_in;
    
    float* Lab_cc_in = inputs.Lab_cc_in;
    float* xy_cc_in = inputs.xy_cc_in;
    point * seedArray_in = inputs.seedArray_in;
    float* ss_in = inputs.ss_in;
    int* label_in = inputs.label_in;
    
    int StepX = inputs.StepX;
    int StepY = inputs.StepY;
    int nRows = inputs.nRows;
    int nCols = inputs.nCols;
    int seedNum = inputs.seedNum;
    unsigned long next = inputs.next;
    int Kpm_nbr = inputs.Kpm_nbr;
    int * label_out = inputs.label_out;
    float * dist_out = inputs.dist_out;
    int pw = inputs.pw;
    int iterationNum = inputs.iterationNum;
    
    float* centerL1=new float[seedNum];
    float* centerL2=new float[seedNum];
    float* centera1=new float[seedNum];
    float* centera2=new float[seedNum];
    float* centerb1=new float[seedNum];
    float* centerb2=new float[seedNum];
    float* centerx1=new float[seedNum];
    float* centerx2=new float[seedNum];
    float* centery1=new float[seedNum];
    float* centery2=new float[seedNum];
    float* WSum=new float[seedNum];
    int* clusterSize=new int[seedNum];
    float* Lab_cc = new float[seedNum];
    float* xy_cc  = new float[seedNum];
    
    float *ss = (float *) calloc(seedNum,sizeof(float));
    int * label = (int *) calloc(nCols*nRows,sizeof(int));
    point *seedArray = new point[seedNum];
    
    for (int i=0; i<seedNum; i++) {
        seedArray[i].x = seedArray_in[i].x;
        seedArray[i].y = seedArray_in[i].y;
        centerL1[i] = centerL1_in[i];
        centerL2[i] = centerL2_in[i];
        centera1[i] = centera1_in[i];
        centera2[i] = centera2_in[i];
        centerb1[i] = centerb1_in[i];
        centerb2[i] = centerb2_in[i];
        centerx1[i] = centerx1_in[i];
        centerx2[i] = centerx2_in[i];
        centery1[i] = centery1_in[i];
        centery2[i] = centery2_in[i];
        Lab_cc[i] = Lab_cc_in[i];
        xy_cc[i] = xy_cc_in[i];
        ss[i] = ss_in[i];
    }
    
    for (int i=0; i<nRows*nCols; i++)
        label[i] = label_in[i];
    
    
    int   * nnf  = (int *)   calloc(nCols*nRows*2,sizeof(int));
    float * nnfd = (float *) calloc(nCols*nRows,sizeof(float));
    float * nnfd_pw = (float *) calloc(nCols*nRows,sizeof(float));
    
    
    
    int * label_t1 = (int *) malloc(nCols*nRows*sizeof(int));
    for (int i=0; i<nRows*nCols; i++)
        label_t1[i] = label[i];
    
    float * sp_var = (float *)calloc(seedNum, sizeof(float));
    for (int i=0; i<seedNum; i++)
        sp_var[i] = 1;
    float * mean_sp = (float *)calloc(seedNum*3, sizeof(float));
    
    int ha = nRows;
    int wa = nCols;
    int hb = ha;
    int wb = wa;
    
    // Initialize with random nearest neighbor field (NNF)
    // Effective width and height (possible upper left corners of patches)
    int aew = wa, aeh = ha;
    int bew = wb-pw, beh = hb-pw;
    
    int size_a = ha*wa;
    
    int pos, pos_shift;
    int ax, ay, bx, by, ax_off;
    int iter = 0;
    int ystart, yend, ychange, xstart, xend, xchange;
    int xbest, ybest;
    float distbest;
    float *distbest_ptr = &distbest;
    int *xbest_ptr = &xbest;
    int *ybest_ptr = &ybest;
    int rand_val;
    int xp = 0, yp = 0, xmin, xmax, ymin, ymax;
    
    int mag_init = (StepX+StepY)/2;
    
    int thread_nbr = 0;
    int off_nnf = thread_nbr*size_a;
    
    
    // RANDOM INITIALIZATION
    for (ax = 0; ax < aew; ax++) {
        
        ax_off = ha*ax;
        
        for (ay = 0; ay < aeh; ay++) {
            
            pos = ax_off+ay;
            
            xmin = MAX(ax-mag_init, pw);
            xmax = MIN(ax+mag_init+1, bew-1);
            
            ymin = MAX(ay-mag_init, pw);
            ymax = MIN(ay+mag_init+1, beh-1);
            
            
            int loop = 1;
            
            while (loop) {
                
                //random match
                rand_val = pm_rand(&next);
                bx = xmin+rand_val%(xmax-xmin);
                rand_val = pm_rand(&next);
                by = ymin+rand_val%(ymax-ymin);
                
                if ((abs(ax-bx)>omega) || (abs(ay-by)>omega) || (loop > 1)) {
                    loop = 0;
                    
                    //Map init
                    nnf[pos + off_nnf*2] = bx;
                    nnf[pos + off_nnf*2 + size_a] = by;
                    
                    int lab = label[by + bx*hb];
                    *distbest_ptr = FLT_MAX;
                    float dist_k = improve_guess_2D_prop(ax, ay, bx, by, ha, wa, hb, wb, pw, xbest_ptr, ybest_ptr, distbest_ptr,
                            L1_1, L2_1, a1_1, a2_1, b1_1, b2_1, x1_1, x2_1, y1_1, y2_1,  Lab_2[pos], xy_2[pos],
                            centerL1[lab], centerL2[lab], centera1[lab], centera2[lab], centerb1[lab], centerb2[lab], centerx1[lab], centerx2[lab], centery1[lab], centery2[lab], Lab_cc[lab], xy_cc[lab],
                            seedArray[lab].x, seedArray[lab].y,  mag_init, ss[lab], sp_var[lab], 1, nnfd_pw, 0, 0);
                    nnfd[pos + off_nnf] = dist_k;
                    
                    
                }
                else {
                    loop += 1;
                }
                
            }
            
            
        }
    }
    
    while (iter < iterationNum) {
        
        // In each iteration, improve the NNF, by looping in scanline or
        // reverse-scanline order
        ystart = 1; yend = aeh; ychange = 1;
        xstart = 1; xend = aew; xchange = 1;
        
        if (iter % 2 == 1) {
            xstart = xend-1; xend = -1; xchange = -1;
            ystart = yend-1; yend = -1; ychange = -1;
        }
        
        for (ay = ystart; ay != yend; ay += ychange) {
            
            for (ax = xstart; ax != xend; ax += xchange) {
                
                //Current position in a
                pos = ha*ax+ay;
                
                
                int minx = MAX(ax-mag_init,pw);
                int miny = MAX(ay-mag_init,pw);
                int maxx = MIN(ax+mag_init,bew-1);
                int maxy = MIN(ay+mag_init,beh-1);
                
                
                // Current (best) guess
                xbest = nnf[pos + off_nnf*2];
                ybest = nnf[pos + off_nnf*2 + size_a];
                distbest = nnfd[pos + off_nnf];
                distbest_ptr = &distbest;
                xbest_ptr = &xbest;
                ybest_ptr = &ybest;
                
                // PROPAGATION: Improve current guess by trying instead
                // correspondences from left and above (below and right on
                // odd iterations)
                
                //XSHIFT
                if ((unsigned) (ax - xchange-pw) < (unsigned) aew-pw) {
                    pos_shift = ha*(ax-xchange)+ay;
                    xp = nnf[pos_shift + off_nnf*2]  + xchange;
                    yp = nnf[pos_shift + off_nnf*2 + size_a];
                    if ((xp !=  *xbest_ptr)  || (yp !=  *ybest_ptr)) {
                        if ((xp-pw > minx) && (xp+pw<maxx) && (yp-pw > miny) && (yp+pw<maxy) ) {
                            if ( (abs(yp-ay)>omega) || (abs(xp-ax)> omega) ) {
                                
                                int lab = label[yp+xp*hb];
                                
                                improve_guess_2D_prop(ax, ay, xp, yp, ha, wa, hb, wb, pw, xbest_ptr, ybest_ptr, distbest_ptr,
                                        L1_1, L2_1, a1_1, a2_1, b1_1, b2_1, x1_1, x2_1, y1_1, y2_1, Lab_2[pos], xy_2[pos],
                                        centerL1[lab], centerL2[lab], centera1[lab], centera2[lab], centerb1[lab], centerb2[lab], centerx1[lab], centerx2[lab], centery1[lab], centery2[lab],Lab_cc[lab], xy_cc[lab],
                                        seedArray[lab].x, seedArray[lab].y, mag_init, ss[lab], sp_var[lab],  2,
                                        nnfd_pw, xchange, 0);
                                
                            }
                            
                            
                            
                            
                            
                        }
                    }
                    
                }
                
                
                
                //YSHIFT
                if ((unsigned) (ay - ychange-pw) < (unsigned) aeh-pw) {
                    pos_shift = ha*ax+ay -ychange;
                    xp = nnf[pos_shift + off_nnf*2];
                    yp = nnf[pos_shift + off_nnf*2 + size_a] + ychange;
                    if ((xp !=  *xbest_ptr)  || (yp !=  *ybest_ptr)) {
                        if ((xp-pw > minx) && (xp+pw<maxx) && (yp-pw > miny) && (yp+pw<maxy) ) {
                            if ( (abs(yp-ay)>omega) || (abs(xp-ax)> omega) ){
                                
                                int lab = label[yp+xp*hb];
                                
                                improve_guess_2D_prop(ax, ay, xp, yp, ha, wa, hb, wb, pw, xbest_ptr, ybest_ptr, distbest_ptr,
                                        L1_1, L2_1, a1_1, a2_1, b1_1, b2_1, x1_1, x2_1, y1_1, y2_1,  Lab_2[pos],xy_2[pos],
                                        centerL1[lab], centerL2[lab], centera1[lab], centera2[lab], centerb1[lab], centerb2[lab], centerx1[lab], centerx2[lab], centery1[lab], centery2[lab],Lab_cc[lab], xy_cc[lab],
                                        seedArray[lab].x, seedArray[lab].y, mag_init, ss[lab], sp_var[lab], 3,
                                        nnfd_pw, 0, ychange);
                                
                            }
                        }
                    }
                }
                
                // RANDOM SEARCH: Improve current guess by searching in
                // boxes of exponentially decreasing size around the
                // current position
                
                // Sampling window
                for (int mag = mag_init; mag >= 1; mag /= 2) {
                    
                    xbest = *xbest_ptr;
                    ybest = *ybest_ptr;
                    
                    xmin = MAX(xbest-mag, minx);
                    xmax = MIN(xbest+mag, maxx);
                    if(xmin == xmax) continue;
                    
                    ymin = MAX(ybest-mag, miny);
                    ymax = MIN(ybest+mag, maxy);
                    if(ymin == ymax) continue;
                    
                    //Random match
                    xp = (int) xmin+pm_rand(&next)%MAX(xmax-xmin, 1);
                    yp = (int) ymin+pm_rand(&next)%MAX(ymax-ymin, 1);
                    
                    if ((xp-pw > minx) && (xp+pw<maxx) && (yp-pw > miny) && (yp+pw<maxy) ) {
                        if ( (abs(yp-ay)>omega) || (abs(xp-ax)> omega) ){
                            
                            int lab = label[yp+xp*hb];
                            
                            improve_guess_2D_prop(ax, ay, xp, yp, ha, wa, hb, wb, pw, xbest_ptr, ybest_ptr, distbest_ptr,
                                    L1_1, L2_1, a1_1, a2_1, b1_1, b2_1, x1_1, x2_1, y1_1, y2_1, Lab_2[pos], xy_2[pos],
                                    centerL1[lab], centerL2[lab], centera1[lab], centera2[lab], centerb1[lab], centerb2[lab], centerx1[lab], centerx2[lab], centery1[lab], centery2[lab],Lab_cc[lab], xy_cc[lab],
                                    seedArray[lab].x, seedArray[lab].y, mag_init, ss[lab], sp_var[lab],  1,
                                    nnfd_pw, 0, 0);
                            
                        }
                    }
                }
                
                //Map updates
                nnf[pos  + off_nnf*2]           = *xbest_ptr;
                nnf[pos  + off_nnf*2 + size_a]  = *ybest_ptr;
                nnfd[pos + off_nnf]             = *distbest_ptr;
                
                label_t1[pos] = label[*ybest_ptr + *xbest_ptr*hb];
                
                
            }
            
            
            
        }
        
        //Recopy label map
        for (int i=0; i<nRows*nCols; i++)
            label[i] = label_t1[i];
        
        //Update clusters
        if (iter<iterationNum-1) {
            
            for(int i=0;i<seedNum;i++)            {
                centerL1[i]=0;
                centerL2[i]=0;
                centera1[i]=0;
                centera2[i]=0;
                centerb1[i]=0;
                centerb2[i]=0;
                centerx1[i]=0;
                centerx2[i]=0;
                centery1[i]=0;
                centery2[i]=0;
                WSum[i]=0;
                clusterSize[i]=0;
                seedArray[i].x=0;
                seedArray[i].y=0;
                
                sp_var[i] = 0;
                mean_sp[i] = 0;
                mean_sp[i+seedNum] = 0;
                mean_sp[i+seedNum*2] = 0;
                
            }
            
            
            
            for(int i=0;i<nCols;i++)            {
                for(int j=0;j<nRows;j++)                {
                    int L=label[i*nRows+j];
                    float Weight=W[i][j];
                    
                    centerL1[L]+=Weight*L1[i][j];
                    centerL2[L]+=Weight*L2[i][j];
                    centera1[L]+=Weight*a1[i][j];
                    centera2[L]+=Weight*a2[i][j];
                    centerb1[L]+=Weight*b1[i][j];
                    centerb2[L]+=Weight*b2[i][j];
                    centerx1[L]+=Weight*x1[i][j];
                    centerx2[L]+=Weight*x2[i][j];
                    centery1[L]+=Weight*y1[i][j];
                    centery2[L]+=Weight*y2[i][j];
                    
                    clusterSize[L]++;
                    WSum[L]+=Weight;
                    seedArray[L].x+=i;
                    seedArray[L].y+=j;
                    
                    //compute variance of each SP
                    mean_sp[L] += R[i*nRows+j];
                    mean_sp[L+seedNum] += G[i*nRows+j];
                    mean_sp[L+2*seedNum] += B[i*nRows+j];
                    
                    sp_var[L] += R[i*nRows+j]*R[i*nRows+j] + G[i*nRows+j]*G[i*nRows+j] + B[i*nRows+j]*B[i*nRows+j];
                    
                    
                }
            }
            for(int i=0;i<seedNum;i++)            {
                WSum[i]=(WSum[i]==0)?1:WSum[i];
                clusterSize[i]=(clusterSize[i]==0)?1:clusterSize[i];
            }
            for(int i=0;i<seedNum;i++)            {
                centerL1[i]/=WSum[i];
                centerL2[i]/=WSum[i];
                centera1[i]/=WSum[i];
                centera2[i]/=WSum[i];
                centerb1[i]/=WSum[i];
                centerb2[i]/=WSum[i];
                centerx1[i]/=WSum[i];
                centerx2[i]/=WSum[i];
                centery1[i]/=WSum[i];
                centery2[i]/=WSum[i];
                seedArray[i].x/=clusterSize[i];
                seedArray[i].y/=clusterSize[i];
                
                mean_sp[i] /= clusterSize[i];
                mean_sp[i+seedNum] /= clusterSize[i];
                mean_sp[i+2*seedNum] /= clusterSize[i];
                if (sp_var[i] == 0)
                    sp_var[i] = FLT_MAX;
                else
                    sp_var[i] /= clusterSize[i];
                
                
                Lab_cc[i] = centerL1[i]*centerL1[i] + centerL2[i]*centerL2[i] + centera1[i]*centera1[i] +
                        centera2[i]*centera2[i] + centerb1[i]*centerb1[i] + centerb2[i]*centerb2[i];
                
                xy_cc[i]  = centerx1[i]*centerx1[i] + centerx2[i]*centerx2[i] + centery1[i]*centery1[i] + centery2[i]*centery2[i];
                
            }
            
            
            for( int k = 0; k < seedNum; k++ ) {
                sp_var[k] = sp_var[k] - mean_sp[k]*mean_sp[k] - mean_sp[k+seedNum]*mean_sp[k+seedNum] - mean_sp[k+seedNum*2]*mean_sp[k+seedNum*2];
                sp_var[k] = sqrt(sp_var[k])/3;
            }
                
                
                float threshold_var = 25.0; //beta
                for( int k = 0; k < seedNum; k++ ) {
                    sp_var[k] = (float) exp(sp_var[k]/threshold_var);  // m ^ if var sp ^
                    
                
            }
            
            //Init lab map
            for(int i=0;i<seedNum;i++) {
                int x=seedArray[i].x;
                int y=seedArray[i].y;
                
                int minX=(x-(StepX)<=0)?0:x-StepX;
                int minY=(y-(StepY)<=0)?0:y-StepY;
                ss[i] =  (centerx1[i] - x1_1[minX*nRows+minY])*(centerx1[i] - x1_1[minX*nRows+minY]) +
                        (centery1[i] - y1_1[minX*nRows+minY])*(centery1[i] - y1_1[minX*nRows+minY]) +
                        (centerx2[i] - x2_1[minX*nRows+minY])*(centerx2[i] - x2_1[minX*nRows+minY]) +
                        (centery2[i] - y2_1[minX*nRows+minY])*(centery2[i] - y2_1[minX*nRows+minY]);
            }
            
            
            //Recompute distance
            for (int i=0; i<ha*wa; i++)
                nnfd[i] = FLT_MAX;
            
        }
        
        iter += 1;
        
    }
    
    
    //recopy label out
    for (int i=0; i<nRows*nCols; i++) {
        label_out[i+nRows*nCols*Kpm_nbr] = label[i];
        dist_out[i+nRows*nCols*Kpm_nbr] = nnfd[i];
    }
    
    free(label);
    free(label_t1);
    free(nnf);
    free(nnfd);
    free(nnfd_pw);
    
    free(sp_var);
    free(ss);
    free(mean_sp);
    
    //Clear Memory
    delete []centerL1;
    delete []centerL2;
    delete []centera1;
    delete []centera2;
    delete []centerb1;
    delete []centerb2;
    delete []centerx1;
    delete []centerx2;
    delete []centery1;
    delete []centery2;
    delete []Lab_cc;
    delete []xy_cc;
    delete []seedArray;
    delete []clusterSize;
    delete []WSum;
    
    
    pthread_exit(0);
    
    
}





void DoSuperpixel(
        float** L1, float** L2, float** a1, float** a2, float** b1, float** b2,
        float** x1, float** x2, float** y1, float** y2, float** W, int* label,
        point* seedArray, int seedNum,
        int nCols, int nRows, int StepX, int StepY,
        int iterationNum,
        int thresholdCoef,
        unsigned char* R,
        unsigned char* G,
        unsigned char* B,
        int Kpm,
        int patch_w)
{
    
    
    float* L1_1=new float[nCols*nRows]();
    float* L2_1=new float[nCols*nRows]();
    float* a1_1=new float[nCols*nRows]();
    float* a2_1=new float[nCols*nRows]();
    float* b1_1=new float[nCols*nRows]();
    float* b2_1=new float[nCols*nRows]();
    float* x1_1=new float[nCols*nRows]();
    float* y1_1=new float[nCols*nRows]();
    float* x2_1=new float[nCols*nRows]();
    float* y2_1=new float[nCols*nRows]();
    
    float* Lab_2=new float[nCols*nRows]();
    float* xy_2=new float[nCols*nRows]();
    
    float* centerL1=new float[seedNum];
    float* centerL2=new float[seedNum];
    float* centera1=new float[seedNum];
    float* centera2=new float[seedNum];
    float* centerb1=new float[seedNum];
    float* centerb2=new float[seedNum];
    float* centerx1=new float[seedNum];
    float* centerx2=new float[seedNum];
    float* centery1=new float[seedNum];
    float* centery2=new float[seedNum];
    float* WSum=new float[seedNum];
    int* clusterSize=new int[seedNum];
    float* Lab_cc = new float[seedNum];
    float* xy_cc  = new float[seedNum];
    
    //BILATERAL
    float sigma_s = FLT_MAX;
    float sigma_c = 40;
    float *g_s = (float*) calloc(max(nRows,nCols),sizeof(float));
    float *g_r = (float*) calloc(256,sizeof(float));
    lookup_table2(g_s,g_r,max(nRows,nCols),sigma_s,sigma_c);
    
    int pw = 2;
    for(int i=0;i<nCols;i++) {
        for(int j=0;j<nRows;j++) {
            float count = 0;
            
            int pos_i = i*nRows+j;
            
            label[pos_i] = 0;
            
            unsigned char valr = R[pos_i];
            unsigned char valg = G[pos_i];
            unsigned char valb = B[pos_i];
            unsigned char val = (valr + valg + valb)/3;
            
            x1_1[pos_i] = (float) (x1[i][j]);
            y1_1[pos_i] = (float) y1[i][j];
            x2_1[pos_i] = (float) (x2[i][j]);
            y2_1[pos_i] = (float) (y2[i][j]);
            
            xy_2[pos_i] = (float) x1[i][j]*x1[i][j] +  x2[i][j]*x2[i][j] +  y1[i][j]*y1[i][j] +  y2[i][j]*y2[i][j];
            
            
            for (int dx=-pw; dx<=pw; dx++) {
                for (int dy=-pw; dy<=pw; dy++) {
                    if ((i+dx<nCols)&&(i+dx>=0)&(j+dy>=0)&(j+dy<nRows)){
                        
                        int pos_d = (i+dx)*nRows+(j+dy);
                        unsigned char val2r = R[pos_d];
                        unsigned char val2g = G[pos_d];
                        unsigned char val2b = B[pos_d];
                        unsigned char val2 = (val2r + val2g + val2b)/3;
                        float d_s = g_s[abs(dx)+abs(dy)];
                        float d_rr = g_r[abs(val-val2)];
                        float kk = d_rr*d_s;
                        L1_1[pos_i] += (float) (L1[i+dx][j+dy]*kk);
                        L2_1[pos_i] += (float) (L2[i+dx][j+dy]*kk);
                        a1_1[pos_i] += (float) (a1[i+dx][j+dy]*kk);
                        b1_1[pos_i] += (float) (b1[i+dx][j+dy]*kk);
                        a2_1[pos_i] += (float) (a2[i+dx][j+dy]*kk);
                        b2_1[pos_i] += (float) (b2[i+dx][j+dy]*kk);
                        Lab_2[pos_i] += (float) (L1[i+dx][j+dy]*L1[i+dx][j+dy] +
                                L2[i+dx][j+dy]*L2[i+dx][j+dy] +
                                a1[i+dx][j+dy]*a1[i+dx][j+dy] +
                                a2[i+dx][j+dy]*a2[i+dx][j+dy] +
                                b1[i+dx][j+dy]*b1[i+dx][j+dy] +
                                b2[i+dx][j+dy]*b2[i+dx][j+dy])*kk;
                        
                        
                        count += (float) kk;
                    }
                }
            }
            L1_1[pos_i]  /= count;
            L2_1[pos_i]  /= count;
            a1_1[pos_i]  /= count;
            a2_1[pos_i]  /= count;
            b1_1[pos_i]  /= count;
            b2_1[pos_i]  /= count;
            
            Lab_2[pos_i] /= count;
            
            
        }
    }
    
    
    //Initialization
    for(int i=0;i<seedNum;i++)
    {
        centerL1[i]=0;
        centerL2[i]=0;
        centera1[i]=0;
        centera2[i]=0;
        centerb1[i]=0;
        centerb2[i]=0;
        centerx1[i]=0;
        centerx2[i]=0;
        centery1[i]=0;
        centery2[i]=0;
        int x=seedArray[i].x;
        int y=seedArray[i].y;
        int minX=(x-StepX/4<=0)?0:x-StepX/4;
        int minY=(y-StepY/4<=0)?0:y-StepY/4;
        int maxX=(x+StepX/4>=nCols-1)?nCols-1:x+StepX/4;
        int maxY=(y+StepY/4>=nRows-1)?nRows-1:y+StepY/4;
        int Count=0;
        for(int j=minX;j<=maxX;j++)
            for(int k=minY;k<=maxY;k++)
            {
                Count++;
                
                centerL1[i]+=L1[j][k];
                centerL2[i]+=L2[j][k];
                centera1[i]+=a1[j][k];
                centera2[i]+=a2[j][k];
                centerb1[i]+=b1[j][k];
                centerb2[i]+=b2[j][k];
                centerx1[i]+=x1[j][k];
                centerx2[i]+=x2[j][k];
                centery1[i]+=y1[j][k];
                centery2[i]+=y2[j][k];
                
                label[k+j*nRows] = i;
                
            }
        centerL1[i]/=Count;
        centerL2[i]/=Count;
        centera1[i]/=Count;
        centera2[i]/=Count;
        centerb1[i]/=Count;
        centerb2[i]/=Count;
        centerx1[i]/=Count;
        centerx2[i]/=Count;
        centery1[i]/=Count;
        centery2[i]/=Count;
        
        Lab_cc[i] = centerL1[i]*centerL1[i] + centerL2[i]*centerL2[i] + centera1[i]*centera1[i] +
                centera2[i]*centera2[i] + centerb1[i]*centerb1[i] + centerb2[i]*centerb2[i];
        
        xy_cc[i]  = centerx1[i]*centerx1[i] + centerx2[i]*centerx2[i] + centery1[i]*centery1[i] + centery2[i]*centery2[i];
        
    }
    
    float *ss = (float *) calloc(seedNum,sizeof(float));
    
    //Init lab map
    for(int i=0;i<seedNum;i++) {
        int x=seedArray[i].x;
        int y=seedArray[i].y;
        int minX=(x-StepX/2<=0)?0:x-StepX/2;
        int minY=(y-StepY/2<=0)?0:y-StepY/2;
        int maxX=(x+StepX/2>=nCols-1)?nCols-1:x+StepX/2;
        int maxY=(y+StepY/2>=nRows-1)?nRows-1:y+StepY/2;
        for(int j=minX;j<=maxX;j++) {
            for(int k=minY;k<=maxY;k++) {
                label[k+j*nRows] = i;
            }
        }
        minX=(x-(StepX)<=0)?0:x-StepX;
        minY=(y-(StepY)<=0)?0:y-StepY;
        maxX=(x+(StepX)>=nCols-1)?nCols-1:x+StepX;
        maxY=(y+(StepY)>=nRows-1)?nRows-1:y+StepY;
        ss[i] = (centerx1[i] - x1_1[minX*nRows+minY])*(centerx1[i] - x1_1[minX*nRows+minY]) +
                (centery1[i] - y1_1[minX*nRows+minY])*(centery1[i] - y1_1[minX*nRows+minY]) +
                (centerx2[i] - x2_1[minX*nRows+minY])*(centerx2[i] - x2_1[minX*nRows+minY]) +
                (centery2[i] - y2_1[minX*nRows+minY])*(centery2[i] - y2_1[minX*nRows+minY]);
    }
    
    
    int * label_Kpm = (int*) calloc(nRows*nCols*Kpm,sizeof(int));
    float * dist_Kpm = (float*) calloc(nRows*nCols*Kpm,sizeof(float));
    
    //Thread stuff
    int thread_nbr = Kpm;
    
    //Thread argument structures
    pthread_t*   thread_list = (pthread_t*) calloc(thread_nbr, sizeof(pthread_t));
    pm_struct* thread_args = (pm_struct*)calloc(thread_nbr, sizeof(pm_struct));
    
    unsigned long next = 1789;
    
    //Launching of the THREADS
    for (int i=0; i < thread_nbr; i++) {
        //Thread arguments
        
        thread_args[i].next = next;
        pm_rand(&next);
        thread_args[i].Kpm_nbr = i;
        
        thread_args[i].R = R;
        thread_args[i].G = G;
        thread_args[i].B = B;
        
        thread_args[i].L1 = L1;
        thread_args[i].L2 = L2;
        thread_args[i].a1 = a1;
        thread_args[i].a2 = a2;
        thread_args[i].b1 = b1;
        thread_args[i].b2 = b2;
        thread_args[i].x1 = x1;
        thread_args[i].x2 = x2;
        thread_args[i].y1 = y1;
        thread_args[i].y2 = y2;
        thread_args[i].L1_1 = L1_1;
        thread_args[i].L2_1 = L2_1;
        thread_args[i].a1_1 = a1_1;
        thread_args[i].a2_1 = a2_1;
        thread_args[i].b1_1 = b1_1;
        thread_args[i].b2_1 = b2_1;
        thread_args[i].Lab_2 = Lab_2;
        
        thread_args[i].x1_1 = x1_1;
        thread_args[i].x2_1 = x2_1;
        thread_args[i].y1_1 = y1_1;
        thread_args[i].y2_1 = y2_1;
        thread_args[i].xy_2 = xy_2;
        
        thread_args[i].dist_out = dist_Kpm;
        thread_args[i].label_out = label_Kpm;
        thread_args[i].label_in = label;
        
        thread_args[i].nRows = nRows;
        thread_args[i].nCols = nCols;
        
        thread_args[i].seedArray_in = seedArray;
        thread_args[i].seedNum = seedNum;
        thread_args[i].StepX = StepX;
        thread_args[i].StepY = StepY;
        
        thread_args[i].centerL1_in = centerL1;
        thread_args[i].centerL2_in = centerL2;
        thread_args[i].centera1_in = centera1;
        thread_args[i].centera2_in = centera2;
        thread_args[i].centerb1_in = centerb1;
        thread_args[i].centerb2_in = centerb2;
        thread_args[i].centerx1_in = centerx1;
        thread_args[i].centerx2_in = centerx2;
        thread_args[i].centery1_in = centery1;
        thread_args[i].centery2_in = centery2;
        
        thread_args[i].Lab_cc_in = Lab_cc;
        thread_args[i].xy_cc_in = xy_cc;
        thread_args[i].ss_in = ss;
        thread_args[i].W = W;
        thread_args[i].iterationNum = iterationNum;
        
        thread_args[i].pw = patch_w;
        
        if (pthread_create(&thread_list[i], NULL, pm_core, &thread_args[i]))
            printf("Error creating a thread!\n");
        
    }
    
    /*Wait for all threads to end*/
    for (int j=0; j<thread_nbr; j++) {
        pthread_join(thread_list[j],NULL);
    }
    
    int ha = nRows;
    int wa = nCols;
    
    for (int i=0; i<ha*wa; i++)
        label[i] = label_Kpm[i];
    
    
    //refaire sans grosse structure
    //fusion of label maps
    float * sum_w = (float *) calloc(ha*wa*Kpm,sizeof(float));
    int * diff_labels_Kpm = (int *) malloc(ha*wa*Kpm*sizeof(int));
    for (int i=0; i<ha*wa*Kpm; i++) {
        diff_labels_Kpm[i] = -1;
    }
    
    for (int i=0; i<ha*wa; i++) {
        
        for (int j=0; j<Kpm; j++) {
            
            int lab = label_Kpm[i+ha*wa*j];
            
            int k = 0;
            while (k<Kpm) {
                if (diff_labels_Kpm[i+ha*wa*k] == lab) {
                    sum_w[i+ha*wa*k] += 1;
                    k = Kpm;
                }
                else {
                    if (diff_labels_Kpm[i+ha*wa*k] == -1) {
                        diff_labels_Kpm[i+ha*wa*k] = lab;
                        sum_w[i+ha*wa*k] += 1;
                        k = Kpm;
                    }
                }
                k += 1;
            }
        }
        int max = 0;
        for (int k=0; k<Kpm; k++) {
            if ((sum_w[i+ha*wa*k])>max) {
                max = sum_w[i+ha*wa*k];
                label[i] = diff_labels_Kpm[i+ha*wa*k];
            }
        }
    }
    
    
    //EnforceConnection
    int threshold=(nCols*nRows)/(seedNum*thresholdCoef);
    preEnforceConnectivity(label,nCols,nRows);
    EnforceConnectivity(L1_1,L2_1,a1_1,a2_1,b1_1,b2_1,x1_1,x2_1,y1_1,y2_1,W,label,threshold,nCols,nRows);
    
    
    free(diff_labels_Kpm);
    free(sum_w);
    free(label_Kpm);
    free(dist_Kpm);
    
    free(ss);
    free(thread_list);
    free(thread_args);
    
    free(g_r);
    free(g_s);
    
    delete []centerL1;
    delete []centerL2;
    delete []centera1;
    delete []centera2;
    delete []centerb1;
    delete []centerb2;
    delete []centerx1;
    delete []centerx2;
    delete []centery1;
    delete []centery2;
    delete []Lab_cc;
    delete []xy_cc;
    
    delete []L1_1;
    delete []L2_1;
    delete []a1_1;
    delete []a2_1;
    delete []b1_1;
    delete []b2_1;
    delete []y1_1;
    delete []x1_1;
    delete []y2_1;
    delete []x2_1;
    delete []Lab_2;
    delete []xy_2;
    delete []WSum;
    delete []clusterSize;
    
    return;
}

#endif
