


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "mex.h"
#include "matrix.h"

#ifndef MAX
#define MAX(a, b) ((a)>(b)?(a):(b))
#endif

#ifndef MIN
#define MIN(a, b) ((a)<(b)?(a):(b))
#endif

#define EPSILON 0.0000001

#ifndef enc_type
#define enc_type unsigned char
#endif

#define omega 4      //5
#define mag_init 10  //15
#define nbr_iter 3   //3




int
        pm_rand(unsigned long* next)
{
    (*next) = (*next) * 1103515245 + 12345;
    return (unsigned int)((*next)/65536) % 32768;
}


// Measure distance between entire 2 patches,
// terminating early if we exceed a cutoff distance
float dist_2D(float *l, float *a, float *b, float *l2, float * a2, float * b2, int ax, int ay, int bx, int by,
        int ha, int hb, int size_a, int size_b, int pw, float cutoff,
        int cy, int cx, float beta) {
    
    float ans = 0;
    int dx, dy;
    float ld = 0, ad = 0, bd = 0, ld2 = 0, ad2 = 0, bd2 = 0;
    
    float count = 0;
    for (dx = -pw; dx <= pw; dx++) {
        for (dy = -pw; dy <= pw; dy++) {
            
            ld = l[(ax+dx)*ha+ay+dy] - l[(bx+dx)*hb+by+dy];
            ad = a[(ax+dx)*ha+ay+dy] - a[(bx+dx)*hb+by+dy];
            bd = b[(ax+dx)*ha+ay+dy] - b[(bx+dx)*hb+by+dy];
            ld2 = l2[(ax+dx)*ha+ay+dy] - l2[(bx+dx)*hb+by+dy];
            ad2 = a2[(ax+dx)*ha+ay+dy] - a2[(bx+dx)*hb+by+dy];
            bd2 = b2[(ax+dx)*ha+ay+dy] - b2[(bx+dx)*hb+by+dy];
            
            ans += ld*ld + ad*ad + bd*bd + ld2*ld2 + ad2*ad2 + bd2*bd2;
            
            count += 1;
            
        }
        //if (ans >= cutoff) { return cutoff; }
    }
    ans /= count;
    
    
    //pondération distance au centre
    float sigma = 10;
//     float dist_s_c = 2 - exp(- ((cx-bx)*(cx-bx) + (cy-by)*(cy-by))/(sigma*sigma));
    float dist_s_c = 1 + (1 - exp(- ((cx-bx)*(cx-bx) + (cy-by)*(cy-by))/(sigma*sigma)))*5;
    
    
    ans += dist_s_c*beta; //1000000;
    
// //     
//     if ((ay == 215) && (ax == 520)) {
//         printf("dist %f, distworst %f, %d, %d, %d, %d, %d, %f\n", ans, cutoff, pw, ax, bx, ay, by, dist_s_c);
// //         mexEvalString("drawnow");
//     }
    
    
    return ans;
}




// Compare two distance and change the optimal ptr if the new distance is
// lower than the previous best
void improve_guess_2D_queue(float *l, float *a, float *b, float *l2, float * a2, float * b2, int ax, int ay,
        int *xbest_ptr, int *ybest_ptr, float *distbest_ptr,
        int bx, int by, int ha, int hb, int pw,
        int size_a, int size_b, int * nnf, float *nnfd, int * max_dist, int Kpm, float off,
        int cy, int cx, float* max_dist_pm, int * min_dist, float beta) {
    
    //Distance computation
    float dist = dist_2D(l, a, b, l2, a2, b2, ax, ay, bx, by, ha, hb, size_a, size_b, pw, *distbest_ptr, cy, cx, beta);
    dist += off;
    
 
    
    //Comparison
    if (dist < *distbest_ptr) {
//         *distbest_ptr = dist;
//         *xbest_ptr = bx;
//         *ybest_ptr = by;
        
        int pos = ay+ha*ax;
        
        max_dist_pm[pos] += dist - *distbest_ptr;
        
        //Map updates
        nnf[pos + max_dist[pos]*size_a*2]           = bx;
        nnf[pos + max_dist[pos]*size_a*2 + size_a]  = by;
        nnfd[pos + max_dist[pos]*size_a]            = dist;
        
        
        if (dist < nnfd[pos + min_dist[pos]*size_a]) {
                min_dist[pos] = max_dist[pos];
        }
        
        //Find max pos to replace
        float max_dist_k = 0;
        for (int k=0; k<Kpm; k++) {
            if (nnfd[pos + k*size_a] > max_dist_k) {
                max_dist_k = nnfd[pos + k*size_a];
                max_dist[pos] = k;
            }
        }
        
        *xbest_ptr = nnf[pos + max_dist[pos]*size_a*2];
        *ybest_ptr = nnf[pos + max_dist[pos]*size_a*2 + size_a];
        *distbest_ptr = nnfd[pos + max_dist[pos]*size_a];
        
    }
    
}



void  patchmatch_queue(float *l, float * a, float * b, float *l2, float * a2, float * b2,
        int miny, int maxy, int minx, int maxx, int Kpm,
        int ha, int wa, int * klabels, int lab, float * nnfd, int pw,
        int cy, int cx, float * max_dist_pm, int *nnf, int * max_dist, float beta, 
        unsigned long next) {
    
    
    int * min_dist = (int *) calloc(ha*wa,sizeof(int));
    
    int hb = ha;
    int wb = wa;
    
    // Initialize with random nearest neighbor field (NNF)
    // Effective width and height (possible upper left corners of patches)
    int aew = maxx-pw, aeh = maxy-pw;
    int bew = maxx-pw, beh = maxy-pw;
    
    int size_a = ha*wa;
    int size_b = hb*wb;
    
    int pos, pos_shift;
    int ax, ay, bx, by, ax_off;
    int iter = 0;
    int ystart, yend, ychange, xstart, xend, xchange;
    int xbest, ybest;
    float distbest;
    int distbest_temp;
    float *distbest_ptr;
    int *xbest_ptr;
    int *ybest_ptr;
    int rand_val;
    int xp, yp, xmin, xmax, ymin, ymax;
    int i;
  
    
    float offset = 999;
    float off = 0;
    
    // RANDOM INITIALIZATION
    for (ax = minx+pw; ax < aew; ax++) {
        
        ax_off = ha*ax;
        
        for (ay = miny+pw; ay < aeh; ay++) {
            
            pos = ax_off+ay;
            
            max_dist_pm[pos] = 0;
            
            
            xmin = MAX(cx-mag_init, pw);
            xmax = MIN(cx+mag_init+1, bew-1);
            
            ymin = MAX(cy-mag_init, pw);
            ymax = MIN(cy+mag_init+1, beh-1);
            
            
            float max_dist_k = -1;
            float min_dist_k = 999;
            for (int k=0; k<Kpm; k++) {
                
                int loop = 1;
//
                while (loop) {
                    
                    
                    //random match
                    rand_val = pm_rand(&next);
                    bx = xmin+rand_val%(xmax-xmin);
                    rand_val = pm_rand(&next);
                    by = ymin+rand_val%(ymax-ymin);
                    
//                      if ((lab == 1) && (klabels[by+ha*bx] != 0)) {
//                     printf("lab %d\n", klabels[by+ha*bx]);
//                 mexEvalString("drawnow");
//
                    if ((klabels[by+ha*bx] == lab) || (loop > 10)){//>= 0) //
                        loop = 0;
//                     off = 0;
//                 else
//                     off = offset;
//                         loop = 0;
//                 }
                        
                        //Map init
                        nnf[pos + k*size_a*2] = bx;
                        nnf[pos + k*size_a*2 + size_a] = by;
                        float dist_k = (float) dist_2D(l, a, b, l2, a2, b2, ax, ay, bx, by, ha, hb, size_a, size_b, pw, FLT_MAX, cy, cx, beta);
                        dist_k += off;
                        nnfd[pos + k*size_a] = dist_k;
                        
                        max_dist_pm[pos] += dist_k;
                        
                        if (max_dist_pm[pos] < 0) {
                            printf("%f\n", max_dist_pm[pos]);
                            mexEvalString("drawnow");
                        }
                        
                        //find max dist
                        if (dist_k > max_dist_k) {
                            max_dist_k = dist_k;
                            max_dist[pos] = k;
                        }
                        if (dist_k < min_dist_k) {
                            min_dist_k = dist_k;
                            min_dist[pos] = k;
                        }
                        
                        
                        //////
//                         max_dist_pm[pos] = 0.1;
                        /////
                        
                        
                        
                    }
                    else {
                       loop += 1;   
                    }
                    
                }
            }
            
            
        }
    }
    
//     printf("out init\n");
//                 mexEvalString("drawnow");
    
    
    while (iter < nbr_iter) {
        
//             printf("iter %d\n", iter);
//                 mexEvalString("drawnow");
        
        // In each iteration, improve the NNF, by looping in scanline or
        // reverse-scanline order
        ystart = miny+1+pw; yend = aeh; ychange = 1;
        xstart = minx+1+pw; xend = aew; xchange = 1;
        
        if (iter % 2 == 1) {
            xstart = xend-1; xend = minx+pw-1; xchange = -1;
            ystart = yend-1; yend = miny+pw-1; ychange = -1;
        }
        
        for (ay = ystart; ay != yend; ay += ychange) {
            
//                         printf("ay %d\n", ay);
//                 mexEvalString("drawnow");
            
            
            for (ax = xstart; ax != xend; ax += xchange) {
                
                //Current position in a
                pos = ha*ax+ay;
                
                // Current (best) guess
                xbest = nnf[pos + min_dist[pos]*size_a*2];
                ybest = nnf[pos + min_dist[pos]*size_a*2 + size_a];
                distbest = nnfd[pos + max_dist[pos]*size_a];
                distbest_temp = distbest;
                distbest_ptr = &distbest;
                xbest_ptr = &xbest;
                ybest_ptr = &ybest;
                
                // PROPAGATION: Improve current guess by trying instead
                // correspondences from left and above (below and right on
                // odd iterations)
                
                //XSHIFT
                if ((unsigned) (ax - xchange-pw) < (unsigned) aew-pw) {
                    pos_shift = ha*(ax-xchange)+ay;
                    xp = nnf[pos_shift]  + xchange;
                    yp = nnf[pos_shift + size_a];
                    if ((unsigned) xp-pw < (unsigned) bew-pw) {
                        if ((xp !=  *xbest_ptr)  || (yp !=  *ybest_ptr)) {
                            if (klabels[yp + xp*hb] == lab) { //>= 0) {
//                                 off = 0;
//                             else
//                                 off = offset;
                                if ( (abs(yp-ay)>=omega) || (abs(xp-ax)>= omega) )
                                    improve_guess_2D_queue(l, a, b, l2, a2, b2, ax, ay, xbest_ptr, ybest_ptr,
                                            distbest_ptr, xp, yp, ha,  hb, pw, size_a, size_b, nnf, nnfd, max_dist, Kpm, off, cy, cx, max_dist_pm, min_dist, beta);
                            }
                        }
                    }
                }
                
                
                //YSHIFT
                if ((unsigned) (ay - ychange-pw) < (unsigned) aeh-pw) {
                    pos_shift = ha*ax+ay -ychange;
                    xp = nnf[pos_shift];
                    yp = nnf[pos_shift + size_a] + ychange;
                    if ((unsigned) yp-pw < (unsigned) beh-pw) {
                        if ((xp !=  *xbest_ptr)  || (yp !=  *ybest_ptr)) {
                            
                            if (klabels[yp + xp*hb] == lab)  { //>= 0) {
//                                 off = 0;
//                             else
//                                 off = offset;
                                if ( (abs(yp-ay)>=omega) || (abs(xp-ax)>= omega) )
                                    improve_guess_2D_queue(l, a, b, l2, a2, b2, ax, ay, xbest_ptr, ybest_ptr,
                                            distbest_ptr, xp, yp, ha, hb, pw, size_a, size_b, nnf, nnfd, max_dist, Kpm, off, cy, cx, max_dist_pm, min_dist, beta);
                            }
                        }
                    }
                }
                
                // RANDOM SEARCH: Improve current guess by searching in
                // boxes of exponentially decreasing size around the
                // current position
                
                // Sampling window
                for (int mag = mag_init; mag >= 1; mag /= 2) {
                    
//                     xbest = *xbest_ptr;
//                     ybest = *ybest_ptr;
                    
//                     xmin = MAX(xbest-mag, minx+pw);
//                     xmax = MIN(xbest+mag+1, bew-1);
                    xmin = MAX(cx-mag, minx+pw);
                    xmax = MIN(cx+mag+1, bew-1);
                    if(xmin == xmax) continue;
                    
                    ymin = MAX(cy-mag, miny+pw);
                    ymax = MIN(cy+mag+1, beh-1);
                    if(ymin == ymax) continue;
                    
                    //Random match
                    xp = (int) xmin+pm_rand(&next)%MAX(xmax-xmin, 1);
                    yp = (int) ymin+pm_rand(&next)%MAX(ymax-ymin, 1);
                    
                    if (((unsigned) yp-pw < (unsigned) beh)
                    && ((unsigned) xp-pw < (unsigned) bew)) {
                        
                        if (klabels[yp + xp*hb] == lab) {  //>= 0) {
//                             off = 0;
//                         else
//                             off = offset;
                            
                            if ( (abs(yp-ay)>=omega) || (abs(xp-ax)>= omega) )
                                improve_guess_2D_queue(l, a, b, l2, a2, b2, ax, ay, xbest_ptr, ybest_ptr,
                                        distbest_ptr, xp, yp, ha, hb, pw, size_a, size_b, nnf, nnfd, max_dist, Kpm, off, cy, cx, max_dist_pm, min_dist, beta);
                        }
                    }
                }
                
//                 //Map updates
//                 nnf[pos + max_dist[pos]*size_a*2]           = *xbest_ptr;
//                 nnf[pos + max_dist[pos]*size_a*2 + size_a]  = *ybest_ptr;
//                 nnfd[pos + max_dist[pos]*size_a]         = *distbest_ptr;
//
//                 //Find max pos to replace
//                 float max_dist_k = 0;
//                 for (int k=0; k<Kpm; k++) {
//                     if (nnfd[pos + k*size_a] > max_dist_k) {
//                          max_dist_k = nnfd[pos + k*size_a];
//                          max_dist[pos] = k;
//                     }
//                 }
                
//                 if ((ay==20) && (ax==20)) {
//                     printf("%d, %d, %f, %d\n", *xbest_ptr, *ybest_ptr, *distbest_ptr, max_dist[ay+ax*ha]);
//                     mexEvalString("drawnow");
//                     for (int pp=0; pp<Kpm; pp++) {
//                         printf("PM in %f\n",nnfd[ay+ax*ha + pp*ha*wa]);
//                         mexEvalString("drawnow");
//                     }
//                 }
                
            }
            
            
            
        }
        
        
        iter++;
    }
    
    
    free(min_dist);
    
    
}



// Compare two distance and change the optimal ptr if the new distance is
// lower than the previous best
void improve_guess_2D(float *l, float *a, float *b, float *l2, float * a2, float * b2, int ax, int ay,
        int *xbest_ptr, int *ybest_ptr, float *distbest_ptr,
        int bx, int by, int ha, int hb, int pw,
        int size_a, int size_b, int * nnf, float *nnfd, int Kpm, float off,
        int cy, int cx, float* max_dist_pm, float beta) {
    
    //Distance computation
    float dist = dist_2D(l, a, b, l2, a2, b2, ax, ay, bx, by, ha, hb, size_a, size_b, pw, *distbest_ptr, cy, cx, beta);
    dist += off;
    
    
    //Comparison
    if (dist < *distbest_ptr) {
        *distbest_ptr = dist;
        *xbest_ptr = bx;
        *ybest_ptr = by;
        
    }
    
}





typedef struct{
    float *l;
    float *a;
    float *b;
    float *l2;
    float *a2;
    float *b2;
    int miny;
    int minx;
    int maxy;
    int maxx;
    int Kpm;
    int ha;
    int wa;
    int * klabels;
    int lab;
    int pw;
    float * max_dist_pm;
    int thread_nbr;
    int *nnf;
    float * nnfd;
    unsigned long next;
    float beta; 
    int miny_i;
    int minx_i;
    int maxy_i;
    int maxx_i;
}pm_struct;




void*
        pm_core(void *arg)
{
    pm_struct inputs = *(pm_struct*) arg;
    float *l     = inputs.l;
    float *a     = inputs.a;
    float *b     = inputs.b;
    float *l2    = inputs.l2;
    float *a2    = inputs.a2;
    float *b2    = inputs.b2;
    int miny     = inputs.miny;
    int minx     = inputs.minx;
    int maxy     = inputs.maxy;
    int maxx     = inputs.maxx;
    int miny_i     = inputs.miny_i;
    int minx_i     = inputs.minx_i;
    int maxy_i     = inputs.maxy_i;
    int maxx_i     = inputs.maxx_i;
    int Kpm      = inputs.Kpm;
    int ha       = inputs.ha;
    int wa       = inputs.wa;
    int * klabels  = inputs.klabels;
    int lab      = inputs.lab;
    int pw       = inputs.pw;
    float * max_dist_pm = inputs.max_dist_pm;
    int thread_nbr = inputs.thread_nbr;
    int * nnf = inputs.nnf;
    float * nnfd = inputs.nnfd;
    unsigned long next = inputs.next;
    float beta = inputs.beta;
    
    int cx = minx_i + (maxx_i-minx_i)/2;
    int cy = miny_i + (maxy_i-miny_i)/2;
    
    int hb = ha;
    int wb = wa;
    
    // Initialize with random nearest neighbor field (NNF)
    // Effective width and height (possible upper left corners of patches)
    int aew = maxx-pw, aeh = maxy-pw;
    int bew = maxx_i-pw, beh = maxy_i-pw;
    
    int size_a = ha*wa;
    int size_b = hb*wb;
    
    int pos, pos_shift;
    int ax, ay, bx, by, ax_off;
    int iter = 0;
    int ystart, yend, ychange, xstart, xend, xchange;
    int xbest, ybest;
    float distbest;
    int distbest_temp;
    float *distbest_ptr;
    int *xbest_ptr;
    int *ybest_ptr;
    int rand_val;
    int xp, yp, xmin, xmax, ymin, ymax;
    int i;
    
    
    int off_nnf = thread_nbr*size_a;
    
    
    float offset = 999;
    float off = 0;
    
    // RANDOM INITIALIZATION
    for (ax = minx+pw; ax < aew; ax++) {
        
        ax_off = ha*ax;
        
        for (ay = miny+pw; ay < aeh; ay++) {
            
            pos = ax_off+ay;
            
            xmin = MAX(cx-mag_init, minx_i);
            xmax = MIN(cx+mag_init+1, bew-1);
            
            ymin = MAX(cy-mag_init, miny_i);
            ymax = MIN(cy+mag_init+1, beh-1);
            
            
            int loop = 1;
//
            while (loop) {
                
                //random match
                rand_val = pm_rand(&next);
                bx = xmin+rand_val%(xmax-xmin);
                rand_val = pm_rand(&next);
                by = ymin+rand_val%(ymax-ymin);
                
//                      if ((lab == 1) && (klabels[by+ha*bx] != 0)) {
//                     printf("lab %d\n", klabels[by+ha*bx]);
//                 mexEvalString("drawnow");
//
                if ((klabels[by+ha*bx] == lab) || (loop > 10)) {//>= 0) //
                    loop = 0;
//                     off = 0;
//                 else
//                     off = offset;
//                         loop = 0;
//                 }
                    
                    //Map init
                    nnf[pos + off_nnf*2] = bx;
                    nnf[pos + off_nnf*2 + size_a] = by;
                    float dist_k = (float) dist_2D(l, a, b, l2, a2, b2, ax, ay, bx, by, ha, hb, size_a, size_b, pw, FLT_MAX, cy, cx, beta);
                    dist_k += off;
                    nnfd[pos + off_nnf] = dist_k;
                    
                    
                }
                else {
                    loop += 1;
                }
                
            }
            
            
        }
    }
    
//     printf("out init\n");
//                 mexEvalString("drawnow");
    
    
    while (iter <  nbr_iter) {
        
//             printf("iter %d\n", iter);
//                 mexEvalString("drawnow");
        
        // In each iteration, improve the NNF, by looping in scanline or
        // reverse-scanline order
        ystart = miny+1+pw; yend = aeh; ychange = 1;
        xstart = minx+1+pw; xend = aew; xchange = 1;
        
        if (iter % 2 == 1) {
            xstart = xend-1; xend = minx+pw-1; xchange = -1;
            ystart = yend-1; yend = miny+pw-1; ychange = -1;
        }
        
        for (ay = ystart; ay != yend; ay += ychange) {
            
//                         printf("ay %d\n", ay);
//                 mexEvalString("drawnow");
            
            
            for (ax = xstart; ax != xend; ax += xchange) {
                
                //Current position in a
                pos = ha*ax+ay;
                
                // Current (best) guess
                xbest = nnf[pos + off_nnf*2];
                ybest = nnf[pos + off_nnf*2 + size_a];
                distbest = nnfd[pos + off_nnf];
                distbest_temp = distbest;
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
//                     if ((unsigned) xp-pw < (unsigned) bew-pw) {
                    if ((xp-pw > minx_i) && (xp+pw<bew-1)) {
                        if ((xp !=  *xbest_ptr)  || (yp !=  *ybest_ptr)) {
                            if (klabels[yp + xp*hb] == lab) { //>= 0) {
//                                 off = 0;
//                             else
//                                 off = offset;
                                if ( (abs(yp-ay)>=omega) || (abs(xp-ax)>= omega) )
                                    improve_guess_2D(l, a, b, l2, a2, b2, ax, ay, xbest_ptr, ybest_ptr,
                                            distbest_ptr, xp, yp, ha,  hb, pw, size_a, size_b, nnf, nnfd, Kpm, off, cy, cx, max_dist_pm, beta);
                            }
                        }
                    }
                }
                
                
                //YSHIFT
                if ((unsigned) (ay - ychange-pw) < (unsigned) aeh-pw) {
                    pos_shift = ha*ax+ay -ychange;
                    xp = nnf[pos_shift + off_nnf*2];
                    yp = nnf[pos_shift + off_nnf*2 + size_a] + ychange;
//                     if ((unsigned) yp-pw < (unsigned) beh-pw) {
                     if ((yp-pw > miny_i) && (yp+pw<beh-1)) {
                        if ((xp !=  *xbest_ptr)  || (yp !=  *ybest_ptr)) {
                            
                            if (klabels[yp + xp*hb] == lab)  { //>= 0) {
//                                 off = 0;
//                             else
//                                 off = offset;
                                if ( (abs(yp-ay)>=omega) || (abs(xp-ax)>= omega) )
                                    improve_guess_2D(l, a, b, l2, a2, b2, ax, ay, xbest_ptr, ybest_ptr,
                                            distbest_ptr, xp, yp, ha, hb, pw, size_a, size_b, nnf, nnfd, Kpm, off, cy, cx, max_dist_pm, beta);
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
                    
                    xmin = MAX(xbest-mag, minx_i);
                    xmax = MIN(xbest+mag+1, bew-1);
//                     xmin = MAX(cx-mag, minx+pw);
//                     xmax = MIN(cx+mag+1, bew-1);
                    if(xmin == xmax) continue;
                    
                    ymin = MAX(ybest-mag, miny_i);
                    ymax = MIN(ybest+mag+1, beh-1);
                    if(ymin == ymax) continue;
                    
                    //Random match
                    xp = (int) xmin+pm_rand(&next)%MAX(xmax-xmin, 1);
                    yp = (int) ymin+pm_rand(&next)%MAX(ymax-ymin, 1);
                    
//                     if (((unsigned) yp-pw < (unsigned) beh) && ((unsigned) xp-pw < (unsigned) bew)) {
                    
                    
                     if ((xp-pw > minx_i) && (xp+pw<bew-1) && (yp-pw > miny_i) && (yp+pw<beh-1)) {
                    
                        if (klabels[yp + xp*hb] == lab) {  //>= 0) {
//                             off = 0;
//                         else
//                             off = offset;
                            
                            if ( (abs(yp-ay)>=omega) || (abs(xp-ax)>= omega) )
                                improve_guess_2D(l, a, b, l2, a2, b2, ax, ay, xbest_ptr, ybest_ptr,
                                        distbest_ptr, xp, yp, ha, hb, pw, size_a, size_b, nnf, nnfd, Kpm, off, cy, cx, max_dist_pm, beta);
                        }
                    }
                }
                
                //Map updates
                nnf[pos  + off_nnf*2]           = *xbest_ptr;
                nnf[pos  + off_nnf*2 + size_a]  = *ybest_ptr;
                nnfd[pos + off_nnf]             = *distbest_ptr;
//
//                 //Find max pos to replace
//                 float max_dist_k = 0;
//                 for (int k=0; k<Kpm; k++) {
//                     if (nnfd[pos + k*size_a] > max_dist_k) {
//                          max_dist_k = nnfd[pos + k*size_a];
//                          max_dist[pos] = k;
//                     }
//                 }
                
//                 if ((ay==20) && (ax==20)) {
//                     printf("%d, %d, %f, %d\n", *xbest_ptr, *ybest_ptr, *distbest_ptr, max_dist[ay+ax*ha]);
//                     mexEvalString("drawnow");
//                     for (int pp=0; pp<Kpm; pp++) {
//                         printf("PM in %f\n",nnfd[ay+ax*ha + pp*ha*wa]);
//                         mexEvalString("drawnow");
//                     }
//                 }
                
            }
            
            
            
        }
        
        
        iter++;
    }
    
    
    
    pthread_exit(0);
    
}



void  patchmatch(float *l, float * a, float * b, float *l2, float * a2, float * b2,
        int miny, int maxy, int minx, int maxx, 
        int ha, int wa, int Kpm, int * klabels, int lab, float * nnfd, int pw,
        float * max_dist_pm, int *nnf, int * max_dist, int use_mthread, int use_DIR, float beta) {
    
    int cx = minx + (maxx-minx)/2;
    int cy = miny + (maxy-miny)/2;
    
        
    unsigned long next = 1789;
    srand(time(NULL));
    
    if (use_mthread == 0) {
//         printf("hello\n");
//         mexEvalString("drawnow");
           patchmatch_queue(l, a, b, l2, a2, b2, miny, maxy, minx, maxx, Kpm, ha, wa, klabels, lab, nnfd, pw,
                            cy, cx, max_dist_pm, nnf, max_dist, beta, next);
    }
    else {
    
    
    int thread_nbr = Kpm;
    
    //Thread argument structures
    pthread_t* thread_list = (pthread_t*) calloc(thread_nbr, sizeof(pthread_t));
    pm_struct* thread_args = (pm_struct*)calloc(thread_nbr, sizeof(pm_struct));
    
    int minyi,minxi,maxyi,maxxi;
    

    //Launching of the THREADS
    for (int i=0; i < thread_nbr; i++) {
        
        //Thread arguments
        thread_args[i].l = l;
        thread_args[i].a = a;
        thread_args[i].b = b;
        thread_args[i].l2 = l2;
        thread_args[i].a2 = a2;
        thread_args[i].b2 = b2;
        
        thread_args[i].Kpm = Kpm;
        thread_args[i].ha = ha;
        thread_args[i].wa = wa;
        thread_args[i].klabels = klabels;
        thread_args[i].lab = lab;
        thread_args[i].pw = pw;
        thread_args[i].max_dist_pm = max_dist_pm;
        thread_args[i].nnf = nnf;
        thread_args[i].nnfd = nnfd;
        thread_args[i].thread_nbr = i;
        thread_args[i].beta = beta;
        
//         next = pm_rand(&next);
//         thread_args[i].next = next;
        thread_args[i].next = i;
    
        if (use_DIR==0) {
            minyi = MAX(miny-pw,pw);
            maxyi = MIN(maxy+pw,ha);
            minxi = MAX(minx-pw,pw);
            maxxi = MIN(maxx+pw,wa);
        }
        
        else {
            switch ((int) i/4) {
                case 0: //top left search
                    minyi = MAX(miny-pw,pw);
                    maxyi = cy+pw;
                    minxi = MAX(minx-pw,pw);
                    maxxi = cx+pw;
                    break;
                case 1: //top right search
                    minyi = MAX(miny-pw,pw);
                    maxyi = cy+pw;
                    minxi = cx-pw;
                    maxxi = MIN(maxx+pw,wa);
                    break;
                case 2: //bottom left search
                    minyi = cy-pw;
                    maxyi = MIN(maxy+pw,ha);
                    minxi = MAX(minx-pw,pw);
                    maxxi = cx+pw;
                    break;
                case 3: //bottom right search
                    minyi = cy-pw;
                    maxyi = MIN(maxy+pw,ha);
                    minxi = cx-pw;
                    maxxi = MIN(maxx+pw,wa);
                    break;
                default:
                    printf("problem\n");
            }
            
            
        }
        
        thread_args[i].miny = miny;
        thread_args[i].minx = minx;
        thread_args[i].maxy = maxy;
        thread_args[i].maxx = maxx;
        
        thread_args[i].miny_i = minyi;
        thread_args[i].minx_i = minxi;
        thread_args[i].maxy_i = maxyi;
        thread_args[i].maxx_i = maxxi;
        
        
        if (pthread_create(&thread_list[i], NULL, pm_core, &thread_args[i]))
            printf("Error creating a thread!\n");
        
    }
    
    /*Wait for all threads to end*/
    for (int j=0; j<thread_nbr; j++) {
        pthread_join(thread_list[j],NULL);
    }
    
    
    for (int y=miny; y<maxy; y++) {
        for (int x=minx; x<maxx; x++) {
            max_dist_pm[y+x*ha] = 0;
            for (int k=0; k<Kpm; k++)
                max_dist_pm[y+x*ha] += nnfd[y+x*ha + k*ha*wa];
            max_dist_pm[y+x*ha] /= Kpm;
        }
    }
    
        
    
    free(thread_args);
    free(thread_list);
    
    
    }

    
    
    
}









