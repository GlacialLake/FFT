#include "fft.h"
#include<math.h>
#include<string.h>
#include<stdio.h>
#define PI 3.141592653589793

void fft(double *re_src, double *im_src, double *re_dst, double *im_dst, int N)
{
    int b = 1;
    while(!(N & b)){
        b <<= 1;
    }
    if(b != N){
        return;
    }

    memcpy(re_dst, re_src, N * sizeof(double));
    memcpy(im_dst, im_src, N * sizeof(double));

    int M = (int)((double)log(N) / log(2));

    //change the indexes
    int NV2 = N >> 1;
    int new_ord = 0;
    int bit;
    double swap_tmp;
    for(int i = 0; i < N - 1; ++i){
        if(i < new_ord){
           swap_tmp = re_dst[new_ord];
           re_dst[new_ord] = re_dst[i];
           re_dst[i] = swap_tmp;
           swap_tmp = im_dst[new_ord];
           im_dst[new_ord] = im_dst[i];
           im_dst[i] = swap_tmp;
        }
        bit = NV2;
        while(bit <= new_ord){
            new_ord -= bit;
            bit >>= 1;
        }
        new_ord += bit;
    }

    int prev_dft_len;
    double re_coeff, im_coeff;
    double re_coeff_step, im_coeff_step;
    double re_diff, im_diff;
    double re_new_coeff, im_new_coeff;
    int bf_up, bf_down;
    
    int dft_len = 1;
    for(int l = 0; l < M; ++l){
        dft_len <<= 1;
        prev_dft_len = dft_len >> 1;

        //the coeff
        re_coeff = 1;
        im_coeff = 0;

        //the increase scale of the coeffs on the current layer
        re_coeff_step = cos(PI / prev_dft_len);
        im_coeff_step = -sin(PI / prev_dft_len);

        //one group refers to the butterflies on the current layer which has the same coeff
        for(int group = 0; group < prev_dft_len; ++group){

            //bf_up refers to the top fork's index of the current butterfly
            //bf_down refers to the bottom fork's index of the current butterfly
            for(bf_up = group; bf_up < N; bf_up += dft_len){
                bf_down = bf_up + prev_dft_len;
                re_diff = re_dst[bf_down] * re_coeff - im_dst[bf_down] * im_coeff;
                im_diff = re_dst[bf_down] * im_coeff + im_dst[bf_down] * re_coeff;
                re_dst[bf_down] = re_dst[bf_up] - re_diff;
                im_dst[bf_down] = im_dst[bf_up] - im_diff;
                re_dst[bf_up] += re_diff;
                im_dst[bf_up] += im_diff;
            }

            //update the coeffs on the current layer
            re_new_coeff = re_coeff * re_coeff_step - im_coeff * im_coeff_step;
            im_new_coeff = re_coeff * im_coeff_step + im_coeff * re_coeff_step;
            re_coeff = re_new_coeff;
            im_coeff = im_new_coeff;
        }
    }

}

void ifft(double *re_src, double *im_src, double *re_dst, double *im_dst, int N)
{
    int b = 1;
    while(!(N & b)){
        b <<= 1;
    }
    if(b != N){
        return;
    }

    double *ptr, *ptr1, *ptr2;
    int n;
    ptr = im_src;
    n = N;
    while(n--){
        *ptr = -*ptr;
        ++ptr;
    }

    //use the previous function as a tool
    fft(re_src, im_src, re_dst, im_dst, N);
    
    ptr = re_dst;
    ptr1 = im_dst;
    ptr2 = im_src;
    n = N;
    while(n--){
        *ptr = (*ptr) / N;
        *ptr1 = (-*ptr1) / N;
        
        //restore the origin array
        *ptr2 = -*ptr2;
        
        ++ptr;
        ++ptr1;
        ++ptr2;
    }
}
