#include<stdio.h>
#include "fft.h"
int main(void)
{
    double re_src[] = {1, 2, 3.6, -1, 0, 2, 3, 4};
    double im_src[] = {2, 4, -1, 0, 6.7, 2, 0, -1};
    double re_dst[8];
    double im_dst[8];
    double re_recovered[8];
    double im_recovered[8];
    fft(re_src, im_src, re_dst, im_dst, 8);
    ifft(re_dst, im_dst, re_recovered, im_recovered, 8);
    for(int i = 0; i < 8; ++i){
        printf("RE::%.5f %.5f\n", re_src[i], re_recovered[i]);
        printf("IM::%.5f %.5f\n", im_src[i], im_recovered[i]);
    }
    return 0;
}
