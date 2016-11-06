#ifndef FFT_H
#define FFT_H
void fft(double *re_src, double *im_src, double *re_dst, double *im_dst, int N);
void ifft(double *re_src, double *im_src, double *re_dst, double *im_dst, int N);
#endif
