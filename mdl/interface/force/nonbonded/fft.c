/* ZFFT http://fy.chalmers.se/~tfylb/ZFFT.html*/
/* Copyright (c) 2000 by Lennart Bengtsson */
/* Permission to use, copy, modify, distribute, and sell this software and its */
/* documentation for any purpose is hereby granted without fee, provided that */
/* the above copyright notice appear in all copies and that both the copyright */
/* notice and this permission notice appear in supporting documentation. */
/* No representations are made about the suitability of this software for */
/* any purpose.  It is provided "as is" without express or implied warranty. */
/* Note also that this kind of optimization violates the fortran standard */
/* so the compiled code should be tested carefully. So far ZFFT has been */
/* tested on SUN, SGI and IBM servers. */
#ifndef HAVE_FFT

#ifdef __cplusplus
extern "C" {
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct {
  double re;
  double im;
} zomplex;

#ifndef min
#define min(a,b)               ( ((a) < (b)) ? (a) : (b) )
#endif
#ifndef max
#define max(a,b)               ( ((a) > (b)) ? (a) : (b) )
#endif

/* Table of constant values */
static  int const1 = 1;
static  int const16384 = 16384;
static  int const0 = 0;
static  int constTrue = 1;
static  int constFalse = 0;

int radix2f_(double *a1,double * a2, double *b1,double * b2,int * b,int *ldb, double *c2,double * d2)
{
  /* System generated locals */
  int b1_dim2, b1_offset, b2_dim2, b2_offset, i__1;

    /* Local variables */
  int i__;
  double i1, i2, j2, r1, r2, s2;

    /* Parameter adjustments */
  a2 -= 3;
  a1 -= 3;
  b2_dim2 = *ldb;
  b2_offset = 1 + 2 * (1 + b2_dim2 * 1);
  b2 -= b2_offset;
  b1_dim2 = *ldb;
  b1_offset = 1 + 2 * (1 + b1_dim2 * 1);
  b1 -= b1_offset;

    /* Function Body */
  i__1 = *b;
  for (i__ = 1; i__ <= i__1; ++i__) {
    r2 = a2[(i__ << 1) + 1];
    i2 = a2[(i__ << 1) + 2];
    s2 = r2 * *c2 - i2 * *d2;
    j2 = i2 * *c2 + r2 * *d2;
    r1 = a1[(i__ << 1) + 1];
    i1 = a1[(i__ << 1) + 2];
    b1[((i__ * b1_dim2 + 1) << 1) + 1] = r1 + s2;
    b1[((i__ * b1_dim2 + 1) << 1) + 2] = i1 + j2;
    b2[((i__ * b2_dim2 + 1) << 1) + 1] = r1 - s2;
    b2[((i__ * b2_dim2 + 1) << 1) + 2] = i1 - j2;
  }
  return 0;
} /* radix2f_ */

int radix3f_(double *a1, double *a2, double *a3, double *b1, double *b2, double *b3, int *b, int *ldb, double *c2, double *d2, double *c3, double *d3)
{
  /* Initialized data */

  double sin60 = .86602540378443864676;
  double half = .5;

  /* System generated locals */
  int b1_dim2, b1_offset, b2_dim2, b2_offset, b3_dim2, b3_offset, i__1;

  /* Local variables */
  int i__;
  double i1, i2, i3, j2, j3, r1, r2, r3, t1, t2, t3, t4, s2, s3;

    /* Parameter adjustments */
  a3 -= 3;
  a2 -= 3;
  a1 -= 3;
  b3_dim2 = *ldb;
  b3_offset = 1 + 2 * (1 + b3_dim2 * 1);
  b3 -= b3_offset;
  b2_dim2 = *ldb;
  b2_offset = 1 + 2 * (1 + b2_dim2 * 1);
  b2 -= b2_offset;
  b1_dim2 = *ldb;
  b1_offset = 1 + 2 * (1 + b1_dim2 * 1);
  b1 -= b1_offset;

    /* Function Body */
  i__1 = *b;
  for (i__ = 1; i__ <= i__1; ++i__) {
    r3 = a2[(i__ << 1) + 1];
    i3 = a2[(i__ << 1) + 2];
    s3 = r3 * *c2 - i3 * *d2;
    j3 = i3 * *c2 + r3 * *d2;
    r2 = a3[(i__ << 1) + 1];
    i2 = a3[(i__ << 1) + 2];
    s2 = r2 * *c3 - i2 * *d3;
    j2 = i2 * *c3 + r2 * *d3;
    t1 = s2 + s3;
    t2 = j2 + j3;
    t3 = s2 - s3;
    t4 = j2 - j3;
    r1 = a1[(i__ << 1) + 1];
    i1 = a1[(i__ << 1) + 2];
    b1[((i__ * b1_dim2 + 1) << 1) + 1] = r1 + t1;
    b1[((i__ * b1_dim2 + 1) << 1) + 2] = i1 + t2;
    b2[((i__ * b2_dim2 + 1) << 1) + 1] = r1 - half * t1 + sin60 * t4;
    b2[((i__ * b2_dim2 + 1) << 1) + 2] = i1 - half * t2 - sin60 * t3;
    b3[((i__ * b3_dim2 + 1) << 1) + 1] = r1 - half * t1 - sin60 * t4;
    b3[((i__ * b3_dim2 + 1) << 1) + 2] = i1 - half * t2 + sin60 * t3;
  }
  return 0;
} /* radix3f_ */

int radix4f_(double *a1, double *a2, double *a3, double *a4, double *b1, double *b2, double *b3, double *b4, int *b, int *ldb, double *c2, double *d2, double *c3, double *d3, double *c4, double *d4)
{
  /* System generated locals */
  int b1_dim2, b1_offset, b2_dim2, b2_offset, b3_dim2, b3_offset, 
    b4_dim2, b4_offset, i__1;

  /* Local variables */
  int i__;
  double i0, i1, r0, r1, r2, r3, i2, i3, s1, j1, s2, j2, s3, j3, 
    t1, t2, t3, t4, t5, t6, t7, t8;

  /* Parameter adjustments */
  a4 -= 3;
  a3 -= 3;
  a2 -= 3;
  a1 -= 3;
  b4_dim2 = *ldb;
  b4_offset = 1 + 2 * (1 + b4_dim2 * 1);
  b4 -= b4_offset;
  b3_dim2 = *ldb;
  b3_offset = 1 + 2 * (1 + b3_dim2 * 1);
  b3 -= b3_offset;
  b2_dim2 = *ldb;
  b2_offset = 1 + 2 * (1 + b2_dim2 * 1);
  b2 -= b2_offset;
  b1_dim2 = *ldb;
  b1_offset = 1 + 2 * (1 + b1_dim2 * 1);
  b1 -= b1_offset;

    /* Function Body */
  i__1 = *b;
  for (i__ = 1; i__ <= i__1; ++i__) {
    r1 = a2[(i__ << 1) + 1];
    i1 = a2[(i__ << 1) + 2];
    s1 = r1 * *c2 - i1 * *d2;
    j1 = i1 * *c2 + r1 * *d2;
    r2 = a3[(i__ << 1) + 1];
    i2 = a3[(i__ << 1) + 2];
    s2 = r2 * *c3 - i2 * *d3;
    j2 = i2 * *c3 + r2 * *d3;
    r3 = a4[(i__ << 1) + 1];
    i3 = a4[(i__ << 1) + 2];
    s3 = r3 * *c4 - i3 * *d4;
    j3 = i3 * *c4 + r3 * *d4;
    r0 = a1[(i__ << 1) + 1];
    i0 = a1[(i__ << 1) + 2];
    t1 = r0 + s2;
    t2 = s1 + s3;
    t3 = i0 + j2;
    t4 = j1 + j3;
    t5 = r0 - s2;
    t6 = j1 - j3;
    t7 = i0 - j2;
    t8 = s1 - s3;
    b1[((i__ * b1_dim2 + 1) << 1) + 1] = t1 + t2;
    b1[((i__ * b1_dim2 + 1) << 1) + 2] = t3 + t4;
    b2[((i__ * b2_dim2 + 1) << 1) + 1] = t5 - t6;
    b2[((i__ * b2_dim2 + 1) << 1) + 2] = t7 + t8;
    b3[((i__ * b3_dim2 + 1) << 1) + 1] = t1 - t2;
    b3[((i__ * b3_dim2 + 1) << 1) + 2] = t3 - t4;
    b4[((i__ * b4_dim2 + 1) << 1) + 1] = t5 + t6;
    b4[((i__ * b4_dim2 + 1) << 1) + 2] = t7 - t8;
  }
  return 0;
} /* radix4f_ */

int radix5f_(double *a1,double * a2,double * a3,double * a4,double * a5,double * b1,double * b2,double * b3,double * b4,double * b5,int *b, int *ldb, double *c2,double * d2,double * c3,double * d3,double * c4,double * d4,double * c5,double * d5)
{
  /* Initialized data */

  double cos36 = .8090169943749474241;
  double sin36 = .58778525229247312917;
  double cos72 = .3090169943749474241;
  double sin72 = .95105651629515357212;

  /* System generated locals */
  int b1_dim2, b1_offset, b2_dim2, b2_offset, b3_dim2, b3_offset, 
    b4_dim2, b4_offset, b5_dim2, b5_offset, i__1;

    /* Local variables */
  int i__;
  double r1, r2, r3, r4, r5, i1, i2, i3, i4, i5, s2, s3, s4, s5, 
    j2, j3, j4, j5, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12,
    t13, t14, t15, t16;

    /* Parameter adjustments */
  a5 -= 3;
  a4 -= 3;
  a3 -= 3;
  a2 -= 3;
  a1 -= 3;
  b5_dim2 = *ldb;
  b5_offset = 1 + 2 * (1 + b5_dim2 * 1);
  b5 -= b5_offset;
  b4_dim2 = *ldb;
  b4_offset = 1 + 2 * (1 + b4_dim2 * 1);
  b4 -= b4_offset;
  b3_dim2 = *ldb;
  b3_offset = 1 + 2 * (1 + b3_dim2 * 1);
  b3 -= b3_offset;
  b2_dim2 = *ldb;
  b2_offset = 1 + 2 * (1 + b2_dim2 * 1);
  b2 -= b2_offset;
  b1_dim2 = *ldb;
  b1_offset = 1 + 2 * (1 + b1_dim2 * 1);
  b1 -= b1_offset;

    /* Function Body */
  i__1 = *b;
  for (i__ = 1; i__ <= i__1; ++i__) {
    r5 = a2[(i__ << 1) + 1];
    i5 = a2[(i__ << 1) + 2];
    s5 = r5 * *c2 - i5 * *d2;
    j5 = i5 * *c2 + r5 * *d2;
    r4 = a3[(i__ << 1) + 1];
    i4 = a3[(i__ << 1) + 2];
    s4 = r4 * *c3 - i4 * *d3;
    j4 = i4 * *c3 + r4 * *d3;
    r3 = a4[(i__ << 1) + 1];
    i3 = a4[(i__ << 1) + 2];
    s3 = r3 * *c4 - i3 * *d4;
    j3 = i3 * *c4 + r3 * *d4;
    r2 = a5[(i__ << 1) + 1];
    i2 = a5[(i__ << 1) + 2];
    s2 = r2 * *c5 - i2 * *d5;
    j2 = i2 * *c5 + r2 * *d5;
    t1 = s2 + s5;
    t2 = s3 + s4;
    t3 = s2 - s5;
    t4 = s3 - s4;
    t5 = j2 + j5;
    t6 = j3 + j4;
    t7 = j2 - j5;
    t8 = j3 - j4;
    r1 = a1[(i__ << 1) + 1];
    i1 = a1[(i__ << 1) + 2];
    t9 = t1 * cos72 - t2 * cos36;
    t10 = t5 * cos72 - t6 * cos36;
    t11 = t2 * cos72 - t1 * cos36;
    t12 = t6 * cos72 - t5 * cos36;
    t13 = t7 * sin72 + t8 * sin36;
    t14 = t3 * sin72 + t4 * sin36;
    t15 = t8 * sin72 - t7 * sin36;
    t16 = t4 * sin72 - t3 * sin36;
    b1[((i__ * b1_dim2 + 1) << 1) + 1] = r1 + t1 + t2;
    b1[((i__ * b1_dim2 + 1) << 1) + 2] = i1 + t5 + t6;
    b2[((i__ * b2_dim2 + 1) << 1) + 1] = r1 + t9 + t13;
    b2[((i__ * b2_dim2 + 1) << 1) + 2] = i1 + t10 - t14;
    b3[((i__ * b3_dim2 + 1) << 1) + 1] = r1 + t11 - t15;
    b3[((i__ * b3_dim2 + 1) << 1) + 2] = i1 + t12 + t16;
    b4[((i__ * b4_dim2 + 1) << 1) + 1] = r1 + t11 + t15;
    b4[((i__ * b4_dim2 + 1) << 1) + 2] = i1 + t12 - t16;
    b5[((i__ * b5_dim2 + 1) << 1) + 1] = r1 + t9 - t13;
    b5[((i__ * b5_dim2 + 1) << 1) + 2] = i1 + t10 + t14;
  }
  return 0;
} /* radix5f_ */

int radix7f_(double *a1,double * a2,double * a3, double *a4,double * a5,double * a6,double * a7,double * b1,double * b2,double * b3,double * b4,double * b5, double *b6,double * b7,int * b,int * ldb)
{
  /* Initialized data */

  double c1 = .623489801858733530525;
  double s1 = .78183148246802980870844;
  double c2 = .2225209339563144042889;
  double s2 = .97492791218182360701813;
  double c3 = .9009688679024191262361;
  double s3 = .43388373911755812047577;

  /* System generated locals */
  int b1_dim2, b1_offset, b2_dim2, b2_offset, b3_dim2, b3_offset, 
    b4_dim2, b4_offset, b5_dim2, b5_offset, b6_dim2, b6_offset, 
    b7_dim2, b7_offset, i__1;

  /* Local variables */
  int i__;
  double d1, p1, p2, p3, p4, t1, t2, t3, t4, t5, t6, t7, p5, p6, 
    p7, r1, r2, d2, r3, d3, i1, j1, i2, j2, i3, j3, t11, t21, t31, 
    t12, t22, t32, t13, t23, t33, t14, t24, t34;

    /* Parameter adjustments */
  a7 -= 3;
  a6 -= 3;
  a5 -= 3;
  a4 -= 3;
  a3 -= 3;
  a2 -= 3;
  a1 -= 3;
  b7_dim2 = *ldb;
  b7_offset = 1 + 2 * (1 + b7_dim2 * 1);
  b7 -= b7_offset;
  b6_dim2 = *ldb;
  b6_offset = 1 + 2 * (1 + b6_dim2 * 1);
  b6 -= b6_offset;
  b5_dim2 = *ldb;
  b5_offset = 1 + 2 * (1 + b5_dim2 * 1);
  b5 -= b5_offset;
  b4_dim2 = *ldb;
  b4_offset = 1 + 2 * (1 + b4_dim2 * 1);
  b4 -= b4_offset;
  b3_dim2 = *ldb;
  b3_offset = 1 + 2 * (1 + b3_dim2 * 1);
  b3 -= b3_offset;
  b2_dim2 = *ldb;
  b2_offset = 1 + 2 * (1 + b2_dim2 * 1);
  b2 -= b2_offset;
  b1_dim2 = *ldb;
  b1_offset = 1 + 2 * (1 + b1_dim2 * 1);
  b1 -= b1_offset;

    /* Function Body */
  i__1 = *b;
  for (i__ = 1; i__ <= i__1; ++i__) {
    t1 = a2[(i__ << 1) + 1];
    t2 = a7[(i__ << 1) + 1];
    r1 = t1 + t2;
    d1 = t1 - t2;
    t11 = r1 * c1;
    t21 = r1 * c2;
    t31 = r1 * c3;
    t12 = d1 * s1;
    t22 = d1 * s2;
    t32 = d1 * s3;
    t3 = a3[(i__ << 1) + 1];
    t4 = a6[(i__ << 1) + 1];
    r2 = t3 + t4;
    d2 = t3 - t4;
    t11 -= r2 * c2;
    t21 += r2 * c3;
    t31 -= r2 * c1;
    t12 += d2 * s2;
    t22 -= d2 * s3;
    t32 -= d2 * s1;
    t5 = a4[(i__ << 1) + 1];
    t6 = a5[(i__ << 1) + 1];
    r3 = t5 + t6;
    d3 = t5 - t6;
    t7 = a1[(i__ << 1) + 1];
    b1[((i__ * b1_dim2 + 1) << 1) + 1] = t7 + r1 + r2 + r3;
    t11 -= r3 * c3;
    t21 -= r3 * c1;
    t31 += r3 * c2;
    t12 += d3 * s3;
    t22 -= d3 * s1;
    t32 += d3 * s2;
    p1 = a2[(i__ << 1) + 2];
    p2 = a7[(i__ << 1) + 2];
    i1 = p1 + p2;
    j1 = p1 - p2;
    t13 = i1 * c1;
    t23 = i1 * c2;
    t33 = i1 * c3;
    t14 = j1 * s1;
    t24 = j1 * s2;
    t34 = j1 * s3;
    p3 = a3[(i__ << 1) + 2];
    p4 = a6[(i__ << 1) + 2];
    i2 = p3 + p4;
    j2 = p3 - p4;
    t13 -= i2 * c2;
    t23 += i2 * c3;
    t33 -= i2 * c1;
    t14 += j2 * s2;
    t24 -= j2 * s3;
    t34 -= j2 * s1;
    p5 = a4[(i__ << 1) + 2];
    p6 = a5[(i__ << 1) + 2];
    i3 = p5 + p6;
    j3 = p5 - p6;
    p7 = a1[(i__ << 1) + 2];
    b1[((i__ * b1_dim2 + 1) << 1) + 2] = p7 + i1 + i2 + i3;
    t13 -= i3 * c3;
    t23 -= i3 * c1;
    t33 += i3 * c2;
    t14 += j3 * s3;
    t24 -= j3 * s1;
    t34 += j3 * s2;
    t11 = t7 + t11;
    t21 = t7 - t21;
    t31 = t7 - t31;
    t13 = p7 + t13;
    t23 = p7 - t23;
    t33 = p7 - t33;
    b2[((i__ * b2_dim2 + 1) << 1) + 1] = t11 - t14;
    b2[((i__ * b2_dim2 + 1) << 1) + 2] = t13 + t12;
    b7[((i__ * b7_dim2 + 1) << 1) + 1] = t11 + t14;
    b7[((i__ * b7_dim2 + 1) << 1) + 2] = t13 - t12;
    b3[((i__ * b3_dim2 + 1) << 1) + 1] = t21 - t24;
    b3[((i__ * b3_dim2 + 1) << 1) + 2] = t23 + t22;
    b6[((i__ * b6_dim2 + 1) << 1) + 1] = t21 + t24;
    b6[((i__ * b6_dim2 + 1) << 1) + 2] = t23 - t22;
    b4[((i__ * b4_dim2 + 1) << 1) + 1] = t31 - t34;
    b4[((i__ * b4_dim2 + 1) << 1) + 2] = t33 + t32;
    b5[((i__ * b5_dim2 + 1) << 1) + 1] = t31 + t34;
    b5[((i__ * b5_dim2 + 1) << 1) + 2] = t33 - t32;
  }
  return 0;
} /* radix7f_ */

int radix2b_(double *a1, double *a2, double *b1, double *b2,int * b, int *ldb, double *c2, double *d2)
{
  /* System generated locals */
  int b1_dim2, b1_offset, b2_dim2, b2_offset, i__1;

    /* Local variables */
  int i__;
  double i1, i2, j2, r1, r2, s2;
  /* Parameter adjustments */
  a2 -= 3;
  a1 -= 3;
  b2_dim2 = *ldb;
  b2_offset = 1 + 2 * (1 + b2_dim2 * 1);
  b2 -= b2_offset;
  b1_dim2 = *ldb;
  b1_offset = 1 + 2 * (1 + b1_dim2 * 1);
  b1 -= b1_offset;

  /* Function Body */
  i__1 = *b;
  for (i__ = 1; i__ <= i__1; ++i__) {
    r2 = a2[(i__ << 1) + 1];
    i2 = a2[(i__ << 1) + 2];
    s2 = r2 * *c2 + i2 * *d2;
    j2 = i2 * *c2 - r2 * *d2;
    r1 = a1[(i__ << 1) + 1];
    i1 = a1[(i__ << 1) + 2];
    b1[((i__ * b1_dim2 + 1) << 1) + 1] = r1 + s2;
    b1[((i__ * b1_dim2 + 1) << 1) + 2] = i1 + j2;
    b2[((i__ * b2_dim2 + 1) << 1) + 1] = r1 - s2;
    b2[((i__ * b2_dim2 + 1) << 1) + 2] = i1 - j2;
  }     
  return 0;
} /* radix2b_ */

int radix3b_(double *a1,double * a2, double *a3,double * b1,double * b2,double * b3, int *b, int *ldb, double *c2, double *d2, double *c3, double *d3)
{
  /* Initialized data */

  double sin60 = .86602540378443864676;
  double half = .5;

  /* System generated locals */
  int b1_dim2, b1_offset, b2_dim2, b2_offset, b3_dim2, b3_offset, i__1;

  /* Local variables */
  int i__;
  double i1, i2, i3, j2, j3, r1, r2, r3, t1, t2, t3, t4, s2, s3;

    /* Parameter adjustments */
  a3 -= 3;
  a2 -= 3;
  a1 -= 3;
  b3_dim2 = *ldb;
  b3_offset = 1 + 2 * (1 + b3_dim2 * 1);
  b3 -= b3_offset;
  b2_dim2 = *ldb;
  b2_offset = 1 + 2 * (1 + b2_dim2 * 1);
  b2 -= b2_offset;
  b1_dim2 = *ldb;
  b1_offset = 1 + 2 * (1 + b1_dim2 * 1);
  b1 -= b1_offset;

    /* Function Body */
  i__1 = *b;
  for (i__ = 1; i__ <= i__1; ++i__) {
    r2 = a2[(i__ << 1) + 1];
    i2 = a2[(i__ << 1) + 2];
    s2 = r2 * *c2 + i2 * *d2;
    j2 = i2 * *c2 - r2 * *d2;
    r3 = a3[(i__ << 1) + 1];
    i3 = a3[(i__ << 1) + 2];
    s3 = r3 * *c3 + i3 * *d3;
    j3 = i3 * *c3 - r3 * *d3;
    t1 = s2 + s3;
    t2 = j2 + j3;
    t3 = s2 - s3;
    t4 = j2 - j3;
    r1 = a1[(i__ << 1) + 1];
    i1 = a1[(i__ << 1) + 2];
    b1[((i__ * b1_dim2 + 1) << 1) + 1] = r1 + t1;
    b1[((i__ * b1_dim2 + 1) << 1) + 2] = i1 + t2;
    b2[((i__ * b2_dim2 + 1) << 1) + 1] = r1 - half * t1 + sin60 * t4;
    b2[((i__ * b2_dim2 + 1) << 1) + 2] = i1 - half * t2 - sin60 * t3;
    b3[((i__ * b3_dim2 + 1) << 1) + 1] = r1 - half * t1 - sin60 * t4;
    b3[((i__ * b3_dim2 + 1) << 1) + 2] = i1 - half * t2 + sin60 * t3;
  }
  return 0;
} /* radix3b_ */

int radix4b_(double *a1, double *a2, double *a3, double *a4, double *b1, double *b2, double *b3, double *b4, int *b, int *ldb, double *c2, double *d2, double *c3, double *d3, double *c4, double *d4)
{
  /* System generated locals */
  int b1_dim2, b1_offset, b2_dim2, b2_offset, b3_dim2, b3_offset, 
    b4_dim2, b4_offset, i__1;

  /* Local variables */
  int i__;
  double i0, i1, r0, r1, r2, r3, i2, i3, s1, j1, s2, j2, s3, j3, 
    t1, t2, t3, t4, t5, t6, t7, t8;

  /* Parameter adjustments */
  a4 -= 3;
  a3 -= 3;
  a2 -= 3;
  a1 -= 3;
  b4_dim2 = *ldb;
  b4_offset = 1 + 2 * (1 + b4_dim2 * 1);
  b4 -= b4_offset;
  b3_dim2 = *ldb;
  b3_offset = 1 + 2 * (1 + b3_dim2 * 1);
  b3 -= b3_offset;
  b2_dim2 = *ldb;
  b2_offset = 1 + 2 * (1 + b2_dim2 * 1);
  b2 -= b2_offset;
  b1_dim2 = *ldb;
  b1_offset = 1 + 2 * (1 + b1_dim2 * 1);
  b1 -= b1_offset;

    /* Function Body */
  i__1 = *b;
  for (i__ = 1; i__ <= i__1; ++i__) {
    r3 = a2[(i__ << 1) + 1];
    i3 = a2[(i__ << 1) + 2];
    s3 = r3 * *c2 + i3 * *d2;
    j3 = i3 * *c2 - r3 * *d2;
    r2 = a3[(i__ << 1) + 1];
    i2 = a3[(i__ << 1) + 2];
    s2 = r2 * *c3 + i2 * *d3;
    j2 = i2 * *c3 - r2 * *d3;
    r1 = a4[(i__ << 1) + 1];
    i1 = a4[(i__ << 1) + 2];
    s1 = r1 * *c4 + i1 * *d4;
    j1 = i1 * *c4 - r1 * *d4;
    r0 = a1[(i__ << 1) + 1];
    i0 = a1[(i__ << 1) + 2];
    t1 = r0 + s2;
    t2 = s1 + s3;
    t3 = i0 + j2;
    t4 = j1 + j3;
    t5 = r0 - s2;
    t6 = j1 - j3;
    t7 = i0 - j2;
    t8 = s1 - s3;
    b1[((i__ * b1_dim2 + 1) << 1) + 1] = t1 + t2;
    b1[((i__ * b1_dim2 + 1) << 1) + 2] = t3 + t4;
    b2[((i__ * b2_dim2 + 1) << 1) + 1] = t5 - t6;
    b2[((i__ * b2_dim2 + 1) << 1) + 2] = t7 + t8;
    b3[((i__ * b3_dim2 + 1) << 1) + 1] = t1 - t2;
    b3[((i__ * b3_dim2 + 1) << 1) + 2] = t3 - t4;
    b4[((i__ * b4_dim2 + 1) << 1) + 1] = t5 + t6;
    b4[((i__ * b4_dim2 + 1) << 1) + 2] = t7 - t8;
  }
  return 0;
} /* radix4b_ */

int radix5b_(double *a1,double * a2,double * a3,double * a4,double * a5,double * b1,double * b2,double * b3,double * b4,double * b5,int * b,int *ldb, double *c2,double * d2,double * c3,double * d3,double * c4,double * d4, double *c5, double *d5)
{
  /* Initialized data */

  double cos36 = .8090169943749474241;
  double sin36 = .58778525229247312917;
  double cos72 = .3090169943749474241;
  double sin72 = .95105651629515357212;

  /* System generated locals */
  int b1_dim2, b1_offset, b2_dim2, b2_offset, b3_dim2, b3_offset, 
    b4_dim2, b4_offset, b5_dim2, b5_offset, i__1;

    /* Local variables */
  int i__;
  double r1, r2, r3, r4, r5, i1, i2, i3, i4, i5, s2, s3, s4, s5, 
    j2, j3, j4, j5, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12,
    t13, t14, t15, t16;

    /* Parameter adjustments */
  a5 -= 3;
  a4 -= 3;
  a3 -= 3;
  a2 -= 3;
  a1 -= 3;
  b5_dim2 = *ldb;
  b5_offset = 1 + 2 * (1 + b5_dim2 * 1);
  b5 -= b5_offset;
  b4_dim2 = *ldb;
  b4_offset = 1 + 2 * (1 + b4_dim2 * 1);
  b4 -= b4_offset;
  b3_dim2 = *ldb;
  b3_offset = 1 + 2 * (1 + b3_dim2 * 1);
  b3 -= b3_offset;
  b2_dim2 = *ldb;
  b2_offset = 1 + 2 * (1 + b2_dim2 * 1);
  b2 -= b2_offset;
  b1_dim2 = *ldb;
  b1_offset = 1 + 2 * (1 + b1_dim2 * 1);
  b1 -= b1_offset;

    /* Function Body */
  i__1 = *b;
  for (i__ = 1; i__ <= i__1; ++i__) {
    r2 = a2[(i__ << 1) + 1];
    i2 = a2[(i__ << 1) + 2];
    s2 = r2 * *c2 + i2 * *d2;
    j2 = i2 * *c2 - r2 * *d2;
    r3 = a3[(i__ << 1) + 1];
    i3 = a3[(i__ << 1) + 2];
    s3 = r3 * *c3 + i3 * *d3;
    j3 = i3 * *c3 - r3 * *d3;
    r4 = a4[(i__ << 1) + 1];
    i4 = a4[(i__ << 1) + 2];
    s4 = r4 * *c4 + i4 * *d4;
    j4 = i4 * *c4 - r4 * *d4;
    r5 = a5[(i__ << 1) + 1];
    i5 = a5[(i__ << 1) + 2];
    s5 = r5 * *c5 + i5 * *d5;
    j5 = i5 * *c5 - r5 * *d5;
    t1 = s2 + s5;
    t2 = s3 + s4;
    t3 = s2 - s5;
    t4 = s3 - s4;
    t5 = j2 + j5;
    t6 = j3 + j4;
    t7 = j2 - j5;
    t8 = j3 - j4;
    r1 = a1[(i__ << 1) + 1];
    i1 = a1[(i__ << 1) + 2];
    t9 = t1 * cos72 - t2 * cos36;
    t10 = t5 * cos72 - t6 * cos36;
    t11 = t2 * cos72 - t1 * cos36;
    t12 = t6 * cos72 - t5 * cos36;
    t13 = t7 * sin72 + t8 * sin36;
    t14 = t3 * sin72 + t4 * sin36;
    t15 = t8 * sin72 - t7 * sin36;
    t16 = t4 * sin72 - t3 * sin36;
    b1[((i__ * b1_dim2 + 1) << 1) + 1] = r1 + t1 + t2;
    b1[((i__ * b1_dim2 + 1) << 1) + 2] = i1 + t5 + t6;
    b2[((i__ * b2_dim2 + 1) << 1) + 1] = r1 + t9 + t13;
    b2[((i__ * b2_dim2 + 1) << 1) + 2] = i1 + t10 - t14;
    b3[((i__ * b3_dim2 + 1) << 1) + 1] = r1 + t11 - t15;
    b3[((i__ * b3_dim2 + 1) << 1) + 2] = i1 + t12 + t16;
    b4[((i__ * b4_dim2 + 1) << 1) + 1] = r1 + t11 + t15;
    b4[((i__ * b4_dim2 + 1) << 1) + 2] = i1 + t12 - t16;
    b5[((i__ * b5_dim2 + 1) << 1) + 1] = r1 + t9 - t13;
    b5[((i__ * b5_dim2 + 1) << 1) + 2] = i1 + t10 + t14;
  }
  return 0;
} /* radix5b_ */

int radix7b_(double *a1, double *a2,double * a3,double * a4,double * a5,double * a6,double * a7,double * b1,double * b2,double * b3,double * b4,double * b5, double *b6, double *b7, int *b, int *ldb)
{
  /* Initialized data */

  double c1 = .623489801858733530525;
  double s1 = .78183148246802980870844;
  double c2 = .2225209339563144042889;
  double s2 = .97492791218182360701813;
  double c3 = .9009688679024191262361;
  double s3 = .43388373911755812047577;

  /* System generated locals */
  int b1_dim2, b1_offset, b2_dim2, b2_offset, b3_dim2, b3_offset, 
    b4_dim2, b4_offset, b5_dim2, b5_offset, b6_dim2, b6_offset, 
    b7_dim2, b7_offset, i__1;

  /* Local variables */
  int i__;
  double d1, p1, p2, p3, p4, t1, t2, t3, t4, t5, t6, t7, p5, p6, 
    p7, r1, r2, d2, r3, d3, i1, j1, i2, j2, i3, j3, t11, t21, t31, 
    t12, t22, t32, t13, t23, t33, t14, t24, t34;

    /* Parameter adjustments */
  a7 -= 3;
  a6 -= 3;
  a5 -= 3;
  a4 -= 3;
  a3 -= 3;
  a2 -= 3;
  a1 -= 3;
  b7_dim2 = *ldb;
  b7_offset = 1 + 2 * (1 + b7_dim2 * 1);
  b7 -= b7_offset;
  b6_dim2 = *ldb;
  b6_offset = 1 + 2 * (1 + b6_dim2 * 1);
  b6 -= b6_offset;
  b5_dim2 = *ldb;
  b5_offset = 1 + 2 * (1 + b5_dim2 * 1);
  b5 -= b5_offset;
  b4_dim2 = *ldb;
  b4_offset = 1 + 2 * (1 + b4_dim2 * 1);
  b4 -= b4_offset;
  b3_dim2 = *ldb;
  b3_offset = 1 + 2 * (1 + b3_dim2 * 1);
  b3 -= b3_offset;
  b2_dim2 = *ldb;
  b2_offset = 1 + 2 * (1 + b2_dim2 * 1);
  b2 -= b2_offset;
  b1_dim2 = *ldb;
  b1_offset = 1 + 2 * (1 + b1_dim2 * 1);
  b1 -= b1_offset;

    /* Function Body */
  i__1 = *b;
  for (i__ = 1; i__ <= i__1; ++i__) {
    t1 = a7[(i__ << 1) + 1];
    t2 = a2[(i__ << 1) + 1];
    r1 = t1 + t2;
    d1 = t1 - t2;
    t11 = r1 * c1;
    t21 = r1 * c2;
    t31 = r1 * c3;
    t12 = d1 * s1;
    t22 = d1 * s2;
    t32 = d1 * s3;
    t3 = a6[(i__ << 1) + 1];
    t4 = a3[(i__ << 1) + 1];
    r2 = t3 + t4;
    d2 = t3 - t4;
    t11 -= r2 * c2;
    t21 += r2 * c3;
    t31 -= r2 * c1;
    t12 += d2 * s2;
    t22 -= d2 * s3;
    t32 -= d2 * s1;
    t5 = a5[(i__ << 1) + 1];
    t6 = a4[(i__ << 1) + 1];
    r3 = t5 + t6;
    d3 = t5 - t6;
    t7 = a1[(i__ << 1) + 1];
    b1[((i__ * b1_dim2 + 1) << 1) + 1] = t7 + r1 + r2 + r3;
    t11 -= r3 * c3;
    t21 -= r3 * c1;
    t31 += r3 * c2;
    t12 += d3 * s3;
    t22 -= d3 * s1;
    t32 += d3 * s2;
    p1 = a7[(i__ << 1) + 2];
    p2 = a2[(i__ << 1) + 2];
    i1 = p1 + p2;
    j1 = p1 - p2;
    t13 = i1 * c1;
    t23 = i1 * c2;
    t33 = i1 * c3;
    t14 = j1 * s1;
    t24 = j1 * s2;
    t34 = j1 * s3;
    p3 = a6[(i__ << 1) + 2];
    p4 = a3[(i__ << 1) + 2];
    i2 = p3 + p4;
    j2 = p3 - p4;
    t13 -= i2 * c2;
    t23 += i2 * c3;
    t33 -= i2 * c1;
    t14 += j2 * s2;
    t24 -= j2 * s3;
    t34 -= j2 * s1;
    p5 = a5[(i__ << 1) + 2];
    p6 = a4[(i__ << 1) + 2];
    i3 = p5 + p6;
    j3 = p5 - p6;
    p7 = a1[(i__ << 1) + 2];
    b1[((i__ * b1_dim2 + 1) << 1) + 2] = p7 + i1 + i2 + i3;
    t13 -= i3 * c3;
    t23 -= i3 * c1;
    t33 += i3 * c2;
    t14 += j3 * s3;
    t24 -= j3 * s1;
    t34 += j3 * s2;
    t11 = t7 + t11;
    t21 = t7 - t21;
    t31 = t7 - t31;
    t13 = p7 + t13;
    t23 = p7 - t23;
    t33 = p7 - t33;
    b2[((i__ * b2_dim2 + 1) << 1) + 1] = t11 - t14;
    b2[((i__ * b2_dim2 + 1) << 1) + 2] = t13 + t12;
    b7[((i__ * b7_dim2 + 1) << 1) + 1] = t11 + t14;
    b7[((i__ * b7_dim2 + 1) << 1) + 2] = t13 - t12;
    b3[((i__ * b3_dim2 + 1) << 1) + 1] = t21 - t24;
    b3[((i__ * b3_dim2 + 1) << 1) + 2] = t23 + t22;
    b6[((i__ * b6_dim2 + 1) << 1) + 1] = t21 + t24;
    b6[((i__ * b6_dim2 + 1) << 1) + 2] = t23 - t22;
    b4[((i__ * b4_dim2 + 1) << 1) + 1] = t31 - t34;
    b4[((i__ * b4_dim2 + 1) << 1) + 2] = t33 + t32;
    b5[((i__ * b5_dim2 + 1) << 1) + 1] = t31 + t34;
    b5[((i__ * b5_dim2 + 1) << 1) + 2] = t33 - t32;
  }
  return 0;
} /* radix7b_ */

int gcd_(int* n1,int* n2)
{
  /* System generated locals */
  int ret_val;

    /* Local variables */
  int i__, j, n;

  i__ = *n1;
  j = *n2;
 L10:
  n = min(i__,j);
  ret_val = max(i__,j);
  if (n == 0) {
    return ret_val;
  }
  i__ = n;
  j = ret_val - n;
  goto L10;
} /* gcd_ */

void fftfact_(int *n,int *nf,int *f)
{

  /* Local variables */
  int i__, n7;

  /* Parameter adjustments */
  --f;

  /* Function Body */
  i__ = *n;
  *nf = 0;
  n7 = 0;
 L10:
  if (i__ == 1) {
    if (n7 > 1) {
      printf("Only one factor of 7 allowed in fft dimension.\n");
      exit(-1);
    }
    return;
  }
  if (i__ % 7 == 0) {
    ++(*nf);
    f[*nf] = 7;
    i__ /= 7;
    ++n7;
    goto L10;
  }
  if (i__ % 5 == 0) {
    ++(*nf);
    f[*nf] = 5;
    i__ /= 5;
    goto L10;
  }
  if (i__ % 4 == 0) {
    ++(*nf);
    f[*nf] = 4;
    i__ /= 4;
    goto L10;
  }
  if (i__ % 3 == 0) {
    ++(*nf);
    f[*nf] = 3;
    i__ /= 3;
    goto L10;
  }
  if (i__ % 2 == 0) {
    ++(*nf);
    f[*nf] = 2;
    i__ /= 2;
    goto L10;
  }
  printf("Illegal fft size %i\n",(*n));
  exit(-1);
} /* fftfact_ */

void fftcoeff_(int *n,zomplex * coeff)
{
  /* System generated locals */
  int i__1, i__2;
  double d__1, d__2;
  zomplex z__1;


    /* Local variables */
  int i__;
  double angle;

    /* Parameter adjustments */
  --coeff;

  /* Function Body */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    angle = (double) (i__ - 1) / (double) (*n) * 2. * 
      3.141592653589793238462643;
    i__2 = i__;
    d__1 = cos(angle);
    d__2 = sin(angle);
    z__1.re = d__1, z__1.im = d__2;
    coeff[i__2].re = z__1.re, coeff[i__2].im = z__1.im;
    /* L10: */
  }
} /* fftcoeff_ */

void zfftcopyin_(int *n,int * b,int * bdim,zomplex * a,int *  inc ,int * lda,int * rindex,zomplex * w,int * lpref)
{
  /* System generated locals */
  int w_dim1, w_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
  int i__, j, ia;

  /*    if(lpref){ */
  /*      if(*lda <= *inc) */
  /*        prefetch(a,lda,inc,b,n); */
  /*      else */
  /*        prefetch(a,inc,lda,n,b); */
  /*    } */
  /* Copy data */
  /* Parameter adjustments */
  --rindex;
  w_dim1 = *bdim;
  w_offset = 1 + w_dim1 * 1;
  w -= w_offset;
  --a;

  /* Function Body */
  i__1 = *n;
  for (j = 1; j <= i__1; ++j) {
    ia = rindex[j];
    i__2 = *b;
    for (i__ = 1; i__ <= i__2; ++i__) {
      i__3 = i__ + j * w_dim1;
      i__4 = ia;
      w[i__3].re = a[i__4].re, w[i__3].im = a[i__4].im;
      ia += *lda;
    }
  }
} /* zfftcopyin_ */

int zfftrev_(int *n, int *inc, int *nf, int *f, int * nmax , int *rindex,int * ipadok)
{
  int imax[10], irev, i__, s, fi, ii[10], stride[10];

  /* OK up to 100000 elements */
  /* Parameter adjustments */
  --f;
  --rindex;

  /* Function Body */
  s = 1;
  for (fi = *nf; fi >= 1; --fi) {
    ii[fi - 1] = 1;
    stride[fi - 1] = s;
    s *= f[fi];
    imax[fi - 1] = s;
  }
  i__ = 1;
  irev = 1;
 L20:
  rindex[i__] = (irev - 1) * *inc + 1;
  if (i__ == *n) {
    goto L100;
  }
  ++i__;
  fi = 1;
 L10:
  irev += stride[fi - 1];
  ii[fi - 1] += stride[fi - 1];
  if (ii[fi - 1] > imax[fi - 1]) {
    ii[fi - 1] = 1;
    irev -= imax[fi - 1];
    ++fi;
    goto L10;
  }
  goto L20;
 L100:
  if (*n * *inc / f[*nf] % 1024 == 0) {
    *ipadok = 0;
  } else {
    *ipadok = 1;
  }
  return 0;
} /* zfftrev_ */

/* Row (inc=1) transform code on bitreversed vectors */
/* n      length of transform (sub) vector */
/* b      # of transforms in this block */
/* bdim   # of rows in w, b.le.bdim */
/* w      bdim*n matrix containing data elements */
/* nf     # of factors in n */
/* f      vector of factors */
/* coeff  exp[2 pi i k/n], k=0..n-1] */
/* a      pointer to output vector */
/* inc    from caller, stride within transform */
/* lda    from caller, stride between transforms measured in elements */
/* ipadok if zero, don't copy data back to a here because doing so would */
/*        trash the L1 cache */
/* nrem   product of remaining factors (n unless blocking is used) */
/* irow   row #, max(irow)*nrem = total length */
int zfftm1rowf_(int *n, int *b, int *bdim,zomplex * w, int *nf, int *f, double *coeff, zomplex *a, int *inc, int *lda, 
		int *ipadok, int *nprev, int *irow, int *nrem)
{
  /* System generated locals */
  int i__1, i__2, i__3, i__4;

    /* Local variables */
  int aind, ntot, i__, j, p, index, ff, ic, nt, ninner, stride;
  int astride, bstride;


    /* Parameter adjustments */
  --coeff;
  --w;
  --f;
  --a;

    /* Function Body */
  stride = *bdim;
  ninner = 1;
  bstride = *bdim;
  ntot = *n * *bdim;
  nt = *n;
  i__1 = *nf - *ipadok;
  for (p = 1; p <= i__1; ++p) {
    ff = f[p];
    bstride *= ff;
    nt /= ff;
    i__2 = ntot;
    i__3 = bstride;
    for (i__ = 1; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__3) {
      index = i__;
      ic = (*irow - 1) * nt * *nrem;
      i__4 = ninner;
      for (j = 1; j <= i__4; ++j) {
      
	if (ff == 2) {
	  radix2f_((double*)(&w[index]), (double*)(&w[index + stride]), (double*)(&w[index]), (double*)(&w[index + stride]), b, &const1, &coeff[(ic << 1) + 1], &coeff[(ic << 1) + 2]);
	} else if (ff == 3) {
	  radix3f_((double*)(&w[index]), (double*)(&w[index + stride]), (double*)(&w[index + (stride << 1)]), (double*)(&w[index]), (double*)(&w[index + stride]), (double*)(&w[index + (stride << 1)]), b, &const1, &coeff[(ic << 1) + 1], &coeff[(ic << 1) + 2], &coeff[((ic + ic) <<  1) + 1], &coeff[ic * 4 + 2]);
	} else if (ff == 4) {
	  radix4f_((double*)(&w[index]), (double*)(&w[index + stride]), (double*)(&w[index + (stride << 1)]), (double*)(&w[index + stride * 3]), (double*)(&w[index]),  (double*)(&w[index + stride]), (double*)(&w[index + (stride << 1)]), (double*)(&w[index + stride * 3]), b, &const1, &coeff[(ic << 1) + 1], &coeff[(ic << 1) + 2], &coeff[ic * 4 + 1], &coeff[ic * 4 + 2], &coeff[ic * 6 + 1], &coeff[ic * 6 + 2]);
	} else if (ff == 5) {
	  radix5f_((double*)(&w[index]), (double*)(&w[index + stride]), (double*)(&w[index + (stride << 1)]), (double*)(&w[index + stride * 3]), (double*)(&w[index + (stride << 2)]), (double*)(&w[index]), (double*)(&w[index + stride]), (double*)(&w[index + (stride << 1)]), (double*)(&w[index + stride * 3]), (double*)(&w[index + (stride << 2)]), b, &const1, &coeff[(ic << 1) + 1], &coeff[(ic << 1) + 2], &coeff[ic * 4 + 1], &coeff[ic * 4 + 2], &coeff[ic *  6 + 1], &coeff[ic * 6 + 2], &coeff[(ic << 3) + 1], &coeff[(ic << 3) + 2]);
	} else if (ff == 7) {
	  radix7f_((double*)(&w[index]), (double*)(&w[index + stride]), (double*)(&w[index + (stride << 1)]), (double*)(&w[index + stride * 3]), (double*)(&w[index + (stride << 2)]), (double*)(&w[index + stride * 5]), (double*)(&w[index + stride * 6]), (double*)(&w[index]), (double*)(&w[index + stride]), (double*)(&w[index + (stride << 1)]), (double*)(&w[index + stride * 3]), (double*)(&w[index + (stride << 2)]), (double*)(&w[index + stride * 5]),  (double*)(&w[index + stride * 6]), b, &const1);
	}
	index += *bdim;
	ic += nt * (*nprev * *nrem);
      }
    }
    stride *= ff;
    ninner *= ff;
    /* L10: */
  }
  if (*ipadok == 0) {
    return 0;
  }
  ff = f[p];
  astride = *inc * ninner;
  index = 1;
  aind = 1;
  ic = (*irow - 1) * *nrem;
  i__1 = ninner;
  for (j = 1; j <= i__1; ++j) {
    if (ff == 2) {
      radix2f_((double*)(&w[index]), (double*)(&w[index + stride]), (double*)(&a[aind]), (double*)(&a[aind + astride]), b, lda, &coeff[(ic << 1) + 1], &coeff[(ic << 1) + 2]);
    } else if (ff == 3) {
      radix3f_((double*)(&w[index]), (double*)(&w[index + stride]), (double*)(&w[index + (stride << 1)]), (double*)(&a[aind]), (double*)(&a[aind + astride]), (double*)(&a[aind + (astride << 1)]),b, lda, &coeff[(ic << 1) + 1], &coeff[(ic << 1) + 2], & coeff[ic * 4 + 1], &coeff[ic * 4 + 2]);
    } else if (ff == 4) {
      radix4f_((double*)(&w[index]), (double*)(&w[index + stride]), (double*)(&w[index + (stride << 1)]), (double*)(&w[index + stride * 3]), (double*)(&a[aind]), (double*)(&a[aind + astride]), (double*)(&a[aind + (astride << 1)]),(double*)(&a[aind + astride * 3]), b, lda, & coeff[(ic << 1) + 1], &coeff[(ic << 1) + 2], &coeff[ic * 4 + 1], &coeff[ic * 4 + 2], &coeff[ic * 6 + 1], &coeff[ic * 6 + 2]);
    } else if (ff == 5) {
      radix5f_((double*)(&w[index]), (double*)(&w[index + stride]), (double*)(&w[index + (stride << 1)]), (double*)(&w[index + stride * 3]), (double*)(&w[index + (stride << 2)]), (double*)(&a[aind]), (double*)(&a[aind + astride]), (double*)(&a[aind + (astride << 1)]),(double*)(&a[aind + astride * 3]), (double*)(&a[aind + (astride << 2)]), b, lda, & coeff[(ic << 1) + 1], &coeff[(ic << 1) + 2], &coeff[ic * 4 + 1], &coeff[ic * 4 + 2], &coeff[ic * 6 + 1], &coeff[ic * 6 + 2], &coeff[(ic << 3) + 1], &coeff[(ic << 3) + 2]);
    } else if (ff == 7) {
      radix7f_((double*)(&w[index]), (double*)(&w[index + stride]), (double*)(&w[index + (stride << 1)]), (double*)(&w[index + stride * 3]), (double*)(&w[index + (stride << 2)]), (double*)(&w[index + stride * 5]), (double*)(&w[index + stride * 6]), (double*)(&a[aind]), (double*)(&a[aind + astride]), (double*)(&a[aind + (astride << 1)]),(double*)(&a[aind + astride * 3]), (double*)(&a[aind + (astride << 2)]), (double*)(&a[aind + astride * 5]), (double*)(&a[aind + astride * 6]), b, lda);
    }
    index += *bdim;
    ic += *nrem;
    aind += *inc;
  }
  return 0;
} /* zfftm1rowf_ */

int zfftm1rowb_(int *n, int *b, int *bdim,zomplex *w, int *nf, int *f, double *coeff, zomplex *a, int *inc, int *lda, 
		int *ipadok, int *nprev, int *irow, int *nrem)
{
  /* System generated locals */
  int i__1, i__2, i__3, i__4;

    /* Local variables */
  int aind, ntot, i__, j, p, index, ff, ic, nt, ninner, stride;
  int astride, bstride;

    /* Parameter adjustments */
  --coeff;
  --w;
  --f;
  --a;

    /* Function Body */
  stride = *bdim;
  ninner = 1;
  bstride = *bdim;
  ntot = *n * *bdim;
  nt = *n;
  i__1 = *nf - *ipadok;
  for (p = 1; p <= i__1; ++p) {
    ff = f[p];
    bstride *= ff;
    nt /= ff;
    i__2 = ntot;
    i__3 = bstride;
    for (i__ = 1; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__3) {
      index = i__;
      ic = (*irow - 1) * nt * *nrem;
      i__4 = ninner;
      for (j = 1; j <= i__4; ++j) {
	if (ff == 2) {
	  radix2b_((double*)(&w[index]), (double*)(&w[index + stride]), (double*)(&w[index]), (double*)(&w[index + stride]), b, &const1, &coeff[(ic << 1) + 1],  &coeff[(ic << 1) + 2]);
	} else if (ff == 3) {
	  radix3b_((double*)(&w[index]), (double*)(&w[index + stride]), (double*)(&w[index + (stride << 1)]), (double*)(&w[index]), (double*)(&w[index + stride]), (double*)(&w[index + (stride << 1)]), b, &const1, &coeff[(ic << 1) + 1], &coeff[(ic << 1) + 2], &coeff[ic * 4 + 1], &coeff[ic * 4 + 2]);
	} else if (ff == 4) {
	  radix4b_((double*)(&w[index]), (double*)(&w[index + stride]), (double*)(&w[index + (stride << 1)]), (double*)(&w[index + stride * 3]), (double*)(&w[index]),  (double*)(&w[index + stride]), (double*)(&w[index + (stride << 1)]), (double*)(&w[index + stride * 3]), b, &const1, &coeff[(ic << 1) + 1], &coeff[(ic << 1) + 2], &coeff[ic * 4  + 1], &coeff[ic * 4 + 2], &coeff[ic * 6 + 1], &coeff[ic * 6 + 2]);
	} else if (ff == 5) {
	  radix5b_((double*)(&w[index]), (double*)(&w[index + stride]), (double*)(&w[index + (stride << 1)]), (double*)(&w[index + stride * 3]), (double*)(&w[index + (stride << 2)]), (double*)(&w[index]), (double*)(&w[index + stride]), (double*)(&w[index + (stride << 1)]), (double*)(&w[index + stride * 3]), (double*)(&w[index + (stride << 2)]), b, &const1, &coeff[(ic << 1) + 1], &coeff[(ic << 1) + 2], &coeff[ic * 4 + 1], &coeff[ic * 4 + 2], &coeff[ic * 6 + 1], &coeff[ic * 6 + 2], &coeff[(ic << 3) + 1], &coeff[(ic << 3) + 2]);
	} else if (ff == 7) {
	  radix7b_((double*)(&w[index]), (double*)(&w[index + stride]), (double*)(&w[index + (stride << 1)]), (double*)(&w[index + stride * 3]), (double*)(&w[index + (stride << 2)]), (double*)(&w[index + stride * 5]), (double*)(&w[index + stride * 6]), (double*)(&w[index]), (double*)(&w[index + stride]), (double*)(&w[index + (stride << 1)]), (double*)(&w[index + stride * 3]), (double*)(&w[index + (stride << 2)]), (double*)(&w[index + stride * 5]), (double*)(&w[index + stride * 6]), b, &const1);
	}
	index += *bdim;
	ic += nt * (*nprev * *nrem);
      }
    }
    stride *= ff;
    ninner *= ff;
    /* L10: */
  }
  if (*ipadok == 0) {
    return 0;
  }
  ff = f[p];
  astride = *inc * ninner;
  index = 1;
  aind = 1;
  ic = (*irow - 1) * *nrem;
  i__1 = ninner;
  for (j = 1; j <= i__1; ++j) {
    if (ff == 2) {
      radix2b_((double*)(&w[index]), (double*)(&w[index + stride]), (double*)(&a[aind]), (double*)(&a[aind + astride]), b, lda, &coeff[(ic << 1) + 1], &coeff[(ic << 1)+ 2]);
    } else if (ff == 3) {
      radix3b_((double*)(&w[index]), (double*)(&w[index + stride]), (double*)(&w[index + (stride << 1)]), (double*)(&a[aind]), (double*)(&a[aind + astride]), (double*)(&a[aind + (astride << 1)]), b, lda, &coeff[(ic << 1) + 1], &coeff[(ic << 1) + 2], & coeff[ic * 4 + 1], &coeff[ic * 4 + 2]);
    } else if (ff == 4) {
      radix4b_((double*)(&w[index]), (double*)(&w[index + stride]), (double*)(&w[index + (stride << 1)]), (double*)(&w[index + stride * 3]), (double*)(&a[aind]), (double*)(&a[aind + astride]), (double*)(&a[aind + (astride << 1)]),(double*)(&a[aind + astride * 3]), b, lda, & coeff[(ic << 1) + 1], &coeff[(ic << 1) + 2], &coeff[ic * 4 + 1], &coeff[ic * 4 + 2], &coeff[ic * 6 + 1], &coeff[ic * 6 + 2]);
    } else if (ff == 5) {
      radix5b_((double*)(&w[index]), (double*)(&w[index + stride]), (double*)(&w[index + (stride << 1)]), (double*)(&w[index + stride * 3]), (double*)(&w[index + (stride << 2)]), (double*)(&a[aind]), (double*)(&a[aind + astride]), (double*)(&a[aind + (astride << 1)]),(double*)(&a[aind + astride * 3]), (double*)(&a[aind + (astride << 2)]), b, lda, & coeff[(ic << 1) + 1], &coeff[(ic << 1) + 2], &coeff[ic * 4 + 1], &coeff[ic * 4 + 2], &coeff[ic * 6 + 1], &coeff[ic * 6 + 2], &coeff[(ic << 3) + 1], &coeff[(ic << 3) + 2]);
    } else if (ff == 7) {
      radix7b_((double*)(&w[index]), (double*)(&w[index + stride]), (double*)(&w[index + (stride << 1)]), (double*)(&w[index + stride * 3]), (double*)(&w[index + (stride << 2)]), (double*)(&w[index + stride * 5]), (double*)(&w[index + stride * 6]), (double*)(&a[aind]), (double*)(&a[aind + astride]), (double*)(&a[aind + (astride << 1)]),(double*)(&a[aind + astride * 3]), (double*)(&a[aind + (astride << 2)]), (double*)(&a[aind + astride * 5]), (double*)(&a[aind + astride * 6]), b, lda);
    }
    index += *bdim;
    ic += *nrem;
    aind += *inc;
  } 
  return 0;
} /* zfftm1rowb_ */


/* Multiple fft routine for large vectors (length>=128) */
/* Uses both L2 and L1 cache blocking. */
/* The maximum supported transform length is 8192 */
int zfftm1big_(int *sign,int * n,int * p,zomplex * a,int * inc,int * lda,int * lpref,int * nf,int * f,zomplex *coeff,int *rindex)
{
  /* System generated locals */
  int a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7;


    /* Local variables */
  int bpad, bincache;
  zomplex wbig[65536];
  int b, i__, j, k, m, iwpad, b1, n1, n2, w1, ia, if__, iw;
  int nf1;
  int pad, iiw;


  /* Max length of col transform (size of L1 cache) */
  /* Work area size (size of L2 cache) */
  /* Parameter adjustments */
  --coeff;
  a_dim1 = *lda;
  a_offset = 1 + a_dim1 * 1;
  a -= a_offset;
  --f;
  --rindex;

    /* Function Body */
  n1 = 1;
  n2 = *n;
  if__ = 1;
 L10:
  if (n1 < n2) {
    n1 *= f[if__];
    n2 /= f[if__];
    ++if__;
    goto L10;
  }
  nf1 = if__ - 1;
  /* Determine the # of transforms in each block, b1. b1 is bounded from above */
  /* by the max work area size (the size of the L2 cache) and from below by the */
  /* desire to make efficient use of loop unrolling/pipelining. We want the work */
  /* matrix and a piece of the input matrix to reside in L2 cache at the same */
  /* time, so that we can do the copy back operation without additional memory */
  /* operations. Therefore we try to limit the work area size to half of the */
  /* L2 cache size. */
  b1 = 1024 / n1;
 L30:
  /* The work area contains n2 columns of n1*b1+bpad*b1+pad elements each. */
  /* The padding is chosen so that n1+bpad and 8192/16/b1 are relatively */
  /* prime (8192 is the page size.) */
  bincache = 512 / b1;
  bpad = 0;
 L20:
  i__1 = n1 + bpad;
  if (gcd_(&i__1, &bincache) > 1) {
    ++bpad;
    goto L20;
  }
  /* If 8192/16/b1 is not an int, we have to adjust with some extra padding pad */
  /* Because the columns are roughly 2 pages long, there is a "2" below. */
  pad = (512 - b1 * bincache) << 1;
  w1 = (n1 + bpad) * b1 + pad;
  if (w1 * n2 > 32768) {
    /* If the work area turned out to be too large, first reduce the block size */
    if (b1 > 8) {
      --b1;
      goto L30;
    }
    /* If the work area is still too large, give up. */
    /* Transform lenghts up to 8192 are safe with the current cache parameters. */
    if (w1 * n2 > 65536) {
      printf("Large FFT work area too small.\n");
      exit(-1);
    }
  }
  /* Loop over fft blocks */
  i__1 = *p;
  i__2 = b1;
  for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
    /* Computing MIN */
    i__3 = b1, i__4 = *p - i__ + 1;
    b = min(i__3,i__4);
    /* Loop over columns of the work matrix */
    iw = 1;
    iwpad = 1;
    i__3 = n2;
    for (j = 1; j <= i__3; ++j) {
      /* Fill column with data */
      zfftcopyin_(&n1, &b, &b1, &a[i__ * a_dim1 + 1], inc, lda, &rindex[iw], &wbig[iwpad - 1], lpref);
      /* Do a block of column transforms */
      if (*sign >= 1) {
	zfftm1rowf_(&n1, &b, &b1, &wbig[iwpad - 1], &nf1, &f[1], (double*)(&coeff[1]), &a[a_offset], &const0, &const0, &const0, &const1, & const1, &n2);
      } else {
	zfftm1rowb_(&n1, &b, &b1, &wbig[iwpad - 1], &nf1, &f[1], (double*)(&coeff[1]), &a[a_offset], &const0, &const0, &const0, &const1, & const1, &n2);
      }
      iw += n1;
      iwpad += w1;
    }
    /* Loop over block rows of the work matrix */
    iw = 1;
    i__3 = n1;
    for (j = 1; j <= i__3; ++j) {
      /* Do a block of row transforms */
      if (*sign >= 1) {
	i__4 = *nf - nf1;
	zfftm1rowf_(&n2, &b, &w1, &wbig[iw - 1], &i__4, &f[nf1 + 1], (double*)(&coeff[1]), &a[a_offset], &const0, &const0, &const0, &n1, &j,   &const1);
      } else {
	i__4 = *nf - nf1;
	zfftm1rowb_(&n2, &b, &w1, &wbig[iw - 1], &i__4, &f[nf1 + 1], (double*)(&coeff[1]), &a[a_offset], &const0, &const0, &const0, &n1, &j,   &const1);
      }
      iw += b1;
    }
    /* Copy back data */
    if (*lda <= *inc) {
      if (*lda == 1) {
	ia = 1;
	iiw = 1;
	i__3 = n2;
	for (j = 1; j <= i__3; ++j) {
	  iw = iiw;
	  i__4 = n1;
	  for (k = 1; k <= i__4; ++k) {
	    i__5 = b;
	    for (m = 1; m <= i__5; ++m) {
	      /* An index checker would not like this */
	      i__6 = ia + m - 1 + i__ * a_dim1;
	      i__7 = iw + m - 2;
	      a[i__6].re = wbig[i__7].re, a[i__6].im = wbig[i__7]
		.im;
	    }
	    ia += *inc;
	    iw += b1;
	  }
	  iiw += w1;
	}
      } else {
	ia = 1;
	iiw = 1;
	i__3 = n2;
	for (j = 1; j <= i__3; ++j) {
	  iw = iiw;
	  i__4 = n1;
	  for (k = 1; k <= i__4; ++k) {
	    i__5 = b;
	    for (m = 1; m <= i__5; ++m) {
	      i__6 = ia + (m - 1 + i__) * a_dim1;
	      i__7 = iw + m - 2;
	      a[i__6].re = wbig[i__7].re, a[i__6].im = wbig[i__7]
		.im;
	    }
	    ia += *inc;
	    iw += b1;
	  }
	  iiw += w1;
	}
      }
    } else {
      /* lda.gt.inc case */
      if (*inc == 1) {
	i__3 = b;
	for (m = 1; m <= i__3; ++m) {
	  ia = 1;
	  iw = m;
	  i__4 = n2;
	  for (j = 1; j <= i__4; ++j) {
	    i__5 = n1;
	    for (k = 1; k <= i__5; ++k) {
	      i__6 = ia + k - 1 + (m - 1 + i__) * a_dim1;
	      i__7 = iw + (k - 1) * b1 - 1;
	      a[i__6].re = wbig[i__7].re, a[i__6].im = wbig[i__7]
		.im;
	    }
	    ia += n1;
	    iw += w1;
	  }
	}
      } else {
	i__3 = b;
	for (m = 1; m <= i__3; ++m) {
	  ia = 1;
	  iw = m;
	  i__4 = n2;
	  for (j = 1; j <= i__4; ++j) {
	    i__5 = n1;
	    for (k = 1; k <= i__5; ++k) {
	      i__6 = ia + (k - 1) * *inc + (m - 1 + i__) * 
		a_dim1;
	      i__7 = iw + (k - 1) * b1 - 1;
	      a[i__6].re = wbig[i__7].re, a[i__6].im = wbig[i__7]
		.im;
	    }
	    ia += n1 * *inc;
	    iw += w1;
	  }
	}
      }
    }
    /* next fft block */
  }
  return 0;
} /* zfftm1big_ */

/* copy out routine for the large inc case */
void zfftcopyout_(int *n, int *b, zomplex *a, int *inc, int *lda,zomplex  *w)
{
  /* System generated locals */
  int w_dim1, w_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
  int i__, j, index, ia;

  /* Parameter adjustments */
  w_dim1 = *b;
  w_offset = 1 + w_dim1 * 1;
  w -= w_offset;
  --a;

    /* Function Body */
  ia = 1;
  i__1 = *n;
  for (j = 1; j <= i__1; ++j) {
    index = ia;
    i__2 = *b;
    for (i__ = 1; i__ <= i__2; ++i__) {
      i__3 = index;
      i__4 = i__ + j * w_dim1;
      a[i__3].re = w[i__4].re, a[i__3].im = w[i__4].im;
      index += *lda;
    }
    ia += *inc;
  }
} /* zfftcopyout_ */

int zfftm1di_(int *n,zomplex * coeff)
{
  /* Local variables */
  int nf;


    /* Parameter adjustments */
  --coeff;

  /* Function Body */
  if (*n < 1) {
    printf("zfftm1di: illegal FFT dimension %i\n",(*n));
    exit(-1);
  }
  coeff[1].re = 12345678., coeff[1].im = 0.;
  coeff[2].re = (double) (*n), coeff[2].im = 0.;
  fftfact_(n, &nf, (int*)(&coeff[4]));
  coeff[3].re = (double) nf, coeff[3].im = 0.;
  fftcoeff_(n, &coeff[15]);
  return 0;
} /* zfftm1di_ */

int zfftm1d_(int *sign,int * n,int * p, zomplex* a, int *inc,int * lda, zomplex *coeff, int *lpref)
{
  /* System generated locals */
  int a_dim1, a_offset, i__1, i__2, i__3, i__4;


    /* Local variables */
  int bdim;
  int b;
  int i__;
  zomplex w[1024];
  int nf, ipadok;
  int rindex[16384];

    /* Parameter adjustments */
  --coeff;
  a_dim1 = *lda;
  a_offset = 1 + a_dim1 * 1;
  a -= a_offset;

    /* Function Body */
  if (coeff[1].re != 12345678. || coeff[1].im != 0.) {
    printf("zfftm1d: coeff not initialized.\n");
    exit(-1);
  }
  if ((double) (*n) != coeff[2].re || 0. != coeff[2].im) {
    printf("zfftm1d: coeff initialized with wrong dimension.\n");
    exit(-1);
  }
  nf = (int) coeff[3].re;
  zfftrev_(n, inc, &nf, (int*)(&coeff[4]), &const16384, rindex, &ipadok);
  if (*n <= 128) {
    /* Use L1 blocking only */
    bdim = 1024 / *n;
    i__1 = *p;
    i__2 = bdim;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
      /* Computing MIN */
      i__3 = bdim, i__4 = *p - i__ + 1;
      b = min(i__3,i__4);
      zfftcopyin_(n, &b, &b, &a[i__ * a_dim1 + 1], inc, lda, rindex, w, lpref);
      if (*sign >= 1) {
	zfftm1rowf_(n, &b, &b, w, &nf, (int*)(&coeff[4]), (double*)(&coeff[15]), &a[i__ * a_dim1 + 1], inc, lda, &ipadok, &const1, &const1, &const1);
      } else {
	zfftm1rowb_(n, &b, &b, w, &nf, (int*)(&coeff[4]), (double*)(&coeff[15]), &a[i__ * a_dim1 + 1], inc, lda, &ipadok, &const1, &const1, &const1);
      }
      if (ipadok == 0) {
	zfftcopyout_(n, &b, &a[i__ * a_dim1 + 1], inc, lda, w);
      }
    }
  } else {
    /* Use both L1 and L2 blocking */
    zfftm1big_(sign, n, p, &a[a_offset], inc, lda, lpref, &nf, (int*)(&coeff[4]), 
	       &coeff[15], rindex);
  }
  return 0;
} /* zfftm1d_ */




/* Multidimensional transforms */
int zfft2di_(int *x, int *y,zomplex * coeff)
{

  /* Parameter adjustments */
  --coeff;

  /* Function Body */
  zfftm1di_(x, &coeff[1]);
  zfftm1di_(y, &coeff[*x + 16]);
  return 0;
} /* zfft2di_ */

int zfft2d_(int *sign, int *x, int *y, zomplex *a, int *lda, zomplex *coeff)
{
  /* System generated locals */
  int a_dim1, a_offset;
  int L__1;

  /* Local variables */

    /* Parameter adjustments */
  --coeff;
  a_dim1 = *lda;
  a_offset = 1 + a_dim1 * 1;
  a -= a_offset;

    /* Function Body */
  zfftm1d_(sign, y, x, &a[a_offset], lda, &const1, &coeff[*x + 16], &constTrue);
  L__1 = *x * *y > 32768;
  zfftm1d_(sign, x, y, &a[a_offset], &const1, lda, &coeff[1], &L__1);
  return 0;
} /* zfft2d_ */

int zfft3di_(int *x, int *y, int *z__,zomplex * coeff)
{

  /* Parameter adjustments */
  --coeff;

  /* Function Body */
  zfftm1di_(x, &coeff[1]);
  zfftm1di_(y, &coeff[*x + 16]);
  zfftm1di_(z__, &coeff[*x + *y + 31]);
  return 0;
} /* zfft3di_ */

int zfft3d_(int *sign,int * x, int *y, int *z__,zomplex * a, int *la1, int *la2,zomplex * coeff)
{
  /* System generated locals */
  int a_dim1, a_dim2, a_offset, i__1, i__2;

    /* Local variables */
  int k;
  /* Parameter adjustments */
  --coeff;
  a_dim1 = *la1;
  a_dim2 = *la2;
  a_offset = 1 + a_dim1 * (1 + a_dim2 * 1);
  a -= a_offset;

  /* Function Body */
  i__1 = *z__;
  for (k = 1; k <= i__1; ++k) {
    zfftm1d_(sign, y, x, &a[(k * a_dim2 + 1) * a_dim1 + 1], la1, &const1, &
	     coeff[*x + 16], &constTrue);
    zfftm1d_(sign, x, y, &a[(k * a_dim2 + 1) * a_dim1 + 1], &const1, la1, &
	     coeff[1], &constFalse);
  }
  i__1 = *y;
  for (k = 1; k <= i__1; ++k) {
    i__2 = *la1 * *la2;
    zfftm1d_(sign, z__, x, &a[(k + a_dim2) * a_dim1 + 1], &i__2, &const1, &
	     coeff[*x + *y + 31], &constTrue);
  }
  return 0;
} /* zfft3d_ */


/* SGI-ZFFT wrapper */
zomplex *zfftm1di( int m, zomplex *save)
{
  zomplex *c = 0;
  int i=0;
  c = (zomplex*)malloc((m+15)*sizeof(zomplex));    
  zfftm1di_(&m,c);
  if(save != 0)
    for(i=0;i<m+15;i++)
      save[i] = c[i];
  return c;
}

int zfftm1d( int sign, int m, int n, zomplex *array, int incI, int incJ, zomplex *save)
{
  static int lpref = 0;
  zfftm1d_(&sign,&m,&n,array,&incI,&incJ,save,&lpref);
  return 0;
}

zomplex *zfft2di( int n1, int n2, zomplex *save)
{
  zomplex *c = 0;
  int i=0;
  c = (zomplex*)malloc((n1+n2+30)*sizeof(zomplex));    
  zfft2di_(&n1,&n2,c);
  if(save != 0)
    for(i=0;i<n1+n2+30;i++)
      save[i] = c[i];
  return c;
}

int zfft2d( int sign, int n1, int n2, zomplex *array, int ld, zomplex *save)
{
  zfft2d_(&sign,&n1,&n2,array,&ld,save);
  return 0;
}

zomplex *zfft3di( int n1, int n2, int n3, zomplex *save)
{
  zomplex *c = 0;
  int i=0;
  c = (zomplex*)malloc((n1+n2+n3+45)*sizeof(zomplex));    
  zfft3di_(&n1,&n2,&n3,c);
  if(save != 0)
    for(i=0;i<n1+n2+n3+45;i++)
      save[i] = c[i];
  return c;
}

int zfft3d( int sign, int n1, int n2, int n3, zomplex *array, int ld1, int ld2, zomplex *save)
{
  zfft3d_(&sign,&n1,&n2,&n3,array,&ld1,&ld2,save);
  return 0;
}

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* HAVE_FFT*/
