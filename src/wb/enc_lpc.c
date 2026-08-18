/*
 *===================================================================
 *  3GPP AMR Wideband Floating-point Speech Codec
 *===================================================================
 */
#include <math.h>
#include <memory.h>

#include "typedef.h"
#include "enc_util.h"

#define ORDER        16          /* order of linear prediction filter   */
#define ISF_GAP      128         /* 50                                  */
#define M            16
#define M16k         20          /* Order of LP filter                  */
#define MP1          (M+1)
#define NC16k        (M16k / 2)
#define MU           10923       /* Prediction factor (1.0/3.0) in Q15  */
#define F_MU         (1.0 / 3.0) /* prediction factor                   */
#define N_SURV_MAX   4           /* 4 survivors max                     */

#ifndef PI
#define PI           3.141592654
#endif

/* isp_isf_conversion */
#define SCALE1  (6400.0/PI)

/* chebyshev */
#define  NO_ITER   4    /* no of iterations for tracking the root */
#define  NO_POINTS 100

#define SIZE_BK1        256
#define SIZE_BK2        256
#define SIZE_BK21       64
#define SIZE_BK22       128
#define SIZE_BK23       128
#define SIZE_BK24       32
#define SIZE_BK25       32
#define SIZE_BK21_36b   128
#define SIZE_BK22_36b   128
#define SIZE_BK23_36b   64
#define SIZE_BK_NOISE1  64
#define SIZE_BK_NOISE2  64
#define SIZE_BK_NOISE3  64
#define SIZE_BK_NOISE4  32
#define SIZE_BK_NOISE5  32

extern const Word16 E_ROM_mean_isf[];
extern const Word16 E_ROM_cos[];

extern const Float32 E_ROM_dico1_isf[];
extern const Float32 E_ROM_dico2_isf[];
extern const Float32 E_ROM_dico21_isf[];
extern const Float32 E_ROM_dico22_isf[];
extern const Float32 E_ROM_dico23_isf[];
extern const Float32 E_ROM_dico24_isf[];
extern const Float32 E_ROM_dico25_isf[];
extern const Float32 E_ROM_dico21_isf_36b[];
extern const Float32 E_ROM_dico22_isf_36b[];
extern const Float32 E_ROM_dico23_isf_36b[];
extern const Float32 E_ROM_f_mean_isf[];
extern const Float32 E_ROM_lag_window[];
extern const Float32 E_ROM_grid[];
extern const Float32 E_ROM_f_interpol_frac[];

static void E_LPC_isf_reorder(Word16 *isf, Word16 min_dist, Word16 n)
{
   Word32 i, isf_min;

   isf_min = min_dist;

   for (i = 0; i < n - 1; i++)
   {
      if (isf[i] < isf_min)
      {
         isf[i] = (Word16)isf_min;
      }

      isf_min = isf[i] + min_dist;
   }
}

static void E_LPC_isp_pol_get(Word16 *isp, Word32 *f, Word32 n, Word16 k16)
{
   Word32 i, j, t0, s1, s2;
   Word16 hi, lo;

   s1 = 8388608;
   s2 = 512;

   if(k16)
   {
      s1 >>= 2;
      s2 >>= 2;
   }

   f[0] = s1;              /* f[0] = 1.0; in Q23         */
   f[1] = isp[0] * (-s2);  /* f[1] = -2.0*isp[0] in Q23  */
   f += 2;                 /* Advance f pointer          */
   isp += 2;               /* Advance isp pointer        */

   for(i = 2; i <= n; i++)
   {
      *f = f[ - 2];

      for(j = 1; j < i; j++, f--)
      {
         E_UTIL_l_extract(f[- 1], &hi, &lo);
         t0 = E_UTIL_mpy_32_16(hi, lo, *isp);   /* t0 = f[-1] * isp */
         t0 = (t0 << 1);
         *f = (*f - t0);         /* *f -= t0    */
         *f = (*f + f[ - 2]);    /* *f += f[-2] */
      }

      *f = *f - (*isp * s2);     /* *f -= isp << 8 */
      f += i;     /* Advance f pointer   */
      isp += 2;   /* Advance isp pointer */
   }
}

static void E_LPC_f_isp_pol_get(Float32 isp[], Float32 f[], Word32 n)
{
   Float32 b;
   Word32 i, j;

   f[0] = 1;
   b = (Float32)(-2.0 * *isp);
   f[1] = b;

   for (i = 2; i <= n; i++)
   {
      isp += 2;
      b = (Float32)(-2.0 * *isp);
      f[i] = (Float32)(b * f[i - 1] + 2.0 * f[i - 2]);

      for (j = i - 1; j > 1; j--)
      {
         f[j] += b * f[j - 1] + f[j - 2];
      }

      f[1] += b;
   }
}

void E_LPC_isp_a_conversion(Word16 isp[], Word16 a[], Word16 m)
{
   Word32 f1[NC16k + 1], f2[NC16k];
   Word32 i, j, nc, t0;
   Word16 hi, lo;

   nc = m >> 1;

   if (nc > 8)
   {
      E_LPC_isp_pol_get(&isp[0], f1, nc, 1);

      for (i = 0; i <= nc; i++)
      {
         f1[i] = (f1[i] << 2);
      }
   }
   else
   {
      E_LPC_isp_pol_get(&isp[0], f1, nc, 0);
   }

   if (nc > 8)
   {
      E_LPC_isp_pol_get(&isp[1], f2, nc - 1, 1);

      for (i = 0; i <= nc - 1; i++)
      {
         f2[i] = (f2[i] << 2);
      }
   }
   else
   {
      E_LPC_isp_pol_get(&isp[1], f2, nc - 1, 0);
   }

   for (i = (nc - 1); i > 1; i--)
   {
      f2[i] = f2[i] - f2[i - 2];     /* f2[i] -= f2[i-2]; */
   }

   for (i = 0; i < nc; i++)
   {
      E_UTIL_l_extract(f1[i], &hi, &lo);
      t0 = E_UTIL_mpy_32_16(hi, lo, isp[m - 1]);
      f1[i] = f1[i] + t0;

      E_UTIL_l_extract(f2[i], &hi, &lo);
      t0 = E_UTIL_mpy_32_16(hi, lo, isp[m - 1]);
      f2[i] = f2[i] - t0;
   }

   a[0] = 4096;

   for (i = 1, j = (m - 1); i < nc; i++, j--)
   {
      t0 = f1[i] + f2[i];                    /* f1[i] + f2[i]             */
      a[i] = (Word16)((t0 + 0x800) >> 12);   /* from Q23 to Q12 and * 0.5 */

      t0 = (f1[i] - f2[i]);                  /* f1[i] - f2[i]             */
      a[j] = (Word16)((t0 + 0x800) >> 12);   /* from Q23 to Q12 and * 0.5 */
   }

   E_UTIL_l_extract(f1[nc], &hi, &lo);
   t0 = E_UTIL_mpy_32_16(hi, lo, isp[m - 1]);
   t0 = (f1[nc] + t0);
   a[nc] = (Word16)((t0 + 0x800) >> 12);    /* from Q23 to Q12 and * 0.5 */

   a[m] = (Word16)((isp[m - 1] + 0x4) >> 3); /* from Q15 to Q12          */
}

void E_LPC_f_isp_a_conversion(Float32 *isp, Float32 *a, Word32 m)
{
   Float32 f1[(M16k / 2) + 1], f2[M16k / 2];
   Word32 i, j, nc;

   nc = m / 2;

   E_LPC_f_isp_pol_get(&isp[0], f1, nc);
   E_LPC_f_isp_pol_get(&isp[1], f2, nc-1);

   for (i = (nc - 1); i > 1; i--)
   {
      f2[i] -= f2[i - 2];
   }

   for (i = 0; i < nc; i++)
   {
      f1[i] *= (Float32)(1.0 + isp[m - 1]);
      f2[i] *= (Float32)(1.0 - isp[m - 1]);
   }

   a[0] = 1.0;

   for (i = 1, j = m - 1; i < nc; i++, j--)
   {
      a[i] = (Float32)(0.5 * (f1[i] + f2[i]));
      a[j] = (Float32)(0.5 * (f1[i] - f2[i]));
   }

   a[nc] = (Float32)(0.5 * f1[nc] * (1.0 + isp[m - 1]));
   a[m] = isp[m - 1];
}

void E_LPC_int_isp_find(Word16 isp_old[], Word16 isp_new[],
                        const Word16 frac[], Word16 Az[])
{
   Word32 i, k, fac_old, fac_new, tmp;
   Word16 isp[M];

   for (k = 0; k < 3; k++)
   {
      fac_new = frac[k];
      fac_old = ((32767 - fac_new) + 1);  /* 1.0 - fac_new */

      for (i = 0; i < M; i++)
      {
         tmp = isp_old[i] * fac_old;
         tmp += (isp_new[i] * fac_new);
         isp[i] = (Word16)((tmp + 0x4000) >> 15);
      }

      E_LPC_isp_a_conversion(isp, Az, M);
      Az += MP1;
   }

   E_LPC_isp_a_conversion(isp_new, Az, M);
}

void E_LPC_f_int_isp_find(Float32 isp_old[], Float32 isp_new[], Float32 a[],
                          Word32 nb_subfr, Word32 m)
{
  Float32 isp[M], fnew, fold;
  Float32 *p_a;
  Word32 i, k;

  p_a = a;

  for (k = 0; k < nb_subfr; k++)
  {
    fnew = E_ROM_f_interpol_frac[k];
    fold = (Float32)(1.0 - fnew);

    for (i = 0; i < m; i++)
    {
      isp[i] = isp_old[i] * fold + isp_new[i] * fnew;
    }

    E_LPC_f_isp_a_conversion(isp, p_a, m);
    p_a += (m + 1);
  }
}

void E_LPC_a_weight(Float32 *a, Float32 *ap, Float32 gamma, Word32 m)
{
   Float32 f;
   Word32 i;

   ap[0] = a[0];
   f = gamma;

   for (i = 1; i <= m; i++)
   {
      ap[i] = f*a[i];
      f *= gamma;
   }
}

static void E_LPC_isf_2s3s_decode(Word32 *indice, Word16 *isf_q, Word16 *past_isfq)
{
   Word32 i;
   Word16 tmp;

   for(i = 0; i < 9; i++)
   {
      isf_q[i] = (Word16)((E_ROM_dico1_isf[indice[0] * 9 + i] * 2.56F) + 0.5F);
   }

   for(i = 0; i < 7; i++)
   {
      isf_q[i + 9] = (Word16)((E_ROM_dico2_isf[indice[1] * 7 + i] * 2.56F) + 0.5F);
   }

   for(i = 0; i < 5; i++)
   {
      isf_q[i] = (Word16)(isf_q[i] + (Word16)((E_ROM_dico21_isf_36b[indice[2] * 5 + i] * 2.56F) + 0.5F));
   }

   for(i = 0; i < 4; i++)
   {
      isf_q[i + 5] = (Word16)(isf_q[i + 5] + (Word32)((E_ROM_dico22_isf_36b[indice[3] * 4 + i] * 2.56F) + 0.5F));
   }

   for(i = 0; i < 7; i++)
   {
      isf_q[i + 9] = (Word16)(isf_q[i + 9] + (Word32)((E_ROM_dico23_isf_36b[indice[4] * 7 + i] * 2.56F) + 0.5F));
   }

   for(i = 0; i < ORDER; i++)
   {
      tmp = isf_q[i];
      isf_q[i] = (Word16)(tmp + E_ROM_mean_isf[i]);
      isf_q[i] = (Word16)(isf_q[i] + ((MU * past_isfq[i]) >> 15));
      past_isfq[i] = tmp;
   }

   E_LPC_isf_reorder(isf_q, ISF_GAP, ORDER);
}

void E_LPC_isf_2s5s_decode(Word32 *indice, Word16 *isf_q, Word16 *past_isfq)
{
   Word32 i;
   Word16 tmp;

   for (i = 0; i < 9; i++)
   {
      isf_q[i] = (Word16)((E_ROM_dico1_isf[indice[0] * 9 + i] * 2.56F) + 0.5F);
   }

   for (i = 0; i < 7; i++)
   {
      isf_q[i + 9] =
         (Word16)((E_ROM_dico2_isf[indice[1] * 7 + i] * 2.56F) + 0.5F);
   }

   for (i = 0; i < 3; i++)
   {
      isf_q[i] = (Word16)(isf_q[i] +
         (Word32)((E_ROM_dico21_isf[indice[2] * 3 + i] * 2.56F) + 0.5F));
   }

   for (i = 0; i < 3; i++)
   {
      isf_q[i + 3] = (Word16)(isf_q[i + 3] +
         (Word32)((E_ROM_dico22_isf[indice[3] * 3 + i] * 2.56F) + 0.5F));
   }

   for (i = 0; i < 3; i++)
   {
      isf_q[i + 6] = (Word16)(isf_q[i + 6] +
         (Word32)((E_ROM_dico23_isf[indice[4] * 3 + i] * 2.56F) + 0.5F));
   }

   for (i = 0; i < 3; i++)
   {
      isf_q[i + 9] = (Word16)(isf_q[i + 9] +
         (Word32)((E_ROM_dico24_isf[indice[5] * 3 + i] * 2.56F) + 0.5F));
   }

   for (i = 0; i < 4; i++)
   {
      isf_q[i + 12] = (Word16)(isf_q[i + 12] +
         (Word32)((E_ROM_dico25_isf[indice[6] * 4 + i] * 2.56F) + 0.5F));
   }

   for (i = 0; i < ORDER; i++)
   {
      tmp = isf_q[i];
      isf_q[i] = (Word16)(tmp + E_ROM_mean_isf[i]);
      isf_q[i] = (Word16)(isf_q[i] + ((MU * past_isfq[i]) >> 15));
      past_isfq[i] = tmp;
   }

   E_LPC_isf_reorder(isf_q, ISF_GAP, ORDER);
}

void E_LPC_isf_isp_conversion(Word16 isf[], Word16 isp[], Word16 m)
{
   Word32 i, ind, offset, tmp;

   for (i = 0; i < m - 1; i++)
   {
      isp[i] = isf[i];
   }

   isp[m - 1] = (Word16)(isf[m - 1] << 1);

   for (i = 0; i < m; i++)
   {
      ind = isp[i] >> 7;         /* ind    = b7-b15 of isf[i] */
      offset = isp[i] & 0x007f;  /* offset = b0-b6  of isf[i] */

      tmp = ((E_ROM_cos[ind + 1] - E_ROM_cos[ind]) * offset) << 1;
      isp[i] = (Word16)(E_ROM_cos[ind] + (tmp >> 8));
   }
}

static Float32 E_LPC_chebyshev(Float32 x, Float32 *f, Word32 n)
{
   Float32 b1, b2, b0, x2;
   Word32 i;                        /* for the special case of 10th order  */
                                    /*       filter (n=5)                  */
   x2 = 2.0F * x;                   /* x2 = 2.0*x;                         */
   b2 = f[0];                       /* b2 = f[0];                          */
   b1 = x2 * b2 + f[1];             /* b1 = x2*b2 + f[1];                  */

   for (i = 2; i < n; i++)
   {
      b0 = x2 * b1 - b2 + f[i];     /* b0 = x2 * b1 - 1. + f[2];           */
      b2 = b1;                      /* b2 = x2 * b0 - b1 + f[3];           */
      b1 = b0;                      /* b1 = x2 * b2 - b0 + f[4];           */
   }

   return (x * b1 - b2 + 0.5F * f[n]); /* return (x*b1 - b2 + 0.5*f[5]);   */
}

Word16 E_LPC_isf_sub_vq(Float32 *x, const Float32 *E_ROM_dico, Word32 dim,
                        Word32 E_ROM_dico_size, Float32 *distance)
{
   Float32 dist_min, dist, temp;
   const Float32 *p_E_ROM_dico;
   Word32 i, j, index = 0;

   dist_min = 1.0e30f;
   p_E_ROM_dico = E_ROM_dico;

   for (i = 0; i < E_ROM_dico_size; i++)
   {
      dist = x[0] - *p_E_ROM_dico++;
      dist *= dist;

      for (j = 1; j < dim; j++)
      {
         temp = x[j] - *p_E_ROM_dico++;
         dist += temp * temp;
      }

      if (dist < dist_min)
      {
         dist_min = dist;
         index = i;
      }
   }

   *distance = dist_min;
   memcpy(x, &E_ROM_dico[index * dim], dim * sizeof(Float32));

   return (Word16)index;
}

void E_LPC_lag_wind(Float32 r[], Word32 m)
{
   Word32 i;

   for (i = 0; i < m; i++)
   {
      r[i] *= E_ROM_lag_window[i];
   }
}

void E_LPC_lev_dur(Float32 *a, Float32 *r, Word32 m)
{
   Float32 buf[M];
   Float32 *rc;    /* reflection coefficients  0,...,m-1 */
   Float32 s, at, err;
   Word32 i, j, l;

   rc = &buf[0];

   rc[0] = (-r[1]) / r[0];
   a[0] = 1.0F;
   a[1] = rc[0];
   err = r[0] + r[1] * rc[0];

   for (i = 2; i <= m; i++)
   {
      s = 0.0F;

      for (j = 0; j < i; j++)
      {
         s += r[i - j] * a[j];
      }

      rc[i - 1]= (-s) / (err);

      for (j = 1; j <= (i >> 1); j++)
      {
         l = i - j;
         at = a[j] + rc[i - 1] * a[l];
         a[l] += rc[i - 1] * a[j];
         a[j] = at;
      }

      a[i] = rc[i - 1];
      err += rc[i - 1] * s;

      if (err <= 0.0F)
      {
         err = 0.01F;
      }
   }
}

void E_LPC_a_isp_conversion(Float32 *a, Float32 *isp, Float32 *old_isp,
                            Word32 m)
{
   Float32 f1[(M >> 1) + 1], f2[M >> 1];
   Float32 *pf;
   Float32 xlow, ylow, xhigh, yhigh, xmid, ymid, xint;
   Word32 j, i, nf, ip, order, nc;

   nc = m >> 1;

   for (i=0; i < nc; i++)
   {
      f1[i] = a[i] + a[m - i];
      f2[i] = a[i] - a[m - i];
   }

   f1[nc] = 2.0F * a[nc];

   for (i = 2; i < nc; i++)
   {
      f2[i] += f2[i - 2];
   }

   nf=0;      /* number of found frequencies */
   ip=0;      /* flag to first polynomial   */

   pf = f1;  /* start with F1(z) */
   order = nc;

   xlow = E_ROM_grid[0];
   ylow = E_LPC_chebyshev(xlow, pf, nc);

   j = 0;

   while ( (nf < m - 1) && (j < NO_POINTS) )
   {
      j++;
      xhigh = xlow;
      yhigh = ylow;
      xlow = E_ROM_grid[j];
      ylow = E_LPC_chebyshev(xlow, pf, order);

      if (ylow * yhigh <= 0.0F)  /* if sign change new root exists */
      {
         j--;

         for (i = 0; i < NO_ITER; i++)
         {
            xmid = 0.5F * (xlow + xhigh);
            ymid = E_LPC_chebyshev(xmid, pf, order);

            if (ylow * ymid <= 0.0F)
            {
               yhigh = ymid;
               xhigh = xmid;
            }
            else
            {
               ylow = ymid;
               xlow = xmid;
            }
         }

         xint = xlow - ylow * (xhigh - xlow) / (yhigh - ylow);

         isp[nf] = xint;    /* new root */
         nf++;

         ip = 1 - ip;                  /* flag to other polynomial    */
         pf = ip ? f2 : f1;            /* pointer to other polynomial */
         order = ip ? (nc - 1) : nc;   /* order of other polynomial   */

         xlow = xint;
         ylow = E_LPC_chebyshev(xlow, pf, order);
      }
   }

   isp[m - 1] = a[m];

   if (nf < m - 1)
   {
      for(i = 0; i < m; i++)
      {
         isp[i] = old_isp[i];
      }
   }
}

void E_LPC_isp_isf_conversion(Float32 isp[], Float32 isf[], Word32 m)
{
   Word32 i;

   for(i = 0; i < (m - 1); i++)
   {
      isf[i] = (Float32)(acos(isp[i]) * SCALE1);
   }

   isf[m - 1] = (Float32)(acos(isp[m - 1]) * SCALE1 * 0.5F);
}

static void E_LPC_stage1_isf_vq(Float32 *x, const Float32 *E_ROM_dico,
                                Word32 dim, Word32 E_ROM_dico_size,
                                Word32 *index, Word32 surv)
{
   Float32 dist_min[N_SURV_MAX];
   Float32 dist, temp1, temp2;
   const Float32 *p_E_ROM_dico;
   Word32 i, j, k, l;

   for (i = 0; i < surv; i++)
   {
      dist_min[i] = 1.0e30F;
   }

   for (i = 0; i < surv; i++)
   {
      index[i] = i;
   }

   p_E_ROM_dico = E_ROM_dico;

   for (i = 0; i < E_ROM_dico_size; i++)
   {
      dist = x[0] - *p_E_ROM_dico++;
      dist *= dist;

      for (j = 1; j < dim; j += 2)
      {
         temp1 = x[j] - *p_E_ROM_dico++;
         temp2 = x[j + 1] - *p_E_ROM_dico++;
         dist += temp1 * temp1 + temp2 * temp2;
      }

      for (k = 0; k < surv; k++)
      {
         if (dist < dist_min[k])
         {
            for (l = surv - 1; l > k; l--)
            {
               dist_min[l] = dist_min[l - 1];
               index[l] = index[l - 1];
            }

            dist_min[k] = dist;
            index[k] = i;

            break;
         }
      }
   }

}

void E_LPC_isf_2s3s_quantise(Float32 *isf1, Word16 *isf_q, Word16 *past_isfq,
                             Word32 *indice, Word32 nb_surv)
{
   Float32 isf[ORDER], isf_stage2[ORDER];
   Float32 temp, min_err, distance;
   Word32 surv1[N_SURV_MAX];     /* indices of survivors from 1st stage */
   Word32 tmp_ind[5];
   Word32 i, k;

   for (i = 0; i < ORDER; i++)
   {
      isf[i] = (Float32)((isf1[i] - E_ROM_f_mean_isf[i]) -
         F_MU * past_isfq[i] * 0.390625F);
   }

   E_LPC_stage1_isf_vq(&isf[0], E_ROM_dico1_isf, 9, SIZE_BK1, surv1, nb_surv);

   distance = 1.0e30F;

   for (k = 0; k < nb_surv; k++)
   {
      for (i = 0; i < 9; i++)
      {
         isf_stage2[i] = isf[i] - E_ROM_dico1_isf[i + surv1[k] * 9];
      }

      tmp_ind[0] = E_LPC_isf_sub_vq(&isf_stage2[0], E_ROM_dico21_isf_36b, 5,
         SIZE_BK21_36b, &min_err);
      temp = min_err;
      tmp_ind[1] = E_LPC_isf_sub_vq(&isf_stage2[5], E_ROM_dico22_isf_36b, 4,
         SIZE_BK22_36b, &min_err);
      temp += min_err;

      if (temp < distance)
      {
         distance = temp;
         indice[0] = surv1[k];

         for (i = 0; i < 2; i++)
         {
            indice[i + 2] = tmp_ind[i];
         }
      }
   }

   E_LPC_stage1_isf_vq(&isf[9], E_ROM_dico2_isf, 7, SIZE_BK2, surv1, nb_surv);

   distance = 1.0e30F;

   for (k = 0; k < nb_surv; k++)
   {
      for (i = 0; i < 7; i++)
      {
         isf_stage2[i] = isf[9 + i] - E_ROM_dico2_isf[i + surv1[k] * 7];
      }

      tmp_ind[0] = E_LPC_isf_sub_vq(&isf_stage2[0], E_ROM_dico23_isf_36b, 7,
         SIZE_BK23_36b, &min_err);
      temp = min_err;

      if (temp < distance)
      {
         distance = temp;
         indice[1] = surv1[k];
         indice[4]= tmp_ind[0];
      }
   }

   E_LPC_isf_2s3s_decode(indice, isf_q, past_isfq);
}

void E_LPC_isf_2s5s_quantise(Float32 *isf1, Word16 *isf_q, Word16 *past_isfq,
                             Word32 *indice, Word32 nb_surv)
{
   Float32 isf[ORDER], isf_stage2[ORDER];
   Float32 temp, min_err, distance;
   Word32 surv1[N_SURV_MAX];     /* indices of survivors from 1st stage */
   Word32 tmp_ind[5];
   Word32 i, k;

   for (i=0; i<ORDER; i++)
   {
      isf[i] = (Float32)((isf1[i] - E_ROM_f_mean_isf[i]) -
         F_MU * past_isfq[i] * 0.390625F);
   }

   E_LPC_stage1_isf_vq(&isf[0], E_ROM_dico1_isf, 9, SIZE_BK1, surv1, nb_surv);

   distance = 1.0e30F;

   for (k = 0; k < nb_surv; k++)
   {
      for (i = 0; i < 9; i++)
      {
         isf_stage2[i] = isf[i] - E_ROM_dico1_isf[i + surv1[k] * 9];
      }

      tmp_ind[0] = E_LPC_isf_sub_vq(&isf_stage2[0], E_ROM_dico21_isf, 3,
         SIZE_BK21, &min_err);
      temp = min_err;
      tmp_ind[1] = E_LPC_isf_sub_vq(&isf_stage2[3], E_ROM_dico22_isf, 3,
         SIZE_BK22, &min_err);
      temp += min_err;
      tmp_ind[2] = E_LPC_isf_sub_vq(&isf_stage2[6], E_ROM_dico23_isf, 3,
         SIZE_BK23, &min_err);
      temp += min_err;

      if (temp < distance)
      {
         distance = temp;
         indice[0] = surv1[k];

         for (i = 0; i < 3; i++)
         {
            indice[i + 2] = tmp_ind[i];
         }
      }
   }

   E_LPC_stage1_isf_vq(&isf[9], E_ROM_dico2_isf, 7, SIZE_BK2, surv1, nb_surv);

   distance = 1.0e30F;

   for (k=0; k<nb_surv; k++)
   {
      for (i = 0; i < 7; i++)
      {
         isf_stage2[i] = isf[9+i] - E_ROM_dico2_isf[i+surv1[k]*7];
      }

      tmp_ind[0] = E_LPC_isf_sub_vq(&isf_stage2[0], E_ROM_dico24_isf, 3,
         SIZE_BK24, &min_err);
      temp = min_err;
      tmp_ind[1] = E_LPC_isf_sub_vq(&isf_stage2[3], E_ROM_dico25_isf, 4,
         SIZE_BK25, &min_err);
      temp += min_err;

      if (temp < distance)
      {
         distance = temp;
         indice[1] = surv1[k];

         for (i = 0; i < 2; i++)
         {
            indice[i + 5]= tmp_ind[i];
         }
      }
   }

   E_LPC_isf_2s5s_decode(indice, isf_q, past_isfq);
}
