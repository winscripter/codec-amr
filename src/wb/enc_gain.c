/*
 *===================================================================
 *  3GPP AMR Wideband Floating-point Speech Codec
 *===================================================================
 */
#include <math.h>
#include <memory.h>
#include "typedef.h"
#include "enc_util.h"

#define L_FRAME         256   /* Frame size                                */
#define L_SUBFR         64    /* Subframe size                             */
#define HP_ORDER        3
#define L_INTERPOL1     4
#define L_INTERPOL2     16
#define PIT_SHARP       27853 /* pitch sharpening factor = 0.85 Q15        */
#define F_PIT_SHARP     0.85F /* pitch sharpening factor                   */
#define PIT_MIN         34    /* Minimum pitch lag with resolution 1/4     */
#define UP_SAMP         4
#define DIST_ISF_MAX    120
#define DIST_ISF_THRES  60
#define GAIN_PIT_THRES  0.9F
#define GAIN_PIT_MIN    0.6F

extern const Float32 E_ROM_corrweight[];
extern const Float32 E_ROM_inter4_1[];
extern const Word16 E_ROM_inter4_2[];

void E_GAIN_clip_init(Float32 mem[])
{
   mem[0] = DIST_ISF_MAX;
   mem[1] = GAIN_PIT_MIN;
}

Word32 E_GAIN_clip_test(Float32 mem[])
{
   Word32 clip;

   clip = 0;
   if ((mem[0] < DIST_ISF_THRES) && (mem[1] > GAIN_PIT_THRES))
   {
      clip = 1;
   }

   return (clip);
}

void E_GAIN_clip_isf_test(Float32 isf[], Float32 mem[])
{
   Word32 i;
   Float32 dist, dist_min;

   dist_min = isf[1] - isf[0];

   for (i = 2; i < M - 1; i++)
   {
      dist = isf[i] - isf[i-1];

      if (dist < dist_min)
      {
         dist_min = dist;
      }
   }

   dist = 0.8F * mem[0] + 0.2F * dist_min;

   if (dist > DIST_ISF_MAX)
   {
      dist = DIST_ISF_MAX;
   }

   mem[0] = dist;
}

void E_GAIN_clip_pit_test(Float32 gain_pit, Float32 mem[])
{
   Float32 gain;

   gain = 0.9F * mem[1] + 0.1F * gain_pit;

   if (gain < GAIN_PIT_MIN)
   {
      gain = GAIN_PIT_MIN;
   }

   mem[1] = gain;
}

void E_GAIN_lp_decim2(Float32 x[], Word32 l, Float32 *mem)
{
   Float32 x_buf[L_FRAME + 3];
   Float32 temp;
   Word32 i, j;

   memcpy(x_buf, mem, 3 * sizeof(Float32));
   memcpy(&x_buf[3], x, l * sizeof(Float32));

   for (i = 0; i < 3; i++)
   {
      mem[i] = ((x[l - 3 + i] > 1e-10) | (x[l - 3 + i] < -1e-10)) ? x[l - 3 + i] : 0;
   }

   for (i = 0, j = 0; i < l; i += 2, j++)
   {
      temp = x_buf[i] * 0.13F;
      temp += x_buf[i + 1] * 0.23F;
      temp += x_buf[i + 2] * 0.28F;
      temp += x_buf[i + 3] * 0.23F;
      temp += x_buf[i + 4] * 0.13F;
      x[j] = temp;
   }
}

Word32 E_GAIN_open_loop_search(Float32 *wsp, Word32 L_min, Word32 L_max,
                           Word32 nFrame, Word32 L_0, Float32 *gain,
                           Float32 *hp_wsp_mem, Float32 hp_old_wsp[],
                           UWord8 weight_flg)
{
   Word32  i, j, k, L = 0;
   Float32  o, R0, R1, R2, R0_max = -1.0e23f;
   const Float32 *ww, *we;
   Float32 *data_a, *data_b, *hp_wsp, *p, *p1;

   ww = &E_ROM_corrweight[198];
   we = &E_ROM_corrweight[98 + L_max - L_0];

   for (i = L_max; i > L_min; i--)
   {
      p  = &wsp[0];
      p1 = &wsp[-i];

      R0 = 0.0;

      for (j = 0; j < nFrame; j += 2)
      {
         R0 += p[j] * p1[j];
         R0 += p[j + 1] * p1[j + 1];
      }

      R0 *= *ww--;

      if ((L_0 > 0) & (weight_flg == 1))
      {
         R0 *= *we--;
      }

      if (R0 >= R0_max)
      {
         R0_max = R0;
         L = i;
      }
   }

   data_a = hp_wsp_mem;
   data_b = hp_wsp_mem + HP_ORDER;
   hp_wsp = hp_old_wsp + L_max;

   for (k = 0; k < nFrame; k++)
   {
      data_b[0] = data_b[1];
      data_b[1] = data_b[2];
      data_b[2] = data_b[3];

      data_b[HP_ORDER] = wsp[k];

      o = data_b[0] * 0.83787057505665F;
      o += data_b[1] * -2.50975570071058F;
      o += data_b[2] * 2.50975570071058F;
      o += data_b[3] * -0.83787057505665F;

      o -= data_a[0] * -2.64436711600664F;
      o -= data_a[1] * 2.35087386625360F;
      o -= data_a[2] * -0.70001156927424F;

      data_a[2] = data_a[1];
      data_a[1] = data_a[0];

      data_a[0] = o;
      hp_wsp[k] = o;
   }

   p  = &hp_wsp[0];
   p1 = &hp_wsp[-L];

   R0 = 0.0F;
   R1 = 0.0F;
   R2 = 0.0F;

   for (j = 0; j < nFrame; j++)
   {
      R1 += p1[j] * p1[j];
      R2 += p[j] * p[j];
      R0 += p[j] * p1[j];
   }

   *gain = (Float32)(R0 / (sqrt(R1 * R2) + 1e-5));

   memcpy(hp_old_wsp, &hp_old_wsp[nFrame], L_max * sizeof(Float32));

   return(L);
}

static void E_GAIN_sort(Word32 n, Word32 *ra)
{
   Word32 l, j, ir, i, rra;

   l = (n >> 1) + 1;
   ir = n;
   for (;;)
   {
      if (l > 1)
      {
         rra = ra[--l];
      }
      else
      {
         rra = ra[ir];
         ra[ir] = ra[1];

         if (--ir == 1)
         {
            ra[1] = rra;
            return;
         }
      }

      i = l;
      j = l << 1;
      while (j <= ir)
      {
         if (j < ir && ra[j] < ra[j + 1])
         {
            ++j;
         }

         if (rra < ra[j])
         {
            ra[i] = ra[j];
            j += (i = j);
         }
         else
         {
            j = ir + 1;
         }
      }

      ra[i] = rra;
   }
}

Word32 E_GAIN_olag_median(Word32 prev_ol_lag, Word32 old_ol_lag[5])
{
   Word32 tmp[6] = {0};
   Word32 i;

   for (i = 4; i > 0; i--)
   {
      old_ol_lag[i] = old_ol_lag[i-1];
   }

   old_ol_lag[0] = prev_ol_lag;

   for (i = 0; i < 5; i++)
   {
      tmp[i+1] = old_ol_lag[i];
   }

   E_GAIN_sort(5, tmp);

   return tmp[3];
}

static void E_GAIN_norm_corr(Float32 exc[], Float32 xn[], Float32 h[],
                             Word32 t_min, Word32 t_max, Float32 corr_norm[])
{
   Float32 excf[L_SUBFR];  /* filtered past excitation */
   Float32 alp, ps, norm;
   Word32 t, j, k;

   k = - t_min;
   E_UTIL_f_convolve(&exc[k], h, excf);

   for (t = t_min; t <= t_max; t++)
   {
      ps = 0.0F;
      alp = 0.01F;

      for (j = 0; j < L_SUBFR; j++)
      {
         ps += xn[j] * excf[j];
         alp += excf[j] * excf[j];
      }

      norm = (Float32)(1.0F / sqrt(alp));
      corr_norm[t] = ps * norm;

      if (t != t_max)
      {
         k--;

         for (j = L_SUBFR - 1; j > 0; j--)
         {
            excf[j] = excf[j - 1] + exc[k] * h[j];
         }

         excf[0] = exc[k];
      }
   }
}

static Float32 E_GAIN_norm_corr_interpolate(Float32 *x, Word32 frac)
{
   Float32 s, *x1, *x2;
   const Float32 *c1, *c2;

   if (frac < 0)
   {
      frac += 4;
      x--;
   }

   x1 = &x[0];
   x2 = &x[1];
   c1 = &E_ROM_inter4_1[frac];
   c2 = &E_ROM_inter4_1[4 - frac];

   s = x1[0] * c1[0] + x2[0] * c2[0];
   s += x1[-1] * c1[4] + x2[1] * c2[4];
   s += x1[-2] * c1[8] + x2[2] * c2[8];
   s += x1[-3] * c1[12] + x2[3] * c2[12];

   return s;
}

Word32 E_GAIN_closed_loop_search(Float32 exc[], Float32 xn[], Float32 h[],
                             Word32 t0_min, Word32 t0_max, Word32 *pit_frac,
                             Word32 i_subfr, Word32 t0_fr2, Word32 t0_fr1)
{
   Float32 corr_v[15 + 2 * L_INTERPOL1 + 1];
   Float32 cor_max, max, temp;
   Float32 *corr;
   Word32 i, fraction, step;
   Word32 t0, t_min, t_max;

   t_min = t0_min - L_INTERPOL1;
   t_max = t0_max + L_INTERPOL1;
   corr = &corr_v[-t_min];
   E_GAIN_norm_corr(exc, xn, h, t_min, t_max, corr);
   max = corr[t0_min];
   t0  = t0_min;

   for(i = t0_min + 1; i <= t0_max; i++)
   {
      if( corr[i] > max)
      {
         max = corr[i];
         t0 = i;
      }
   }

   if((i_subfr == 0) & (t0 >= t0_fr1))
   {
      *pit_frac = 0;
      return(t0);
   }

   step = 1;                /* 1/4 subsample resolution */
   fraction = -3;
   if (((i_subfr == 0) & (t0 >= t0_fr2)) | (t0_fr2 == PIT_MIN))
   {
      step = 2;              /* 1/2 subsample resolution */
      fraction = -2;
   }

   if (t0 == t0_min)
   {
      fraction = 0;
   }

   cor_max = E_GAIN_norm_corr_interpolate(&corr[t0], fraction);

   for (i = (fraction + step); i <= 3; i += step)
   {
      temp = E_GAIN_norm_corr_interpolate(&corr[t0], i);

      if (temp > cor_max)
      {
         cor_max = temp;
         fraction = i;
      }
   }

   if (fraction < 0)
   {
      fraction += 4;
      t0 -= 1;
   }

   *pit_frac = fraction;

   return (t0);
}

void E_GAIN_adaptive_codebook_excitation(Word16 exc[], Word16 T0, Word32 frac, Word16 L_subfr)
{
   Word32 i, j, k, L_sum;
   Word16 *x;

   x = &exc[-T0];

   frac = -(frac);

   if (frac < 0)
   {
      frac = (frac + UP_SAMP);
      x--;
   }

   x = x - L_INTERPOL2 + 1;

   for (j = 0; j < L_subfr; j++)
   {
      L_sum = 0L;

      for (i = 0, k = ((UP_SAMP - 1) - frac); i < 2 * L_INTERPOL2; i++, k += UP_SAMP)
      {
         L_sum = L_sum + (x[i] * E_ROM_inter4_2[k]);
      }

      L_sum = (L_sum + 0x2000) >> 14;
      exc[j] = E_UTIL_saturate(L_sum);
      x++;
   }
}

void E_GAIN_pitch_sharpening(Word16 *x, Word16 pit_lag)
{
   Word32 L_tmp, i;

   for (i = pit_lag; i < L_SUBFR; i++)
   {
      L_tmp = x[i] << 15;
      L_tmp += x[i - pit_lag] * PIT_SHARP;
      x[i] = (Word16)((L_tmp + 0x4000) >> 15);
   }
}

void E_GAIN_f_pitch_sharpening(Float32 *x, Word32 pit_lag)
{
   Word32 i;

   for (i = pit_lag; i < L_SUBFR; i++)
   {
      x[i] += x[i - pit_lag] * F_PIT_SHARP;
   }
}

Word32 E_GAIN_voice_factor(Word16 exc[], Word16 Q_exc, Word16 gain_pit,
                          Word16 code[], Word16 gain_code)
{

   Word32 i, L_tmp, tmp, exp, ener1, exp1, ener2, exp2;

   ener1 = E_UTIL_dot_product12(exc, exc, L_SUBFR, &exp1) >> 16;
   exp1 = exp1 - (Q_exc + Q_exc);
   L_tmp = (gain_pit * gain_pit) << 1;
   exp = E_UTIL_norm_l(L_tmp);
   tmp = (L_tmp << exp) >> 16;
   ener1 = (ener1 * tmp) >> 15;
   exp1 = (exp1 - exp) - 10;        /* 10 -> gain_pit Q14 to Q9 */

   ener2 = E_UTIL_dot_product12(code, code, L_SUBFR, &exp2) >> 16;

   exp = E_UTIL_norm_s(gain_code);
   tmp = gain_code << exp;
   tmp = (tmp * tmp) >> 15;
   ener2 = (ener2 * tmp) >> 15;
   exp2 = exp2 - (exp + exp);

   i = exp1 - exp2;

   if (i >= 0)
   {
      ener1 = ener1 >> 1;
      ener2 = ener2 >> (i + 1);
   }
   else
   {
      i = 1 - i;
      if (i < 32)
      {
         ener1 = ener1 >> i;
      }
      else
      {
         ener1 = 0;
      }
      ener2 = ener2 >> 1;
   }

   tmp = ener1 - ener2;
   ener1 = (ener1 + ener2) + 1;

   tmp = (tmp << 15) / ener1;

   return (tmp);
}
