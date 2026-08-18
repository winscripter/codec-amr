/*
 *===================================================================
 *  3GPP AMR Wideband Floating-point Speech Codec
 *===================================================================
 */
#include <math.h>
#include <memory.h>
#include "typedef.h"
#include "enc_main.h"
#include "enc_lpc.h"

#ifdef WIN32
#pragma warning( disable : 4310)
#endif

#define MAX_16 (Word16)0x7FFF
#define MIN_16 (Word16)0x8000
#define MAX_31 (Word32)0x3FFFFFFF
#define MIN_31 (Word32)0xC0000000
#define L_FRAME16k   320     /* Frame size at 16kHz         */
#define L_SUBFR16k   80      /* Subframe size at 16kHz      */
#define L_SUBFR      64      /* Subframe size               */
#define M16k         20      /* Order of LP filter          */
#define L_WINDOW     384     /* window size in LP analysis  */
#define PREEMPH_FAC  0.68F   /* preemphasis factor          */

extern const Word16 E_ROM_pow2[];
extern const Word16 E_ROM_log2[];
extern const Word16 E_ROM_isqrt[];
extern const Float32 E_ROM_fir_6k_7k[];
extern const Float32 E_ROM_hp_gain[];
extern const Float32 E_ROM_fir_ipol[];
extern const Float32 E_ROM_hamming_cos[];

Word16 E_UTIL_random(Word16 *seed)
{
  *seed = (Word16) (*seed * 31821L + 13849L);

  return(*seed);
}

Word16 E_UTIL_saturate(Word32 inp)
{
   Word16 out;

   if ((inp < MAX_16) & (inp > MIN_16))
   {
      out = (Word16)inp;
   }
   else
   {
      if (inp > 0)
      {
         out = MAX_16;
      }
      else
      {
         out = MIN_16;
      }
   }

   return(out);
}

Word32 E_UTIL_saturate_31(Word32 inp)
{
   Word32 out;

   if ((inp < MAX_31) & (inp > MIN_31))
   {
      out = inp;
   }
   else
   {
      if (inp > 0)
      {
         out = MAX_31;
      }
      else
      {
         out = MIN_31;
      }
   }

   return(out);
}

Word16 E_UTIL_norm_s (Word16 var1)
{
   Word16 var_out;

   if (var1 == 0)
   {
      var_out = 0;
   }
   else
   {
      if (var1 == -1)
      {
         var_out = 15;
      }
      else
      {
         if (var1 < 0)
         {
            var1 = (Word16)~var1;
         }
         for (var_out = 0; var1 < 0x4000; var_out++)
         {
            var1 <<= 1;
         }
      }
   }

   return (var_out);
}

Word16 E_UTIL_norm_l (Word32 L_var1)
{
   Word16 var_out;

   if (L_var1 == 0)
   {
      var_out = 0;
   }
   else
   {
      if (L_var1 == (Word32) 0xffffffffL)
      {
         var_out = 31;
      }
      else
      {
         if (L_var1 < 0)
         {
            L_var1 = ~L_var1;
         }
         for (var_out = 0; L_var1 < (Word32) 0x40000000L; var_out++)
         {
            L_var1 <<= 1;
         }
      }
   }

   return (var_out);
}

void E_UTIL_l_extract(Word32 L_32, Word16 *hi, Word16 *lo)
{
   *hi = (Word16)(L_32 >> 16);
   *lo = (Word16)((L_32 >> 1) - ((*hi * 16384) << 1));
   return;
}

Word32 E_UTIL_mpy_32_16 (Word16 hi, Word16 lo, Word16 n)
{
   Word32 L_32;

   L_32 = (hi * n) << 1;
   L_32 = L_32 + (((lo * n) >> 15) << 1);

   return (L_32);
}

Word32 E_UTIL_pow2(Word16 exponant, Word16 fraction)
{
   Word32 L_x, tmp, i, exp;
   Word16 a;

   L_x = fraction * 32;          /* L_x = fraction<<6             */
   i = L_x >> 15;                /* Extract b10-b16 of fraction   */
   a = (Word16)(L_x);            /* Extract b0-b9   of fraction   */
   a = (Word16)(a & (Word16)0x7fff);
   L_x = E_ROM_pow2[i] << 16;    /* table[i] << 16                */
   tmp = E_ROM_pow2[i] - E_ROM_pow2[i + 1];  /* table[i] - table[i+1] */
   L_x = L_x - ((tmp * a) << 1); /* L_x -= tmp*a*2                */
   exp = 30 - exponant;
   L_x = (L_x + (1 << (exp - 1))) >> exp;

   return(L_x);
}

static void E_UTIL_normalised_log2(Word32 L_x, Word16 exp, Word16 *exponent,
                          Word16 *fraction)
{
   Word32 i, a, tmp;
   Word32 L_y;

   if (L_x <= 0)
   {
      *exponent = 0;
      *fraction = 0;
      return;
   }

   *exponent = (Word16)(30 - exp);

   L_x = L_x >> 10;
   i = L_x >> 15;         /* Extract b25-b31               */
   a = L_x;               /* Extract b10-b24 of fraction   */
   a = a & 0x00007fff;
   i = i - 32;
   L_y = E_ROM_log2[i] << 16;               /* table[i] << 16        */
   tmp = E_ROM_log2[i] - E_ROM_log2[i + 1]; /* table[i] - table[i+1] */
   L_y = L_y - ((tmp * a) << 1);            /* L_y -= tmp*a*2        */
   *fraction = (Word16)(L_y >> 16);
}

void E_UTIL_log2_32 (Word32 L_x, Word16 *exponent, Word16 *fraction)
{
   Word16 exp;

   exp = E_UTIL_norm_l(L_x);
   E_UTIL_normalised_log2((L_x << exp), exp, exponent, fraction);
}

static Float32 E_UTIL_interpol(Float32 *x, Word32 frac, Word32 up_samp, Word32 nb_coef)
{
   Word32 i;
   Float32 s;
   Float32 *x1, *x2;
   const Float32 *c1, *c2;

   x1 = &x[0];
   x2 = &x[1];
   c1 = &E_ROM_fir_ipol[frac];
   c2 = &E_ROM_fir_ipol[up_samp - frac];

   s = 0.0;

   for(i = 0; i < nb_coef; i++)
   {
      s += x1[-i] * c1[up_samp * i] + x2[i] * c2[up_samp * i];
   }

   return s;
}

void E_UTIL_hp50_12k8(Float32 signal[], Word32 lg, Float32 mem[])
{
   Word32 i;
   Float32 x0, x1, x2, y0, y1, y2;

   y1 = mem[0];
   y2 = mem[1];
   x0 = mem[2];
   x1 = mem[3];

   for(i = 0; i < lg; i++)
   {
      x2 = x1;
      x1 = x0;
      x0 = signal[i];

      y0 = y1 * 1.978881836F + y2 * -0.979125977F + x0 * 0.989501953F + x1 * -1.979003906F + x2 * 0.989501953F;

      signal[i] = y0;
      y2 = y1;
      y1 = y0;
   }

   mem[0] = ((y1 > 1e-10) | (y1 < -1e-10)) ? y1 : 0;
   mem[1] = ((y2 > 1e-10) | (y2 < -1e-10)) ? y2 : 0;
   mem[2] = ((x0 > 1e-10) | (x0 < -1e-10)) ? x0 : 0;
   mem[3] = ((x1 > 1e-10) | (x1 < -1e-10)) ? x1 : 0;
}

static void E_UTIL_hp400_12k8(Float32 signal[], Word32 lg, Float32 mem[])
{
   Word32 i;

   Float32 x0, x1, x2;
   Float32 y0, y1, y2;

   y1 = mem[0];
   y2 = mem[1];
   x0 = mem[2];
   x1 = mem[3];

   for(i = 0; i < lg; i++)
   {
      x2 = x1;
      x1 = x0;
      x0 = signal[i];
      y0 = y1 * 1.787109375F + y2 * -0.864257812F + x0 * 0.893554687F + x1 * - 1.787109375F + x2 * 0.893554687F;
      signal[i] = y0;
      y2 = y1;
      y1 = y0;
   }

   mem[0] = y1;
   mem[1] = y2;
   mem[2] = x0;
   mem[3] = x1;
}

static void E_UTIL_bp_6k_7k(Float32 signal[], Word32 lg, Float32 mem[])
{
   Float32 x[L_SUBFR16k + 30];
   Float32 s0, s1, s2, s3;
   Float32 *px;
   Word32 i, j;

   memcpy(x, mem, 30 * sizeof(Float32));
   memcpy(x + 30, signal, lg * sizeof(Float32));

   px = x;

   for(i = 0; i < lg; i++)
   {
      s0 = 0;
      s1 = px[0] * E_ROM_fir_6k_7k[0];
      s2 = px[1] * E_ROM_fir_6k_7k[1];
      s3 = px[2] * E_ROM_fir_6k_7k[2];

      for(j = 3; j < 31; j += 4)
      {
         s0 += px[j] * E_ROM_fir_6k_7k[j];
         s1 += px[j + 1] * E_ROM_fir_6k_7k[j + 1];
         s2 += px[j + 2] * E_ROM_fir_6k_7k[j + 2];
         s3 += px[j + 3] * E_ROM_fir_6k_7k[j + 3];
      }

      px++;

      signal[i] = (Float32)((s0 + s1 + s2 + s3) * 0.25F);   /* gain of coef = 4.0 */
   }

   memcpy(mem, x + lg, 30 * sizeof(Float32));
}

void E_UTIL_preemph(Word16 x[], Word16 mu, Word32 lg, Word16 *mem)
{
   Word32 i, L_tmp;
   Word16 temp;

   temp = x[lg - 1];

   for (i = lg - 1; i > 0; i--)
   {
      L_tmp = x[i] << 15;
      L_tmp -= x[i - 1] * mu;
      x[i] = (Word16)((L_tmp + 0x4000) >> 15);
   }

   L_tmp = (x[0] << 15);
   L_tmp -= *mem * mu;
   x[0] = (Word16)((L_tmp + 0x4000) >> 15);

   *mem = temp;
}

void E_UTIL_f_preemph(Float32 *signal, Float32 mu, Word32 L, Float32 *mem)
{
   Word32 i;
   Float32 temp;

   temp = signal[L - 1];

   for (i = L - 1; i > 0; i--)
   {
      signal[i] = signal[i] - mu * signal[i - 1];
   }

   signal[0] -= mu * (*mem);
   *mem = temp;
}

void E_UTIL_deemph(Float32 *signal, Float32 mu, Word32 L, Float32 *mem)
{
   Word32 i;

   signal[0] = signal[0] + mu * (*mem);

   for (i = 1; i < L; i++)
   {
      signal[i] = signal[i] + mu * signal[i - 1];
   }

   *mem = signal[L - 1];

   if ((*mem < 1e-10) & (*mem > -1e-10))
   {
      *mem = 0;
   }
}

void E_UTIL_synthesis(Float32 a[], Float32 x[], Float32 y[], Word32 l,
                      Float32 mem[], Word32 update_m)
{

   Float32 buf[L_FRAME16k + M16k];     /* temporary synthesis buffer */
   Float32 s;
   Float32 *yy;
   Word32 i, j;

   memcpy(buf, mem, M * sizeof(Float32));
   yy = &buf[M];

   for (i = 0; i < l; i++)
   {
      s = x[i];

      for (j = 1; j <= M; j += 4)
      {
         s -= a[j] * yy[i - j];
         s -= a[j + 1] * yy[i - (j + 1)];
         s -= a[j + 2] * yy[i - (j + 2)];
         s -= a[j + 3] * yy[i - (j + 3)];
      }

      yy[i] = s;
      y[i] = s;
   }

   if (update_m)
   {
      memcpy(mem, &yy[l - M], M * sizeof(Float32));
   }
}

static void E_UTIL_down_samp(Float32 *res, Float32 *res_d, Word32 L_frame_d)
{
   Word32 i, j, frac;
   Float32 pos, fac;

   fac = 0.8F;
   pos = 0;

   for(i = 0; i < L_frame_d; i++)
   {
      j = (Word32)pos;    /* j = (Word32)( (Float32)i * inc); */
      frac = (Word32)(((pos - (Float32)j) * 4) + 0.5);
      res_d[i] = fac * E_UTIL_interpol(&res[j], frac, 4, 15);
      pos += 1.25F;
   }
}

void E_UTIL_decim_12k8(Float32 sig16k[], Word32 lg, Float32 sig12k8[], Float32 mem[])
{
   Float32 signal[(2 * 15) + L_FRAME16k];

   memcpy(signal, mem, 2 * 15 * sizeof(Float32));
   memcpy(&signal[2 * 15], sig16k, lg * sizeof(Float32));
   E_UTIL_down_samp(signal + 15, sig12k8, lg * 4 / 5);
   memcpy(mem, &signal[lg], 2 * 15 * sizeof(Float32));
}

void E_UTIL_residu(Float32 *a, Float32 *x, Float32 *y, Word32 l)
{
   Float32 s;
   Word32 i;

   for (i = 0; i < l; i++)
   {
      s = x[i];
      s += a[1] * x[i - 1];
      s += a[2] * x[i - 2];
      s += a[3] * x[i - 3];
      s += a[4] * x[i - 4];
      s += a[5] * x[i - 5];
      s += a[6] * x[i - 6];
      s += a[7] * x[i - 7];
      s += a[8] * x[i - 8];
      s += a[9] * x[i - 9];
      s += a[10] * x[i - 10];
      s += a[11] * x[i - 11];
      s += a[12] * x[i - 12];
      s += a[13] * x[i - 13];
      s += a[14] * x[i - 14];
      s += a[15] * x[i - 15];
      s += a[16] * x[i - 16];
      y[i] = s;
   }
}

void E_UTIL_convolve(Word16 x[], Word16 q, Float32 h[], Float32 y[])
{
   Float32 fx[L_SUBFR];
   Float32 temp, scale;
   Word32 i, n;

   scale = (Float32)pow(2, -q);

   for (i = 0; i < L_SUBFR; i++)
   {
      fx[i] = (Float32)(scale * x[i]);
   }

   for (n = 0; n < L_SUBFR; n += 2)
   {
      temp = 0.0;
      for (i = 0; i <= n; i++)
      {
         temp += (Float32)(fx[i] * h[n - i]);
      }
      y[n] = temp;

      temp = 0.0;
      for (i = 0; i <= (n + 1); i += 2)
      {
         temp += (Float32)(fx[i] * h[(n + 1) - i]);
         temp += (Float32)(fx[i + 1] * h[n - i]);
      }
      y[n + 1] = temp;

   }
}

void E_UTIL_f_convolve(Float32 x[], Float32 h[], Float32 y[])
{
   Float32 temp;
   Word32 i, n;

   for (n = 0; n < L_SUBFR; n += 2)
   {
      temp = 0.0;

      for (i = 0; i <= n; i++)
      {
         temp += x[i] * h[n - i];
      }

      y[n] = temp;

      temp = 0.0;

      for (i = 0; i <= (n + 1); i += 2)
      {
         temp += x[i] * h[(n + 1) - i];
         temp += x[i + 1] * h[n - i];
      }

      y[n + 1] = temp;
   }
}

void E_UTIL_signal_up_scale(Word16 x[], Word16 exp)
{
   Word32 i;
   Word32  tmp;

   for (i = 0; i < (PIT_MAX + L_INTERPOL + L_SUBFR); i++)
   {
      tmp = x[i] << exp;
      x[i] = E_UTIL_saturate(tmp);
   }
}

void E_UTIL_signal_down_scale(Word16 x[], Word32 lg, Word16 exp)
{
   Word32 i, tmp;

   for (i = 0; i < lg; i++)
   {
      tmp = x[i] << 16;
      tmp = tmp >> exp;
      x[i] = (Word16)((tmp + 0x8000) >> 16);
   }
}

Word32 E_UTIL_dot_product12(Word16 x[], Word16 y[], Word32 lg, Word32 *exp)
{
   Word32 i, sft, L_sum, L_sum1, L_sum2, L_sum3, L_sum4;

   L_sum1 = 0L;
   L_sum2 = 0L;
   L_sum3 = 0L;
   L_sum4 = 0L;

   for (i = 0; i < lg; i += 4)
   {
      L_sum1 += x[i] * y[i];
      L_sum2 += x[i + 1] * y[i + 1];
      L_sum3 += x[i + 2] * y[i + 2];
      L_sum4 += x[i + 3] * y[i + 3];
   }

   L_sum1 = E_UTIL_saturate_31(L_sum1);
   L_sum2 = E_UTIL_saturate_31(L_sum2);
   L_sum3 = E_UTIL_saturate_31(L_sum3);
   L_sum4 = E_UTIL_saturate_31(L_sum4);
   L_sum1 += L_sum3;
   L_sum2 += L_sum4;
   L_sum1 = E_UTIL_saturate_31(L_sum1);
   L_sum2 = E_UTIL_saturate_31(L_sum2);
   L_sum = L_sum1 + L_sum2;
   L_sum = (E_UTIL_saturate_31(L_sum) << 1) + 1;

   sft = E_UTIL_norm_l(L_sum);
   L_sum = (L_sum << sft);

   *exp = (30 - sft);

   return (L_sum);
}

void E_UTIL_normalised_inverse_sqrt(Word32 *frac, Word16 *exp)
{
   Word32 i, a, tmp;

   if (*frac <= (Word32) 0)
   {
      *exp = 0;
      *frac = 0x7fffffffL;
      return;
   }

   if ((Word16) (*exp & 1) == 1)  /* If exponant odd -> shift right */
   {
      *frac = (*frac >> 1);
   }

   *exp = (Word16)(-((*exp - 1) >> 1));

   *frac = (*frac >> 9);
   i = *frac >> 16;                    /* Extract b25-b31 */
   *frac = (*frac >> 1);
   a = (Word16)*frac;                  /* Extract b10-b24 */
   a = a & 0x00007fff;

   i = i - 16;

   *frac = E_ROM_isqrt[i] << 16;                /* table[i] << 16         */
   tmp = E_ROM_isqrt[i] - E_ROM_isqrt[i + 1];   /* table[i] - table[i+1]) */

   *frac = *frac - ((tmp * a) << 1);            /* frac -=  tmp*a*2       */
}

Word32 E_UTIL_enc_synthesis(Float32 Aq[], Float32 exc[], Float32 synth16k[], Coder_State *st)
{
   Float32 synth[L_SUBFR];
   Float32 HF[L_SUBFR16k];   /* High Frequency vector      */
   Float32 Ap[M + 1];
   Float32 HF_SP[L_SUBFR16k];   /* High Frequency vector (from original signal) */
   Float32 HP_est_gain, HP_calc_gain, HP_corr_gain, fac, tmp, ener, dist_min;
   Float32 dist, gain2;
   Word32 i, hp_gain_ind = 0;

   E_UTIL_synthesis(Aq, exc, synth, L_SUBFR, st->mem_syn2, 1);
   E_UTIL_deemph(synth, PREEMPH_FAC, L_SUBFR, &(st->mem_deemph));
   E_UTIL_hp50_12k8(synth, L_SUBFR, st->mem_sig_out);

   memcpy(HF_SP, synth16k, L_SUBFR16k * sizeof(Float32));

   for(i = 0; i < L_SUBFR16k; i++)
   {
      HF[i] = (Float32)E_UTIL_random(&(st->mem_seed));
   }

   ener = 0.01F;
   tmp = 0.01F;

   for(i = 0; i < L_SUBFR; i++)
   {
      ener += exc[i] * exc[i];
   }

   for(i = 0; i < L_SUBFR16k; i++)
   {
      tmp += HF[i] * HF[i];
   }

   tmp = (Float32)(sqrt(ener / tmp));

   for(i = 0; i < L_SUBFR16k; i++)
   {
      HF[i] *= tmp;
   }

   E_UTIL_hp400_12k8(synth, L_SUBFR, st->mem_hp400);
   ener = 0.001f;
   tmp = 0.001f;

   for(i = 1; i < L_SUBFR; i++)
   {
      ener += synth[i] * synth[i];
      tmp += synth[i] * synth[i - 1];
   }

   fac = tmp / ener;

   HP_est_gain = 1.0F - fac;
   gain2 = (1.0F - fac) * 1.25F;

   if(st->mem_vad_hist)
   {
      HP_est_gain = gain2;
   }

   if(HP_est_gain < 0.1)
   {
      HP_est_gain = 0.1f;
   }

   if(HP_est_gain > 1.0)
   {
      HP_est_gain = 1.0f;
   }

   E_LPC_a_weight(Aq, Ap, 0.6f, M);
   E_UTIL_synthesis(Ap, HF, HF, L_SUBFR16k, st->mem_syn_hf, 1);
   E_UTIL_bp_6k_7k(HF, L_SUBFR16k, st->mem_hf);
   E_UTIL_bp_6k_7k(HF_SP, L_SUBFR16k, st->mem_hf2);

   ener = 0.001F;
   tmp = 0.001F;

   for(i = 0; i < L_SUBFR16k; i++)
   {
      ener += HF_SP[i] * HF_SP[i];
      tmp += HF[i] * HF[i];
   }

   HP_calc_gain = (Float32)sqrt(ener /tmp);
   st->mem_gain_alpha *= st->dtx_encSt->mem_dtx_hangover_count / 7;

   if(st->dtx_encSt->mem_dtx_hangover_count > 6)
   {
      st->mem_gain_alpha = 1.0F;
   }

   HP_corr_gain = (HP_calc_gain * st->mem_gain_alpha) + ((1.0F - st->mem_gain_alpha) * HP_est_gain);

   dist_min = 100000.0F;

   for(i = 0; i < 16; i++)
   {
      dist = (HP_corr_gain - E_ROM_hp_gain[i]) * (HP_corr_gain - E_ROM_hp_gain[i]);

      if(dist_min > dist)
      {
         dist_min = dist;
         hp_gain_ind = i;
      }
   }

   HP_corr_gain = (Float32)E_ROM_hp_gain[hp_gain_ind];

   return(hp_gain_ind);
}

void E_UTIL_autocorr(Float32 *x, Float32 *r)
{
   Float32 t[L_WINDOW + M];
   Word32 i, j;

   for (i = 0; i < L_WINDOW; i += 4)
   {
      t[i] = x[i] * E_ROM_hamming_cos[i];
      t[i + 1] = x[i + 1] * E_ROM_hamming_cos[i + 1];
      t[i + 2] = x[i + 2] * E_ROM_hamming_cos[i + 2];
      t[i + 3] = x[i + 3] * E_ROM_hamming_cos[i + 3];
   }

   memset(&t[L_WINDOW], 0, M * sizeof(Float32));
   memset(r, 0, (M + 1) * sizeof(Float32));

   for (j = 0; j < L_WINDOW; j++)
   {
      r[0] += t[j] * t[j];
      r[1] += t[j] * t[j + 1];
      r[2] += t[j] * t[j + 2];
      r[3] += t[j] * t[j + 3];
      r[4] += t[j] * t[j + 4];
      r[5] += t[j] * t[j + 5];
      r[6] += t[j] * t[j + 6];
      r[7] += t[j] * t[j + 7];
      r[8] += t[j] * t[j + 8];
      r[9] += t[j] * t[j + 9];
      r[10] += t[j] * t[j + 10];
      r[11] += t[j] * t[j + 11];
      r[12] += t[j] * t[j + 12];
      r[13] += t[j] * t[j + 13];
      r[14] += t[j] * t[j + 14];
      r[15] += t[j] * t[j + 15];
      r[16] += t[j] * t[j + 16];
   }

   if (r[0] < 1.0F)
   {
      r[0] = 1.0F;
   }
}
