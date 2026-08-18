/*
 *===================================================================
 *  3GPP AMR Wideband Floating-point Speech Codec
 *===================================================================
 */
#include <math.h>
#include <memory.h>
#include "typedef.h"
#include "dec_main.h"
#include "dec_lpc.h"

#define MAX_16       (Word16)0x7FFF
#define MIN_16       (Word16)0x8000
#define L_SUBFR      64       /* Subframe size                    */
#define L_SUBFR16k   80       /* Subframe size at 16kHz           */
#define M16k         20       /* Order of LP filter               */
#define PREEMPH_FAC  22282    /* preemphasis factor (0.68 in Q15) */
#define FAC4         4
#define FAC5         5
#define UP_FAC       20480    /* 5/4 in Q14                       */
#define INV_FAC5     6554     /* 1/5 in Q15                       */
#define NB_COEF_UP   12
#define L_FIR        31
#define MODE_7k      0
#define MODE_24k     8


extern const Word16 D_ROM_pow2[];
extern const Word16 D_ROM_isqrt[];
extern const Word16 D_ROM_log2[];
extern const Word16 D_ROM_fir_up[];
extern const Word16 D_ROM_fir_6k_7k[];
extern const Word16 D_ROM_fir_7k[];
extern const Word16 D_ROM_hp_gain[];

#ifdef WIN32
#pragma warning( disable : 4310)
#endif

Word16 D_UTIL_random(Word16 *seed)
{
   *seed = (Word16)(*seed * 31821L + 13849L);
   return(*seed);
}

Word32 D_UTIL_pow2(Word16 exponant, Word16 fraction)
{
	Word32 L_x, tmp, i, exp;
	Word16 a;

	L_x = fraction * 32;          /* L_x = fraction<<6             */
	i = L_x >> 15;                /* Extract b10-b16 of fraction   */
	a = (Word16)(L_x);            /* Extract b0-b9   of fraction   */
	a = (Word16)(a & (Word16)0x7fff);
	L_x = D_ROM_pow2[i] << 16;    /* table[i] << 16                */
	tmp = D_ROM_pow2[i] - D_ROM_pow2[i + 1];  /* table[i] - table[i+1] */
	tmp = L_x - ((tmp * a) << 1); /* L_x -= tmp*a*2                */
	exp = 30 - exponant;
	if (exp <= 31)
	{
		L_x = tmp >> exp;

		if ((1 << (exp - 1)) & tmp)
		{
			L_x++;
		}
	}
	else
	{
		L_x = 0;
	}

	return(L_x);
}

Word16 D_UTIL_norm_l(Word32 L_var1)
{
   Word16 var_out;

   if(L_var1 == 0)
   {
      var_out = 0;
   }
   else
   {
      if(L_var1 == (Word32)0xffffffffL)
      {
         var_out = 31;
      }
      else
      {
         if(L_var1 < 0)
         {
            L_var1 = ~L_var1;
         }

         for(var_out = 0; L_var1 < (Word32)0x40000000L; var_out++)
         {
            L_var1 <<= 1;
         }
      }
   }

   return(var_out);
}

Word16 D_UTIL_norm_s(Word16 var1)
{
   Word16 var_out;

   if(var1 == 0)
   {
      var_out = 0;
   }
   else
   {
      if(var1 == -1)
      {
         var_out = 15;
      }
      else
      {
         if(var1 < 0)
         {
            var1 = (Word16)~var1;
         }

         for(var_out = 0; var1 < 0x4000; var_out++)
         {
            var1 <<= 1;
         }
      }
   }
   return(var_out);
}

Word32 D_UTIL_dot_product12(Word16 x[], Word16 y[], Word16 lg, Word16 *exp)
{
   Word32 sum, i, sft;

   sum = 0L;

   for(i = 0; i < lg; i++)
   {
      sum += x[i] * y[i];
   }
   sum = (sum << 1) + 1;

   sft = D_UTIL_norm_l(sum);
   sum = sum << sft;
   *exp = (Word16)(30 - sft);   /* exponent = 0..30 */

   return(sum);
}

void D_UTIL_normalised_inverse_sqrt(Word32 *frac, Word16 *exp)
{
   Word32 i, tmp;
   Word16 a;

   if(*frac <= (Word32)0)
   {
      *exp = 0;
      *frac = 0x7fffffffL;
      return;
   }

   if((*exp & 0x1) == 1)   /* If exponant odd -> shift right */
   {
      *frac = *frac >> 1;
   }
   *exp = (Word16)(-((*exp - 1) >> 1));
   *frac = *frac >> 9;
   i = *frac >>16;      /* Extract b25-b31   */
   *frac = *frac >> 1;
   a = (Word16)(*frac); /* Extract b10-b24   */
   a = (Word16)(a & (Word16)0x7fff);
   i = i - 16;
   *frac = D_ROM_isqrt[i] << 16; /* table[i] << 16    */
   tmp = D_ROM_isqrt[i] - D_ROM_isqrt[i + 1];   /* table[i] - table[i+1]) */
   *frac = *frac - ((tmp * a) << 1);   /* frac -=  tmp*a*2  */
}

Word32 D_UTIL_inverse_sqrt(Word32 L_x)
{
   Word32 L_y;
   Word16 exp;

   exp = D_UTIL_norm_l(L_x);
   L_x = (L_x << exp);   /* L_x is normalized */
   exp = (Word16)(31 - exp);
   D_UTIL_normalised_inverse_sqrt(&L_x, &exp);

   if(exp < 0)
   {
      L_y = (L_x >> -exp);   /* denormalization   */
   }
   else
   {
      L_y = (L_x << exp);   /* denormalization   */
   }

   return(L_y);
}

static void D_UTIL_normalised_log2(Word32 L_x, Word16 exp, Word16 *exponent,
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
   L_y = D_ROM_log2[i] << 16;               /* table[i] << 16        */
   tmp = D_ROM_log2[i] - D_ROM_log2[i + 1]; /* table[i] - table[i+1] */
   L_y = L_y - ((tmp * a) << 1);            /* L_y -= tmp*a*2        */
   *fraction = (Word16)(L_y >> 16);
}

void D_UTIL_log2(Word32 L_x, Word16 *exponent, Word16 *fraction)
{
   Word16 exp;

   exp = D_UTIL_norm_l(L_x);
   D_UTIL_normalised_log2((L_x <<exp), exp, exponent, fraction);
}

void D_UTIL_l_extract(Word32 L_32, Word16 *hi, Word16 *lo)
{
   *hi = (Word16)(L_32 >> 16);
   *lo = (Word16)((L_32 >> 1) - (*hi * 32768));
}

Word32 D_UTIL_mpy_32_16(Word16 hi, Word16 lo, Word16 n)
{
   Word32 L_32;

   L_32 = hi * n;
   L_32 += (lo * n) >> 15;

   return(L_32 << 1);
}

Word32 D_UTIL_mpy_32(Word16 hi1, Word16 lo1, Word16 hi2, Word16 lo2)
{
   Word32 L_32;

   L_32 = hi1 * hi2;
   L_32 += (hi1 * lo2) >> 15;
   L_32 += (lo1 * hi2) >> 15;

   return(L_32 << 1);
}

Word16 D_UTIL_saturate(Word32 inp)
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

void D_UTIL_signal_up_scale(Word16 x[], Word16 lg, Word16 exp)
{
    Word32 i, tmp;

    for (i = 0; i < lg; i++)
    {
       tmp = x[i] << exp;
       x[i] = D_UTIL_saturate(tmp);
    }
}

void D_UTIL_signal_down_scale(Word16 x[], Word16 lg, Word16 exp)
{
   Word32 i, tmp;

   for(i = 0; i < lg; i++)
   {
      tmp = x[i] << 16;
      tmp = tmp >> exp;
      x[i] = (Word16)((tmp + 0x8000) >> 16);
   }
}

static void D_UTIL_deemph_32(Word16 x_hi[], Word16 x_lo[], Word16 y[], Word16 mu, Word16 L, Word16 *mem)
{
   Word32 i, fac;
   Word32 tmp;

   fac = mu >> 1;
   tmp = (x_hi[0] << 12) + x_lo[0];
   tmp = (tmp << 6) + (*mem * fac);
   tmp = (tmp + 0x2000) >> 14;
   y[0] = D_UTIL_saturate(tmp);

   for(i = 1; i < L; i++)
   {
      tmp = (x_hi[i] << 12) + x_lo[i];
      tmp = (tmp << 6) + (y[i - 1] * fac);
      tmp = (tmp + 0x2000) >> 14;
      y[i] = D_UTIL_saturate(tmp);
   }

   *mem = y[L - 1];
}

static void D_UTIL_synthesis_32(Word16 a[], Word16 m, Word16 exc[],
                                Word16 Qnew, Word16 sig_hi[], Word16 sig_lo[],
                                Word16 lg)
{
   Word32 i, j, a0, s;
   Word32 tmp, tmp2;

   s = D_UTIL_norm_s((Word16)a[0]) - 2;
   a0 = a[0] >> (4 + Qnew);   /* input / 16 and >>Qnew */

   for(i = 0; i < lg; i++)
   {
      tmp = 0;

      for(j = 1; j <= m; j++)
      {
         tmp -= sig_lo[i - j] * a[j];
      }

      tmp = tmp >> (15 - 4);
      tmp2 = exc[i] * a0;

      for(j = 1; j <= m; j++)
      {
         tmp2 -= sig_hi[i - j] * a[j];
      }

      tmp += tmp2 << 1;
      tmp <<= s;

      sig_hi[i] = (Word16)(tmp >> 13);
      sig_lo[i] = (Word16)((tmp  >> 1) - (sig_hi[i] * 4096));
   }
}

static void D_UTIL_hp50_12k8(Word16 signal[], Word16 lg, Word16 mem[])
{
   Word32 i, L_tmp;
   Word16 y2_hi, y2_lo, y1_hi, y1_lo, x0, x1, x2;

   y2_hi = mem[0];
   y2_lo = mem[1];
   y1_hi = mem[2];
   y1_lo = mem[3];
   x0 = mem[4];
   x1 = mem[5];

   for(i = 0; i < lg; i++)
   {
      x2 = x1;
      x1 = x0;
      x0 = signal[i];
      L_tmp = 8192L;   /* rounding to maximise precision */
      L_tmp = L_tmp + (y1_lo * 16211);
      L_tmp = L_tmp + (y2_lo * (-8021));
      L_tmp = L_tmp >> 14;
      L_tmp = L_tmp + (y1_hi * 32422);
      L_tmp = L_tmp + (y2_hi * (-16042));
      L_tmp = L_tmp + (x0 * 8106);
      L_tmp = L_tmp + (x1 * (-16212));
      L_tmp = L_tmp + (x2 * 8106);
      L_tmp = L_tmp << 2;  /* coeff Q11 --> Q14 */
      y2_hi = y1_hi;
      y2_lo = y1_lo;
      D_UTIL_l_extract(L_tmp, &y1_hi, &y1_lo);
      L_tmp = (L_tmp + 0x4000) >> 15;   /* coeff Q14 --> Q15 with saturation */
      signal[i] = D_UTIL_saturate(L_tmp);

   }
   mem[0] = y2_hi;
   mem[1] = y2_lo;
   mem[2] = y1_hi;
   mem[3] = y1_lo;
   mem[4] = x0;
   mem[5] = x1;
}

Word16 D_UTIL_interpol(Word16 *x, Word16 const *fir, Word16 frac,
                       Word16 resol, Word16 nb_coef)
{
   Word32 i, k;
   Word32 sum;

   x = x - nb_coef + 1;
   sum = 0L;

   for(i = 0, k = ((resol - 1) - frac); i < 2 * nb_coef; i++,
      k = (Word16)(k + resol))
   {
      sum = sum + (x[i] * fir[k]);
   }

   if((sum < 536846336) & (sum > -536879104))
   {
      sum = (sum + 0x2000) >> 14;
   }
   else if(sum > 536846336)
   {
      sum = 32767;
   }
   else
   {
      sum = -32768;
   }

   return((Word16)sum);   /* saturation can occur here */
}

static void D_UTIL_up_samp(Word16 *sig_d, Word16 *sig_u, Word16 L_frame)
{
   Word32 pos, i, j;
   Word16 frac;

   pos = 0;   /* position with 1/5 resolution */

   for(j = 0; j < L_frame; j++)
   {
      i = (pos * INV_FAC5) >> 15;   /* integer part = pos * 1/5 */
      frac = (Word16)(pos - ((i << 2) + i));   /* frac = pos - (pos/5)*5   */
      sig_u[j] = D_UTIL_interpol(&sig_d[i], D_ROM_fir_up, frac, FAC5, NB_COEF_UP);
      pos = pos + FAC4;   /* position + 4/5 */
   }
}

static void D_UTIL_oversamp_16k(Word16 sig12k8[], Word16 lg, Word16 sig16k[], Word16 mem[])
{
   Word16 lg_up;
   Word16 signal[L_SUBFR + (2 * NB_COEF_UP)];

   memcpy(signal, mem, (2 * NB_COEF_UP) * sizeof(Word16));
   memcpy(signal + (2 * NB_COEF_UP), sig12k8, lg * sizeof(Word16));
   lg_up = (Word16)(((lg * UP_FAC) >> 15) << 1);
   D_UTIL_up_samp(signal + NB_COEF_UP, sig16k, lg_up);
   memcpy(mem, signal + lg, (2 * NB_COEF_UP) * sizeof(Word16));
}

void D_UTIL_hp400_12k8(Word16 signal[], Word16 lg, Word16 mem[])
{
   Word32 i, L_tmp;
   Word16 y2_hi, y2_lo, y1_hi, y1_lo, x0, x1, x2;

   y2_hi = mem[0];
   y2_lo = mem[1];
   y1_hi = mem[2];
   y1_lo = mem[3];
   x0 = mem[4];
   x1 = mem[5];

   for(i = 0; i < lg; i++)
   {
      x2 = x1;
      x1 = x0;
      x0 = signal[i];
      L_tmp = 8192L + (y1_lo * 29280);
      L_tmp = L_tmp + (y2_lo * (-14160));
      L_tmp = (L_tmp >> 14);
      L_tmp = L_tmp + (y1_hi * 58560);
      L_tmp = L_tmp + (y2_hi * (-28320));
      L_tmp = L_tmp + (x0 * 1830);
      L_tmp = L_tmp + (x1 * (-3660));
      L_tmp = L_tmp + (x2 * 1830);
      L_tmp = (L_tmp << 1);   /* coeff Q12 --> Q13 */
      y2_hi = y1_hi;
      y2_lo = y1_lo;
      D_UTIL_l_extract(L_tmp, &y1_hi, &y1_lo);

      signal[i] = (Word16)((L_tmp + 0x8000) >> 16);
   }
   mem[0] = y2_hi;
   mem[1] = y2_lo;
   mem[2] = y1_hi;
   mem[3] = y1_lo;
   mem[4] = x0;
   mem[5] = x1;
}

static void D_UTIL_synthesis(Word16 a[], Word16 m, Word16 x[], Word16 y[], Word16 lg, Word16 mem[], Word16 update)
{
   Word32 i, j, tmp, s;
   Word16 y_buf[L_SUBFR16k + M16k], a0;
   Word16 *yy;

   yy = &y_buf[m];

   s = D_UTIL_norm_s(a[0]) - 2;
   memcpy(y_buf, mem, m * sizeof(Word16));

   a0 = (Word16)(a[0] >> 1);   /* input / 2 */

   for(i = 0; i < lg; i++)
   {
      tmp = x[i] * a0;

      for(j = 1; j <= m; j++)
      {
         tmp -= a[j] * yy[i - j];
      }
      tmp <<= s;

      y[i] = yy[i] = (Word16)((tmp + 0x800) >> 12);
   }

   if(update)
   {
      memcpy(mem, &yy[lg - m], m * sizeof(Word16));
   }
}

void D_UTIL_bp_6k_7k(Word16 signal[], Word16 lg, Word16 mem[])
{
   Word32 x[L_SUBFR16k + (L_FIR - 1)];
   Word32 i, j, tmp;

   for(i = 0; i < (L_FIR - 1); i++)
   {
      x[i] = (Word16)mem[i];   /* gain of filter = 4 */
   }

   for(i = 0; i < lg; i++)
   {
      x[i + L_FIR - 1] = signal[i] >> 2;   /* gain of filter = 4 */
   }

   for(i = 0; i < lg; i++)
   {
      tmp = 0;

      for(j = 0; j < L_FIR; j++)
      {
         tmp += x[i + j] * D_ROM_fir_6k_7k[j];
      }

      signal[i] = (Word16)((tmp + 0x4000) >> 15);
   }

   for(i = 0; i < (L_FIR - 1); i++)
   {
      mem[i] = (Word16)x[lg + i];
   }
}

static void D_UTIL_hp_7k(Word16 signal[], Word16 lg, Word16 mem[])
{

   Word32 i, j, tmp;
   Word16 x[L_SUBFR16k + (L_FIR - 1)];

   memcpy(x, mem, (L_FIR - 1) * sizeof(Word16));
   memcpy(&x[L_FIR - 1], signal, lg * sizeof(Word16));

   for(i = 0; i < lg; i++)
   {
      tmp = 0;

      for(j = 0; j < L_FIR; j++)
      {
         tmp += x[i + j] * D_ROM_fir_7k[j];
      }

      signal[i] = (Word16)((tmp + 0x4000) >> 15);
   }

   memcpy(mem, x + lg, (L_FIR - 1) * sizeof(Word16));
}

void D_UTIL_dec_synthesis(Word16 Aq[], Word16 exc[], Word16 Q_new,
                          Word16 synth16k[], Word16 prms, Word16 HfIsf[],
                          Word16 mode, Word16 newDTXState, Word16 bfi,
                          Decoder_State *st)
{
   Word32 tmp, i;
   Word16 exp;
   Word16 ener, exp_ener;
   Word32 fac;
   Word16 synth_hi[M + L_SUBFR], synth_lo[M + L_SUBFR];
   Word16 synth[L_SUBFR];
   Word16 HF[L_SUBFR16k];   /* High Frequency vector      */
   Word16 Ap[M16k + 1];
   Word16 HfA[M16k + 1];
   Word16 HF_corr_gain;
   Word16 HF_gain_ind;
   Word32 gain1, gain2;
   Word16 weight1, weight2;

   memcpy(synth_hi, st->mem_syn_hi, M * sizeof(Word16));
   memcpy(synth_lo, st->mem_syn_lo, M * sizeof(Word16));
   D_UTIL_synthesis_32(Aq, M, exc, Q_new, synth_hi + M, synth_lo + M, L_SUBFR);
   memcpy(st->mem_syn_hi, synth_hi + L_SUBFR, M * sizeof(Word16));
   memcpy(st->mem_syn_lo, synth_lo + L_SUBFR, M * sizeof(Word16));
   D_UTIL_deemph_32(synth_hi + M, synth_lo + M, synth, PREEMPH_FAC, L_SUBFR,
      &(st->mem_deemph));
   D_UTIL_hp50_12k8(synth, L_SUBFR, st->mem_sig_out);
   D_UTIL_oversamp_16k(synth, L_SUBFR, synth16k, st->mem_oversamp);

   for(i = 0; i < L_SUBFR16k; i++)
   {
      HF[i] = (Word16)(D_UTIL_random(&(st->mem_seed2)) >> 3);
   }

   D_UTIL_signal_down_scale(exc, L_SUBFR, 3);
   Q_new = (Word16)(Q_new - 3);
   ener = (Word16)(D_UTIL_dot_product12(exc, exc, L_SUBFR, &exp_ener) >> 16);
   exp_ener = (Word16)(exp_ener - (Q_new << 1));

   tmp = (Word16)(D_UTIL_dot_product12(HF, HF, L_SUBFR16k, &exp) >> 16);

   if(tmp > ener)
   {
      tmp = tmp >> 1;   /* Be sure tmp < ener */
      exp = (Word16)(exp + 1);
   }

   tmp = (tmp << 15) / ener;

   if(tmp > 32767)
   {
      tmp = 32767;
   }

   tmp = tmp << 16;   /* result is normalized */
   exp = (Word16)(exp - exp_ener);
   D_UTIL_normalised_inverse_sqrt(&tmp, &exp);

   if(exp >= 0)
   {
      tmp = tmp >> (15 - exp);
   }
   else
   {
      tmp = tmp >> (-exp);
      tmp = tmp >> 15;
   }

   if(tmp > 0x7FFF)
   {
      tmp = 0x7FFF;
   }

   for(i = 0; i < L_SUBFR16k; i++)
   {
      HF[i] = (Word16)((HF[i] * tmp) >> 15);
   }

   D_UTIL_hp400_12k8(synth, L_SUBFR, st->mem_hp400);
   tmp = 0L;

   for(i = 0; i < L_SUBFR; i++)
   {
      tmp = tmp + (synth[i] * synth[i]);
   }

   tmp = (tmp << 1) + 1;
   exp = D_UTIL_norm_l(tmp);
   ener = (Word16)((tmp << exp) >> 16);   /* ener = r[0] */
   tmp = 0L;

   for(i = 1; i < L_SUBFR; i++)
   {
      tmp = tmp + (synth[i] * synth[i - 1]);
   }

   tmp = (tmp << 1) + 1;
   tmp = (tmp << exp) >> 16;   /* tmp = r[1] */

   if(tmp > 0)
   {
      fac = ((tmp << 15) / ener);

      if(fac > 32767)
      {
         fac = 32767;
      }
   }
   else
   {
      fac = 0;
   }

   gain1 = (32767 - fac);
   gain2 = ((32767 - fac) * 20480) >> 15;
   gain2 = (gain2 << 1);

   if(gain2 > 32767)
      gain2 = 32767;

   if(st->mem_vad_hist > 0)
   {
      weight1 = 0;
      weight2 = 32767;
   }
   else
   {
      weight1 = 32767;
      weight2 = 0;
   }

   tmp = (weight1 * gain1) >> 15;
   tmp = tmp + ((weight2 * gain2) >> 15);

   if(tmp != 0)
   {
      tmp = tmp + 1;
   }

   if(tmp < 3277)
   {
      tmp = 3277;
   }

   if((mode == MODE_24k) & (bfi == 0))
   {
      HF_gain_ind = prms;
      HF_corr_gain = D_ROM_hp_gain[HF_gain_ind];

      for(i = 0; i < L_SUBFR16k; i++)
      {
         HF[i] = (Word16)(((HF[i] * HF_corr_gain) >> 15) << 1);
      }
   }
   else
   {
      for(i = 0; i < L_SUBFR16k; i++)
      {
         HF[i] = (Word16)((HF[i] * tmp) >> 15);
      }
   }

   if((mode <= MODE_7k) & (newDTXState == SPEECH))
   {
      D_LPC_isf_extrapolation(HfIsf);
      D_LPC_isp_a_conversion(HfIsf, HfA, 0, M16k);
      D_LPC_a_weight(HfA, Ap, 29491, M16k);   /* fac=0.9 */
      D_UTIL_synthesis(Ap, M16k, HF, HF, L_SUBFR16k, st->mem_syn_hf, 1);
   }
   else
   {
      D_LPC_a_weight(Aq, Ap, 19661, M);   /* fac=0.6 */
      D_UTIL_synthesis(Ap, M, HF, HF, L_SUBFR16k, st->mem_syn_hf + (M16k - M), 1);
   }

   D_UTIL_bp_6k_7k(HF, L_SUBFR16k, st->mem_hf);

   if(mode == MODE_24k)
   {
      D_UTIL_hp_7k(HF, L_SUBFR16k, st->mem_hf3);
   }

   for(i = 0; i < L_SUBFR16k; i++)
   {
      tmp = (synth16k[i] + HF[i]);
      synth16k[i] = D_UTIL_saturate(tmp);
   }
}

void D_UTIL_preemph(Word16 x[], Word16 mu, Word16 lg, Word16 *mem)
{
   Word32 i, L_tmp;
   Word16 temp;

   temp = x[lg - 1];

   for(i = lg - 1; i > 0; i--)
   {
      L_tmp = x[i] << 15;
      L_tmp = L_tmp - (x[i - 1] * mu);
      x[i] = (Word16)((L_tmp + 0x4000) >> 15);
   }

   L_tmp = x[0] << 15;
   L_tmp = L_tmp - (*mem * mu);
   x[0] = (Word16)((L_tmp + 0x4000) >> 15);
   *mem = temp;
}
