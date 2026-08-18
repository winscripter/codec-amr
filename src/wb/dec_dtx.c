/*
 *===================================================================
 *  3GPP AMR Wideband Floating-point Speech Codec
 *===================================================================
 */
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include "typedef.h"
#include "dec_dtx.h"
#include "dec_lpc.h"
#include "dec_util.h"

#define MAX_31                      (Word32)0x3FFFFFFF
#define L_FRAME                     256   /* Frame size                          */
#define RX_SPEECH_LOST              2
#define RX_SPEECH_BAD               3
#define RX_SID_FIRST                4
#define RX_SID_UPDATE               5
#define RX_SID_BAD                  6
#define RX_NO_DATA                  7
#define ISF_GAP                     128   /* 50                                  */
#define D_DTX_MAX_EMPTY_THRESH      50
#define GAIN_FACTOR                 75
#define ISF_FACTOR_LOW              256
#define ISF_FACTOR_STEP             2
#define ISF_DITH_GAP                448
#define D_DTX_HANG_CONST            7     /* yields eight frames of SP HANGOVER  */
#define D_DTX_ELAPSED_FRAMES_THRESH (24 + 7 - 1)
#define RANDOM_INITSEED             21845 /* own random init value               */

int D_DTX_reset(D_DTX_State *st, const Word16 *isf_init)
{
   Word32 i;

   if(st == (D_DTX_State*)NULL)
   {
      return(-1);
   }

   st->mem_since_last_sid = 0;
   st->mem_true_sid_period_inv = (1 << 13);
   st->mem_log_en = 3500;
   st->mem_log_en_prev = 3500;
   st->mem_cng_seed = RANDOM_INITSEED;
   st->mem_hist_ptr = 0;

   /* Init isf_hist[] and decoder log frame energy */
   memcpy(st->mem_isf, isf_init, M * sizeof(Word16));
   memcpy(st->mem_isf_prev, isf_init, M * sizeof(Word16));

   for(i = 0; i < D_DTX_HIST_SIZE; i++)
   {
      memcpy(&st->mem_isf_buf[i * M], isf_init, M * sizeof(Word16));
      st->mem_log_en_buf[i] = 3500;
   }
   st->mem_dtx_hangover_count = D_DTX_HANG_CONST;
   st->mem_dec_ana_elapsed_count = 127;
   st->mem_sid_frame = 0;
   st->mem_valid_data = 0;
   st->mem_dtx_hangover_added = 0;
   st->mem_dtx_global_state = SPEECH;
   st->mem_data_updated = 0;
   st->mem_dither_seed = RANDOM_INITSEED;
   st->mem_cn_dith = 0;

   return(0);
}

int D_DTX_init(D_DTX_State **st, const Word16 *isf_init)
{
   D_DTX_State *s;

   if(st == (D_DTX_State**)NULL)
   {
      return(-1);
   }

   *st = NULL;

   if((s = (D_DTX_State*)malloc(sizeof(D_DTX_State))) == NULL)
   {
      return(-1);
   }

   D_DTX_reset(s, isf_init);
   *st = s;

   return(0);
}

void D_DTX_exit(D_DTX_State **st)
{
   if(st == NULL || *st == NULL)
   {
      return;
   }

   free(*st);
   *st = NULL;

   return;
}

UWord8 D_DTX_rx_handler(D_DTX_State *st, UWord8 frame_type)
{
   UWord8 newState;
   UWord8 encState;

   if((frame_type == RX_SID_FIRST) | (frame_type == RX_SID_UPDATE) |
      (frame_type == RX_SID_BAD) | (((st->mem_dtx_global_state == DTX) |
      (st->mem_dtx_global_state == D_DTX_MUTE)) & ((frame_type == RX_NO_DATA) |
      (frame_type == RX_SPEECH_BAD) | (frame_type == RX_SPEECH_LOST))))
   {
      newState = DTX;

      if((st->mem_dtx_global_state == D_DTX_MUTE) &
         ((frame_type == RX_SID_BAD) | (frame_type == RX_SID_FIRST) |
         (frame_type == RX_SPEECH_LOST) | (frame_type == RX_NO_DATA)))
      {
         newState = D_DTX_MUTE;
      }

      st->mem_since_last_sid = D_UTIL_saturate(st->mem_since_last_sid + 1);

      if(st->mem_since_last_sid > D_DTX_MAX_EMPTY_THRESH)
      {
         newState = D_DTX_MUTE;
      }
   }
   else
   {
      newState = SPEECH;
      st->mem_since_last_sid = 0;
   }

   if((st->mem_data_updated == 0) & (frame_type == RX_SID_UPDATE))
   {
      st->mem_dec_ana_elapsed_count = 0;
   }

   st->mem_dec_ana_elapsed_count++;

   if(st->mem_dec_ana_elapsed_count > 127)
   {
      st->mem_dec_ana_elapsed_count = 127;
   }

   st->mem_dtx_hangover_added = 0;

   if((frame_type == RX_SID_FIRST) | (frame_type == RX_SID_UPDATE) |
      (frame_type == RX_SID_BAD) | (frame_type == RX_NO_DATA))
   {
      encState = DTX;
   }
   else
   {
      encState = SPEECH;
   }

   if(encState == SPEECH)
   {
      st->mem_dtx_hangover_count = D_DTX_HANG_CONST;
   }
   else
   {
      if(st->mem_dec_ana_elapsed_count > D_DTX_ELAPSED_FRAMES_THRESH)
      {
         st->mem_dtx_hangover_added = 1;
         st->mem_dec_ana_elapsed_count = 0;
         st->mem_dtx_hangover_count = 0;
      }
      else if(st->mem_dtx_hangover_count == 0)
      {
         st->mem_dec_ana_elapsed_count = 0;
      }
      else
      {
         st->mem_dtx_hangover_count--;
      }
   }

   if(newState != SPEECH)
   {
      st->mem_sid_frame = 0;
      st->mem_valid_data = 0;

      if(frame_type == RX_SID_FIRST)
      {
         st->mem_sid_frame = 1;
      }
      else if(frame_type == RX_SID_UPDATE)
      {
         st->mem_sid_frame = 1;
         st->mem_valid_data = 1;
      }
      else if(frame_type == RX_SID_BAD)
      {
         st->mem_sid_frame = 1;
         st->mem_dtx_hangover_added = 0;
      }
   }

   return newState;
}

static void D_DTX_cn_dithering(Word16 isf[M], Word32 *L_log_en_int,
                               Word16 *dither_seed)
{
   Word32 temp, temp1, i, dither_fac, rand_dith,rand_dith2;

   rand_dith = D_UTIL_random(dither_seed) >> 1;
   rand_dith2 = D_UTIL_random(dither_seed) >>1;
   rand_dith = rand_dith + rand_dith2;
   *L_log_en_int = *L_log_en_int + ((rand_dith * GAIN_FACTOR) << 1);

   if(*L_log_en_int < 0)
   {
      *L_log_en_int = 0;
   }

   dither_fac = ISF_FACTOR_LOW;
   rand_dith = D_UTIL_random(dither_seed) >> 1;
   rand_dith2 = D_UTIL_random(dither_seed) >> 1;
   rand_dith = rand_dith + rand_dith2;
   temp = isf[0] + (((rand_dith * dither_fac) + 0x4000) >> 15);

   if(temp < ISF_GAP)
   {
      isf[0] = ISF_GAP;
   }
   else
   {
      isf[0] = (Word16)temp;
   }

   for(i = 1; i < M - 1; i++)
   {
      dither_fac = dither_fac + ISF_FACTOR_STEP;
      rand_dith = D_UTIL_random(dither_seed) >> 1;
      rand_dith2 = D_UTIL_random(dither_seed) >> 1;
      rand_dith = rand_dith + rand_dith2;
      temp = isf[i] + (((rand_dith * dither_fac) + 0x4000) >> 15);
      temp1 = temp - isf[i - 1];

      /* Make sure that isf spacing remains at least ISF_DITH_GAP Hz */
      if(temp1 < ISF_DITH_GAP)
      {
         isf[i] = (Word16)(isf[i - 1] + ISF_DITH_GAP);
      }
      else
      {
         isf[i] = (Word16)temp;
      }
   }

   if(isf[M - 2] > 16384)
   {
      isf[M - 2] = 16384;
   }
}

void D_DTX_exe(D_DTX_State *st, Word16 *exc2, Word16 new_state, Word16 isf[], Word16 **prms)
{

   Word32 i, j, L_tmp, ptr;
   Word32 exp0, int_fac;
   Word32 gain;
   Word32 L_isf[M], L_log_en_int, level32, ener32;
   Word16 log_en_index;
   Word16 tmp_int_length;
   Word16 exp, log_en_int_e, log_en_int_m, level;

   if((st->mem_dtx_hangover_added != 0) & (st->mem_sid_frame != 0))
   {
      ptr = st->mem_hist_ptr + 1;

      if(ptr == D_DTX_HIST_SIZE)
      {
         ptr = 0;
      }

      memcpy(&st->mem_isf_buf[ptr * M], &st->mem_isf_buf[st->mem_hist_ptr * M], M * sizeof(Word16));

      st->mem_log_en_buf[ptr] = st->mem_log_en_buf[st->mem_hist_ptr];

      st->mem_log_en = 0;
      memset(L_isf, 0, M * sizeof(Word32));

      for(i = 0; i < D_DTX_HIST_SIZE; i++)
      {
         st->mem_log_en = (Word16)(st->mem_log_en + st->mem_log_en_buf[i]);

         for(j = 0; j < M; j++)
         {
            L_isf[j] = L_isf[j] + st->mem_isf_buf[i * M + j];
         }
      }

      st->mem_log_en = (Word16)(st->mem_log_en >> 1);
      st->mem_log_en = (Word16)(st->mem_log_en + 1024);

      if(st->mem_log_en < 0)
      {
         st->mem_log_en = 0;
      }

      for(j = 0; j < M; j++)
      {
         st->mem_isf[j] = (Word16)(L_isf[j]>>3);
      }
   }

   if(st->mem_sid_frame != 0)
   {
      memcpy(st->mem_isf_prev, st->mem_isf, M * sizeof(Word16));
      st->mem_log_en_prev = st->mem_log_en;

      if(st->mem_valid_data != 0) /* new data available (no CRC) */
      {
         tmp_int_length = st->mem_since_last_sid;

         if(tmp_int_length > 32)
         {
            tmp_int_length = 32;
         }

         if(tmp_int_length >= 2)
         {
            st->mem_true_sid_period_inv =
               (Word16)(0x2000000 / (tmp_int_length << 10));
         }
         else
         {
            st->mem_true_sid_period_inv = 1 << 14;   /* 0.5 it Q15 */
         }

         D_LPC_isf_noise_d(*prms, st->mem_isf);
         (*prms) += 5;
         log_en_index = *(*prms)++;

         st->mem_cn_dith = *(*prms)++;
         st->mem_log_en = (Word16)(log_en_index << (15 - 6));
         st->mem_log_en = (Word16)((st->mem_log_en * 12483) >> 15);

         if((st->mem_data_updated == 0) ||
            (st->mem_dtx_global_state == SPEECH))
         {
            memcpy(st->mem_isf_prev, st->mem_isf, M * sizeof(Word16));
            st->mem_log_en_prev = st->mem_log_en;
         }
      }
   }

   if((st->mem_sid_frame != 0) && (st->mem_valid_data != 0))
   {
      st->mem_since_last_sid = 0;
   }

   if(st->mem_since_last_sid < 32)
   {
      int_fac = st->mem_since_last_sid << 10;   /* Q10 */
   }
   else
   {
      int_fac = 32767;
   }

   int_fac = (int_fac * st->mem_true_sid_period_inv) >> 15;

   if(int_fac > 1024)
   {
      int_fac = 1024;
   }
   int_fac = int_fac << 4;   /* Q10 -> Q14 */
   L_log_en_int = (int_fac * st->mem_log_en) << 1;   /* Q14 * Q9 -> Q24 */

   for(i = 0; i < M; i++)
   {
      isf[i] = (Word16)((int_fac * st->mem_isf[i]) >> 15);
   }
   int_fac = 16384 - int_fac;
   L_log_en_int = L_log_en_int + ((int_fac * st->mem_log_en_prev) << 1);

   for(i = 0; i < M; i++)
   {
      L_tmp = isf[i] + ((int_fac * st->mem_isf_prev[i]) >> 15);
      isf[i] = (Word16)(L_tmp << 1);   /* Q14 -> Q15 */
   }

   if(st->mem_cn_dith != 0)
   {
      D_DTX_cn_dithering(isf, &L_log_en_int, &st->mem_dither_seed);
   }

   L_log_en_int = (L_log_en_int >> 9); /* Q25 -> Q16 */
   log_en_int_e = (Word16)((L_log_en_int)>>16);
   log_en_int_m = (Word16)((L_log_en_int - (log_en_int_e << 16)) >> 1);
   log_en_int_e = (Word16)(log_en_int_e + (16 - 1));

   level32 = D_UTIL_pow2(log_en_int_e, log_en_int_m);   /* Q16 */
   exp0 = D_UTIL_norm_l(level32);
   level32 = (level32 << exp0);
   exp0 = (15 - exp0);
   level = (Word16)(level32 >> 16);

   for(i = 0; i < L_FRAME; i++)
   {
      exc2[i] = (Word16)((D_UTIL_random(&(st->mem_cng_seed)) >> 4));
   }

   ener32 = D_UTIL_dot_product12(exc2, exc2, L_FRAME, &exp);
   D_UTIL_normalised_inverse_sqrt(&ener32, &exp);
   gain = ener32 >>16;
   gain = (level * gain) >> 15;
   exp = (Word16)(exp0 + exp  + 4);

   if(exp >= 0)
   {
      for(i = 0; i < L_FRAME; i++)
      {
         L_tmp = (exc2[i] * gain) >> 15;   /* Q0 * Q15 */
         exc2[i] = (Word16)(L_tmp << exp);
      }
   }
   else
   {
      exp = (Word16)-exp;

      for(i = 0; i < L_FRAME; i++)
      {
         L_tmp = (exc2[i] * gain) >> 15;   /* Q0 * Q15 */
         exc2[i] = (Word16)(L_tmp >> exp);
      }
   }

   if(new_state == D_DTX_MUTE)
   {
      tmp_int_length = st->mem_since_last_sid;

      if(tmp_int_length > 32)
      {
         tmp_int_length = 32;
      }

      st->mem_true_sid_period_inv = D_UTIL_saturate((0x02000000 / (tmp_int_length << 10)));
      st->mem_since_last_sid = 0;
      st->mem_log_en_prev = st->mem_log_en;
      st->mem_log_en = D_UTIL_saturate(st->mem_log_en - 64);
   }

   if((st->mem_sid_frame != 0) && ((st->mem_valid_data != 0) ||
      ((st->mem_valid_data == 0) && (st->mem_dtx_hangover_added) != 0)))
   {
      st->mem_since_last_sid = 0;
      st->mem_data_updated = 1;
   }
}

void D_DTX_activity_update(D_DTX_State *st, Word16 isf[], Word16 exc[])
{

   Word32 L_frame_en, log_en;
   Word32 i;
   Word16 log_en_e, log_en_m;

   st->mem_hist_ptr = (Word16)(st->mem_hist_ptr + 1);

   if(st->mem_hist_ptr == D_DTX_HIST_SIZE)
   {
      st->mem_hist_ptr = 0;
   }

   memcpy(&st->mem_isf_buf[st->mem_hist_ptr * M], isf, M * sizeof(Word16));

   L_frame_en = 0;

   for(i = 0; i < L_FRAME; i++)
   {
      L_frame_en = L_frame_en + (exc[i] * exc[i]);
      if (L_frame_en > MAX_31)
      {
         L_frame_en = MAX_31;
         break;
      }
   }

   D_UTIL_log2(L_frame_en, &log_en_e, &log_en_m);
   log_en = log_en_e << 7;   /* Q7 */
   log_en = log_en + (log_en_m >> (15 - 7));
   log_en = log_en - 1024;
   st->mem_log_en_buf[st->mem_hist_ptr] = (Word16)log_en;
}
