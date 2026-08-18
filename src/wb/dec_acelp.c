/*
 *===================================================================
 *  3GPP AMR Wideband Floating-point Speech Codec
 *===================================================================
 */
#include <memory.h>
#include "typedef.h"
#include "dec_util.h"

#define L_SUBFR      64    /* Subframe size              */
#define PRED_ORDER   4
#define MEAN_ENER    30    /* average innovation energy  */
extern const Word16 D_ROM_ph_imp_low[];
extern const Word16 D_ROM_ph_imp_mid[];

static void D_ACELP_add_pulse(Word32 pos[], Word32 nb_pulse,
                              Word32 track, Word16 code[])
{
   Word32 i, k;

   for(k = 0; k < nb_pulse; k++) {
      i = ((pos[k] & (16 - 1)) << 2) + track;

      if((pos[k] & 16) == 0)
      {
         code[i] = (Word16)(code[i] + 512);
      }
      else
      {
         code[i] = (Word16)(code[i] - 512);
      }
   }
}

static void D_ACELP_decode_1p_N1(Word32 index, Word32 N,
                                 Word32 offset, Word32 pos[])
{
   Word32 i, pos1, mask;

   mask = ((1 << N) - 1);

   pos1 = ((index & mask) + offset);
   i = ((index >> N) & 1);

   if(i == 1)
   {
      pos1 += 16;
   }

   pos[0] = pos1;
}

static void D_ACELP_decode_2p_2N1(Word32 index, Word32 N,
                                  Word32 offset, Word32 pos[])
{
   Word32 i, pos1, pos2;
   Word32 mask;

   mask = ((1 << N) - 1);

   pos1 = (((index >> N) & mask) + offset);
   i = (index >> (2 * N)) & 1;
   pos2 = ((index & mask) + offset);

   if((pos2 - pos1) < 0)
   {
      if(i == 1)
      {
         pos1 += 16;
      }
      else
      {
         pos2 += 16;
      }
   }
   else
   {
      if(i == 1)
      {
         pos1 += 16;
         pos2 += 16;
      }
   }

   pos[0] = pos1;
   pos[1] = pos2;
}

static void D_ACELP_decode_3p_3N1(Word32 index, Word32 N,
                                  Word32 offset, Word32 pos[])
{
   Word32 j, mask, idx;

   mask = ((1 << ((2 * N) - 1)) - 1);
   idx = index & mask;
   j = offset;

   if(((index >> ((2 * N) - 1)) & 1) == 1)
   {
      j += (1 << (N - 1));
   }

   D_ACELP_decode_2p_2N1(idx, N - 1, j, pos);
   mask = ((1 << (N + 1)) - 1);
   idx = (index >> (2 * N)) & mask;
   D_ACELP_decode_1p_N1(idx, N, offset, pos + 2);
}

static void D_ACELP_decode_4p_4N1(Word32 index, Word32 N,
                                  Word32 offset, Word32 pos[])
{
   Word32 j, mask, idx;

   mask = ((1 << ((2 * N) - 1)) - 1);
   idx = index & mask;
   j = offset;

   if(((index >> ((2 * N) - 1)) & 1) == 1)
   {
      j += (1 << (N - 1));
   }

   D_ACELP_decode_2p_2N1(idx, N - 1, j, pos);
   mask = ((1 << ((2 * N) + 1)) - 1);
   idx = (index >> (2 * N)) & mask;
   D_ACELP_decode_2p_2N1(idx, N, offset, pos + 2);
}

static void D_ACELP_decode_4p_4N(Word32 index, Word32 N,
                                 Word32 offset, Word32 pos[])
{
   Word32 j, n_1;

   n_1 = N - 1;
   j = offset + (1 << n_1);

   switch((index >> ((4 * N) - 2)) & 3)
   {
   case 0:
      if(((index >> ((4 * n_1) + 1)) & 1) == 0)
      {
         D_ACELP_decode_4p_4N1(index, n_1, offset, pos);
      }
      else
      {
         D_ACELP_decode_4p_4N1(index, n_1, j, pos);
      }
      break;

   case 1:
      D_ACELP_decode_1p_N1((index >> ((3 * n_1) + 1)), n_1, offset, pos);
      D_ACELP_decode_3p_3N1(index, n_1, j, pos + 1);
      break;

   case 2:
      D_ACELP_decode_2p_2N1((index >> ((2 * n_1) + 1)), n_1, offset, pos);
      D_ACELP_decode_2p_2N1(index, n_1, j, pos + 2);
      break;

   case 3:
      D_ACELP_decode_3p_3N1((index >> (n_1 + 1)), n_1, offset, pos);
      D_ACELP_decode_1p_N1(index, n_1, j, pos + 3);
      break;
   }
}

static void D_ACELP_decode_5p_5N(Word32 index, Word32 N,
                                 Word32 offset, Word32 pos[])
{
   Word32 j, n_1;
   Word32 idx;

   /*
    * Decode 5 pulses with 5*N bits
    */
   n_1 = N - 1;
   j = offset + (1 << n_1);
   idx = (index >> ((2 * N) + 1));

   if(((index >> ((5 * N) - 1)) & 1) == 0)
   {
      D_ACELP_decode_3p_3N1(idx, n_1, offset, pos);
      D_ACELP_decode_2p_2N1(index, N, offset, pos + 3);
   }
   else
   {
      D_ACELP_decode_3p_3N1(idx, n_1, j, pos);
      D_ACELP_decode_2p_2N1(index, N, offset, pos + 3);
   }
}

static void D_ACELP_decode_6p_6N_2(Word32 index, Word32 N,
                                   Word32 offset, Word32 pos[])
{
   Word32 j, n_1, offsetA, offsetB;

   n_1 = N - 1;
   j = offset + (1 << n_1);
   offsetA = offsetB = j;

   if(((index >> ((6 * N) - 5)) & 1) == 0)
   {
      offsetA = offset;
   }
   else
   {
      offsetB = offset;
   }

   switch((index >> ((6 * N) - 4)) & 3)
   {
      case 0:
         D_ACELP_decode_5p_5N(index >> N, n_1, offsetA, pos);
         D_ACELP_decode_1p_N1(index, n_1, offsetA, pos + 5);
         break;

      case 1:
         D_ACELP_decode_5p_5N(index >> N, n_1, offsetA, pos);
         D_ACELP_decode_1p_N1(index, n_1, offsetB, pos + 5);
         break;

      case 2:
         D_ACELP_decode_4p_4N(index >> ((2 * n_1) + 1), n_1, offsetA, pos);
         D_ACELP_decode_2p_2N1(index, n_1, offsetB, pos + 4);
         break;

      case 3:
         D_ACELP_decode_3p_3N1(index >> ((3 * n_1) + 1), n_1, offset, pos);
         D_ACELP_decode_3p_3N1(index, n_1, j, pos + 3);
         break;
   }
}

void D_ACELP_decode_2t(Word16 index, Word16 code[])
{
   Word32 i0, i1;

   memset(code, 0, 64 * sizeof(Word16));

   i0 = (index >> 5) & 0x0000003E;
   i1 = ((index & 0x0000001F) << 1) + 1;

   if(((index >> 6) & 32) == 0)
   {
      code[i0] = 512;
   }
   else
   {
      code[i0] = -512;
   }

   if((index & 32) == 0)
   {
      code[i1] = 512;
   }
   else
   {
      code[i1] = -512;
   }
}

void D_ACELP_decode_4t(Word16 index[], Word16 nbbits, Word16 code[])
{
   Word32 k, L_index, pos[6];

   memset(code, 0, 64 * sizeof(Word16));

   if(nbbits == 20)
   {
      for(k = 0; k < 4; k++)
      {
         L_index = index[k];
         D_ACELP_decode_1p_N1(L_index, 4, 0, pos);
         D_ACELP_add_pulse(pos, 1, k, code);
      }
   }
   else if(nbbits == 36)
   {
      for(k = 0; k < 4; k++)
      {
         L_index = index[k];
         D_ACELP_decode_2p_2N1(L_index, 4, 0, pos);
         D_ACELP_add_pulse(pos, 2, k, code);
      }
   }
   else if(nbbits == 44)
   {
      for(k = 0; k < 4 - 2; k++)
      {
         L_index = index[k];
         D_ACELP_decode_3p_3N1(L_index, 4, 0, pos);
         D_ACELP_add_pulse(pos, 3, k, code);
      }

      for(k = 2; k < 4; k++)
      {
         L_index = index[k];
         D_ACELP_decode_2p_2N1(L_index, 4, 0, pos);
         D_ACELP_add_pulse(pos, 2, k, code);
      }
   }
   else if(nbbits == 52)
   {
      for(k = 0; k < 4; k++)
      {
         L_index = index[k];
         D_ACELP_decode_3p_3N1(L_index, 4, 0, pos);
         D_ACELP_add_pulse(pos, 3, k, code);
      }
   }
   else if(nbbits == 64)
   {
      for(k = 0; k < 4; k++)
      {
         L_index = ((index[k] << 14) + index[k + 4]);
         D_ACELP_decode_4p_4N(L_index, 4, 0, pos);
         D_ACELP_add_pulse(pos, 4, k, code);
      }
   }
   else if(nbbits == 72)
   {
      for(k = 0; k < 4 - 2; k++)
      {
         L_index = ((index[k] << 10) + index[k + 4]);
         D_ACELP_decode_5p_5N(L_index, 4, 0, pos);
         D_ACELP_add_pulse(pos, 5, k, code);
      }

      for(k = 2; k < 4; k++)
      {
         L_index = ((index[k] << 14) + index[k + 4]);
         D_ACELP_decode_4p_4N(L_index, 4, 0, pos);
         D_ACELP_add_pulse(pos, 4, k, code);
      }
   }
   else if(nbbits == 88)
   {
      for(k = 0; k < 4; k++)
      {
         L_index = ((index[k] << 11) + index[k + 4]);
         D_ACELP_decode_6p_6N_2(L_index, 4, 0, pos);
         D_ACELP_add_pulse(pos, 6, k, code);
      }
   }
}

void D_ACELP_phase_dispersion(Word16 gain_code, Word16 gain_pit, Word16 code[], Word16 mode, Word16 disp_mem[])
{
   Word32 code2[2 * L_SUBFR] = {0};
   Word32 i, j, state;
   Word16 *prev_gain_pit, *prev_gain_code, *prev_state;

   prev_state = disp_mem;
   prev_gain_code = disp_mem + 1;
   prev_gain_pit = disp_mem + 2;

   if(gain_pit < 9830)   /* 0.6 in Q14 */
   {
      state = 0;
   }
   else if(gain_pit < 14746)   /* 0.9 in Q14 */
   {
      state = 1;
   }
   else
   {
      state = 2;
   }

   for(i = 5; i > 0; i--)
   {
      prev_gain_pit[i] = prev_gain_pit[i - 1];
   }
   prev_gain_pit[0] = gain_pit;

   if((gain_code - *prev_gain_code) > (*prev_gain_code << 1))
   {
      /* onset */
      if(state < 2)
      {
         state = state + 1;
      }
   }
   else
   {
      j = 0;

      for(i = 0; i < 6; i++)
      {
         if(prev_gain_pit[i] < 9830)   /* 0.6 in Q14 */
            j = (j + 1);
      }

      if(j > 2)
      {
         state = 0;
      }

      if((state - *prev_state) > 1)
      {
         state = state - 1;
      }
   }
   *prev_gain_code = gain_code;
   *prev_state = (Word16)state;

   state = state + mode;   /* level of dispersion */

   if(state == 0)
   {
      for(i = 0; i < L_SUBFR; i++)
      {
         if(code[i] != 0)
         {
            for(j = 0; j < L_SUBFR; j++)
            {
               code2[i + j] = code2[i + j] +
                  (((code[i] * D_ROM_ph_imp_low[j]) + 0x4000) >> 15);
            }
         }
      }
   }
   else if(state == 1)
   {
      for(i = 0; i < L_SUBFR; i++)
      {
         if(code[i] != 0)
         {
            for(j = 0; j < L_SUBFR; j++)
            {
               code2[i + j] = code2[i + j] +
                  (((code[i] * D_ROM_ph_imp_mid[j]) + 0x4000) >> 15);
            }
         }
      }
   }

   if(state < 2)
   {
      for(i = 0; i < L_SUBFR; i++)
      {
         code[i] = (Word16)(code2[i] + code2[i + L_SUBFR]);
      }
   }
}
