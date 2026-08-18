/*
 * ===================================================================
 *  TS 26.104
 *  REL-5 V5.4.0 2004-03
 *  REL-6 V6.1.0 2004-03
 *  3GPP AMR Floating-point Speech Codec
 * ===================================================================
 *
 */
 
#include "../amr_config.h"
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <float.h>
#include "sp_enc.h"
#include "rom_enc.h"
#include "sp_enc_vad.h"
#include "sp_enc_vad2.h"

static Float64 Dotproduct40( Float32 *x, Float32 *y )
{
   Float64 acc = 0.0;

   for ( Word32 i = 0; i < 40; i++ ) {
      acc += x[i] * y[i];
   }

   return acc;
}

static void Autocorr( Float32 x[], Float32 r[], const Float32 wind[] )
{
   Float32 y[L_WINDOW + M + 1];   /* Windowed signal */

   for ( int i = 0; i < L_WINDOW; i++ ) {
      y[i] = x[i] * wind[i];
   }

   memset( &y[L_WINDOW], 0, 44 );

   for ( int i = 0; i <= M; i++ ) {
      Float64 sum = 0;

      for ( int j = 0; j < L_WINDOW; j += 40 ) {
         sum += Dotproduct40( &y[j], &y[j + i] );
      }

      r[i] = (Float32)sum;
   }
}

static void Levinson( Float32 *old_A, Float32 *r, Float32 *A, Float32 *rc )
{
   Float32 rct[M];

   rct[0] = ( -r[1] ) / r[0];
   A[0] = 1.0F;
   A[1] = rct[0];
   Float32 err = r[0] + r[1] * rct[0];

   if ( err <= 0.0 )
      err = 0.01F;

   for ( int i = 2; i <= M; i++ ) {
      Float32 sum = 0.0F;
      for ( int j = 0; j < i; j++ )
         sum += r[i - j] * A[j];

      rct[i - 1] = ( -sum ) / ( err );

      for ( int j = 1; j <= ( i / 2 ); j++ ) {
         l = i - j;
         Float32 at = A[j] + rct[i - 1] *A[l];
         A[l] += rct[i - 1] *A[j];
         A[j] = at;
      }

      A[i] = rct[i - 1];
      err += rct[i - 1] *sum;

      if ( err <= 0.0 )
         err = 0.01F;
   }

   memcpy( rc, rct, 4 * sizeof( Float32 ) );
   memcpy( old_A, A, MP1 * sizeof( Float32 ) );
}

static void lpc( Float32 *old_A, Float32 x[], Float32 x_12k2[], Float32 a[], enum ModeNB mode )
{
   Float32 r[MP1], rc[4];

   if ( mode == MR122 ) {
      Autocorr( x_12k2, r, window_160_80 );

      for ( int i = 1; i <= M; i++ ) {
         r[i] = r[i] * lag_wind[i - 1];
      }

      r[0] *= 1.0001F;
      if ( r[0] < 1.0F ) r[0] = 1.0F;

      Levinson( old_A, r, &a[MP1], rc );
      Autocorr( x_12k2, r, window_232_8 );

      for ( int i = 1; i <= M; i++ ) {
         r[i] = r[i] * lag_wind[i - 1];
      }

      r[0] *= 1.0001F;

      if ( r[0] < 1.0F ) r[0] = 1.0F;

      Levinson( old_A, r, &a[MP1 * 3], rc );
   }
   else {
      Autocorr( x, r, window_200_40 );
      for ( int i = 1; i <= M; i++ ) {
         r[i] = r[i] * lag_wind[i - 1];
      }

      r[0] *= 1.0001F;
      if ( r[0] < 1.0F ) r[0] = 1.0F;

      Levinson( old_A, r, &a[MP1 * 3], rc );
   }
}

static Float32 Chebps( Float32 x, Float32 f[] )
{
   Float32 x2 = 2.0F * x;
   Float32 b2 = 1.0F;
   Float32 b1 = x2 + f[1];

   for ( int i = 2; i < 5; i++ ) {
      Float32 b0 = x2 * b1 - b2 + f[i];
      b2 = b1;
      b1 = b0;
   }

   return( x * b1 - b2 + f[5] );
}

static void Az_lsp( Float32 a[], Float32 lsp[], Float32 old_lsp[] )
{
   Float32 xlow, ylow, xhigh, yhigh, xmid, ymid, xint;
   Float32 y;
   Float32 *coef;
   Float32 f1[6], f2[6];

   f1[0] = 1.0F;
   f2[0] = 1.0F;

   for ( Word32 i = 0; i < ( NC ); i++ ) {
      f1[i + 1] = a[i + 1] +a[M - i] - f1[i];
      f2[i + 1] = a[i + 1] -a[M - i] + f2[i];
   }
   f1[NC] *= 0.5F;
   f2[NC] *= 0.5F;

   Word32 nf = 0;   /* number of found frequencies */
   Word32 ip = 0;   /* indicator for f1 or f2 */
   coef = f1;
   xlow = grid[0];
   ylow = Chebps( xlow, coef );
   Word32 j = 0;

   while ( ( nf < M ) && ( j < 60 ) ) {
      j++;
      xhigh = xlow;
      yhigh = ylow;
      xlow = grid[j];
      ylow = Chebps( xlow, coef );

      if ( ylow * yhigh <= 0 ) {
         for ( Word32 i = 0; i < 4; i++ ) {
            xmid = ( xlow + xhigh ) * 0.5F;
            ymid = Chebps( xmid, coef );

            if ( ylow * ymid <= 0.0F ) {
               yhigh = ymid;
               xhigh = xmid;
            }
            else {
               ylow = ymid;
               xlow = xmid;
            }
         }

         y = yhigh - ylow;

         if ( y == 0 ) {
            xint = xlow;
         }
         else {
            y = ( xhigh - xlow ) / ( yhigh - ylow );
            xint = xlow - ylow * y;
         }
         lsp[nf] = xint;
         xlow = xint;
         nf++;

         if ( ip == 0 ) {
            ip = 1;
            coef = f2;
         }
         else {
            ip = 0;
            coef = f1;
         }
         ylow = Chebps( xlow, coef );
      }
   }

   if ( nf < M ) {
      memcpy( lsp, old_lsp, M <<2 );
   }
}

static void Lsp_Az( Float32 lsp[], Float32 a[] )
{
   Float32 f1[6], f2[6];

   Get_lsp_pol( &lsp[0], f1 );
   Get_lsp_pol( &lsp[1], f2 );

   for ( Word32 i = 5; i > 0; i-- ) {
      f1[i] += f1[i - 1];
      f2[i] -= f2[i - 1];
   }
   a[0] = 1;

   for ( Word32 i = 1, j = 10; i <= 5; i++, j-- ) {
      a[i] = ( Float32 )( ( f1[i] + f2[i] ) * 0.5F );
      a[j] = ( Float32 )( ( f1[i] - f2[i] ) * 0.5F );
   }
}

static void Int_lpc_1and3_2( Float32 lsp_old[], Float32 lsp_mid[], Float32 lsp_new[], Float32 az[] )
{
   Float32 lsp[M];

   for ( Word32 i = 0; i < M; i += 2 ) {
      lsp[i] = ( lsp_mid[i] + lsp_old[i] ) * 0.5F;
      lsp[i + 1] = ( lsp_mid[i + 1] +lsp_old[i+1] ) * 0.5F;
   }

   Lsp_Az( lsp, az );
   az += MP1 * 2;

   for ( Word32 i = 0; i < M; i += 2 ) {
      lsp[i] = ( lsp_mid[i] + lsp_new[i] ) * 0.5F;
      lsp[i + 1] = ( lsp_mid[i + 1] +lsp_new[i+1] ) * 0.5F;
   }

   Lsp_Az( lsp, az );
}

static void Lsp_lsf( Float32 lsp[], Float32 lsf[] )
{
   for ( Word32 i = 0; i < M; i++ ) {
      lsf[i] = ( Float32 )( acos( lsp[i] )*SCALE_LSP_FREQ );
   }
}

static void Lsf_wt( Float32 *lsf, Float32 *wf )
{
   wf[0] = lsf[1];

   for ( Word32 i = 1; i < 9; i++ ) {
      wf[i] = lsf[i + 1] -lsf[i - 1];
   }
   wf[9] = 4000.0F - lsf[8];

   for ( Word32 i = 0; i < 10; i++ ) {
      Float32 temp = wf[i] < 450.0F
         ? 3.347F - SLOPE1_WGHT_LSF * wf[i]
         : temp = 1.8F - SLOPE2_WGHT_LSF * ( wf[i] - 450.0F );
      wf[i] = temp * temp;
   }
}

static Word16 Vq_subvec( Float32 *lsf_r1, Float32 *lsf_r2, const Float32 *dico, Float32 *wf1, Float32 *wf2, Word16 dico_size )
{
   Word32 index = 0;
   Float64 dist_min = DBL_MAX;

   for ( Word32 i = 0; i < dico_size; i++ ) {
      Float64 temp = lsf_r1[0] - dico[i * 4 + 0];
      Float64 dist = temp * temp * wf1[0];

      temp = lsf_r1[1] - dico[i * 4 + 1];
      dist += temp * temp * wf1[1];

      temp = lsf_r2[0] - dico[i * 4 + 2];
      dist += temp * temp * wf2[0];

      temp = lsf_r2[1] - dico[i * 4 + 3];
      dist += temp * temp * wf2[1];

      if ( dist < dist_min ) {
         dist_min = dist;
         index = i;
      }
   }

   lsf_r1[0] = dico[index * 4 + 0];
   lsf_r1[1] = dico[index * 4 + 1];
   lsf_r2[0] = dico[index * 4 + 2];
   lsf_r2[1] = dico[index * 4 + 3];

   return ( Word16 )index;
}

static Word16 Vq_subvec_s( Float32 *lsf_r1, Float32 *lsf_r2, const Float32 *dico, Float32 *wf1, Float32 *wf2, Word16 dico_size )
{
   Word32 i, index = 0;
   Word16 sign = 0;

   Float64 dist_min = DBL_MAX;

   for ( i = 0; i < dico_size; i++ ) {
      Float64 temp1 = lsf_r1[0] - dico[i * 4 + 0];
      Float64 temp2 = lsf_r1[0] + dico[i * 4 + 0];
      Float64 dist1 = temp1 * temp1 * wf1[0];
      Float64 dist2 = temp2 * temp2 * wf1[0];

      temp1 = lsf_r1[1] - dico[i * 4 + 1];
      temp2 = lsf_r1[1] + dico[i * 4 + 1];
      dist1 += temp1 * temp1 * wf1[1];
      dist2 += temp2 * temp2 * wf1[1];

      temp1 = lsf_r2[0] - dico[i * 4 + 2];
      temp2 = lsf_r2[0] + dico[i * 4 + 2];
      dist1 += temp1 * temp1 * wf2[0];
      dist2 += temp2 * temp2 * wf2[0];

      temp1 = lsf_r2[1] - dico[i * 4 + 3];
      temp2 = lsf_r2[1] + dico[i * 4 + 3];
      dist1 += temp1 * temp1 * wf2[1];
      dist2 += temp2 * temp2 * wf2[1];

      if ( dist1 < dist_min ) {
         dist_min = dist1;
         index = i;
         sign = 0;
      }

      if ( dist2 < dist_min ) {
         dist_min = dist2;
         index = i;
         sign = 1;
      }
   }

   if ( sign == 0 ) {
      lsf_r1[0] = dico[index * 4 + 0];
      lsf_r1[1] = dico[index * 4 + 1];
      lsf_r2[0] = dico[index * 4 + 2];
      lsf_r2[1] = dico[index * 4 + 3];
   }
   else {
      lsf_r1[0] = -dico[index * 4 + 0];
      lsf_r1[1] = -dico[index * 4 + 1];
      lsf_r2[0] = -dico[index * 4 + 2];
      lsf_r2[1] = -dico[index * 4 + 3];
   }

   index = index << 1;
   index = index + sign;

   return ( Word16 )index;
}

static void Reorder_lsf( Float32 *lsf, Float32 min_dist )
{
   Float32 lsf_min = min_dist;

   for ( Word32 i = 0; i < M; i++ ) {
      if ( lsf[i] < lsf_min ) {
         lsf[i] = lsf_min;
      }
      lsf_min = lsf[i] + min_dist;
   }
}

static void Lsf_lsp( Float32 lsf[], Float32 lsp[] )
{
   for ( Word32 i = 0; i < M; i++ ) {
      lsp[i] = ( Float32 )cos( SCALE_FREQ_LSP * lsf[i] );
   }
}

static Word16 Vq_subvec3( Float32 *lsf_r1, const Float32 *dico, Float32 *wf1, Word16 dico_size, Word32 use_half )
{
   Float64 dist, dist_min;
   Float32 temp;
   Word32 i, index = 0;
   Word32 dico_index;

   dist_min = FLT_MAX;

   if ( use_half == 0 ) {
      for ( i = 0; i < dico_size; i++ ) {
         dico_index = 3 * i;

         temp = lsf_r1[0] - dico[dico_index + 0];
         temp *= wf1[0];
         dist = temp * temp;

         temp = lsf_r1[1] - dico[dico_index + 1];
         temp *= wf1[1];
         dist += temp * temp;

         temp = lsf_r1[2] - dico[dico_index + 2];
         temp *= wf1[2];
         dist += temp * temp;

         if ( dist < dist_min ) {
            dist_min = dist;
            index = i;
         }
      }

      dico_index = 3 * index;
   }
   else {
      for ( i = 0; i < dico_size; i++ ) {
         dico_index = 6 * i;

         temp = lsf_r1[0] - dico[dico_index + 0];
         temp *= wf1[0];
         dist = temp * temp;

         temp = lsf_r1[1] - dico[dico_index + 1];
         temp *= wf1[1];
         dist += temp * temp;

         temp = lsf_r1[2] - dico[dico_index + 2];
         temp *= wf1[2];
         dist += temp * temp;

         if ( dist < dist_min ) {
            dist_min = dist;
            index = i;
         }
      }

      dico_index = 6 * index;
   }

   lsf_r1[0] = dico[dico_index + 0];
   lsf_r1[1] = dico[dico_index + 1];
   lsf_r1[2] = dico[dico_index + 2];

   return ( Word16 )index;
}

static Word16 Vq_subvec4( Float32 *lsf_r1, const Float32 *dico, Float32 *wf1, Word16 dico_size )
{
   Float64 dist, dist_min;
   Float32 temp;
   Word32 i, index = 0;
   Word32 dico_index;

   dist_min = FLT_MAX;

   for ( i = 0; i < dico_size; i++ ) {
      dico_index = 4 * i;

      temp = lsf_r1[0] - dico[dico_index + 0];
      temp *= wf1[0];
      dist = temp * temp;

      temp = lsf_r1[1] - dico[dico_index + 1];
      temp *= wf1[1];
      dist += temp * temp;

      temp = lsf_r1[2] - dico[dico_index + 2];
      temp *= wf1[2];
      dist += temp * temp;

      temp = lsf_r1[3] - dico[dico_index + 3];
      temp *= wf1[3];
      dist += temp * temp;

      if ( dist < dist_min ) {
         dist_min = dist;
         index = i;
      }
   }

   dico_index = 4 * index;

   lsf_r1[0] = dico[dico_index + 0];
   lsf_r1[1] = dico[dico_index + 1];
   lsf_r1[2] = dico[dico_index + 2];
   lsf_r1[3] = dico[dico_index + 3];

   return ( Word16 )index;
}

static void Q_plsf_3( enum ModeNB mode, Float32 *past_rq, Float32 *lsp1, Float32 *lsp1_q, Word16 *indice, Word32 *pred_init_i )
{
   Float32 lsf1[M], wf1[M], lsf_p[M], lsf_r1[M];
   Float32 lsf1_q[M];
   Float32 pred_init_err;
   Float32 min_pred_init_err;
   Float32 temp_r1[M];
   Float32 temp_p[M];

   Lsp_lsf( lsp1, lsf1 );
   Lsf_wt( lsf1, wf1 );

   if ( mode != MRDTX ) {
      for ( Word32 i = 0; i < M; i++ ) {
         lsf_p[i] = mean_lsf_3[i] + past_rq[i] * pred_fac[i];
         lsf_r1[i] = lsf1[i] - lsf_p[i];
      }
   }
   else {
      *pred_init_i = 0;
      min_pred_init_err = FLT_MAX;

      for ( Word32 j = 0; j < PAST_RQ_INIT_SIZE; j++ ) {
         pred_init_err = 0;

         for ( Word32 i = 0; i < M; i++ ) {
            temp_p[i] = mean_lsf_3[i] + past_rq_init[j * M + i];
            temp_r1[i] = lsf1[i] - temp_p[i];
            pred_init_err += temp_r1[i] * temp_r1[i];
         }

         if ( pred_init_err < min_pred_init_err ) {
            min_pred_init_err = pred_init_err;
            memcpy( lsf_r1, temp_r1, M <<2 );
            memcpy( lsf_p, temp_p, M <<2 );
            memcpy( past_rq, &past_rq_init[j * M], M <<2 );
            *pred_init_i = j;
         }
      }
   }

   if ( ( mode == MR475 ) || ( mode == MR515 ) ) {
      indice[0] = Vq_subvec3( &lsf_r1[0], dico1_lsf_3, &wf1[0], DICO1_SIZE_3, 0 );
      indice[1] = Vq_subvec3( &lsf_r1[3], dico2_lsf_3, &wf1[3], DICO2_SIZE_3 /2, 1 );
      indice[2] = Vq_subvec4( &lsf_r1[6], mr515_3_lsf, &wf1[6], MR515_3_SIZE );
   }
   else if ( mode == MR795 ) {
      indice[0] = Vq_subvec3( &lsf_r1[0], mr795_1_lsf, &wf1[0], MR795_1_SIZE, 0 );
      indice[1] = Vq_subvec3( &lsf_r1[3], dico2_lsf_3, &wf1[3], DICO2_SIZE_3, 0 );
      indice[2] = Vq_subvec4( &lsf_r1[6], dico3_lsf_3, &wf1[6], DICO3_SIZE_3 );
   }
   else { /* MR59, MR67, MR74, MR102 , MRDTX */
      indice[0] = Vq_subvec3( &lsf_r1[0], dico1_lsf_3, &wf1[0], DICO1_SIZE_3, 0 );
      indice[1] = Vq_subvec3( &lsf_r1[3], dico2_lsf_3, &wf1[3], DICO2_SIZE_3, 0 );
      indice[2] = Vq_subvec4( &lsf_r1[6], dico3_lsf_3, &wf1[6], DICO3_SIZE_3 );
   }

   for ( i = 0; i < M; i++ ) {
      lsf1_q[i] = lsf_r1[i] + lsf_p[i];
      past_rq[i] = lsf_r1[i];
   }

   Reorder_lsf( lsf1_q, 50.0F );
   Lsf_lsp( lsf1_q, lsp1_q );
}

static void Q_plsf_5( Float32 *past_rq, Float32 *lsp1, Float32 *lsp2, Float32 *lsp1_q, Float32 *lsp2_q, Word16 *indice )
{
   Float32 lsf1[M], lsf2[M], wf1[M], wf2[M], lsf_p[M], lsf_r1[M], lsf_r2[M];
   Float32 lsf1_q[M], lsf2_q[M];

   Lsp_lsf( lsp1, lsf1 );
   Lsp_lsf( lsp2, lsf2 );
   Lsf_wt( lsf1, wf1 );
   Lsf_wt( lsf2, wf2 );

   for ( Word32 i = 0; i < M; i++ ) {
      lsf_p[i] = mean_lsf_5[i] + past_rq[i] * 0.65F;
      lsf_r1[i] = lsf1[i] - lsf_p[i];
      lsf_r2[i] = lsf2[i] - lsf_p[i];
   }

   indice[0] = Vq_subvec( &lsf_r1[0], &lsf_r2[0], dico1_lsf_5, &wf1[0], &wf2[0], DICO1_SIZE_5 );
   indice[1] = Vq_subvec( &lsf_r1[2], &lsf_r2[2], dico2_lsf_5, &wf1[2], &wf2[2], DICO2_SIZE_5 );
   indice[2] = Vq_subvec_s( &lsf_r1[4], &lsf_r2[4], dico3_lsf_5, &wf1[4], &wf2[4], DICO3_SIZE_5 );
   indice[3] = Vq_subvec( &lsf_r1[6], &lsf_r2[6], dico4_lsf_5, &wf1[6], &wf2[6], DICO4_SIZE_5 );
   indice[4] = Vq_subvec( &lsf_r1[8], &lsf_r2[8], dico5_lsf_5, &wf1[8], &wf2[8], DICO5_SIZE_5 );

   for ( Word32 i = 0; i < M; i++ ) {
      lsf1_q[i] = lsf_r1[i] + lsf_p[i];
      lsf2_q[i] = lsf_r2[i] + lsf_p[i];
      past_rq[i] = lsf_r2[i];
   }

   Reorder_lsf( lsf1_q, 50.0F );
   Reorder_lsf( lsf2_q, 50.0F );
   Lsf_lsp( lsf1_q, lsp1_q );
   Lsf_lsp( lsf2_q, lsp2_q );
}

static void Int_lpc_1and3( Float32 lsp_old[], Float32 lsp_mid[], Float32 lsp_new[], Float32 az[] )
{
   Float32 lsp[M];

   for ( Word32 i = 0; i < M; i++ ) {
      lsp[i] = ( lsp_mid[i] + lsp_old[i] ) * 0.5F;
   }

   Lsp_Az( lsp, az );
   az += MP1;

   Lsp_Az( lsp_mid, az );
   az += MP1;

   for ( Word32 i = 0; i < M; i++ ) {
      lsp[i] = ( lsp_mid[i] + lsp_new[i] ) * 0.5F;
   }

   Lsp_Az( lsp, az );
   az += MP1;

   Lsp_Az( lsp_new, az );
}

static void Int_lpc_1to3_2( Float32 lsp_old[], Float32 lsp_new[], Float32 az[] )
{
   Float32 lsp[M];

   for ( Word32 i = 0; i < M; i += 2 ) {
      lsp[i] = lsp_new[i] * 0.25F + lsp_old[i] * 0.75F;
      lsp[i + 1] = lsp_new[i + 1] *0.25F + lsp_old[i + 1] *0.75F;
   }

   Lsp_Az( lsp, az );
   az += MP1;

   for ( Word32 i = 0; i < M; i += 2 ) {
      lsp[i] = ( lsp_old[i] + lsp_new[i] ) * 0.5F;
      lsp[i + 1] = ( lsp_old[i + 1] +lsp_new[i+1] )*0.5F;
   }

   Lsp_Az( lsp, az );
   az += MP1;

   for ( Word32 i = 0; i < M; i += 2 ) {
      lsp[i] = lsp_old[i] * 0.25F + lsp_new[i] * 0.75F;
      lsp[i + 1] = lsp_old[i + 1] *0.25F + lsp_new[i + 1] *0.75F;
   }

   Lsp_Az( lsp, az );
}

static void Int_lpc_1to3( Float32 lsp_old[], Float32 lsp_new[], Float32 az[] )
{
   Float32 lsp[M];

   for ( Word32 i = 0; i < M; i++ ) {
      lsp[i] = lsp_new[i] * 0.25F + lsp_old[i] * 0.75F;
   }

   Lsp_Az( lsp, az );
   az += MP1;

   for ( Word32 i = 0; i < M; i++ ) {
      lsp[i] = ( lsp_old[i] + lsp_new[i] ) * 0.5F;
   }

   Lsp_Az( lsp, az );
   az += MP1;

   for ( Word32 i = 0; i < M; i++ ) {
      lsp[i] = lsp_old[i] * 0.25F + lsp_new[i] * 0.75F;
   }

   Lsp_Az( lsp, az );
   az += MP1;

   Lsp_Az( lsp_new, az );
}

static void lsp( enum ModeNB req_mode, enum ModeNB used_mode, Float32 *lsp_old,
      Float32 *lsp_old_q, Float32 *past_rq, Float32 az[], Float32 azQ[], Float32
      lsp_new[], Word16 **anap )
{
   Float32 lsp_new_q[M];   /* LSPs at 4th subframe */
   Float32 lsp_mid[M], lsp_mid_q[M];   /* LSPs at 2nd subframe */
   Word32 pred_init_i;   /* init index for MA prediction in DTX mode */

   if ( req_mode == MR122 ) {
      Az_lsp( &az[MP1], lsp_mid, lsp_old );
      Az_lsp( &az[MP1 * 3], lsp_new, lsp_mid );

      Int_lpc_1and3_2( lsp_old, lsp_mid, lsp_new, az );

      if ( used_mode != MRDTX ) {
         Q_plsf_5( past_rq, lsp_mid, lsp_new, lsp_mid_q, lsp_new_q, *anap );
         Int_lpc_1and3( lsp_old_q, lsp_mid_q, lsp_new_q, azQ );
         ( *anap ) += 5;
      }
   }
   else {
      Az_lsp( &az[MP1 * 3], lsp_new, lsp_old );
      Int_lpc_1to3_2( lsp_old, lsp_new, az );

      if ( used_mode != MRDTX ) {
         Q_plsf_3( req_mode, past_rq, lsp_new, lsp_new_q, *anap, &pred_init_i );
         Int_lpc_1to3( lsp_old_q, lsp_new_q, azQ );
         ( *anap ) += 3;
      }
   }

   memcpy( lsp_old, lsp_new, M <<2 );
   memcpy( lsp_old_q, lsp_new_q, M <<2 );
}

static Word16 check_lsp( Word16 *count, Float32 *lsp )
{
   Float32 dist_min1 = FLT_MAX;

   for ( Word32 i = 3; i < 8; i++ ) {
      Float32 dist = lsp[i] - lsp[i + 1];

      if ( dist < dist_min1 ) {
         dist_min1 = dist;
      }
   }
   Float32 dist_min2 = FLT_MAX;

   for ( Word32 i = 1; i < 3; i++ ) {
      Float32 dist = lsp[i] - lsp[i + 1];

      if ( dist < dist_min2 ) {
         dist_min2 = dist;
      }
   }

   Float32 dist_th = lsp[1] > 0.98F ? 0.018F : lsp[1] > 0.93F ? 0.024F : 0.034F;
   *count = (dist_min1 < 0.046F || dist_min2 < dist_th) ? *count + 1 : 0;

   if ( *count >= 12 ) {
      *count = 12;
      return 1;
   }
   else {
      return 0;
   }
}

static void Weight_Ai( Float32 a[], const Float32 fac[], Float32 a_exp[] )
{
   a_exp[0] = a[0];
   for ( Word32 i = 1; i <= M; i++ ) {
      a_exp[i] = a[i] * fac[i - 1];
   }
}

static void Residu( Float32 a[], Float32 x[], Float32 y[] )
{
   for ( Word32 i = 0; i < L_SUBFR; i++ ) {
      Float32 s = x[i] * a[0];
      for ( Word32 j = 1; j <= 10; j++ ) {
         s += x[i - j] * a[j];
      }
      y[i] = s;
   }
}

static void Syn_filt( Float32 a[], Float32 x[], Float32 y[], Float32 mem[], Word16 update )
{
   Float64 tmp[50];

   for ( Word32 i = 0; i < M; i++ ) {
      tmp[i] = mem[i];
   }

   for ( Word32 i = 0; i < L_SUBFR; i++ ) {
      Float64 sum = x[i] * a[0];

      for ( Word32 j = 1; j <= M; j++ ) {
         sum -= a[j] * tmp[M + i - j];
      }

      tmp[M + i] = sum;
      y[i] = ( Float32 )sum;
   }

   if ( update != 0 ) {
      for ( Word32 i = 0; i < M; i++ ) {
         mem[i] = y[L_SUBFR - M + i];
      }
   }
}

static Word32 pre_big( enum ModeNB mode, const Float32 gamma1[], const Float32
      gamma1_12k2[], const Float32 gamma2[], Float32 A_t[], Word16 frame_offset,
      Float32 speech[], Float32 mem_w[], Float32 wsp[] )
{
   Float32 Ap1[MP1], Ap2[MP1];
   Float32* g1 = mode <= MR795 ? gamma1 : gamma1_12k2;
   Word32 offset = frame_offset > 0 ? MP1 << 1 : 0;

   for ( Word32 i = 0; i < 2; i++ ) {
      Weight_Ai( &A_t[offset], g1, Ap1 );
      Weight_Ai( &A_t[offset], gamma2, Ap2 );
      Residu( Ap1, &speech[frame_offset], &wsp[frame_offset] );
      Syn_filt( Ap2, &wsp[frame_offset], &wsp[frame_offset], mem_w, 1 );
      offset += MP1;
      frame_offset += L_SUBFR;
   }

   return 0;
}

static void comp_corr( Float32 sig[], Word32 L_frame, Word32 lag_max, Word32 lag_min, Float32 corr[] )
{
   for ( Word32 i = lag_max; i >= lag_min; i-- ) {
      Float64 T0 = 0.0;

      for ( Word32 j = 0; j < L_frame; j += 40 ) {
         T0 += Dotproduct40( &sig[j], &sig[j - i] );
      }

      corr[-i] = ( Float32 )T0;
   }
}

static Word32 gmed_n( Word32 ind[], Word32 n )
{
   Word32 ix = 0;
   Word32 tmp[9];
   Word32 tmp2[9];

   for ( Word32 i = 0; i < n; i++ ) {
      tmp2[i] = ind[i];
   }

   for ( Word32 i = 0; i < n; i++ ) {
      Word32 max = -32767;

      for ( j = 0; j < n; j++ ) {
         if ( tmp2[j] >= max ) {
            max = tmp2[j];
            ix = j;
         }
      }
      tmp2[ix] = -32768;
      tmp[i] = ix;
   }

   Word32 medianIndex = tmp[( n >>1 )];

   return( ind[medianIndex] );
}

static void subframePreProc( enum ModeNB mode, const Float32 gamma1[], const
      Float32 gamma1_12k2[], const Float32 gamma2[], Float32 *A, Float32 *Aq,
      Float32 *speech, Float32 *mem_err, Float32 *mem_w0, Float32 *zero, Float32
      ai_zero[], Float32 *exc, Float32 h1[], Float32 xn[], Float32 res2[],
      Float32 error[] )
{
   Float32 Ap1[MP1];
   Float32 Ap2[MP1];
   Float32 *g1 = ( mode == MR122 ) || ( mode == MR102 ) ? gamma1_12k2 : gamma1;

   Weight_Ai( A, g1, Ap1 );
   Weight_Ai( A, gamma2, Ap2 );
   memcpy( ai_zero, Ap1, MP1 <<2 );
   Syn_filt( Aq, ai_zero, h1, zero, 0 );
   Syn_filt( Ap2, h1, h1, zero, 0 );
   Residu( Aq, speech, res2 );
   memcpy( exc, res2, L_SUBFR <<2 );
   Syn_filt( Aq, exc, error, mem_err, 0 );
   Residu( Ap1, error, xn );
   Syn_filt( Ap2, xn, xn, mem_w0, 0 );
}

static void getRange( Word32 T0, Word16 delta_low, Word16 delta_range, Word16 pitmin, Word16 pitmax, Word32 *T0_min, Word32 *T0_max )
{
   *T0_min = T0 - delta_low;

   if ( *T0_min < pitmin ) {
      *T0_min = pitmin;
   }

   *T0_max = *T0_min + delta_range;

   if ( *T0_max > pitmax ) {
      *T0_max = pitmax;
      *T0_min = *T0_max - delta_range;
   }
}

static void Norm_Corr( Float32 exc[], Float32 xn[], Float32 h[], Word32 t_min,
      Word32 t_max, Float32 corr_norm[] )
{
   Float32 exc_temp[L_SUBFR];

   Word32 k = -t_min;
   Float32 *p_exc = &exc[ - t_min];

   for ( Word32 j = 0; j < L_SUBFR; j++ ) {
      Float32 sum = 0;

      for ( Word32 i = 0; i <= j; i++ ) {
         sum += p_exc[i] * h[j - i];
      }
      exc_temp[j] = sum;
   }

   for ( Word32 i = t_min; i <= t_max; i++ ) {
      Float32 norm = (Float32)Dotproduct40( exc_temp, exc_temp );
      norm = norm == 0 ? 1.0 : ( Float32 )( 1.0 / ( sqrt( norm ) ) );

      Float32 corr = (Float32)Dotproduct40( xn, exc_temp );
      corr_norm[i] = corr * norm;

      if ( i != t_max ) {
         k--;
         for ( Word32 j = L_SUBFR - 1; j > 0; j-- ) {
            exc_temp[j] = exc_temp[j - 1] + exc[k] * h[j];
         }
         exc_temp[0] = exc[k];
      }
   }
}

static Float32 Interpol_3or6( Float32 *x, Word32 frac, Word16 flag3 )
{
   Float32 *x1, *x2;
   const Float32 *c1, *c2;

   if ( flag3 != 0 ) {
      frac <<= 1;
   }

   if ( frac < 0 ) {
      frac += UP_SAMP_MAX;
      x--;
   }
   x1 = &x[0];
   x2 = &x[1];
   c1 = &b24[frac];
   c2 = &b24[UP_SAMP_MAX - frac];
   Float32 s = 0;

   for ( Word32 i = 0, k = 0; i < L_INTER_SRCH; i++, k += UP_SAMP_MAX ) {
      s += x1[ - i] * c1[k];
      s += x2[i] * c2[k];
   }

   return s;
}

static void searchFrac( Word32 *lag, Word32 *frac, Word16 last_frac, Float32 corr[], Word16 flag3 )
{
   Float32 max = Interpol_3or6( &corr[ * lag], *frac, flag3 );

   for ( Word32 i = *frac + 1; i <= last_frac; i++ ) {
      Float32 corr_int = Interpol_3or6( &corr[ * lag], i, flag3 );

      if ( corr_int > max ) {
         max = corr_int;
         *frac = i;
      }
   }

   if ( flag3 == 0 ) {
      if ( *frac == -3 ) {
         *frac = 3;
         *lag -= 1;
      }
   }
   else {
      if ( *frac == -2 ) {
         *frac = 1;
         *lag -= 1;
      }

      if ( *frac == 2 ) {
         *frac = -1;
         *lag += 1;
      }
   }
}

static Word32 Enc_lag3( Word32 T0, Word32 T0_frac, Word32 T0_prev, Word32 T0_min, Word32 T0_max, Word16 delta_flag, Word16 flag4 )
{
   Word32 index, i, tmp_ind, uplag, tmp_lag;

   if ( delta_flag == 0 ) {
      index = T0 <= 85
         ? T0 * 3 - 58 + T0_frac
         : T0 + 112;
   }
   else {
      if ( flag4 == 0 ) {
         index = 3 * ( T0 - T0_min ) + 2 + T0_frac;
      } else {
         tmp_lag = T0_prev;

         if ( ( tmp_lag - T0_min ) > 5 ) tmp_lag = T0_min + 5;
         if ( ( T0_max - tmp_lag ) > 4 ) tmp_lag = T0_max - 4;

         uplag = T0 + T0 + T0 + T0_frac;
         i = tmp_lag - 2;
         tmp_ind = i + i + i;

         if ( tmp_ind >= uplag ) {
            index = ( T0 - tmp_lag ) + 5;
         }
         else {
            i = tmp_lag + 1;
            i = i + i + i;
            index = i > uplag ? ( uplag - tmp_ind ) + 3 : ( T0 - tmp_lag ) + 11;
         }
      }
   }
   return index;
}

static Word32 Enc_lag6( Word32 T0, Word32 T0_frac, Word32 T0_min, Word16 delta_flag )
{
   if ( delta_flag == 0 ) {
      if ( T0 <= 94 ) {
         return T0 * 6 - 105 + T0_frac;
      }
      else {
         return T0 + 368;
      }
   }
   else {
      return 6 * ( T0 - T0_min ) + 3 + T0_frac;
   }
}

static Word32 Pitch_fr( Word32 *T0_prev_subframe, enum ModeNB mode, Word32 T_op[],
      Float32 exc[], Float32 xn[], Float32 h[], Word16 i_subfr, Word32 *pit_frac
      , Word16 *resu3, Word32 *ana_index )
{
   Float32 corr_v[40];
   Float32 max;
   Float32 *corr;
   Word32 i, t_min, t_max, T0_min, T0_max;
   Word32 lag, frac, tmp_lag;
   Word16 max_frac_lag, flag3, flag4, last_frac;
   Word16 delta_int_low, delta_int_range, delta_frc_low, delta_frc_range;
   Word16 pit_min;
   Word16 frame_offset;
   Word16 delta_search;

   max_frac_lag = mode_dep_parm[mode].max_frac_lag;
   flag3 = mode_dep_parm[mode].flag3;
   frac = mode_dep_parm[mode].first_frac;
   last_frac = mode_dep_parm[mode].last_frac;
   delta_int_low = mode_dep_parm[mode].delta_int_low;
   delta_int_range = mode_dep_parm[mode].delta_int_range;
   delta_frc_low = mode_dep_parm[mode].delta_frc_low;
   delta_frc_range = mode_dep_parm[mode].delta_frc_range;
   pit_min = mode_dep_parm[mode].pit_min;

   delta_search = 1;

   if ( ( i_subfr == 0 ) || ( i_subfr == L_FRAME_BY2 ) ) {
      if ( ( ( mode != MR475 ) && ( mode != MR515 ) ) || ( i_subfr != L_FRAME_BY2 ) ) {
         delta_search = 0;
         frame_offset = 1;

         if ( i_subfr == 0 )
            frame_offset = 0;

         getRange( T_op[frame_offset], delta_int_low, delta_int_range, pit_min, PIT_MAX, &T0_min, &T0_max );
      }
      else {
         getRange( *T0_prev_subframe, delta_frc_low, delta_frc_range, pit_min, PIT_MAX, &T0_min, &T0_max );
      }
   }
   else {
      getRange( *T0_prev_subframe, delta_frc_low, delta_frc_range, pit_min, PIT_MAX, &T0_min, &T0_max );
   }

   t_min = T0_min - L_INTER_SRCH;
   t_max = T0_max + L_INTER_SRCH;
   corr = &corr_v[ - t_min];
   Norm_Corr( exc, xn, h, t_min, t_max, corr );

   max = corr[T0_min];
   lag = T0_min;

   for ( i = T0_min + 1; i <= T0_max; i++ ) {
      if ( corr[i] >= max ) {
         max = corr[i];
         lag = i;
      }
   }

   if ( ( delta_search == 0 ) && ( lag > max_frac_lag ) ) {
      frac = 0;
   }
   else {
      if ( ( delta_search != 0 ) && ( ( mode == MR475 ) || ( mode == MR515 ) || ( mode == MR59 ) || ( mode == MR67 ) ) ) {
         tmp_lag = *T0_prev_subframe;

         if ( ( tmp_lag - T0_min ) > 5 ) tmp_lag = T0_min + 5;
         if ( ( T0_max - tmp_lag ) > 4 ) tmp_lag = T0_max - 4;

         if ( ( lag == tmp_lag ) || ( lag == ( tmp_lag - 1 ) ) ) {
            searchFrac( &lag, &frac, last_frac, corr, flag3 );
         } else if ( lag == ( tmp_lag - 2 ) ) {
            frac = 0;
            searchFrac( &lag, &frac, last_frac, corr, flag3 );
         } else if ( lag == ( tmp_lag + 1 ) ) {
            last_frac = 0;
            searchFrac( &lag, &frac, last_frac, corr, flag3 );
         }
         else {
            frac = 0;
         }
      }
      else
         searchFrac( &lag, &frac, last_frac, corr, flag3 );
   }

   if ( flag3 != 0 ) {
      flag4 = 0;

      if ( ( mode == MR475 ) || ( mode == MR515 ) || ( mode == MR59 ) || ( mode == MR67 ) ) {
         flag4 = 1;
      }

      *ana_index = Enc_lag3( lag, frac, *T0_prev_subframe, T0_min, T0_max, delta_search, flag4 );
   }
   else {
      *ana_index = Enc_lag6( lag, frac, T0_min, delta_search );
   }

   *T0_prev_subframe = lag;
   *resu3 = flag3;
   *pit_frac = frac;
   return( lag );
}

static void Pred_lt_3or6( Float32 exc[], Word32 T0, Word32 frac, Word16 flag3 )
{
   Float32 s;
   Float32 *x0;
   const Float32 *c1;
   const Float32 *c2;
   Word32 i, k;

   x0 = &exc[-T0];

   frac = -frac;

   if ( flag3 != 0 ) {
      frac <<= 1;
   }

   if ( frac < 0 ) {
      frac += UP_SAMP_MAX;
      x0--;
   }

   c1 = &b60[frac];
   c2 = &b60[UP_SAMP_MAX - frac];

   for ( i = 0; i < L_SUBFR; i++ ) {
      s = x0[0] * c1[0] + x0[1] * c2[0];

      for ( k = 1; k < 10; k++ ) {
         s += x0[-k] * c1[k * 6];
         s += x0[k + 1] * c2[k * 6];
      }

      exc[i] = ( Float32 )floor( s + 0.5F );

      x0++;
   }
}

static void Pred_lt_3or6_fixed( Word32 exc[], Word32 T0, Word32 frac, Word32 flag3 )
{
   Word32 s;
   Word32 i, j;
   Word32 *x0;
   const Word32 *c1, *c2;

   x0 = &exc[-T0];
   frac = -frac;

   if ( flag3 != 0 ) {
      frac <<= 1;
   }

   if ( frac < 0 ) {
      frac += 6;
      x0--;
   }

   c1 = &inter6[frac];
   c2 = &inter6[6 - frac];

   for ( i = 0; i < 40; i++ ) {
      s = 0;

      for ( j = 0; j < 10; j++ ) {
         s += x0[-j] * c1[j * 6];
         s += x0[j + 1] * c2[j * 6];
      }

      exc[i] = ( s + 0x4000 ) >> 15;

      x0++;
   }
}

static Float32 G_pitch( Float32 xn[], Float32 y1[], Float32 gCoeff[] )
{
   Float32 sum = (Float32)Dotproduct40( y1, y1 );
   sum += 0.01F;
   gCoeff[0] = sum;
   sum = (Float32)Dotproduct40( xn, y1 );
   gCoeff[1] = sum;
   Float32 gain = ( Float32 )( gCoeff[1] / gCoeff[0] );

   if ( gain < 0.0 )
      gain = 0.0F;

   if ( gain > 1.2 )
      gain = 1.2F;

   return( gain );
}

static Word16 check_gp_clipping( Float32 *gp, Float32 g_pitch )
{
   Float32 sum = g_pitch;

   for ( Word32 i = 0; i < N_FRAME; i++ ) {
      sum += gp[i];
   }

   return sum > 7.6F ? 1 : 0;
}

static Word16 q_gain_pitch( enum ModeNB mode, Float32 gp_limit, Float32 *gain, Float32 gain_cand[], Word32 gain_cind[] )
{
   Float32 err_min, err;
   Word32 i, index;

   err_min = ( Float32 )fabs( *gain - qua_gain_pitch[0] );
   index = 0;

   for ( i = 1; i < NB_QUA_PITCH; i++ ) {
      if ( qua_gain_pitch[i] <= gp_limit ) {
         err = ( Float32 )fabs( *gain - qua_gain_pitch[i] );

         if ( err < err_min ) {
            err_min = err;
            index = i;
         }
      }
   }

   if ( mode == MR795 ) {
      Word32 ii;

      if ( index == 0 ) {
         ii = index;
      }
      else {
         ii = index - 1;

         if ( index == ( NB_QUA_PITCH - 1 ) || ( qua_gain_pitch[index + 1] >
               gp_limit ) ) {
            ii = index - 2;
         }
      }

      for ( i = 0; i < 3; i++ ) {
         gain_cind[i] = ii;
         gain_cand[i] = qua_gain_pitch[ii];
         ii++;
      }
      *gain = qua_gain_pitch[index];
   }
   else {
      *gain = qua_gain_pitch_MR122[index];
   }
   return( Word16 )index;
}

static void cl_ltp( Word32 *T0_prev_subframe, Float32 *gp, enum ModeNB mode,
      Word16 frame_offset, Word32 T_op[], Float32 *h1, Float32 *exc, Float32
      res2[], Float32 xn[], Word16 lsp_flag, Float32 xn2[], Float32 y1[], Word32
      *T0, Word32 *T0_frac, Float32 *gain_pit, Float32 gCoeff[], Word16 **anap,
      Float32 *gp_limit )
{
   Float32 s;
   Word32 i, n;
   Word16 gpc_flag, resu3;   /* flag for upsample resolution */

   Word32 exc_tmp[314];
   Word32 *exc_tmp_p;

   exc_tmp_p = exc_tmp + PIT_MAX + L_INTERPOL;

   *T0 = Pitch_fr( T0_prev_subframe, mode, T_op, exc, xn, h1, frame_offset, T0_frac, &resu3, &i );
   *( *anap )++ = ( Word16 )i;

   for (i = -(PIT_MAX + L_INTERPOL); i < 40; i++)
      exc_tmp_p[i] = (Word32)exc[i];

   Pred_lt_3or6_fixed( exc_tmp_p, *T0, *T0_frac, resu3 );

   for (i = -(PIT_MAX + L_INTERPOL); i < 40; i++)
      exc[i] = (Float32)exc_tmp_p[i];

   for ( n = 0; n < L_SUBFR; n++ ) {
      s = 0;

      for ( i = 0; i <= n; i++ ) {
         s += exc[i] * h1[n - i];
      }
      y1[n] = s;
   }

   *gain_pit = G_pitch( xn, y1, gCoeff );
   gpc_flag = 0;
   *gp_limit = 2.0F;

   if ( ( lsp_flag != 0 ) && ( *gain_pit > 0.95F ) ) {
      gpc_flag = check_gp_clipping( gp, *gain_pit );
   }

   if ( ( mode == MR475 ) || ( mode == MR515 ) ) {
      if ( *gain_pit > 0.85 ) {
         *gain_pit = 0.85F;
      }

      if ( gpc_flag != 0 )
         *gp_limit = GP_CLIP;
   }
   else {
      if ( gpc_flag != 0 ) {
         *gp_limit = GP_CLIP;
         *gain_pit = GP_CLIP;
      }

      if ( mode == MR122 ) {
         *( *anap )++ = q_gain_pitch( MR122, *gp_limit, gain_pit, NULL, NULL );
      }
   }

   for ( i = 0; i < L_SUBFR; i++ ) {
      xn2[i] = xn[i] - y1[i] * *gain_pit;
      res2[i] = res2[i] - exc[i] * *gain_pit;
   }
}

static Float32 DotProduct( Float32 *x, Float32 *y, Word32 len )
{
   Float32 acc = 0.0F;

   for ( Word32 i = 0; i < len; i++ )
      acc += x[i] * y[i];

   return( acc );
}

static void cor_h_x( Float32 h[], Float32 x[], Float32 dn[] )
{
   dn[0] = (Float32)Dotproduct40( h, x );

   for ( Word32 i = 1; i < L_CODE; i++ )
      dn[i] = (Float32)DotProduct( h, &x[i], L_CODE - i );
}

static void set_sign( Float32 dn[], Float32 sign[], Float32 dn2[], Word16 n )
{
   Word32 i, j, k, pos = 0;

   for ( i = 0; i < L_CODE; i++ ) {
      Float32 val = dn[i];

      if ( val >= 0 ) {
         sign[i] = 1.0F;
      }
      else {
         sign[i] = -1.0F;
         val = -val;
      }

      dn[i] = val;
      dn2[i] = val;
   }

   for ( i = 0; i < NB_TRACK; i++ ) {
      for ( k = 0; k < ( 8 - n ); k++ ) {
         Float32 min = FLT_MAX;

         for ( j = i; j < L_CODE; j += STEP ) {
            if ( dn2[j] >= 0 ) {
               val = dn2[j] - min;

               if ( val < 0 ) {
                  min = dn2[j];
                  pos = j;
               }
            }
         }
         dn2[pos] = -1.0F;
      }
   }
}

static void cor_h( Float32 h[], Float32 sign[], Float32 rr[][L_CODE] )
{
   Float32 sum;
   Word32 ii, i, j, k;

   rr[0][0] = ( Float32 )Dotproduct40( h, h );

   sum = 0.0F;

   for ( k = 0; k < 13; k++ ) {
      sum += h[k] * h[k];
      rr[39 - k][39 - k] = sum;
   }

   for ( ii = 1; ii <= 10; ii++ ) {
      j = L_CODE - 1;
      i = j - ii;
      sum = 0.0F;

      for ( k = 0; k < L_CODE - ii; k++ ) {
         sum += h[k] * h[k + ii];

         rr[i][j] = sum * sign[i] * sign[j];
         rr[j][i] = rr[i][j];

         i--;
         j--;
      }
   }
}

static void search_2i40_9bits( Word16 subNr, Float32 dn[], Float32 rr[][L_CODE], Word32 codvec[] )
{
   Float32 ps0, ps1, psk, alp, alp0, alp1, alpk, sq, sq1;
   Word32 i0, i1, ix, i;
   Word16 ipos[2];

   psk = -1;
   alpk = 1;

   for ( i = 0; i < 2; i++ ) {
      codvec[i] = i;
   }

   for ( Word16 track1 = 0; track1 < 2; track1++ ) {
      ipos[0] = startPos[( subNr << 1 ) + ( track1 << 3 )];
      ipos[1] = startPos[( subNr << 1 ) + 1 + ( track1 << 3 )];

      for ( i0 = ipos[0]; i0 < L_CODE; i0 += STEP ) {
         ps0 = dn[i0];
         alp0 = rr[i0][i0];

         sq = -1;
         alp = 1;
         ix = ipos[1];

         for ( i1 = ipos[1]; i1 < L_CODE; i1 += STEP ) {
            ps1 = ps0 + dn[i1];
            alp1 = alp0 + rr[i1][i1] + 2.0F * rr[i0][i1];
            sq1 = ps1 * ps1;

            if ( ( alp * sq1 ) > ( sq * alp1 ) ) {
               sq = sq1;
               alp = alp1;
               ix = i1;
            }
         }

         if ( ( alpk * sq ) > ( psk * alp ) ) {
            psk = sq;
            alpk = alp;
            codvec[0] = i0;
            codvec[1] = ix;
         }
      }
   }
}

static void build_code_2i40_9bits( Word16 subNr, Word32 codvec[], Float32 dn_sign[], Float32 cod[], Float32 h[], Float32 y[], Word16 *anap )
{
   Float32 s;
   Float32 *p0, *p1;
   Word32 _sign[2];
   Word32 i, j, k, track, index, indx = 0, rsign = 0;
   Word8 first, *pt;

   pt = &trackTable[subNr + ( subNr << 2 )];
   memset( cod, 0, 160 );

   for ( k = 0; k < 2; k++ ) {
      i = codvec[k];
      j = ( Word32 )dn_sign[i];
      index = i / 5;
      track = i % 5;
      first = pt[track];

      if ( first == 0 ) {
         track = k == 0 ? 0 : 1;
         index = k == 0 ? index : index << 3;
      } else {
         track = k == 0 ? 0 : 1;
         index = k == 0 ? index + 64 : index << 3;
      }

      if ( j > 0 ) {
         cod[i] = 0.9998779296875F;
         _sign[k] = 1;
         rsign = rsign + ( 1 << track );
      }
      else {
         cod[i] = -1;
         _sign[k] = -1;
      }
      indx = indx + index;
   }
   p0 = h - codvec[0];
   p1 = h - codvec[1];

   for ( i = 0; i < L_CODE; i++ ) {
      s = *p0++ * _sign[0];
      s += *p1++ * _sign[1];
      y[i] = s;
   }
   anap[0] = ( Word16 )indx;
   anap[1] = ( Word16 )rsign;
}

static void code_2i40_9bits( Word16 subNr, Float32 x[], Float32 h[], Word32 T0, Float32 pitch_sharp, Float32 code[], Float32 y[], Word16 *anap )
{
   Float32 rr[L_CODE][L_CODE];
   Float32 dn[L_CODE], dn_sign[L_CODE], dn2[L_CODE];
   Word32 codvec[2];

   if ( ( T0 < L_CODE ) && ( pitch_sharp != 0.0F ) ) {
      for ( Word32 i = T0; i < L_CODE; i++ ) {
         h[i] += h[i - T0] * pitch_sharp;
      }
      cor_h_x( h, x, dn );
      set_sign( dn, dn_sign, dn2, 8 );
      cor_h( h, dn_sign, rr );
      search_2i40_9bits( subNr, dn, rr, codvec );
      build_code_2i40_9bits( subNr, codvec, dn_sign, code, h, y, anap );

      if ( ( T0 < L_CODE ) && ( pitch_sharp != 0.0F ) )
         for ( i = T0; i < L_CODE; i++ ) {
            code[i] += code[i - T0] * pitch_sharp;
         }
   }
}

static void search_2i40_11bits( Float32 dn[], Float32 rr[][L_CODE], Word32 codvec[] )
{
   Float64 alpk, alp, alp0, alp1;
   Float32 psk, ps0, ps1, sq, sq1;
   Word32 i, i0, i1, ix = 0;
   Word16 ipos[2];
   Word16 track1, track2;

   psk = -1;
   alpk = 1;

   for ( i = 0; i < 2; i++ ) {
      codvec[i] = i;
   }

   for ( track1 = 0; track1 < 2; track1++ ) {
      for ( track2 = 0; track2 < 4; track2++ ) {
         ipos[0] = startPos1[track1];
         ipos[1] = startPos2[track2];

         for ( i0 = ipos[0]; i0 < L_CODE; i0 += STEP ) {
            ps0 = dn[i0];
            alp0 = rr[i0][i0] * 0.25F;

            sq = -1;
            alp = 1;
            ix = ipos[1];

            for ( i1 = ipos[1]; i1 < L_CODE; i1 += STEP ) {
               ps1 = ps0 + dn[i1];

               alp1 = alp0 + rr[i1][i1] * 0.25F;
               alp1 += rr[i0][i1] * 0.5F;
               sq1 = ps1 * ps1;

               if ( ( alp * sq1 ) > ( sq * alp1 ) ) {
                  sq = sq1;
                  alp = alp1;
                  ix = i1;
               }
            }

            if ( ( alpk * sq ) > ( psk * alp ) ) {
               psk = sq;
               alpk = alp;
               codvec[0] = i0;
               codvec[1] = ix;
            }
         }
      }
   }
}

static void build_code_2i40_11bits( Word32 codvec[], Float32 dn_sign[], Float32 cod[], Float32 h[], Float32 y[], Word16 *anap )
{
   Float64 s;
   Float32 *p0, *p1;
   Word32 _sign[2];
   Word32 i, j, k, track, index, indx = 0, rsign = 0;

   memset( cod, 0, 160 );

   for ( k = 0; k < 2; k++ ) {
      i = codvec[k];
      j = ( Word16 )dn_sign[i];
      index = i / 5;
      track = i % 5;

      if ( track == 0 ) {
         track = 1;
         index = index << 6;
      }
      else if ( track == 1 ) {
         if ( k == 0 ) {
            track = 0;
            index = index << 1;
         }
         else {
            track = 1;
            index = ( index << 6 ) + 16;
         }
      }
      else if ( track == 2 ) {
         track = 1;
         index = ( index << 6 ) + 32;
      }
      else if ( track == 3 ) {
         track = 0;
         index = ( index << 1 ) + 1;
      }
      else if ( track == 4 ) {
         track = 1;
         index = ( index << 6 ) + 48;
      }

      if ( j > 0 ) {
         cod[i] = 0.9998779296875F;
         _sign[k] = 1;
         rsign = rsign + ( 1 << track );
      }
      else {
         cod[i] = -1;
         _sign[k] = -1;
      }
      indx = indx + index;
   }
   p0 = h - codvec[0];
   p1 = h - codvec[1];

   for ( i = 0; i < L_CODE; i++ ) {
      s = *p0++ * _sign[0];
      s += *p1++ * _sign[1];
      y[i] = ( Float32 )s;
   }
   anap[0] = ( Word16 )indx;
   anap[1] = ( Word16 )rsign;
}

static void code_2i40_11bits( Float32 x[], Float32 h[], Word32 T0, Float32 pitch_sharp, Float32 code[], Float32 y[], Word16 *anap )
{
   Float32 rr[L_CODE][L_CODE];
   Float32 dn[L_CODE], dn2[L_CODE], dn_sign[L_CODE];
   Word32 codvec[2];

   if ( ( T0 < L_CODE ) && ( pitch_sharp != 0.0F ) ) {
      for ( Word32 i = T0; i < L_CODE; i++ ) {
         h[i] = h[i] + ( h[i - T0] * pitch_sharp );
      }
   }

   cor_h_x( h, x, dn );
   set_sign( dn, dn_sign, dn2, 8 );
   cor_h( h, dn_sign, rr );
   search_2i40_11bits( dn, rr, codvec );
   build_code_2i40_11bits( codvec, dn_sign, code, h, y, anap );

   if ( ( T0 < L_CODE ) && ( pitch_sharp != 0.0F ) ) {
      for ( i = T0; i < L_CODE; i++ ) {
         code[i] = code[i] + code[i - T0] * pitch_sharp;
      }
   }
}

static void search_3i40( Float32 dn[], Float32 dn2[], Float32 rr[][L_CODE], Word32 codvec[] )
{
   Float32 psk, ps0, ps1, sq, sq1, alpk, alp, alp0, alp1, ps = 0.0F;
   Float32 *rr2, *rr1, *rr0, *pdn, *pdn_max;
   Word32 ipos[3];
   Word32 i0, i1, i2, ix, i, pos, track1, track2;

   psk = -1.0F;
   alpk = 1.0F;

   for ( track1 = 1; track1 < 4; track1 += 2 ) {
      for ( track2 = 2; track2 < 5; track2 += 2 ) {
         ipos[0] = 0;
         ipos[1] = track1;
         ipos[2] = track2;

         for ( i = 0; i < 3; i++ ) {
            for ( i0 = ipos[0]; i0 < L_CODE; i0 += STEP ) {
               if ( dn2[i0] >= 0 ) {
                  ps0 = dn[i0];
                  alp0 = rr[i0][i0];
                  sq = -1.0F;
                  alp = 1.0F;
                  ps = 0.0F;
                  ix = ipos[1];
                  i1 = ipos[1];
                  rr1 = &rr[i1][i1];
                  rr0 = &rr[i0][i1];
                  pdn = &dn[i1];
                  pdn_max = &dn[L_CODE];

                  do {
                     ps1 = ps0 + *pdn;
                     alp1 = alp0 + *rr1 + 2.0F * *rr0;
                     sq1 = ps1 * ps1;

                     if ( ( alp * sq1 ) > ( sq * alp1 ) ) {
                        sq = sq1;
                        ps = ps1;
                        alp = alp1;
                        ix = ( Word16 )( pdn - dn );
                     }
                     pdn += STEP;
                     rr1 += ( 40 * STEP + STEP );
                     rr0 += STEP;
                  } while ( pdn < pdn_max );
                  i1 = ix;

                  ps0 = ps;
                  alp0 = alp;
                  sq = -1.0F;
                  alp = 1.0F;
                  ps = 0.0F;
                  ix = ipos[2];
                  i2 = ipos[2];
                  rr2 = &rr[i2][i2];
                  rr1 = &rr[i1][i2];
                  rr0 = &rr[i0][i2];
                  pdn = &dn[i2];

                  do {
                     ps1 = ps0 + *pdn;
                     alp1 = alp0 + *rr2 + 2.0F * ( *rr1 + *rr0 );
                     sq1 = ps1 * ps1;

                     if ( ( alp * sq1 ) > ( sq * alp1 ) ) {
                        sq = sq1;
                        ps = ps1;
                        alp = alp1;
                        ix = ( Word16 )( pdn - dn );
                     }
                     pdn += STEP;
                     rr2 += ( 40 * STEP + STEP );
                     rr1 += STEP;
                     rr0 += STEP;
                  } while ( pdn < pdn_max );
                  i2 = ix;

                  if ( ( alpk * sq ) > ( psk * alp ) ) {
                     psk = sq;
                     alpk = alp;
                     codvec[0] = i0;
                     codvec[1] = i1;
                     codvec[2] = i2;
                  }
               }
            }

            pos = ipos[2];
            ipos[2] = ipos[1];
            ipos[1] = ipos[0];
            ipos[0] = pos;
         }
      }
   }
}

static void build_code_3i40_14bits( Word32 codvec[], Float32 dn_sign[], Float32 cod[], Float32 h[], Float32 y[], Word16 *anap )
{
   Float64 s;
   Float32 _sign[3];
   Float32 *p0, *p1, *p2;
   Word32 i, j, k, track, index, indx, rsign;

   memset( cod, 0, 160 );
   indx = 0;
   rsign = 0;

   for ( k = 0; k < 3; k++ ) {
      i = codvec[k];
      j = ( Word16 )dn_sign[i];
      index = i / 5;
      track = i % 5;

      if ( track == 1 )
         index = index << 4;
      else if ( track == 2 ) {
         track = 2;
         index = index << 8;
      }
      else if ( track == 3 ) {
         track = 1;
         index = ( index << 4 ) + 8;
      }
      else if ( track == 4 ) {
         track = 2;
         index = ( index << 8 ) + 128;
      }

      if ( j > 0 ) {
         cod[i] = 0.9998779296875F;
         _sign[k] = 1.0F;
         rsign = rsign + ( 1 << track );
      }
      else {
         cod[i] = -1.0F;
         _sign[k] = -1.0F;
      }
      indx = indx + index;
   }
   p0 = h - codvec[0];
   p1 = h - codvec[1];
   p2 = h - codvec[2];

   for ( i = 0; i < L_CODE; i++ ) {
      s = *p0++ * _sign[0];
      s += *p1++ * _sign[1];
      s += *p2++ * _sign[2];
      y[i] = ( Float32 )s;
   }
   anap[0] = ( Word16 )indx;
   anap[1] = ( Word16 )rsign;
}

static void code_3i40_14bits( Float32 x[], Float32 h[], Word32 T0, Float32 pitch_sharp, Float32 code[], Float32 y[], Word16 *anap )
{
   Float32 rr[L_CODE][L_CODE];
   Float32 dn[L_CODE], dn2[L_CODE], dn_sign[L_CODE];
   Word32 codvec[3];

   if ( ( T0 < L_CODE ) && ( pitch_sharp != 0 ) ) {
      for ( Word32 i = T0; i < L_CODE; i++ ) {
         h[i] = h[i] + ( h[i - T0] * pitch_sharp );
      }
   }

   cor_h_x( h, x, dn );
   set_sign( dn, dn_sign, dn2, 6 );
   cor_h( h, dn_sign, rr );
   search_3i40( dn, dn2, rr, codvec );
   build_code_3i40_14bits( codvec, dn_sign, code, h, y, anap );

   if ( ( T0 < L_CODE ) && ( pitch_sharp != 0 ) ) {
      for ( Word32 i = T0; i < L_CODE; i++ ) {
         code[i] = code[i] + ( code[i - T0] * pitch_sharp );
      }
   }
}

static void search_4i40( Float32 dn[], Float32 dn2[], Float32 rr[][L_CODE], Word32 codvec[] )
{
   Float64 alpk, alp, alp0, alp1;
   Float32 ps, psk, ps0, ps1, sq, sq1;
   Word32 ipos[4];
   Word32 i0, i1, i2, i3, ix, i, pos, track;

   psk = -1;
   alpk = 1;

   for ( i = 0; i < 4; i++ ) {
      codvec[i] = i;
   }

   for ( track = 3; track < 5; track++ ) {
      ipos[0] = 0;
      ipos[1] = 1;
      ipos[2] = 2;
      ipos[3] = track;

      for ( i = 0; i < 4; i++ ) {
         for ( i0 = ipos[0]; i0 < L_CODE; i0 += STEP ) {
            if ( dn2[i0] >= 0 ) {
               ps0 = dn[i0];
               alp0 = rr[i0][i0] * 0.25F;
               sq = -1;
               alp = 1;
               ps = 0;
               ix = ipos[1];

               for ( i1 = ipos[1]; i1 < L_CODE; i1 += STEP ) {
                  ps1 = ps0 + dn[i1];
                  alp1 = alp0 + rr[i1][i1] * 0.25F;
                  alp1 = alp1 + rr[i0][i1] * 0.5F;
                  sq1 = ps1 * ps1;

                  if ( ( alp * sq1 ) > ( sq * alp1 ) ) {
                     sq = sq1;
                     ps = ps1;
                     alp = alp1;
                     ix = i1;
                  }
               }

               i1 = ix;
               ps0 = ps;
               alp0 = alp * 0.25F;
               sq = -1;
               alp = 1;
               ps = 0;
               ix = ipos[2];

               for ( i2 = ipos[2]; i2 < L_CODE; i2 += STEP ) {
                  ps1 = ps0 + dn[i2];
                  alp1 = alp0 + rr[i2][i2] * 0.0625F;
                  alp1 += rr[i1][i2] * 0.125F;
                  alp1 += rr[i0][i2] * 0.125F;
                  sq1 = ps1 * ps1;

                  if ( ( alp * sq1 ) > ( sq * alp1 ) ) {
                     sq = sq1;
                     ps = ps1;
                     alp = alp1;
                     ix = i2;
                  }
               }

               i2 = ix;

               ps0 = ps;
               alp0 = alp;
               sq = -1;
               alp = 1;
               ps = 0;
               ix = ipos[3];

               for ( i3 = ipos[3]; i3 < L_CODE; i3 += STEP ) {
                  ps1 = ps0 + dn[i3];
                  alp1 = alp0 + rr[i3][i3] * 0.0625F;
                  alp1 += rr[i2][i3] * 0.125F;
                  alp1 += rr[i1][i3] * 0.125F;
                  alp1 += rr[i0][i3] * 0.125F;
                  sq1 = ps1 * ps1;

                  if ( ( alp * sq1 ) > ( sq * alp1 ) ) {
                     sq = sq1;
                     ps = ps1;
                     alp = alp1;
                     ix = i3;
                  }
               }

               if ( ( alpk * sq ) > ( psk * alp ) ) {
                  psk = sq;
                  alpk = alp;
                  codvec[0] = i0;
                  codvec[1] = i1;
                  codvec[2] = i2;
                  codvec[3] = ix;
               }
            }
         }

         pos = ipos[3];
         ipos[3] = ipos[2];
         ipos[2] = ipos[1];
         ipos[1] = ipos[0];
         ipos[0] = pos;
      }
   }
}

static void build_code_4i40( Word32 codvec[], Float32 dn_sign[], Float32 cod[], Float32 h[], Float32 y[], Word16 *anap )
{
   Float64 s;
   Float32 *p0, *p1, *p2, *p3;
   Word32 _sign[4];
   Word32 i, j, k, track, index, indx = 0, rsign = 0;

   memset( cod, 0, 160 );

   for ( k = 0; k < 4; k++ ) {
      i = codvec[k];
      j = ( Word16 )dn_sign[i];
      index = i / 5;
      track = i % 5;
      index = gray[index];

      if ( track == 1 )
         index = index << 3;
      else if ( track == 2 ) {
         index = index << 6;
      } else if ( track == 3 ) {
         index = index << 10;
      } else if ( track == 4 ) {
         track = 3;
         index = ( index << 10 ) + 512;
      }

      if ( j > 0 ) {
         cod[i] = 0.9998779296875F;
         _sign[k] = 1;
         rsign = rsign + ( 1 << track );
      }
      else {
         cod[i] = -1;
         _sign[k] = -1;
      }
      indx = indx + index;
   }
   p0 = h - codvec[0];
   p1 = h - codvec[1];
   p2 = h - codvec[2];
   p3 = h - codvec[3];

   for ( i = 0; i < L_CODE; i++ ) {
      s = *p0++ * _sign[0];
      s += *p1++ * _sign[1];
      s += *p2++ * _sign[2];
      s += *p3++ * _sign[3];
      y[i] = ( Float32 )( s );
   }
   anap[0] = ( Word16 )indx;
   anap[1] = ( Word16 )rsign;
}

static void code_4i40_17bits( Float32 x[], Float32 h[], Word32 T0, Float32 pitch_sharp, Float32 code[], Float32 y[], Word16 *anap )
{
   Float32 rr[L_CODE][L_CODE];
   Float32 dn[L_CODE], dn2[L_CODE], dn_sign[L_CODE];
   Word32 codvec[4];

   if ( ( T0 < L_CODE ) && ( pitch_sharp != 0 ) ) {
      for ( Word32 i = T0; i < L_CODE; i++ ) {
         h[i] = h[i] + h[i - T0] * pitch_sharp;
      }
   }

   cor_h_x( h, x, dn );
   set_sign( dn, dn_sign, dn2, 4 );
   cor_h( h, dn_sign, rr );
   search_4i40( dn, dn2, rr, codvec );
   build_code_4i40( codvec, dn_sign, code, h, y, anap );

   if ( ( T0 < L_CODE ) && ( pitch_sharp != 0 ) ) {
      for ( Word32 i = T0; i < L_CODE; i++ ) {
         code[i] = code[i] + ( code[i - T0] * pitch_sharp );
      }
   }
}

static void set_sign12k2( Float32 dn[], Float32 cn[], Float32 sign[], Word32 pos_max[], Word16 nb_track, Word32 ipos[], Word16 step )
{
   Float32 b[L_CODE];
   Float32 val, cor, k_cn, k_dn, max, max_of_all, sum;
   Word32 i, j, pos = 0;

   sum = 0.01F;
   sum += (Float32)Dotproduct40( cn, cn );
   k_cn = ( Float32 )( 1 / sqrt( sum ) );
   sum = 0.01F;
   sum += (Float32)Dotproduct40( dn, dn );
   k_dn = ( Float32 )( 1 / sqrt( sum ) );

   for ( i = 0; i < L_CODE; i++ ) {
      val = dn[i];
      cor = ( k_cn * cn[i] ) + ( k_dn * val );
      sign[i] = 1;

      if ( cor < 0 ) {
         sign[i] = -1;
         cor = -cor;
         val = -val;
      }

      dn[i] = val;
      b[i] = cor;
   }
   max_of_all = -1;

   for ( i = 0; i < nb_track; i++ ) {
      max = -1;

      for ( j = i; j < L_CODE; j += step ) {
         cor = b[j];
         val = cor - max;

         if ( val > 0 ) {
            max = cor;
            pos = j;
         }
      }

      pos_max[i] = pos;
      val = max - max_of_all;

      if ( val > 0 ) {
         max_of_all = max;
         ipos[0] = i;
      }
   }

   pos = ipos[0];
   ipos[nb_track] = pos;

   for ( i = 1; i < nb_track; i++ ) {
      pos++;

      if ( pos >= nb_track ) {
         pos = 0;
      }
      ipos[i] = pos;
      ipos[i + nb_track] = pos;
   }
}

static void search_8i40( Float32 dn[], Float32 rr[][L_CODE], Word32 ipos[], Word32 pos_max[], Word32 codvec[] )
{
   Float32 psk, ps, ps0, ps1, ps2;
   Float32 sq, sq2;
   Float32 alpk, alp, alp0, alp1, alp2;
   Word32 i0, i1, i2, i3, i4, i5, i6, i7;
   Word32 ia, ib;
   Word32 i, j, k, pos;

   psk = -1.0F;
   alpk = 1.0F;

   for ( i = 0; i < 8; i++ ) {
      codvec[i] = i;
   }

   i0 = pos_max[ipos[0]];

   for ( i = 1; i < 5; i++ ) {
      i1 = pos_max[ipos[1]];
      ps0 = dn[i0] + dn[i1];
      alp0 = rr[i0][i0] + rr[i1][i1] + 2.0F * rr[i0][i1];
      sq = -1.0F;
      alp = 1.0F;
      ps = 0.0F;
      ia = ipos[2];
      ib = ipos[3];

      for ( i2 = ipos[2]; i2 < L_CODE; i2 += 4 ) {
         ps1 = ps0 + dn[i2];
         alp1 = alp0 + rr[i2][i2] + 2.0F * ( rr[i0][i2] + rr[i1][i2] );
         for ( i3 = ipos[3]; i3 < L_CODE; i3 += 4 ) {
            ps2 = ps1 + dn[i3];
            sq2 = ps2 * ps2;
            alp2 = alp1 + rr[i3][i3] + 2.0F * ( rr[i0][i3] + rr[i1][i3] + rr[i2][i3] );
            if ( alp * sq2 > sq * alp2 ) {
               sq = sq2;
               ps = ps2;
               alp = alp2;
               ia = i2;
               ib = i3;
            }
         }
      }

      i2 = ia;
      i3 = ib;
      ps0 = ps;
      alp0 = alp;
      sq = -1.0F;
      alp = 1.0F;
      ps = 0.0F;
      ia = ipos[4];
      ib = ipos[5];

      for ( i4 = ipos[4]; i4 < L_CODE; i4 += 4 ) {
         ps1 = ps0 + dn[i4];

         alp1 = alp0 + rr[i4][i4] + 2.0F * ( rr[i0][i4] + rr[i1][i4] + rr[i2][i4] + rr[i3][i4] );

         for ( i5 = ipos[5]; i5 < L_CODE; i5 += 4 ) {

            ps2 = ps1 + dn[i5];
            sq2 = ps2 * ps2;

            alp2 = alp1 + rr[i5][i5] + 2.0F * ( rr[i0][i5] + rr[i1][i5] + rr[i2][i5] + rr[i3][i5] + rr[i4][i5] );

            if ( alp * sq2 > sq * alp2 ) {
               sq = sq2;
               ps = ps2;
               alp = alp2;
               ia = i4;
               ib = i5;
            }
         }
      }

      i4 = ia;
      i5 = ib;
      ps0 = ps;
      alp0 = alp;
      sq = -1.0F;
      alp = 1.0F;
      ps = 0.0F;

      ia = ipos[6];
      ib = ipos[7];

      for ( i6 = ipos[6]; i6 < L_CODE; i6 += 4 ) {

         ps1 = ps0 + dn[i6];

         alp1 = alp0
              + rr[i6][i6]
              + 2.0F * ( rr[i0][i6]
                        + rr[i1][i6]
                        + rr[i2][i6]
                        + rr[i3][i6]
                        + rr[i4][i6]
                        + rr[i5][i6] );

         for ( i7 = ipos[7]; i7 < L_CODE; i7 += 4 ) {

            ps2 = ps1 + dn[i7];
            sq2 = ps2 * ps2;

            alp2 = alp1
                 + rr[i7][i7]
                 + 2.0F * ( rr[i0][i7]
                           + rr[i1][i7]
                           + rr[i2][i7]
                           + rr[i3][i7]
                           + rr[i4][i7]
                           + rr[i5][i7]
                           + rr[i6][i7] );

            if ( alp * sq2 > sq * alp2 ) {
               sq = sq2;
               ps = ps2;
               alp = alp2;
               ia = i6;
               ib = i7;
            }
         }
      }

      if ( alpk * sq > psk * alp ) {
         psk = sq;
         alpk = alp;
         codvec[0] = i0;
         codvec[1] = i1;
         codvec[2] = i2;
         codvec[3] = i3;
         codvec[4] = i4;
         codvec[5] = i5;
         codvec[6] = ia;
         codvec[7] = ib;
      }

      pos = ipos[1];

      for ( j = 1, k = 2; k < 8; j++, k++ ) {
         ipos[j] = ipos[k];
      }

      ipos[7] = pos;
   }
}

static void build_code_8i40_31bits( Word32 codvec[], Float32 dn_sign[], Float32 cod[], Float32 h[], Float32 y[], Word32 sign_indx[], Word32 pos_indx[] )
{
   Float64 s;
   Float32 *p0, *p1, *p2, *p3, *p4, *p5, *p6, *p7;
   Word32 sign[8];
   Word32 i, j, k, track, sign_index, pos_index;

   memset( cod, 0, L_CODE <<2 );

   for ( i = 0; i < NB_TRACK_MR102; i++ ) {
      pos_indx[i] = -1;
      sign_indx[i] = -1;
   }

   for ( k = 0; k < 8; k++ ) {
      i = codvec[k];
      j = ( Word32 )dn_sign[i];
      pos_index = i >> 2;
      track = i & 3;

      if ( j > 0 ) {
         cod[i] = cod[i] + 0.99987792968750F;
         sign[k] = 1;
         sign_index = 0;
      } else {
         cod[i] = cod[i] - 0.99987792968750F;
         sign[k] = -1;
         sign_index = 1;
      }

      if ( pos_indx[track] < 0 ) {
         pos_indx[track] = pos_index;
         sign_indx[track] = sign_index;
      } else {
         if ( ( ( sign_index ^ sign_indx[track] ) & 1 ) == 0 ) {
            if ( pos_indx[track] <= pos_index ) {
               pos_indx[track + NB_TRACK_MR102] = pos_index;
            } else {
               pos_indx[track + NB_TRACK_MR102] = pos_indx[track];
               pos_indx[track] = pos_index;
               sign_indx[track] = sign_index;
            }
         }
         else {
            if ( pos_indx[track] <= pos_index ) {   /*swap*/
               pos_indx[track + NB_TRACK_MR102] = pos_indx[track];
               pos_indx[track] = pos_index;
               sign_indx[track] = sign_index;
            } else {
               pos_indx[track + NB_TRACK_MR102] = pos_index;
            }
         }
      }
   }

   p0 = h - codvec[0];
   p1 = h - codvec[1];
   p2 = h - codvec[2];
   p3 = h - codvec[3];
   p4 = h - codvec[4];
   p5 = h - codvec[5];
   p6 = h - codvec[6];
   p7 = h - codvec[7];

   for ( i = 0; i < L_CODE; i++ ) {
      s = 0.0F;

      for ( j = 0; j < 8; j++ ) {
         s += p[j][i] * sign[j];
      }

      y[i] = ( Float32 )s;
   }
}

static Word16 compress10( Word32 pos_indxA, Word32 pos_indxB, Word32 pos_indxC )
{
   Word32 ia = pos_indxA >> 1;
   Word32 ib = ( ( pos_indxB >> 1 ) * 5 );
   Word32 ic = ( ( pos_indxC >> 1 ) * 25 );
   Word32 indx = ( ia + ( ib + ic ) ) << 3;
   ia = pos_indxA & 1;
   ib = ( pos_indxB & 1 ) << 1;
   ic = ( pos_indxC & 1 ) << 2;
   indx = indx + ( ia + ( ib + ic ) );
   return( Word16 )indx;
}

static void compress_code( Word32 sign_indx[], Word32 pos_indx[], Word16 indx[])
{
   for ( Word32 i = 0; i < NB_TRACK_MR102; i++ ) {
      indx[i] = ( Word16 )sign_indx[i];
   }

   indx[NB_TRACK_MR102] = compress10( pos_indx[0], pos_indx[4], pos_indx[1] );
   indx[NB_TRACK_MR102 + 1] = compress10( pos_indx[2], pos_indx[6], pos_indx[5] );

   Word32 ib = ( pos_indx[7] >> 1 ) & 1;
   Word32 ia = (ib == 1) ? 4 - ( pos_indx[3] >> 1 ) : pos_indx[3] >> 1
   ib = ( ( pos_indx[7] >> 1 ) * 5 );
   ib = ( ( ia + ib ) << 5 ) + 12;
   Word32 ic = ( ( ib * 1311 ) >> 15 ) << 2;
   ia = pos_indx[3] & 1;
   ib = ( pos_indx[7] & 1 ) << 1;
   indx[NB_TRACK_MR102 + 2] = ( Word16 )( ia + ( ib + ic ) );
}

static void code_8i40_31bits( Float32 x[], Float32 cn[], Float32 h[],
                             Word32 T0, Float32 pitch_sharp, Float32 code[],
                             Float32 y[], Word16 anap[] )
{
   Float32 rr[L_CODE][L_CODE];
   Float32 dn[L_CODE], sign[L_CODE];
   Word32 ipos[8], pos_max[NB_TRACK_MR102], codvec[8], linear_signs[NB_TRACK_MR102], linear_codewords[8];
   Word32 i;

   if ( pitch_sharp > 1.0F )
      pitch_sharp = 1.0F;

   if ( pitch_sharp != 0 ) {
      for ( i = T0; i < L_SUBFR; i++ ) {
         h[i] += h[i - T0] * pitch_sharp;
      }
   }

   cor_h_x( h, x, dn );
   set_sign12k2( dn, cn, sign, pos_max, NB_TRACK_MR102, ipos, STEP_MR102 );
   cor_h( h, sign, rr );
   search_8i40( dn, rr, ipos, pos_max, codvec );
   build_code_8i40_31bits( codvec, sign, code, h, y, linear_signs, linear_codewords );
   compress_code( linear_signs, linear_codewords, anap );

   if ( pitch_sharp != 0 ) {
      for ( i = T0; i < L_SUBFR; i++ ) {
         code[i] += code[i - T0] * pitch_sharp;
      }
   }
}

static void search_pair(
      Float32 dn[],
      Float32 rr[][L_CODE],
      const Word32 selected[],
      Word32 n_selected,
      Word32 track_a,
      Word32 track_b,
      Word32 *best_a,
      Word32 *best_b,
      Float32 *best_ps,
      Float32 *best_alp )
{
   Float32 best_sq = -1.0F;
   Float32 ps_best = 0.0F;
   Float32 alp_best = 1.0F;

   for ( Word32 a = track_a; a < L_CODE; a += STEP ) {
      for ( Word32 b = track_b; b < L_CODE; b += STEP ) {
         Float32 ps = dn[a] + dn[b];
         Float32 alp = rr[a][a] + rr[b][b] + 2.0F * rr[a][b];

         for ( Word32 k = 0; k < n_selected; k++ ) {
            Word32 i = selected[k];

            ps += dn[i];

            alp += rr[i][i]
                 + 2.0F * rr[i][a]
                 + 2.0F * rr[i][b];
         }

         Float32 sq = ps * ps;

         if ( ( alp_best * sq ) > ( best_sq * alp ) ) {
            best_sq = sq;
            ps_best = ps;
            alp_best = alp;

            *best_a = a;
            *best_b = b;
         }
      }
   }

   *best_ps = ps_best;
   *best_alp = alp_best;
}

static void search_10i40(
      Float32 dn[],
      Float32 rr[][L_CODE],
      Word32 ipos[],
      Word32 pos_max[],
      Word32 codvec[] )
{
   Float32 psk = -1.0F;
   Float32 alpk = 1.0F;

   Word32 selected[10];
   Word32 i0;
   Word32 i1;
   Word32 pos;

   for ( Word32 i = 0; i < 10; i++ ) {
      codvec[i] = i;
      selected[i] = i;
   }

   i0 = pos_max[ipos[0]];

   for ( Word32 iteration = 0; iteration < 4; iteration++ ) {
      Float32 ps;
      Float32 alp;
      Word32 a;
      Word32 b;

      i1 = pos_max[ipos[1]];

      selected[0] = i0;
      selected[1] = i1;

      search_pair(dn, rr, selected, 2, ipos[2], ipos[3], &a, &b, &ps, &alp );
      selected[2] = a;
      selected[3] = b;
      search_pair(dn, rr, selected, 4, ipos[4], ipos[5], &a, &b, &ps, &alp );
      selected[4] = a;
      selected[5] = b;
      search_pair(dn, rr, selected, 6, ipos[6], ipos[7], &a, &b, &ps, &alp );
      selected[6] = a;
      selected[7] = b;
      search_pair(dn, rr, selected, 8, ipos[8], ipos[9], &a, &b, &ps, &alp );

      selected[8] = a;
      selected[9] = b;

      {
         Float32 sq = ps * ps;

         if ( ( alpk * sq ) > ( psk * alp ) ) {
            psk = sq;
            alpk = alp;

            for ( Word32 k = 0; k < 10; k++ ) {
               codvec[k] = selected[k];
            }
         }
      }

      pos = ipos[1];

      for ( Word32 k = 1; k < 9; k++ ) {
         ipos[k] = ipos[k + 1];
      }

      ipos[9] = pos;
   }
}

static void build_code_10i40_35bits( Word32 codvec[], Float32 dn_sign[], Float32 cod[], Float32 h[], Float32 y[], Word16 indx[] )
{
   Word32 i, j, k, track, index, sign[10];
   Float32 *p0, *p1, *p2, *p3, *p4, *p5, *p6, *p7, *p8, *p9;
   Float64 s;

   memset( cod, 0, 160 );
   memset( y, 0, 160 );

   for ( i = 0; i < NB_TRACK; i++ ) {
      indx[i] = -1;
   }

   for ( k = 0; k < 10; k++ ) {
      i = codvec[k];
      j = ( Word16 )dn_sign[i];
      index = ( Word16 )( i / 5 );
      track = ( Word16 )( i % 5 );

      if ( j > 0 ) {
         cod[i] = cod[i] + 1;
         sign[k] = 1;
      } else {
         cod[i] = cod[i] - 1;
         sign[k] = -1;
         index = index + 8;
      }

      if ( indx[track] < 0 ) {
         indx[track] = ( Word16 )index;
      } else {
         if ( ( ( index ^ indx[track] ) & 8 ) == 0 ) {
            if ( indx[track] <= index ) {
               indx[track + 5] = ( Word16 )index;
            } else {
               indx[track + 5] = ( Word16 )indx[track];
               indx[track] = ( Word16 )index;
            }
         } else {
            if ( ( indx[track] & 7 ) <= ( index & 7 ) ) {
               indx[track + 5] = ( Word16 )indx[track];
               indx[track] = ( Word16 )index;
            } else {
               indx[track + 5] = ( Word16 )index;
            }
         }
      }
   }
   p0 = h - codvec[0];
   p1 = h - codvec[1];
   p2 = h - codvec[2];
   p3 = h - codvec[3];
   p4 = h - codvec[4];
   p5 = h - codvec[5];
   p6 = h - codvec[6];
   p7 = h - codvec[7];
   p8 = h - codvec[8];
   p9 = h - codvec[9];

   for ( i = 0; i < L_CODE; i++ ) {
      s = 0.0F;

      for ( j = 0; j < 10; j++ ) {
         s += p[j][i] * sign[j];
      }

      y[i] = s;
   }
}

static void q_p( Word16 *ind, Word32 n )
{
   Word16 tmp = *ind;
   *ind = n < 5
      ? ( Word16 )( ( tmp & 0x8 ) | gray[tmp & 0x7] )
      : gray[tmp & 7];
}

static void code_10i40_35bits( Float32 x[], Float32 cn[], Float32 h[],
    Word32 T0, Float32 gain_pit, Float32 code[],
    Float32 y[], Word16 anap[] )
{
   Float32 rr[L_CODE][L_CODE];
   Float32 dn[L_CODE], sign[L_CODE];
   Word32 ipos[10], pos_max[NB_TRACK], codvec[10];

   if ( gain_pit > 1.0F )
      gain_pit = 1.0F;

   if ( gain_pit != 0 ) {
      for ( Word32 i = T0; i < L_SUBFR; i++ ) {
         h[i] += h[i - T0] * gain_pit;
      }
   }
   
   cor_h_x( h, x, dn );
   set_sign12k2( dn, cn, sign, pos_max, NB_TRACK, ipos, STEP );
   cor_h( h, sign, rr );
   search_10i40( dn,  rr, ipos, pos_max, codvec );
   build_code_10i40_35bits( codvec, sign, code, h, y, anap );

   for ( Word32 i = 0; i < 10; i++ ) {
      q_p( &anap[i], i );
   }

   if ( gain_pit != 0 ) {
      for ( Word32 i = T0; i < L_SUBFR; i++ ) {
         code[i] += code[i - T0] * gain_pit;
      }
   }
}

static void cbsearch( enum ModeNB mode, Word16 subnr, Float32 x[],
                     Float32 h[], Word32 T0, Float32 pitch_sharp,
                     Float32 gain_pit, Float32 code[], Float32 y[],
                     Float32 *res2, Word16 **anap )
{
   switch (mode){
   case MR475:
   case MR515:
      code_2i40_9bits( subnr, x, h, T0, pitch_sharp, code, y, *anap );
      ( *anap ) += 2;
      break;
   case MR59:
      code_2i40_11bits( x, h, T0, pitch_sharp, code, y, *anap );
      ( *anap ) += 2;
      break;
   case MR67:
      code_3i40_14bits( x, h, T0, pitch_sharp, code, y, *anap );
      ( *anap ) += 2;
      break;
   case MR74:
   case MR795:
      code_4i40_17bits( x, h, T0, pitch_sharp, code, y, *anap );
      ( *anap ) += 2;
      break;
   case MR102:
      code_8i40_31bits( x, res2, h, T0, pitch_sharp, code, y, *anap );
      *anap += 7;
      break;
   default:
      code_10i40_35bits( x, res2, h, T0, gain_pit, code, y, *anap );
      *anap += 10;
   }
}

static void Log2_norm( Word32 x, Word32 exp, Word32 *exponent, Word32 *fraction )
{
   if ( x <= 0 ) {
      *exponent = 0;
      *fraction = 0;
      return;
   }

   Word32 i = x >> 25;
   i = i - 32;
   Word32 a = x >> 9;
   a = a & 0xFFFE;   /* 2a */
   Word32 y = ( log2_table[i] << 16 ) - a * ( log2_table[i] - log2_table[i + 1] );
   *fraction = y >> 16;
   *exponent = 30 - exp;
}

static void Log2( Word32 x, Word32 *exponent, Word32 *fraction )
{
   int exp;
   frexp( ( Float64 )x, &exp );
   exp = 31 - exp;
   Log2_norm( x <<exp, exp, exponent, fraction );
}

static Word32 Pow2( Word32 exponent, Word32 fraction )
{
   Word32 i, a, tmp, x, exp;
   i = fraction >> 10;
   a = ( fraction << 5 ) & 0x7fff;
   x = pow2_table[i] << 16;
   tmp = pow2_table[i] - pow2_table[i + 1];
   x -= ( tmp * a ) << 1;

   if ( exponent >= -1 ) {
      exp = ( 30 - exponent );

      if ( ( x & ( ( Word32 )1 << ( exp - 1 ) ) ) != 0 ) {
         return ( x >> exp ) + 1;
      } else
         return x >> exp;
   }
   else
      return 0;
}

static void gc_pred( Word32 *past_qua_en, enum ModeNB mode, Float32 *code, Word32 *gcode0_exp, Word32 *gcode0_fra, Float32 *en )
{
   Float64 ener_code;
   Word32 exp, frac, ener, ener_tmp, tmp;
   int exp_code;

   ener_code = Dotproduct40( code, code );

   if ( mode == MR122 ) {
      ener = (Word32)(ener_code * 33554432);
      ener = ( ( ener + 0x00008000L ) >> 16 ) * 52428;

      Log2( ener, &exp, &frac );
      ener = ( ( exp - 30 ) << 16 ) + ( frac << 1 );

      ener_tmp = 44 * qua_gain_code_MR122[past_qua_en[0]];
      ener_tmp += 37 * qua_gain_code_MR122[past_qua_en[1]];
      ener_tmp += 22 * qua_gain_code_MR122[past_qua_en[2]];
      ener_tmp += 12 * qua_gain_code_MR122[past_qua_en[3]];

      ener_tmp = ener_tmp << 1;
      ener_tmp += 783741L;

      ener = ( ener_tmp - ener ) >> 1;
      *gcode0_exp = ener >> 16;
      *gcode0_fra = ( ener >> 1 ) - ( *gcode0_exp << 15 );
   }
   else {
      ener = (Word32)(ener_code * 134217728);
      if (ener < 0)
         ener = 0x7fffffff;

      frexp( ( Float64 )ener, &exp_code );
      exp_code = 31 - exp_code;
      ener <<= exp_code;

      Log2_norm( ener, exp_code, &exp, &frac );

      tmp = ( exp * ( -49320 ) ) + ( ( ( frac * ( -24660 ) ) >> 15 ) << 1 );

      if ( mode == MR102 ) {
         tmp += 2134784;
      } else if ( mode == MR795 ) {
         tmp += 2183936;
         *en = (Float32)ener_code;
      } else if ( mode == MR74 ) {
         tmp += 2085632;
      } else if ( mode == MR67 ) {
         tmp += 2065152;   /* Q14 */
      } else /* MR59, MR515, MR475 */ {
         tmp += 2134784;
      }

      tmp = tmp << 9;

      tmp += 5571 * qua_gain_code[past_qua_en[0]];
      tmp += 4751 * qua_gain_code[past_qua_en[1]];
      tmp += 2785 * qua_gain_code[past_qua_en[2]];
      tmp += 1556 * qua_gain_code[past_qua_en[3]];

      tmp = tmp >> 15;   /* Q8  */

      if ( mode == MR74 ) {
         tmp = tmp * 10878;
      }
      else {
         tmp = tmp * 10886;
      }
      tmp = tmp >> 9;   /* -> Q15 */

      *gcode0_exp = tmp >> 15;
      *gcode0_fra = tmp - ( *gcode0_exp * 32768 );
   }
}

static void calc_filt_energies( enum ModeNB mode, Float32 xn[], Float32 xn2[],
      Float32 y1[], Float32 y2[], Float32 gCoeff[], Float32 coeff[], Float32 *
      cod_gain )
{
   Float32 sum, ener_init = 0.01F;

   if ( ( mode == MR795 ) || ( mode == MR475 ) )
      ener_init = 0;
   coeff[0] = gCoeff[0];
   coeff[1] = -2.0F * gCoeff[1];

   sum = (Float32)Dotproduct40( y2, y2 );
   sum += ener_init;
   coeff[2] = sum;

   sum = (Float32)Dotproduct40( xn, y2 );
   sum += ener_init;
   coeff[3] = -2.0F * sum;

   sum = (Float32)Dotproduct40( y1, y2 );
   sum += ener_init;
   coeff[4] = 2.0F * sum;

   if ( ( mode == MR475 ) || ( mode == MR795 ) ) {
      sum = (Float32)Dotproduct40( xn2, y2 );
      *cod_gain = sum <= 0 ? 0 : sum / coeff[2];
   }
}

static void MR475_update_unq_pred( Word32 *past_qua_en, Float32 gcode0, Float32 cod_gain )
{
   Float32 qua_ener, pred_err_fact;
   Word32 i, index, energy, max, s;

   if ( cod_gain <= 0 ) {
      qua_ener = -32.0F;
   }
   else {
      pred_err_fact = gcode0 != 0 ? cod_gain / gcode0 : 10f;
      qua_ener = pred_err_fact < 0.0251189F ? -32.0F
         : pred_err_fact > 7.8125F ? 17.8558F
         : ( Float32 )( 20.0F*log10( pred_err_fact ) );
   }
   energy = (Word32)(qua_ener * 1024 + 0.5F);
   max = abs(energy - qua_gain_code[0]);
   index = 0;

   for ( i = 1; i < NB_QUA_CODE+VQ_SIZE_HIGHRATES+VQ_SIZE_LOWRATES+MR475_VQ_SIZE*2+3; i++ )
   {
      s = abs(energy - qua_gain_code[i]);
      if (s < max){
         max = s;
         index = i;
         if (s == 0) {
            break;
         }
      }
   }

   for ( i = 3; i > 0; i-- ) {
      past_qua_en[i] = past_qua_en[i - 1];
   }
   past_qua_en[0] = index;
}

static Word16 MR475_gain_quant( Word32 *past_qua_en, Word32 sf0_gcode0_exp, Word32
                               sf0_gcode0_fra, Float32 sf0_coeff[], Float32 sf0_target_en,
                               Float32 sf1_code_nosharp[], Word32 sf1_gcode0_exp, Word32
                               sf1_gcode0_fra, Float32 sf1_coeff[], Float32 sf1_target_en,
                               Float32 gp_limit, Float32 *sf0_gain_pit, Float32
                               *sf0_gain_cod, Float32 *sf1_gain_pit, Float32 *sf1_gain_cod )
{
   Float32 temp, temp2, g_pitch, g2_pitch, g_code, g2_code, g_pit_cod, dist_min, sf0_gcode0, sf1_gcode0;
   const Float32 *p;
   Word32 i, tmp, g_code_tmp, gcode0, index = 0;

   sf0_gcode0 = (Float32)Pow2(sf0_gcode0_exp, sf0_gcode0_fra);
   sf1_gcode0 = (Float32)Pow2(sf1_gcode0_exp, sf1_gcode0_fra);

   if ( ( sf0_target_en * 2.0F ) < sf1_target_en ) {
      for ( i = 0; i < 5; i++ )
         sf0_coeff[i] *= 2f;
   }
   else if ( sf0_target_en > ( sf1_target_en * 4.0F ) ) {
      for ( i = 0; i < 5; i++ )
         sf1_coeff[i] *= 2f;
   }

   dist_min = FLT_MAX;
   p = &table_gain_MR475[0];

   for ( i = 0; i < MR475_VQ_SIZE; i++ ) {
      g_pitch = *p++;
      g_code = *p++;
      g_code *= sf0_gcode0;
      g2_pitch = g_pitch * g_pitch;
      g2_code = g_code * g_code;
      g_pit_cod = g_code * g_pitch;
      temp = sf0_coeff[0] * g2_pitch;
      temp += sf0_coeff[1] * g_pitch;
      temp += sf0_coeff[2] * g2_code;
      temp += sf0_coeff[3] * g_code;
      temp += sf0_coeff[4] * g_pit_cod;
      temp2 = g_pitch - gp_limit;
      g_pitch = *p++;
      g_code = *p++;

      if ( temp2 <= 0 && ( g_pitch <= gp_limit ) ) {
         g_code *= sf1_gcode0;
         g2_pitch = g_pitch * g_pitch;
         g2_code = g_code * g_code;
         g_pit_cod = g_code * g_pitch;
         temp += sf1_coeff[0] * g2_pitch;
         temp += sf1_coeff[1] * g_pitch;
         temp += sf1_coeff[2] * g2_code;
         temp += sf1_coeff[3] * g_code;
         temp += sf1_coeff[4] * g_pit_cod;

         if ( temp < dist_min ) {
            dist_min = temp;
            index = i;
         }
      }
   }

   tmp = index << 2;
   p = &table_gain_MR475[tmp];
   *sf0_gain_pit = *p++;
   g_code_tmp = (Word32)(*p++ * 4096 + 0.5F);

   gcode0 = Pow2( 14, sf0_gcode0_fra );
   if ( sf0_gcode0_exp < 11 ) {
      *sf0_gain_cod = (Float32)(( g_code_tmp * gcode0 ) >> ( 25 - sf0_gcode0_exp ));
   } else {
      i = ( ( g_code_tmp * gcode0 ) << ( sf0_gcode0_exp - 9 ) );
      *sf0_gain_cod = ( i >> ( sf0_gcode0_exp - 9 ) ) != ( g_code_tmp * gcode0 ) ? (Float32)0x7FFF : (Float32)(i >> 16);
   }

   *sf0_gain_cod *= 0.5F;

   for ( i = 3; i > 0; i-- ) {
      past_qua_en[i] = past_qua_en[i - 1];
   }
   past_qua_en[0] = NB_QUA_CODE+VQ_SIZE_HIGHRATES+VQ_SIZE_LOWRATES + (index << 1);

   gc_pred( past_qua_en, MR475, sf1_code_nosharp, &sf1_gcode0_exp, &sf1_gcode0_fra, &sf0_gcode0 );

   tmp += 2;
   p = &table_gain_MR475[tmp];
   *sf1_gain_pit = *p++;
   g_code_tmp = (Word32)(*p++ * 4096 + 0.5F);

   gcode0 = Pow2( 14, sf1_gcode0_fra );
   if ( sf1_gcode0_exp < 11 ) {
      *sf1_gain_cod = (Float32)(( g_code_tmp * gcode0 ) >> ( 25 - sf1_gcode0_exp ));
   }
   else {
      i = ( ( g_code_tmp * gcode0 ) << ( sf1_gcode0_exp - 9 ) );
      *sf1_gain_cod = ( i >> ( sf1_gcode0_exp - 9 ) ) != ( g_code_tmp * gcode0 ) ? (Float32)0x7FFF : (Float32)(i >> 16);
   }

   *sf1_gain_cod *= 0.5F;

   for ( i = 3; i > 0; i-- ) {
      past_qua_en[i] = past_qua_en[i - 1];
   }
   past_qua_en[0] = NB_QUA_CODE+VQ_SIZE_HIGHRATES+VQ_SIZE_LOWRATES + (index << 1) + 1;

   return( Word16 )index;
}

static Word16 q_gain_code( Float32 gcode0, Float32 *gain, Word32 *qua_ener_index)
{
   Float64 err_min, err;
   const Float32 *p;
   Word32 i, index;


   p = &gain_factor[0];

   err_min = fabs( *gain - ( gcode0 * *p++ ) );
   index = 0;

   for ( i = 1; i < NB_QUA_CODE; i++ ) {
      err = fabs( *gain - ( gcode0 * *p++ ) );

      if ( err < err_min ) {
         err_min = err;
         index = i;
      }
   }
   p = &gain_factor[index];
   *gain = (Float32)floor(gcode0 * *p);
   *qua_ener_index = index;

   return( Word16 )index;
}

static void MR795_gain_code_quant3( Word32 gcode0_exp, Word32 gcode0_fra, Float32 g_pitch_cand[],
      Word32 g_pitch_cind[], Float32 coeff[], Float32 *gain_pit, Word32 *
      gain_pit_ind, Float32 *gain_cod, Word32 *gain_cod_ind, Word32 *qua_ener_index )
{
   Float32 gcode0, dist_min, g_pitch, g2_pitch, g_code, g2_code, g_pit_cod, tmp0, tmp;
   const Float32 *p;
   Word32 i, j, cod_ind, pit_ind, g_code0, g_code_tmp;

   gcode0 = (Float32)Pow2( gcode0_exp, gcode0_fra);
   
   dist_min = FLT_MAX;
   cod_ind = 0;
   pit_ind = 0;

   for ( j = 0; j < 3; j++ ) {
      g_pitch = g_pitch_cand[j];
      g2_pitch = g_pitch * g_pitch;
      tmp0 = coeff[0] * g2_pitch;
      tmp0 += coeff[1] * g_pitch;
      p = &gain_factor[0];

      for ( i = 0; i < NB_QUA_CODE; i++ ) {
         g_code = *p++;
         g_code = g_code * gcode0;
         g2_code = g_code * g_code;
         g_pit_cod = g_code * g_pitch;
         tmp = tmp0 + coeff[2] * g2_code;
         tmp += coeff[3] * g_code;
         tmp += coeff[4] * g_pit_cod;

         if ( tmp < dist_min ) {
            dist_min = tmp;
            cod_ind = i;
            pit_ind = j;
         }
      }
   }

   p = &gain_factor[cod_ind];
   g_code_tmp = (Word32)(2048 * *p);
   *qua_ener_index = cod_ind;

   g_code0 = Pow2( 14, gcode0_fra);
   i = ( g_code_tmp * g_code0 ) << 1;
   gcode0_exp = 9 - gcode0_exp;
   i = gcode0_exp > 0 ? i >> gcode0_exp : i << -gcode0_exp;
   *gain_cod = (Float32)(i >> 16);
   if (*gain_cod > 32767)
      *gain_cod = 32767;

   *gain_cod *= 0.5F;

   *gain_cod_ind = ( Word16 )cod_ind;
   *gain_pit = g_pitch_cand[pit_ind];
   *gain_pit_ind = g_pitch_cind[pit_ind];
}

static void calc_unfilt_energies( Float32 res[], Float32 exc[], Float32 code[], Float32 gain_pit, Float32 en[], Float32 *ltpg )
{
   Float32 sum, pred_gain;
   Word32 i;

   en[0] = en[0] < 200 ? 0 : (Float32)Dotproduct40( res, res );
   en[1] = (Float32)Dotproduct40( exc, exc );
   en[2] = (Float32)Dotproduct40( exc, code );
   en[3] = 0;

   for ( i = 0; i < L_SUBFR; i++ ) {
      sum = res[i] - ( exc[i] * gain_pit );
      en[3] += sum * sum;
   }

   if ( en[3] > 0 && en[0] != 0 ) {
      pred_gain = en[0] / en[3];
      *ltpg = ( Float32 )( log10( pred_gain ) / log10( 2 ) );
   } else {
      *ltpg = 0;
   }
}

static Float32 gmed_n_f( Float32 ind[], Word16 n )
{
   Word32 medianIndex;
   Word32 i, j, ix = 0;
   Word32 tmp[9];
   Float32 tmp2[9];
   Float32 max;

   for ( i = 0; i < n; i++ ) {
      tmp2[i] = ind[i];
   }

   for ( i = 0; i < n; i++ ) {
      max = -FLT_MAX;

      for ( j = 0; j < n; j++ ) {
         if ( tmp2[j] >= max ) {
            max = tmp2[j];
            ix = j;
         }
      }
      tmp2[ix] = -FLT_MAX;
      tmp[i] = ix;
   }
   medianIndex = tmp[n >> 1];
   return( ind[medianIndex] );
}

static void gain_adapt( Float32 *prev_gc, Word16 *onset, Float32 *ltpg_mem, Float32 *prev_alpha, Float32 ltpg, Float32 gain_cod, Float32 *alpha )
{
   Float32 result, filt;   /* alpha factor, median-filtered LTP coding gain */
   Word32 i;
   Word16 adapt;

   if ( ltpg <= 0.3321928F /*LTP_GAIN_THR1*/ ) {
      adapt = 0;
   } else {
      if ( ltpg <= 0.6643856 /*LTP_GAIN_THR2*/ ) {
         adapt = 1;
      } else {
         adapt = 2;
      }
   }

   if ( ( gain_cod > 2.0F * *prev_gc ) && ( gain_cod > 100 ) ) {
      *onset = 8;
   } else {
      if ( *onset != 0 ) {
         (*onset)--;
      }
   }

   if ( ( *onset != 0 ) && ( adapt < 2 ) ) {
      adapt++;
   }
   ltpg_mem[0] = ltpg;
   filt = gmed_n_f( ltpg_mem, 5 );

   if ( adapt == 0 ) {
      if ( filt > 0.66443 ) {
         result = 0;
      } else {
         result = filt < 0 ? 0.5 : ( Float32 )( 0.5-0.75257499*filt );
      }
   } else {
      result = 0;
   }

   if ( *prev_alpha == 0 ) {
      result = 0.5F * result;
   }

   *alpha = result;
   *prev_alpha = result;
   *prev_gc = gain_cod;

   for ( i = LTPG_MEM_SIZE - 1; i > 0; i-- ) {
      ltpg_mem[i] = ltpg_mem[i - 1];
   }
}

static Word16 MR795_gain_code_quant_mod( Float32 gain_pit, Word32 gcode0_exp, Word32 gcode0_fra,
      Float32 en[], Float32 alpha, Float32 gain_cod_unq, Float32 *gain_cod,
      Word32 *qua_ener_index )
{
   Float32 coeff[5];
   Float32 gcode0, g2_pitch, g_code, g2_code, d2_code, dist_min, gain_code, tmp;
   const Float32 *p;
   Word32 i, index, g_code_tmp, g_code0;

   gcode0 = (Float32)Pow2(gcode0_exp, gcode0_fra);

   gain_code = *gain_cod;
   g2_pitch = gain_pit * gain_pit;
   coeff[0] = ( Float32 )( sqrt( alpha * en[0] ) );
   coeff[1] = alpha * en[1] * g2_pitch;
   coeff[2] = 2.0F * alpha * en[2] * gain_pit;
   coeff[3] = alpha * en[3];
   coeff[4] = ( 1.0F - alpha ) * en[3];

   dist_min = FLT_MAX;
   index = 0;
   p = &gain_factor[0];

   for ( i = 0; i < NB_QUA_CODE; i++ ) {
      g_code = *p++;
      g_code = g_code * gcode0;
      if ( g_code >= ( 2.0F * gain_code ) )
         break;
      g2_code = g_code * g_code;
      d2_code = g_code - gain_cod_unq;
      d2_code = d2_code * d2_code;
      tmp = coeff[1] + coeff[2] * g_code;
      tmp += coeff[3] * g2_code;
      tmp = ( Float32 )sqrt( tmp );
      tmp = tmp - coeff[0];
      tmp = tmp * tmp;
      tmp += coeff[4] * d2_code;
      if ( tmp < dist_min ) {
         dist_min = tmp;
         index = i;
      }
   }

   p = &gain_factor[index];
   g_code_tmp = (Word32)(2048 * *p);
   *qua_ener_index = index;

   g_code0 = Pow2( 14, gcode0_fra);
   i = ( g_code_tmp * g_code0 ) << 1;
   gcode0_exp = 9 - gcode0_exp;
   i = gcode0_exp > 0 ? i >> gcode0_exp : i << ( -gcode0_exp );

   *gain_cod = (Float32)(i >> 16);
   if (*gain_cod > 32767)
      *gain_cod = 32767;

   *gain_cod *= 0.5F;
   return( Word16 )index;
}

static void MR795_gain_quant( Float32 *prev_gc, Word16 *onset, Float32 *ltpg_mem
      , Float32 *prev_alpha, Float32 res[], Float32 exc[], Float32 code[],
      Float32 coeff[], Float32 code_en, Word32 gcode0_exp, Word32 gcode0_fra, Float32 cod_gain,
      Float32 gp_limit, Float32 *gain_pit, Float32 *gain_cod, Word32 *qua_ener_index,
      Word16 **anap )
{
   Float32 en[4], g_pitch_cand[3];
   Float32 ltpg, alpha, gain_cod_unq;   /* code gain (unq.) */
   Word32 g_pitch_cind[3];   /* pitch gain indices */
   Word32 gain_pit_index, gain_cod_index;

   gain_pit_index = q_gain_pitch( MR795, gp_limit, gain_pit, g_pitch_cand, g_pitch_cind );
   MR795_gain_code_quant3( gcode0_exp, gcode0_fra, g_pitch_cand, g_pitch_cind, coeff, gain_pit,
         &gain_pit_index, gain_cod, &gain_cod_index, qua_ener_index );

   calc_unfilt_energies( res, exc, code, *gain_pit, en, &ltpg );
   gain_adapt( prev_gc, onset, ltpg_mem, prev_alpha, ltpg, *gain_cod, &alpha );
   if ( ( en[0] != 0 ) && ( alpha > 0 ) ) {
      en[3] = code_en;
      gain_cod_unq = cod_gain;
      gain_cod_index = MR795_gain_code_quant_mod( *gain_pit, gcode0_exp, gcode0_fra, en, alpha, gain_cod_unq, gain_cod, qua_ener_index );
   }
   *( *anap )++ = ( Word16 )gain_pit_index;
   *( *anap )++ = ( Word16 )gain_cod_index;
}

static Word16 Qua_gain( enum ModeNB mode, Word32 gcode0_exp, Word32 gcode0_fra, Float32 coeff[], Float32
      gp_limit, Float32 *gain_pit, Float32 *gain_cod, Word32 *qua_ener_index)
{
   Float32 g_pitch, g2_pitch, g_code, g2_code, g_pit_cod, tmp, dist_min, gcode0;
   const Float32 *table_gain, *p;
   Word32 i, index = 0, gcode_0, g_code_tmp;
   Word16 table_len;

   gcode0 = (Float32)Pow2( gcode0_exp, gcode0_fra );

   if ( ( mode == MR102 ) || ( mode == MR74 ) || ( mode == MR67 ) ) {
      table_len = VQ_SIZE_HIGHRATES;
      table_gain = table_highrates;
      *qua_ener_index = NB_QUA_CODE;
   }
   else {
      table_len = VQ_SIZE_LOWRATES;
      table_gain = table_lowrates;
      *qua_ener_index = NB_QUA_CODE + VQ_SIZE_HIGHRATES;
   }

   dist_min = FLT_MAX;
   p = &table_gain[0];

   for ( i = 0; i < table_len; i++ ) {
      g_pitch = *p++;
      g_code = *p++;

      if ( g_pitch <= gp_limit ) {
         g_code *= gcode0;
         g2_pitch = g_pitch * g_pitch;
         g2_code = g_code * g_code;
         g_pit_cod = g_code * g_pitch;
         tmp = coeff[0] * g2_pitch;
         tmp += coeff[1] * g_pitch;
         tmp += coeff[2] * g2_code;
         tmp += coeff[3] * g_code;
         tmp += coeff[4] * g_pit_cod;

         if ( tmp < dist_min ) {
            dist_min = tmp;
            index = i;
         }
      }
   }

   p = &table_gain[index << 1];
   *gain_pit = *p++;
   g_code_tmp = (Word32)(4096 * *p);

   gcode_0 = Pow2( 14, gcode0_fra );
   if ( gcode0_exp < 11 ) {
      *gain_cod = (Float32)((g_code_tmp * gcode_0) >> ( 25 - gcode0_exp ));
   }
   else {
      i = ( ( g_code_tmp * gcode_0) << ( gcode0_exp - 9 ) );
      *gain_cod = ( i >> ( gcode0_exp - 9 ) ) != ( g_code_tmp * gcode_0) ? 0x7FFF : (Float32)(i >> 16);
   }
   *gain_cod = *gain_cod * 0.5F;
   *qua_ener_index += index;

   return( Word16 )index;
}

static void gainQuant( enum ModeNB mode, Word32 even_subframe, Word32 *
      past_qua_en, Word32 *past_qua_en_unq, Float32 *sf0_coeff, Float32 *
      sf0_target_en, Word32 *sf0_gcode0_exp, Word32 *sf0_gcode0_fra,Word16 **gain_idx_ptr, Float32 *
      sf0_gain_pit, Float32 *sf0_gain_cod, Float32 *res, Float32 *exc, Float32
      code[], Float32 xn[], Float32 xn2[], Float32 y1[], Float32 y2[], Float32
      gCoeff[], Float32 gp_limit, Float32 *gain_pit, Float32 *gain_cod, Float32
      *prev_gc, Word16 *onset, Float32 *ltpg_mem, Float32 *prev_alpha, Word16 **
      anap )
{
   Float32 coeff[5];
   Float32 gcode0, cod_gain, en = 0;
   Word32 i, exp, frac, qua_ener_index;

   if ( mode == MR475 ) {
      if ( even_subframe != 0 ) {
         *gain_idx_ptr = ( *anap )++;
         past_qua_en_unq[0] = past_qua_en[0];
         past_qua_en_unq[1] = past_qua_en[1];
         past_qua_en_unq[2] = past_qua_en[2];
         past_qua_en_unq[3] = past_qua_en[3];

         gc_pred( past_qua_en, mode, code, sf0_gcode0_exp, sf0_gcode0_fra, &en );
         gcode0 = (Float32)Pow2(*sf0_gcode0_exp, *sf0_gcode0_fra);

         calc_filt_energies( mode, xn, xn2, y1, y2, gCoeff, sf0_coeff, &cod_gain);
         *gain_cod = cod_gain;
         *sf0_target_en = (Float32)Dotproduct40( xn, xn );

         MR475_update_unq_pred( past_qua_en_unq, gcode0, cod_gain );
      } else {
         gc_pred( past_qua_en_unq, mode, code, &exp, &frac, &en );

         calc_filt_energies( mode, xn, xn2, y1, y2, gCoeff, coeff, &cod_gain );
         en = (Float32)Dotproduct40( xn, xn );

         **gain_idx_ptr = MR475_gain_quant( past_qua_en, *sf0_gcode0_exp, *sf0_gcode0_fra, sf0_coeff,
               *sf0_target_en, code, exp, frac, coeff, en, gp_limit, sf0_gain_pit,
               sf0_gain_cod, gain_pit, gain_cod );
      }
   } else {
      gc_pred( past_qua_en, mode, code, &exp, &frac, &en );

      if ( mode == MR122 ) {
         gcode0 = (Float32)Pow2( exp, frac );
         if (gcode0 > 2047.9375F) gcode0 = 2047.9375F;

         *gain_cod = (Float32)(Dotproduct40( xn2, y2 ) / ( Dotproduct40( y2, y2 )+ 0.01F ));

         if ( *gain_cod < 0 )
            *gain_cod = 0.0F;
         *( *anap )++ = q_gain_code( gcode0, gain_cod,&qua_ener_index);
      } else {
         calc_filt_energies( mode, xn, xn2, y1, y2, gCoeff, coeff, &cod_gain );

         if ( mode == MR795 ) {
            MR795_gain_quant( prev_gc, onset, ltpg_mem, prev_alpha, res, exc,
                  code, coeff, en, exp, frac , cod_gain, gp_limit, gain_pit,
                  gain_cod, &qua_ener_index, anap );
         } else {
            *( *anap )++ = Qua_gain( mode, exp, frac, coeff, gp_limit, gain_pit,
                  gain_cod, &qua_ener_index);
         }
      }

      for ( i = 3; i > 0; i-- ) {
         past_qua_en[i] = past_qua_en[i - 1];
      }
      past_qua_en[0] = qua_ener_index;
   }
}

static void subframePostProc( Float32 *speech, Word16 i_subfr, Float32 gain_pit,
      Float32 gain_code, Float32 *a_q, Float32 synth[], Float32 xn[], Float32
      code[], Float32 y1[], Float32 y2[], Float32 *mem_syn, Float32 *mem_err,
      Float32 *mem_w0, Float32 *exc, Float32 *sharp )
{
   Word32 i, j;

   *sharp = gain_pit;
   if ( *sharp > 0.794556F ) {
      *sharp = 0.794556F;
   }

   for ( i = 0; i < L_SUBFR; i += 4 ) {
      exc[i + i_subfr] = (Float32)floor((gain_pit * exc[i + i_subfr] + gain_code * code[i]) + 0.5F);
      exc[i + i_subfr + 1] = (Float32)floor((gain_pit * exc[i + i_subfr + 1] + gain_code * code[i + 1]) + 0.5F);
      exc[i + i_subfr + 2] = (Float32)floor((gain_pit * exc[i + i_subfr + 2] + gain_code * code[i + 2]) + 0.5F);
      exc[i + i_subfr + 3] = (Float32)floor((gain_pit * exc[i + i_subfr + 3] + gain_code * code[i + 3]) + 0.5F);
   }

   Syn_filt( a_q, &exc[i_subfr], &synth[i_subfr], mem_syn, 1 );

   for ( i = L_SUBFR - M, j = 0; i < L_SUBFR; i++, j++ ) {
      mem_err[j] = speech[i_subfr + i] - synth[i_subfr + i];
      mem_w0[j] = xn[i] - y1[i] * gain_pit - y2[i] * gain_code;
   }
}

static void Convolve( Float32 x[], Float32 h[], Float32 y[] )
{
   for ( Word32 n = 0; n < L_SUBFR; n++ ) {
      Float32 s = 0.0F;

      for ( Word32 i = 0; i <= n; i++ ) {
         s += x[i] * h[n - i];
      }
      y[n] = s;
   }
}

static Word16 tx_dtx_handler( Word16 vad_flag, Word16 *decAnaElapsedCount, Word16 *dtxHangoverCount, enum ModeNB *used_mode )
{
   Word16 compute_new_sid_possible;

   *decAnaElapsedCount += 1;
   compute_new_sid_possible = 0;

   if ( vad_flag != 0 ) {
      *dtxHangoverCount = DTX_HANG_CONST;
   } else {
      if ( *dtxHangoverCount == 0 ) {
         *decAnaElapsedCount = 0;
         *used_mode = MRDTX;
         compute_new_sid_possible = 1;
      } else {
         *dtxHangoverCount -= 1;
         if ( ( *decAnaElapsedCount + *dtxHangoverCount ) < DTX_ELAPSED_FRAMES_THRESH ) {
            *used_mode = MRDTX;
         }
      }
   }
   return compute_new_sid_possible;
}

static void dtx_buffer( Word16 *hist_ptr, Float32 *lsp_hist, Float32 lsp_new[], Float32 speech[], Float32 *log_en_hist )
{
   Float64 frame_en;

   *hist_ptr += 1;

   if ( *hist_ptr == DTX_HIST_SIZE ) {
      *hist_ptr = 0;
   }

   memcpy( &lsp_hist[ * hist_ptr * M], lsp_new, sizeof( Float32 )*M );
   frame_en = Dotproduct40( speech, speech );
   frame_en += Dotproduct40( &speech[40], &speech[40] );
   frame_en += Dotproduct40( &speech[80], &speech[80] );
   frame_en += Dotproduct40( &speech[120], &speech[120] );

   if ( frame_en > 1 ) {
      log_en_hist[ * hist_ptr] = ( Float32 )( log10( frame_en * 0.00625F )*
            1.660964F );
   } else {
      log_en_hist[ * hist_ptr] = -3.660965F;
   }
}

static Word32 dtx_enc( Word16 *log_en_index, Float32 log_en_hist[], Float32
      lsp_hist[], Word16 *lsp_index, Word32 *init_lsf_vq_index, Word16
      compute_sid_flag, Float32 past_rq[], Word32 *past_qua_en, Word16 **anap )
{
   Float32 log_en, lsf[M], lsp[M], lsp_q[M];
   Word32 i, j;


   if ( ( compute_sid_flag != 0 ) ) {
      log_en = 0;
      memset( lsp, 0, sizeof( Float32 )*M );

      for ( i = 0; i < DTX_HIST_SIZE; i++ ) {
         log_en += log_en_hist[i];

         for ( j = 0; j < M; j++ ) {
            lsp[j] += lsp_hist[i * M + j];
         }
      }
      log_en = log_en * 0.125F;

      for ( j = 0; j < M; j++ ) {
         lsp[j] = lsp[j] * 0.125F;
      }

      log_en = log_en + 2.5F;
      *log_en_index = ( Word16 )( ( log_en * 4 ) + 0.5F );   /* 6 bits */
      if ( *log_en_index > 63 ) {
         *log_en_index = 63;
      }

      if ( *log_en_index < 0 ) {
         *log_en_index = 0;
      }

      if (*log_en_index > 46){
         past_qua_en[0] = NB_QUA_CODE+VQ_SIZE_HIGHRATES+VQ_SIZE_LOWRATES+(MR475_VQ_SIZE*2) + 46;
         past_qua_en[1] = NB_QUA_CODE+VQ_SIZE_HIGHRATES+VQ_SIZE_LOWRATES+(MR475_VQ_SIZE*2) + 46;
         past_qua_en[2] = NB_QUA_CODE+VQ_SIZE_HIGHRATES+VQ_SIZE_LOWRATES+(MR475_VQ_SIZE*2) + 46;
         past_qua_en[3] = NB_QUA_CODE+VQ_SIZE_HIGHRATES+VQ_SIZE_LOWRATES+(MR475_VQ_SIZE*2) + 46;
      }
      else {
         past_qua_en[0] = NB_QUA_CODE+VQ_SIZE_HIGHRATES+VQ_SIZE_LOWRATES+(MR475_VQ_SIZE*2) + *log_en_index;
         past_qua_en[1] = NB_QUA_CODE+VQ_SIZE_HIGHRATES+VQ_SIZE_LOWRATES+(MR475_VQ_SIZE*2) + *log_en_index;
         past_qua_en[2] = NB_QUA_CODE+VQ_SIZE_HIGHRATES+VQ_SIZE_LOWRATES+(MR475_VQ_SIZE*2) + *log_en_index;
         past_qua_en[3] = NB_QUA_CODE+VQ_SIZE_HIGHRATES+VQ_SIZE_LOWRATES+(MR475_VQ_SIZE*2) + *log_en_index;
      }

      Lsp_lsf( lsp, lsf );
      Reorder_lsf( lsf, 0.00625F );
      Lsf_lsp( lsf, lsp );
      Q_plsf_3( MRDTX, past_rq, lsp, lsp_q, lsp_index, init_lsf_vq_index );
   }

   *( *anap )++ = ( Word16 )*init_lsf_vq_index;
   *( *anap )++ = lsp_index[0];
   *( *anap )++ = lsp_index[1];
   *( *anap )++ = lsp_index[2];
   *( *anap )++ = *log_en_index;

   return 0;
}

static void cod_amr( cod_amrState *st, enum ModeNB mode, Float32 new_speech[], Word16 ana[], enum ModeNB *used_mode, Float32 synth[] )
{
   Float32 A_t[( MP1 ) * 4];   /* A(z) unquantized for the 4 subframes */
   Float32 Aq_t[( MP1 ) * 4];   /* A(z)   quantized for the 4 subframes */
   Float32 *A, *Aq;   /* Pointer on Aq_t */
   Float32 lsp_new[M];
   Float32 xn[L_SUBFR];   /* Target vector for pitch search */
   Float32 xn2[L_SUBFR];   /* Target vector for codebook search */
   Float32 code[L_SUBFR];   /* Fixed codebook excitation */
   Float32 y1[L_SUBFR];   /* Filtered adaptive excitation */
   Float32 y2[L_SUBFR];   /* Filtered fixed codebook excitation */
   Float32 gCoeff[3];   /* Correlations between xn, y1, & y2: */
   Float32 res[L_SUBFR];   /* Short term (LPC) prediction residual */
   Float32 res2[L_SUBFR];   /* Long term (LTP) prediction residual */
   Float32 xn_sf0[L_SUBFR];   /* Target vector for pitch search */
   Float32 y2_sf0[L_SUBFR];   /* Filtered codebook innovation */
   Float32 code_sf0[L_SUBFR];   /* Fixed codebook excitation */
   Float32 h1_sf0[L_SUBFR];   /* The impulse response of sf0 */
   Float32 mem_syn_save[M];   /* Filter memory */
   Float32 mem_w0_save[M];   /* Filter memory */
   Float32 mem_err_save[M];   /* Filter memory */
   Float32 sharp_save = 0;   /* Sharpening */
   Float32 gain_pit_sf0;   /* Quantized pitch gain for sf0 */
   Float32 gain_code_sf0;   /* Quantized codebook gain for sf0 */
   Word16 i_subfr_sf0 = 0;   /* Position in exc[] for sf0 */
   Float32 gain_pit, gain_code;
   Float32 gp_limit;   /* pitch gain limit value */
   Word32 T0_sf0 = 0;   /* Integer pitch lag of sf0 */
   Word32 T0_frac_sf0 = 0;   /* Fractional pitch lag of sf0 */
   Word32 T0, T0_frac;
   Word32 T_op[2];
   Word32 evenSubfr;
   Word32 i;
   Word16 i_subfr, subfrNr;
   Word16 lsp_flag = 0;   /* indicates resonance in LPC filter */
   Word16 compute_sid_flag;
   Word16 vad_flag;

   memcpy( st->new_speech, new_speech, L_FRAME <<2 );

   if ( st->dtx ) {
#ifdef VAD2
     vad_flag = vad2 (st->vadSt, st->new_speech);
     vad_flag = vad2 (st->vadSt, st->new_speech+80) || vad_flag;
#else
      vad_flag = vad( st->vadSt, st->new_speech );
#endif
      if ( *used_mode < 0 )
         vad_flag = 1;
      *used_mode = mode;

      compute_sid_flag = tx_dtx_handler( vad_flag, &st->dtxEncSt->decAnaElapsedCount, &st->dtxEncSt->dtxHangoverCount, used_mode );
   }
   else {
      compute_sid_flag = 0;
      *used_mode = mode;
   }

   lpc( st->lpcSt->LevinsonSt->old_A, st->p_window, st->p_window_12k2, A_t, mode);
   lsp( mode, *used_mode, st->lspSt->lsp_old, st->lspSt->lsp_old_q, st->lspSt->
         qSt->past_rq, A_t, Aq_t, lsp_new, &ana );

   dtx_buffer( &st->dtxEncSt->hist_ptr, st->dtxEncSt->lsp_hist, lsp_new, st->
         new_speech, st->dtxEncSt->log_en_hist );

   if ( *used_mode == MRDTX ) {
      dtx_enc( &st->dtxEncSt->log_en_index, st->dtxEncSt->log_en_hist, st->
            dtxEncSt->lsp_hist, st->dtxEncSt->lsp_index, &st->dtxEncSt->
            init_lsf_vq_index, compute_sid_flag, &st->lspSt->qSt->past_rq[0], st
            ->gainQuantSt->gc_predSt->past_qua_en, &ana );
      memset( st->old_exc, 0, ( PIT_MAX + L_INTERPOL )<<2 );
      memset( st->mem_w0, 0, M <<2 );
      memset( st->mem_err, 0, M <<2 );
      memset( st->zero, 0, L_SUBFR <<2 );
      memset( st->hvec, 0, L_SUBFR <<2 );
      memset( st->lspSt->qSt->past_rq, 0, M <<2 );
      memcpy( st->lspSt->lsp_old, lsp_new, M <<2 );
      memcpy( st->lspSt->lsp_old_q, lsp_new, M <<2 );

      st->clLtpSt->pitchSt->T0_prev_subframe = 0;
      st->sharp = 0;
   }
   else {
      lsp_flag = check_lsp( &st->tonStabSt->count, st->lspSt->lsp_old );
   }

#ifdef VAD2
   if (st->dtx) {
      st->vadSt->Rmax = 0.0;
      st->vadSt->R0 = 0.0;
   }
#endif

   for ( subfrNr = 0, i_subfr = 0; subfrNr < 2; subfrNr++, i_subfr +=
         L_FRAME_BY2 ) {
      pre_big( mode, gamma1, gamma1_12k2, gamma2, A_t, i_subfr, st->speech, st->
            mem_w, st->wsp );

      if ( ( mode != MR475 ) && ( mode != MR515 ) ) {
         ol_ltp( mode, st->vadSt, &st->wsp[i_subfr], &T_op[subfrNr], st->
               ol_gain_flg, &st->pitchOLWghtSt->old_T0_med, &st->pitchOLWghtSt->
               wght_flg, &st->pitchOLWghtSt->ada_w, st->old_lags, st->dtx,
               subfrNr );
      }
   }

   if ( ( mode == MR475 ) || ( mode == MR515 ) ) {
      ol_ltp( mode, st->vadSt, &st->wsp[0], &T_op[0], st->ol_gain_flg, &st->
            pitchOLWghtSt->old_T0_med, &st->pitchOLWghtSt->wght_flg, &st->
            pitchOLWghtSt->ada_w, st->old_lags, st->dtx, 1 );
      T_op[1] = T_op[0];
   }

#ifdef VAD2
   if (st->dtx) {
      LTP_flag_update(st->vadSt, mode);
   }
#endif

#ifndef VAD2
   if ( st->dtx ) {
      vad_pitch_detection( st->vadSt, T_op );
   }
#endif

   if ( *used_mode == MRDTX ) {
      goto the_end;
   }

   A = A_t;

   Aq = Aq_t;
   evenSubfr = 0;
   subfrNr = -1;

   for ( i_subfr = 0; i_subfr < L_FRAME; i_subfr += L_SUBFR ) {
      subfrNr += 1;
      evenSubfr = 1 - evenSubfr;

      if ( ( evenSubfr != 0 ) && ( *used_mode == MR475 ) ) {
         memcpy( mem_syn_save, st->mem_syn, M <<2 );
         memcpy( mem_w0_save, st->mem_w0, M <<2 );
         memcpy( mem_err_save, st->mem_err, M <<2 );
         sharp_save = st->sharp;
      }

      if ( *used_mode != MR475 ) {
         subframePreProc( *used_mode, gamma1, gamma1_12k2, gamma2, A, Aq, &st->
               speech[i_subfr], st->mem_err, st->mem_w0, st->zero, st->ai_zero,
               &st->exc[i_subfr], st->h1, xn, res, st->error );
      }

      else {
         subframePreProc( *used_mode, gamma1, gamma1_12k2, gamma2, A, Aq, &st->
               speech[i_subfr], st->mem_err, mem_w0_save, st->zero, st->ai_zero,
               &st->exc[i_subfr], st->h1, xn, res, st->error );

         if ( evenSubfr != 0 ) {
            memcpy( h1_sf0, st->h1, L_SUBFR <<2 );
         }
      }

      memcpy( res2, res, L_SUBFR <<2 );

      cl_ltp( &st->clLtpSt->pitchSt->T0_prev_subframe, st->tonStabSt->gp, *
            used_mode, i_subfr, T_op, st->h1, &st->exc[i_subfr], res2, xn,
            lsp_flag, xn2, y1, &T0, &T0_frac, &gain_pit, gCoeff, &ana, &gp_limit
            );

      if ( ( subfrNr == 0 ) && ( st->ol_gain_flg[0] > 0 ) ) {
         st->old_lags[1] = T0;
      }

      if ( ( subfrNr == 3 ) && ( st->ol_gain_flg[1] > 0 ) ) {
         st->old_lags[0] = T0;
      }

      cbsearch( *used_mode, subfrNr, xn2, st->h1, T0, st->sharp, gain_pit, code,
            y2, res2, &ana );

      gainQuant( *used_mode, evenSubfr, st->gainQuantSt->gc_predSt->past_qua_en,
            st->gainQuantSt->gc_predUncSt->past_qua_en, st->gainQuantSt->
            sf0_coeff, &st->gainQuantSt->sf0_target_en, &st->gainQuantSt->
            sf0_gcode0_exp, &st->gainQuantSt->
            sf0_gcode0_fra, &st->gainQuantSt->gain_idx_ptr, &gain_pit_sf0, &
            gain_code_sf0, res, &st->exc[i_subfr], code, xn, xn2, y1, y2, gCoeff
            , gp_limit, &gain_pit, &gain_code, &st->gainQuantSt->adaptSt->
            prev_gc, &st->gainQuantSt->adaptSt->onset, st->gainQuantSt->adaptSt
            ->ltpg_mem, &st->gainQuantSt->adaptSt->prev_alpha, &ana );

      for ( i = 0; i < N_FRAME - 1; i++ ) {
         st->tonStabSt->gp[i] = st->tonStabSt->gp[i + 1];
      }
      st->tonStabSt->gp[N_FRAME - 1] = gain_pit;

      if ( *used_mode != MR475 ) {
         subframePostProc( st->speech, i_subfr, gain_pit, gain_code, Aq, synth,
               xn, code, y1, y2, st->mem_syn, st->mem_err, st->mem_w0, st->exc,
               &st->sharp );
      }
      else {
         if ( evenSubfr != 0 ) {
            i_subfr_sf0 = i_subfr;
            memcpy( xn_sf0, xn, L_SUBFR <<2 );
            memcpy( y2_sf0, y2, L_SUBFR <<2 );
            memcpy( code_sf0, code, L_SUBFR <<2 );
            T0_sf0 = T0;
            T0_frac_sf0 = T0_frac;
            subframePostProc( st->speech, i_subfr, gain_pit, gain_code, Aq,
                  synth, xn, code, y1, y2, mem_syn_save, st->mem_err,
                  mem_w0_save, st->exc, &st->sharp );
            st->sharp = sharp_save;
         } else {
            memcpy( st->mem_err, mem_err_save, M <<2 );

            Pred_lt_3or6( &st->exc[i_subfr_sf0], T0_sf0, T0_frac_sf0, 1 );
            Convolve( &st->exc[i_subfr_sf0], h1_sf0, y1 );
            Aq -= MP1;
            subframePostProc( st->speech, i_subfr_sf0, gain_pit_sf0,
                  gain_code_sf0, Aq, synth, xn_sf0, code_sf0, y1, y2_sf0, st->
                  mem_syn, st->mem_err, st->mem_w0, st->exc, &sharp_save );

            Aq += MP1;

            subframePreProc( *used_mode, gamma1, gamma1_12k2, gamma2, A, Aq, &st
                  ->speech[i_subfr], st->mem_err, st->mem_w0, st->zero, st->
                  ai_zero, &st->exc[i_subfr], st->h1, xn, res, st->error );

            Pred_lt_3or6( &st->exc[i_subfr], T0, T0_frac, 1 );
            Convolve( &st->exc[i_subfr], st->h1, y1 );
            subframePostProc( st->speech, i_subfr, gain_pit, gain_code, Aq,
                  synth, xn, code, y1, y2, st->mem_syn, st->mem_err, st->mem_w0,
                  st->exc, &st->sharp );
         }
      }
      A += MP1;
      Aq += MP1;
   }
the_end:

   for ( i = 0; i < PIT_MAX; i++ ) {
      st->old_wsp[i] = st->old_wsp[L_FRAME + i];
   }

   for ( i = 0; i < PIT_MAX + L_INTERPOL; i++ ) {
      st->old_exc[i] = st->old_exc[L_FRAME + i];
   }

   for ( i = 0; i < L_TOTAL - L_FRAME; i++ ) {
      st->old_speech[i] = st->old_speech[L_FRAME + i];
   }
}

static Word32 Pre_Process_reset( Pre_ProcessState *state )
{
   if ( state == ( Pre_ProcessState * )NULL ) {
      fprintf( stderr, "Pre_Process_reset: invalid parameter\n" );
      return-1;
   }
   state->y2 = 0;
   state->y1 = 0;
   state->x0 = 0;
   state->x1 = 0;
   return 0;
}

static void Pre_Process_exit( Pre_ProcessState **state )
{
   if ( state == NULL || *state == NULL )
      return;

   free( *state );
   *state = NULL;
   return;
}

static Word32 Pre_Process_init( Pre_ProcessState **state )
{
   Pre_ProcessState * s;

   if ( state == ( Pre_ProcessState * * )NULL ) {
      fprintf( stderr, "Pre_Process_init: invalid parameter\n" );
      return-1;
   }
   *state = NULL;

   if ( ( s = ( Pre_ProcessState * ) malloc( sizeof( Pre_ProcessState ) ) ) == NULL ) {
      fprintf( stderr, "Pre_Process_init: can not malloc state structure\n" );
      return-1;
   }
   Pre_Process_reset( s );
   *state = s;
   return 0;
}

static void Pre_Process( Float32 *y2, Float32 *y1, Float32 *x0, Float32*x1, Word16 *speech, Float32 *f_speech )
{
   for ( Word32 i = 0; i < 160; i++ ) {
      Float32 x2 = *x1;
      *x1 = *x0;
      *x0 = speech[i];
      Float32 tmp = ( Float32 )( 0.4636230465* *x0 - 0.92724705 * *x1 + 0.4636234515 * x2 + 1.906005859 * *y1 - 0.911376953 * *y2 );
      f_speech[i] = tmp;
      *y2 = *y1;
      *y1 = tmp;
   }

   if ( ( fabs( *y1 )+fabs( *y2 ) ) < 0.0000000001 )
      *y2 = *y1 = 0;
}

static void cod_amr_reset( cod_amrState *s, Word32 dtx )
{
   Word32 i;

   s->dtx = dtx;
   s->clLtpSt->pitchSt->T0_prev_subframe = 0;

   memset( s->lspSt->qSt->past_rq, 0, sizeof( Float32 )*M );
   memcpy( s->lspSt->lsp_old, lsp_init_data, sizeof( lsp_init_data ) );
   memcpy( s->lspSt->lsp_old_q, lsp_init_data, sizeof( lsp_init_data ) );

   for ( i = 0; i < NPRED; i++ ) {
      s->gainQuantSt->gc_predSt->past_qua_en[i] = NB_QUA_CODE+VQ_SIZE_HIGHRATES+VQ_SIZE_LOWRATES+MR475_VQ_SIZE*2+DTX_VQ_SIZE;
      s->gainQuantSt->gc_predUncSt->past_qua_en[i] = NB_QUA_CODE+VQ_SIZE_HIGHRATES+VQ_SIZE_LOWRATES+MR475_VQ_SIZE*2+DTX_VQ_SIZE;
   }

   s->gainQuantSt->adaptSt->onset = 0;
   s->gainQuantSt->adaptSt->prev_alpha = 0.0F;
   s->gainQuantSt->adaptSt->prev_gc = 0.0F;
   memset( s->gainQuantSt->adaptSt->ltpg_mem, 0, sizeof( Float32 )*LTPG_MEM_SIZE );
   s->gainQuantSt->sf0_gcode0_exp = 0;
   s->gainQuantSt->sf0_gcode0_fra = 0;
   s->gainQuantSt->sf0_target_en = 0.0F;
   memset( s->gainQuantSt->sf0_coeff, 0, sizeof( Float32 )*5 );
   s->gainQuantSt->gain_idx_ptr = NULL;
   s->pitchOLWghtSt->old_T0_med = 40;
   s->pitchOLWghtSt->ada_w = 0.0F;
   s->pitchOLWghtSt->wght_flg = 0;
   s->tonStabSt->count = 0;
   memset( s->tonStabSt->gp, 0, sizeof( Float32 )*N_FRAME );
   s->lpcSt->LevinsonSt->old_A[0] = 1.0F;
   memset( &s->lpcSt->LevinsonSt->old_A[1], 0, sizeof( Float32 )*M );

#ifdef VAD2
   s->vadSt->pre_emp_mem = 0.0;
   s->vadSt->update_cnt = 0;
   s->vadSt->hyster_cnt = 0;
   s->vadSt->last_update_cnt = 0;
   for ( i = 0; i < NUM_CHAN; i++ ) {
     s->vadSt->ch_enrg_long_db[i] = 0.0;
     s->vadSt->ch_enrg[i] = 0.0;
     s->vadSt->ch_noise[i] = 0.0;
   }
   s->vadSt->Lframe_cnt = 0L;
   s->vadSt->tsnr = 0.0;
   s->vadSt->hangover = 0;
   s->vadSt->burstcount = 0;
   s->vadSt->fupdate_flag = 0;
   s->vadSt->negSNRvar = 0.0;
   s->vadSt->negSNRbias = 0.0;
   s->vadSt->R0 = 0.0;
   s->vadSt->Rmax = 0.0;
   s->vadSt->LTP_flag = 0;
#else
   s->vadSt->oldlag_count = 0;
   s->vadSt->oldlag = 0;
   s->vadSt->pitch = 0;
   s->vadSt->tone = 0;
   s->vadSt->complex_high = 0;
   s->vadSt->complex_low = 0;
   s->vadSt->complex_hang_timer = 0;
   s->vadSt->vadreg = 0;
   s->vadSt->burst_count = 0;
   s->vadSt->hang_count = 0;
   s->vadSt->complex_hang_count = 0;

   for ( i = 0; i < 3; i++ ) {
      s->vadSt->a_data5[i][0] = 0;
      s->vadSt->a_data5[i][1] = 0;
   }

   for ( i = 0; i < 5; i++ ) {
      s->vadSt->a_data3[i] = 0;
   }

   for ( i = 0; i < COMPLEN; i++ ) {
      s->vadSt->bckr_est[i] = NOISE_INIT;
      s->vadSt->old_level[i] = NOISE_INIT;
      s->vadSt->ave_level[i] = NOISE_INIT;
      s->vadSt->sub_level[i] = 0;
   }
   s->vadSt->best_corr_hp = CVAD_LOWPOW_RESET;
   s->vadSt->speech_vad_decision = 0;
   s->vadSt->complex_warning = 0;
   s->vadSt->sp_burst_count = 0;
   s->vadSt->corr_hp_fast = CVAD_LOWPOW_RESET;
#endif

   s->dtxEncSt->hist_ptr = 0;
   s->dtxEncSt->log_en_index = 0;
   s->dtxEncSt->init_lsf_vq_index = 0;
   s->dtxEncSt->lsp_index[0] = 0;
   s->dtxEncSt->lsp_index[1] = 0;
   s->dtxEncSt->lsp_index[2] = 0;

   for ( i = 0; i < DTX_HIST_SIZE; i++ ) {
      memcpy( &s->dtxEncSt->lsp_hist[i * M], lsp_init_data, sizeof( Float32 )*M );
   }
   memset( s->dtxEncSt->log_en_hist, 0, M * sizeof( Float32 ) );
   s->dtxEncSt->dtxHangoverCount = DTX_HANG_CONST;
   s->dtxEncSt->decAnaElapsedCount = DTX_ELAPSED_FRAMES_THRESH;
   s->new_speech = s->old_speech + L_TOTAL - L_FRAME;
   s->speech = s->new_speech - L_NEXT;
   s->p_window = s->old_speech + L_TOTAL - L_WINDOW;
   s->p_window_12k2 = s->p_window - L_NEXT;
   s->wsp = s->old_wsp + PIT_MAX;
   s->exc = s->old_exc + PIT_MAX + L_INTERPOL;
   s->zero = s->ai_zero + MP1;
   s->error = s->mem_err + M;
   s->h1 = &s->hvec[L_SUBFR];

   memset( s->old_speech, 0, sizeof( Float32 )*L_TOTAL );
   memset( s->old_exc, 0, sizeof( Float32 )*( PIT_MAX + L_INTERPOL ) );
   memset( s->old_wsp, 0, sizeof( Float32 )*PIT_MAX );
   memset( s->mem_syn, 0, sizeof( Float32 )*M );
   memset( s->mem_w, 0, sizeof( Float32 )*M );
   memset( s->mem_w0, 0, sizeof( Float32 )*M );
   memset( s->mem_err, 0, sizeof( Float32 )*M );
   memset( s->ai_zero, 0, sizeof( Float32 )*L_SUBFR );
   memset( s->hvec, 0, sizeof( Float32 )*L_SUBFR );

   for ( i = 0; i < 5; i++ ) {
      s->old_lags[i] = 40;
   }
   s->sharp = 0.0F;
}

static Word32 cod_amr_init( cod_amrState **state, Word32 dtx )
{
   cod_amrState * s;

   if ( ( s = ( cod_amrState * ) malloc( sizeof( cod_amrState ) ) ) == NULL ) {
      return-1;
   }

   if ( ( s->clLtpSt = ( clLtpState * ) malloc( sizeof( clLtpState ) ) ) == NULL
         ) {
      return-1;
   }

   if ( ( s->clLtpSt->pitchSt = ( Pitch_frState * ) malloc( sizeof(
         Pitch_frState ) ) ) == NULL ) {
      return-1;
   }

   if ( ( s->lspSt = ( lspState * ) malloc( sizeof( lspState ) ) ) == NULL ) {
      return-1;
   }

   if ( ( s->lspSt->qSt = ( Q_plsfState * ) malloc( sizeof( Q_plsfState ) ) ) ==
         NULL ) {
      return-1;
   }

   if ( ( s->gainQuantSt = ( gainQuantState * ) malloc( sizeof( gainQuantState ) ) ) == NULL ) {
      return-1;
   }

   if ( ( s->gainQuantSt->gc_predSt = ( gc_predState * ) malloc( sizeof( gc_predState ) ) ) == NULL ) {
      return-1;
   }

   if ( ( s->gainQuantSt->gc_predUncSt = ( gc_predState * ) malloc( sizeof( gc_predState ) ) ) == NULL ) {
      return-1;
   }

   if ( ( s->gainQuantSt->adaptSt = ( gain_adaptState * ) malloc( sizeof( gain_adaptState ) ) ) == NULL ) {
      return-1;
   }

   if ( ( s->pitchOLWghtSt = ( pitchOLWghtState * ) malloc( sizeof( pitchOLWghtState ) ) ) == NULL ) {
      return-1;
   }

   if ( ( s->tonStabSt = ( tonStabState * ) malloc( sizeof( tonStabState ) ) ) == NULL ) {
      return-1;
   }

   if ( ( s->lpcSt = ( lpcState * ) malloc( sizeof( lpcState ) ) ) == NULL ) {
      return-1;
   }

   if ( ( s->lpcSt->LevinsonSt = ( LevinsonState * ) malloc( sizeof( LevinsonState ) ) ) == NULL ) {
      return-1;
   }

   if ( ( s->vadSt = ( vadState * ) malloc( sizeof( vadState ) ) ) == NULL ) {
      return-1;
   }

   if ( ( s->dtxEncSt = ( dtx_encState * ) malloc( sizeof( dtx_encState ) ) ) == NULL ) {
      return-1;
   }

   cod_amr_reset( s, dtx );
   *state = s;
   return 0;
}

static void cod_amr_exit( cod_amrState **state )
{
   if ( state == NULL || *state == NULL )
      return;

   free( ( *state )->vadSt );
   free( ( *state )->gainQuantSt->gc_predSt );
   free( ( *state )->gainQuantSt->gc_predUncSt );
   free( ( *state )->gainQuantSt->adaptSt );
   free( ( *state )->clLtpSt->pitchSt );
   free( ( *state )->lspSt->qSt );
   free( ( *state )->lpcSt->LevinsonSt );
   free( ( *state )->lpcSt );
   free( ( *state )->lspSt );
   free( ( *state )->clLtpSt );
   free( ( *state )->gainQuantSt );
   free( ( *state )->pitchOLWghtSt );
   free( ( *state )->tonStabSt );
   free( ( *state )->dtxEncSt );
   free( *state );
   *state = NULL;
}

void * Speech_Encode_Frame_init( int dtx )
{
   Speech_Encode_FrameState * s;

   if ( ( s = ( Speech_Encode_FrameState * ) malloc( sizeof( Speech_Encode_FrameState ) ) ) == NULL ) {
      fprintf( stderr, "Speech_Encode_Frame_init: can not malloc state structure\n" );
      return NULL;
   }
   s->pre_state = NULL;
   s->cod_amr_state = NULL;
   s->dtx = dtx;

   if ( Pre_Process_init( &s->pre_state ) || cod_amr_init( &s->cod_amr_state,
         dtx ) ) {
      Speech_Encode_Frame_exit( ( void ** )( &s ) );
      return NULL;
   }
   return s;
}

int Speech_Encode_Frame_reset( void *st, int dtx )
{
   Speech_Encode_FrameState * state;
   state = ( Speech_Encode_FrameState * )st;

   if ( ( Speech_Encode_FrameState * )state == NULL ) {
      fprintf( stderr, "Speech_Encode_Frame_reset: invalid parameter\n" );
      return-1;
   }
   Pre_Process_reset( state->pre_state );
   cod_amr_reset( state->cod_amr_state, dtx );
   return 0;
}

void Speech_Encode_Frame_exit( void **st )
{
   if ( ( Speech_Encode_FrameState * )( *st ) == NULL )
      return;
   Pre_Process_exit( &( ( ( Speech_Encode_FrameState * )( *st ) )->pre_state ) );
   cod_amr_exit( &( ( ( Speech_Encode_FrameState * )( *st ) )->cod_amr_state ) );

   free( *st );
   *st = NULL;
}

void Speech_Encode_Frame( void *st, enum ModeNB mode, Word16 *new_speech, Word16 *
      prm, enum ModeNB *used_mode )
{
   Float32 syn[L_FRAME];   /* Buffer for synthesis speech */
   Float32 speech[160];

   Speech_Encode_FrameState * state;
   state = ( Speech_Encode_FrameState * )st;

   for ( Word32 i = 0; i < 160; i++ ) {
      new_speech[i] = ( Word16 )( new_speech[i] & 0xfff8 );
   }

   Pre_Process( &state->pre_state->y2, &state->pre_state->y1, &state->pre_state->x0, &state->pre_state->x1, new_speech, speech );
   cod_amr( state->cod_amr_state, mode, speech, prm, used_mode, syn );
}
