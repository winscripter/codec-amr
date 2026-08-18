#ifndef VAD2
static void vad_tone_detection( vadState *st, Float32 T0, Float32 t1 )
{
   if ( ( t1 > 0 ) && ( T0 > t1 * TONE_THR ) ) {
      st->tone = st->tone | 0x00004000;
   }
}
#endif

#ifdef VAD2
static Word16 Lag_max( Float32 corr[], Float32 sig[], Word16 L_frame,
		       Word32 lag_max, Word32 lag_min, Float32 *cor_max,
		       Word32 dtx, Float32 *rmax, Float32 *r0 )
#else
static Word16 Lag_max( vadState *vadSt, Float32 corr[], Float32 sig[], Word16
      L_frame, Word32 lag_max, Word32 lag_min, Float32 *cor_max, Word32 dtx )
#endif
{
   Float32 max, T0;
   Float32 *p;
   Word32 i, j, p_max;

   max = -FLT_MAX;
   p_max = lag_max;

   for ( i = lag_max, j = ( PIT_MAX - lag_max - 1 ); i >= lag_min; i--, j-- ) {
      if ( corr[ - i] >= max ) {
         max = corr[ - i];
         p_max = i;
      }
   }

   T0 = 0.0F;
   p = &sig[ - p_max];

   for ( i = 0; i < L_frame; i++, p++ ) {
      T0 += *p * *p;
   }

   if ( dtx ) {
#ifdef VAD2
     *rmax = max;
     *r0 = T0;
#else
     vad_tone_detection( vadSt, max, T0 );
#endif
   }

   if ( T0 > 0.0F )
      T0 = 1.0F / ( Float32 )sqrt( T0 );
   else
      T0 = 0.0F;

   max *= T0;
   *cor_max = max;
   return( ( Word16 )p_max );
}

#ifndef VAD2
static void hp_max( Float32 corr[], Float32 sig[], Word32 L_frame, Word32
      lag_max, Word32 lag_min, Float32 *cor_hp_max )
{
   Float32 T0, t1, max;
   Float32 *p, *p1;
   Word32 i;

   max = -FLT_MAX;
   T0 = 0;

   for ( i = lag_max - 1; i > lag_min; i-- ) {
      T0 = ( ( corr[ - i] * 2 ) - corr[ - i-1] )-corr[ - i + 1];
      T0 = ( Float32 )fabs( T0 );

      if ( T0 >= max ) {
         max = T0;
      }
   }

   p = sig;
   p1 = &sig[0];
   T0 = 0;

   for ( i = 0; i < L_frame; i++, p++, p1++ ) {
      T0 += *p * *p1;
   }
   p = sig;
   p1 = &sig[ - 1];
   t1 = 0;

   for ( i = 0; i < L_frame; i++, p++, p1++ ) {
      t1 += *p * *p1;
   }

   T0 = T0 - t1;
   T0 = ( Float32 )fabs( T0 );

   if ( T0 != 0 ) {
      *cor_hp_max = max / T0;
   }
   else {
      *cor_hp_max = 0;
   }
}
#endif

#ifndef VAD2
static void vad_tone_detection_update( vadState *st, Word16 one_lag_per_frame )
{
   st->tone = st->tone >> 1;

   if ( one_lag_per_frame != 0 ) {
      st->tone = st->tone >> 1;
      st->tone = st->tone | 0x00002000;
   }
}
#endif

static Word32 Pitch_ol( enum ModeNB mode, vadState *vadSt, Float32 signal[],
      Word32 pit_min, Word32 pit_max, Word16 L_frame, Word32 dtx, Word16 idx )
{
   Float32 corr[PIT_MAX + 1];
   Float32 max1, max2, max3, p_max1, p_max2, p_max3;
   Float32 *corr_ptr;
   Word32 i, j;
#ifdef VAD2
   Float32 r01, r02, r03;
   Float32 rmax1, rmax2, rmax3;
#else
   Float32 corr_hp_max;
#endif

#ifndef VAD2
   if ( dtx ) {
      if ( ( mode == MR475 ) || ( mode == MR515 ) ) {
         vad_tone_detection_update( vadSt, 1 );
      }
      else {
         vad_tone_detection_update( vadSt, 0 );
      }
   }
#endif

   corr_ptr = &corr[pit_max];
   comp_corr( signal, L_frame, pit_max, pit_min, corr_ptr );

#ifdef VAD2
   j = pit_min << 2;
   p_max1 = Lag_max( corr_ptr, signal, L_frame, pit_max, j, &max1, dtx, &rmax1, &r01 );
   i = j - 1;
   j = pit_min << 1;
   p_max2 = Lag_max( corr_ptr, signal, L_frame, i, j, &max2, dtx, &rmax2, &r02 );
   i = j - 1;
   p_max3 = Lag_max( corr_ptr, signal, L_frame, i, pit_min, &max3, dtx, &rmax3, &r03 );
#else
   j = pit_min << 2;
   p_max1 = Lag_max( vadSt, corr_ptr, signal, L_frame, pit_max, j, &max1, dtx );
   i = j - 1;
   j = pit_min << 1;
   p_max2 = Lag_max( vadSt, corr_ptr, signal, L_frame, i, j, &max2, dtx );
   i = j - 1;
   p_max3 = Lag_max( vadSt, corr_ptr, signal, L_frame, i, pit_min, &max3, dtx );

   if ( dtx ) {
      if ( idx == 1 ) {
         hp_max( corr_ptr, signal, L_frame, pit_max, pit_min, &corr_hp_max );
         vadSt->best_corr_hp = corr_hp_max * 0.5F;
      }
   }
#endif

   if ( ( max1 * 0.85F ) < max2 ) {
      max1 = max2;
      p_max1 = p_max2;
#ifdef VAD2
      if (dtx) {
	rmax1 = rmax2;
	r01 = r02;
      }
#endif
   }

   if ( ( max1 * 0.85F ) < max3 ) {
      p_max1 = p_max3;
#ifdef VAD2
      if (dtx) {
	rmax1 = rmax3;
	r01 = r03;
      }
#endif
   }
#ifdef VAD2
   if (dtx) {
     vadSt->Rmax += rmax1;   /* Save max correlation */
     vadSt->R0   += r01;     /* Save max energy */
   }
#endif
   return( Word32 )p_max1;
}

static Word32 Lag_max_wght( vadState *vadSt, Float32 corr[], Float32 signal[],
      Word32 old_lag, Word32 *cor_max, Word32 wght_flg, Float32 *gain_flg,
      Word32 dtx )
{
   Float32 t0, t1, max;
   Float32 *psignal, *p1signal;
   const Float32 *ww, *we;
   Word32 i, j, p_max;


   ww = &corrweight[250];
   we = &corrweight[266 - old_lag];
   max = -FLT_MAX;
   p_max = PIT_MAX;

   if ( wght_flg > 0 ) {
      for ( i = PIT_MAX; i >= PIT_MIN; i-- ) {
         t0 = corr[ - i] * *ww--;
         t0 *= *we--;

         if ( t0 >= max ) {
            max = t0;
            p_max = i;
         }
      }
   }
   else {
      for ( i = PIT_MAX; i >= PIT_MIN; i-- ) {
         t0 = corr[ - i] * *ww--;

         if ( t0 >= max ) {
            max = t0;
            p_max = i;
         }
      }
   }
   psignal = &signal[0];
   p1signal = &signal[ - p_max];
   t0 = 0;
   t1 = 0;

   for ( j = 0; j < L_FRAME_BY2; j++, psignal++, p1signal++ ) {
      t0 += *psignal * *p1signal;
      t1 += *p1signal * *p1signal;
   }

   if ( dtx ) {
#ifdef VAD2
       vadSt->Rmax += t0;   /* Save max correlation */
       vadSt->R0   += t1;   /* Save max energy */
#else
      vad_tone_detection_update( vadSt, 0 );
      vad_tone_detection( vadSt, t0, t1 );
#endif
   }

   *gain_flg = t0 - ( t1 * 0.4F );
   *cor_max = 0;
   return( p_max );
}

static Word32 Pitch_ol_wgh( Word32 *old_T0_med, Word16 *wght_flg, Float32 *ada_w,
      vadState *vadSt, Float32 signal[], Word32 old_lags[], Float32 ol_gain_flg[],
      Word16 idx, Word32 dtx )
{
   Float32 corr[PIT_MAX + 1];
#ifndef VAD2
   Float32 corr_hp_max;
#endif
   Float32 *corrPtr;
   Word32 i, max1, p_max1;

   corrPtr = &corr[PIT_MAX];
   comp_corr( signal, L_FRAME_BY2, PIT_MAX, PIT_MIN, corrPtr );
   p_max1 = Lag_max_wght( vadSt, corrPtr, signal, *old_T0_med,
         &max1, *wght_flg, &ol_gain_flg[idx], dtx );

   if ( ol_gain_flg[idx] > 0 ) {
      for ( i = 4; i > 0; i-- ) {
         old_lags[i] = old_lags[i - 1];
      }
      old_lags[0] = p_max1;
      *old_T0_med = gmed_n( old_lags, 5 );
      *ada_w = 1;
   }
   else {
      *old_T0_med = p_max1;
      *ada_w = *ada_w * 0.9F;
   }

   if ( *ada_w < 0.3 ) {
      *wght_flg = 0;
   }
   else {
      *wght_flg = 1;
   }

#ifndef VAD2
   if ( dtx ) {
      if ( idx == 1 ) {
         hp_max( corrPtr, signal, L_FRAME_BY2, PIT_MAX, PIT_MIN, &corr_hp_max );
         vadSt->best_corr_hp = corr_hp_max * 0.5F;
      }
   }
#endif
   return( p_max1 );
}

static void ol_ltp( enum ModeNB mode, vadState *vadSt, Float32 wsp[], Word32 *T_op
      , Float32 ol_gain_flg[], Word32 *old_T0_med, Word16 *wght_flg, Float32 *ada_w
      , Word32 *old_lags, Word32 dtx, Word16 idx )
{
   if ( mode != MR102 ) {
      ol_gain_flg[0] = 0;
      ol_gain_flg[1] = 0;
   }

   if ( ( mode == MR475 ) || ( mode == MR515 ) ) {
      *T_op = Pitch_ol( mode, vadSt, wsp, PIT_MIN, PIT_MAX, L_FRAME, dtx, idx );
   }
   else {
      if ( mode <= MR795 ) {
         *T_op = Pitch_ol( mode, vadSt, wsp, PIT_MIN, PIT_MAX, L_FRAME_BY2, dtx, idx );
      }
      else if ( mode == MR102 ) {
         *T_op = Pitch_ol_wgh( old_T0_med, wght_flg, ada_w, vadSt, wsp, old_lags, ol_gain_flg, idx, dtx );
      }
      else {
         *T_op = Pitch_ol( mode, vadSt, wsp, PIT_MIN_MR122, PIT_MAX, L_FRAME_BY2, dtx, idx );
      }
   }
}

#ifndef VAD2
static void complex_estimate_adapt( vadState *st, Word16 low_power )
{
   Float32 alpha;

   if ( st->best_corr_hp < st->corr_hp_fast ) {
      alpha = st->corr_hp_fast < CVAD_THRESH_ADAPT_HIGH ? CVAD_ADAPT_FAST : CVAD_ADAPT_REALLY_FAST;
   } else {
      alpha = st->corr_hp_fast < CVAD_THRESH_ADAPT_HIGH ? CVAD_ADAPT_FAST : CVAD_ADAPT_SLOW;
   }

   st->corr_hp_fast = st->corr_hp_fast - alpha * st->corr_hp_fast + alpha * st-> best_corr_hp;

   if ( st->corr_hp_fast < CVAD_MIN_CORR ) {
      st->corr_hp_fast = CVAD_MIN_CORR;
   }

   if ( low_power != 0 ) {
      st->corr_hp_fast = CVAD_MIN_CORR;
   }
}
#endif

#ifndef VAD2
static Word32 complex_vad( vadState *st, Word16 low_power )
{
   st->complex_high = st->complex_high >> 1;
   st->complex_low = st->complex_low >> 1;

   if ( low_power == 0 ) {
      if ( st->corr_hp_fast > CVAD_THRESH_ADAPT_HIGH ) {
         st->complex_high = st->complex_high | 0x00004000;
      }

      if ( st->corr_hp_fast > CVAD_THRESH_ADAPT_LOW ) {
         st->complex_low = st->complex_low | 0x00004000;
      }
   }

   if ( st->corr_hp_fast > CVAD_THRESH_HANG ) {
      st->complex_hang_timer += 1;
   } else {
      st->complex_hang_timer = 0;
   }
   return( Word16 )( ( ( st->complex_high & 0x00007f80 ) == 0x00007f80 ) || ( (
         st->complex_low & 0x00007fff ) == 0x00007fff ) );
}
#endif

#ifndef VAD2
static void update_cntrl( vadState *st, Float32 level[] )
{
   Float32 stat_rat, num, denom;
   Float32 alpha;
   Word32 i;

   if ( st->complex_warning != 0 ) {
      if ( st->stat_count < CAD_MIN_STAT_COUNT ) {
         st->stat_count = CAD_MIN_STAT_COUNT;
      }
   }

   if ( ( ( st->pitch & 0x6000 ) == 0x6000 ) || ( ( st->tone & 0x00007c00 ) == 0x7c00 ) ) {
      st->stat_count = STAT_COUNT;
   } else {
      if ( ( st->vadreg & 0x7f80 ) == 0 ) {
         st->stat_count = STAT_COUNT;
      } else {
         stat_rat = 0;

         for ( i = 0; i < COMPLEN; i++ ) {
            if ( level[i] > st->ave_level[i] ) {
               num = level[i];
               denom = st->ave_level[i];
            } else {
               num = st->ave_level[i];
               denom = level[i];
            }

            if ( num < STAT_THR_LEVEL ) {
               num = STAT_THR_LEVEL;
            }

            if ( denom < STAT_THR_LEVEL ) {
               denom = STAT_THR_LEVEL;
            }
            stat_rat += num / denom * 64;
         }

         if ( stat_rat > STAT_THR ) {
            st->stat_count = STAT_COUNT;
         } else {
            if ( ( st->vadreg & 0x4000 ) != 0 ) {
               if ( st->stat_count != 0 ) {
                  st->stat_count -= 1;
               }
            }
         }
      }
   }

   alpha = ALPHA4;

   if ( st->stat_count == STAT_COUNT ) {
      alpha = 1.0F;
   } else if ( ( st->vadreg & 0x4000 ) == 0 ) {
      alpha = ALPHA5;
   }

   for ( i = 0; i < COMPLEN; i++ ) {
      st->ave_level[i] += alpha * ( level[i] - st->ave_level[i] );
   }
}
#endif

#ifndef VAD2
static void noise_estimate_update( vadState *st, Float32 level[] )
{
   Float32 alpha_up, alpha_down, bckr_add;
   Word32 i;

   update_cntrl( st, level );
   bckr_add = 2;

   if ( ( ( 0x7800 & st->vadreg ) == 0 ) && ( ( st->pitch & 0x7800 ) == 0 ) && ( st->complex_hang_count == 0 ) ) {
      alpha_up = ALPHA_UP1;
      alpha_down = ALPHA_DOWN1;
   }
   else {
      if ( ( st->stat_count == 0 ) && ( st->complex_hang_count == 0 ) ) {
         alpha_up = ALPHA_UP2;
         alpha_down = ALPHA_DOWN2;
      } else {
         alpha_up = 0;
         alpha_down = ALPHA3;
         bckr_add = 0;
      }
   }

   for ( i = 0; i < COMPLEN; i++ ) {
      Float32 temp;
      temp = st->old_level[i] - st->bckr_est[i];

      if ( temp < 0 ) {
         st->bckr_est[i] = ( -2 + ( st->bckr_est[i] + ( alpha_down * temp ) ) );
         if ( st->bckr_est[i] < NOISE_MIN ) {
            st->bckr_est[i] = NOISE_MIN;
         }
      }
      else {
         st->bckr_est[i] = ( bckr_add + ( st->bckr_est[i] + ( alpha_up * temp ) ) );

         if ( st->bckr_est[i] > NOISE_MAX ) {
            st->bckr_est[i] = NOISE_MAX;
         }
      }
   }

   for ( i = 0; i < COMPLEN; i++ ) {
      st->old_level[i] = level[i];
   }
}
#endif

#ifndef VAD2
static Word16 hangover_addition( vadState *st, Float32 noise_level, Word16
      low_power )
{
   Word16 hang_len, burst_len;

   if ( noise_level > HANG_NOISE_THR ) {
      burst_len = BURST_LEN_HIGH_NOISE;
      hang_len = HANG_LEN_HIGH_NOISE;
   } else {
      burst_len = BURST_LEN_LOW_NOISE;
      hang_len = HANG_LEN_LOW_NOISE;
   }

   if ( low_power != 0 ) {
      st->burst_count = 0;
      st->hang_count = 0;
      st->complex_hang_count = 0;
      st->complex_hang_timer = 0;
      return 0;
   }

   if ( st->complex_hang_timer > CVAD_HANG_LIMIT ) {
      if ( st->complex_hang_count < CVAD_HANG_LENGTH ) {
         st->complex_hang_count = CVAD_HANG_LENGTH;
      }
   }

   if ( st->complex_hang_count != 0 ) {
      st->burst_count = BURST_LEN_HIGH_NOISE;
      st->complex_hang_count -= 1;
      return 1;
   }
   else {
      if ( ( ( st->vadreg & 0x3ff0 ) == 0 ) && ( st->corr_hp_fast > CVAD_THRESH_IN_NOISE ) ) {
         return 1;
      }
   }

   if ( ( st->vadreg & 0x4000 ) != 0 ) {
      st->burst_count += 1;

      if ( st->burst_count >= burst_len ) {
         st->hang_count = hang_len;
      }
      return 1;
   } else {
      st->burst_count = 0;

      if ( st->hang_count > 0 ) {
         st->hang_count -= 1;
         return 1;
      }
   }
   return 0;
}
#endif

#ifndef VAD2
static Word16 vad_decision( vadState *st, Float32 level[COMPLEN], Float32
      pow_sum )
{
   Float32 snr_sum, temp, vad_thr, noise_level;
   Word32 i;
   Word16 low_power_flag;

   snr_sum = 0;

   for ( i = 0; i < COMPLEN; i++ ) {
      temp = level[i] / st->bckr_est[i];
      snr_sum += temp * temp;
   }
   snr_sum = snr_sum * 56.8889F;

   noise_level = st->bckr_est[0] + st->bckr_est[1] + st->bckr_est[2] + st->bckr_est[3] + st->bckr_est[4] + st->bckr_est[5] + st->bckr_est[6] + st->bckr_est[7] + st->bckr_est[8];
   noise_level = noise_level * 0.111111F;

   vad_thr = VAD_SLOPE * ( noise_level - VAD_P1 ) + VAD_THR_HIGH;

   if ( vad_thr < VAD_THR_LOW ) {
      vad_thr = VAD_THR_LOW;
   }

   st->vadreg >>= 1;

   if ( snr_sum > vad_thr ) {
      st->vadreg = st->vadreg | 0x4000;
   }

   if ( pow_sum < VAD_POW_LOW ) {
      low_power_flag = 1;
   } else {
      low_power_flag = 0;
   }

   complex_estimate_adapt( st, low_power_flag );
   st->complex_warning = complex_vad( st, low_power_flag );
   noise_estimate_update( st, level );

   st->speech_vad_decision = hangover_addition( st, noise_level, low_power_flag);
   return( st->speech_vad_decision );
}
#endif

#ifndef VAD2
static Float32 level_calculation( Float32 data[], Float32 *sub_level, Word16 count1, Word16 count2, Word16 ind_m, Word16 ind_a, Word16 scale )
{
   Float32 level, temp1;
   Word32 i;

   temp1 = 0;

   for ( i = count1; i < count2; i++ ) {
      temp1 += ( Float32 )fabs( data[ind_m * i + ind_a] );
   }
   level = temp1 + *sub_level;
   *sub_level = temp1;

   for ( i = 0; i < count1; i++ ) {
      level += ( Float32 )fabs( data[ind_m * i + ind_a] );
   }
   return( scale * level );
}
#endif

#ifndef VAD2
static void filter3( Float32 *in0, Float32 *in1, Float32 *data )
{
   Float32 temp1 = *in1 - ( COEFF3 * *data );
   Float32 temp2 = *data + ( COEFF3 * temp1 );
   *data = temp1;
   *in1 = ( *in0 - temp2 ) * 0.5F;
   *in0 = ( *in0 + temp2 ) * 0.5F;
}
#endif

#ifndef VAD2
static void filter5( Float32 *in0, Float32 *in1, Float32 data[] )
{
   Float32 temp0 = *in0 - ( COEFF5_1 * data[0] );
   Float32 temp1 = data[0] + ( COEFF5_1 * temp0 );
   data[0] = temp0;
   temp0 = *in1 - ( COEFF5_2 * data[1] );
   Float32 temp2 = data[1] + ( COEFF5_2 * temp0 );
   data[1] = temp0;
   *in0 = ( temp1 + temp2 ) * 0.5F;
   *in1 = ( temp1 - temp2 ) * 0.5F;
}
#endif

#ifndef VAD2
static void first_filter_stage( Float32 in[], Float32 out[], Float32 data[] )
{
   Float32 temp0, temp1, temp2, temp3;
   Float32 data0, data1;
   Word32 i;

   data0 = data[0];
   data1 = data[1];

   for ( i = 0; i < L_SUBFR; i++ ) {
      temp0 = ( in[4 * i + 0] * 0.25F ) - ( COEFF5_1 * data0 );
      temp1 = data0 + ( COEFF5_1 * temp0 );
      temp3 = ( in[4*i+1]*0.25F )-( COEFF5_2 * data1 );
      temp2 = data1 + ( COEFF5_2 * temp3 );
      out[4 * i + 0] = temp1 + temp2;
      out[4 * i + 1] = temp1 - temp2;
      data0 = ( in[4 * i + 2] * 0.25F ) - ( COEFF5_1 * temp0 );
      temp1 = temp0 + ( COEFF5_1 * data0 );
      data1 = ( in[4 * i + 3] * 0.25F ) - ( COEFF5_2 * temp3 );
      temp2 = temp3 + ( COEFF5_2 * data1 );
      out[4 * i + 2] = temp1 + temp2;
      out[4 * i + 3] = temp1 - temp2;
   }
   data[0] = data0;
   data[1] = data1;
}
#endif

#ifndef VAD2
static void filter_bank( vadState *st, Float32 in[], Float32 level[] )
{
   Word32 i;
   Float32 tmp_buf[FRAME_LEN];

   first_filter_stage( in, tmp_buf, st->a_data5[0] );

   for ( i = 0; i < FRAME_LEN / 4; i++ ) {
      filter5( &tmp_buf[4 * i], &tmp_buf[4 * i + 2], st->a_data5[1] );
      filter5( &tmp_buf[4 * i +1], &tmp_buf[4 * i + 3], st->a_data5[2] );
   }

   for ( i = 0; i < FRAME_LEN / 8; i++ ) {
      filter3( &tmp_buf[8 * i + 0], &tmp_buf[8 * i + 4], &st->a_data3[0] );
      filter3( &tmp_buf[8 * i + 2], &tmp_buf[8 * i + 6], &st->a_data3[1] );
      filter3( &tmp_buf[8 * i + 3], &tmp_buf[8 * i + 7], &st->a_data3[4] );
   }

   for ( i = 0; i < FRAME_LEN / 16; i++ ) {
      filter3( &tmp_buf[16 * i + 0], &tmp_buf[16 * i + 8], &st->a_data3[2] );
      filter3( &tmp_buf[16 * i + 4], &tmp_buf[16 * i + 12], &st->a_data3[3] );
   }

   level[8] = level_calculation( tmp_buf, &st->sub_level[8], FRAME_LEN /4 - 8, FRAME_LEN /4, 4, 1, 1 );
   level[7] = level_calculation( tmp_buf, &st->sub_level[7], FRAME_LEN /8 - 4, FRAME_LEN /8, 8, 7, 2 );
   level[6] = level_calculation( tmp_buf, &st->sub_level[6], FRAME_LEN /8 - 4, FRAME_LEN /8, 8, 3, 2 );
   level[5] = level_calculation( tmp_buf, &st->sub_level[5], FRAME_LEN /8 - 4, FRAME_LEN /8, 8, 2, 2 );
   level[4] = level_calculation( tmp_buf, &st->sub_level[4], FRAME_LEN /8 - 4, FRAME_LEN /8, 8, 6, 2 );
   level[3] = level_calculation( tmp_buf, &st->sub_level[3], FRAME_LEN /16 - 2, FRAME_LEN /16, 16, 4, 2 );
   level[2] = level_calculation( tmp_buf, &st->sub_level[2], FRAME_LEN /16 - 2, FRAME_LEN /16, 16, 12, 2 );
   level[1] = level_calculation( tmp_buf, &st->sub_level[1], FRAME_LEN /16 - 2, FRAME_LEN /16, 16, 8, 2 );
   level[0] = level_calculation( tmp_buf, &st->sub_level[0], FRAME_LEN /16 - 2, FRAME_LEN /16, 16, 0, 2 );
}
#endif
