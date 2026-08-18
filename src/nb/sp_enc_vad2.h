#ifndef VAD2
static Word16 vad( vadState *st, Float32 in_buf[] )
{
   Float32 level[COMPLEN];
   Float32 pow_sum;
   Word32 i;

   pow_sum = 0L;

   for ( i = -40; i < 120; i += 8 ) {
      pow_sum += in_buf[i] * in_buf[i];
      pow_sum += in_buf[i + 1] *in_buf[i + 1];
      pow_sum += in_buf[i + 2] * in_buf[i + 2];
      pow_sum += in_buf[i + 3] * in_buf[i + 3];
      pow_sum += in_buf[i + 4] * in_buf[i + 4];
      pow_sum += in_buf[i + 5] * in_buf[i + 5];
      pow_sum += in_buf[i + 6] * in_buf[i + 6];
      pow_sum += in_buf[i + 7] * in_buf[i + 7];
   }

   if ( pow_sum < POW_PITCH_THR ) {
      st->pitch = ( Word16 )( st->pitch & 0x3fff );
   }

   if ( pow_sum < POW_COMPLEX_THR ) {
      st->complex_low = ( Word16 )( st->complex_low & 0x3fff );
   }

   filter_bank( st, in_buf, level );
   return( vad_decision( st, level, pow_sum ) );
}
#endif

#ifndef VAD2
static void vad_pitch_detection( vadState *st, Word32 T_op[] )
{
   Word32 lagcount, i;

   lagcount = 0;

   for ( i = 0; i < 2; i++ ) {
      if ( abs( st->oldlag - T_op[i] ) < LTHRESH ) {
         lagcount += 1;
      }

      st->oldlag = T_op[i];
   }

   st->pitch = st->pitch >> 1;

   if ( ( st->oldlag_count + lagcount ) >= NTHRESH ) {
      st->pitch = st->pitch | 0x4000;
   }

   st->oldlag_count = lagcount;
}
#endif


#ifdef VAD2

int vad2 (vadState *st, Float32 *farray_ptr)
{
  static int	ch_tbl [NUM_CHAN][2] = {

    { 2,  3},
    { 4,  5},
    { 6,  7},
    { 8,  9},
    {10, 11},
    {12, 13},
    {14, 16},
    {17, 19},
    {20, 22},
    {23, 26},
    {27, 30},
    {31, 35},
    {36, 41},
    {42, 48},
    {49, 55},
    {56, 63}

  };

  static int	vm_tbl [90] = {
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    3, 3, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 7, 7, 7,
    8, 8, 9, 9, 10, 10, 11, 12, 12, 13, 13, 14, 15,
    15, 16, 17, 17, 18, 19, 20, 20, 21, 22, 23, 24,
    24, 25, 26, 27, 28, 28, 29, 30, 31, 32, 33, 34,
    35, 36, 37, 37, 38, 39, 40, 41, 42, 43, 44, 45,
    46, 47, 48, 49, 50, 50, 50, 50, 50, 50, 50, 50,
    50, 50
  };

  static Word16 hangover_table[20] =
  {
    30, 30, 30, 30, 30, 30, 28, 26, 24, 22, 20, 18, 16, 14, 12, 10, 8, 8, 8, 8
  };

  static Word16 burstcount_table[20] =
  {
    8, 8, 8, 8, 8, 8, 8, 8, 7, 6, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4
  };

  static Word16 vm_threshold_table[20] =
  {
    34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 40, 51, 71, 100, 139, 191, 257, 337, 432
  };

  float		data_buffer [FFT_LEN1], enrg, snr;
  float		tne, tce, ftmp;
  int		ch_snr [NUM_CHAN];
  int		i, j, j1, j2;
  int		vm_sum;
  int		update_flag;
  float		ch_enrg_dev;
  float		ch_enrg_db [NUM_CHAN];
  float		alpha;
  float		peak, avg, peak2avg;
  int		sine_wave_flag;
  float		tce_db, tne_db;
  float		xt;
  int	tsnrq;
  int	ivad;

  st->Lframe_cnt++;

  for (i = 0; i < DELAY0; i++)
    data_buffer [i] = 0.0;

  data_buffer [DELAY0] = *farray_ptr + PRE_EMP_FAC1 * st->pre_emp_mem;

  for (i = DELAY0+1, j = 1; i < DELAY0+FRM_LEN1; i++, j++)
    data_buffer [i] = *(farray_ptr + j) + PRE_EMP_FAC1 *
      *(farray_ptr + j - 1);

  st->pre_emp_mem = *(farray_ptr + FRM_LEN1 - 1);

  for (i = DELAY0+FRM_LEN1; i < FFT_LEN1; i++)
    data_buffer [i] = 0.0;
  real_fft (data_buffer, +1);
  alpha = (st->Lframe_cnt == 1) ? 1.0 : CEE_SM_FAC1;
  for (i = LO_CHAN; i <= HI_CHAN; i++)
    {
      enrg = 0.0;
      j1 = ch_tbl [i][0], j2 = ch_tbl [i][1];
      for (j = j1; j <= j2; j++)
	enrg += square(data_buffer [2*j]) + square(data_buffer [2*j+1]);
      enrg /= (float) (j2 - j1 + 1);
      st->ch_enrg [i] = (1 - alpha) * st->ch_enrg [i] + alpha * enrg;
      if (st->ch_enrg [i] < MIN_CHAN_ENRG) st->ch_enrg [i] = MIN_CHAN_ENRG;
    }

  tce = 0.0;
  for (i = LO_CHAN; i <= HI_CHAN; i++)
    tce += st->ch_enrg [i];

  peak = avg = 0.;
  for (i = LO_CHAN; i <= HI_CHAN; i++) {
    if (i >= SINE_START_CHAN && st->ch_enrg [i] > peak)
      peak = st->ch_enrg [i];
    avg += st->ch_enrg [i];
  }
  avg /= HI_CHAN - LO_CHAN + 1;
  peak2avg = (avg < 1./NORM_ENRG) ? 0. : 10.*log10 (peak/avg);

  if (peak2avg > 10.)
    sine_wave_flag = TRUE;
  else
    sine_wave_flag = FALSE;

  if (st->Lframe_cnt <= INIT_FRAMES) {
    if (sine_wave_flag == TRUE) {
      for (i = LO_CHAN; i <= HI_CHAN; i++)
        st->ch_noise [i] = INE;
    }
    else {
      for (i = LO_CHAN; i <= HI_CHAN; i++)
        st->ch_noise [i] = max(st->ch_enrg [i], INE);
    }
  }

  for (i = LO_CHAN; i <= HI_CHAN; i++) {
    snr = 10.0 * log10 ((double)st->ch_enrg [i] / st->ch_noise [i]);
    if (snr < 0.0) snr = 0.0;
    ch_snr [i] = (snr + 0.1875) / 0.375;
  }

  vm_sum = 0;
  for (i = LO_CHAN; i <= HI_CHAN; i++) {
    j = min(ch_snr[i],89);
    vm_sum += vm_tbl [j];
  }

  if (st->Lframe_cnt <= INIT_FRAMES  || st->fupdate_flag == TRUE ) {
#if NORM_ENERG==4
    tce_db = 49.918;
#elif NORM_ENERG==1
    tce_db = 55.938;
#else
    tce_db = (96. - 22. - 10*log10 (FFT_LEN1/2) - 10.*log10 (NORM_ENRG));
#endif

    st->negSNRvar = 0.0;
    st->negSNRbias = 0.0;

    tne = 0.0;
    for (i = LO_CHAN; i <= HI_CHAN; i++)
      tne += st->ch_noise [i];

    tne_db = 10 * log10 (tne);

    xt = tce_db - tne_db;
    st->tsnr = xt;

  }
  else {
    xt = 0;
    for (i=LO_CHAN; i<=HI_CHAN; i++)
      xt += st->ch_enrg[i]/st->ch_noise[i];
    xt = 10*log10(xt/NUM_CHAN);

    if (xt > st->tsnr)
      st->tsnr = 0.9*st->tsnr + 0.1*xt;
    else if (xt > 0.625*st->tsnr)
      st->tsnr = 0.998*st->tsnr + 0.002*xt;
  }

  tsnrq = (int)(st->tsnr/3.);
  tsnrq = min(19, max(0, tsnrq));

  if (xt < 0) {
    st->negSNRvar = min (0.99*st->negSNRvar + 0.01*xt*xt, 4.0);
    st->negSNRbias = max (12.0*(st->negSNRvar - 0.65), 0.0);
  }

  if (vm_sum > vm_threshold_table[tsnrq] + st->negSNRbias) {
    ivad = 1;
    if (++st->burstcount > burstcount_table[tsnrq]) {
      st->hangover = hangover_table[tsnrq];
    }
  } else {
    st->burstcount = 0;
    if (--st->hangover <= 0) {
      ivad = 0;
      st->hangover = 0;
    } else {
      ivad = 1;
    }
  }

  for (i = LO_CHAN; i <= HI_CHAN; i++)
    ch_enrg_db [i] = 10.*log10( st->ch_enrg [i] );

  ch_enrg_dev = 0.;
  if (st->Lframe_cnt == 1)
    for (i = LO_CHAN; i <= HI_CHAN; i++)
      st->ch_enrg_long_db [i] = ch_enrg_db [i];
  else
    for (i = LO_CHAN; i <= HI_CHAN; i++)
      ch_enrg_dev += fabs( st->ch_enrg_long_db [i] - ch_enrg_db [i] );

  ftmp = st->tsnr - xt;
  if (ftmp <= 0.0 || st->tsnr <= 0.0)
    alpha = HIGH_ALPHA1;
  else if (ftmp > st->tsnr)
    alpha = LOW_ALPHA1;
  else
    alpha = HIGH_ALPHA1 - (ALPHA_RANGE1 * ftmp / st->tsnr);

  for (i = LO_CHAN; i <= HI_CHAN; i++) {
    st->ch_enrg_long_db[i] = alpha*st->ch_enrg_long_db[i] + (1.-alpha)*ch_enrg_db[i];
  }

  update_flag = FALSE;
  st->fupdate_flag = FALSE;
  if ((vm_sum <= UPDATE_THLD) ||
      (st->Lframe_cnt <= INIT_FRAMES && sine_wave_flag == FALSE)) {
    update_flag = TRUE;
    st->update_cnt = 0;
  }
  else if (tce > NOISE_FLOOR && ch_enrg_dev < DEV_THLD1 &&
           sine_wave_flag == FALSE && st->LTP_flag == FALSE) {
    st->update_cnt++;
    if (st->update_cnt >= UPDATE_CNT_THLD1) {
      update_flag = TRUE;
      st->fupdate_flag = TRUE;
    }
  }

  if ( st->update_cnt == st->last_update_cnt )
    st->hyster_cnt++;
  else
    st->hyster_cnt = 0;
  st->last_update_cnt = st->update_cnt;

  if ( st->hyster_cnt > HYSTER_CNT_THLD1 )
    st->update_cnt = 0;

  /* Update the channel noise estimates */
  if (update_flag == TRUE) {
    for (i = LO_CHAN; i <= HI_CHAN; i++) {
      st->ch_noise [i] = (1.0 - CNE_SM_FAC1) * st->ch_noise [i] +
	CNE_SM_FAC1 * st->ch_enrg [i];
      if (st->ch_noise [i] < MIN_CHAN_ENRG) st->ch_noise [i] = MIN_CHAN_ENRG;
    }
  }

  return (ivad);

}

static double	phs_tbl [SIZE];		/* holds the complex sinusoids */

void		real_fft (float *farray_ptr, int isign)
{

  float		ftmp1_real, ftmp1_imag, ftmp2_real, ftmp2_imag;
  int		i, j;
  static int	first = TRUE;

  void		cmplx_fft (float *, int);
  void		fill_tbl ();

  if (first == TRUE) {
    fill_tbl ();
    first = FALSE;
  }

  if (isign == 1) {
    cmplx_fft (farray_ptr, isign);
    ftmp1_real = *farray_ptr;
    ftmp2_real = *(farray_ptr + 1);
    *farray_ptr = ftmp1_real + ftmp2_real;
    *(farray_ptr + 1) = ftmp1_real - ftmp2_real;

    for (i = 2, j = SIZE - i; i <= SIZE_BY_TWO; i = i + 2, j = SIZE - i) {
      ftmp1_real = *(farray_ptr + i) + *(farray_ptr + j);
      ftmp1_imag = *(farray_ptr + i + 1) - *(farray_ptr + j + 1);
      ftmp2_real = *(farray_ptr + i + 1) + *(farray_ptr + j + 1);
      ftmp2_imag = *(farray_ptr + j) - *(farray_ptr + i);

      *(farray_ptr + i) = (ftmp1_real + phs_tbl [i] * ftmp2_real - phs_tbl [i + 1] * ftmp2_imag) / 2.0;
      *(farray_ptr + i + 1) = (ftmp1_imag + phs_tbl [i] * ftmp2_imag + phs_tbl [i + 1] * ftmp2_real) / 2.0;
      *(farray_ptr + j) = (ftmp1_real + phs_tbl [j] * ftmp2_real + phs_tbl [j + 1] * ftmp2_imag) / 2.0;
      *(farray_ptr + j + 1) = (-ftmp1_imag - phs_tbl [j] * ftmp2_imag + phs_tbl [j + 1] * ftmp2_real) / 2.0;
    }
  }
  else {
    ftmp1_real = *farray_ptr;
    ftmp2_real = *(farray_ptr + 1);
    *farray_ptr = (ftmp1_real + ftmp2_real) / 2.0;
    *(farray_ptr + 1) = (ftmp1_real - ftmp2_real) / 2.0;

    for (i = 2, j = SIZE - i; i <= SIZE_BY_TWO; i = i + 2, j = SIZE - i) {
      ftmp1_real = *(farray_ptr + i) + *(farray_ptr + j);
      ftmp1_imag = *(farray_ptr + i + 1) - *(farray_ptr + j + 1);
      ftmp2_real = -(*(farray_ptr + i + 1) + *(farray_ptr + j + 1));
      ftmp2_imag = -(*(farray_ptr + j) - *(farray_ptr + i));

      *(farray_ptr + i) = (ftmp1_real + phs_tbl [i] * ftmp2_real + phs_tbl [i + 1] * ftmp2_imag) / 2.0;
      *(farray_ptr + i + 1) = (ftmp1_imag + phs_tbl [i] * ftmp2_imag - phs_tbl [i + 1] * ftmp2_real) / 2.0;
      *(farray_ptr + j) = (ftmp1_real + phs_tbl [j] * ftmp2_real - phs_tbl [j + 1] * ftmp2_imag) / 2.0;
      *(farray_ptr + j + 1) = (-ftmp1_imag - phs_tbl [j] * ftmp2_imag - phs_tbl [j + 1] * ftmp2_real) / 2.0;
    }
    cmplx_fft (farray_ptr, isign);
  }
}

void		cmplx_fft (float *farray_ptr, int isign)
{
  int		i, j, k, ii, jj, kk, ji, kj;
  float		ftmp, ftmp_real, ftmp_imag;

  for (i = 0, j = 0; i < SIZE-2; i = i + 2) {
    if (j > i) {
      ftmp = *(farray_ptr+i);
      *(farray_ptr+i) = *(farray_ptr+j);
      *(farray_ptr+j) = ftmp;

      ftmp = *(farray_ptr+i+1);
      *(farray_ptr+i+1) = *(farray_ptr+j+1);
      *(farray_ptr+j+1) = ftmp;
    }
    k = SIZE_BY_TWO;
    while (j >= k) {
      j -= k;
      k >>= 1;
    }
    j += k;
  }

  if (isign == 1) {
    for (i = 0; i < NUM_STAGE; i++) {		/* i is stage counter */
      jj = (2 << i);				/* FFT size */
      kk = (jj << 1);				/* 2 * FFT size */
      ii = SIZE / jj;				/* 2 * number of FFT's */
      for (j = 0; j < jj; j = j + 2) {		/* j is sample counter */
        ji = j * ii;				/* ji is phase table index */
        for (k = j; k < SIZE; k = k + kk) {	/* k is butterfly top */
          kj = k + jj;				/* kj is butterfly bottom */

          ftmp_real = *(farray_ptr + kj) * phs_tbl [ji] -
	    *(farray_ptr + kj + 1) * phs_tbl [ji + 1];

          ftmp_imag = *(farray_ptr + kj + 1) * phs_tbl [ji] +
	    *(farray_ptr + kj) * phs_tbl [ji + 1];

          *(farray_ptr + kj) = (*(farray_ptr + k) - ftmp_real) / 2.0;
          *(farray_ptr + kj + 1) = (*(farray_ptr + k + 1) - ftmp_imag) / 2.0;

          *(farray_ptr + k) = (*(farray_ptr + k) + ftmp_real) / 2.0;
          *(farray_ptr + k + 1) = (*(farray_ptr + k + 1) + ftmp_imag) / 2.0;
        }
      }
    }
  }
  else {
    for (i = 0; i < NUM_STAGE; i++) {		/* i is stage counter */
      jj = (2 << i);				/* FFT size */
      kk = (jj << 1);				/* 2 * FFT size */
      ii = SIZE / jj;				/* 2 * number of FFT's */
      for (j = 0; j < jj; j = j + 2) {		/* j is sample counter */
        ji = j * ii;				/* ji is phase table index */
        for (k = j; k < SIZE; k = k + kk) {	/* k is butterfly top */
          kj = k + jj;				/* kj is butterfly bottom */

          ftmp_real = *(farray_ptr + kj) * phs_tbl [ji] +
	    *(farray_ptr + kj + 1) * phs_tbl [ji + 1];

          ftmp_imag = *(farray_ptr + kj + 1) * phs_tbl [ji] -
	    *(farray_ptr + kj) * phs_tbl [ji + 1];

          *(farray_ptr + kj) = *(farray_ptr + k) - ftmp_real;
          *(farray_ptr + kj + 1) = *(farray_ptr + k + 1) - ftmp_imag;

          *(farray_ptr + k) = *(farray_ptr + k) + ftmp_real;
          *(farray_ptr + k + 1) = *(farray_ptr + k + 1) + ftmp_imag;
        }
      }
    }
  }
}


void		fill_tbl ()
{
  int		i;
  double	delta_f, theta;

  delta_f = - PI / (double) SIZE_BY_TWO;
  for (i = 0; i < SIZE_BY_TWO; i++) {
    theta = delta_f * (double) i;
    phs_tbl[2*i] = cos(theta);
    phs_tbl[2*i+1] = sin(theta);
  }
  return;
}		/* end fill_tbl () */

void LTP_flag_update (vadState * st, Word16 mode)
{
  Float32 thresh;

  if ((mode == MR475) || (mode == MR515))
    thresh = 0.55;
  else if (mode == MR102)
    thresh = 0.60;
  else
    thresh = 0.65;

  if (st->Rmax  > thresh*st->R0) st->LTP_flag = TRUE;
  else st->LTP_flag = FALSE;

  return;
}

#endif
