#ifndef _SP_ENC_H
#define _SP_ENC_H

#include "typedef.h"
#include "sp_runtime.h"

typedef struct
{
   Float32 y2;
   Float32 y1;
   Float32 x0;
   Float32 x1;
}
Pre_ProcessState;

#ifdef VAD2

#define	FRM_LEN1		80
#define	DELAY0			24
#define	FFT_LEN1		128

#define	UPDATE_CNT_THLD1	50

#define	INIT_FRAMES		4

#define	CNE_SM_FAC1		0.1
#define	CEE_SM_FAC1		0.55

#define	HYSTER_CNT_THLD1	6
#define	HIGH_ALPHA1		0.9
#define	LOW_ALPHA1		0.7
#define	ALPHA_RANGE1		(HIGH_ALPHA1-LOW_ALPHA1)

#define NORM_ENRG		(4.0)
#define	MIN_CHAN_ENRG		(0.0625 / NORM_ENRG)
#define	INE			(16.0 / NORM_ENRG)
#define	NOISE_FLOOR		(1.0 / NORM_ENRG)

#define	PRE_EMP_FAC1		(-0.8)

#define	NUM_CHAN		16
#define	LO_CHAN			0
#define	HI_CHAN			15
#define	UPDATE_THLD		35

#define	SINE_START_CHAN		2
#define	P2A_THRESH		10.0
#define	DEV_THLD1		28.0

#define	SIZE			128
#define	SIZE_BY_TWO		64
#define	NUM_STAGE		6

#define	PI			3.141592653589793

#define	TRUE			1
#define	FALSE			0

#define	min(a,b)		((a)<(b)?(a):(b))
#define	max(a,b)		((a)>(b)?(a):(b))
#define	square(a)		((a)*(a))

typedef struct
{
  Float32 pre_emp_mem;
  Word16  update_cnt;
  Word16  hyster_cnt;
  Word16  last_update_cnt;
  Float32 ch_enrg_long_db[NUM_CHAN];
  Word32  Lframe_cnt;
  Float32 ch_enrg[NUM_CHAN];
  Float32 ch_noise[NUM_CHAN];
  Float32 tsnr;
  Word16  hangover;
  Word16  burstcount;
  Word16  fupdate_flag;
  Float32 negSNRvar;
  Float32 negSNRbias;
  Float32 R0;
  Float32 Rmax;
  Word16  LTP_flag;
}
vadState;
#else
typedef struct
{
   Float32 bckr_est[COMPLEN];   /* background noise estimate */
   Float32 ave_level[COMPLEN];
   Float32 old_level[COMPLEN];   /* input levels of the previous frame */
   Float32 sub_level[COMPLEN];
   Float32 a_data5[3][2];   /* memory for the filter bank */
   Float32 a_data3[5];   /* memory for the filter bank */
   Float32 best_corr_hp;   /* FIP filtered value */
   Float32 corr_hp_fast;   /* filtered value */
   Word32 vadreg;   /* flags for intermediate VAD decisions */
   Word32 pitch;   /* flags for pitch detection */
   Word32 oldlag_count, oldlag;   /* variables for pitch detection */
   Word32 complex_high;   /* flags for complex detection */
   Word32 complex_low;   /* flags for complex detection */
   Word32 complex_warning;   /* complex background warning */
   Word32 tone;   /* flags for tone detection */
   Word16 burst_count;   /* counts length of a speech burst */
   Word16 hang_count;   /* hangover counter */
   Word16 stat_count;   /* stationary counter */
   Word16 complex_hang_count;   /* complex hangover counter, used by VAD */
   Word16 complex_hang_timer;   /* hangover initiator, used by CAD */
   Word16 speech_vad_decision;   /* final decision */
   Word16 sp_burst_count;
}
vadState;
#endif
#define DTX_HIST_SIZE 8
#define DTX_ELAPSED_FRAMES_THRESH (24 + 7 -1)
#define DTX_HANG_CONST 7   /* yields eight frames of SP HANGOVER */
typedef struct
{
   Float32 lsp_hist[M * DTX_HIST_SIZE];
   Float32 log_en_hist[DTX_HIST_SIZE];
   Word32 init_lsf_vq_index;
   Word16 hist_ptr;
   Word16 log_en_index;
   Word16 lsp_index[3];
   Word16 dtxHangoverCount;
   Word16 decAnaElapsedCount;
}
dtx_encState;

typedef struct
{
   Float32 gp[N_FRAME];
   Word16 count;
}
tonStabState;

typedef struct
{
   Word32 past_qua_en[4];
}
gc_predState;

typedef struct
{
   Float32 prev_alpha;
   Float32 prev_gc;
   Float32 ltpg_mem[LTPG_MEM_SIZE];
   Word16 onset;
}
gain_adaptState;

typedef struct
{
   Float32 sf0_target_en;
   Float32 sf0_coeff[5];
   Word32 sf0_gcode0_exp;
   Word32 sf0_gcode0_fra;
   Word16 *gain_idx_ptr;
   gc_predState * gc_predSt;
   gc_predState * gc_predUncSt;
   gain_adaptState * adaptSt;
}
gainQuantState;
typedef struct
{
   Word32 T0_prev_subframe;
}
Pitch_frState;
typedef struct
{
   Pitch_frState * pitchSt;
}
clLtpState;
typedef struct
{
   Float32 ada_w;
   Word32 old_T0_med;
   Word16 wght_flg;
}
pitchOLWghtState;
typedef struct
{
   Float32 past_rq[M];
}
Q_plsfState;
typedef struct
{
   Float32 lsp_old[M];
   Float32 lsp_old_q[M];
   Q_plsfState * qSt;
}
lspState;
typedef struct
{
   Float32 old_A[M + 1];
}
LevinsonState;
typedef struct
{
   LevinsonState * LevinsonSt;
}
lpcState;
typedef struct
{
   Float32 old_speech[L_TOTAL];
   Float32 *speech, *p_window, *p_window_12k2;
   Float32 *new_speech;   /* Global variable */
   Float32 old_wsp[L_FRAME + PIT_MAX];
   Float32 *wsp;
   Word32 old_lags[5];
   Float32 ol_gain_flg[2];
   Float32 old_exc[L_FRAME + PIT_MAX + L_INTERPOL];
   Float32 *exc;
   Float32 ai_zero[L_SUBFR + MP1];
   Float32 *zero;
   Float32 *h1;
   Float32 hvec[L_SUBFR * 2];
   lpcState * lpcSt;
   lspState * lspSt;
   clLtpState * clLtpSt;
   gainQuantState * gainQuantSt;
   pitchOLWghtState * pitchOLWghtSt;
   tonStabState * tonStabSt;
   vadState * vadSt;
   Word32 dtx;
   dtx_encState * dtxEncSt;
   Float32 mem_syn[M], mem_w0[M], mem_w[M];
   Float32 mem_err[M + L_SUBFR], *error;
   Float32 sharp;
}
cod_amrState;
typedef struct
{
   cod_amrState * cod_amr_state;
   Pre_ProcessState * pre_state;
   Word32 dtx;
}Speech_Encode_FrameState;

#ifndef ModeNBDefined

enum ModeNB { MR475 = 0,
            MR515,
            MR59,
            MR67,
            MR74,
            MR795,
            MR102,
            MR122,
            MRDTX
};
#define ModeNBDefined
#endif

void *Speech_Encode_Frame_init (int dtx);
int Speech_Encode_Frame_reset(void *st, int dtx);
void Speech_Encode_Frame_exit (void **st);
void Speech_Encode_Frame (void *st, enum ModeNB mode, short *newSpeech,
                   short *prm, enum ModeNB *usedMode);

#endif
