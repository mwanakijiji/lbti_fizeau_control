/* LBTI Path Length Corrector.
 * Read phase camera and issue corrections to Fast Path-Length Corrector PZTs, or NAC PZTs if testing.
 *
 * Main program thread, sometimes referred to as "background", handles all INDI messages, disk IO and
 * other relatively slow periodic tasks via callbacks. Two more threads perform image processing, one for
 * each beam. The thread for beam 2 just performs image processing on beam 2. The thread for beam 1,
 * sometimes referred to as the "real time" thread (for deprecated historical reasons), performs image
 * processing on beam 1, then combines its results with beam 2 and computes and issues the PZT commands.
 * A fourth thread reads OVMS messages.
 * 
 * Both image processing threads run all the time. The concept of "loop closed" just determines whether
 * commands are sent to the beam steering PZTs.
 * 
 * The camera image consists of two light beams, processing for each is the same:
 *   1. extract beam based on operator-defined red (beam1) or green (beam2) circle (see GUI)
 *   2. pad to larger background with mean to increase FFT resolution.
 *   3. compute complex FFT: mag used to compute tip and tilt corrections, arg for phase correction.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdarg.h>
#include <dirent.h>
#include <float.h>
#include <math.h>
#include <errno.h>
#include <ctype.h>
#include <time.h>
#include <fcntl.h>
#include <unistd.h>
#include <pthread.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <netdb.h>
#include <arpa/inet.h>
#include <fftw3.h>

#include "ip.h"
#include "fq.h"
#include "zlib.h"
#include "astro.h"
#include "indidevapi.h"
#include "configfile.h"
#include "indiplcprops.h"
#include "../mkfilter/biquad.h"

/* should have been in netdb.h ? */
extern int h_errno;

/* compile-time config setup */
#define	EXITSLEEP	2		/* seconds to wait before exiting */
static char plc_cfgfn[] = "plc.cfg";    /* name of our cfg file */

/* variables from config file
 */
static char *pcam_host;			/* phase camera imaging server host */
static int pcam_port;			/* phase camera imaging server port */
static char *pzt_host;			/* path length correcter server host */
static int pzt_port;			/* path length corrector server port */
static int trace;			/* whether to log debug info */
static int nopchw;			/* no phasecam hw */
static int nofpchw;			/* no FPC hw */
static int noovms;			/* no OVMS */
static double max_pl_corr;		/* max phase correction in one loop, degrees */
static double max_tip_corr;		/* max tip correction in one loop, mas */
static double max_tip_err;		/* max tip error, mas */
static double max_tilt_corr;		/* max tilt correction in one loop, mas */
static double max_tilt_err;		/* max tilt error, mas */
static char *darks_dir;			/* dir for darks */
static char *flats_dir;			/* dir for flats */
static char *badpix_dir;		/* dir for bad pixel maps */
static int nflatimg;			/* n images per flat */
static int badmapsecs;			/* n seconds per bad pixel map */
static double bad_space_std;		/* bad spacial pixel std dev */
static double bad_time_std;		/* bad temporal pixel std dev */
static double cg_min;			/* threshold between using phase or CG */
static double spc_offload_ms;		/* move spc this often, ms */
static double spc_offload_scale;	/* steps per micron */
static double spc_offload_alpha;	/* exponential smoothing parameter */
static double spc_offload_max;		/* max allowed piston, um */
static char *spc_offload_motnam;	/* SPC motor INDI name */
static char *spc_offload_devnam;	/* SPC INDI device name */
static int sat_pix;			/* smallest pixel value considered saturated, 0..64k */
static double ramp_time;		/* gain ramp time after closing loop, secs */
static double b2sp0_time;		/* time to collect beam2 setpoint baseline phase, secs */
static int beam_pad0;			/* beam magnification factor from config file */
static double rtcolper;			/* default rtdata collection period, seconds */
static double rtpubper;			/* rtdata publishing period, seconds */
static int nautobeam;                   /* n frames history for find beam stats */
static double cg_boxcar_period;		/* peropd of CG average, seconds */
static double wv_boxcar_period;		/* peropd of water vapor average, seconds */
static double plpfc;			/* piston low pass filter cutoff freq, Hz */
static double cg_freeze;		/* seconds to freeze PIDs after CG kick */
static char *camera_host;		/* camera info host */
static int camera_port;			/* camera info port */
static int camera_period;		/* period to poll camera info, ms */
static double lam1_PC;			/* nominal Phasecam beam 1 operating wavelength, microns */
static double lam2_PC;			/* nominal Phasecam beam 2 operating wavelength, microns */
static double lam_NOMIC;		/* nominal NOMIC operating wavelength, microns */
static double peakignore;		/* fraction of mag image to ignore looking for tip/tilt peaks */
static int cpu_plcloop;			/* CPU affinity for plc loop */
static int cpu_beam2;			/* CPU affinity for beam2 loop */
static double opd_bump;			/* OVMS OPD change considered too much */
static double opd_freeze;		/* freeze loop time, secs */
static double fscan_min;		/* minimum sweep position, microns from reference pos */
static double fscan_max;		/* maximum sweep position, microns from reference pos */
static double fscan_step;		/* microns per step */
static int fscan_nsweeps;		/* number of full sweeps to perform, 0 to disable */
static int fscan_settle;		/* microsecs to wait after each step to allow for PZT settling */
static double UBC_PLPGain_SF;		/* gain scale factor */
static double UBC_PLIGain_SF;		/* gain scale factor */
static double UBC_PLDGain_SF;		/* gain scale factor */
static double UBC_TTPGain_SF;		/* gain scale factor */
static double UBC_TTIGain_SF;		/* gain scale factor */
static double UBC_TTDGain_SF;		/* gain scale factor */
static double NAC_PLPGain_SF;		/* gain scale factor */
static double NAC_PLDGain_SF;		/* gain scale factor */
static double NAC_PLIGain_SF;		/* gain scale factor */
static double NAC_TTPGain_SF;		/* gain scale factor */
static double NAC_TTDGain_SF;		/* gain scale factor */
static double NAC_TTIGain_SF;		/* gain scale factor */
static double OVMS_TIP_GAIN_SF;		/* gain scale factor */
static double OVMS_TILT_GAIN_SF;	/* gain scale factor */
static double OVMS_OPD_GAIN_SF;		/* gain scale factor */



/* info to collect closed-loop operating statistics.
 * this is filled in by the real-time loops, the background loop copies occasionally to publish.
 * access is protected with the mutex.
 */
typedef struct {

    // accessor guard
    pthread_mutex_t lock;

    // operational parameters
    int plc_loops;			/* n times through plc thread */
    double t_start;			/* when these samples were started */
    int n_samples;			/* number of samples */	
    int nf_jumps;			/* number of CG frame jumps in this period */

    /******
     * measured values
     */

    // PZT commanded values
    double p_sum, p_sum2;		/* sum(piston) and sum^2(piston) */
    double a_sum, a_sum2;		/* sum(az==tilt) and sum^2(az==tilt) */
    double e_sum, e_sum2;		/* sum(el==tip) and sum^2(el==tip) */

    // beam 1
    double ph1_sum, ph1_sum2;		/* sum(unwrap_phase) and sum^2(unwrap_phase) */
    double cph1_sum, sph1_sum;		/* sum(cos(unwrap_phase)) and sum(sin(unwrap_phase)) */
    double snr1_sum;			/* sum(SNR) */

    // beam 2
    double ph2_sum, ph2_sum2;		/* sum(unwrap_phase) and sum^2(unwrap_phase) */
    double cph2_sum, sph2_sum;		/* sum(cos(unwrap_phase)) and sum(sin(unwrap_phase)) */
    double snr2_sum;			/* sum(SNR) */

    /******
     * corresponding derived values for publishing
     */

    double p_mean, p_sd;		/* piston mean, std dev */
    double a_mean, a_sd;		/* az==tilt mean, std dev */
    double e_mean, e_sd;		/* el==tip mean, std dev */

    double ph1_mean;			/* unw_phase mean */
    double ph1_sd;			/* unw_phase std dev */
    double cph1_mean, sph1_mean;	/* mean cos(unwrap_phase) and sin(unwrap_phase) */
    double snr1_mean;			/* SNR mean */

    double ph2_mean;			/* unw_phase mean */
    double ph2_sd;			/* unw_phase std dev */
    double cph2_mean, sph2_mean;	/* mean cos(unwrap_phase) and sin(unwrap_phase) */
    double snr2_mean;			/* SNR mean */

} LoopStats;

static LoopStats loopstats;		/* stats */
static int plc_period;			/* stat accumulation period, ms, set from EXPTIME from sci camera */


/* The following variables support sending back real-time data samples via the pcamdata INDI property.
 *   The name of the FITS header field and the duration of samples to capture is set by the pcamdatafield
 *   property. The real-time loop just adds new entries to the ring buffers at the position of rtdata_i.
 *   The background thread publishes the buffers guarded with plc_bglock.
 */
static const int MAXRTSAMPS = 10000;	/* maximum samples */
static double rtdata_t[MAXRTSAMPS];	/* sample time, JD */
static double rtdata_y[MAXRTSAMPS];	/* measured samples */
static double rtdata_dt;		/* max time sample period, seconds */
static int rtdata_n;			/* n entries actually in use in each ring buffer */
static int rtdata_i;			/* index of next rtdata_[] buffer entry to use */
static char jd_kw[] = "PCJD";		/* FITS keyword for current time as JD */


/* handy shorthands for some INDI properties.
 * the strange forms using dereferences based on USE_UBC allow the defines to be used as lvalues.
 * N.B. locations are FITS 1-based
 */
#define	SAVEFEVERY	atoi(savefits.tp[EVERY_SVF].text)
#define	SAVESTATS	(savestats.sp[YES_SVST].s == ISS_ON)
#define	SAVEFDIR	savefits.tp[DIR_SVF].text
#define	SPCOFFLOAD	(spcoffload.sp[YES_SPCOL].s == ISS_ON)
#define	SENDFITS	sendfits.np[HZ_SNDF].value
#define	APPLYDARK	(applydark.sp[APPLY_AD].s == ISS_ON)
#define	APPLYFLAT	(applyflat.sp[APPLY_AF].s == ISS_ON)
#define	APPLYBADPIX	(applybadpix.sp[APPLY_AB].s == ISS_ON)
#define	USE_UBC		(useubc.sp[YES_UUBC].s == ISS_ON)
#define AUTOBEAM_ON     (findbeam.np[STATE_PFB].value > 0)
#define SHOWBEAM1	(showbeam1.sp[YES_SBM1].s == ISS_ON)
#define	LOOPCLOSED	(closeloop.sp[YES_CLP].s == ISS_ON)
#define	MINFSNR		(*(USE_UBC ? &ubcsettings.np[MINFSNR_UBCS].value  : &nacsettings.np[MINFSNR_NACS].value))
#define	PLSP		(*(USE_UBC ? &ubcsettings.np[PLSP_UBCS].value     : &nacsettings.np[PLSP_NACS].value))
#define	PLPGAIN		(USE_UBC   ? UBC_PLPGain_SF*ubcsettings.np[PLPGAIN_UBCS].value  : NAC_PLPGain_SF*nacsettings.np[PLPGAIN_NACS].value)
#define	PLDGAIN		(USE_UBC   ? UBC_PLDGain_SF*ubcsettings.np[PLDGAIN_UBCS].value  : NAC_PLDGain_SF*nacsettings.np[PLDGAIN_NACS].value)
#define	PLIGAIN		(USE_UBC   ? UBC_PLIGain_SF*ubcsettings.np[PLIGAIN_UBCS].value  : NAC_PLIGain_SF*nacsettings.np[PLIGAIN_NACS].value)
#define	CGSP		(*(USE_UBC ? &ubcsettings.np[CGSP_UBCS].value     : &nacsettings.np[CGSP_NACS].value))
#define	CGSCALE		(*(USE_UBC ? &ubcsettings.np[CGSC_UBCS].value     : &nacsettings.np[CGSC_NACS].value))
#define	PWVG		(*(USE_UBC ? &ubcsettings.np[PWVGAIN_UBCS].value  : &nacsettings.np[PWVGAIN_NACS].value))
#define	TIPSP		(*(USE_UBC ? &ubcsettings.np[TIPSP_UBCS].value    : &nacsettings.np[TIPSP_NACS].value))
#define	TILTSP		(*(USE_UBC ? &ubcsettings.np[TILTSP_UBCS].value   : &nacsettings.np[TILTSP_NACS].value))
#define	TTPGAIN		(USE_UBC   ? UBC_TTPGain_SF*ubcsettings.np[TTPGAIN_UBCS].value  : NAC_TTPGain_SF*nacsettings.np[TTPGAIN_NACS].value)
#define	TTDGAIN		(USE_UBC   ? UBC_TTDGain_SF*ubcsettings.np[TTDGAIN_UBCS].value  : NAC_TTDGain_SF*nacsettings.np[TTDGAIN_NACS].value)
#define	TTIGAIN		(USE_UBC   ? UBC_TTIGain_SF*ubcsettings.np[TTIGAIN_UBCS].value  : NAC_TTIGain_SF*nacsettings.np[TTIGAIN_NACS].value)



/* local copy of appropriate beam configuration to avoid changes during one plcOneIteration().
 * protect with beam_lock when these change.
 * N.B. beam positions are FITS coordinates
 */
static int BEAM1X, BEAM1Y;		/* beam 1 center */
static int BEAM2X, BEAM2Y;		/* beam 2 center */
static int BEAMR;			/* same radius; width = 2*BEAMR + 1 */
static int beam_pad;			/* current working beam pad, derived from beam_pad0 and BEAMR */
#define	BEAMW	 (BEAMR*2+1)		/* beam width on raw camera image, pixels */
#define	PADIMGSZ (BEAMW*beam_pad)	/* padded beam image size, pixels */
#define	BEAM_P2MAS(p) (440.0*(p)/PADIMGSZ)	/* convert pixels on FFT mag image to mas tip/tilt */
#define	BEAM_MAS2P(t) ((t)*PADIMGSZ/440.0)	/* convert mas tilt/tip to pixels on FFT mag image */
static pthread_mutex_t beam_lock;	/* guard for changing above values */


/* persistent piezo info to support plcLoopThread
 */
static FILE *pzt_f;			/* acromag socket */
static double pzt_tip0;			/* tip when loop closed */
static double pzt_tilt0;		/* tilt when loop closed */
static double pzt_pist0;		/* piston when loop closed */


/* misc variables */
static pthread_mutex_t plc_bglock;	/* used to lock copying from rt -> bg */
static int sf_tid;			/* SENDFITS timer ID */
static double Tip_Now, Tilt_Now, Piston_Now;	/* tip, tilt and piston, asec, asec, microns */
static double Piston_Mean;		/* Piston_Now smoothed with spc_offload_alpha */
static bool init_plloop;		/* see computeCorrections() */
static int spc_offload_tid;		/* spcoffloadCB() timer ID */
static double spc_offload_net;		/* microns since loop closed */


/* variables involved with camera calibrations */
static char dark_fn[128];		/* basename of dark file, if not empty */
static char flat_fn[128];		/* basename of flat file, if not empty */
static char badpix_fn[128];		/* basename of bad pixel map file, if not empty */
static FImage dark_fim;			/* dark image, if DARKISDEF */
static FImage flat_fim;			/* flat image, if FLATISDEF */
static FImage badpix_fim;		/* bad pixel map image, if FLATISDEF */
#define	DARKISDEF()  (dark_fn[0]!=0)	/* whether dark_fim is valid */
#define	DARKUNDEF()  (dark_fn[0]=0)	/* mark dark_fim as invalid */
#define	FLATISDEF()  (flat_fn[0]!=0)	/* whether flat_fim is valid */
#define	FLATUNDEF()  (flat_fn[0]=0)	/* mark flat_fim as invalid */
#define	BADPIXISDEF()  (badpix_fn[0]!=0)/* whether badpix_fim is valid */
#define	BADPIXUNDEF()  (badpix_fn[0]=0)	/* mark badpix_fim as invalid */


/* simple ring buffer for use when running stats are needed.
 * initialize with desired size.
 */
class RingBuffer {
    private:
	int bufsize, head, nbuf;
	double sum, sum2;
	double *buffer;

    public:
	RingBuffer (int size) {
	    bufsize = size;
	    buffer = (double *) malloc (bufsize*sizeof(double));
	    reset();
	}
	~RingBuffer() {
	    free (buffer);
	}
	void resize (int size) {
	    bufsize = size;
	    buffer = (double *) realloc (buffer, bufsize*sizeof(double));
	    reset();
	}
	void reset() {
	    memset (buffer, 0, bufsize*sizeof(double));
	    sum = sum2 = 0;
	    head = 0;
	    nbuf = 0;
	}
	void add (double v) {
	    double old = buffer[head];
	    sum -= old;
	    sum2 -= old*old;
	    buffer[head] = v;
	    sum += v;
	    sum2 += v*v;
	    if (++head >= bufsize)
		head = 0;
	    if (++nbuf > bufsize)
		nbuf = bufsize;
	}
	double rms() {
	    return (sqrt(sum2/nbuf));
	}
	double stddev() {
	    return (sqrt(nbuf*sum2 - sum*sum)/nbuf);
	}
	double mean() {
	    return (sum/nbuf);
	}
};
static RingBuffer *autobeam_x_rb;      /* used for autobeam left-right corrections */
static RingBuffer *autobeam_y_rb;      /* used for autobeam up-down corrections */

/* a simple box-car summing class that records how long it has been since it was reset
 */
class BoxCar {
    private:
	int nused;
	double sum;
	struct timeval tv0;

    public:
	BoxCar() {
	    reset();
	}
	void reset() {
	    nused = 0;
	    sum = 0;
	    gettimeofday (&tv0, NULL);
	}
	void add (double v) {
	    sum += v;
	    nused++;
	}
	int n() {
	    return (nused);
	}
	double age() {				// seconds since reset()
	    struct timeval tv;
	    gettimeofday (&tv, NULL);
	    return ((tv.tv_sec + tv.tv_usec*1e-6) - (tv0.tv_sec + tv0.tv_usec*1e-6));
	}
	double mean(void) {
	    return (nused > 0 ? sum/nused : 0);
	}
};

/* collection of FFT info, images and derived results from processing one beam.
 */
class OneBeam {

    private:

	bool unwphase_known;		// used to init unwrapped phase to first measurement

    public:

	/* mark unwrapped phase as unknown
	 */
	void initUnwrappedPhase() {
	    unwphase_known = false;
	}

	/* update unwrapped phase
	 */
	void updateUnwrappedPhase() {

	    // init phase_unw to phase if first time through
	    if (unwphase_known) {
		phase_unw += phase - last_phase;
		if (phase < last_phase - 180)
		    phase_unw += 360;
		else if (phase > last_phase + 180)
		    phase_unw -= 360;
	    } else {
		phase_unw = phase;
		unwphase_known = true;
	    }

	    // phase history
	    last_phase = phase;
	}

	/* These are used by the FFTW package for computing magnitude and phase of images containing fringes.
	 * They are set up each time we get a new beam location/size or when phasecam image size changes.
	 */
	int padsize;				/* current size */
	fftw_complex *fft_input;		/* FFT beam input array */
	fftw_complex *fft_output;		/* FFT beam output array */
	fftw_plan fft_theplan;			/* FFT "plan" */
	FImage pad_fim;				/* beam after extraction and padding */
	FImage fftmag_fim;			/* FFT magnitude image */
	FImage fftphs_fim;			/* FFT phase image */

	/* pointer to raw camera image
	 * (only way to get it to the beam2 thread)
	 */
	FImage *raw_fimp;

	/* wavelength for this beam, set from config file
	 */
	double lambda;

	/* reference setpoint.
	 * beam1 set from PLSP (from GUI), beam2 is set from initial phase measurement
	 */
	double sp0;

	/* values derived from the FFT
	 */
	CamPix tt_pixel;			/* raw tip/tilt max pixel */
	CamPix phs_pixel;			/* raw corresponding phase pixel */
	double setpoint;			/* net SP around which raw phase is remapped */
	double phase;				/* measured phase, +-180 from SP */
	double last_phase;			/* previous phase */
	double phase_unw;			/* unwrapped phase, degrees from start */
	double tt_snr;				/* SNR of tt_pixel */
	double tip, tilt;			/* beam separation from mag image, mas */

};

static OneBeam beam1, beam2;					/* info for each beam */

/* these are the raw, padded, FFT mag and phase images the background will publish.
 * copied from the selected beam from one or the other beam[12].
 * N.B. protect the copy with plc_bglock
 */
static FImage bg_cam_fim;		/* bg copy of raw camera image */
static FImage bg_pad_fim;		/* bg copy of padded beam */
static FImage bg_fftmag_fim;		/* bg copy of FFT mag */
static FImage bg_fftphs_fim;		/* bg copy of FFT phase */


/* these variables are used to synchronize the thread that processes beam2.
 * the thread loops forever, but each time around waits for b2_run before processing a new image,
 * then when the processing is complete it signals back using b2_done.
 */
static pthread_cond_t b2_run_cv;	/* cond variable to wait for new image to process */
static pthread_mutex_t b2_run_lock;	/* run CV's matching mutex */
static pthread_cond_t b2_done_cv;	/* cond variable to signal processing is complete */
static pthread_mutex_t b2_done_lock;	/* done CV's matching mutex */
static volatile bool b2_run = false;	/* set when ok for b2 thread to run */
static volatile bool b2_done = false;	/* set when measurements are ready */


/* to reduce impact on plcLoopThread(), FITS headers are written to disk by queing them
 * to fitshdr_q and letting them get written leisurely by saveFITSHeaderQCB().
 * access to the q is protected by fitshdrq_lock.
 */
typedef struct {
    FITSRow *var;			/* just like FImage var */
    int nvar;				/* just like FImage nvar */
    double jd;				/* time of this image */
} HDRVar;
static FQ *fitshdr_q;
static pthread_mutex_t fitshdrq_lock;
static const int saveFHms = 250;


/* this structure is used to record OVMS data.
 * the actual collection is performed in a separate thread so the lock must be used when
 * changing or acccess the data.
 */
typedef struct {

    /* the following are set once from the config file
     */
    int port;                           /* UDP port */
    char *group;                        /* multicast group */
    char *interface;                    /* network interface to use */

    /* gains: init from config file, updated via INDI, used by real-time thread.
     * don't bother protecting with a mutex.
     */
    double tip_gain, tilt_gain, opd_gain;

    /* the following are written by one thread and used by another so access must be protected with the mutex
     */
    pthread_mutex_t lock;               /* lock for the following items */
    double timestamp;                   /* UNIX time of the data, secs since 1/1/1970 UTC */
    double horizon;			/* prediction horizon, ms */
    double tip;                         /* raw tip */
    double tilt;                        /* raw tilt */
    double opd;                         /* raw opd */
    bool ok;                            /* whether getting valid packets ok */

    /* these check for sudden bumps, as per opd_bump and opd_freeze
     */
    bool bump;				/* whether currently in a bump condition */
    double bump_t0;			/* when bump started */

} OVMS;

static OVMS ovms;


/* path length vibration filter
 */
static BiQuadFilter vib_bqfilter;
#define VIB_FC          vibfilter.np[FC_PLVF].value                      // center freq, Hz
#define VIB_Q           vibfilter.np[Q_PLVF].value                       // Q ~ FC/BW
#define VIB_FGAIN       vibfilter.np[GAIN_PLVF].value                    // filter gain, dB
#define VIB_VGAIN       (pow (10.0, vibfilter.np[GAIN_PLVF].value/20))   // filter gain as voltage ratio
#define VIB_ZGAIN       vibgain.np[GAIN_PLVG].value                      // piezo gain


/* Simple PID controller.
 * For an intro see http://brettbeauregard.com/blog/2011/04/improving-the-beginners-pid-introduction
 */
class PIDController {

    private:

        double k_p, k_i, k_d;                           // current tuning values
        double i_sum;                                   // current integrator sum
        double on_time;                                 // when controller was turned on, secs since epoch
        double last_now;                                // time of last loop, secs since epoch
        double last_error;                              // last difference from set point
	double last_output;				// last output value
        double max_change;				// max allowed change from one step to the next
        double gain_ramp_time;                          // time to ramp up gains after turn on, seconds
	bool first_time;				// flag for first run
	BiQuadFilter difflp;				// low pass filter for differentiation
	bool use_difflp;				// whether to use

	// do not reset tuning values, ramp time or difflp, let them persist across turnOn's
	void initVariables() {
	    on_time = last_now = 0;
	    i_sum = 0;
	    last_error = 0;
	    max_change = DBL_MAX;
	    first_time = true;
	}

    public:

        // class constructor
        PIDController() {
	    initVariables();
	    k_p = k_i = k_d = 0;
	    gain_ramp_time = 0;
	    use_difflp = false;
        }

        // set desired tuning parameters
        void setTuning (double newkp, double newki, double newkd) {
            k_p = newkp;
            k_i = newki;
            k_d = newkd;

	    /* let next run collect new history */
	    first_time = true;
        }

	// set low pass filter for differentiator
	void setDiffLP (double Fc, double Fs, double Q) {
	    difflp.defineLowpassFilter (Fc, Fs, Q);
	    difflp.initRun();
	    use_difflp = true;
	}

        // set maximum step change, in units of the control variable
        void setMaxChange (double chg) {
            max_change = chg;
        }

        // set time over which gains are ramped up after turn on
        void setGainRampTime (double grt) {
            gain_ramp_time = grt;
        }

        // turn on the controller, given now and initial output value
        void turnOn (double now, double output) { 
	    initVariables();
            on_time = now;
            i_sum = output;
	    use_difflp = false;
        }

        // next iteration of controller, given time and set point - measurement
        double run (double now, double error) {

	    /* initialize history values on first call */
	    if (first_time) {
		first_time = false;
		last_error = error;
		last_now = now;
		last_output = i_sum;
		return (last_output);
	    }

            /* ramp up gains gradually over ramp_time */
            double gramp;
            if (gain_ramp_time > 0) {
                gramp = (now - on_time)/gain_ramp_time;
                if (gramp > 1)
                    gramp = 1;
            } else
                gramp = 1;
            double kp = k_p * gramp;
            double kd = k_d * gramp;
            double ki = k_i * gramp;

            /* error variables */
            double dt = (now - last_now);
            double de = (error - last_error) / dt;
	    if (use_difflp)
		de = difflp.run(de);
            i_sum += (ki * error) * dt;

            /* compute PID Output */
            double output = kp * error + i_sum + kd * de;
	    if (output < last_output - max_change)
		output = last_output - max_change;
	    else if (output > last_output + max_change)
		output = last_output + max_change;

            /* maintain state */
            last_error = error;
            last_now = now;
	    last_output = output;

            return (output);
        }
};

/* PID controllers for each driven variable */
static PIDController piston_pid;		// runs in units of microns
static PIDController tip_pid;			// runs in units of arc secs
static PIDController tilt_pid;			// runs in units of arc secs


/* local functions */
static void initOnce (void);
static void readConfig(void);
static int openSocket (const char *host, int port, bool log_success, char whynot[]);
static void noNagle (int fd);
static void *plcLoopThread (void *not_used);
static void plcOpenLoop ();
static int plcCloseLoop ();
static void startBeam2Thread(void);
static void *beam2Thread (void *not_used);
static void plcOneIteration(void);
static void publishLoopInfo (LoopStats *ls);
static void loopStatusCB (void *);
static void cameraLoopPeriodCB(void *);
static void saveFITSHeaderQCB (void *);
static void uptimeCB (void *);
static void sendFITSCB(void *not_used);
static void spcOffloadCB(void *not_used);
static void readNextLocalFITS (FImage *fip);
static void jdFilename (double jd, char *full);
static int takeDark (char whynot[]);
static int saveDark (FImage *fip, char whynot[]);
static void applyDark (FImage *fip);
static int findDark(char whynot[]);
static int takeFlat (char whynot[]);
static int saveFlat (FImage *fip, char whynot[]);
static void applyFlat (FImage *fip);
static int findFlat(char whynot[]);
static int takeBadPixMap (char whynot[]);
static int saveBadPixMap (FImage *fip, char whynot[]);
static void applyBadPixMap (FImage *fip);
static int findBadPixMap(char whynot[]);
static void codeSaturatedPixels (FImage *fip);
static int openFITS (char name[], FImage *fip, char whynot[]);
static void darkFileName(char *path);
static void flatFileName(char *path);
static void badpixFileName(char *path);
static void loopstatsUpdate (double tnow, double pist, double az, double el);
static int calibrateDACs (FILE *fp, double *tip, double *tilt, double *piston);
static int scan4Fringes(void);
static void setTTP (FILE *fp, double tip, double tilt, double piston);
static void setPupil (FILE *fp, double tip, double tilt, double piston);
static int imagePipeline (FImage *fip, double *newtip, double *newtilt, double *newpl);
static void processOneBeam (OneBeam *bp, int beam_x, int beam_y);
static void padFITS (FImage *dst, FImage *src);
static void copy2BG (void);
static void publishImage (FImage *fip, IBLOBVectorProperty *bvp, int eleidx);
static void publishRTData(void);
static void publishHeader (FImage *fip);
static void addFITSNow (FImage *fip, const char *key, const char *comment);
static void forceBeamFit (FImage *fip, int& fits_x, int& fits_y, int& r);
static void forceBothBeamsFit (FImage *fip);
static void readNextImage (FImage *fip);
static void collectRTData (FImage *fip, const char *name);
static int readTO (int fd, char buf[], int maxlen, int msto);
static int readPCam (char buf[], int max, int blocking);
static void pcamReadOneImage (FImage *fip);
static void prepWorking (OneBeam *bp);
static void setupFFTW(OneBeam *bp, int w, int h);
static void setupFITS (FImage *fip, int w, int h);
static void extractBeam (FImage *out, FImage *in, int beam_x, int beam_y);
static int computeCorrections (OneBeam *b1p, OneBeam *b2p, double *tipc, double *tiltc, double *plc);
static double pin (double v, double center);
static void FFTPeakRHP (OneBeam *bp, double *cx, double *cy);
static CamPix imageValue (FImage *fip, int x, int y);
static void contrastMoment (FImage *fip, double *xp);
static void convertFFTW2FITS (fftw_complex *in, FImage *mag, FImage *phs);
static void convertFITS2FFTW (FImage *fip, fftw_complex *out);
static void saveFITS (FImage *fip);
static void saveFITSHeader (FImage *fip);
static void insertIncreasing (double *cp, int n, double newc);
static int quadrantShiftIndex (int i, int sz);
static void quadrantShiftFFTW (fftw_complex *in, int sz);
static int socketMsg (FILE *socket, char buf[], int maxlen, int msto, const char *fmt, ...);
static void cleanStr (char buf[], int l);
static void startPLCThead (void);
static void startOVMSThread(void);
static void setupVibFilterBiQuad(void);
static void accumulateAutoBeamFrames (void);
static void findBeam (FImage *fip, int *minpix, int *maxpix, double *center_x, double *center_y);
static double jdNow (void);
static double secsNow (void);
static double secs2JD (double secs);
static void publishCLShutdown (const char *fmt, ...);
static void timeTrace (bool init, const char *fmt, ...);
static void insureDir (char *dir);
static int mkpath (char *path, int mode);
static void setBeamPad(void);
static void noJitter(void);
static void bye (const char *fmt, ...);

/* send client definitions of all or specified properties */
void
ISGetProperties (char const *dev, char const *nam)
{

        if (dev && strcmp (PLCDEV, dev))
	    return;

	initOnce();

	if (!nam || !strcmp (nam, plcversion.name))
	    IDDefNumber (&plcversion, NULL);
	if (!nam || !strcmp (nam, closeloop.name))
	    IDDefSwitch (&closeloop, NULL);
	if (!nam || !strcmp (nam, spcoffload.name))
	    IDDefSwitch (&spcoffload, NULL);
	if (!nam || !strcmp (nam, useubc.name))
	    IDDefSwitch (&useubc, NULL);
	if (!nam || !strcmp (nam, takedark.name))
	    IDDefSwitch (&takedark, NULL);
	if (!nam || !strcmp (nam, applydark.name))
	    IDDefSwitch (&applydark, NULL);
	if (!nam || !strcmp (nam, takeflat.name))
	    IDDefSwitch (&takeflat, NULL);
	if (!nam || !strcmp (nam, applyflat.name))
	    IDDefSwitch (&applyflat, NULL);
	if (!nam || !strcmp (nam, takebadpix.name))
	    IDDefSwitch (&takebadpix, NULL);
	if (!nam || !strcmp (nam, applybadpix.name))
	    IDDefSwitch (&applybadpix, NULL);
	if (!nam || !strcmp (nam, loopinfo.name))
	    IDDefNumber (&loopinfo, NULL);
	if (!nam || !strcmp (nam, ubcsettings.name))
	    IDDefNumber (&ubcsettings, NULL);
	if (!nam || !strcmp (nam, nacsettings.name))
	    IDDefNumber (&nacsettings, NULL);
	if (!nam || !strcmp (nam, plsetpoint.name))
	    IDDefNumber (&plsetpoint, NULL);
	if (!nam || !strcmp (nam, pcamheader.name))
	    IDDefNumber (&pcamheader, NULL);
	if (!nam || !strcmp (nam, savefits.name))
	    IDDefText (&savefits, NULL);
	if (!nam || !strcmp (nam, savestats.name))
	    IDDefSwitch (&savestats, NULL);
	if (!nam || !strcmp (nam, sendfits.name))
	    IDDefNumber (&sendfits, NULL);
	if (!nam || !strcmp (nam, showbeam1.name))
	    IDDefSwitch (&showbeam1, NULL);
	if (!nam || !strcmp (nam, pcamimage.name))
	    IDDefBLOB (&pcamimage, NULL);
	if (!nam || !strcmp (nam, pcambeamimage.name))
	    IDDefBLOB (&pcambeamimage, NULL);
	if (!nam || !strcmp (nam, pcamfftmimage.name))
	    IDDefBLOB (&pcamfftmimage, NULL);
	if (!nam || !strcmp (nam, pcamfftpimage.name))
	    IDDefBLOB (&pcamfftpimage, NULL);
	if (!nam || !strcmp (nam, pcamdatafield.name))
	    IDDefText (&pcamdatafield, NULL);
	if (!nam || !strcmp (nam, pcamdata.name))
	    IDDefBLOB (&pcamdata, NULL);
	if (!nam || !strcmp (nam, findbeam.name))
	    IDDefNumber (&findbeam, NULL);
	if (!nam || !strcmp (nam, ovmsgains.name))
	    IDDefNumber (&ovmsgains, NULL);
	if (!nam || !strcmp (nam, vibfilter.name))
	    IDDefNumber (&vibfilter, NULL);
	if (!nam || !strcmp (nam, vibgain.name))
	    IDDefNumber (&vibgain, NULL);
	if (!nam || !strcmp (nam, ovmsgains.name))
	    IDDefNumber (&ovmsgains, NULL);
}

void
ISNewNumber (const char *dev, const char *name, double *doubles, char *names[],
int n)
{
        if (!IUCrackNumber (&plsetpoint, dev, name, doubles, names, n)) {

	    /* special property to set setpoint without baggage of full-on nacsettings or ubcsettings.
	     * N.B. setting from here is also mirrored in the them.
	     */

	    double sp = plsetpoint.np[PLSP_PLSP].value;

	    // paranoid
	    if (isnan(sp)) {
		plsetpoint.s = IPS_ALERT;
		IDSetNumber (&plsetpoint, "UBC PLSetpoint was attempted with NaN");
		return;
	    }

	    plsetpoint.s = IPS_OK;
	    if (plsetpoint.np[FORNAC_PLSP].value) {
		if (!USE_UBC) {
		    pthread_mutex_lock (&beam_lock);
		    beam1.sp0 = sp;
		    pthread_mutex_unlock (&beam_lock);
		}
		IDSetNumber (&plsetpoint, "New NAC PLSetpoint %g", sp);
		nacsettings.s = IPS_OK;
		nacsettings.np[PLSP_NACS].value = sp;
		IDSetNumber (&nacsettings, NULL);
	    } else {
		if (USE_UBC) {
		    pthread_mutex_lock (&beam_lock);
		    beam1.sp0 = sp;
		    pthread_mutex_unlock (&beam_lock);
		}
		IDSetNumber (&plsetpoint, "New UBC PLSetpoint %g", sp);
		ubcsettings.s = IPS_OK;
		ubcsettings.np[PLSP_UBCS].value = sp;
		IDSetNumber (&ubcsettings, NULL);
	    }

	}

	else

        if (!IUCrackNumber (&ubcsettings, dev, name, doubles, names, n)) {

	    /* N.B. don't check beams or FFTW here, let imagePipeline do it in case it is in
	     * a separate thread
	     */

	    /* protect BEAM changes from imagePipeline if currently using UBC */
	    if (USE_UBC) {
		pthread_mutex_lock (&beam_lock);
		BEAM1X = (int)(ubcsettings.np[BEAM1X_UBCS].value);
		BEAM1Y = (int)(ubcsettings.np[BEAM1Y_UBCS].value);
		BEAM2X = (int)(ubcsettings.np[BEAM2X_UBCS].value);
		BEAM2Y = (int)(ubcsettings.np[BEAM2Y_UBCS].value);
		BEAMR  = (int)(ubcsettings.np[BEAMR_UBCS].value);
		setBeamPad();
		beam1.sp0 = PLSP;
		pthread_mutex_unlock (&beam_lock);
	    }

	    /* update tuning on the fly */
	    tip_pid.setTuning (TTPGAIN, TTIGAIN, TTDGAIN);
	    tilt_pid.setTuning (TTPGAIN, TTIGAIN, TTDGAIN);
	    piston_pid.setTuning (PLPGAIN, PLIGAIN, PLDGAIN);

	    /* ack */
	    ubcsettings.s = IPS_OK;
	    IDSetNumber (&ubcsettings, NULL);

	    /* mirror to plsetpoint but don't publish to avoid circular loop */
	    plsetpoint.np[PLSP_PLSP].value = ubcsettings.np[PLSP_UBCS].value;
	    plsetpoint.np[FORNAC_PLSP].value = 0;
	}

        else if (!IUCrackNumber (&nacsettings, dev, name, doubles, names, n)) {

	    /* N.B. don't check beams or FFTW here, let imagePipeline do it in case it is in
	     * a separate thread
	     */

	    /* protect BEAM changes from imagePipeline if currently using NAC */
	    if (!USE_UBC) {
		pthread_mutex_lock (&beam_lock);
		BEAM1X = (int)(nacsettings.np[BEAM1X_NACS].value);
		BEAM1Y = (int)(nacsettings.np[BEAM1Y_NACS].value);
		BEAM2X = (int)(nacsettings.np[BEAM2X_NACS].value);
		BEAM2Y = (int)(nacsettings.np[BEAM2Y_NACS].value);
		BEAMR  = (int)(nacsettings.np[BEAMR_NACS].value);
		setBeamPad();
		beam1.sp0 = PLSP;
		pthread_mutex_unlock (&beam_lock);
	    }

	    /* update tuning on the fly */
	    tip_pid.setTuning (TTPGAIN, TTIGAIN, TTDGAIN);
	    tilt_pid.setTuning (TTPGAIN, TTIGAIN, TTDGAIN);
	    piston_pid.setTuning (PLPGAIN, PLIGAIN, PLDGAIN);

	    /* ack */
	    nacsettings.s = IPS_OK;
	    IDSetNumber (&nacsettings, NULL);

	    /* mirror to plsetpoint but don't publish to avoid circular loop */
	    plsetpoint.np[PLSP_PLSP].value = nacsettings.np[PLSP_NACS].value;
	    plsetpoint.np[FORNAC_PLSP].value = 1;
	}

	else if (!IUCrackNumber (&sendfits, dev, name, doubles, names, n)) {

	    /* new sendfits rate */

	    /* stop any existing timer */
	    if (sf_tid) {
		IERmTimer (sf_tid);
		sf_tid = 0;
	    }

            /* always want at least 1 */
            if (SENDFITS < 1)
                SENDFITS = 1;

	    /* start new loop */
            sf_tid = IEAddTimer ((int)(1000.0/SENDFITS), sendFITSCB, NULL);

	    /* ack */
	    sendfits.s = IPS_OK;
            IDSetNumber (&sendfits, "Sending FITS and data at %g Hz", SENDFITS);
	}

	else if (!IUCrackNumber (&findbeam, dev, name, doubles, names, n)) {

	    if (LOOPCLOSED) {

                // control real-time beam positioning
                if (AUTOBEAM_ON) {
                    // init stats
                    autobeam_x_rb->reset();
                    autobeam_y_rb->reset();
                    findbeam.np[BEAMX_PFB].value = -1;
                    findbeam.np[BEAMY_PFB].value = -1;
                    findbeam.s = IPS_BUSY;
                    IDSetNumber (&findbeam, "Auto beam active");
                } else {
                    findbeam.s = IPS_OK;
                    IDSetNumber (&findbeam, "Auto beam off");
                }

	    } else {

                if (AUTOBEAM_ON) {
                    // just find beam one time, adjust FPC and SPC
                    findbeam.s = IPS_BUSY;
                    findbeam.np[BEAMX_PFB].value = -1;
                    findbeam.np[BEAMY_PFB].value = -1;
                    IDSetNumber (&findbeam, "Finding beam centroid");

                    accumulateAutoBeamFrames();
                    double x = autobeam_x_rb->mean();
                    double y = autobeam_y_rb->mean();

                    // TODO: move FPC and SPC

                    // report results
                    findbeam.np[STATE_PFB].value = 0;           // one-shot
                    findbeam.np[BEAMX_PFB].value = x + 1;	// report as FITS
                    findbeam.np[BEAMY_PFB].value = y + 1;
                    findbeam.s = IPS_OK;
                    IDSetNumber (&findbeam, "Beam center at %g %g", x+1, y+1);
                } else {
                    findbeam.s = IPS_OK;
                    IDSetNumber (&findbeam, "Find beam not active");
                }
	    }
	}

	else if (!IUCrackNumber (&ovmsgains, dev, name, doubles, names, n)) {

	    /* sanity checks
	     */
	    if (noovms) {
		ovmsgains.s = IPS_ALERT;
		IDSetNumber (&ovmsgains, "OVMS is not enabled -- see noovms in plc.cfg");
		return;
	    }
	    if (!ovms.ok) {
		ovmsgains.s = IPS_ALERT;
		IDSetNumber (&ovmsgains, "OVMS is not reporting properly");
		return;
	    }

	    /* copy in 
	     */
	    ovms.tip_gain = OVMS_TIP_GAIN_SF*ovmsgains.np[TIP_PLOG].value;
	    ovms.tilt_gain = OVMS_TILT_GAIN_SF*ovmsgains.np[TILT_PLOG].value;
	    ovms.opd_gain = OVMS_OPD_GAIN_SF*ovmsgains.np[OPD_PLOG].value;

	    /* ack
	     */
	    ovmsgains.s = IPS_OK;
	    IDSetNumber (&ovmsgains, "Installed new OVMS gains: tip=%g tilt=%g opd=%g\n", 
		ovms.tip_gain, ovms.tilt_gain, ovms.opd_gain);

	}

	else if (!IUCrackNumber (&vibfilter, dev, name, doubles, names, n)) {

	    /* define filter with new values.
	     * N.B. we assume loop has been closed long enough to get good value for RATE_PLCLI.
	     */
	    setupVibFilterBiQuad();

            /* ack
             */
            vibfilter.s = IPS_OK;
            IDSetNumber (&vibfilter, "Installed new vibration filter: Fc %g, Q %g, Gain %g dB, Gain %g volts",
                        VIB_FC, VIB_Q, VIB_FGAIN, VIB_VGAIN);
	}

	else if (!IUCrackNumber (&vibgain, dev, name, doubles, names, n)) {

	    /* it's sufficient to having updated the property, just ack
	     */
	    vibgain.s = IPS_OK;
	    IDSetNumber (&vibgain, "Vibration piezo gain set to %g", VIB_ZGAIN);
	}

}

void
ISNewText (const char *dev, const char *name, char *texts[], char *names[],
int n)
{
	if (IUCrackText (&savefits, dev, name, texts, names, n) == 0) {

	    /* just ack, updating variable is enough */
	    savefits.s = IPS_OK;
	    IDSetText (&savefits, "Saving every %d images to %s", SAVEFEVERY, SAVEFDIR);
	    insureDir (SAVEFDIR);
	}

	else if (IUCrackText (&pcamdatafield, dev, name, texts, names, n) == 0) {

	    double maxsecs = atoi (pcamdatafield.tp[MAXSECS_PCDF].text);

	    if (maxsecs == 0 || pcamdatafield.tp[FIELD_PCDF].text[0] == '\0') {

		pcamdatafield.s = IPS_OK;
		IDSetText (&pcamdatafield, "Real-time data collection disabled");
		rtdata_dt = 0;

	    } else if (maxsecs < 0) {

		/* restart sample rtdata at the default interval */

		rtdata_dt = rtcolper;
		rtdata_n = 0;
		rtdata_i = 0;

		char buf[32];
		sprintf (buf, "%g", rtdata_dt);
		IUSaveText (&pcamdatafield.tp[MAXSECS_PCDF], buf);

		pcamdatafield.s = IPS_OK;
		IDSetText (&pcamdatafield, "Now recording %s with %.2f seconds history (RTCOLPER default)",
		    pcamdatafield.tp[FIELD_PCDF].text, rtdata_dt);

	    } else {

		/* restart sample rtdata at the specified interval */

		rtdata_dt = maxsecs;
		rtdata_n = 0;
		rtdata_i = 0;

		pcamdatafield.s = IPS_OK;
		IDSetText (&pcamdatafield, "Now recording %s with %.2f seconds history",
		    pcamdatafield.tp[FIELD_PCDF].text, rtdata_dt);
	    }
	}
}

void
ISNewSwitch (const char *dev, const char *name, ISState *states, char *names[], int n)
{
	/* ignore if not ours */
	if (strcmp (dev, PLCDEV))
	    return;

	if (!strcmp (name, closeloop.name) && n == 1) {

	    /* closed-loop control on/off */

            ISwitch *sp = IUFindSwitch (&closeloop, names[0]);
            if (!sp)
                return; /* bogus */
	    bool want_on = states[0] == ISS_ON;

	    if (want_on) {
		// N.B. take care to open pzt_f before LOOPCLOSED reports true
		if (plcCloseLoop() == 0) {
		    sp->s = states[0];
		    closeloop.s = IPS_OK;
		    IDSetSwitch (&closeloop, "Path Length Loop is now running on %s",
				USE_UBC ? "sky" : "NAC");
		}
	    } else {
		// N.B. take care to insure LOOPCLOSED reports false before closing pzt_f
		sp->s = states[0];
		plcOpenLoop();
		closeloop.s = IPS_OK;
		IDSetSwitch (&closeloop, "Path Length Loop is now open");
	    }
	}

	else if (!strcmp (name, takedark.name) && n == 1) {

	    ISwitch *sp = IUFindSwitch (&takedark, names[0]);

	    if (sp == &takedark.sp[TAKE_TD]) {
		sp->s = states[0];
		if (sp->s == ISS_ON) {
		    char msg[1024];

		    if (takeDark(msg) < 0) {
			takedark.s = IPS_ALERT;
			takedark.sp[TAKE_TD].s = ISS_OFF;
			IDSetSwitch (&takedark, "New dark failed: %s", msg);

			DARKUNDEF();
			applydark.sp[APPLY_AD].s = ISS_OFF;
			IDSetSwitch (&applydark, "New dark failed: %s", msg);
		    } else {
			takedark.s = IPS_OK;
			takedark.sp[TAKE_TD].s = ISS_OFF;
			IDSetSwitch (&takedark, "New dark : %s", dark_fn);
		    }

		} else {
		    /* just ack */
		    takedark.s = IPS_OK;
		    IDSetSwitch (&takedark, "No dark taken");
		}
	    }
	}

	else if (!strcmp (name, applydark.name) && n == 1) {

	    ISwitch *sp = IUFindSwitch (&applydark, names[0]);
	    char whynot[1024];

	    if (sp == &applydark.sp[APPLY_AD]) {
		sp->s = states[0];
		if (sp->s == ISS_ON) {
		    if (findDark(whynot) == 0) {
			applydark.s = IPS_OK;
			IDSetSwitch (&applydark, "Darking is On using %s", dark_fn);
		    } else {
			sp->s = ISS_OFF;
			applydark.s = IPS_ALERT;
			IDSetSwitch (&applydark, "No suitable darks found: %s", whynot);
		    }
		} else {
		    applydark.s = IPS_OK;
		    IDSetSwitch (&applydark, "Darking is Off");
		}
	    }
	}

	else if (!strcmp (name, takeflat.name) && n == 1) {

	    ISwitch *sp = IUFindSwitch (&takeflat, names[0]);

	    if (sp == &takeflat.sp[TAKE_TF]) {
		sp->s = states[0];
		if (sp->s == ISS_ON) {
		    char msg[1024];

		    if (takeFlat(msg) < 0) {
			takeflat.s = IPS_ALERT;
			takeflat.sp[TAKE_TF].s = ISS_OFF;
			IDSetSwitch (&takeflat, "New flat failed: %s", msg);

			FLATUNDEF();
			applyflat.sp[APPLY_AF].s = ISS_OFF;
			IDSetSwitch (&applyflat, "New flat failed: %s", msg);
		    } else {
			takeflat.s = IPS_OK;
			takeflat.sp[TAKE_TF].s = ISS_OFF;
			IDSetSwitch (&takeflat, "New flat from %d frames: %s", nflatimg, flat_fn);
		    }

		} else {
		    /* just ack */
		    takeflat.s = IPS_OK;
		    IDSetSwitch (&takeflat, "No flat taken");
		}
	    }
	}

	else if (!strcmp (name, applyflat.name) && n == 1) {

	    ISwitch *sp = IUFindSwitch (&applyflat, names[0]);
	    char whynot[1024];

	    if (sp == &applyflat.sp[APPLY_AF]) {
		sp->s = states[0];
		if (sp->s == ISS_ON) {
		    if (findFlat(whynot) == 0) {
			applyflat.s = IPS_OK;
			IDSetSwitch (&applyflat, "Flating is On using %s", flat_fn);
		    } else {
			sp->s = ISS_OFF;
			applyflat.s = IPS_ALERT;
			IDSetSwitch (&applyflat, "No suitable flats found: %s", whynot);
		    }
		} else {
		    applyflat.s = IPS_OK;
		    IDSetSwitch (&applyflat, "Flating is Off");
		}
	    }
	}

	else if (!strcmp (name, takebadpix.name) && n == 1) {

	    ISwitch *sp = IUFindSwitch (&takebadpix, names[0]);

	    takebadpix.s = IPS_BUSY;
	    IDSetSwitch (&takebadpix, NULL);

	    if (sp == &takebadpix.sp[TAKE_TB]) {
		sp->s = states[0];
		if (sp->s == ISS_ON) {
		    char msg[1024];

		    if (takeBadPixMap(msg) < 0) {
			takebadpix.s = IPS_ALERT;
			takebadpix.sp[TAKE_TB].s = ISS_OFF;
			IDSetSwitch (&takebadpix, "New badpix failed: %s", msg);

			BADPIXUNDEF();
			applybadpix.sp[APPLY_AB].s = ISS_OFF;
			IDSetSwitch (&applybadpix, "New badpix failed: %s", msg);
		    } else {
			int nbadframes;
			getIntFITS (&badpix_fim, "PCBNFRAM", &nbadframes);
			takebadpix.s = IPS_OK;
			takebadpix.sp[TAKE_TB].s = ISS_OFF;
			IDSetSwitch (&takebadpix, "New bad pixel map from %d frames: %s",
						    nbadframes, badpix_fn);
		    }

		} else {
		    /* just ack */
		    takebadpix.s = IPS_OK;
		    IDSetSwitch (&takebadpix, "No badpix taken");
		}
	    }
	}

	else if (!strcmp (name, applybadpix.name) && n == 1) {

	    ISwitch *sp = IUFindSwitch (&applybadpix, names[0]);
	    char whynot[1024];

	    if (sp == &applybadpix.sp[APPLY_AB]) {
		sp->s = states[0];
		if (sp->s == ISS_ON) {
		    if (findBadPixMap(whynot) == 0) {
			applybadpix.s = IPS_OK;
			IDSetSwitch (&applybadpix, "Now scrubbing bad pixels using %s", badpix_fn);
		    } else {
			sp->s = ISS_OFF;
			applybadpix.s = IPS_ALERT;
			IDSetSwitch (&applybadpix, "No suitable badpix map found: %s", whynot);
		    }
		} else {
		    applybadpix.s = IPS_OK;
		    IDSetSwitch (&applybadpix, "Scrubbing bad pixels is now Off");
		}
	    }
	}

	else if (!strcmp (name, useubc.name) && n == 1) {

	    /* enough to just update state */
	    ISwitch *sp = IUFindSwitch (&useubc, names[0]);
	    sp->s = states[0];

	    /* N.B. don't check beams or FFTW here, let imagePipeline do it in case it is in
	     * a separate thread
	     */

	    /* ack */
	    useubc.s = IPS_OK;
	    IDSetSwitch (&useubc, NULL);
	}

	else if (!strcmp (name, savestats.name) && n == 1) {

	    /* enough to just update state */
	    ISwitch *sp = IUFindSwitch (&savestats, names[0]);
	    sp->s = states[0];

	    /* ack */
	    savestats.s = IPS_OK;
	    if (sp->s == ISS_ON)
		IDSetSwitch (&savestats, "Saving FITS stats files to %s", SAVEFDIR);
	    else
		IDSetSwitch (&savestats, "Saving FITS stats files now Off");

	}

	else if (!strcmp (name, showbeam1.name) && n == 1) {

	    /* enough to just update state */
	    ISwitch *sp = IUFindSwitch (&showbeam1, names[0]);
	    sp->s = states[0];

	    /* ack */
	    showbeam1.s = IPS_OK;
	    if (sp->s == ISS_ON)
		IDSetSwitch (&showbeam1, "Now showing beam 1");
	    else
		IDSetSwitch (&showbeam1, "Now showing beam 2");

	}

	else if (!strcmp (name, spcoffload.name) && n == 1) {

	    /* update state */
	    ISwitch *sp = IUFindSwitch (&spcoffload, names[0]);
	    sp->s = states[0];

	    /* make sure any existing timer is off */
	    if (spc_offload_tid) {
		IERmTimer(spc_offload_tid);
		spc_offload_tid = 0;
	    }

	    /* start fresh if want on */
	    spcoffload.s = IPS_OK;
	    if (SPCOFFLOAD) {
		/* want on: ACK and start trimer */
		IDSetSwitch (&spcoffload, "SPC Offloading is now On");
		spcOffloadCB(NULL);
	    } else {
		/* want off: just ACK, timer is already off */
		IDSetSwitch (&spcoffload, "SPC Offloading is now Off");
	    }

	}
}

void
ISNewBLOB (const char *dev, const char *name, int sizes[],
    int blobsizes[], char *blobs[], char *formats[], char *names[], int n)
{
}

/* indiserver is sending us a message from a snooped device */
void
ISSnoopDevice (XMLEle *root)
{
	// if (!IUSnoopNumber (root, &envnow)) {
}

/* call to perform one-time initializations.
 * harmlessly does nothing if called again.
 */
static void
initOnce()
{
	static int inited = 0;

	/* really, just once */
	if (inited)
	    return;
	inited = 1;


	/* lock us into memory */
	if (mlockall(MCL_CURRENT|MCL_FUTURE) < 0)
	    IDLog ("Can not lock into memory: %s\n", strerror(errno));

	/* read config file */
	readConfig();

	/* insure directories exist */
	insureDir (darks_dir);
	insureDir (flats_dir);
	insureDir (badpix_dir);
	insureDir (SAVEFDIR);

	/* init file sharing lock */
	pthread_mutex_init (&plc_bglock, NULL);

	/* init loopstats lock */
	pthread_mutex_init (&loopstats.lock, NULL);

	/* init beam changes lock */
	pthread_mutex_init (&beam_lock, NULL);

	/* init fits header q and lock */
	pthread_mutex_init (&fitshdrq_lock, NULL);
	fitshdr_q = newFQ(100);

	/* start main loop */
	startPLCThead();

	/* start OVMS thread */
	if (!noovms)
	    startOVMSThread();

	/* start beam2 thread */
	startBeam2Thread();

	/* init beam info */
	if (USE_UBC) {
	    BEAM1X = (int)(ubcsettings.np[BEAM1X_UBCS].value);
	    BEAM1Y = (int)(ubcsettings.np[BEAM1Y_UBCS].value);
	    BEAM2X = (int)(ubcsettings.np[BEAM2X_UBCS].value);
	    BEAM2Y = (int)(ubcsettings.np[BEAM2Y_UBCS].value);
	    BEAMR  = (int)(ubcsettings.np[BEAMR_UBCS].value);
	} else {
	    BEAM1X = (int)(nacsettings.np[BEAM1X_NACS].value);
	    BEAM1Y = (int)(nacsettings.np[BEAM1Y_NACS].value);
	    BEAM2X = (int)(nacsettings.np[BEAM2X_NACS].value);
	    BEAM2Y = (int)(nacsettings.np[BEAM2Y_NACS].value);
	    BEAMR  = (int)(nacsettings.np[BEAMR_NACS].value);
	}
	setBeamPad();
	beam1.sp0 = PLSP;

	/* start various background tasks */
	IEAddTimer (plc_period, loopStatusCB, NULL);
	IEAddTimer (5000, uptimeCB, NULL);
	IEAddTimer (saveFHms, saveFITSHeaderQCB, NULL);
	if (camera_period > 0)
	    IEAddTimer (camera_period, cameraLoopPeriodCB, NULL);

	/* start sendfits loop immediately if active */
	if (SENDFITS > 0)
	    sf_tid = IEAddTimer ((int)(1000.0/SENDFITS), sendFITSCB, NULL);

	/* clear client messages */
	closeloop.s = IPS_OK;
	IDSetSwitch (&closeloop, "");
}

/* read config file.
 */
static void
readConfig()
{
	CFValues *cvp = cfLoad (plc_cfgfn, 1, NULL);

	pcam_host = cfMValue (cvp, "pcamimghost", "localhost");
	pcam_port = cfIntValue (cvp, "pcamimgport", 6001);
	pzt_host = cfMValue (cvp, "pzthost", "localhost");
	pzt_port = cfIntValue (cvp, "pztport", 8630);

	ubcsettings.np[BEAM1X_UBCS].value          = cfDblValue (cvp, "UBC_Beam1_x", 100.0);
	ubcsettings.np[BEAM1Y_UBCS].value          = cfDblValue (cvp, "UBC_Beam1_y", 100.0);
	ubcsettings.np[BEAM2X_UBCS].value          = cfDblValue (cvp, "UBC_Beam2_x", 100.0);
	ubcsettings.np[BEAM2Y_UBCS].value          = cfDblValue (cvp, "UBC_Beam2_y", 100.0);
	ubcsettings.np[BEAMR_UBCS].value           = cfDblValue (cvp, "UBC_Beam_r", 0.0);
	ubcsettings.np[MINFSNR_UBCS].value         = cfDblValue (cvp, "UBC_MinFTSNR", 0.0);
	ubcsettings.np[PLSP_UBCS].value            = cfDblValue (cvp, "UBC_PLSetpoint", 0.0);
	ubcsettings.np[PWVGAIN_UBCS].value         = cfDblValue (cvp, "UBC_PLPWVGain", 0.0);
	ubcsettings.np[PLPGAIN_UBCS].value         = cfDblValue (cvp, "UBC_PLPGain", 0.0);
	ubcsettings.np[PLDGAIN_UBCS].value         = cfDblValue (cvp, "UBC_PLDGain", 0.0);
	ubcsettings.np[PLIGAIN_UBCS].value         = cfDblValue (cvp, "UBC_PLIGain", 0.0);
	ubcsettings.np[CGSP_UBCS].value            = cfDblValue (cvp, "UBC_CGSP", 0.0);
	ubcsettings.np[CGSC_UBCS].value            = cfDblValue (cvp, "UBC_CGScale", 0.0);
	ubcsettings.np[TIPSP_UBCS].value           = cfDblValue (cvp, "UBC_TipSetpoint", 0.0);
	ubcsettings.np[TILTSP_UBCS].value          = cfDblValue (cvp, "UBC_TiltSetpoint", 0.0);
	ubcsettings.np[TTPGAIN_UBCS].value         = cfDblValue (cvp, "UBC_TTPGain", 0.0);
	ubcsettings.np[TTDGAIN_UBCS].value         = cfDblValue (cvp, "UBC_TTDGain", 0.0);
	ubcsettings.np[TTIGAIN_UBCS].value         = cfDblValue (cvp, "UBC_TTIGain", 0.0);

	nacsettings.np[BEAM1X_NACS].value          = cfDblValue (cvp, "NAC_Beam1_x", 100.0);
	nacsettings.np[BEAM1Y_NACS].value          = cfDblValue (cvp, "NAC_Beam1_y", 100.0);
	nacsettings.np[BEAM2X_NACS].value          = cfDblValue (cvp, "NAC_Beam2_x", 100.0);
	nacsettings.np[BEAM2Y_NACS].value          = cfDblValue (cvp, "NAC_Beam2_y", 100.0);
	nacsettings.np[BEAMR_NACS].value           = cfDblValue (cvp, "NAC_Beam_r", 0.0);
	nacsettings.np[MINFSNR_NACS].value         = cfDblValue (cvp, "NAC_MinFTSNR", 0.0);
	nacsettings.np[PLSP_NACS].value            = cfDblValue (cvp, "NAC_PLSetpoint", 0.0);
	nacsettings.np[PWVGAIN_NACS].value         = cfDblValue (cvp, "NAC_PLPWVGain", 0.0);
	nacsettings.np[PLPGAIN_NACS].value         = cfDblValue (cvp, "NAC_PLPGain", 0.0);
	nacsettings.np[PLDGAIN_NACS].value         = cfDblValue (cvp, "NAC_PLDGain", 0.0);
	nacsettings.np[PLIGAIN_NACS].value         = cfDblValue (cvp, "NAC_PLIGain", 0.0);
	nacsettings.np[CGSP_UBCS].value            = cfDblValue (cvp, "NAC_CGSP", 0.0);
	nacsettings.np[CGSC_UBCS].value            = cfDblValue (cvp, "NAC_CGScale", 0.0);
	nacsettings.np[TIPSP_NACS].value           = cfDblValue (cvp, "NAC_TipSetpoint", 0.0);
	nacsettings.np[TILTSP_NACS].value          = cfDblValue (cvp, "NAC_TiltSetpoint", 0.0);
	nacsettings.np[TTPGAIN_NACS].value         = cfDblValue (cvp, "NAC_TTPGain", 0.0);
	nacsettings.np[TTDGAIN_NACS].value         = cfDblValue (cvp, "NAC_TTDGain", 0.0);
	nacsettings.np[TTIGAIN_NACS].value         = cfDblValue (cvp, "NAC_TTIGain", 0.0);

	max_pl_corr = cfDblValue (cvp, "MAX_PL_CORR", 30.0);
	max_tip_corr = cfDblValue (cvp, "MAX_TIP_CORR", 30.0);
	max_tip_err = cfDblValue (cvp, "MAX_TIP_ERR", 30.0);
	max_tilt_corr = cfDblValue (cvp, "MAX_TILT_CORR", 30.0);
	max_tilt_err = cfDblValue (cvp, "MAX_TILT_ERR", 30.0);
	sat_pix = cfIntValue (cvp, "SAT_PIX", 65535);
	bad_space_std = cfDblValue (cvp, "BAD_SPACE_PIX_STD", 5.0);
	bad_time_std = cfDblValue (cvp, "BAD_TIME_PIX_STD", 5.0);
	cg_min = cfDblValue (cvp, "CG_MIN", 0.3);
	cg_freeze = cfDblValue (cvp, "CG_FREEZE", 10.0)/1000.0;	// ms to secs
	cg_boxcar_period = cfDblValue (cvp, "CG_BOXCAR_PERIOD", 1.0);
	spc_offload_ms = cfDblValue (cvp, "SPC_OFFLOAD_TIME", 30.0) * 1000.0;
	spc_offload_scale = cfDblValue (cvp, "SPC_OFFLOAD_SCALE", 580.0);
	spc_offload_alpha = cfDblValue (cvp, "SPC_OFFLOAD_ALPHA", 0.1);
	spc_offload_max = cfDblValue (cvp, "SPC_OFFLOAD_MAX", 50.0);
	spc_offload_devnam = cfMValue (cvp, "SPC_OFFLOAD_MOTNAM", "Ubcs.SPC_Trans");
	char *dot = strchr (spc_offload_devnam, '.');
	if (!dot)
	    bye ("SPC_OFFLOAD_MOTNAM format must be device.name: %s", spc_offload_devnam);
	spc_offload_motnam = dot+1;
	*dot = '\0';

	IUSaveText (&savefits.tp[EVERY_SVF], cfValue (cvp, "savefits_every", "0"));
	IUSaveText (&savefits.tp[DIR_SVF], cfValue (cvp, "savefits_dir", "/tmp"));
	darks_dir = cfMValue (cvp, "darks_dir", "/tmp");
	(void) mkdir (darks_dir, 0755);
	flats_dir = cfMValue (cvp, "flats_dir", "/tmp");
	(void) mkdir (flats_dir, 0755);
	nflatimg = cfIntValue (cvp, "NFLATIMGS", 10);
	badpix_dir = cfMValue (cvp, "badpix_dir", "/tmp");
	(void) mkdir (badpix_dir, 0755);
	badmapsecs = cfIntValue (cvp, "BADMAPSECS", 5);
	ramp_time = cfDblValue (cvp, "RAMP_TIME", 10.0);
	b2sp0_time = cfDblValue (cvp, "B2SP0_TIME", 1.0);
	wv_boxcar_period = cfDblValue (cvp, "WV_BOXCAR_PERIOD", 1.0);

	plc_period = cfDblValue (cvp, "LOOP_STATS_PERIOD", 27.0);
	SENDFITS = cfDblValue (cvp, "sendfits", 0.0);
	trace = cfIntValue (cvp, "trace", 0);
	nopchw = cfIntValue (cvp, "nopchw", 1);
	nofpchw = cfIntValue (cvp, "nofpchw", 1);
	noovms = cfIntValue (cvp, "noovms", 1);
	plpfc = cfDblValue (cvp, "PLPFC", 20.0);

	beam_pad0 = cfIntValue (cvp, "BEAM_PAD", 1);

	rtcolper = cfDblValue (cvp, "RTCOLPER", 10.0);
	rtpubper = cfDblValue (cvp, "RTPUBPER", 0.5);

        ovmsgains.np[TIP_PLOG].value = ovms.tip_gain = cfDblValue (cvp, "OVMS_TIP_GAIN", 0);
        ovmsgains.np[TILT_PLOG].value = ovms.tilt_gain = cfDblValue (cvp, "OVMS_TILT_GAIN", 0);
        ovmsgains.np[OPD_PLOG].value = ovms.opd_gain = cfDblValue (cvp, "OVMS_OPD_GAIN", 0);
        OVMS_TIP_GAIN_SF = cfDblValue (cvp, "OVMS_TIP_GAIN_SF", 1);
        OVMS_TILT_GAIN_SF = cfDblValue (cvp, "OVMS_TILT_GAIN_SF", 1);
        OVMS_OPD_GAIN_SF = cfDblValue (cvp, "OVMS_OPD_GAIN_SF", 1);
        ovms.port = cfIntValue (cvp, "OVMS_PORT", 1);
        ovms.group = cfMValue (cvp, "OVMS_GROUP", "undefined");
        ovms.interface = cfMValue (cvp, "OVMS_NET_IF", "undefined");
	opd_bump = cfDblValue (cvp, "OVMS_OPD_BUMP", 1000.0);
	opd_freeze = cfDblValue (cvp, "OVMS_OPD_FREEZE", 0.0);

	camera_host = cfMValue (cvp, "CAMERA_HOST", "");
	camera_port = cfIntValue (cvp, "CAMERA_PORT", 11123);
	camera_period = cfIntValue (cvp, "CAMERA_POLL", 1000);

	VIB_FC = cfDblValue (cvp, "VIB_FILTER_FC", 12.0);
	VIB_Q = cfDblValue (cvp, "VIB_FILTER_Q", 1.0);
	VIB_FGAIN = cfDblValue (cvp, "VIB_FILTER_GAIN", 50.0);
	VIB_ZGAIN = cfDblValue (cvp, "UBC_VIB_Gain", 0.0);

        nautobeam = cfIntValue (cvp, "N_AUTOBEAM_FRAMES", 100);
        autobeam_x_rb = new RingBuffer (nautobeam);
        autobeam_y_rb = new RingBuffer (nautobeam);

	beam1.lambda = lam1_PC = cfDblValue (cvp, "LAMBDA1_PC", 2.2);
	beam2.lambda = lam2_PC = cfDblValue (cvp, "LAMBDA2_PC", 2.2);
	lam_NOMIC = cfDblValue (cvp, "LAMBDA_NOMIC", 11.0);

	cpu_plcloop = cfIntValue (cvp, "CPU_PLCLOOP", 0);
	cpu_beam2 = cfIntValue (cvp, "CPU_BEAM2", 1);

	peakignore = cfDblValue (cvp, "PEAKIGNORE", 0.5);

	fscan_min = cfDblValue (cvp, "FSCAN_MIN", 0.0);
	fscan_max = cfDblValue (cvp, "FSCAN_MAX", 0.0);
	fscan_step = cfDblValue (cvp, "FSCAN_STEP", 0.0);
	fscan_nsweeps = cfIntValue (cvp, "FSCAN_NSWEEPS", 0);
	fscan_settle = cfIntValue (cvp, "FSCAN_SETTLE", 0) * 1000;	// ms to us

	UBC_PLPGain_SF = cfDblValue (cvp, "UBC_PLPGain_SF", 1.0);
	UBC_PLIGain_SF = cfDblValue (cvp, "UBC_PLIGain_SF", 1.0);
	UBC_PLDGain_SF = cfDblValue (cvp, "UBC_PLDGain_SF", 1.0);
	UBC_TTPGain_SF = cfDblValue (cvp, "UBC_TTPGain_SF", 1.0);
	UBC_TTIGain_SF = cfDblValue (cvp, "UBC_TTIGain_SF", 1.0);
	UBC_TTDGain_SF = cfDblValue (cvp, "UBC_TTDGain_SF", 1.0);
	NAC_PLPGain_SF = cfDblValue (cvp, "NAC_PLPGain_SF", 1.0);
	NAC_PLDGain_SF = cfDblValue (cvp, "NAC_PLDGain_SF", 1.0);
	NAC_PLIGain_SF = cfDblValue (cvp, "NAC_PLIGain_SF", 1.0);
	NAC_TTPGain_SF = cfDblValue (cvp, "NAC_TTPGain_SF", 1.0);
	NAC_TTDGain_SF = cfDblValue (cvp, "NAC_TTDGain_SF", 1.0);
	NAC_TTIGain_SF = cfDblValue (cvp, "NAC_TTIGain_SF", 1.0);

	cfFree (cvp);

}

/* disable the Nagle algorithm on the given TCP socket.
 * this has the effect of sending data immediately without attempting to collect small packets
 * together to reduce overhead.
 */
static void
noNagle (int socket)
{
        int on = 1;
        if (setsockopt (socket, IPPROTO_TCP, TCP_NODELAY, &on, sizeof(on)) < 0)
            bye ("Can not disable nagle algorithm\n");
}

/* return a socket descriptor connected to the given host/port.
 * if trouble return -1 with brief reason in whynot[].
 */
static int
openSocket (const char *host, int port, bool log_success, char whynot[])
{
	struct sockaddr_in serv_addr;
	struct hostent *hp;
	int sockfd;

	if (trace)
	    IDLog ("Opening %s:%d\n", host, port);

	/* lookup host address */
	hp = gethostbyname (host);
	if (!hp) {
	    sprintf (whynot, "gethostbyname(%s): %s", host, strerror(errno));
	    return (-1);
	}

	/* create a socket */
	(void) memset ((char *)&serv_addr, 0, sizeof(serv_addr));
	serv_addr.sin_family = AF_INET;
	serv_addr.sin_addr.s_addr = ((struct in_addr *)(hp->h_addr_list[0]))->s_addr;
	serv_addr.sin_port = htons(port);
	if ((sockfd = socket (AF_INET, SOCK_STREAM, 0)) < 0) {
	    sprintf (whynot, "socket(%s,%d): %s", host, port, strerror(errno));
	    return (-1);
	}

	/* connect */
	if (connect (sockfd,(struct sockaddr *)&serv_addr,sizeof(serv_addr)) < 0) {
	    sprintf (whynot, "connect(%s,%d): %s", host, port, strerror(errno));
	    return (-1);
	}

	/* done */
	if (log_success)
	    IDLog ("%s:%d opened successfully on fd %d\n", host, port, sockfd);
	return (sockfd);
}


/* start the "real-time" thread.
 * exit if real trouble
 */
static void
startPLCThead ()
{
	pthread_attr_t attr;
	pthread_t plc_t;
	int pret;

	/* create new detached thread which starts loop immediately */
	pthread_attr_init (&attr);
	pthread_attr_setdetachstate (&attr, PTHREAD_CREATE_DETACHED);

	/* crank up priority */
	if (pthread_attr_setschedpolicy (&attr, SCHED_RR))
	    IDLog ("INFO: Can not set scheduling policy for plcLoop to SCHED_RR\n");
	else {
	    IDLog ("INFO: Successfully set scheduling policy for plcLoop to SCHED_RR\n");
	    struct sched_param shp;
	    int max = sched_get_priority_max(SCHED_RR);
	    shp.sched_priority = max;
	    if (pthread_attr_setschedparam (&attr, &shp) != 0)
		IDLog ("INFO: Can not set scheduling priority for plcLoop to %d\n", max);
	    else
		IDLog ("INFO: Successfully set scheduling priority for plcLoop to %d\n", max);
	}

	pret = pthread_create (&plc_t, &attr, plcLoopThread, NULL);
	if (pret) {
	    /* error creating thread */
	    closeloop.sp[YES_CLP].s = ISS_OFF;
	    closeloop.s = IPS_ALERT;
	    IDSetSwitch (&closeloop, "Error creating thread: %s", strerror(pret));
	    bye ("pthread_create(): %s\n", strerror(pret));
	}

	/* stay on cpu_plcloop */
	cpu_set_t cpuset;
	CPU_ZERO (&cpuset);
	CPU_SET (cpu_plcloop, &cpuset);
	if (pthread_setaffinity_np(plc_t, sizeof(cpu_set_t), &cpuset))
	    IDLog ("INFO: Can not lock plcloop thread to cpu %d\n", cpu_plcloop);
	else
	    IDLog ("INFO: Successfully locked plcloop thread to cpu %d\n", cpu_plcloop);
}

/* this function is the separate "closed loop" thread to run the PLC algorithm.
 * we run forever unless some fatal error.
 */
static void *
plcLoopThread (void *not_used)
{
	/* attempt to set us to high priority */
	if (setpriority (PRIO_PROCESS, 0, -19) < 0)
	    IDLog ("Warning: Can not set high process priority: %s\n", strerror(errno));

        /* ok, loop forever */
        while (1) {
	    plcOneIteration();
	    loopstats.plc_loops++;
	}

	/* lint */
	return (NULL);
}


/* call this function when "closing" the loop.
 * return 0 if ok, else -1 with closeloop already reported
 */
static int
plcCloseLoop ()
{
	/* init loopstats */
	pthread_mutex_lock (&loopstats.lock);
	memset (&loopstats, 0, sizeof(loopstats));
	pthread_mutex_unlock (&loopstats.lock);

	/* open acromag command socket */
	if (!nofpchw) {
	    if (pzt_f)
		fclose (pzt_f);
	    char ynot[1024];
	    int fd = openSocket (pzt_host, pzt_port, true, ynot);
	    if (fd < 0)
		bye ("Error opening ACROMAG server %s:%d: %s", pzt_host, pzt_port, ynot);
            noNagle (fd);
	    pzt_f = fdopen (fd, "r+");
	    if (!pzt_f)
		bye ("Can not create r+ FILE from ACROMAG server socket %d", fd);
	    setbuf (pzt_f, NULL);
	}

	/* read initial DAC values */
	if (!nofpchw) {
	    if (calibrateDACs (pzt_f, &Tip_Now, &Tilt_Now, &Piston_Now) < 0)
		return (-1);	// already published reason

	    // save and restore after loop is opened
	    pzt_tip0 = Tip_Now;
	    pzt_tilt0 = Tilt_Now;
	    pzt_pist0 = Piston_Now;
	}

	/* move to fringes or abort if can't find */
	if (!nofpchw && MINFSNR >= 0 && fscan_nsweeps > 0 && scan4Fringes() < 0)
	    return (-1);		// plcOpenLoop() already restored starting pos and closed pzt_f
	
	/* flag to tell computeCorrections() to restart */
	init_plloop = true;

	/* init beam setpoints */
	beam1.sp0 = beam1.setpoint = PLSP;
	beam2.sp0 = beam2.setpoint = 0;

	/* restart unwrapped phase stats */
	beam1.initUnwrappedPhase();
	beam2.initUnwrappedPhase();

	return (0);
}


/* call this function to "open" the loop.
 */
static void
plcOpenLoop ()
{

	/* restore starting position */
	if (pzt_f)
	    setTTP (pzt_f, pzt_tip0, pzt_tilt0, pzt_pist0);

	/* reset and publish null stats */
	pthread_mutex_lock (&loopstats.lock);
	memset (&loopstats, 0, sizeof(loopstats));
	publishLoopInfo (&loopstats);
	pthread_mutex_unlock (&loopstats.lock);

	/* reset GUI setpoint to current beam1 value in range += 180 */
	double b1sp = beam1.setpoint;
	while (b1sp > 180)
	    b1sp -= 360;
	while (b1sp < -180)
	    b1sp += 360;
	if (USE_UBC) {
	    ubcsettings.np[PLSP_UBCS].value = b1sp;
	    IDSetNumber (&ubcsettings, NULL);
	} else {
	    nacsettings.np[PLSP_NACS].value = b1sp;
	    IDSetNumber (&nacsettings, NULL);
	}

        /* close pzt connection */
	if (pzt_f) {
	    fclose (pzt_f);
	    pzt_f = NULL;
	}
}

/* read current values of tip, tilt and piston from the channels that correspond to the
 * current USE_UBC setting
 * return 0 if ok, else -1 with closeloop already reported
 */
static int
calibrateDACs(FILE *fp, double *tip, double *tilt, double *piston)
{
	char buf[1024];

	if (USE_UBC) {
	    if (socketMsg (fp, buf, sizeof(buf), 2000, "read_fpc_tip_tilt_piston\n") < 0) {
		publishCLShutdown ("Error sending read_fpc_tip_tilt_piston: %s\n", buf);
		return (-1);
	    }
	    if (sscanf (buf, "%lf %lf %lf", tip, tilt, piston) != 3) {
		publishCLShutdown ("Error response from read_fpc_tip_tilt_piston: %s\n", buf);
		return (-1);
	    }
	    IDLog ("Inital FPC Tip = %g asec, Tilt = %g asec, Piston = %g um\n", *tip, *tilt, *piston);

	} else {
	    if (socketMsg (fp, buf, sizeof(buf), 2000, "read_field_tip_tilt_piston\n") < 0) {
		publishCLShutdown ("Error sending read_field_tip_tilt_piston: %s\n", buf);
		return (-1);
	    }
	    if (sscanf (buf, "%lf %lf %lf", tip, tilt, piston) != 3) {
		publishCLShutdown ("Error response from read_field_tip_tilt_piston: %s\n", buf);
		return (-1);
	    }
	    IDLog ("Inital Field Tip = %g asec, Tilt = %g asec, Piston = %g um\n", *tip, *tilt, *piston);
	}

	return (0);
}

/* atomically write the given TTP values to the appropriate piezo family, depending on USE_UBC.
 * tip and tilt in arc secs, piston in microns.
 * N.B. must be done by same thread as setPupil
 * TODO: known issue: if alternate here and setPupil socket err msgs can be delayed enough to get reported 
 *   by subsequenct call
 */
static void
setTTP (FILE *fp, double tip, double tilt, double piston)
{
	char buf[1024];
	int nr;

	if (USE_UBC) {
	    if (socketMsg (fp, buf, sizeof(buf), 0, "write_fpc_tip_tilt_piston %g %g %g\n",
			tip, tilt, piston) < 0) {
		publishCLShutdown ("write_fpc_tip_tilt_piston(%g,%g,%g): %s", tip, tilt, piston, buf);
		return;
	    }
	    if ((nr = readTO (fileno(fp), buf, sizeof(buf), 0)) > 0) {
		publishCLShutdown ("write_fpc_tip_tilt_piston(%g,%g,%g): %.*s", tip, tilt, piston, nr, buf);
		return;
	    }

	} else {

	    if (socketMsg (fp, buf, sizeof(buf), 0, "write_field_tip_tilt_piston %g %g %g\n",
			tip, tilt, piston) < 0) {
		publishCLShutdown ("write_field_tip_tilt_piston(%g,%g,%g): %s", tip, tilt, piston, buf);
		return;
	    }
	    if ((nr = readTO (fileno(fp), buf, sizeof(buf), 0)) > 0) {
		publishCLShutdown ("write_field_tip_tilt_piston(%g,%g,%g): %.*s", tip, tilt, piston, nr, buf);
		return;
	    }
	}
}

/* scan the PZTs in piston end-to-end, searching for the location which creates the largest snr.
 * if this snr is at least MINFSNR, move the PZTs to that location and return 0.
 * if not found publish why with publishCLShutdown() which calls plcOpenLoop() which restores
 *   starting positions and return -1.
 * N.B. we assume !nofpchw
 * N.B. we use pzt_tip0 and pzt_tilt0 while scanning
 */
static int
scan4Fringes()
{
	double max_snr = 0;				// max SNR found so far
	double max_at;					// piston @ max_snr
	int seqn = 0, prev_seqn;			// insure a fresh image each time
	double scan_tip = pzt_tip0;			// scanning tip 
	double scan_tilt = pzt_tilt0;			// scanning tip 
	double snr;					// current SNR

	IDSetSwitch (&closeloop, "Begin scanning for fringes ...");

	// perform all sweeps
	for (int i = 0; i < fscan_nsweeps; i++) {

	    // perform one sweep
	    for (double piston = fscan_min; piston <= fscan_max; piston += fscan_step) {

		// move PZT
		setTTP (pzt_f, scan_tip, scan_tilt, piston);	// closes pzt_f if any acromag problems
		if (!pzt_f)
		    return (-1);
		usleep (fscan_settle);

		// get snr from next _unique_ frame being collected by the imaging threads
		do {
		    // wait for stable bg_cam_fim copy
		    pthread_mutex_lock (&plc_bglock);

		    // get this frame's sequence number and SNR
		    prev_seqn = seqn;
		    if (getIntFITS (&bg_cam_fim, "PCSEQN", &seqn) < 0) {
			publishCLShutdown ("Can not find PCSEQN");
			return (-1);
		    }
		    if (getRealFITS (&bg_cam_fim, "PCMFSNR", &snr) < 0) {
			publishCLShutdown ("Can not find PCMFSNR");
			return (-1);
		    }

		    // done with bg_cam_fim, ok to let real time clobber it again
		    pthread_mutex_unlock (&plc_bglock);


		} while (prev_seqn == seqn);

		IDLog ("SNR Scan: piston %g SNR %g\n", piston, snr);

		// record position if this is higher
		if (snr > max_snr) {
		    max_snr = snr;
		    max_at = piston;
		}
	    }
	}


	// decide
	if (max_snr >= MINFSNR) {

	    // found a qualifying position: move there coming at the same direction to avoid hystersis
	    setTTP (pzt_f, scan_tip, scan_tilt, fscan_min);
	    usleep (2*fscan_settle);
	    setTTP (pzt_f, scan_tip, scan_tilt, max_at);
	    usleep (fscan_settle);

	    // record as new current
	    Tip_Now = scan_tip;
	    Tilt_Now = scan_tilt;
	    Piston_Now = max_at;

	    IDSetSwitch (&closeloop, "Scanning successful: best SNR %g at %g microns from reference",
			    max_snr, max_at);
	    return (0);

	} else {
	    // search failed
	    publishCLShutdown ("Fringe scanning failed: best SNR found was %g but min is %g\n",
	    	max_snr, MINFSNR);
	    return (-1);
	}
}

/* issue the given absolute values to the Pupil mirror
 * return 0 if ok, else -1 with reason already publshed via publishCLShutdown
 * N.B. must be done by same thread as setTTP
 * TODO: known issue: if alternate here and setTTP socket err msgs can be delayed enough to get reported 
 *   by subsequenct call
 */
static void
setPupil (FILE *fp, double tip, double tilt, double piston)
{
	char buf[1024];
	int nr;

	if (socketMsg (fp, buf, sizeof(buf), 0, "write_pupil_tip_tilt_piston %g %g %g\n",
		    tip, tilt, piston) < 0) {
	    publishCLShutdown ("write_pupil_tip_tilt_piston(%g,%g,%g): %s", tip, tilt, piston, buf);
	    return;
	}
	if ((nr = readTO (fileno(fp), buf, sizeof(buf), 0)) > 0) {
	    publishCLShutdown ("write_pupil_tip_tilt_piston(%g,%g,%g): %.*s", tip, tilt, piston, nr, buf);
	    return;
	}
}

/* one iteration of the main PLC loop
 * if trouble open the loop,
 */
static void
plcOneIteration()
{
	FImage rt_fim;
	double cmd_tip, cmd_tilt, cmd_pl;

	noJitter();
	timeTrace (true, "start PLC iteration\n");

	/* clean stack garbage */
	initFImage (&rt_fim);

	/* read next image */
	readNextImage (&rt_fim);
	noJitter();

	/* calculate stats and corrections, install if closed */
	if (imagePipeline (&rt_fim, &cmd_tip, &cmd_tilt, &cmd_pl) == 0 && !nofpchw) {

	    if (trace)
		IDLog ("Commanding acromag\n");

	    /* add this to log precise moment when command to send to acromag
	    { struct timeval tv;
	      gettimeofday (&tv, NULL);
	      IDLog ("Sending Acromag command at %d.%06d\n", tv.tv_sec, tv.tv_usec);
	    }
	    */

	    // engage!
	    if (pzt_f && LOOPCLOSED) {
		// update global record
		Tip_Now = cmd_tip;
		Tilt_Now = cmd_tilt;
		Piston_Now = cmd_pl;
		Piston_Mean = spc_offload_alpha*Piston_Now + (1-spc_offload_alpha)*Piston_Mean;

		// go
		setTTP (pzt_f, Tip_Now, Tilt_Now, Piston_Now);
	    }

	}
	noJitter();

	/* copy images and rtdata buffers for INDI "background" thread but don't wait if busy */
	if (!pthread_mutex_trylock (&plc_bglock)) {

	    /* copy selected images to b1_*_fim */
	    copy2BG();

	    /* add desired field, if set, to rtdata ring buffers */
	    collectRTData (&rt_fim, pcamdatafield.tp[FIELD_PCDF].text);

	    /* unlock */
	    pthread_mutex_unlock (&plc_bglock);
	}

	/* save image and/or stats if desired */
	if (SAVEFEVERY > 0 && (loopstats.plc_loops % SAVEFEVERY) == 0)
	    saveFITS (&rt_fim);
	if (SAVESTATS)
	    saveFITSHeader (&rt_fim);

	/* finished with rt_fim */
	resetFImage (&rt_fim);
}

/* add pcamdatafield from the given FITS header to rtdata ring buffers
 */
static void
collectRTData (FImage *fip, const char *name)
{

	// skip if no name or sample period
	if (!name || name[0] == '\0' || rtdata_dt <= 0)
	    return;

	double jdvalue, yvalue;
	bool ok = true;

	// make sure we can get the desired field
	if (getRealFITS (fip, name, &yvalue) < 0) {
	    int ivalue;
	    if (getIntFITS (fip, name, &ivalue) < 0)
		ok = false;
	    else
		yvalue = ivalue;
	}
	if (!ok) {
	    pcamdatafield.s = IPS_ALERT;
	    IDSetText (&pcamdatafield, "Field %s not found", pcamdatafield.tp[FIELD_PCDF].text);
	    rtdata_dt = 0;
	    return;
	}

	// and its time stamp
	if (ok && getRealFITS (fip, jd_kw, &jdvalue) < 0)
	    ok = false;
	if (!ok) {
	    pcamdatafield.s = IPS_ALERT;
	    IDSetText (&pcamdatafield, "Field %s not found", jd_kw);
	    rtdata_dt = 0;
	    return;
	}

	// add to rtdata if fields found ok
	if (ok) {
	    rtdata_t[rtdata_i] = jdvalue;
	    rtdata_y[rtdata_i] = yvalue;
	    if (rtdata_n < MAXRTSAMPS)
		rtdata_n++;
	    if (++rtdata_i >= MAXRTSAMPS)
		rtdata_i = 0;
	}
}


/* given a raw phasecam image, break into two beams and process each to find new commanded
 * corrections. return 0 if corrections are valid to be applied, else -1.
 */
static int
imagePipeline (FImage *fip, double *cmd_tip, double *cmd_tilt, double *cmd_pl)
{
	// apply bad pixel map if switched on
	if (APPLYBADPIX)
	    applyBadPixMap (fip);

	// apply dark if switched on
	if (APPLYDARK)
	    applyDark (fip);

	// apply flat if switched on
	if (APPLYFLAT)
	    applyFlat (fip);

	// insure unchanging beam parameters
	pthread_mutex_lock (&beam_lock);

	// capture fip as raw image for both beams
	beam1.raw_fimp = fip;
	beam2.raw_fimp = fip;

	// insure working files and FFT are proper size for beam
	prepWorking(&beam1);
	prepWorking(&beam2);

	// signal beam2 thread to start processing
	pthread_mutex_lock (&b2_run_lock);
	b2_run = true;
	pthread_cond_signal (&b2_run_cv);
	pthread_mutex_unlock (&b2_run_lock);

	// meanwhile we process beam1 here
	processOneBeam (&beam1, BEAM1X, BEAM1Y);

	// wait for beam2 to finish via cv
	pthread_mutex_lock (&b2_done_lock);
	while (!b2_done)
	    pthread_cond_wait (&b2_done_cv, &b2_done_lock);
	b2_done = false;
	pthread_mutex_unlock (&b2_done_lock);


	// update beam1 unwrapped phase by itself
	beam1.updateUnwrappedPhase();

	// update beam2 unwrapped phase pinned around current scaled beam1 offset
	beam2.phase_unw += beam2.phase - beam2.last_phase;
	beam2.last_phase = beam2.phase;

	// pin unwrapped beam2 phase to scaled beam1 phase
	double center = beam1.phase_unw + beam2.setpoint
		+ (beam1.phase_unw - beam1.setpoint) * ((beam1.lambda - beam2.lambda)/beam2.lambda);
	beam2.phase_unw = pin (beam2.phase_unw, center);


	// compute corrections
	int cret = computeCorrections (&beam1, &beam2, cmd_tip, cmd_tilt, cmd_pl);

	// ok to allow beam parameters to change again
	pthread_mutex_unlock (&beam_lock);

	return (cret);
}

/* this thread runs forever and is used to compute beam2 info.
 */
static void *
beam2Thread (void *dummy)
{
	/* attempt to set us to high priority */
	if (setpriority (PRIO_PROCESS, 0, -19) < 0)
	    IDLog ("Warning: Can not set high process priority: %s\n", strerror(errno));

	/* forever */
	while (true) {

	    // wait to be told to start
	    pthread_mutex_lock (&b2_run_lock);
	    while (!b2_run)
		pthread_cond_wait (&b2_run_cv, &b2_run_lock);
	    b2_run = false;
	    pthread_mutex_unlock (&b2_run_lock);

	    // process beam2
	    processOneBeam (&beam2, BEAM2X, BEAM2Y);

	    // signal we are finished
	    pthread_mutex_lock (&b2_done_lock);
	    b2_done = true;
	    pthread_cond_signal (&b2_done_cv);
	    pthread_mutex_unlock (&b2_done_lock);
	}

	bye ("beam2Thread ended unexpectedly\n");
	return (NULL);
}

/* given one beam descriptor and the beam location on the raw image,
 * extract the beam from raw_fimp, pad it, perform the FFTs and derive the measured phase, tip and tilt.
 */
static void
processOneBeam (OneBeam *bp, int beam_x, int beam_y)
{
	// form padded image
	FImage justbeam_fim;
	initFImage (&justbeam_fim);
	setupFITS (&justbeam_fim, BEAMW, BEAMW);
	extractBeam (&justbeam_fim, bp->raw_fimp, beam_x, beam_y);
	padFITS (&bp->pad_fim, &justbeam_fim);
	resetFImage (&justbeam_fim);

	// prep for FFT
	convertFITS2FFTW (&bp->pad_fim, bp->fft_input);
	quadrantShiftFFTW (bp->fft_input, bp->pad_fim.sw);

	// compute the 2D FFT from fft_input to fft_output
	fftw_execute (bp->fft_theplan);

	// convert FFT to mag and phase images
	convertFFTW2FITS (bp->fft_output, &bp->fftmag_fim, &bp->fftphs_fim);

	// find centroided peak of region on right half of mag image to find tip and tilt
	double cx, cy;
	FFTPeakRHP (bp, &cx, &cy);

	// find value on phase image at same location as tip/tilt mag peak 
	bp->phs_pixel = imageValue (&bp->fftphs_fim, (int)(cx+0.5), (int)(cy+0.5));

	// map unsigned pixel 0..65k to phase -180..+180 degrees
	bp->phase = 360.0*bp->phs_pixel/MAXCAMPIX - 180.0;

	// convert t/t to mas between beams to be independent of pixel scale
	bp->tilt = BEAM_P2MAS(2*cx - bp->fftmag_fim.sw);
	bp->tip = BEAM_P2MAS(2*cy - bp->fftmag_fim.sh);

	// remap phase to setpoint+-180 to maximize available range before wrapping
	// bp->phase = pin (bp->phase, bp->setpoint);

	// normally we could perform unwrapping here too but the algorithms are different per beam.
}

/* copy src to the center of dst, padding with first pixel in src.
 * N.B. we assume each dimension of dst is at least as large as src
 */
static void
padFITS (FImage *dst, FImage *src)
{
	CamPix *sp = (CamPix*)src->image;
	CamPix *dp = (CamPix*)dst->image;
	int wrap = dst->sw - src->sw;
	int ndpix = dst->sw * dst->sh;
	CamPix pad = *sp;

	/* first copy pad pixel everywhere in dst */
	for (int i = 0; i < ndpix; i++)
	    *dp++ = pad;

	/* then copy src into center of dst */
	dp = (CamPix*)dst->image;
	dp += ((dst->sh - src->sh)/2)*dst->sw + (dst->sw - src->sw)/2;
	for (int i = 0; i < src->sh; i++) {
	    for (int j = 0; j < src->sw; j++)
		*dp++ = *sp++;
	    dp += wrap;
	}
}

/* read one image into fip */
static void
readNextImage (FImage *fip)
{
	if (nopchw) {
	    /* use any local FITS files found */
	    readNextLocalFITS (fip);
	} else {
	    /* read real camera image */
	    pcamReadOneImage (fip);
	}

	/* always insure beams are within a quadrant */
	forceBothBeamsFit (fip);
}

static void
quadrantShiftFFTW (fftw_complex *in, int sz)
{
	int i, npix;

	npix = sz*sz;
	for (i = 0; i < npix/2; i++) {
	    int v = quadrantShiftIndex (i, sz);
	    double tmp = in[i][0];
	    in[i][0] = in[v][0];
	    in[v][0] = tmp;
	}
	
}

/* update loopstats from both beams and the given commanded positions.
 * guard with loopstats.lock so bg loop can copy and publish without effecting our thread speed much.
 * FYI: tilt == az, tip == el
 */
static void
loopstatsUpdate (double tnow, double cmd_pist, double cmd_az, double cmd_el)
{
	double elapsedms;

	// don't let bg read loopstats while we are updating
	pthread_mutex_lock (&loopstats.lock);

	// see if it is time to compute and publish stats
	elapsedms = (tnow - loopstats.t_start)*1000.0;
	if (elapsedms >= plc_period) {

	    // compute stats, if at least a few
	    if (loopstats.n_samples > 2) {
		loopstats.p_mean = loopstats.p_sum/loopstats.n_samples;
		loopstats.p_sd = sqrt(loopstats.n_samples*loopstats.p_sum2
					- loopstats.p_sum*loopstats.p_sum)/loopstats.n_samples;

		loopstats.a_mean = loopstats.a_sum/loopstats.n_samples;
		loopstats.a_sd = sqrt(loopstats.n_samples*loopstats.a_sum2
					- loopstats.a_sum*loopstats.a_sum)/loopstats.n_samples;

		loopstats.e_mean = loopstats.e_sum/loopstats.n_samples;
		loopstats.e_sd = sqrt(loopstats.n_samples*loopstats.e_sum2
					- loopstats.e_sum*loopstats.e_sum)/loopstats.n_samples;


		loopstats.ph1_mean = loopstats.ph1_sum/loopstats.n_samples;
		loopstats.ph1_sd = sqrt(loopstats.n_samples*loopstats.ph1_sum2
					- loopstats.ph1_sum*loopstats.ph1_sum)/loopstats.n_samples;
		loopstats.cph1_mean = loopstats.cph1_sum/loopstats.n_samples;
		loopstats.sph1_mean = loopstats.sph1_sum/loopstats.n_samples;
		loopstats.snr1_mean = loopstats.snr1_sum/loopstats.n_samples;

		loopstats.ph2_mean = loopstats.ph2_sum/loopstats.n_samples;
		loopstats.ph2_sd = sqrt(loopstats.n_samples*loopstats.ph2_sum2
					- loopstats.ph2_sum*loopstats.ph2_sum)/loopstats.n_samples;
		loopstats.cph2_mean = loopstats.cph2_sum/loopstats.n_samples;
		loopstats.sph2_mean = loopstats.sph2_sum/loopstats.n_samples;
		loopstats.snr2_mean = loopstats.snr2_sum/loopstats.n_samples;

	    }

	    // reset accumulators now
	    loopstats.t_start = tnow;
	    loopstats.p_sum = 0.0;
	    loopstats.p_sum2 = 0.0;
	    loopstats.a_sum = 0.0;
	    loopstats.a_sum2 = 0.0;
	    loopstats.e_sum = 0.0;
	    loopstats.e_sum2 = 0.0;
	    loopstats.ph1_sum = 0.0;
	    loopstats.ph1_sum2 = 0.0;
	    loopstats.cph1_sum = 0.0;
	    loopstats.sph1_sum = 0.0;
	    loopstats.snr1_sum = 0.0;
	    loopstats.ph2_sum = 0.0;
	    loopstats.ph2_sum2 = 0.0;
	    loopstats.cph2_sum = 0.0;
	    loopstats.sph2_sum = 0.0;
	    loopstats.snr2_sum = 0.0;
	    loopstats.nf_jumps = 0;
	    loopstats.n_samples = 0;
	}

	// accumulate new values

	loopstats.p_sum += cmd_pist;
	loopstats.p_sum2 += cmd_pist*cmd_pist;
	loopstats.a_sum += cmd_az;
	loopstats.a_sum2 += cmd_az*cmd_az;
	loopstats.e_sum += cmd_el;
	loopstats.e_sum2 += cmd_el*cmd_el;

	loopstats.ph1_sum += beam1.phase_unw;
	loopstats.ph1_sum2 += beam1.phase_unw*beam1.phase_unw;
	loopstats.cph1_sum += cos(beam1.phase_unw*(beam1.lambda/lam_NOMIC)*(M_PI/180.0));
	loopstats.sph1_sum += sin(beam1.phase_unw*(beam1.lambda/lam_NOMIC)*(M_PI/180.0));
	loopstats.snr1_sum += beam1.tt_snr;

	loopstats.ph2_sum += beam2.phase_unw;
	loopstats.ph2_sum2 += beam2.phase_unw*beam2.phase_unw;
	loopstats.cph2_sum += cos(beam2.phase_unw*(beam2.lambda/lam_NOMIC)*(M_PI/180.0));
	loopstats.sph2_sum += sin(beam2.phase_unw*(beam2.lambda/lam_NOMIC)*(M_PI/180.0));
	loopstats.snr2_sum += beam2.tt_snr;

	// TODO loopstats.nf_jumps += cg_jump;
	loopstats.n_samples++;

	// release loopstats to bg
	pthread_mutex_unlock (&loopstats.lock);
}

/* extract the beam at x, y from raw camera image in to out.
 * set area surrounding beam to beam's mean to avoid transition ring in FFT magnitude.
 * also set center nine pixels of beam to mean to map out known UBC obscuration.
 */
static void
extractBeam (FImage *out, FImage *in, int beam_x, int beam_y)
{
	CamPix *inp = (CamPix *)in->image;
	CamPix *outp = (CamPix *)out->image;
	int inw = in->sw;
	int x, y, x0, y0;
	double sum;
	int nsum;
	CamPix mean;

	/* assert size
	 */
	if (out->sw != BEAMW)
	    bye ("extractBeam() image size error: want %d found %d\n", BEAMW, out->sw);

	/* copy beam rectangle from in to out.
	 * N.B. BEAMX/Y are FITS 1-based
	 */
	x0 = (beam_x-1) - BEAMR;
	y0 = (beam_y-1) - BEAMR;
	for (x = 0; x < BEAMW; x++)
	    for (y = 0; y < BEAMW; y++)
		outp[y*BEAMW + x] = inp[(y+y0)*inw + (x+x0)];

	/* find mean inside the beam
	 */
	sum = 0;
	nsum = 0;
	CamPix *in_beamcenter = (CamPix *) out->image;
	in_beamcenter += (BEAMW*BEAMR) + BEAMR;
	for (y = -BEAMR; y <= BEAMR; y++) {
	    for (x = -BEAMR; x <= BEAMR; x++) {
		if (x*x + y*y <= BEAMR*BEAMR) {
		    sum += in_beamcenter[y*BEAMW + x];
		    nsum++;
		}
	    }
	}
	mean = (CamPix) (sum/nsum);

	/* set everything outside the beam to the mean
	 */
	for (y = -BEAMR; y <= BEAMR; y++) {
	    for (x = -BEAMR; x <= BEAMR; x++) {
		if (x*x + y*y > BEAMR*BEAMR) {
		    outp[(y+BEAMR)*BEAMW + (x+BEAMR)] = mean;
		}
	    }
	}

	/* also set center 9 pixels of beam to mean to map out center obstruction
	 */
	if (1 || USE_UBC) {					// do it everywhere for Mar 2016 on-sky NAC
	    outp[BEAMW*(BEAMR-1) + (BEAMR-1)] = mean;
	    outp[BEAMW*(BEAMR-1) + (BEAMR  )] = mean;
	    outp[BEAMW*(BEAMR-1) + (BEAMR+1)] = mean;
	    outp[BEAMW*(BEAMR  ) + (BEAMR-1)] = mean;
	    outp[BEAMW*(BEAMR  ) + (BEAMR  )] = mean;
	    outp[BEAMW*(BEAMR  ) + (BEAMR+1)] = mean;
	    outp[BEAMW*(BEAMR+1) + (BEAMR-1)] = mean;
	    outp[BEAMW*(BEAMR+1) + (BEAMR  )] = mean;
	    outp[BEAMW*(BEAMR+1) + (BEAMR+1)] = mean;
	}
}

/* convert the pixels in fip to the real portion of fftw array out.
 * N.B. we do no size checking!
 */
static void
convertFITS2FFTW (FImage *fip, fftw_complex *out)
{
	CamPix *in = (CamPix *)fip->image;
	int i, n = fip->sw * fip->sh;

	for (i = 0; i < n; i++) {
	    out[i][0] = (double)(*in++);
	    out[i][1] = 0.0;
	}
}

/* convert the given complex array into magnitude and phase FImages.
 * images are arranged with freq [0,0] in center.
 * N.B. we scale from min to second-highest max to avoid spike at freq [0,0]
 * N.B. *in is converted from real/img to mag/phase IN PLACE.
 * N.B. we do no bounds checking.
 */
static void
convertFFTW2FITS (fftw_complex *in, FImage *mag, FImage *phs)
{
	double mmin, mmax, msmax, mrange;		/* magnitude range */
	int w = mag->sw, h = mag->sh;
	int i, npix = w * h;
	CamPix *magp = (CamPix*)mag->image;
	CamPix *phsp = (CamPix*)phs->image;

	/* convert real/img to magnitude/phase IN PLACE.
	 * also find min, max and second max in order to set range of mag pixels.
	 */ 
	mmin = DBL_MAX;
	mmax = msmax = DBL_MIN;
	for (i = 0; i < npix; i++) {
	    double real = in[i][0];
	    double img =  in[i][1];
	    double mag = sqrt(real*real + img*img);
	    in[i][0] = mag;
	    if (mag < mmin) 
		mmin = mag;
	    if (mag > mmax) {
		msmax = mmax;
		mmax = mag;
	    } else if (mag > msmax)
		msmax = mag;
	    in[i][1] = (180/M_PI)*atan2(img, real);	/* degrees -180 .. 180 */
	}
	mrange = msmax - mmin;
	if (mrange == 0)
	    mrange = 1;

	/* convert to FITS and rearrange quadrants to put freq [0,0] in center as per fftshift.
	 */
	for (i = 0; i < npix; i++) {
	    int v;
	    double t;

	    v = quadrantShiftIndex (i, w);

	    /* mag file -- mapped mmin .. msmax => 0 .. 65535 */
	    t = in[v][0];
	    if (t > msmax) t = msmax;
	    *magp++ = (CamPix)(65535.0*(t - mmin)/mrange);

	    /* phase file -- in range 0 .. 65535 */
	    *phsp++ = (CamPix)(65535.0*(in[v][1] + 180.0)/360.0);
	}
}

/* given an index and square image size (width or height),
 * return the index that performs the FFT quandrant shift.
 * N.B. assumes square image sz x sz
 */
static int
quadrantShiftIndex (int i, int sz)
{
	int sz2 = sz/2;
	int x = i%sz;
	int y = i/sz;

	// return (x == sz2-1 && y == sz2-1 ? 0 : sz*sz/2);

	if ((sz%2) == 0) {
	    // even
	    if (x >= sz2) {
		x -= sz2;
		if (y >= sz2)
		    y -= sz2;
		else
		    y += sz2;
	    } else {
		x += sz2;
		if (y >= sz2)
		    y -= sz2;
		else
		    y += sz2;
	    }
	} else {
	    // odd
	    if (x >= sz2) {
		x -= sz2;
		if (y >= sz2)
		    y -= sz2;
		else
		    y += (sz2 + 1);
	    } else {
		x += sz2 + 1;
		if (y >= sz2)
		    y -= sz2;
		else
		    y += (sz2 + 1);
	    }
	}

	return (y*sz + x);
}

/* given both beams compute the commanded tip, tilt and path length values.
 * N.B. add lots of header fields to b1p->raw_fimp ONLY.
 * return 0 if ok, -1 if image does not pass safety or quality threshold.
 */
static int
computeCorrections (OneBeam *b1p, OneBeam *b2p, double *cmd_tipp, double *cmd_tiltp, double *cmd_plp)
{
        static const double RATE_SETUP = 2;     /* seconds to wait for accurate loop rate estimate */
	static double last_tnow;		/* time of last iteration */
        static double t_closed;			/* time when loop was closed */
	static double last_used_cmd_pl;		/* use this when loop is on hold from bad stats */
	static double last_used_cmd_tip;	/* use this when loop is on hold from bad stats */
	static double last_used_cmd_tilt;	/* use this when loop is on hold from bad stats */
	static double cgx;			/* boxcar horizontal contrast gradient offset */
	static double wvx;			/* boxcar water vapor */
	static double cg_err;			/* CG control terms */
        static bool vib_ready;                  /* flag set when vib filter ready after closing loop */
        static bool difflp_ready;               /* flag set when diff lowpasas ready after closing loop */
	static BoxCar cg_boxcar;		/* form mean of several CG values */
	static BoxCar wv_boxcar;		/* form mean of several WV terms */
	static double cg_jump;			/* current amount of CG jump, microns */
	static double wv_term, wv_err;		/* mean water vapor correction before and after PWVG */
	static double t_cg_kick;		/* time when CG kick occured */
	static double p_phase;			/* pseudo-phase from water vapor, degrees */
	static double phase_gd;			/* phase group delay offset */
	static bool b2sp0_init;			/* whether beam2 sp0 has been set */
	static BoxCar b2sp0_boxcar;		/* used to average beam2 phase for sp0 */

	double tnow = secsNow();		/* now, seconds since epoch */
	double pl_err = 0;			/* phase correction this loop, degrees */
	double tip_err = 0, tilt_err = 0;	/* tip/tilt correction this loop, mas */
        double cgx_raw = 0;                     /* raw contrast gradient value */
        double wvx_raw = 0;                     /* raw water vapor */
	double cmd_tip=0, cmd_tilt=0, cmd_pl=0;	/* commanded values */
	double running_time = 0;		/* handy seconds loop has been running */
	bool cg_use = false;			/* performance switch */

	/* init_plloop is set the first time we are called after the loop is closed.
	 * use it to know to initialize all our local state. Then...
	 * ... MUST reset init_plloop before returning
	 */
	if (init_plloop) {
            last_tnow = t_closed = tnow;
	    last_used_cmd_pl = Piston_Mean = Piston_Now;
	    last_used_cmd_tip = Tip_Now;
	    last_used_cmd_tilt = Tilt_Now;
	    cgx = 0;
	    wvx = 0;
	    cg_err = 0;
            vib_ready = false;
            difflp_ready = false;
	    cg_boxcar.reset();
	    wv_boxcar.reset();
	    cg_jump = 0;
	    wv_term = wv_err = 0;
	    t_cg_kick = 0;
	    p_phase = 0;
	    phase_gd = 0;
	    b2sp0_init = false;
	    b2sp0_boxcar.reset();
	    spc_offload_net = 0;

	    // init PID loops
	    tip_pid.turnOn (tnow, Tip_Now);
	    tip_pid.setTuning (TTPGAIN, TTIGAIN, TTDGAIN);
	    tip_pid.setGainRampTime (ramp_time);
	    tip_pid.setMaxChange (max_tip_corr);

	    tilt_pid.turnOn (tnow, Tilt_Now);
	    tilt_pid.setTuning (TTPGAIN, TTIGAIN, TTDGAIN);
	    tilt_pid.setGainRampTime (ramp_time);
	    tilt_pid.setMaxChange (max_tilt_corr);

	    piston_pid.turnOn (tnow, Piston_Now);
	    piston_pid.setTuning (PLPGAIN, PLIGAIN, PLDGAIN);
	    piston_pid.setGainRampTime (ramp_time);
	    piston_pid.setMaxChange (max_pl_corr/360.0*lam1_PC);		// loop works in um
	}

	/* seconds loop has been running */
	running_time = LOOPCLOSED && t_closed > 0 ? (tnow - t_closed) : 0;

	/* set b2 baseline setpoint from initial estimate of b1-b2 difference scaled by lamda */
	if (!b2sp0_init) {
	    // don't start collecting until PIDs are fully engaged
	    if (running_time > ramp_time) {
		// find beam2 baseline SP from unwrapped phase differences adjusted for wavelength
		double phase_diff = b2p->phase_unw - b1p->phase_unw
		    - (b1p->phase_unw-b1p->setpoint) * ((b1p->lambda - b2p->lambda)/b2p->lambda);
		b2sp0_boxcar.add(phase_diff);
		if (running_time > ramp_time + b2sp0_time) {
		    b2p->sp0 = b2sp0_boxcar.mean();
		    b2sp0_init = true;
		}
	    }
	}

	/* grab a stable snapshot of OVMS data */
	OVMS ovms_copy;
	pthread_mutex_lock (&ovms.lock);
	memcpy (&ovms_copy, &ovms, sizeof(ovms_copy));
	pthread_mutex_unlock (&ovms.lock);

	/* compute TT errors using beam 1 */
	tip_err = (TIPSP - b1p->tip);				// mas
	tilt_err = (TILTSP - b1p->tilt);			// mas

	/* set up quality thresholds based on beam 1 */
	bool ttok = fabs(tip_err) < max_tip_err && fabs(tilt_err) < max_tilt_err;

	/* only update control loops if quality metrics are adequate
	 * and we are not in a hold waiting for OVMS bump, CG kick to complete or OPD freeze
	 */
	bool run_loop = !ovms_copy.bump
			&& ttok
			&& tnow > t_cg_kick + cg_freeze;
	if (run_loop) {
	    /* find instantaneous:
	     *   offset between the two unwrapped phases
	     *   pseudo phase;
	     *   water vapor
	     */
	    // phase_gd = pin (b2p->phase - b1p->phase, 0.0);
	    phase_gd = b2p->phase_unw - b1p->phase_unw;
	    p_phase = phase_gd * (b2p->lambda / (b1p->lambda - b2p->lambda) );
	    // wvx_raw = pin (b1p->phase - p_phase, 0.0);
	    wvx_raw = b1p->phase_unw - b1p->setpoint - p_phase;

	    /* average the instantaneous water vapor term */
	    wv_boxcar.add(wvx_raw);

	    /* measure and average contrast gradient within the padded image of beam 1 */
	    contrastMoment (&beam1.pad_fim, &cgx_raw);
	    cg_boxcar.add(cgx_raw);

	    /* process CG after each averaging period
	     */
	    if (cg_boxcar.age() >= cg_boxcar_period) {

		/* get mean CG then reset boxcar */
		cgx = cg_boxcar.mean();
		cg_boxcar.reset();

		/* compute CG correction and whether to use it for phase jump */
		cg_err = CGSP - cgx;
		cg_use = fabs(cg_err) > cg_min && fabs(CGSCALE) > 0;

		/* jump a full wave if CG above threshold
		 * and record time to hold off PIDs a few frames
		 */
		if (cg_use) {
                    cg_jump += ((cg_err < 0) ? lam1_PC: -lam1_PC) * CGSCALE;	// microns
		    t_cg_kick = tnow;
		}
	    }

	    /* process WV after each averaging period
	     */
	    if (wv_boxcar.age() >= wv_boxcar_period) {

		/* get mean PWV term then reset boxcar */
		wv_term = wv_boxcar.mean();
		wv_boxcar.reset();

		/* update set point correction due to water vapor */
		wv_err = PWVG*wv_term;
	    }

	    /* apply WV term to b1p */
	    b1p->setpoint = b1p->sp0 + wv_err;

	    /* update beam2 setpoint to follow beam1 if baseline sp is ready */
	    if (b2sp0_init)
		b2p->setpoint = b2p->sp0 + (b1p->setpoint - b1p->sp0) * (b1p->lambda/b2p->lambda);

	    /* net path length error is based solely on beam1 */
	    pl_err = b1p->setpoint - b1p->phase_unw;

	    /* inject error into PID to compute absolute commmanded path length */
	    cmd_pl = piston_pid.run (tnow, pl_err/360.0*lam1_PC);		// degrees to microns

	    /* compute commanded tip and tilt control entirely based on beam1 */
	    cmd_tip = tip_pid.run (tnow, tip_err*0.001);	// want asec from mas
	    cmd_tilt = tilt_pid.run (tnow, tilt_err*0.001);	// want asec from mas

	    /* record values for reuse in case next frame is rejected */
	    last_used_cmd_pl = cmd_pl;
	    last_used_cmd_tip = cmd_tip;
	    last_used_cmd_tilt = cmd_tilt;

	} else {

	    /* this image does not qualify, start with last-good positions again
	     */
	    cmd_pl = last_used_cmd_pl;
	    cmd_tip = last_used_cmd_tip;
	    cmd_tilt = last_used_cmd_tilt;

	}

	/* ramp up Vib filter gains to match PID controller.
	 * this also gives filters a few iterations to fill history before becoming effective.
	 */
	double gramp;
	if (ramp_time > 0) {
	    gramp = running_time/ramp_time;
	    if (gramp > 1)
		gramp = 1;
	} else
	    gramp = 1;

	/* add vibration filter feedforward after waiting
         * to get good estimate of loop rate RATE_PLCLI.
	 */
        double vib_corr = 0;
        if (running_time > RATE_SETUP) {
            if (!vib_ready) {
                setupVibFilterBiQuad();
                vib_ready = true;
            } else {
                vib_corr = VIB_ZGAIN * vib_bqfilter.run(b1p->phase_unw) / VIB_VGAIN * gramp;
                cmd_pl += vib_corr;
            }
        }


	/* add OVMS feedforward only if enabled.
	 * apply to Pupil if NAC, else add to FPC commanded
	 */
	if (!noovms && ovms_copy.ok && !ovms_copy.bump) {
	    double o_tip = ovms_copy.tip * ovms_copy.tip_gain;
	    double o_tilt = ovms_copy.tilt * ovms_copy.tilt_gain;
	    double o_opd = ovms_copy.opd * ovms_copy.opd_gain;
	    if (USE_UBC) {
		cmd_tip += o_tip;
		cmd_tilt += o_tilt;
		cmd_pl += o_opd;
	    } else {
		if (pzt_f && LOOPCLOSED)
		    setPupil (pzt_f, o_tip, o_tilt, o_opd);
	    }
	}

	/* define differentiator low pass filter after we have waited long enough to have a good
	 * estimate of RATE_PLCLI
	 */
        if (!difflp_ready && running_time > RATE_SETUP) {
	    piston_pid.setDiffLP (plpfc, loopinfo.np[RATE_PLCLI].value, 0.7); // 0.7 = flat
	    difflp_ready = true;
	}

        /* add autobeam corrections if enabled
         */
        double autobeam_x = 0, autobeam_y = 0;
        if (AUTOBEAM_ON) {
            int minp, maxp;
            findBeam (b1p->raw_fimp, &minp, &maxp, &autobeam_x, &autobeam_y);
            autobeam_x_rb->add(autobeam_x);
            autobeam_y_rb->add(autobeam_y);

            // TODO: move SPC
            // double meancx = autobeam_x_rb->mean();
            // double meancy = autobeam_y_rb->mean();
        }

	/* keep loopstats honest, even if image does not qualify for a move update */
	loopstatsUpdate (tnow, cmd_pl, cmd_tilt, cmd_tip);

	/* history */
	double loopdt = (tnow - last_tnow)*1000.0;	// time since last loop, ms
	last_tnow = tnow;

	/* reset init flag */
	init_plloop = false;

	/* store in beam1 raw header */
	FImage *fip = b1p->raw_fimp;

	setCommentFITS (fip, "COMMENT", "--------- Operational data --------------------------------------");
	addRealFITS (fip, "SWVERSN",   plcversion.np[DRIVER_PVER].value, 16, "INDI Driver version number");
	addRealFITS (fip, jd_kw,       secs2JD(tnow), 16, "JD now");
	addRealFITS (fip, "PCB1X",     BEAM1X, 6, "Beam 1 center X");
	addRealFITS (fip, "PCB1Y",     BEAM1Y, 6, "Beam 1 center Y");
	addRealFITS (fip, "PCB2X",     BEAM2X, 6, "Beam 2 center X");
	addRealFITS (fip, "PCB2Y",     BEAM2Y, 6, "Beam 2 center Y");
	addRealFITS (fip, "PCBRAD",    BEAMR, 6, "Beam radius");
	addRealFITS (fip, "PCB1LAM",   b1p->lambda, 6, "Beam 1 wavelength, nm");
	addIntFITS  (fip, "PCRWMAG1",  b1p->tt_pixel, "Raw FFT Magnitude pixel value, beam 1, 0..65k");
	addIntFITS  (fip, "PCRWPHS1",  b1p->phs_pixel, "Raw FFT Phase pixel value @ Mag, beam 1, 0..65k");
	addRealFITS (fip, "PCB2LAM",   b2p->lambda, 6, "Beam 2 wavelength, nm");
	addIntFITS  (fip, "PCRWMAG2",  b2p->tt_pixel, "Raw FFT Magnitude pixel value, beam 2, 0..65k");
	addIntFITS  (fip, "PCRWPHS2",  b2p->phs_pixel, "Raw FFT Phase pixel value @ Mag, beam 2, 0..65k");
	addIntFITS  (fip, "PCCLOSED",  LOOPCLOSED, "Whether loop is ostensibly closed");
	addIntFITS  (fip, "PCRUNTHS",  run_loop, "Whether loop ran this frame");
	addIntFITS  (fip, "PCLSTATP",  plc_period, "Closed loop stats collection period, ms");
	addRealFITS (fip, "PCIMGDT",   loopdt, 6, "Time since previous image, msec");
	addRealFITS (fip, "PCLOOPHZ",  loopinfo.np[RATE_PLCLI].value, 6, "Closed loop rate, Hz");
	addRealFITS (fip, "PCRUNTIM",  running_time, 6, "Elapsed time since loop was closed, secs");
	addIntFITS  (fip, "PCBMPAD",   beam_pad, "Extracted beam padding factor before FFT");
	addIntFITS  (fip, "PCFFTSZ",   PADIMGSZ, "FFT size");
	addStringFITS (fip, "PCSRC",   USE_UBC ? "UBC" : "NAC", "Source");

	setCommentFITS (fip, "COMMENT", "--------- Path length information -------------------------------");
	addRealFITS (fip, "PCPLSP1",   b1p->setpoint, 6, "Path length setpoint, beam 1, degrees");
	addRealFITS (fip, "PCPLSP01",  b1p->sp0, 6, "Baseline path length setpoint, beam 1, degrees");
	addRealFITS (fip, "PCFTPHS1",  b1p->phase, 6, "Raw Phase value, beam 1, degrees, PCPLSP1 +-180");
	addRealFITS (fip, "PCUNWPH1",  b1p->phase_unw, 6, "Unwrapped phase, beam 1, degrees");
	addRealFITS (fip, "PCUNPMN1",  loopstats.ph1_mean, 6, "Unwrapped phase mean, beam 1, degrees");
	addRealFITS (fip, "PCUNPSD1",  loopstats.ph1_sd, 6, "Unwrapped phase StdDev, beam 1, degs");
	addRealFITS (fip, "PCPLSP2",   b2p->setpoint, 6, "Path length setpoint, beam 2, degrees");
	addRealFITS (fip, "PCPLSP02",  b2p->sp0, 6, "Baseline path length setpoint, beam 2, degrees");
	addRealFITS (fip, "PCFTPHS2",  b2p->phase, 6, "Raw Phase value, beam 2, degrees, PCPLSP2 +-180");
	addRealFITS (fip, "PCUNWPH2",  b2p->phase_unw, 6, "Unwrapped phase, beam 2, degrees");
	addRealFITS (fip, "PCUNPMN2",  loopstats.ph2_mean, 6, "Unwrapped phase mean, beam 2, degrees");
	addRealFITS (fip, "PCUNPSD2",  loopstats.ph2_sd, 6, "Unwrapped phase StdDev, beam 2, degs");
	addRealFITS (fip, "PCPERR",    pl_err, 6, "Path length error, degrees");

	addRealFITS (fip, "PCPLPG",    (USE_UBC
					    ? ubcsettings.np[PLPGAIN_UBCS].value
					    : nacsettings.np[PLPGAIN_NACS].value),
				       6, "Path length proportional gain");

	addRealFITS (fip, "PCPLDG",    (USE_UBC
					    ? ubcsettings.np[PLDGAIN_UBCS].value
					    : nacsettings.np[PLDGAIN_NACS].value),
				       6, "Path length derivative gain");

	addRealFITS (fip, "PCPLIG",    (USE_UBC
					    ? ubcsettings.np[PLIGAIN_UBCS].value
					    : nacsettings.np[PLIGAIN_NACS].value),
				       6, "Path length integral gain");

	addRealFITS (fip, "PCPLDFC",   difflp_ready ? plpfc : 0, 6,
						"Path length derivative filter cutoff, Hz, or 0 if off");
	addRealFITS (fip, "PCNETOFL",  spc_offload_net, 6, "Net offloading this run, microns");


	setCommentFITS (fip, "COMMENT", "--------- Contrast Gradient information -------------------------");
	addRealFITS (fip, "PCCGVAL",   cgx_raw, 6, "Raw Contrast Gradient value this frame");
	addRealFITS (fip, "PCCGVALA",  cgx, 6, "PCCGVAL boxcar average");
	addRealFITS (fip, "PCCGVALP",  cg_boxcar_period, 6, "PCCGVAL boxcar period, sec");
	addRealFITS (fip, "PCCGSP",    CGSP, 6, "Contrast gradient setpoint");
	addRealFITS (fip, "PCCGERR",   cg_err, 6, "Contrast Gradient difference from setpoint");
	addRealFITS (fip, "PCCGMIN",   cg_min, 6, "Contrast Gradient threshold");
	addRealFITS (fip, "PCCGSCAL",  CGSCALE, 6, "Contrast Gradient jump scale multiplier");
	addRealFITS (fip, "PCCGJUMP",  cg_jump, 6, "Net CG jump correction, microns");

	setCommentFITS (fip, "COMMENT", "--------- Water Vapor information -------------------------------");
	addRealFITS (fip, "PCPWVG",    PWVG, 6, "Water vapor correction gain factor");
	addRealFITS (fip, "PCPWVGD",   phase_gd, 6, "Group delay: PCFTPHS2 - PCFTPHS1, degs");
	addRealFITS (fip, "PCPWVPP",   p_phase, 6, "Pseudo phase: PCPWVGD * PCB2LAM/(PCB1LAM-PCB2LAM");
	addRealFITS (fip, "PCPWVRAW",  wvx_raw, 6, "Water vapor term this frame: PCFTPHS1 - PCPWVPP");
	addRealFITS (fip, "PCPWVT",    wv_term, 6, "PCPWVRAW boxcar average");
	addRealFITS (fip, "PCPWVTP",   wv_boxcar_period, 6, "Water vapor boxcar period, secs");
	addRealFITS (fip, "PCPWVER",   wv_err, 6, "Water vapor set point change: PCPWVT * PCPWVG");

	setCommentFITS (fip, "COMMENT", "--------- Tip information ---------------------------------------");
	addRealFITS (fip, "PCFTMGY1",  b1p->tip, 6, "FFT measured tip (el), beam 1, mas");
	addRealFITS (fip, "PCFTMGY2",  b2p->tip, 6, "FFT measured tip (el), beam 2, mas");
	addRealFITS (fip, "PCTIPSP",   TIPSP, 6, "Tip setpoint, mas");
	addRealFITS (fip, "PCTIPER",   tip_err, 6, "Tip error, mas");
	addRealFITS (fip, "PCTIPSTD",  loopstats.e_sd, 6, "Tip error std dev, mas");

	addRealFITS (fip, "PCTTPG",    (USE_UBC
					    ? ubcsettings.np[TTPGAIN_UBCS].value
					    : nacsettings.np[TTPGAIN_NACS].value),
					6, "Tip and Tilt propotional gain");

	addRealFITS (fip, "PCTTIG",    (USE_UBC
					    ? ubcsettings.np[TTIGAIN_UBCS].value
					    : nacsettings.np[TTIGAIN_NACS].value),
					6, "Tip and Tilt integral gain");

	addRealFITS (fip, "PCTTDG",    (USE_UBC
					    ? ubcsettings.np[TTDGAIN_UBCS].value
					    : nacsettings.np[TTDGAIN_NACS].value),
					6, "Tip and Tilt derivative gain");


	setCommentFITS (fip, "COMMENT", "--------- Tilt information --------------------------------------");
	addRealFITS (fip, "PCFTMGX1",  b1p->tilt, 6, "FFT measured tilt (az), beam 1, mas");
	addRealFITS (fip, "PCFTMGX2",  b2p->tilt, 6, "FFT measured tilt (az), beam 2, mas");
	addRealFITS (fip, "PCTLTSP",   TILTSP, 6, "Tilt setpoint, mas");
	addRealFITS (fip, "PCTLTER",   tilt_err, 6, "Tilt error, mas");
	addRealFITS (fip, "PCTLTSTD",  loopstats.a_sd, 6, "Tilt error std dev, mas");

	setCommentFITS (fip, "COMMENT", "--------- Quality metrics ---------------------------------------");
	addIntFITS  (fip, "PCTTOK",    ttok, "Whether PCTIPER < PCMAXTPE and PCTLTER < PCMAXTLE, 0 or 1");
	addRealFITS (fip, "PCMFSNR",   b1p->tt_snr, 6, "FFT Mag SNR");
	addRealFITS (fip, "PCFMINS",   MINFSNR, 6, "Minimum FFT Mag SNR");
	addRealFITS (fip, "PCMAXPLC",  max_pl_corr, 6, "Max allowed PL correction per frame, degrees");
	addRealFITS (fip, "PCBADPS",   bad_space_std, 6, "Bad spacial pixel std dev cutoff");
	addRealFITS (fip, "PCBADPT",   bad_time_std, 6, "Bad temporal pixel std dev cutoff");
	addRealFITS (fip, "PCRAMPDT",  ramp_time, 6, "Gain ramp time, seconds");

	setCommentFITS (fip, "COMMENT", "--------- Commanded values --------------------------------------");
	if (LOOPCLOSED) {
	    addRealFITS (fip, "PCPLCMD",   cmd_pl, 6, "Commanded piston, microns");
	    addRealFITS (fip, "PCTIPCMD",  cmd_tip, 6,  "Commanded tip, arcsec");
	    addRealFITS (fip, "PCTLTCMD",  cmd_tilt, 6, "Commanded tilt, arcsec");
	} else {
	    addRealFITS (fip, "PCPLCMD",   cmd_pl, 6, "Commanded piston if loop was closed, microns");
	    addRealFITS (fip, "PCTIPCMD",  cmd_tip, 6,  "Commanded tip if loop was closed, arcsec");
	    addRealFITS (fip, "PCTLTCMD",  cmd_tilt, 6, "Commanded tilt if loop was closed, arcsec");
	}

        if (AUTOBEAM_ON) {
	    setCommentFITS (fip, "COMMENT", "--------- Autobeam ----------------------------------------------");
            addRealFITS (fip, "PCABIX",   autobeam_x, 6, "Autobeam X this image");
            addRealFITS (fip, "PCABIY",   autobeam_y, 6, "Autobeam Y this image");
            addRealFITS (fip, "PCABMX",   autobeam_x_rb->mean(), 6, "Autobeam X running mean");
            addRealFITS (fip, "PCABMY",   autobeam_y_rb->mean(), 6, "Autobeam Y running beam");
        }

	if (!noovms) {
	    setCommentFITS (fip, "COMMENT", "--------- OVMS --------------------------------------------------");
	    addIntFITS  (fip, "PCOVOK",  ovms_copy.ok, "whether OVMS is working properly");
	    addIntFITS  (fip, "PCOVBMP", ovms_copy.bump, "whether OVMS is detecting a bump");
	    addRealFITS (fip, "PCOVAGE", (tnow - ovms_copy.timestamp)*1e3, 6, "OVMS packet latency, ms");
	    addRealFITS (fip, "PCOVHZN", ovms_copy.horizon, 6, "OVMS prediction horizon, ms");
	    addRealFITS (fip, "PCOVPLG", ovms_copy.opd_gain/OVMS_OPD_GAIN_SF, 6, "OVMS opd gain");
	    addRealFITS (fip, "PCOVPLV", ovms_copy.opd, 6, "OVMS opd value");
	    addRealFITS (fip, "PCOVTPG", ovms_copy.tip_gain/OVMS_TIP_GAIN_SF, 6, "OVMS tip gain");
	    addRealFITS (fip, "PCOVTPV", ovms_copy.tip, 6, "OVMS tip value");
	    addRealFITS (fip, "PCOVTLG", ovms_copy.tilt_gain/OVMS_TILT_GAIN_SF, 6, "OVMS tilt gain");
	    addRealFITS (fip, "PCOVTLV", ovms_copy.tilt, 6, "OVMS tilt value");
	}

        setCommentFITS (fip, "COMMENT", "--------- Vibration filter --------------------------------------");
        addRealFITS (fip, "PCVIBFFC",  VIB_FC, 6, "Vibration filter center frequency, Hz");
        addRealFITS (fip, "PCVIBFQ",   VIB_Q, 6, "Vibration filter Q ~ Fc/BW");
        addRealFITS (fip, "PCVIBFFG",  VIB_FGAIN, 6, "Vibration filter gain, dB");
        addRealFITS (fip, "PCVIBFVG",  VIB_VGAIN, 6, "Vibration filter gain, voltage");
        addRealFITS (fip, "PCVIBFZG",  VIB_ZGAIN, 6, "Vibration filter piezo gain");
        addRealFITS (fip, "PCVCORR",   vib_corr, 6, "Vibration correction, degrees");

	/* always return ok so filters keep running.
	 * if quality metrics not good new values just repeat from previous.
	 */
	if (cmd_tipp)
	    *cmd_tipp = cmd_tip;
	if (cmd_tiltp)
	    *cmd_tiltp = cmd_tilt;
	if (cmd_plp)
	    *cmd_plp = cmd_pl;

	return (0);
}

/* return v after possibly adjusting by +-360 as necessary so it is within center+-180
 */
static double
pin (double v, double center)
{
    	while (v > center+180)
	    v -= 360;
	while (v < center-180)
	    v += 360;
	return (v);
}

/* Find brighest pixel location in right half subimage of FFT mag for the given beam
 */
static void
FFTPeakRHP (OneBeam *bp, double *cxp, double *cyp)
{
	/* handy
	 */
	FImage *fip = &bp->fftmag_fim;
	int w = fip->sw, h = fip->sh;

	/* analysis is confined to columns to the right of this left edge
	 */
        const int left = peakignore*w;

	/* find location of maximum pixel on right side of left boundary
	 */
	int maxx = 0, maxy = 0;
	CamPix maxp = 0;
	CamPix *ip = (CamPix*)fip->image;
	for (int y = 0; y < h; y++) {
	    for (int x = 0; x < w; x++) {
		CamPix p = *ip++;
		if (x >= left) {
		    if (p > maxp) {
			maxp = p;
			maxx = x;
			maxy = y;
		    }
		}
	    }
	}

	/* find sum inside box that surrounds the brightest pixel and stats outside
	 */
	const int HWPSF = 3;
	ip = (CamPix*)fip->image;
	double bkgsum = 0, bkgsum2 = 0, boxsum = 0;
	int ninbox = 0, noutbox = 0;
	for (int y = 0; y < h; y++) {
	    for (int x = 0; x < w; x++) {
		CamPix p = *ip++;
		if (x >= left) {
		    if (x <= maxx+HWPSF && x >= maxx-HWPSF && y <= maxy+HWPSF && y >= maxy-HWPSF) {
			// signal+noise inside box
			boxsum += p;
			ninbox++;
		    } else {
			// background stats outside box
			bkgsum += p;
			bkgsum2 += (double)p*(double)p;
			noutbox++;
		    }
		}
	    }
	}

	/* compute SNR of brightest pixel compared to surrounding noise
	 */
	double bgmean = bkgsum/noutbox;
	double bgstddev = sqrt(noutbox*bkgsum2 - bkgsum*bkgsum)/noutbox;
	bp->tt_pixel = boxsum/ninbox;
	bp->tt_snr = bgstddev > 0 ? (bp->tt_pixel-bgmean)/bgstddev : 0;

#ifdef CENTROID_PEAK
	/* N.B. This algorithm tends to jump from one pixel to the next unless the box radius is > FWHM.
	 */

	/* refine peak location as centroid around the brightest pixel
	 */
	ip = (CamPix *)fip->image;
	ip += maxy*w + maxx;
	double sump = 0, sumxp = 0, sumyp = 0;
	const int R = 7;
	for (int y = -R; y <= R; y++) {
	    for (int x = -R; x <= R; x++) {
		double p = (double) ip[y*w + x];
		sumxp += x*p;
		sumyp += y*p;
		sump += p;
	    }
	}

	/* reported center is location of peak modified by centroid center.
	 * N.B. beware all-black frame, eg, if dark frame is brighter everywhere than current image
	 */
	*cxp = sump ? (maxx + sumxp/sump) : 0;
	*cyp = sump ? (maxy + sumyp/sump) : 0;

#else	// finite difference in each direction
	ip = (CamPix *)fip->image;
	ip += maxy*w + maxx;

	double slope, x_offset, y_offset;
	slope = (double)ip[1] + ip[-1] - 2.0*ip[0];
	x_offset = slope == 0 ? 0 : -0.5*(ip[1] - ip[-1]) / slope;
	slope = (double)ip[w] + ip[-w] - 2.0*ip[0];
	y_offset = slope == 0 ? 0 : -0.5*(ip[w] - ip[-w]) / slope;

	*cxp = maxx + x_offset;
	*cyp = maxy + y_offset;
#endif
}

/* find the value of an image at coordinates [x,y].
 */
static CamPix
imageValue (FImage *fip, int x, int y)
{
	CamPix *ip = (CamPix*)fip->image;
	int w = fip->sw;
	return (ip[y*w + x]);
}

/* given image containing centered beam and surrounded by its mean,
 * return a point that indicates which left-right direction has more
 * contrast as a means to assess location along group.
 * N.B. we assume the corner pixel is the beam mean.
 */
static void
contrastMoment (FImage *fip, double *xp)
{
	CamPix *ip = (CamPix*)fip->image;
	int w = fip->sw;
	int h = fip->sh;
	double sumxp, sump;
	double mean;
	int x, y;

	// mean is corner pixel
	mean = ip[0];

	// move ip to beam center
	ip += (h/2)*w + w/2;

	// find first contrast moment about the mean
	sumxp = sump = 0.0;
	for (y = -BEAMR; y <= BEAMR; y++) {
	    for (x = -BEAMR; x <= BEAMR; x++) {
		if (x*x + y*y <= BEAMR*BEAMR) {
		    double p = fabs(ip[y*w + x] - mean);
		    //xp = sqrt(x*x+y*y)*cos(atan2(y,x)-tip) //x when CG projected along fringes
		    sumxp += x*p;
		    sump += p;
		}
	    }
	}

	// return location WRT beam center.
	*xp = sumxp/sump;
}


/* set up for (re)using fftw_execute to size w and h.
 */
static void
setupFFTW (OneBeam *bp, int w, int h)
{
	int npix = w*h;

	/* create FFTW pieces */
	if (bp->fft_input)
	    fftw_free (bp->fft_input);
	bp->fft_input = (fftw_complex*) fftw_malloc (sizeof(fftw_complex)*npix);

	if (bp->fft_output)
	    fftw_free (bp->fft_output);
	bp->fft_output = (fftw_complex*) fftw_malloc (sizeof(fftw_complex)*npix);

	IDLog ("Creating new FFTW plan for %d x %d\n", w, h);
	if (bp->fft_theplan)
	    fftw_destroy_plan(bp->fft_theplan);
	bp->fft_theplan = fftw_plan_dft_2d (w, h, bp->fft_input, bp->fft_output, FFTW_FORWARD, FFTW_MEASURE);
	if (!bp->fft_theplan)
	    bye ("Could not create FFTW plan for %d x %d\n", w, h);
}

/* set up a valid FITS file of given size, just malloc image but don't init
 */
static void
setupFITS (FImage *fip, int w, int h)
{
	int npix = w*h;
	int nimbytes = npix*sizeof(CamPix);
	int nhdr = FITS_HROWS*FITS_HCOLS;
	char fits[nhdr + nimbytes], *ffp;
	char buf[1024];

	/* start with minimal FITS headers */
	ffp = fits;
	ffp += sprintf (ffp, "SIMPLE  = %20s%-50s", "T",     "");
	ffp += sprintf (ffp, "BITPIX  = %20d%-50s", 16,      " / bits/pix");
	ffp += sprintf (ffp, "NAXIS   = %20d%-50s", 2,       " / n image axes");
	ffp += sprintf (ffp, "NAXIS1  = %20d%-50s", w,       " / n columns");
	ffp += sprintf (ffp, "NAXIS2  = %20d%-50s", h,       " / n rows");
	ffp += sprintf (ffp, "BZERO   = %20d%-50s", 32768,   " / physical zero");
	ffp += sprintf (ffp, "BSCALE  = %20d%-50s", 1,       " / physical scale");
	ffp += sprintf (ffp, "%-79s", "END");
	while (ffp < fits + nhdr)
	    *ffp++ = ' ';

	/* convert fits in-memory image to FImage format -- pixel values don't matter yet */
	resetFImage (fip);
	if (readFITSmem (fits, nhdr+nimbytes, fip, buf) < 0)
	    bye ("Failed to create beam 1 image: %s\n", buf);
}

/* force each beam to fit in a quadrant of fip (ok if different quadrants).
 * if any changes are required, publish new positions and size
 */
static void
forceBothBeamsFit (FImage *fip)
{
	// always publish the first check
	static bool first_time = true;

	bool changed = false;
	int newr = BEAMR;
	int new1x, new1y;
	int new2x, new2y;

	// check beam1
	new1x = BEAM1X;
	new1y = BEAM1Y;
	forceBeamFit (fip, new1x, new1y, newr);
	if (new1x != BEAM1X || new1y != BEAM1Y || newr != BEAMR)
	    changed = true;

	// check beam2
	new2x = BEAM2X;
	new2y = BEAM2Y;
	forceBeamFit (fip, new2x, new2y, newr);
	if (new2x != BEAM2X || new2y != BEAM2Y || newr != BEAMR)
	    changed = true;

	// publish if first time or anything changed
	if (first_time || changed) {

	    // guard from INDI thread updates
	    pthread_mutex_lock (&beam_lock);

	    BEAM1X = new1x;
	    BEAM1Y = new1y;
	    BEAM2X = new2x;
	    BEAM2Y = new2y;
	    BEAMR = newr;

	    if (USE_UBC) {
		ubcsettings.np[BEAM1X_UBCS].value = BEAM1X;
		ubcsettings.np[BEAM1Y_UBCS].value = BEAM1Y;
		ubcsettings.np[BEAM2X_UBCS].value = BEAM2X;
		ubcsettings.np[BEAM2Y_UBCS].value = BEAM2Y;
		ubcsettings.np[BEAMR_UBCS].value  = BEAMR;
		IDSetNumber (&ubcsettings, NULL);
	    } else {
		nacsettings.np[BEAM1X_NACS].value = BEAM1X;
		nacsettings.np[BEAM1Y_NACS].value = BEAM1Y;
		nacsettings.np[BEAM2X_NACS].value = BEAM2X;
		nacsettings.np[BEAM2Y_NACS].value = BEAM2Y;
		nacsettings.np[BEAMR_NACS].value  = BEAMR;
		IDSetNumber (&nacsettings, NULL);
	    }

	    pthread_mutex_unlock (&beam_lock);

	}

	// no longer forced
	first_time = false;
}

/* given a beam center location, force it to fit within one quadrant of fip.
 * N.B. incoming coords are FITS 1-based
 * N.B. we assume fip is square
 */
static void
forceBeamFit (FImage *fip, int& fits_x, int& fits_y, int& r)
{
	/* handy quadrant size */
	int qsz = fip->sw/2;

	/* change to 0-based */
	int x = fits_x - 1;
	int y = fits_y - 1;

	/* anchor which quadrant contains center */
	int startx = x < qsz ? 0 : qsz;
	int starty = y < qsz ? 0 : qsz;
	int endx = startx + qsz - 1;	// last valid index
	int endy = starty + qsz - 1;	//   "

	/* constrain radius to fit inside quadrant */
	if (2*r+1 > qsz)
	    r = (qsz-1)/2;

	/* constrain location to fit inside quadrant */
	if (x - r < startx)
	    x = startx + r;
	else if (x + r > endx)
	    x = endx - r;
	if (y - r < starty)
	    y = starty + r;
	else if (y + r > endy)
	    y = endy - r;

	/* back to fits */
	fits_x = x + 1;
	fits_y = y + 1;
}

/* write the given FITS file to name as per jdFilename().
 * N.B. fip is NO LONGER VALID on return.
 * exit if trouble.
 */
static void
saveFITS (FImage *fip)
{
	char full[1024], msg[1024];
	double jd;
	int fd;

	/* get JD from file if possible */
	if (getRealFITS (fip, jd_kw, &jd) < 0)
	    jd = jdNow();

	/* write file */
	jdFilename (jd, full);
	fd = open (full, O_RDWR|O_CREAT, S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH);
	if (fd < 0)
	    bye ("open(%s): %s\n", full, strerror(errno));
	if (writeFITS (fd, fip, msg, 0) < 0)
	    bye ("write(%s): %s\n", full, msg);
	close (fd);
}

/* q fip header for saving later to fitshdr_q.
 * exit if trouble.
 */
static void
saveFITSHeader (FImage *fip)
{
	pthread_mutex_lock (&fitshdrq_lock);
	int nbytes = fip->nvar * sizeof(FITSRow);
	HDRVar *hvp = (HDRVar *) malloc (sizeof(HDRVar));
	hvp->var = (FITSRow*) memcpy (malloc (nbytes), fip->var, nbytes);
	hvp->nvar = fip->nvar;
	if (getRealFITS (fip, jd_kw, &hvp->jd) < 0)
	    hvp->jd = jdNow();
	pushFQ (fitshdr_q, hvp);
	pthread_mutex_unlock (&fitshdrq_lock);
}

/* called periodically to write any pending fits q entries
 */
static void
saveFITSHeaderQCB (void *not_used)
{
	/* drain q */
	while (1) {
	    pthread_mutex_lock (&fitshdrq_lock);
	    HDRVar *hvp = (HDRVar*) popFQ (fitshdr_q);
	    pthread_mutex_unlock (&fitshdrq_lock);
	    if (!hvp)
		break;

	    char full[1024];
	    double jd = hvp->jd;
	    char *us;
	    FILE *fp;
	    struct stat s;
	    int i, nf;


	    /* start with image file name now */
	    jdFilename (jd, full);

	    /* chop off minute portion of name */
	    us = strrchr (full, '_');
	    if (!us)
		bye ("BUG! file name does not contain underscore: '%s'\n", full);
	    while (--us > full)
		if (*us == '_')
		    break;
	    if (us == full)
		bye ("BUG! can not find minute position in file name: %s\n", full);

	    /* set new suffix */
	    sprintf (us, ".txt");

	    /* if file does not exist, create header line */
	    if (stat (full, &s) < 0) {
		fp = fopen (full, "w");
		if (!fp)
		    bye ("fopen(%s,\"w\"): %s\n", full, strerror(errno));
		fprintf (fp, "# ");
		int fieldn = 1;
		for (nf = i = 0; i < hvp->nvar; i++) {
		    char *row = hvp->var[i];
		    char name[32]; double value;
		    if (sscanf (row, "%8s = %lf", name, &value) == 2) {
			/* field is numeric, so add name to header */
			if (nf++ > 0)
			    fprintf (fp, ",");
			fprintf (fp, "%d-%s", fieldn++, name);
		    }
		}
		fprintf (fp, "\n");
		fclose (fp);
	    }

	    /* open for append */
	    fp = fopen (full, "a");
	    if (!fp)
		bye ("fopen(%s,\"a\"): %s\n", full, strerror(errno));

	    /* add each numeric field */
	    for (nf = i = 0; i < hvp->nvar; i++) {
		char *row = hvp->var[i];
		char name[32]; double value;
		if (sscanf (row, "%8s = %lf", name, &value) == 2) {
		    /* field is numeric, but we want the value as a string to preserve precision */
		    char *eq = strchr (row, '=');
		    char value_str[128];
		    if (sscanf (eq, "= %s", value_str) == 1) {
			if (nf++ > 0)
			    fprintf (fp, ",");
			fprintf (fp, "%s", value_str);
		    }
		}
	    }
	    fprintf (fp, "\n");

	    /* clean up */
	    fclose (fp);
	    free (hvp->var);
	    free (hvp);
	}

	/* repeat */
	IEAddTimer (saveFHms, saveFITSHeaderQCB, NULL);
}

/* this callback is called every plc_period ms to publish closed-loop info.
 */
static void
loopStatusCB(void *not_used)
{
	static double last_t;
	static int last_plc_loops;
	LoopStats lstmp;
	double t;
	double dt;

	/* grab a local copy of the loopstats as quickly as possible so as not to
	 * interfere with the real-time thread
	 */
	pthread_mutex_lock (&loopstats.lock);
	memcpy (&lstmp, &loopstats, sizeof(lstmp));
	pthread_mutex_unlock (&loopstats.lock);

	/* time now */
	t = secsNow();

	/* assume loop has started over if new count is less than old */
	if (loopstats.plc_loops < last_plc_loops) {
	    last_plc_loops = loopstats.plc_loops;
	    last_t = t;
	}

	/* compute plcLoopThread() rate over a little longer period */
	dt = t - last_t;
	if (dt > .5) {
	    double rate = (lstmp.plc_loops - last_plc_loops)/dt;
	    loopinfo.np[RATE_PLCLI].value = rate;
	    last_plc_loops = lstmp.plc_loops;
	    last_t = t;

	    // stop disk if too fast
	    if (rate > 3000) {
		if (SAVEFEVERY > 0) {
		    IUSaveText (&savefits.tp[EVERY_SVF], "0");
		    IDSetText (&savefits, "Too fast -- disabling saving FITS");
		}
		if (SAVESTATS) {
		    savestats.sp[YES_SVST].s = ISS_OFF;
		    IDSetSwitch (&savestats, "Too fast -- disabling saving stats");
		}
	    }
	}

	publishLoopInfo (&lstmp);

	/* now publish info from the FITS header.
	 * N.B. must sync with stable bg_cam_fim
	 */
	pthread_mutex_lock (&plc_bglock);
	publishHeader (&bg_cam_fim);
	pthread_mutex_unlock (&plc_bglock);

	/* repeat */
	IEAddTimer (plc_period, loopStatusCB, NULL);
}

/* publish loopstats to loopinfo INDI property
 */
static void
publishLoopInfo (LoopStats *ls)
{

	loopinfo.np[PISTM_PLCLI].value = ls->p_mean;
	loopinfo.np[PISTS_PLCLI].value = ls->p_sd;
	loopinfo.np[AZM_PLCLI].value = ls->a_mean;
	loopinfo.np[AZS_PLCLI].value = ls->a_sd;
	loopinfo.np[ELM_PLCLI].value = ls->e_mean;
	loopinfo.np[ELS_PLCLI].value = ls->e_sd;

	loopinfo.np[PHMEAN1_PLCLI].value = ls->ph1_mean;
	loopinfo.np[PHSD1_PLCLI].value = ls->ph1_sd;
	loopinfo.np[PHMC1_PLCLI].value = ls->cph1_mean;
	loopinfo.np[PHMS1_PLCLI].value = ls->sph1_mean;
	loopinfo.np[SNR1_PLCLI].value = ls->snr1_mean;
	loopinfo.np[PCLAM1_PLCLI].value = lam1_PC;

	loopinfo.np[PHMEAN2_PLCLI].value = ls->ph2_mean;
	loopinfo.np[PHSD2_PLCLI].value = ls->ph2_sd;
	loopinfo.np[PHMC2_PLCLI].value = ls->cph2_mean;
	loopinfo.np[PHMS2_PLCLI].value = ls->sph2_mean;
	loopinfo.np[SNR2_PLCLI].value = ls->snr2_mean;
	loopinfo.np[PCLAM2_PLCLI].value = lam2_PC;

	loopinfo.np[NOMICLAM_PLCLI].value = lam_NOMIC;
	loopinfo.np[NFJUMPS_PLCLI].value = ls->nf_jumps;
	loopinfo.np[PERIOD_PLCLI].value = plc_period/1000.0;	// want seconds
	loopinfo.np[JD_PLCLI].value = secs2JD(ls->t_start);

	loopinfo.s = IPS_OK;
	IDSetNumber (&loopinfo, NULL);
}

/* callback to query for new EXPTIME from camera.
 */
static void
cameraLoopPeriodCB(void *not_used)
{
	char ynot[1024];
	int sockfd;
	FILE *camera_fp;

	// fresh connection
	sockfd = openSocket (camera_host, camera_port, false, ynot);
	if (sockfd < 0)
	    bye ("Can not connect to EXPTIME camera: %s", ynot);

	// prepare for line-oriented i/o to client
	camera_fp = fdopen (sockfd, "r+");
	if (!camera_fp)
	    bye ("Could not create FILE* from socket %d to %s: %s\n",
		    sockfd, camera_host, strerror(errno));
	setbuf (camera_fp, NULL);          // immediate writes
	// IDLog ("Successfully connected to %s:%d for camera info\n", camera_host, camera_port);

	// ask for exposure time
	fprintf (camera_fp, "EXPTIME?\n");

	// get response
	char buf[1024];
	if (ferror(camera_fp) || feof(camera_fp) || !fgets (buf, sizeof(buf), camera_fp)) {
	    if (feof(camera_fp))
		IDLog ("Connection closed to %s\n", camera_host);
	    else
		IDLog ("Connection error to %s: %s\n", camera_host, strerror(errno));
	    goto repeat;
	}

	// update plc_period if ok
	double new_period;
	if (sscanf (buf, "EXPTIME %lf", &new_period) == 1) {		// seconds
	    int new_plc_period = new_period*1000;			// secs to ms
	    if (new_plc_period != plc_period) {
		const double MIN_EXP = 0.005;
		if (new_plc_period < MIN_EXP) {
		    IDLog ("Received EXPTIME from %s, %d ms, but ignoring because < %g\n",
					    camera_host, new_plc_period, MIN_EXP);
		} else {
		    plc_period = new_plc_period;
		    IDLog ("Received new EXPTIME from %s, setting plc_period to %d ms\n",
				camera_host, plc_period);
		}
	    }
	} else {
	    IDLog ("Bogus response from %s: %s\n", camera_host, buf);
	    goto repeat;
	}

    repeat:

	// again
	fclose (camera_fp);	// also closes sockfd

	IEAddTimer (camera_period, cameraLoopPeriodCB, NULL);
}

/* called every 1000/SENDFITS ms to send PCamImage.
 * send the images the real-time loop is copying into bg_*_fim if available
 */
static void
sendFITSCB(void *not_used)
{
	static int n_loops;

	/* wait for stable bg_*_fim copy */
	pthread_mutex_lock (&plc_bglock);

	/* publish each image */
	publishImage (&bg_cam_fim, &pcamimage, IMG_PCI);
	publishImage (&bg_pad_fim, &pcambeamimage, IMG_PCBI);
	publishImage (&bg_fftmag_fim, &pcamfftmimage, IMG_PCMI);
	publishImage (&bg_fftphs_fim, &pcamfftpimage, IMG_PCPI);

	/* publish the rtdata */
	publishRTData();

	pthread_mutex_unlock (&plc_bglock);

	/* repeat.
	 * TODO: better is to sleep less by the time it took to do this work
	 */
	n_loops++;
	sf_tid = IEAddTimer ((int)(1000.0/SENDFITS), sendFITSCB, NULL);

}

/* publish the given image to the given BLOB @ eleidx in Z compressed format.
 */
static void
publishImage (FImage *fip, IBLOBVectorProperty *bvp, int eleidx)
{
	char ynot[1024];		/* error message */
	char *blob;			/* fip as malloced blob */
	int blobl;			/* size of blob */
        unsigned char *comp;            /* malloced z compressed blob */
	unsigned long ncomp;            /* size of z compressed blob */

	/* convert fip to blob */
	blobl = writeFITSmem (fip, &blob, ynot, 1);
	if (blobl < 0)
	    bye ("Can not make BLOB to publish image: %s\n", ynot);

	/* z compress */
	comp = (unsigned char *) malloc (ncomp = 1.01*(blobl+12));  /* see comments compress.c */
	int c = compress (comp, &ncomp, (unsigned char *)blob, blobl);
	if (c != Z_OK)
	    bye ("Compress returned %d\n", c);

	/* publish */
	strcpy (bvp->bp[eleidx].format, ".fits.z");
	bvp->bp[eleidx].blob = comp;
	bvp->bp[eleidx].bloblen = ncomp;
	bvp->bp[eleidx].size = blobl;
	bvp->s = IPS_OK;
	IDSetBLOB (bvp, NULL);

	/* clean up */
	free (comp);
	free (blob);
}

/* copy the real-time images from the currently selected beam to their bg_*_fim counterparts.
 * N.B. we assume plc_bglock has already been acquired.
 */
static void
copy2BG (void)
{
	OneBeam *bp = SHOWBEAM1 ? &beam1 : &beam2;

	resetFImage (&bg_cam_fim);
	cloneFImage (&bg_cam_fim, bp->raw_fimp, 1);
	codeSaturatedPixels (&bg_cam_fim);

	resetFImage (&bg_pad_fim);
	cloneFImage (&bg_pad_fim, &bp->pad_fim, 1);

	resetFImage (&bg_fftmag_fim);
	cloneFImage (&bg_fftmag_fim, &bp->fftmag_fim, 1);

	resetFImage (&bg_fftphs_fim);
	cloneFImage (&bg_fftphs_fim, &bp->fftphs_fim, 1);
}

/* publish a FITS header
 */
static void
publishHeader (FImage *fip)
{
	double v;

        if (!getRealFITS (fip, jd_kw, &v))        pcamheader.np[PCJD_PCH].value = v;
        if (!getRealFITS (fip, "PCB1X", &v))      pcamheader.np[PCB1X_PCH].value = v;
        if (!getRealFITS (fip, "PCB1Y", &v))      pcamheader.np[PCB1Y_PCH].value = v;
        if (!getRealFITS (fip, "PCBRAD", &v))     pcamheader.np[PCBRAD_PCH].value = v;
        if (!getRealFITS (fip, "PCCLOSED", &v))   pcamheader.np[PCCLOSED_PCH].value = v;
        if (!getRealFITS (fip, "PCFTOK", &v))     pcamheader.np[PCFTOK_PCH].value = v;
        if (!getRealFITS (fip, "PCLOOPHZ", &v))   pcamheader.np[PCLOOPHZ_PCH].value = v;
        if (!getRealFITS (fip, "PCPLSP1", &v))    pcamheader.np[PCPLSP1_PCH].value = v;
        if (!getRealFITS (fip, "PCPLSP01", &v))   pcamheader.np[PCPLSP01_PCH].value = v;
        if (!getRealFITS (fip, "PCPLPG", &v))     pcamheader.np[PCPLPG_PCH].value = v;
        if (!getRealFITS (fip, "PCPLDG", &v))     pcamheader.np[PCPLDG_PCH].value = v;
        if (!getRealFITS (fip, "PCPLIG", &v))     pcamheader.np[PCPLIG_PCH].value = v;
        if (!getRealFITS (fip, "PCTIPSP", &v))    pcamheader.np[PCTIPSP_PCH].value = v;
        if (!getRealFITS (fip, "PCTLTSP", &v))    pcamheader.np[PCTLTSP_PCH].value = v;

	IDSetNumber (&pcamheader, NULL);

}

/* publish rtdata_t/y to pcamdata.
 * N.B. we assume we are called while plc_bglock is acquired
 */
static void
publishRTData()
{
	static float blob[2*MAXRTSAMPS];	/* room for MAXRTSAMPS age/value pairs */
	static double last_pub;			/* time we last published */

	/* don't run until at least one entry is in ring buffers */
	if (rtdata_n <= 0)
	    return;

	/* don't publish any faster than rtpubper
	 */
	double jd0 = jdNow();
	if ((jd0 - last_pub)*SPD < rtpubper)
	    return;
	last_pub = jd0;

	/* collect all samples from the newest back to no more than rtdata_dt seconds old.
	 * samples are stored as age (in seconds) and value pairs, sorted newest first.
	 * since rtdata_i is the next index to be used, we start at (rtdata_i - 1) because that
	 *   is the index of the last entry added (N.B. beware front edge)
	 */

	int nsamples;
	for (nsamples = 0; nsamples < rtdata_n; nsamples++) {
	    int rti = ((rtdata_i - 1) - nsamples + MAXRTSAMPS) % MAXRTSAMPS;
	    double age = (jd0-rtdata_t[rti])*SPD;
	    if (age > rtdata_dt)
		break;
	    blob[2*nsamples] = (float)age;
	    blob[2*nsamples+1] = (float)rtdata_y[rti];
	}
	int blob_size = nsamples * 2 * sizeof(blob[0]);

	// IDLog ("Sending %d samples = %d bytes\n", nsamples, blob_size);
	// IDLog ("  %10.6f %g\n", blob[0], blob[1]);
	// IDLog ("  %10.6f %g\n", blob[2], blob[3]);
	// IDLog ("  %10.6f %g\n", blob[4], blob[5]);

	/* publish blob, if any */
	if (nsamples > 0) {
	    pcamdata.bp[DATA_PCD].blob = blob;
	    pcamdata.bp[DATA_PCD].bloblen = blob_size;
	    pcamdata.bp[DATA_PCD].size = blob_size;	// no compression
	    pcamdata.s = IPS_OK;
	    IDSetBLOB (&pcamdata, NULL);
	}

}


/* set all pixels in fip larger than sat_pix to MAXCAMPIX
 */
static void
codeSaturatedPixels (FImage *fip)
{
	int i;

	for (i = 0; i < fip->totbytes; i += 2) {
	    CamPix p = *(CamPix*)&fip->image[i];
	    if (p >= sat_pix)
		*(CamPix*)&fip->image[i] = MAXCAMPIX;
	}

}

/* fill path[] with full pathname for new image at jd based roughly on ISO 8601.
 * prepare any intermediate dirs as required.
 * return filename.
 */
static void
jdFilename (double jd, char *path)
{
	double jdsecs;
	time_t sec, usec;
	struct tm *tmp;
	int pl;

	/* time from jd */
	jdsecs = (jd - 2440587.5)*SPD;
	sec = (time_t)floor(jdsecs);
	usec = (time_t)floor((jdsecs - (double)sec)*1e6);

	/* break down into components */
	tmp = gmtime (&sec);

	/* insure dir structure SAVEFDIR/UTDate/hour/min */
	pl = 0;
	pl += sprintf (path+pl, "%s", SAVEFDIR);
	pl += strftime (path+pl, 100, "/PLC_%Y_%m_%d", tmp);
	(void) mkdir (path, 0755);
	pl += strftime (path+pl, 100, "/PLC_%Y_%m_%dT%H", tmp);
	(void) mkdir (path, 0755);
	pl += strftime (path+pl, 100, "/PLC_%Y_%m_%dT%H_%M", tmp);
	(void) mkdir (path, 0755);

	/* add file name */
	pl += strftime (path+pl, 100, "/PLC_%Y_%m_%dT%H_%M_%S", tmp);
	sprintf (path+pl, "_%06d.fits", (int)usec);	/* time to 1 usec */
}

/* read next *.fts file found in current directory into *fip.
 */
static void
readNextLocalFITS (FImage *fip)
{
	static int seqn;
	static int skipnfts;;
	struct dirent *dep;
	char name[sizeof(dep->d_name)];
	char msg[1024];
	char *blob;
	FILE *fp;
	DIR *dp;
	int i, n, nblob;

	/* scan for next FITS file */
    RD:
	i = 0;
	dp = opendir (".");
	while ((dep = readdir(dp))) {
	    if (strstr(dep->d_name, ".fts") || strstr(dep->d_name, ".fits")) {
		if (i++ == skipnfts) {
		    strcpy (name, dep->d_name);
		    skipnfts++;
		    break;
		}
	    }
	}
	closedir (dp);
	if (i == 0)
	    bye ("No test image .fts or .fits files found\n");
	if (!dep) {
	    skipnfts = 0;
	    goto RD;
	}

	/* read into blob */
	fp = fopen (name, "r");
	nblob = 0;
	blob = (char *) malloc(4096);		/* malloc seed */
	while ((n = fread (blob+nblob, 1, 4096, fp)) > 0)
	    blob = (char *) realloc (blob, (nblob+=n)+4096);
	fclose (fp);

	/* convert to FImage */
	resetFImage (fip);
	readFITSmem (blob, nblob, fip, msg);

	/* finished with blob */
	free (blob);

	/* add fake sequence number */
	addRealFITS (fip, "PCSEQN", seqn++, 6, "Phasecam Serial number");
}

/* like read(2) but with timeout.
 * return count, or 0 if EOF or no more, or -1 if fatal error
 */
static int
readTO (int fd, char buf[], int maxlen, int msto)
{
	struct timeval tv;
	fd_set fds;
	int r;

	FD_ZERO (&fds);
	FD_SET (fd, &fds);
	tv.tv_sec = msto/1000;
	tv.tv_usec = 1000*(msto%1000);

	r = select (fd+1, &fds, NULL, NULL, &tv);
	if (r <= 0)
	    return (r);
	return (read (fd, buf, maxlen));
}

/* rather like read(2) but dedicated to phasecam image socket on pcam_host:pcam_port.
 * always begin by opening the socket if necessary. then
 * if blocking:
 *   read exactly max bytes unless there is a socket issue (error or EOF). if there is, close
 *   the socket and return -1 such that we may be called again to try more. if there are no
 *   socket issues, all max bytes are in buf[] and we return 0 to indicate success.
 * if !blocking:
 *   keep reading from pcam until the socket would block, using buf[max] again and again if
 *   necessary and thus discarding the data read. return 0 if the next read will block, or -1
 *   if the socket had an error or EOF, in which case we close it so we can be called again.
 * only call bye() if all attempts to open seem hopeless.
 */
static int
readPCam (char buf[], int max, int blocking)
{
	static int pcam_fd;

	/* first insure socket is open */
	if (pcam_fd <= 0) {
	    static const int MAXPCTRIES = 10;
	    char ynot[1024];
	    int i;

	    for (i = 0; i < MAXPCTRIES; i++) {
		pcam_fd = openSocket (pcam_host, pcam_port, true, ynot);
		if (pcam_fd < 0)
		    sleep (1);
		else
		    break;      // got it
	    }
	    if (i == MAXPCTRIES)
		bye ("Failed to open phase cam socket %s:%d after %d tries: %s\n",
				pcam_host, pcam_port, MAXPCTRIES, ynot);
	}

	/* now do the appropriate read */
	if (blocking) {
	    /* blocking read of exactly max bytes */
	    int nread, ntot;

	    for (ntot = 0; ntot < max; ntot += nread) {
		nread = read (pcam_fd, buf+ntot, max-ntot);
		if (nread <= 0) {
		    if (nread < 0)
			IDLog ("Phasecam image socket %s:%d blocking read: %s\n",
					pcam_host, pcam_port, strerror(errno));
		    else
			IDLog ("EOF from Phasecam image socket at %s:%d on fd %d\n",
					pcam_host, pcam_port, pcam_fd);
		    close (pcam_fd);
		    pcam_fd = 0;
		    return (-1);
		}
	    }

	    noJitter();
	    noJitter();
	    noJitter();
	    noJitter();

	    return (0);

	} else {
	    /* non-blocking reads until next would block */
            int nread;
 
            do {
                nread = readTO (pcam_fd, buf, max, 0);
            } while (nread > 0);
 
            if (nread == 0)
                return (0);

	    // must be an error
	    IDLog ("Phasecam image socket %s:%d non-blocking read: %s\n",
	    				pcam_host, pcam_port, strerror(errno));
	    close (pcam_fd);
	    pcam_fd = 0;
	    return (-1);
	}
}

/* read one image from the phase cam into fip.
 * we silently reopen the camera connection here unless it looks really hopeless. We are so kind because
 *   it is normal for phasecam to close image socket whenever it is told to change configurations so we
 *   want to recover from that gracefully.
 * exit if trouble
 */
static void
pcamReadOneImage (FImage *fip)
{
	int nrows, ncols, seqn, nimbytes, r;
	int nhdr = FITS_HROWS*FITS_HCOLS;
	char buf[40960];

    /* come here when socket is closed out from under us and we need to start over
     */
    restart:

	/* look for "SOF" + 3 4-byte ints */
	if (trace)
	    IDLog ("Looking for SOF\n");
	char header[15];		// "SOF" + 3 ints
	do {
	    if (readPCam (header, sizeof(header), 1) < 0)
		goto restart;
	} while (strncmp (header, "SOF", 3));

	/* extract size from header -- any delays from odd alignment still way less than more reads()
	 * N.B. we assume little-endian
	 */
	nrows = *(int*)&header[3];
	ncols = *(int*)&header[7];
	seqn = *(int*)&header[11];

	/* read next whole fresh image */
	nimbytes = nrows * ncols * 2;
	if (nrows <= 0 || nrows > 256 || ncols <= 0 || ncols > 256) {
	    bye ("Preposterous image size reported from phase camera: %d x %d\n", nrows, ncols);
	} else {
	    char image[nhdr + nimbytes], *fits;

	    /* start with FITS header block */
	    fits = image;
	    fits += sprintf (fits, "SIMPLE  = %20s%-50s", "T",   "");
	    fits += sprintf (fits, "BITPIX  = %20d%-50s", 16,    " / bits/pix");
	    fits += sprintf (fits, "NAXIS   = %20d%-50s", 2,     " / n image axes");
	    fits += sprintf (fits, "NAXIS1  = %20d%-50s", ncols, " / n columns");
	    fits += sprintf (fits, "NAXIS2  = %20d%-50s", nrows, " / n rows");
	    fits += sprintf (fits, "BZERO   = %20d%-50s", 32768, " / physical zero");
	    fits += sprintf (fits, "BSCALE  = %20d%-50s", 1,     " / physical scale");
	    fits += sprintf (fits, "PCSEQN  = %20d%-50s", seqn,  " / Phasecam Serial number");
	    fits += sprintf (fits, "%-79s", "END");
	    while (fits < image + nhdr)
		*fits++ = ' ';

	    /* read pixels */
	    if (trace)
		IDLog ("Reading new phase cam image\n");
	    if (readPCam (image + nhdr, nimbytes, 1) < 0)
		goto restart;

	    /* convert camera's unsigned to signed for readFITSmem().
	     * convert intel to FITS byte order.
	     */
	    for (r = 0; r < nimbytes; r += 2) {
		short *sp = (short *)&image[nhdr+r];
		*sp = (short)((int)(*(unsigned short *)sp) - 32768);

		char t = image[nhdr+r];
		image[nhdr+r] = image[nhdr+r+1];
		image[nhdr+r+1] = t;
	    }

	    /* convert in-memory image to FImage format */
	    if (readFITSmem (image, nhdr+nimbytes, fip, buf) < 0)
		bye ("%s\n", buf);
	}
	noJitter();
}

/* prepare the FFT and working images b1/fftmag/fftphs_fim to be PADIMGSZ,
 */
static void
prepWorking (OneBeam *bp)
{
	// avoid needless rework
	if (PADIMGSZ != bp->padsize) {
	    /* setup FFT */
	    setupFFTW (bp, PADIMGSZ, PADIMGSZ);

	    /* create working images */
	    setupFITS (&bp->pad_fim, PADIMGSZ, PADIMGSZ);
	    setupFITS (&bp->fftmag_fim, PADIMGSZ, PADIMGSZ);
	    setupFITS (&bp->fftphs_fim, PADIMGSZ, PADIMGSZ);

	    /* remember */
	    bp->padsize = PADIMGSZ;
	}
}

/* Send a formatted command to the given socket. If trouble, brief reason is in buf[]
 *    and return -1.
 * Then, if msto > 0,
 *   wait up to msto ms for reply. If get a reply ok, it is in buf[] with trailing
 *     whitespace removed but with an EOS added and return 0.
 *   if no reply put reason why not into buf[] and return -1.
 * Return 0 if all ok.
 */
static int
socketMsg (FILE *socket, char buf[], int buflen, int msto, const char *fmt, ...)
{
	va_list va;
	int l;

	/* send */
	va_start (va, fmt);
	l = vfprintf (socket, fmt, va);
	va_end(va);
	if (l < 0) {
	    strcpy (buf, strerror(errno));
	    return (-1);
	}
	fflush (socket);

	/* receive if given a timeout */
	if (msto > 0) {
	    l = readTO (fileno(socket), buf, buflen-1, msto);
	    if (l <= 0) {
		if (l == 0)
		    strcpy (buf, "EOF");
		else
		    strcpy (buf, strerror(errno));
		return (-1);
	    }
	    cleanStr (buf, l);
	}

	/* ok */
	return (0);
}

/* create a dark reference image.
 * if running closed loop, use the next image, else take one here.
 * if ok: return 0 with base file name in dark_fn[] and dark in dark_fim.
 * if err: return -1 with short reason in whynot[].
 */
static int
takeDark (char whynot[])
{
	/* get bg_cam_fim */
	/* wait for stable bg_cam_fim copy */
	pthread_mutex_lock (&plc_bglock);

	/* reset dark */
	resetFImage (&dark_fim);

	/* store locally */
	cloneFImage (&dark_fim, &bg_cam_fim, 1);

	/* done with bg_cam_fim, ok to let real time clobber it again */
	pthread_mutex_unlock (&plc_bglock);

	/* apply any desired additional correcttions */
	if (APPLYBADPIX)
	    applyBadPixMap (&dark_fim);
	if (APPLYFLAT)
	    applyFlat (&dark_fim);

	/* mark it */
	setStringFITS (&dark_fim, "PCISDARK", "YES", "This image is used as a PhaseCam dark");
	addFITSNow (&dark_fim,    "PCDARKTM", "Time this dark frame was taken");

	/* save it */
	return (saveDark (&dark_fim, whynot));
}

/* apply dark_fim to fip.
 * record name in fip.
 */
static void
applyDark (FImage *fip)
{
	if (fip->sw != dark_fim.sw || fip->sh != dark_fim.sh) {
	    applydark.sp[APPLY_AD].s = ISS_OFF;
	    IDSetSwitch (&applydark, "Darks disabled -- wrong size: want %d x %d but found %d x %d",
						    fip->sw, fip->sh, dark_fim.sw, dark_fim.sh);
	    return;
	}

	CamPix *dp = (CamPix*)dark_fim.image;
	CamPix *ip = (CamPix*)fip->image;
	CamPix *ep = ip + fip->sw*fip->sh;
	while (ip < ep) {
	    CamPix p = *ip;
	    CamPix d = *dp++;
	    if (d < p)
		p -= d;
	    else
		p = 0;
	    *ip++ = p;
	}

	setStringFITS (fip, "PCDARK", dark_fn, "Name of Dark frame used");
}

/* look through darks_dir for a suitable dark for bg_cam_fim.
 * return 0 if ok, else -1 with short reason in whynot[].
 */
static int
findDark(char whynot[])
{
	FImage tim;
	struct dirent *ents;
	struct stat s;
	int imgw, imgh;
	int flw, flh;
	char path[1024], newpath[1024];
	long newest;
	DIR *dir;
	int fd;

	/* wait for stable bg_cam_fim copy */
	pthread_mutex_lock (&plc_bglock);
	getNAXIS (&bg_cam_fim, &flw, &flh, whynot);
	pthread_mutex_unlock (&plc_bglock);

	/* find newest dark with matching size */
	dir = opendir (darks_dir);
	if (!dir) {
	    sprintf (whynot, "%s: %s", darks_dir, strerror(errno));
	    return (-1);
	}
	newest = 0;
	DARKUNDEF();
	initFImage (&tim);
	while ((ents = readdir(dir))) {
	    if (strncmp (ents->d_name, "PLC_DARK", 8) == 0) {
		sprintf (path, "%s/%s", darks_dir, ents->d_name);
		if (stat (path, &s) == 0 && s.st_mtime > newest) {
		    fd = openFITS (path, &tim, whynot);
		    if (fd >= 0) {
			getNAXIS (&tim, &imgw, &imgh, whynot);
			if (imgw == flw && imgh == flh) {
			    newest = s.st_mtime;
			    strcpy (newpath, path);
			    strcpy (dark_fn, ents->d_name);
			}
			close (fd);
			resetFImage (&tim);
		    }
		}
	    }
	}
	closedir (dir);

	/* return if none found */
	if (!DARKISDEF()) {
	    sprintf (whynot, "No matching size");
	    return(-1);
	}

	/* read into dark_fim */
	resetFImage (&dark_fim);
	fd = openFITS (newpath, &dark_fim, whynot);
	if (fd < 0) {
	    DARKUNDEF();
	    return(-1);
	}
	close (fd);

	return (0);
}

/* write the given FITS file to be used as a dark to the darks_dir directory.
 * N.B. we NEVER restore pixel byte order after writing so fip is NO LONGER VALID on return.
 * if ok: return 0 with file name in dark_fn[].
 * if err: return -1 with short reason in whynot[].
 */
static int
saveDark (FImage *fip, char whynot[])
{
	char msg[1024];
	int fd;

	darkFileName(whynot);
	fd = open (whynot, O_RDWR|O_CREAT, S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH);
	if (fd < 0) {
	    sprintf (whynot, "open(%s): %s", dark_fn, strerror(errno));
	    return (-1);
	}
	if (writeFITS (fd, fip, msg, 0) < 0) {
	    close (fd);
	    sprintf (whynot, "write(%s): %s", dark_fn, msg);
	    return (-1);
	}
	close (fd);

	/* ok */
	return (0);
}


/* prepare a full dark file name based on current date and time in ISO 8601.
 */
static void
darkFileName (char *path)
{
	struct timeval tv;
	time_t t0;
	struct tm *tmp;

	gettimeofday (&tv, NULL);
	t0 = tv.tv_sec;
	tmp = gmtime (&t0);
	strftime (dark_fn, 100, "PLC_DARK_%Y_%m_%dT%H_%M_%S.fits", tmp);
	sprintf (path, "%s/%s", darks_dir, dark_fn);
}

/* create a flat reference image.
 * if ok: return 0 with base file name in flat_fn[] and flat in flat_fim.
 * if err: return -1 with short reason in whynot[].
 */
static int
takeFlat (char whynot[])
{
	FImage fframes[nflatimg];
	CamPix pixcore[nflatimg];
	CamPix *pptr, *eptr;
	double sum, mean;
	int seqn = 0, last_seqn;
	int npix, i, j;

	/* take each flat frame */
	for (i = 0; i < nflatimg; i++) {
	    FImage *fp = &fframes[i];

	    /* init next frame */
	    initFImage (fp);

	    /* we make sure each new frame we get from the plc is unique by
	     * insuring it has a different PCSEQN
	     */
	    do {
		/* wait for stable bg_cam_fim copy */
		pthread_mutex_lock (&plc_bglock);

		/* store in next position */
		resetFImage (fp);
		cloneFImage (fp, &bg_cam_fim, 1);

		/* done with bg_cam_fim, ok to let real time clobber it again */
		pthread_mutex_unlock (&plc_bglock);

		/* get this frame's sequence number */
		last_seqn = seqn;
		getIntFITS (fp, "PCSEQN", &seqn);

	    } while (last_seqn == seqn);

	    // report progress
	    if ((i%10) == 0 || i == nflatimg-1) {
		takeflat.s = IPS_BUSY;
		IDSetSwitch (&takeflat, "Took flat frame %d of %d", i, nflatimg);
	    }

	    /* apply any additional desired corrections */
	    if (APPLYBADPIX)
		applyBadPixMap (fp);
	    if (APPLYDARK)
		applyDark (fp);
	}

	/* reset flat with a representative frame */
	resetFImage (&flat_fim);
	cloneFImage (&flat_fim, &fframes[0], 0);

	/* set each pixel in flat_fim to median of those pixels in each frame */
	npix = flat_fim.sw * flat_fim.sh;
	for (i = 0; i < npix; i++) {
	    for (j = 0; j < nflatimg; j++)
		pixcore[j] = *(CamPix*)(&fframes[j].image[2*i]);
	    *(CamPix*)(&flat_fim.image[2*i]) = cmedian(pixcore, nflatimg);
	}

	/* finished with temp images */
	for (i = 0; i < nflatimg; i++)
	    resetFImage (&fframes[i]);

	/* mark it */
	setStringFITS (&flat_fim, "PCISFLAT", "YES", "This image is used as a PhaseCam flat");
	addFITSNow (&flat_fim, "PCFLATTM", "Time this flat frame was taken");

	/* record mean */
	sum = 0;
	pptr = (CamPix *) flat_fim.image;
	eptr = pptr + npix;
	while (pptr < eptr)
	    sum += *pptr++;
	mean = sum/npix;
	addRealFITS (&flat_fim, "PCFLMEAN", mean, 6, "Flat mean pixel value");
	IDLog ("Flat mean = %g\n", mean);

	/* save it */
	return (saveFlat (&flat_fim, whynot));
}

static void
addFITSNow (FImage *fip, const char *key, const char *comment)
{
	struct timeval tv;
	struct tm *tmp;
	time_t t0;
	char str[128];

	gettimeofday (&tv, NULL);
	t0 = tv.tv_sec;
	tmp = gmtime (&t0);
	strftime (str, sizeof(str), "%Y_%m_%dT%H_%M_%S", tmp);
	addStringFITS (fip, key, str, comment);
}


/* apply flat_fim to fip.
 * record name in fip.
 */
static void
applyFlat (FImage *fip)
{
	double mean;

	/* assert */
	if (fip->sw != flat_fim.sw || fip->sh != flat_fim.sh) {
	    applyflat.sp[APPLY_AF].s = ISS_OFF;
	    IDSetSwitch (&applyflat, "Flats disabled -- wrong size: want %d x %d have %d x %d",
			fip->sw, fip->sh, flat_fim.sw, flat_fim.sh);
	    return;
	}

	/* get pre-computed mean */
	if (getRealFITS (&flat_fim, "PCFLMEAN", &mean) < 0) {
	    applyflat.sp[APPLY_AF].s = ISS_OFF;
	    IDSetSwitch (&applyflat, "Flats disabled -- flat has no PCFLMEAN fields");
	    return;
	}

	/* apply */
	CamPix *fp = (CamPix*)flat_fim.image;
	CamPix *ip = (CamPix*)fip->image;
	CamPix *ep = ip + fip->sw*fip->sh;
	while (ip < ep) {
	    CamPix f = *fp++;
	    double newpix = (f == 0) ? MAXCAMPIX : *ip*mean/f;
	    if (newpix < 0)
		*ip++ = 0;
	    else if (newpix > MAXCAMPIX)
		*ip++ = MAXCAMPIX;
	    else
		*ip++ = (CamPix)newpix;
	}

	/* record change */
	setStringFITS (fip, "PCFLAT", flat_fn, "Name of Flat frame used");
}

/* look through flats_dir for a suitable flat for bg_cam_fim.
 * return 0 if ok, else -1 with short reason in whynot[].
 */
static int
findFlat(char whynot[])
{
	FImage tim;
	struct dirent *ents;
	struct stat s;
	int imgw, imgh;
	int flw, flh;
	char path[1024], newpath[1024];
	long newest;
	DIR *dir;
	int fd;

	/* wait for stable bg_cam_fim copy */
	pthread_mutex_lock (&plc_bglock);
	getNAXIS (&bg_cam_fim, &flw, &flh, whynot);
	pthread_mutex_unlock (&plc_bglock);

	/* find newest flat with matching size */
	dir = opendir (flats_dir);
	if (!dir) {
	    sprintf (whynot, "%s: %s", flats_dir, strerror(errno));
	    return (-1);
	}
	newest = 0;
	FLATUNDEF();
	initFImage (&tim);
	while ((ents = readdir(dir))) {
	    if (strncmp (ents->d_name, "PLC_FLAT", 8) == 0) {
		sprintf (path, "%s/%s", flats_dir, ents->d_name);
		if (stat (path, &s) == 0 && s.st_mtime > newest) {
		    fd = openFITS (path, &tim, whynot);
		    if (fd >= 0) {
			getNAXIS (&tim, &imgw, &imgh, whynot);
			if (imgw == flw && imgh == flh) {
			    newest = s.st_mtime;
			    strcpy (newpath, path);
			    strcpy (flat_fn, ents->d_name);
			}
			close (fd);
			resetFImage (&tim);
		    }
		}
	    }
	}
	closedir (dir);

	/* return if none found */
	if (!FLATISDEF()) {
	    sprintf (whynot, "No matching size");
	    return(-1);
	}

	/* read into flat_fim */
	resetFImage (&flat_fim);
	fd = openFITS (newpath, &flat_fim, whynot);
	if (fd < 0) {
	    FLATUNDEF();
	    return(-1);
	}
	close (fd);

	return (0);
}

/* write the given FITS file to be used as a flat to the flats_dir directory.
 * N.B. we NEVER restore pixel byte order after writing so fip is NO LONGER VALID on return.
 * if ok: return 0 with file name in flat_fn[].
 * if err: return -1 with short reason in whynot[].
 */
static int
saveFlat (FImage *fip, char whynot[])
{
	char msg[1024];
	int fd;

	flatFileName(whynot);
	fd = open (whynot, O_RDWR|O_CREAT, S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH);
	if (fd < 0) {
	    sprintf (whynot, "open(%s): %s", flat_fn, strerror(errno));
	    return (-1);
	}
	if (writeFITS (fd, fip, msg, 0) < 0) {
	    close (fd);
	    sprintf (whynot, "write(%s): %s", flat_fn, msg);
	    return (-1);
	}
	close (fd);

	/* ok */
	return (0);
}


/* prepare a full flat file name based on current date and time in ISO 8601.
 */
static void
flatFileName (char *path)
{
	struct timeval tv;
	time_t t0;
	struct tm *tmp;

	gettimeofday (&tv, NULL);
	t0 = tv.tv_sec;
	tmp = gmtime (&t0);
	strftime (flat_fn, 100, "PLC_FLAT_%Y_%m_%dT%H_%M_%S.fits", tmp);
	sprintf (path, "%s/%s", flats_dir, flat_fn);
}

/* create a bad pixel map image.
 * the image contains 1 for each bad pixel and 0 for each good pixel.
 * if ok: return 0 with base file name in badpix_fn[] and the map file in badpix_fim.
 * if err: return -1 with short reason in whynot[].
 */
static int
takeBadPixMap (char whynot[])
{
	time_t t0;
	int seqn = 0, last_seqn;
	int w, h;
	int npix, nimages, i;
	int dt, last_dt = 0;

	/* wait for stable bg_cam_fim copy */
	pthread_mutex_lock (&plc_bglock);
	getNAXIS (&bg_cam_fim, &w, &h, whynot);
	pthread_mutex_unlock (&plc_bglock);
	npix = w * h;

	/* accumulate stats for each pixel over time */
	double sum[npix], sum2[npix];
	memset (sum, 0, npix*sizeof(double));
	memset (sum2, 0, npix*sizeof(double));

	/* local temp image */
	FImage tim;
	initFImage (&tim);

	/* take the calibration images */
	t0 = time (NULL);
	for (nimages = 0; (dt = time(NULL) - t0) < badmapsecs; nimages++) {

	    /* we make sure each new frame we get from the plc is unique by
	     * insuring it has a different PCSEQN
	     */
	    do {
		/* wait for stable bg_cam_fim copy */
		pthread_mutex_lock (&plc_bglock);

		/* store locally */
		resetFImage (&tim);
		cloneFImage (&tim, &bg_cam_fim, 1);

		/* done with bg_cam_fim, ok to let real time clobber it again */
		pthread_mutex_unlock (&plc_bglock);

		/* get this frame's sequence number */
		last_seqn = seqn;
		getIntFITS (&tim, "PCSEQN", &seqn);

	    } while (last_seqn == seqn);

	    /* report countdown */
	    if (dt > last_dt) {
		takebadpix.s = IPS_BUSY;
		IDSetSwitch (&takebadpix, "Bad pixel seconds remaining: %d", badmapsecs - dt);
		last_dt = dt;
	    }

	    /* apply any additional desired corrections */
	    if (APPLYDARK)
		applyDark (&tim);
	    if (APPLYFLAT)
		applyFlat (&tim);

	    /* accumulate pixel stats */
	    for (i = 0; i < npix; i++) {
		CamPix p = *(CamPix*)(&tim.image[2*i]);
		sum[i] += p;
		sum2[i] += (double)p*p;
	    }
	}

	/* init badpix_fim */
	resetFImage (&badpix_fim);
	cloneFImage (&badpix_fim, &tim, 0);

	/* find pixels that change too much over time or are markedly different from neighbors */
	int ntimebad = 0, nspacebad = 0;
	for (i = 0; i < npix; i++) {
	    int bad = 0;

	    // check for time variation
	    double std = sqrt(nimages*sum2[i] - sum[i]*sum[i])/nimages;
	    if (std > bad_time_std) {
		bad = 1;
		ntimebad++;
	    }

	    // check for spatial variation, ignoring border */
	    int r = i/w;
	    int c = i%w;

	    if (r > 0 && r < h && c > 0 && c < w) {
		double *sp = &sum[i];
		double s[8];		// surrounding 8 pixels
		s[0] =                  sp[-w-1];
		insertIncreasing (s, 1, sp[-w  ]);
		insertIncreasing (s, 2, sp[-w+1]);
		insertIncreasing (s, 3, sp[  -1]);
		insertIncreasing (s, 4, sp[  +1]);
		insertIncreasing (s, 5, sp[+w-1]);
		insertIncreasing (s, 6, sp[+w  ]);
		insertIncreasing (s, 7, sp[+w+1]);

		double median = (s[3] + s[4])/2;
		double limit = bad_space_std*(s[5] - s[2])/1.349;	// convert quartile to 1 sigma
		if (sp[0] > median + limit || sp[0] < median - limit) {
		    bad = 1;
		    nspacebad++;
		}
	    }

	    *(CamPix*)(&badpix_fim.image[2*i]) = bad;
	}

	/* mark it */
	addStringFITS (&badpix_fim, "PCISBADM", "YES", "This image is used as a PhaseCam bad pixel map");
	addFITSNow (&badpix_fim,    "PCBADPTM", "Time this bad pixel map was taken");
	addIntFITS (&badpix_fim,    "PCBNFRAM", nimages, "N frames used to build bad pixel map");
	addIntFITS (&badpix_fim,    "PCBNSECS", badmapsecs, "N seconds used to build bad pixel map");
	addIntFITS (&badpix_fim,    "PCBNTIME", ntimebad, "N bad pixels over time");
	addIntFITS (&badpix_fim,    "PCBNSPAC", nspacebad, "N bad pixels over image");
	addRealFITS (&badpix_fim,   "PCBADPS",  bad_space_std, 6, "Bad spacial pixel std dev cutoff");
	addRealFITS (&badpix_fim,   "PCBADPT",  bad_time_std, 6, "Bad temporal pixel std dev cutoff");

	/* save it */
	return (saveBadPixMap (&badpix_fim, whynot));
}


/* apply badpix_fim to fip.
 * record name in fip.
 */
static void
applyBadPixMap (FImage *fip)
{
	/* assert */
	if (fip->sw != badpix_fim.sw || fip->sh != badpix_fim.sh) {
	    applybadpix.sp[APPLY_AB].s = ISS_OFF;
	    IDSetSwitch (&applybadpix, "Bad pixels disabled -- wrong size: want %d x %d have %d x %d",
	    		fip->sw, fip->sh, badpix_fim.sw, badpix_fim.sh);
	    return;
	}

	/* replace each bad pixel with mean of 2 pixels below ---
	 * looking only ahead avoids multiple applications
	 */
	int w = fip->sw;
	CamPix *bp = (CamPix*)badpix_fim.image;
	CamPix *ip = (CamPix*)fip->image;
	CamPix *ep = ip + w*(fip->sh-2);		// skip last 2 rows
	while (ip < ep) {
	    if (*bp++ > 0)
		*ip = (CamPix)(0.5*((double)ip[w] + (double)ip[2*w]));
	    ip++;
	}

	/* record change */
	setStringFITS (fip, "PCBADP", badpix_fn, "Name of bad pixel map used");
}

/* look through badpix_dir for a suitable file for bg_cam_fim.
 * return 0 if ok, else -1 with short reason in whynot[].
 */
static int
findBadPixMap(char whynot[])
{
	FImage tim;
	struct dirent *ents;
	struct stat s;
	int imgw, imgh;
	int flw, flh;
	char path[1024], newpath[1024];
	long newest;
	DIR *dir;
	int fd;

	/* wait for stable bg_cam_fim copy */
	pthread_mutex_lock (&plc_bglock);
	getNAXIS (&bg_cam_fim, &flw, &flh, whynot);
	pthread_mutex_unlock (&plc_bglock);

	/* find newest bad pixel map with matching size */
	dir = opendir (badpix_dir);
	if (!dir) {
	    sprintf (whynot, "%s: %s", badpix_dir, strerror(errno));
	    return (-1);
	}
	newest = 0;
	BADPIXUNDEF();
	initFImage (&tim);
	while ((ents = readdir(dir))) {
	    if (strncmp (ents->d_name, "PLC_BADPIX", 10) == 0) {
		sprintf (path, "%s/%s", badpix_dir, ents->d_name);
		if (stat (path, &s) == 0 && s.st_mtime > newest) {
		    fd = openFITS (path, &tim, whynot);
		    if (fd >= 0) {
			getNAXIS (&tim, &imgw, &imgh, whynot);
			if (imgw == flw && imgh == flh) {
			    newest = s.st_mtime;
			    strcpy (newpath, path);
			    strcpy (badpix_fn, ents->d_name);
			}
			close (fd);
			resetFImage (&tim);
		    }
		}
	    }
	}
	closedir (dir);

	/* return if none found */
	if (!BADPIXISDEF()) {
	    sprintf (whynot, "No matching size");
	    return(-1);
	}

	/* read into badpix_fim */
	resetFImage (&badpix_fim);
	fd = openFITS (newpath, &badpix_fim, whynot);
	if (fd < 0) {
	    BADPIXUNDEF();
	    return(-1);
	}
	close (fd);

	return (0);
}

/* write the given FITS file to be used as a badpix map to the badpix_dir directory.
 * N.B. we NEVER restore pixel byte order after writing so fip is NO LONGER VALID on return.
 * if ok: return 0 with file name in badpix_fn[].
 * if err: return -1 with short reason in whynot[].
 */
static int
saveBadPixMap (FImage *fip, char whynot[])
{
	char msg[1024];
	int fd;

	badpixFileName(whynot);
	fd = open (whynot, O_RDWR|O_CREAT, S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH);
	if (fd < 0) {
	    sprintf (whynot, "open(%s): %s", badpix_fn, strerror(errno));
	    return (-1);
	}
	if (writeFITS (fd, fip, msg, 0) < 0) {
	    close (fd);
	    sprintf (whynot, "write(%s): %s", badpix_fn, msg);
	    return (-1);
	}
	close (fd);

	/* ok */
	return (0);
}

/* prepare a full badpixel file name based on current date and time in ISO 8601.
 */
static void
badpixFileName (char *path)
{
	struct timeval tv;
	time_t t0;
	struct tm *tmp;

	gettimeofday (&tv, NULL);
	t0 = tv.tv_sec;
	tmp = gmtime (&t0);
	strftime (badpix_fn, 100, "PLC_BADPIX_%Y_%m_%dT%H_%M_%S.fits", tmp);
	sprintf (path, "%s/%s", badpix_dir, badpix_fn);
}

/* open the named FITS file into fip.
 * return 0 if ok, else -1 with short reason in whynot[]
 */
static int
openFITS (char name[], FImage *fip, char whynot[])
{
	int fd;

	fd = open (name, O_RDONLY);
	if (fd < 0) {
	    sprintf (whynot, "%s: %s", name, strerror(errno));
	    return(-1);
	}
	if (readFITS (fd, fip, whynot) < 0) {
	    close (fd);
	    return (-1);
	}
	return (fd);
}

/* insert newc into existing cp[n] in increasing order.
 * N.B. cp must in fact be of size cp[n+1] to accommodate the new entry.
 */
static void
insertIncreasing (double *cp, int n, double newc)
{
	int i;

	/* find first element larger than newc */
	for (i = 0; i < n && cp[i] < newc; i++)
	    continue;

	/* move larger values to make room to insert newc at cp[i] */
	memmove (&cp[i+1], &cp[i], (n-i)*sizeof(*cp));

	/* insert newc */
	cp[i] = newc;
}

static void
uptimeCB (void *not_used)
{
	static time_t t0;

	if (t0 == 0)
	    t0 = time (NULL);

	plcversion.np[UPTIME_PVER].value = time(NULL) - t0;
	IDSetNumber (&plcversion, NULL);

	IEAddTimer (5000, uptimeCB, NULL);
}

/* take nautobeam frames and accumulate into fresh ring buffer stats.
 */
static void
accumulateAutoBeamFrames()
{

        /* reset ring buffer stats */
        autobeam_x_rb->reset();
        autobeam_y_rb->reset();

        FImage abeamim, *fip = &abeamim;
        initFImage (fip);

        /* compute stats for nautobeam images
         */
	int seqn = 0, last_seqn;
        for (int i = 0; i < nautobeam; i++) {

	    /* make sure each new frame we get from the camera is unique by
	     * insuring it has a different PCSEQN
	     */
	    do {
		/* read next image */
                resetFImage (fip);
		readNextImage (fip);

		/* get this frame's sequence number */
		last_seqn = seqn;
		getIntFITS (fip, "PCSEQN", &seqn);

	    } while (last_seqn == seqn);

            /* apply desired corrections 
             */
	    if (APPLYBADPIX)
		applyBadPixMap (fip);
	    if (APPLYDARK)
		applyDark (fip);
	    if (APPLYFLAT)
		applyFlat (fip);

            /* find beam and accumulate stats
             */
            double x, y;
            int minp, maxp;
            findBeam (fip, &minp, &maxp, &x, &y);
            autobeam_x_rb->add (x);
            autobeam_y_rb->add (y);

            findbeam.s = IPS_BUSY;
            IDSetNumber (&findbeam, "Find Beam image %3d of %3d used pixels from %5d .. %5d found beam at [%7.2f %7.2f]",
                                i, nautobeam, minp, maxp, x, y);
        }

        resetFImage (fip);
}

/* in fip find center of what looks like the light beam in the quadrant of the BEAM* values.
 * x and y values returned are coordinates within fip of beam center, ie NOT FITS. 
 * also report range of pixel values used for centroid, just FYI.
 */
static void
findBeam (FImage *fip, int *minpix, int *maxpix, double *center_x, double *center_y)
{
        CamPix *pix = (CamPix *)fip->image;             // handy

	/* determine quadrant containing beam.
	 * N.B. BEAMX/Y are FITS 1-based
	 */
	int qx, qy, qw, qh;
	qw = fip->sw/2;
	qh = fip->sh/2;
	qx = ((BEAM1X - 1) < qw) ? 0 : qw;
	qy = ((BEAM1Y - 1) < qh) ? 0 : qh;

	/* find histogram in quadrant
	 */
	int hist[NCAMPIX];
        memset (hist, 0, sizeof(hist));
	for (int y = 0; y < qh; y++) {
	    for (int x = 0; x < qw; x++) {
		int p = (int)pix[(y+qy)*fip->sw + (x+qx)];
                hist[p]++;
	    }
	}

#define BEAM_BRIGHEST_N 100
#ifdef BEAM_BRIGHEST_N

	/* find range from simply the brightest N pixels
	 */
	int dim_pix = MAXCAMPIX, bright_pix = MAXCAMPIX;
	for (int n = 0; dim_pix >= 0 && (n += hist[dim_pix]) < BEAM_BRIGHEST_N; dim_pix--)
	    continue;

#else

	/* find range from peaks and lows in histogram
	 */

        const int NBINS = 10;                // number of histogram bins

        /* find lo and hi hist values
         */
        int lo, hi;
        for (lo = 0; lo < NCAMPIX; lo++)
            if (hist[lo] > 0)
                break;
        for (hi = NCAMPIX; --hi >= 0; )
            if (hist[hi] > 0)
                break;

        /* collect this range into NBINS
         */
        int bins[NBINS];
        memset (bins, 0, sizeof(bins));
        int range = hi-lo;
        for (int i = 0; i < range; i++)
            bins[i*NBINS/range] += hist[lo + i];

        /* left side is at first valley -- should be the break between background and everything else
         */
        int lbini = 0;
        for (int i = 1; i < NBINS-2; i++) {     // never want lbini to be right-most bin 
            if (bins[i-1] > bins[i] && bins[i] < bins[i+1]) {
                lbini = i;
                break;
            }
        }

        int dim_pix = lo + lbini*range/NBINS + range/NBINS/2;

        /* continue right until find the next lower bin
         */
        int rbini = lbini + 1;
        for (int i = rbini; i < NBINS-1; i++) {
            if (bins[i] < bins[lbini]) {
                rbini = i;
                break;
            }
        }

        int bright_pix = lo + rbini*range/NBINS + range/NBINS/2;

	/* uncomment to log the histogram if desired
	for (int i = lo; i <= hi; i++)
	    IDLog ("HIST %d %d %g %d\n", i, hist[i], log10(bins[(i-lo)*NBINS/range]),
                        i >= dim_pix && i <= bright_pix);
        for (int i = 0; i < NBINS; i++)
            IDLog ("BINHIST %d %d\n", i, bins[i]);
	 */

#endif // !BEAM_BRIGHEST_N

	/* find centroid within quadrant using only pixels between dim_pix and bright_pix
	 */
	double qsum = 0, qxsum = 0, qysum = 0;
	for (int y = 0; y < qh; y++) {
	    for (int x = 0; x < qw; x++) {
		int p = (int)pix[(y+qy)*fip->sw + (x+qx)];
		if (p >= dim_pix && p <= bright_pix) {
		    qsum += p;
		    qxsum += p * x;
		    qysum += p * y;
		}
	    }
	}

        *minpix = dim_pix;
        *maxpix = bright_pix;

	*center_x = qxsum/qsum + qx;
	*center_y = qysum/qsum + qy;
}

/* timer callback to move SPC to offload error from FPC.
 * do nothing if currently off
 */
static void
spcOffloadCB (void *not_used)
{
	static char indi_host[] = "lbti-web";
	static FILE *indifp;

	/* insure connection to local INDI server in order to send motor command message.
	 * N.B. this assumes all motors are available via lbti-web
	 */
	if (!indifp) {
	    char ynot[1024];
	    int indifd = openSocket (indi_host, 7624, false, ynot);
	    if (indifd < 0) {
		spcoffload.s = IPS_ALERT;
		spcoffload.sp[YES_SPCOL].s = ISS_OFF;
		IDSetSwitch (&spcoffload, "Offload failed to connect to %s: %s", indi_host, ynot);
		return;
	    }
	    indifp = fdopen (indifd, "w");
	}

	if (closeloop.sp[YES_CLP].s == ISS_ON && SPCOFFLOAD) {
	    if (fabs(Piston_Mean) > spc_offload_max) {
		// move translation stage to bring piston back to center (0)
		double steps = spc_offload_scale * Piston_Mean;
		fprintf (indifp, "<newTextVector device='%s' name='%s'>\n",
					spc_offload_devnam, spc_offload_motnam);
		fprintf (indifp, "  <oneText name='command'>%s%.0f</oneText>\n", 
					steps < 0 ? "&lt;" : "&gt;", fabs(steps));	// use the entity!
		fprintf (indifp, "</newTextVector>\n");
		fflush (indifp);
		if (ferror(indifp) || feof(indifp)) {
		    fclose (indifp);
		    indifp = NULL;
		    spcoffload.s = IPS_ALERT;
		    spcoffload.sp[YES_SPCOL].s = ISS_OFF;
		    IDSetSwitch (&spcoffload, "Offload failed: lost connection to %s", indi_host);
		    return;
		}

		// record
		spc_offload_net += steps/spc_offload_scale;
		IDLog ("Piston_Mean at %g: Moved %s.%s by %g. Net offload now %g microns\n",
			Piston_Mean, spc_offload_devnam, spc_offload_motnam, steps/spc_offload_scale,
			spc_offload_net);
	    }
	}

	// again 
	spc_offload_tid = IEAddTimer (spc_offload_ms, spcOffloadCB, NULL);
}

/* remove trailing whitespace from buf[] and add EOS.
 */
static void
cleanStr (char buf[], int l)
{
	while (l > 0 && isspace(buf[l-1]))
	    l--;
	buf[l] = '\0';
}

/* return JD now
 */
static double
jdNow()
{
	return (secs2JD(secsNow()));
}

/* return seconds since the epoch
 */
static double
secsNow()
{
	struct timeval tv;
	gettimeofday (&tv, NULL);
	return (tv.tv_sec+tv.tv_usec/1e6);
}

static double
secs2JD (double secs)
{
	return (2440587.5 + secs/SPD);
}

unsigned bswap32( unsigned  a_int )
{
	return  ( ( a_int << 24 ) & 0xff000000 ) |
	( ( a_int <<  8 ) & 0x00ff0000 ) |
	( ( a_int >>  8 ) & 0x0000ff00 ) |
	( ( a_int >> 24 ) & 0x000000ff );
}

void bswap64( long long int *a_int )
{
	unsigned long long l32, h32;

	h32 = bswap32( ( unsigned int )( *a_int & 0x00000000ffffffffULL ) );
	l32 = bswap32( ( unsigned int )( ( *a_int >> 32 ) & 0x00000000ffffffffULL ) );
	*a_int = ( h32 << 32 ) | l32;
}

/* thread to monitor OVMS channels forever.
 */
static void *
OVMSThread (void *dummy)
{
    start:
	IDLog ("Starting OVMS thread\n");

	/* subscribe to multicast */

	int asock = socket (AF_INET, SOCK_DGRAM, 0);
	int reuse = 1;
	setsockopt (asock, SOL_SOCKET, SO_REUSEADDR, (char *)&reuse, sizeof(reuse));
	struct timeval tv;
	tv.tv_sec = 2;
	tv.tv_usec = 500000;
	setsockopt (asock, SOL_SOCKET, SO_RCVTIMEO, (char *)&tv, sizeof(tv));
	struct sockaddr_in localSock;
	memset (&localSock, 0, sizeof(localSock));
	localSock.sin_family = AF_INET;
	localSock.sin_port = htons (ovms.port);
	localSock.sin_addr.s_addr = INADDR_ANY;
	if (bind (asock, (struct sockaddr *)&localSock, sizeof(localSock)) < 0)
	    bye ("Can not bind to port %d: %s\n", ovms.port, strerror(errno));
	struct ip_mreq group;
	group.imr_multiaddr.s_addr = inet_addr(ovms.group);
	group.imr_interface.s_addr = inet_addr(ovms.interface);
	if (setsockopt (asock, IPPROTO_IP, IP_ADD_MEMBERSHIP, (char *) &group, sizeof(group)) < 0 ) {
	    publishCLShutdown ("Can not join multicast group %s: %s\n", ovms.group, strerror(errno));
	    sleep (10);
	    goto start;
	}

	/* filled in from the latest status packet */
	int numDevices = 0, numChannels = 0;
	int estimationMode = -1, estimationSize = 0, horizon = 0;

	/* read forever, restart if any trouble */
	IDLog ("OVMS thread running\n");

	/* history OPD motions */
	double OPD_prev = 0;
	bool OPD_prev_set = false;

	while (1) {

	    // handy
	    double now = secsNow();

	    // read next packet
	    double databuf[1024];
	    int nread = read (asock, databuf, sizeof(databuf));
	    if (nread < 0) {
		publishCLShutdown ("OVMS socket error: %s\n", strerror(errno));
		close (asock);
		sleep (10);
		goto start;
	    }

	    // extract status or data packet, depending on sign of time stamp (!)
	    bswap64( (long long int*) &databuf[0]);
	    double timestamp = databuf[0];
	    if (timestamp < 0) {

		// byte-swap each status packet entry 
		for (int i = 1; i < 10; i++)
		    bswap64 ((long long int*) &databuf[i]);

		// extract fields -- see OVMS Development Reference Manual, Table 1
		numDevices = (int) databuf[1];
		numChannels = (int) databuf[2];
		estimationMode = (int) databuf[8];
		estimationSize = (int) databuf[9];
		horizon = databuf[10];

	    } else {

		// values to be determined from data packet
		double OPD = 0;
		double TIP = 0;
		double TILT = 0;
		bool ok;

		/* crack depending on mode, if known.
		 * see OVMS User Manual: 509g502e-UserManual.pdf
		 */
		switch (estimationMode) {
		case 1:	{ // use 3 elements beyond channel data

			int start = 1 + numDevices*numChannels;
			for (int i = start; i < start+3; i++)
			    bswap64 ((long long int*) &databuf[i]);

			TIP  = databuf[start+0];
			TILT = databuf[start+1];
			OPD  = databuf[start+2];

			ok = true;
		    }
		    break;

		case 2:	{ // use 5 elements beyond channel data

			int start = 1 + numDevices*numChannels;
			for (int i = start; i < start+5; i++)
			    bswap64 ((long long int*) &databuf[i]);

			double TIP_SX  = databuf[start+0];
			double TILT_SX = databuf[start+1];

			double TIP_DX  = databuf[start+2];
			double TILT_DX = databuf[start+3];

			TIP  = TIP_DX - TIP_SX;
			TILT = TILT_DX - TILT_SX;
			OPD  = databuf[start+4];

			ok = true;
		    }
		    break;

		case 3:	{ // use 18 elements beyond channel data

			int start = 1 + numDevices*numChannels;
			for (int i = start; i < start+18; i++)
			    bswap64 ((long long int*) &databuf[i]);

			double TIP_M1_SX  = databuf[start+0];
			double TILT_M1_SX = databuf[start+1];
			double OPD_M1_SX  = databuf[start+2];
			double TIP_M2_SX  = databuf[start+3];
			double TILT_M2_SX = databuf[start+4];
			double OPD_M2_SX  = databuf[start+5];
			double TIP_M3_SX  = databuf[start+6];
			double TILT_M3_SX = databuf[start+7];
			double OPD_M3_SX  = databuf[start+8];

			double TIP_M1_DX  = databuf[start+9];
			double TILT_M1_DX = databuf[start+10];
			double OPD_M1_DX  = databuf[start+11];
			double TIP_M2_DX  = databuf[start+12];
			double TILT_M2_DX = databuf[start+13];
			double OPD_M2_DX  = databuf[start+14];
			double TIP_M3_DX  = databuf[start+15];
			double TILT_M3_DX = databuf[start+16];
			double OPD_M3_DX  = databuf[start+17];

			TIP  = (-2*TIP_M1_DX  + 2*TIP_M2_DX  - 1.41*TIP_M3_DX)  -
			       (-2*TIP_M1_SX  + 2*TIP_M2_SX  - 1.41*TIP_M3_SX);
			TILT = (-2*TILT_M1_DX + 2*TILT_M2_DX - 1.41*TILT_M3_DX) -
			       (-2*TILT_M1_SX + 2*TILT_M2_SX - 1.41*TILT_M3_SX);
			OPD  = (-2*OPD_M1_DX  + 2*OPD_M2_DX  - 1.41*OPD_M3_DX)  -
			       (-2*OPD_M1_SX  + 2*OPD_M2_SX  - 1.41*OPD_M3_SX);
			
			ok = true;
		    }
		    break;

		default: { // unknown
			ok = false;
		    }
		    break;
		}

		// install into ovms, protected by lock
		pthread_mutex_lock (&ovms.lock);

		    ovms.timestamp = timestamp;
		    ovms.opd = OPD;
		    ovms.tip = TIP;
		    ovms.tilt = TILT;
		    ovms.horizon = horizon;
		    ovms.ok = ok;

		    // detect and record possible bump
		    if (OPD_prev_set && fabs(OPD-OPD_prev) > opd_bump)
			ovms.bump_t0 = now;
		    ovms.bump = now < ovms.bump_t0 + opd_freeze;

		pthread_mutex_unlock (&ovms.lock);

		OPD_prev = OPD;
		OPD_prev_set = true;

	    }

	}

	bye ("OVMS thread unexpectedly quit\n");

	return (NULL);			// lint
}


/* initialize the path length biquad filter.
 * N.B. only valid while loop is running to get good estimate for loop rate.
 */
static void
setupVibFilterBiQuad()
{
	/* gain is value of peak when surrounding floor is 1, in dB
	 */
        vib_bqfilter.definePeakingFilter (VIB_FC, loopinfo.np[RATE_PLCLI].value, VIB_Q, VIB_FGAIN);
        vib_bqfilter.initRun();
}

static void
startOVMSThread(void)
{
        pthread_attr_t attr;
        pthread_t thr;

	/* init lock */
	pthread_mutex_init (&ovms.lock, NULL);

        /* new thread will be detached so we need never join */
        (void) pthread_attr_init (&attr);
        (void) pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);

        /* start the thread that monitors accelerometers forever */
        if (pthread_create (&thr, &attr, OVMSThread, NULL))
            bye ("Can not start OVMS thread: %s\n", strerror(errno));

}

/* start beam2Thread()
 */
static void
startBeam2Thread(void)
{
        pthread_attr_t attr;
        pthread_t thr;

	/* init locks and cvs */
	(void) pthread_mutex_init (&b2_run_lock, NULL);
	(void) pthread_cond_init (&b2_run_cv, NULL);
	(void) pthread_mutex_init (&b2_done_lock, NULL);
	(void) pthread_cond_init (&b2_done_cv, NULL);
	b2_run = false;
	b2_done = false;

        /* new thread will be detached so we need never join */
        (void) pthread_attr_init (&attr);
        (void) pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);

	/* crank up priority */
	if (pthread_attr_setschedpolicy (&attr, SCHED_RR))
	    IDLog ("INFO: Can not set scheduling policy for beam2 to SCHED_RR\n");
	else {
	    IDLog ("INFO: Successfully set scheduling policy for beam2 to SCHED_RR\n");
	    struct sched_param shp;
	    int max = sched_get_priority_max(SCHED_RR);
	    shp.sched_priority = max;
	    if (pthread_attr_setschedparam (&attr, &shp) != 0)
		IDLog ("INFO: Can not set scheduling priority for beam2 to %d\n", max);
	    else
		IDLog ("INFO: Successfully set scheduling priority for beam2 to %d\n", max);
	}

        /* start the thread that computes beam2 stats forever */
        if (pthread_create (&thr, &attr, beam2Thread, NULL))
            bye ("Can not start beam2Thread: %s\n", strerror(errno));

	/* stay on cpu_beam2 */
	cpu_set_t cpuset;
	CPU_ZERO (&cpuset);
	CPU_SET (cpu_beam2, &cpuset);
	if (pthread_setaffinity_np(thr, sizeof(cpu_set_t), &cpuset))
	    IDLog ("INFO: Can not lock beam2 thread to cpu %d\n", cpu_beam2);
	else
	    IDLog ("INFO: Successfully locked beam2 thread to cpu %d\n", cpu_beam2);
}

/* publish why the closed-loop thread is being shut down
 */
static void
publishCLShutdown (const char *fmt, ...)
{
	va_list ap;
	char buf[1024];
	int l;
	
	l = sprintf (buf, "Loop is now open: ");

	va_start (ap, fmt);
	vsprintf (buf+l, fmt, ap);
	va_end (ap);

	closeloop.sp[YES_CLP].s = ISS_OFF;
	closeloop.s = IPS_ALERT;
	IDSetSwitch (&closeloop, buf);

	plcOpenLoop ();
}

/* insure the given directort exist, creating if necessary.
 * exit if fail to create
 */
static void
insureDir (char *dir)
{
	struct stat s;

	if (stat (dir, &s) < 0) {
	    if (mkpath (dir, 0777) < 0)
		bye ("Can not create %s: %s\n", dir, strerror(errno));
	    else
		IDLog ("Created %s\n", dir);
	}
}

/* recursively create path, if necessary.
 * return 0 if ok else -1 with errno
 */
static int
mkpath (char *path, int mode)
{
	static const char slash = '/';
	static const char eos = '\0';
        char *p = path;
        int ok, atend, atslash;

        do {
            atend = (*p == eos);
            atslash = (*p == slash && p > path);
            ok = 1;
            if (atend || atslash) {
                *p = eos;
                ok = (mkdir (path, mode) == 0 || errno == EEXIST);
                if (atslash)
                    *p = slash;
            }
            p++;
        } while (ok && !atend);

        return (ok ? 0 : -1);
}

/* define beam_pad from beam_pad0 and possibly BEAMR.
 * if beam_pad0 is 0, set beam_pad automatically based on BEAMR, else just use beam_pad.
 *   based on known-good case where BEAMR==7 want padding to be 4 to get 60x60 FFT
 */
static void
setBeamPad()
{
    if (beam_pad0)
	beam_pad = beam_pad0;
    else
	beam_pad = (int) fmax (1.0, 2.0*(floor(30.0/BEAMW+0.5)));
    IDLog ("BEAMR = %d beam_pad = %d\n", BEAMR, beam_pad);
}

/* this function seems to do nothing useful but sprinkling it around in strategic locations
 * greatly reduces the real-time thread jitter.
 */
static void
noJitter()
{
	fclose (fopen ("/dev/null", "w"));
}

/* print count and ms since init call, ms since previous call, then the given message.
 * we add \n if none in fmt.
 * N.B. requies #define	TIMETRACE
 */
// #define	TIMETRACE
static void
timeTrace (bool init, const char *fmt, ...)
{
#if defined(TIMETRACE)
	static struct timeval tv_init, tv_last;
	static int n_init;

	if (!tv_init.tv_sec || init) {
	    gettimeofday (&tv_init, NULL);
	    tv_last = tv_init;
	    n_init = 0;
	}

	struct timeval tv;
	gettimeofday (&tv, NULL);
	double ms_init = 1e3*(tv.tv_sec - tv_init.tv_sec) + (1e-3*tv.tv_usec - 1e-3*tv_init.tv_usec);
	double ms_last = 1e3*(tv.tv_sec - tv_last.tv_sec) + (1e-3*tv.tv_usec - 1e-3*tv_last.tv_usec);

	fprintf (stderr, "%3d %7.3f %7.3f : ", n_init, ms_init, ms_last);

	va_list ap;
	va_start (ap, fmt);
	vfprintf (stderr, fmt, ap);
	va_end (ap);

	if (!strchr (fmt, '\n'))
	    fprintf (stderr, "\n");

	tv_last = tv;
	n_init++;
	if (init)
	    tv_init = tv;
#endif // TIMETRACE
}

/* log the given message, sleep for EXITSLEEP, then exit(1).
 * also try to inform user by sending on a Property. We choose closeloop because 
 * probably when we die and loose loop is the worst case situation.
 * the delay prevents us from being restarted so quickly.
 */
static void
bye (const char *fmt, ...)
{
	char msg[2048];
	int l;
	va_list va;

	/* format message */
	va_start (va, fmt);
	l = vsprintf (msg, fmt, va);
	va_end (va);

	/* strip any trailing nl for INDI msg */
	if (msg[l-1] == '\n')
	    msg[l-1] = '\0';

	/* send to any clients */
	closeloop.sp[YES_CLP].s = ISS_OFF;
	closeloop.s = IPS_ALERT;
	IDSetSwitch (&closeloop, msg);

	/* log */
	IDLog ("%s\n", msg);

	/* sleep (to void repart restart cycling if condition persists) */
	sleep (EXITSLEEP);

	/* gone */
	exit(1);
}

/* dummy function for libip */
void
pm_set (int p)
{
}