
// authors: Raphael Hartmann and Christoph Klauer

#define R_NO_REMAP

#include "rts.h"
#include <chrono>

const char *MODEL;
const char *DATA;
int nKERN;
int nPROCS;
int nRESP;
int *CatToResp = 0;
// number of all parameters
int n_all_parameters;
// number of total trials
int datenzahl;
// loglikelihood vector
double *loglik_vec;

namespace ertmpt {

	int log_lik_flag;
	int for_bridge_flag;

	const char *RAUS;
	const char *diagn_tests;
	const char *LOGLIK;

	int NOTHREADS;
	int BURNIN;
	int THIN;
	int SAMPLE_SIZE;
	int IREP;

	double *ConstProb = 0;
	int *CompMinus = 0;
	int *CompPlus = 0;

	double RMAX;

	double *complete_sample = 0;
	double *complete_bridge = 0;

	int n_bridge_parameters;

	double pr_df_sigma_sqr;
	double pr_shape_omega_sqr;
	double pr_rate_omega_sqr;
	double pr_mean_mu_gamma;
	double pr_var_mu_gamma;
	double PRIOR;
	double pr_shape_exp_mu_beta;
	double pr_rate_exp_mu_beta;
	double pr_sf_scale_matrix_SIG;
	double pr_sf_scale_matrix_TAU;
	int pr_df_add_inv_wish;

}


namespace drtmpt {

  const char *RAUS;
  const char *LOGLIK;
  const char *CONTINUE;
  const char *MEANSOUT;
  const char *TESTSOUT;
  const char *RANDOM;
  const char *TMPDIR;

  //Settings for MCMC runs
  int IREP;
  //adpatation phase in phase == 1 and phase == 3; must be multiple of THIN
  int PHASE1;
  //estimation of mass matrix; must be multiple of IREP and should be larger or equal to PHASE1
  int PHASE2;
  //optional thinning must be divisor of IREP
  int THIN;
  //number of threads for parallelization
  int NOTHREADS;
  //SAMPLE_SIZE must be multiple of NOTHREADS * IREP/THIN; This is true because = NOTHREADS * IREP
  int SAMPLE_SIZE;
  //Maximal number of CPU Threads to be used for DIC calculation and initial starting value calculation
  int MAXTHREADS;
  //when to end phase 4 sampling
  double RMAX;

  // sample
  double *complete_sample = 0;

  // whether DIC is computed or not
  bool DIC;
  //whether trialwise log-likelihood values are to be saved when computing DIC
  bool log_lik_flag;
  //intialize with random values (=0) or with personwise maximum-likelihood values (=1)
  int INITIALIZE;

  //prior precision for population means of process-related parameters
  double PRIOR;
  //LKJ priors for the process-related params
  double etat; //shape parameter for LKJ distribution
  double taut; //scale parameter half-Cauchy sds process parameters
  //LKJ priors for the motor-related params
  double etar; //shape parameter for LKJ distribution
  double taur; //scale parameter half-Cauchy sds motor-time parameters
  //degrees of freedom t-distribution of motor times; mu and scale
  int degf;
  double mu_prior;
  double rsd;
  //Gamma-Prior Omega
  double prioralpha;
  double priorbeta;

  //goon = true: Continuation mode; otherwise normal sampling until max(r) <= rmax
  bool goon;
  //how many samples to add if continuing (goon = true);
  //must be multiple of IREP * NOTHREADS/THIN
  int ADDITION;// = 4 * 500;

  //maxtree-depth in Phases 1 - 3
  int maxtreedepth1_3;// = 5;
  //maxtree-depth in Phase 4
  int maxtreedepth4;// = 9;

  //int kernpar;
  int *kern2free = 0;
  int ifree[3];
  bool *comp = 0;
  double *consts = 0;

}



extern "C" {

	SEXP ertmpt_fit(SEXP re, SEXP re2, SEXP re3, SEXP ch, SEXP in, SEXP in2, SEXP in3, SEXP in4, SEXP in5, SEXP bo1, SEXP bo2, SEXP bo3) {

	  using namespace ertmpt;

		RMAX = REAL(re)[0];

		DATA = R_CHAR(STRING_ELT(ch, 0));
		MODEL = R_CHAR(STRING_ELT(ch, 1));
		RAUS = R_CHAR(STRING_ELT(ch, 2));
		diagn_tests = R_CHAR(STRING_ELT(ch, 3));
		LOGLIK = R_CHAR(STRING_ELT(ch, 4));

		NOTHREADS = INTEGER(in)[0];
		BURNIN = INTEGER(in)[1];
		THIN = INTEGER(in)[2];
		SAMPLE_SIZE = INTEGER(in)[3];
		IREP = INTEGER(in)[4];
		nKERN = INTEGER(in)[5];
		nPROCS = INTEGER(in)[6];
		nRESP = INTEGER(in)[7];

		CatToResp = (int *)calloc(nKERN, sizeof(int));
		ConstProb = (double *)calloc(nPROCS, sizeof(double));
		CompMinus = (int *)calloc(nPROCS, sizeof(int));
		CompPlus = (int *)calloc(nPROCS, sizeof(int));
		for (int i = 0; i < nKERN; i++) {
		  CatToResp[i] = INTEGER(in2)[i];
		}
		for (int i = 0; i < nPROCS; i++) {
		  ConstProb[i] = REAL(re2)[i];
		  CompMinus[i] = INTEGER(bo1)[i];
		  CompPlus[i] = INTEGER(bo2)[i];
		}




		log_lik_flag = INTEGER(bo3)[0];
		for_bridge_flag = INTEGER(bo3)[1];

		pr_df_sigma_sqr = REAL(re3)[0];
		pr_shape_omega_sqr = REAL(re3)[1];
		pr_rate_omega_sqr = REAL(re3)[2];
		pr_mean_mu_gamma = REAL(re3)[3];
		pr_var_mu_gamma = REAL(re3)[4];
		PRIOR = REAL(re3)[5];
		pr_shape_exp_mu_beta = REAL(re3)[6];
		pr_rate_exp_mu_beta = REAL(re3)[7];
		pr_sf_scale_matrix_SIG = REAL(re3)[8];
		pr_sf_scale_matrix_TAU = REAL(re3)[9];
		pr_df_add_inv_wish = INTEGER(in3)[0];

		int *k2f = INTEGER(in4);
		int *f2k = INTEGER(in5);


		mainx(k2f, f2k);


		int outCnt = 0, prtCnt = 0;
		SEXP prob = PROTECT(Rf_allocVector(REALSXP, 1));
		outCnt++;
		SEXP pars_samples = PROTECT(Rf_allocMatrix(REALSXP, SAMPLE_SIZE, (n_all_parameters+1)));
		outCnt++;
		SEXP pars_bridge = PROTECT(Rf_allocMatrix(REALSXP, SAMPLE_SIZE, (n_bridge_parameters+1)));
		outCnt++;
		SEXP loglik = PROTECT(Rf_allocMatrix(REALSXP, SAMPLE_SIZE, datenzahl));
		outCnt++;
		SEXP ans = PROTECT(Rf_allocVector(VECSXP, outCnt));
		prtCnt = outCnt + 1;


		double *Rprob = REAL(prob);
		double *Rpars_samples = REAL(pars_samples);
		double *Rpars_bridge = REAL(pars_bridge);
		double *Rloglik = REAL(loglik);
		

		Rprob[0] = 0.44332211;
		for (int i=0; i!=SAMPLE_SIZE; i++) {
			for (int j = 0; j!=(n_all_parameters+1); j++) {
				Rpars_samples[j*SAMPLE_SIZE + i] = complete_sample[i*(n_all_parameters+1) + j];
			}
			if (for_bridge_flag) {
				for (int j = 0; j!=(n_bridge_parameters+1); j++) {
					Rpars_bridge[j*SAMPLE_SIZE + i] = complete_bridge[i*(n_bridge_parameters+1) + j];
				}
			} else {
				for (int j = 0; j!=(n_bridge_parameters+1); j++) {
					Rpars_bridge[j*SAMPLE_SIZE + i] = 0;
				}
			}
			if (log_lik_flag) {
			  for (int j = 0; j!=datenzahl; j++) {
			    Rloglik[j*SAMPLE_SIZE + i] = loglik_vec[j*SAMPLE_SIZE + i];
			  }
			}
		}
		if (complete_sample) free(complete_sample);
		if (complete_bridge) free(complete_bridge);


		SET_VECTOR_ELT(ans,0,prob);
		SET_VECTOR_ELT(ans,1,pars_samples);
		SET_VECTOR_ELT(ans,2,pars_bridge);
		SET_VECTOR_ELT(ans,3,loglik);
		

		SEXP names = PROTECT(Rf_allocVector(STRSXP, outCnt));
		prtCnt++;
		SET_STRING_ELT(names,0,Rf_mkChar("prob"));
		SET_STRING_ELT(names,1,Rf_mkChar("pars_samples"));
		SET_STRING_ELT(names,2,Rf_mkChar("pars_bridge"));
		SET_STRING_ELT(names,3,Rf_mkChar("LogLik"));
		

		Rf_setAttrib(ans,R_NamesSymbol,names);

		// // free variables
		free(CatToResp);
		free(ConstProb);
		free(CompMinus);
		free(CompPlus);
		free(loglik_vec);


		/* Unprotect the ans and names objects */
		UNPROTECT(prtCnt);

		return(ans);
	}

}



extern "C" {

  SEXP drtmpt_fit(SEXP ch1, SEXP in1, SEXP re1, SEXP bo1, SEXP in2, SEXP re2, SEXP in3, SEXP in4, SEXP re3, SEXP in5) {

    // NAMESPACE
    using namespace drtmpt;

    // PATHS
    DATA = R_CHAR(STRING_ELT(ch1, 0));
    MODEL = R_CHAR(STRING_ELT(ch1, 1));
    RAUS = R_CHAR(STRING_ELT(ch1, 2));
    LOGLIK = R_CHAR(STRING_ELT(ch1, 3));
    CONTINUE = R_CHAR(STRING_ELT(ch1, 4));
    MEANSOUT = R_CHAR(STRING_ELT(ch1, 5));
    TESTSOUT = R_CHAR(STRING_ELT(ch1, 6));
    RANDOM = R_CHAR(STRING_ELT(ch1, 7));
    TMPDIR = R_CHAR(STRING_ELT(ch1, 8));


    // SAMPLING SETTINGS
    IREP = INTEGER(in1)[0];
    PHASE1 = INTEGER(in1)[1];
    PHASE2 = INTEGER(in1)[2];
    THIN = INTEGER(in1)[3];
    NOTHREADS = INTEGER(in1)[4];
    SAMPLE_SIZE = INTEGER(in1)[5];//2 * NOTHREADS * IREP;
    MAXTHREADS = INTEGER(in1)[6];
    nKERN = INTEGER(in1)[7];
    nPROCS = INTEGER(in1)[8];
    nRESP = INTEGER(in1)[9];

    CatToResp = (int *)calloc(nKERN, sizeof(int));
    for (int i = 0; i < nKERN; i++) {
      CatToResp[i] = INTEGER(in1)[10+i];
    }

    RMAX = REAL(re1)[0];


    // FLAGS
    DIC = INTEGER(bo1)[0];
    log_lik_flag = INTEGER(bo1)[1];
    INITIALIZE = INTEGER(bo1)[2];


    // PRIORS
    degf = INTEGER(in2)[0];

    PRIOR = REAL(re2)[0];
    etat = REAL(re2)[1];
    taut = REAL(re2)[2];
    etar = REAL(re2)[3];
    taur = REAL(re2)[4];
    mu_prior = REAL(re2)[5];
    rsd = REAL(re2)[6];
    prioralpha = REAL(re2)[7];
    priorbeta = REAL(re2)[8];


    // HMC OPTIONS
    maxtreedepth1_3 = INTEGER(in3)[0];
    maxtreedepth4 = INTEGER(in3)[1];


    // CONTINUATION
    goon = INTEGER(in4)[0];
    ADDITION = INTEGER(in4)[1];


    // CONSTANTS AND EQUALIZATION
    consts = (double*)malloc(nPROCS * 3 * sizeof(double));
    for (int i = 0; i < nPROCS*3; i++) {
      consts[i] = REAL(re3)[i];
    }
    kern2free = (int*)malloc(nPROCS * 3 * sizeof(int));
    comp = (bool*)malloc(nPROCS * 3 * sizeof(bool));
    for (int i = 0; i < nPROCS*3; i++) {
      kern2free[i] = INTEGER(in5)[i];
      comp[i] = (INTEGER(in5)[i+nPROCS*3] == 1);
      if (i < 3) ifree[i] = INTEGER(in5)[i+nPROCS*6];
    }


    main_d();


    int outCnt = 0, prtCnt = 0;
    SEXP pars_samples = PROTECT(Rf_allocMatrix(REALSXP, SAMPLE_SIZE, n_all_parameters));
    outCnt++;
    SEXP loglik = PROTECT(Rf_allocMatrix(REALSXP, SAMPLE_SIZE, datenzahl));
    outCnt++;
    SEXP ans = PROTECT(Rf_allocVector(VECSXP, outCnt));
    prtCnt = outCnt + 1;


    double *Rpars_samples = REAL(pars_samples);
    double *Rloglik = REAL(loglik);


    for (int i=0; i< SAMPLE_SIZE; i++) {
      for (int j = 0; j < n_all_parameters; j++) {
        Rpars_samples[i + j*SAMPLE_SIZE] = complete_sample[i*(n_all_parameters) + j];
      }
      if (log_lik_flag) {
        for (int j = 0; j < datenzahl; j++) {
          Rloglik[i + j*SAMPLE_SIZE] = loglik_vec[i*datenzahl + j];
        }
      }
    }


    if (complete_sample) free(complete_sample);
    free(loglik_vec);


    SET_VECTOR_ELT(ans, 0, pars_samples);
    SET_VECTOR_ELT(ans, 1, loglik);


    SEXP names = PROTECT(Rf_allocVector(STRSXP, outCnt));
    prtCnt++;
    SET_STRING_ELT(names, 0, Rf_mkChar("pars_samples"));
    SET_STRING_ELT(names, 1, Rf_mkChar("loglik"));


    Rf_setAttrib(ans, R_NamesSymbol, names);

    /* Unprotect the ans and names objects */
    UNPROTECT(prtCnt);


    // FREE DYNAMIC VARIABLES
    if (kern2free) free(kern2free);
    if (consts) free(consts);
    if (comp) free(comp);
    if (CatToResp) free(CatToResp);


    return(ans);
  }

}
