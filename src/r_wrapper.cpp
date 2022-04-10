#define R_NO_REMAP
// authors: Raphael Hartmann and Christoph Klauer
// #include <Rcpp.h>
// #include <R.h>
// #include <Rmath.h>
#include "rts.h"
#include "main.h"
#include <Rinternals.h>
#ifdef _OPENMP
#include <omp.h>
#endif


// namespace rtsNS {

	int log_lik_flag;
	int for_bridge_flag;

	const char *DATA;
	const char *MODEL;
	const char *RAUS;
	const char *diagn_tests;
	const char *LOGLIK;

	int NOTHREADS;
	int BURNIN;
	int THIN;
	int SAMPLE_SIZE;
	int IREP;

	int nKERN;
	int nRESP;
	int *CatToResp = 0;

	double *ConstProb = 0;
	int *CompMinus = 0;
	int *CompPlus = 0;

	double RMAX;

	double *complete_sample = 0;
	double *complete_bridge = 0;

	int n_all_parameters;
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

// }

extern "C" {

	SEXP rtmpt_fit(SEXP re, SEXP re2, SEXP re3, SEXP ch, SEXP in, SEXP in2, SEXP in3, SEXP in4, SEXP in5, SEXP bo1, SEXP bo2, SEXP bo3) {

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
		nRESP = INTEGER(in)[6];

		CatToResp = INTEGER(in2);

		ConstProb = REAL(re2);
		CompMinus = INTEGER(bo1);
		CompPlus = INTEGER(bo2);

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
		SEXP ans = PROTECT(Rf_allocVector(VECSXP, outCnt));
		prtCnt = outCnt + 1;


		double *Rprob = REAL(prob);
		double *Rpars_samples = REAL(pars_samples);
		double *Rpars_bridge = REAL(pars_bridge);


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
		}
		if (complete_sample) free(complete_sample);
		if (complete_bridge) free(complete_bridge);


		SET_VECTOR_ELT(ans,0,prob);
		SET_VECTOR_ELT(ans,1,pars_samples);
		SET_VECTOR_ELT(ans,2,pars_bridge);


		SEXP names = PROTECT(Rf_allocVector(STRSXP, outCnt));
		prtCnt++;
		SET_STRING_ELT(names,0,Rf_mkChar("prob"));
		SET_STRING_ELT(names,1,Rf_mkChar("pars_samples"));
		SET_STRING_ELT(names,2,Rf_mkChar("pars_bridge"));


		Rf_setAttrib(ans,R_NamesSymbol,names);

		// // free variables
		// free(CatToResp);
		// free(ConstProb);
		// free(CompMinus);
		// free(CompPlus);


		/* Unprotect the ans and names objects */
		UNPROTECT(prtCnt);

		return(ans);
	}

}
