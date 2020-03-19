#ifndef RTS_H
#define RTS_H
// authors: Christoph Klauer and Raphael Hartmann

//#define _CRT_SECURE_CPP_OVERLOAD_STANDARD_NAMES 1
//#include <string.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <vector>
#include <string>
#include <algorithm>
#include <R.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <cfloat> // f�r DBL_MAX

using std::vector;
using std::setw;
using std::cout;
using std::string;

// #define DATA "../exp3_for_MPT-RT.dat"
//#define DATA "../Daten/broeder_data/Experiment 3 - drift/DM3_zus.txt"
// #define DATA "../Daten/broeder_data/Experiment 2 - threshold/"
//#define DATA "../Daten/broeder_data/Experiment 1 - bias/"
// #define MODEL "../PD_C_overrides_A__guess.info"
//#define MODEL "../Daten/broeder_data/Experiment 3 - drift/2HTMDI.info"
// #define MODEL "../Daten/broeder_data/Experiment 1 - bias/2HTMDI.info"


//#define RAUS "../Daten/broeder_data/Experiment 1 - bias/samples/raus.DI_resp_tukey"
//#define RAUS "../Daten/broeder_data/Experiment 2 - threshold/results/raus.DG_resp_tukey"
//#define RAUS "../Daten/broeder_data/Experiment 3 - drift/results/raus.DI_resp_tukey"
//#define RAUS "../Daten/dube/exp2/raus_dube2_DG_2_responses_tukey"



// #define RAUS "raus"
//#define LOGLIK "log_likelihood"
// #define log_lik_flag false
// #define for_bridge_flag false
const double sqr2pi = 2.506628274631000502415765e0;

//#define diagnose
#define generate_and_save

#ifdef generate_and_save
const bool generate_or_diagnose = true;
const bool save_diagnose=true;
#endif


#ifdef diagnose
const bool generate_or_diagnose = false;
const bool save_diagnose=true;
#endif

#define DEBUG false

//#define PRIOR 1.00

//#define adf 2
// #define IREP 2200
// #define BURNIN 200 // (must be smaller than IREP)
// #define THIN 11 // (must be a divisor of IREP)
// #define NOTHREADS 4
// #define SAMPLE_SIZE 20000 //(a multiple of IREP*NOTHREADS/THIN)
// #define RMAX 1.05

#define T_rng gsl_rng_ranlxd1

#define LNNORM_MAX_X 38.0
#define LNNORM_MIN_X -1.00e9

struct trial {
	int person, tree, category, item,group,rt;
};


struct pfadinfo {
int a;
vector<int> r;
vector<int> pfad_par;
vector<int> pm;
};


struct point {
	double x;
	double h;
	double dh;
};


extern int log_lik_flag;
extern int for_bridge_flag;

extern const char *DATA;
extern const char *MODEL;
extern const char *RAUS;
extern const char *diagn_tests;
extern const char *LOGLIK;

extern int IREP;
extern int BURNIN;
extern int THIN;
extern int NOTHREADS;
extern int SAMPLE_SIZE;
extern double RMAX;

extern int nKERN;
extern int nRESP;

extern int *CatToResp;

extern double *ConstProb;
extern int *CompMinus;
extern int *CompPlus;

extern double *complete_sample;
extern double *complete_bridge;

//Globale Variablen

extern int *cat2tree;
extern int kernpar;
extern int kerncat;
extern int indi;
extern int zweig;
extern int *branch;
extern int nodemax;
extern int *ar;
extern int *nodes_per_tree;
extern int *tree_and_node2par;
extern bool *comp;
extern int ilamfree, ifree,ipred;

extern int *ndrin,*drin;
extern int n_all_parameters;
extern int *nppr;
extern int n_bridge_parameters;

extern int RMAX_reached;
extern bool BURNIN_flag;

extern int igroup;
extern int *t2group;
extern int ireps;
extern int *cat2resp;
extern int respno;
extern int alphaoff;
extern int sigalphaoff;
extern int restparsno;
extern int *free2kern;
extern int *kern2free;
extern double *consts;

extern double pr_df_sigma_sqr;
extern double pr_shape_omega_sqr;
extern double pr_rate_omega_sqr;
extern double pr_mean_mu_gamma;
extern double pr_var_mu_gamma;
extern double PRIOR;
extern double pr_shape_exp_mu_beta;
extern double pr_rate_exp_mu_beta;
extern double pr_sf_scale_matrix_SIG;
extern double pr_sf_scale_matrix_TAU;
extern int pr_df_add_inv_wish;

extern int *pfad_index;
extern vector<pfadinfo> path_info;

double lnnorm(double x);

void make_idaten(vector<trial> daten,int *idaten);

//void by_individuals(vector<trial> daten,int kerntree,double *beta, double *g2,double *likeli,gsl_rng *rst);

double truncnorm(double b, gsl_rng *rst);

double onenorm(gsl_rng *rst);

void make_pij_for_individual(double *x, double *pij, double *pj);

void make_mu(double *mu, double *lams, double *beta,  int *nnodes, double *z, gsl_rng *rst);

void make_lams(double *mu, double *lams, double *beta,  int *nnodes, double *z, gsl_rng *rst);


void bayesreg(int n,double *mean,double *sigma, double *out, gsl_rng *rst);

void invwis(int cases, int nvar, double *xx, double *ssig, double *sigi, double eps, gsl_rng *rst);

double truncexp(double lambda, double upper, gsl_rng *rst);

double equation(int t, int ip, double *mu, double *lams, double *beta);

double oneuni(gsl_rng *rst);

double ars(double step, double &scale,double totallow, double n, double xp, double *beta, double *sigi, double *lambdas,double *lams,int tt, int iz, double start, gsl_rng *rst,
	void gamma_prior(double scale, double norm, double n, double alpha, double p, double *beta, double *sigi, double *lambdas,double *lams,int t,  int iz, bool deriv, point &h));


void r_statistic(int ido, int n_all_parameters, int istream, int iter, double *parmon, double *xwbr, double &rmax);

double rexp(double x);

void diagnosis(vector<trial> daten,int *idaten, int kerntree, gsl_rng *rst);

//void make_index(bool *lambda_comp, int *index);

int make_path_for_one_trial(int branchno, double *pij, double p, gsl_rng *rst);

double oneexp(double lambda, gsl_rng *rst);

void lies(vector<trial> &daten);

void model_design(int kerntree,int *ar, int *branch, int *nodes_per_par, int *nodes_per_tree, int *tree_and_node2par);

void set_ns(vector<trial> daten, int &indi, int &kerntree, int &kerncat,int &igroup, int &ntot);

//void ddf(vector <trial> &daten, int kerntree, double *g2);

double elogdiff(double xa, double xb);

void make_parameters_for_all(double *mu, double *lams, double *beta, double *x_for_all);

void make_pij_for_one_trial(trial one,double *x_for_all,  double *pij, double &pj);

void make_zs_one_trial(trial one,int itrial, int ipath, double *mu, double *lams, double *beta,int *nz_position, double *z,gsl_rng *rst) ;

double malpha(int t, int r,double* restpars,double* slams);

void make_rtau(double *restpar, double *taui, double *slams, gsl_rng *rst);

void make_rmu(vector<trial> daten,double *factor, double *rest,double *restpar, double* slams, gsl_rng *rst);

void make_ralpha(vector<trial> daten,double *factor, double *rest, double *restpar,double *slams, double *taui, gsl_rng *rst);

void make_rsigalpha(vector<trial> daten, double* factor, double *rest, double *restpar, double *slams, bool xflag, gsl_rng *rst);

void make_rsig(vector<trial> daten, double *rest, double *restpar, gsl_rng *rst);

void make_slams(vector<trial> daten,double* factor, double *rest, double *restpars, double *slams, gsl_rng *rst);

double double_truncnorm(double lower, double upper, gsl_rng *rst);

double logsum(double xa, double xb);

double logexgaussian(double lam, double mu, double sd, double t);

void gibbs_times_new(vector<trial> daten,int *nnodes,int nz, int *nz_position,double *beta,int ntau, int *ntau_position,
	gsl_rng *rst1, gsl_rng *rst2, gsl_rng *rst3, gsl_rng *rst4, gsl_rng *rst5, gsl_rng *rst6, gsl_rng *rst7, gsl_rng *rst8, gsl_rng *rst9, gsl_rng *rst10, gsl_rng *rst11, gsl_rng *rst12, gsl_rng *rst13, gsl_rng *rst14, gsl_rng *rst15, gsl_rng *rst16,
	double *lambdas, double* restpars);

void make_tij_for_one_trial_new(trial one,double *rhos, double *lambdas,double *restpars,double *pij);

void tby_individuals(vector<trial> daten, int kerntree,double *beta,double *lambdas, double *restpars, gsl_rng *rst);

double logdiff(double xa, double xb);

void extract_pfadinfo(int *pfad_index, vector<pfadinfo> &path_info);

void loggammagaussian(int n, double lam, double mu, double sd, double t, double& hplus, double& hminus);

double  logf_tij(int a, vector<int> r, double *lams, double *loglams, double mu, double sd, double t);

int gsl_linalg_tri_lower_invert(gsl_matrix * T);


#endif
