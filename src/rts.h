#ifndef RTS_H
#define RTS_H
// authors: Christoph Klauer and Raphael Hartmann
#define _CRT_SECURE_CPP_OVERLOAD_STANDARD_NAMES 1


//#define _CRT_SECURE_CPP_OVERLOAD_STANDARD_NAMES 1
#include <algorithm>
#include <cfloat> // for DBL_MAX
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stack>
#include <string>
#include <thread>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <gsl/gsl_blas.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>

#include <R.h>
#include <Rinternals.h>

//model file path
extern const char *MODEL;
// data file path
extern const char *DATA;
//number of categories
extern int kerncat;
//total number of trials
extern int datenzahl;
//maximum number of paths per category
extern int zweig;
//number of diffusion-model parameters per person
extern int kernpar;
//maximum number of nodes per subtree
extern int nodemax;
//number of persons by group
extern int* ng;
//number of persons
extern int indi;
//index of group by person
extern int *t2group;
//number of groups
extern int igroup;
//total number of model parameters
extern int n_all_parameters;
//maps categories on responses
extern int *cat2resp;
//number of responses
extern int respno;
//loglikelihood vector
extern double *loglik_vec;
//number of process parameters
extern int nKERN;
//number of responses
extern int nRESP;
//maps categories on responses (from R)
extern int *CatToResp;


//trial information inputted
struct trial {
  int person, tree, category, item, group, rt;
};


#define T_rng gsl_rng_ranlxd1
//for lnnorm.cpp
#define LNNORM_MAX_X 38.0
#define LNNORM_MIN_X -1.00e9
#define M_LN_SQRT_PI	0.572364942924700087071713675677
#define M_PISQ 9.86960440108935861883449099987615113531
#define  accuracy -27.63102
#define sqr2pi 2.506628274631000502415765e0


namespace ertmpt {
  
  int mainx(int *k2f, int *f2k); //int argc, char *argv[]
  
  // #define RAUS "raus"
  //#define LOGLIK "log_likelihood"
  // #define log_lik_flag false
  // #define for_bridge_flag false
  
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
  

  struct pfadinfo {
  int a;
  std::vector<int> r;
  std::vector<int> pfad_par;
  std::vector<int> pm;
  };
  
  
  struct point {
  	double x;
  	double h;
  	double dh;
  };
  
  
  extern int log_lik_flag;
  extern int for_bridge_flag;
  
  extern const char *RAUS;
  extern const char *diagn_tests;
  extern const char *LOGLIK;
  
  extern int IREP;
  extern int BURNIN;
  extern int THIN;
  extern int NOTHREADS;
  extern int SAMPLE_SIZE;
  extern double RMAX;
  
  extern double *ConstProb;
  extern int *CompMinus;
  extern int *CompPlus;
  
  extern double *complete_sample;
  extern double *complete_bridge;
  
  //Globale Variablen
  
  extern int *cat2tree;
  extern int *branch;
  extern int *ar;
  extern int *nodes_per_tree;
  extern int *tree_and_node2par;
  extern bool *comp;
  extern int ilamfree, ifree,ipred;
  
  extern int *ndrin,*drin;
  extern int *nppr;
  extern int n_bridge_parameters;
  
  extern int RMAX_reached;
  extern bool BURNIN_flag;
  
  extern int ireps;
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
  extern std::vector<pfadinfo> path_info;
  
  void make_idaten(std::vector<trial> daten,int *idaten);
  
  //void by_individuals(std::vector<trial> daten,int kerntree,double *beta, double *g2,double *likeli,gsl_rng *rst);
  
  void make_pij_for_individual(double *x, double *pij, double *pj);
  
  void make_mu(double *mu, double *lams, double *beta,  int *nnodes, double *z, gsl_rng *rst);
  
  void make_lams(double *mu, double *lams, double *beta,  int *nnodes, double *z, gsl_rng *rst);
  
  
  void bayesreg(int n,double *mean,double *sigma, double *out, gsl_rng *rst);
  
  void invwis(int cases, int nvar, double *xx, double *ssig, double *sigi, double eps, gsl_rng *rst);
  
  double truncexp(double lambda, double upper, gsl_rng *rst);
  
  double equation(int t, int ip, double *mu, double *lams, double *beta);
  
  double ars(double step, double &scale,double totallow, double n, double xp, double *beta, double *sigi, double *lambdas,double *lams,int tt, int iz, double start, gsl_rng *rst,
  	void gamma_prior(double scale, double norm, double n, double alpha, double p, double *beta, double *sigi, double *lambdas,double *lams,int t,  int iz, bool deriv, point &h));
  
  
  void r_statistic(int ido, int n_all_parameters, int istream, int iter, double *parmon, double *xwbr, double &rmax);
  
  double rexp(double x);
  
  void diagnosis(std::vector<trial> daten,int *idaten, int kerntree, gsl_rng *rst);
  
  //void make_index(bool *lambda_comp, int *index);
  
  int make_path_for_one_trial(int branchno, double *pij, double p, gsl_rng *rst);
  
  double oneexp(double lambda, gsl_rng *rst);
  
  void model_design(int kerntree,int *ar, int *branch, int *nodes_per_par, int *nodes_per_tree, int *tree_and_node2par);
  
  //void ddf(std::vector <trial> &daten, int kerntree, double *g2);
  
  double elogdiff(double xa, double xb);
  
  void make_parameters_for_all(double *mu, double *lams, double *beta, double *x_for_all);
  
  void make_pij_for_one_trial(trial one,double *x_for_all,  double *pij, double &pj);
  
  void make_zs_one_trial(trial one,int itrial, int ipath, double *mu, double *lams, double *beta,int *nz_position, double *z,gsl_rng *rst) ;
  
  double malpha(int t, int r,double* restpars,double* slams);
  
  void make_rtau(double *restpar, double *taui, double *slams, gsl_rng *rst);
  
  void make_rmu(std::vector<trial> daten,double *factor, double *rest,double *restpar, double* slams, gsl_rng *rst);
  
  void make_ralpha(std::vector<trial> daten,double *factor, double *rest, double *restpar,double *slams, double *taui, gsl_rng *rst);
  
  void make_rsigalpha(std::vector<trial> daten, double* factor, double *rest, double *restpar, double *slams, bool xflag, gsl_rng *rst);
  
  void make_rsig(std::vector<trial> daten, double *rest, double *restpar, gsl_rng *rst);
  
  void make_slams(std::vector<trial> daten,double* factor, double *rest, double *restpars, double *slams, gsl_rng *rst);
  
  double double_truncnorm(double lower, double upper, gsl_rng *rst);
  
  double logexgaussian(double lam, double mu, double sd, double t);
  
  void gibbs_times_new(std::vector<trial> daten,int *nnodes,int nz, int *nz_position,double *beta,int ntau, int *ntau_position,
  	gsl_rng *rst1, gsl_rng *rst2, gsl_rng *rst3, gsl_rng *rst4, gsl_rng *rst5, gsl_rng *rst6, gsl_rng *rst7, gsl_rng *rst8, gsl_rng *rst9, gsl_rng *rst10, gsl_rng *rst11, gsl_rng *rst12, gsl_rng *rst13, gsl_rng *rst14, gsl_rng *rst15, gsl_rng *rst16,
  	double *lambdas, double* restpars);
  
  void make_tij_for_one_trial_new(trial one,double *rhos, double *lambdas,double *restpars,double *pij);
  
  void tby_individuals(std::vector<trial> daten, int kerntree,double *beta,double *lambdas, double *restpars, gsl_rng *rst);
  
  double logdiff(double xa, double xb);
  
  void extract_pfadinfo(int *pfad_index, std::vector<pfadinfo> &path_info);
  
  void loggammagaussian(int n, double lam, double mu, double sd, double t, double& hplus, double& hminus);
  
  double  logf_tij(int a, std::vector<int> r, double *lams, double *loglams, double mu, double sd, double t);
  
  int gsl_linalg_tri_lower_invert(gsl_matrix * T);

}


namespace drtmpt {
  
  extern bool PROG_BAR_FLAG;
  
  
  const bool generate_or_diagnose = true;
  const bool save_diagnose = true;
  //goon = true: Continuation mode; otherwise normal sampling until max(r) <= rmax
  extern bool goon;
  //how many samples to add if continuing (goon = true);
  //must be multiple of IREP * NOTHREADS/THIN
  extern int ADDITION;
  
  /*
  #ifdef diagnose
   const bool generate_or_diagnose = false;
   const bool save_diagnose = true;
   const bool goon = false;
   const int ADDITION = 0;
  #endif
   */
  
  // sample
  extern double *complete_sample;
  
  // output file path
  extern const char *RAUS;
  //input file path for last sample for restarting sampling
  extern const char *CONTINUE;
  //output file path for trialwise log-likelihood values
  extern const char *LOGLIK;
  extern const char *MEANSOUT;
  extern const char *TESTSOUT;
  extern const char *RANDOM;
  extern const char *TMPDIR;
  
  //whether trialwise log-likelihood values are to be saved when computing DIC
  extern bool log_lik_flag;
  // whether DIC is computed or not
  extern bool DIC;
  //CPUs to be used in computing DIC
  extern int DIC_CPUs;
  
  //intialize with random values (=0) or with personwise maximum-likelihood values (=1)
  extern int INITIALIZE;
  //CPUs to be used if INITIALIZE = 1
  extern int INIT_CPUs;
  
  
  
  //prior precision for population means of process-related parameters
  extern double PRIOR;
  
  //maxtree-depth in Phases 1 - 3
  extern int maxtreedepth1_3;
  //maxtree-depth in Phase 4
  extern int maxtreedepth4;
  
  // adaptive epsilon
  const double gamplus = 0.05;
  const double t0plus = 10.0;
  const double kapplus = 0.75;
  extern double muplus;
  
  //LKJ priors for the process-related params
  extern double etat; //shape parameter for LKJ distribution
  extern double taut; //scale parameter half-Cauchy sds process parameters
  //LKJ priors for the motor-related params
  extern double etar; //shape parameter for LKJ distribution
  extern double taur; //scale parameter half-Cauchy sds motor-time parameters
  
  //Gamma-Prior Omega
  extern double prioralpha;
  extern double priorbeta;
  
  //degrees of freedom t-distribution of motor times; mu and scale
  extern int degf;
  extern double mu_prior;
  extern double rsd;
  
  //inverse chi-square degrees of freedom of distribution of individual motor-time variance
  const double priordf = 2.0;
  
  // Huang-Wand dfs
  const int huang = 2;
  
  //Settings for MCMC runs
  extern int IREP;
  //adpatation phase in phase == 1 and phase == 3; must be multiple of THIN
  extern int PHASE1;
  //estimation of mass matrix; must be multiple of IREP and should be larger or equal to PHASE1
  extern int PHASE2;
  //optional thinning must be divisor of IREP
  extern int THIN;
  extern int NOTHREADS;
  //SAMPLE_SIZE must be multiple of NOTHREADS * IREP/THIN; This is true because = NOTHREADS * IREP
  extern int SAMPLE_SIZE;
  //when to end phase 4 sampling
  extern double RMAX;
  //Maximal number of CPU Threads to be used for DIC calculation and initial starting value calculation
  extern int MAXTHREADS;
  
  extern int RMAX_reached;
  
  //how many categories in each subtree
  #define dNKS(T,IT) nks[T*kerntree+IT]
  //branch probabilities
  #define dPIJ(I,J) pij[I*zweig+J]
  //from tree it and category j within subtree to category number
  #define dTREE2CAT(IT,J) tree2cat[IT*kerncat+J]
  
  //population means of diffusion parameters (on real line)
  #define dMAVW(IG, Type, IP) mavw[IG * 3 * ifreemax + Type * ifreemax + IP]
  //index in hampar
  #define dmapMAVW(IG, Type, IP) mapmavw[IG * 3 * ifreemax + Type * ifreemax + IP]
  
  //individual deviations of diffusion parameters (on real line)
  #define dAVW(T, Type, IP) avw[T * 3 * ifreemax + Type * ifreemax + IP]
  //index in hampar
  #define dmapAVW(T, Type, IP) mapavw[T * 3 * ifreemax + Type * ifreemax + IP]
  //transformed individual diffusion parameters after logistic transformation of model parameters
  #define dTAVW(T, Type, IP) tavw[T * 3 * ifreemax + Type * ifreemax + IP]
  
  //components for computing r statistics
  #define dXWBR(T,I) xwbr[(T-1)*n_all_parameters+I]
  #define dPARMON(I,J) parmon[(I-1)*n_all_parameters+J]
  
  //is node R in path J ending in category I (=1 outcome +; = -1 outcome -) or not = 0?
  #define dAR(I,J,R) ar[I*zweig*nodemax + J*nodemax + R]
  //how many interior nodes in path K ending in category J
  #define dNDRIN(J,K) ndrin[J*zweig+K]
  //index of node no. X in path K ending in category J
  #define dDRIN(J,K,X) drin[J*zweig*nodemax+K*nodemax + X]
  //how many interior edges (n,o) on paths to category C
  #define dNCDRIN(C) ncdrin[C]
  //index of node (PM=0) and output (PM=1) on INth edge on paths to C)
  #define dCDRIN(C,IN,PM) cdrin[C*2*(nodemax*2) +IN*2 + PM]
  
  //diffusion-model parameters attachted to nodes
  #define dTREE_AND_NODE2PAR(T,N,Type) tree_and_node2par[(T)*nodemax*3+(N)*3+Type]
  //index of combination of diffusion-model parameter attached to node
  #define dTREE_AND_NODE2MAP(T, N) tree_and_node2map[(T)*nodemax + (N)]
  
  //how many nodes with parameter-combination J for person I
  #define dNNODES(I,J) nnodes[I*no_patterns+J]
  
  //index of process-completion time of trial X, node N, outcome PM
  #define dTAU_BY_NODE(X,N,PM) tau_by_node[X*2*nodemax+N*2 + PM]
  
  //frequency data for person T by category I
  #define dIDATEN(T,I) idaten[T*kerncat + I]
  
  //derivatives of log-likelihood according to population-mean diffusion-model parameters
  #define dDHWIEN(IG,Type,IP) dhwien[(IG) * 3 * ifreemax + (Type) * ifreemax + IP]
  //derivatives for personwise deviations
  #define dDWIEN(T,Type,IP) dwien[(T) * 3 * ifreemax + (Type) * ifreemax + IP]
  
  //how many diffusion-model combinations M occur in trials for person T outside this person's current paths
  #define dNIPS(T,PM,M) nips[T*2*no_patterns  + PM*no_patterns + M]
  
  //inverse of personwise variance/covariance matrix diffusion-model parameters
  #define dSIGI(I,J) sigi[I*icompg + J]
  //inverse of personwise variance/covarince matrix motor-time times
  #define dGAMI(I,J) gami[I*respno + J]
  //personwise variance/covariance matrix diffusion-model parameters
  #define dSIG(I,J) sig[I*icompg + J]
  //personwise variance/covariance matrix motor times
  #define dGAM(I,J) gam[I*respno + J]
  
  //index of parameter (type,ipar) on scale from 1 to ifreeg; encodes parameters set equal
  #define dKERN2FREE(Type,Ipar) kern2free[Type*kernpar + Ipar]
  //index of parameter (type,ipar) on scale from 1 to icompg (=-1 for constant parameters)
  #define dFREE2COMP(Type,Ipar) free2comp[Type*kernpar + Ipar]
  
  //constant parameter values
  #define dCONSTS(Type,IP) consts[IP*3 + Type]
  //is parameter constant (=false) or not (=true)
  #define dCOMP(Type,IP) comp[IP*3 + Type]
  
  //some performance statistics
  //#define MONITOR(I,IP) monitor[I*2*10 + IP]
  
  //map parameter combination on index
  #define dMAP(IA,IV,IW) map[IA*ifree[1]*ifree[2] + IV*ifree[2] + IW]
  //index of parameter combination N's ifree numbers by type
  #define dCOMB(N,TYPE) comb[N*3 + TYPE]
  
  //store posterior samples
  #define dSAMPLE(I,IP) sample[(I)*(n_all_parameters)+IP]
  
  //posterior sample variance-covariance interim results per thread
  #define dSUPERSIG(ITHREAD,I,J) supersig[ITHREAD * n_all_parameters * n_all_parameters  + I * n_all_parameters + J]
  
  //points stored for adaptive rejection sampler
  struct point {
    double x;
    double h;
    double dh;
  };
  
  //for upper/lower slope functions adaptive rejection sampler
  struct piece {
    double z;
    double slope;
    double absc;
    double center;
  };
  
  //save adaptive rejection sampler results by parameter combination
  struct ars_archiv {
    std::vector<std::vector<point>> hstore;
    std::vector <std::vector<piece>> lowerstore;
    std::vector <std::vector<piece>> upperstore;
    std::vector<double> startstore;
    std::vector<double> scalestore;
    std::vector<double> normstore;
    std::vector<std::vector<double>> sstore;
    
  };
  
  //parameters of logistic transformation of diffusion-model parameters
  struct transform {
    double loc;
    double scale;
    double a;
    double b;
    double range;
  };
  
  //for convolution of process-completion times and motor times densities
  struct my_params {
    int pfadlength;
    double* a;
    double* v;
    double* w;
    int* low_or_up;
    double mu;
    double sig;
    double rt_rest;
  };
  
  //parameter store for NUTS sampler
  struct Theta {
    double* loglambda;
    double* tavw;
    double* tlams;
    gsl_vector* hampar;
  };
  
  //node in tree built by NUTS sampler
  struct Node {
    int status; // status = 0: left leaf; status = 1: right leaf; status = 2: interior
    int level, index;
    struct Node* left, * right;
  };
  
  //information saved at leaves of tree built by NUTS sampler
  struct store {
    int n, s;
    gsl_vector* p;
    Theta* theta, * thetadash;
    int na;
    double a;
  };
  
  
  //Global variables
  
  //Tree structures in NUTS sampler
  extern Node* trees[13];
  
  //variance/covariance matrix of posterior distribution
  extern double* supersig;
  //number of parameters in NUTS sampler
  extern int nhamil;
  extern int phase; //phase = 1,2, 3, or 4
  
  //parameters of logistic transformation per type
  extern transform avwtrans[3];
  
  //which tree is a category in
  extern int *cat2tree;
  //how many paths by category
  extern int *branch;
  //see above under define AR
  extern int *ar;
  //how many nodes per subtree
  extern int *nodes_per_tree;
  //see above under define TREE_AND_NODE2PAR
  extern int *tree_and_node2par;
  //dito
  extern int* tree_and_node2map;
  //dito
  extern bool* comp;
  //how many parameters (after accounting for equality constraints) of each type; sum thereof; maximum per type
  extern int ifree[3], ifreeg, ifreemax;
  //how many parameters (after accounting for equality constraints and constants) of each type; sum thereof
  extern int icomp[3], icompg;
  //how many trials by person
  extern int* n_per_subj;
  
  //offsets in index numbers (hampar) for individual deviations diffusion-model parameters; mean motor times; individual deviations motor times; diffusion-model deviations variance-covariance matrix-related parameters
  extern int iavwoff, irmuoff, ilamoff, isigoff;
  
  //see above under define NNODES
  extern int* nnodes;
  //see above under define; pfadmax[c] = maximum number of interior nodes on a path ending in category c
  extern int *ndrin,*drin, *cdrin, *ncdrin, *pfadmax;
  //number of trials with response r by person
  extern int *nppr;
  //see above under define TAU_BY_NODE
  extern int* tau_by_node;
  //minimum size of Gibbs-Hamiltonian cycles before interim report
  extern int ireps;
  //see above under define KERN2FREE
  extern int* kern2free;
  //see above under define FREE2COMP
  extern int *free2comp;
  //constant value by parameter; -1 for parameters to be estimated
  extern double *consts;
  //some performance statistics
  extern double *monitor;
  //see above under define MAP
  extern int* map;
  //see above under define COMB
  extern int* comb;
  //number of different diffusion-model parameter combinations
  extern int no_patterns;
  //see above under define MAPMAVW
  extern int* mapmavw;
  //see above under define MAPAVW
  extern int* mapavw;
  //minimum rt associated with person, parameter combination, and outcome
  extern std::vector<double> rtmins;
  //batch size for storing posterior samples
  extern int sample_size;
  //total number of taus
  extern int ntau;
  
  
  int main_d();
  //mass matrix hamiltonian
  extern gsl_matrix* supsig;
  //Cholesky factor of inverse of mass matrix
  extern gsl_matrix* sigisqrt;
  //log density diffusion model
  double dwiener_d(double q, double a, double v, double w, double eps);
  //subtract exp(xb) from exp(xa) (>exp(xb)) on log scale
  double logdiff(double xa, double xb);
  //log defective cdf diffusion first-passage time at lower boundary
  double logFjoint_lower(double q, double a, double v, double w);
  //derivative of log-diffusion density by threshold
  double dadwiener_d(double q, double a, double vn, double wn, double d);
  //derivative of log-diffusion density by start point
  double dwdwiener_d(double q, double a, double vn, double wn, double d);
  //set up model_design from file
  void model_design(int kerntree, int *ar, int *branch, int *nodes_per_tree, int *tree_and_node2par);
  //store results of last hamiltonian-gibbs cycle between batches
  void push(int ithread, int n_value_store, int n_all_parameters, gsl_vector* hampar, double* tavw, double* tlams, double* ai, double* loglambdas, double* bi, double* alltaus, double* rest,
            int trialno, int* paths, int* nips, double liknorm[6],double activeeps, double epsm, double Hobjective, double* valuestore, double* parmon, double* parmonstore);
  //retrieve results of last hamiltonian-gibbs cycle between batches
  void pop(int ithread, int n_value_store, int n_all_parameters, gsl_vector* hampar, double* tavw, double* tlams, double* ai, double* loglambdas, double* bi, double* alltaus, double* rest,
           int trialno, int* paths, int* nips, double liknorm[6], double& activeeps, double& epsm, double& Hobjective, double* valuestore, double* parmon, double* parmonstore);
  //sample path (augmented data)
  void make_path(trial one,  int* nips, int itrial, int& path, gsl_vector* hampar,  double *tavw, double* tlams, double* explambda, double *alltaus, double *rest, ars_archiv& ars_store, gsl_rng *rst);
  //show interim results
  void on_screen3(int n_all_parameters, double *xwbr, double *parmon, double *consts, double rmax, int imax, int irun);
  //initialize: flag = 0 random; flag = 1 using max. lik. personwise parameter estimates
  void initialize(int flag, const std::vector<trial> & daten, double xeps, double* parmonstore, int n_value_store, double* valuestore, gsl_rng* rst1, gsl_rng* rst2, gsl_rng* rst3, gsl_rng* rst4);
  //compute r-statistcs
  void r_statistic(int ido, int n_all_parameters, int istream, int iter, double *parmon, double *xwbr, double &rmax, int &imax);
  //the sampler
  void gibbs_times_new(const std::vector<trial> & daten, gsl_rng *rst1, gsl_rng *rst2, gsl_rng *rst3, gsl_rng *rst4);
  //expected mean of diffusion process completion time at threshold pm
  double exp_mean(int pm, double a, double v, double w);
  //log probability of diffusion process ending at threshold pm
  double logprob_upperbound(int pm, double a, double v, double w);
  //helper function for derivatives by threshold and drift rate of logprob_upperbound
  double davlogprob_upperbound(int pm, double a, double v, double w);
  //derivatives by threshold of logprob_upperbound
  double dalogprob_upperbound(int pm, double a, double v, double w, double dav);
  //derivatives by drift rate of logprob_upperbound
  double dvlogprob_upperbound(int pm,  double a, double v, double w, double dav);
  //derivatives by start point of logprob_upperbound
  double dwlogprob_upperbound(int pm,  double a, double v, double w);
  //variance of diffusion-model completion times at lower threshold
  double lower_bound_var(double a, double vn, double wn);
  //type (threshold, drift, start point) of parameter ip on ifree scale
  int is_type(int ip);
  //index of parameter ip on ifree[type] scale
  int ind_type(int type, int ip);
  //Gibbs-sampling of variance-covariance matrix of diffusion-model parameters (phase <=2)
  void sample_sig(gsl_vector* hampar, double* sig, double* sigi, gsl_matrix* cx, double* ai, gsl_rng* rst);
  //Gibbs-sampling of variance-covariance matrix of motor time means (phase <=2)
  void make_rgam(gsl_vector* hampar, double* gam, double* gami, gsl_matrix* cr,  double* bi, gsl_rng* rst);
  //Gibbs-sampling of residual variability Omega^2 (phase <=2)
  void make_romega(gsl_vector* hampar, double* explambdas, double& omega, gsl_rng* rst);
  //NUTs sampler for phase <=2
  bool hnuts(double* scale, int* nips, gsl_vector* hampar, double* tavw, double* tlams, double* sig, double* sigi, gsl_matrix* Ltminusx,
             const std::vector<trial>& daten, double* rscale, double* sl, double* rest, double* loglambdas, double* gam, double* gami, gsl_matrix* Ltminusr,
             double omega,
             double* alltaus, double& liknorm1, double& liknorm2, double& activeeps, double& epsm, double& Hobjective, int m, gsl_rng* rst);
  //compute nips
  void make_nips(const std::vector<trial> & daten, int* paths, int* nips);
  //logistic transformation
  double logit(transform e, double  x);
  //derivative thereof
  double dlogit(transform e, double  x);
  //inverse-logistic transformation
  double invlogit(transform par, double u);
  //log t density (proportional to; without Gamma factors)
  double log_tdist_pdf(double mu, double sig, int ddf, double x);
  //set parameters for logistic transformation
  transform prep_transform(double a, double b, double loc, double scale);
  //derivative of log diffusion density by t
  double dtdwiener_d(double q, double a, double vn, double wn, double& d);
  //adaptive rejection sampler for sampling diffusion-model completion times
  double make_rwiener(int t, int m, int pm, ars_archiv& ars_store, double bound, double a, double v, double w,  gsl_rng* rst);
  //set up adaptive rejection sampler
  void initialize_ars(int t, double* tavw, ars_archiv& ars_store);
  //compute summaries of results, dic, posterior model checks, posterior contrasts, etc.
  void diagnosis(const std::vector<trial> & daten, int* idaten, int kerntree, gsl_rng* rst);
  //store results for continuation in a separate run
  void push_continue(int n_value_store, int irun, double* valuestore, double* parmonstore, gsl_rng* rst1, gsl_rng* rst2, gsl_rng* rst3, gsl_rng* rst4);
  //retrieve results for continuation in a separate run
  void pop_continue(int n_value_store, int& irun, double* valuestore, double* parmonstore, gsl_rng* rst1, gsl_rng* rst2, gsl_rng* rst3, gsl_rng* rst4);
  //compute convolution of densities
  void convolution2(const std::vector<double> &rts, int pfadlength, int* low_or_up, double* a, double* v, double* w, double mu, double sig, std::vector<double>& pbranch);
  //maximum likelihood estimation of personwise diffusion-model parameters
  void tby_individuals(const std::vector<trial> &daten, double* avw, double* lambdas, gsl_rng* rst);
  //sample diffusion-model completion times by inverse cum. distr. function method
  double rwiener_diag(int pm, double bound, double a, double v, double w, gsl_rng* rst);
  //compute derivatives of likelihood by diffusion-model parameters
  void dhudwien(int* nips, gsl_vector* hampar, double* tavw, double* sigi, double* alltaus, double* dstore, gsl_vector* dhampar);
  //compute path and category probabilities for one individual
  void make_pij_for_individual(double* x, double* pij, double* pj);
  //create new Theta
  struct Theta* newTheta();
  //copy theta
  void thetacopy(Theta*& thetadest, Theta* thetasrc);
  //delete and free memory Theta
  void remove_Theta(Theta*& theta);
  //copy vecto
  void pcopy(gsl_vector* pdest, gsl_vector* psrc);
  //make tree structures for NUTs algorithm
  struct Node* make_tree2(int j);
  //provide storage for NUTs algorith
  struct store newstore();
  //delete and clear memory for Node
  void removenode(struct Node*& node);
  //update estimate of mass matrix
  void make_supersigs(int anz, double* parmonstore, gsl_matrix* supsig, gsl_matrix* sigisqrt);
  //NUTs sampler for phase >=3
  bool hnuts2(int* nips, gsl_vector* hampar, double* tavw, double* tlams,
              const std::vector<trial> & daten, double* rest, double* loglambdas,
              double* alltaus, double liknorm[6], double& activeeps, double& epsm, double& Hobjective, int m, bool save, gsl_rng* rst);
  //computes log(1-exp(z))
  double log1pem1(double z);
  //make variance-covariance matrix of diffusion-model parameters from model parameters in hampar
  void make_sigs(int flag, gsl_vector* hampar, double* sig);
  //the reverse direction of make_sigs
  void make_hampar_from_sig(int flag, double* sig, gsl_vector* hampar);
  //derivatives of likelihood for variance-covariance parameters
  void dhudext(gsl_vector* hampar, double* explambda, const std::vector<double>& zt, const std::vector<double>& zr, gsl_matrix* wt, gsl_matrix* wr,
               double etat, double etar, gsl_vector* dhampar);
  //compute contribution to likelihood from LKJ prior
  double joint_likeli4(int flag, gsl_vector* hampar, const std::vector<double> &z, gsl_matrix* w, double eta,
                       double taut, double liknorm4);
  //compute contribution to likelihood from Omega^2 prior
  double joint_likeli5(gsl_vector* hampar, double* explambda, double liknorm6);
  //LKJ transform helper function
  void from_y_to_z(int flag, gsl_vector* hampar, std::vector<double>& z);
  //LKJ transform helper function
  void from_z_to_w(int flag, const std::vector<double> &z, gsl_matrix* w);
  //sets up tree structure for NUTs sampler
  void preptrees(Node* trees[13]);
  //deletes and clears memory of tree structures for NUTs sampler
  void removetrees(Node* trees[13]);
  //make diffusion-model parameters from transformed parameters
  void make_tavwtlams(int flag, gsl_vector* hampar, std::vector<double>& z, gsl_matrix* w, double* tpars);
  //helper function to compute variance-covariance matrix from LKJ parameters
  void from_w_to_sig_sigi(int flag, gsl_vector* hampar, gsl_matrix* w, double* sig);

}

//read data from file
void lies(std::vector<trial> &daten);
//add exp(xa) and exp(xb) on log scale
double logsum(double xa, double xb);
//sets constants
void set_ns(const std::vector<trial> & daten, int& indi, int& kerntree, int& kerncat, int& igroup);
//log cdf standard normal
double lnnorm(double x);
//sample one uniform random variable
double oneuni(gsl_rng *rst);
//sample one standard normal random variable
double onenorm(gsl_rng *rst);
//truncnorm computes normal variate with sd = 1, mean = b, and value >= 0
double truncnorm(double b, gsl_rng *rst);


#endif
