// authors: Christoph Klauer and Raphael Hartmann

#include "rts.h"
#include <ctime>
#include <chrono>

int kerncat;
int kernpar;
int zweig;
int nodemax;
int *ng;
int indi;
int *t2group=0;
int igroup;
int respno;
int *cat2resp = 0;



//set some constants
void set_ns(const std::vector<trial> & daten, int &indi, int &kerntree, int &kerncat, int &igroup) {
  indi = 0; kerntree = 0; kerncat = 0; igroup = 0;
  for (int i = 0; i != static_cast<int>(daten.size()); i++) {
    trial one = daten[i];
    indi = std::max(indi, one.person);
    igroup = std::max(igroup, one.group);
  }
  indi++; igroup++;
  //ntot = static_cast<int>(daten.size());
  std::ifstream info(MODEL);
  info >> zweig;
  info >> kernpar;
  info >> nodemax;
  info >> kerntree;
  info >> kerncat;
  info.close();
}

//from categories to trees
void set_cat2tree(std::vector<trial> & daten, int *cat2tree)
{
  std::ifstream info(MODEL); int schrott;
  for (int j = 0; j != 5; j++) info >> schrott;
  for (int j = 0; j != kerncat; j++) { info >> cat2tree[j]; cat2tree[j]--; }

  for (int i = 0; i != datenzahl; i++) {
    daten[i].tree = cat2tree[daten[i].category];
  }
  info.close();
}

//assign persons to groups t2group; number of persons per group ng
void set_t2group(const std::vector<trial> & daten, int* t2group, int* ng) {

  for (int i = 0; i != datenzahl; i++) {
    trial one = daten[i];
    t2group[one.person] = one.group;
  }

  for (int t = 0; t != indi; t++) ng[t2group[t]]++;
}


namespace ertmpt {

  // Globale Variablen
  int *cat2tree=0;
  int *ar=0;
  int *branch=0;
  int *nodes_per_tree=0;
  int *tree_and_node2par=0;
  bool *comp=0;
  int ifree, ilamfree;
  int ipred;
  int *ndrin=0, *drin=0;
  // int n_all_parameters;
  int *nppr=0;
  // int n_bridge_parameters;

  int RMAX_reached ;
  bool BURNIN_flag;

  int ireps;
  int alphaoff;
  int sigalphaoff;
  int restparsno;
  int *free2kern=0;
  int *kern2free=0;
  double *consts=0;

  int *pfad_index = 0;
  std::vector<pfadinfo> path_info;


  void make_nodes_by_ind(int *idaten, int kerntree, int *nodes_per_par, int &nz, int *nnodes, int &ntau) {
  #define NNODES(I,J) nnodes[I*kernpar+J]
  #define NKS(T,I) nks[T*kerntree+I]
  #define IDATEN(T,I) idaten[T*kerncat + I]
  #define NODES_PER_PAR(I,J) nodes_per_par[I*kernpar + J]
  #define TREE_AND_NODE2PAR(I,J) tree_and_node2par[I*nodemax+J]

  	int *nks = 0; nks = (int *)malloc(indi*kerntree * sizeof(int));
  	for (int i = 0; i != kerntree * indi; i++) nks[i] = 0;
  	for (int i = 0; i != kerncat; ++i) for (int t = 0; t != indi; t++) NKS(t, cat2tree[i]) += IDATEN(t, i);
  	for (int t = 0; t != indi; t++) for (int ip = 0; ip != kernpar; ip++) {
  		NNODES(t, ip) = 0;
  		for (int k = 0; k != kerntree; k++) NNODES(t, ip) += NODES_PER_PAR(k, ip)*NKS(t, k);
  	}
  	nz = ntau = 0;
  	for (int ip = 0; ip != kernpar; ip++) {
  		if (comp[ip]) for (int t = 0; t != indi; t++) nz += NNODES(t, ip);
  		if (comp[ip + kernpar]) for (int t = 0; t != indi; t++) ntau += NNODES(t, ip);
  		if (comp[ip + 2 * kernpar]) for (int t = 0; t != indi; t++) ntau += NNODES(t, ip);
  	}
  	if (nks) free(nks);
  }

  void make_positions(std::vector<trial> daten, int *nnodes, int *nz_position, int *ntau_position) {
  #define NZ_POSITION(X,J) nz_position[X*nodemax+J]
  #define NTAU_POSITION(X,J,PM) ntau_position[2*X*nodemax+2*J+PM]
  	int *boffset = 0; boffset = (int *)malloc(indi*kernpar * sizeof(int));
  	int *loffset = 0; loffset = (int *)malloc(indi*kernpar * sizeof(int));
  	int *btemp = 0; btemp = (int *)malloc(indi*kernpar * sizeof(int));
  	int *ltemp = 0; ltemp = (int *)malloc(indi*kernpar * sizeof(int));
  #define BOFFSET(T,I) boffset[T*kernpar+I]
  #define LOFFSET(T,I) loffset[T*kernpar+I]
  #define BTEMP(T,I) btemp[T*kernpar+I]
  #define LTEMP(T,I) ltemp[T*kernpar+I]

  	int trialno = static_cast<int>(daten.size());
  	for (int i = 0; i != indi * kernpar; i++) boffset[i] = loffset[i] = btemp[i] = ltemp[i] = 0;
  	int jj = 0;
  	for (int ip = 0; ip != kernpar; ip++) if (comp[ip])
  		for (int t = 0; t != indi; t++) {
  			BOFFSET(t, ip) = jj;
  			jj += NNODES(t, ip);
  		}
  	jj = 0;
  	for (int ip = 0; ip != kernpar; ip++)
  		for (int t = 0; t != indi; t++) if ((comp[ip + kernpar]) || (comp[ip + 2 * kernpar])) {
  			LOFFSET(t, ip) = jj;
  			if (comp[ip + kernpar]) jj += NNODES(t, ip);
  			if (comp[2 * kernpar + ip]) jj += NNODES(t, ip);
  		}

  	for (int i = 0; i != nodemax * trialno; i++) nz_position[i] = -1;
  	for (int i = 0; i != 2 * nodemax * trialno; i++) ntau_position[i] = -1;

  	for (int i = 0; i != trialno; i++) {
  		int itree = daten[i].tree, t = daten[i].person;
  		for (int r = 0; r != nodes_per_tree[itree]; r++) {
  			int ip = TREE_AND_NODE2PAR(itree, r);
  			if (comp[ip]) { NZ_POSITION(i, r) = BTEMP(t, ip) + BOFFSET(t, ip); BTEMP(t, ip)++; }
  			if (comp[ip + kernpar]) { NTAU_POSITION(i, r, 0) = LTEMP(t, ip) + LOFFSET(t, ip); LTEMP(t, ip)++; }
  			if (comp[ip + 2 * kernpar]) { NTAU_POSITION(i, r, 1) = LTEMP(t, ip) + LOFFSET(t, ip); LTEMP(t, ip)++; }
  		}
  	}
  	if (DEBUG) for (int t = 0; t != indi; t++) for (int i = 0; i != kernpar; i++) if ((comp[i]) && (BTEMP(t, i) != NNODES(t, i))) Rprintf("PROBLEM%12d%12d\n", t, i);
  	for (int t = 0; t != indi; t++) for (int i = 0; i != kernpar; i++) {
  		if (((comp[i + kernpar]) && (comp[i + 2 * kernpar])) && (LTEMP(t, i) != 2 * NNODES(t, i))) Rprintf("L_PROBLEM%12d%12d\n", t, i);
  		if (((comp[i + kernpar]) && !(comp[i + 2 * kernpar])) && (LTEMP(t, i) != NNODES(t, i))) Rprintf("L_PROBLEM%12d%12d\n", t, i);
  		if ((!(comp[i + kernpar]) && (comp[i + 2 * kernpar])) && (LTEMP(t, i) != NNODES(t, i))) Rprintf("L_PROBLEM%12d%12d\n", t, i);
  		if ((!(comp[i + kernpar]) && !(comp[i + 2 * kernpar])) && (LTEMP(t, i) != 0)) Rprintf("L_PROBLEM%12d%12d\n", t, i);
  	}
  	if (boffset) free(boffset);
  	if (loffset) free(loffset);
  	if (btemp) free(btemp);
  	if (ltemp) free(ltemp);
  }


  void r_statistic(int ido, int n_all_parameters, int istream, int iter, double *parmon, double *xwbr, double &rmax) {

  	if (ido == 1) for (int i = 0; i != 3 * n_all_parameters; i++) xwbr[i] = 0;
  	double r = 1.0 / (istream + 1);
  #define XWBR(T,I) xwbr[(T-1)*n_all_parameters+I]
  #define PARMON(I,J) parmon[(I-1)*n_all_parameters+J]

  	for (int i = 0; i != n_all_parameters; i++) {
  		XWBR(2, i) += gsl_pow_2(PARMON(1, i) - XWBR(3, i))*(1 - r);
  		XWBR(3, i) += (PARMON(1, i) - XWBR(3, i))*r;
  		XWBR(1, i) += (PARMON(2, i) - XWBR(1, i))*r;
  	}
  	if (ido == 3) {
  		rmax = 0.0;
  		for (int i = 0; i != n_all_parameters; i++) {
  			XWBR(3, i) = sqrt((XWBR(1, i) / iter + XWBR(2, i) / istream) / XWBR(1, i)*(iter - 1));
  			if (XWBR(3, i) > rmax) rmax = XWBR(3, i);
  			if (XWBR(1, i) == 0.0) Rprintf("XWBR(1,i)=0, i= %d\n", i);
  		}
  	}
  }

  // MAIN
  int mainx(int *k2f, int *f2k) {
  	//std::cout << "L_Test" << std::setw(12) << 19.2 << std::setw(12) << 3 << std::endl;
  	//printf ("L_Test %12F%12d \n", 15, 3);
  	// Rprintf ("L_TEST%12.2f%12d", 14.2, 3); Rprintf ("L_TEST%12.2f\n", 14.2);

    // setting some variables
    ipred=0;
    ireps=IREP;
    RMAX_reached = 0;
    BURNIN_flag = true;

  	std::vector<trial> daten;
  	int exit_status = 0;

  	//	NagError fail;	INIT_FAIL(fail);
  		/* Choose the base generator */
  		/* Random Seed */

  	gsl_rng *rst1;
  	long int seed = std::time(0); seed = abs(seed * seed);
  	if (DEBUG) Rprintf("%d\n", (seed <= 0));
  	rst1 = gsl_rng_alloc(T_rng);
  	gsl_rng_set(rst1, seed);
  	long int n = gsl_rng_max(rst1) / 2;


  	gsl_rng *rst2;
  	seed = gsl_rng_uniform_int(rst1, n) + 1;
  	if (DEBUG) Rprintf("%d\n", (seed <= 0));
  	rst2 = gsl_rng_alloc(T_rng);
  	gsl_rng_set(rst2, seed);


  	gsl_rng *rst3; rst3 = gsl_rng_alloc(T_rng);
  	if (NOTHREADS >= 3) {
  		seed = gsl_rng_uniform_int(rst1, n) + 1;
  		if (DEBUG) Rprintf("%d\n", (seed <= 0));
  		gsl_rng_set(rst3, seed);
  	}


  	gsl_rng *rst4; rst4 = gsl_rng_alloc(T_rng);
  	if (NOTHREADS >= 4) {
  		seed = gsl_rng_uniform_int(rst1, n) + 1;
  		if (DEBUG) Rprintf("%d\n", (seed <= 0));
  		gsl_rng_set(rst4, seed);
  	}


  	gsl_rng *rst5; rst5 = gsl_rng_alloc(T_rng);
  	if (NOTHREADS >= 5) {
  		seed = gsl_rng_uniform_int(rst1, n) + 1;
  		if (DEBUG) Rprintf("%d\n", (seed <= 0));
  		gsl_rng_set(rst5, seed);
  	}


  	gsl_rng *rst6; rst6 = gsl_rng_alloc(T_rng);
  	if (NOTHREADS >= 6) {
  		seed = gsl_rng_uniform_int(rst1, n) + 1;
  		if (DEBUG) Rprintf("%d\n", (seed <= 0));
  		gsl_rng_set(rst6, seed);
  	}


  	gsl_rng *rst7; rst7 = gsl_rng_alloc(T_rng);
  	if (NOTHREADS >= 7) {
  		seed = gsl_rng_uniform_int(rst1, n) + 1;
  		if (DEBUG) Rprintf("%d\n", (seed <= 0));
  		gsl_rng_set(rst7, seed);
  	}


  	gsl_rng *rst8; rst8 = gsl_rng_alloc(T_rng);
  	if (NOTHREADS >= 8) {
  		seed = gsl_rng_uniform_int(rst1, n) + 1;
  		if (DEBUG) Rprintf("%d\n", (seed <= 0));
  		gsl_rng_set(rst8, seed);
  	}


  	gsl_rng *rst9; rst9 = gsl_rng_alloc(T_rng);
  	if (NOTHREADS >= 9) {
  		seed = gsl_rng_uniform_int(rst1, n) + 1;
  		if (DEBUG) Rprintf("%d\n", (seed <= 0));
  		gsl_rng_set(rst9, seed);
  	}


  	gsl_rng *rst10; rst10 = gsl_rng_alloc(T_rng);
  	if (NOTHREADS >= 10) {
  		seed = gsl_rng_uniform_int(rst1, n) + 1;
  		if (DEBUG) Rprintf("%d\n", (seed <= 0));
  		gsl_rng_set(rst10, seed);
  	}


  	gsl_rng *rst11; rst11 = gsl_rng_alloc(T_rng);
  	if (NOTHREADS >= 11) {
  		seed = gsl_rng_uniform_int(rst1, n) + 1;
  		if (DEBUG) Rprintf("%d\n", (seed <= 0));
  		gsl_rng_set(rst11, seed);
  	}


  	gsl_rng *rst12; rst12 = gsl_rng_alloc(T_rng);
  	if (NOTHREADS >= 12) {
  		seed = gsl_rng_uniform_int(rst1, n) + 1;
  		if (DEBUG) Rprintf("%d\n", (seed <= 0));
  		gsl_rng_set(rst12, seed);
  	}


  	gsl_rng *rst13; rst13 = gsl_rng_alloc(T_rng);
  	if (NOTHREADS >= 13) {
  		seed = gsl_rng_uniform_int(rst1, n) + 1;
  		if (DEBUG) Rprintf("%d\n", (seed <= 0));
  		gsl_rng_set(rst13, seed);
  	}


  	gsl_rng *rst14; rst14 = gsl_rng_alloc(T_rng);
  	if (NOTHREADS >= 14) {
  		seed = gsl_rng_uniform_int(rst1, n) + 1;
  		if (DEBUG) Rprintf("%d\n", (seed <= 0));
  		gsl_rng_set(rst14, seed);
  	}


  	gsl_rng *rst15; rst15 = gsl_rng_alloc(T_rng);
  	if (NOTHREADS >= 15) {
  		seed = gsl_rng_uniform_int(rst1, n) + 1;
  		if (DEBUG) Rprintf("%d\n", (seed <= 0));
  		gsl_rng_set(rst15, seed);
  	}


  	gsl_rng *rst16; rst16 = gsl_rng_alloc(T_rng);
  	if (NOTHREADS >= 16) {
  		seed = gsl_rng_uniform_int(rst1, n) + 1;
  		if (DEBUG) Rprintf("%d\n", (seed <= 0));
  		gsl_rng_set(rst16, seed);
  	}


  	gsl_rng *rst;   rst = gsl_rng_alloc(T_rng);
  	gsl_rng_memcpy(rst, rst1);
  	// std::cout << gsl_rng_size(rst) << std::endl;
  	if (DEBUG) Rprintf("%d\n", static_cast<int>(gsl_rng_size(rst)));
  	// std::cout << gsl_rng_name(rst) << std::endl;
  	if (DEBUG) Rprintf("%s\n", gsl_rng_name(rst));
  	// std::cout << gsl_rng_uniform(rst) << std::endl;
  	if (DEBUG) Rprintf("%f\n", gsl_rng_uniform(rst));
  	// char xx; std::cin >> xx;



  	// Data
  	int kerntree;
  	lies(daten);
  	datenzahl = static_cast<int>(daten.size());
  	
  	
  	//REPEAT:
  	set_ns(daten, indi, kerntree, kerncat, igroup);
  	cat2tree = (int *)malloc(kerncat * sizeof(int));
  	set_cat2tree(daten, cat2tree);
  	t2group = (int *)malloc(indi * sizeof(int));
  	ng = (int*)calloc(igroup, sizeof(int));
  	set_t2group(daten, t2group, ng);
  	int *idaten = 0; idaten = (int *)malloc(indi*kerncat * sizeof(int));
  	make_idaten(daten, idaten);
  	//	std::cout<< std::endl;
  	//	for (int t=0;t!=indi;t++) {for (int j=0;j!=kerncat;j++) std::cout << std::setw(4) << IDATEN(t,j); std::cout<< std::endl;}

		// Model Design
  	// in lies gesetzt	zweig=2; kernpar = 2; nodemax=2;
  	//a = (int *)malloc(kerncat*zweig*kernpar * sizeof(int));
  	ar = (int *)malloc(kerncat*zweig*nodemax * sizeof(int));
  	//    #define A(I,J,K) a[I*zweig*kernpar + J*kernpar + K]
  	//if (!(b = (int *)malloc(kerncat*zweig*kernpar * sizeof(int)))) { printf("Allocation failure\n");	exit_status = -1; }
  	branch = (int *)malloc(kerncat * sizeof(int));
  	//    #define B(I,J,K) b[I*zweig*kernpar + J*kernpar + K]
  	int *nodes_per_par = 0;  nodes_per_par = (int *)malloc(kerntree*kernpar * sizeof(int));
  	//    #define NODES(I,J) nodes_per_par[I*kernpar + J]
  	nodes_per_tree = (int *)malloc(kerntree * sizeof(int));
  	tree_and_node2par = (int *)malloc(kerntree*nodemax * sizeof(int));
  	drin = (int *)malloc(kerncat*zweig*nodemax * sizeof(int));
  	ndrin = (int *)malloc(kerncat*zweig * sizeof(int));
  	pfad_index = (int *)malloc(kerncat*zweig * sizeof(int));

  	// Model specifications
  		// Parameter: beta_comp yes/no
  	comp = (bool *)malloc(3 * kernpar * sizeof(bool));
  	consts = (double *)malloc(kernpar * sizeof(double));

  	model_design(kerntree, ar, branch, nodes_per_par, nodes_per_tree, tree_and_node2par);

  	ifree = 0;
    ilamfree = 0;
    int ilamfree1 = 0;
    int ilamfree2 = 0;
    for (int ip = 0; ip < 3*kernpar; ip++) {
      if (ip < kernpar) {
        if (k2f[ip] == ifree) {
          ifree++;
        }
      } else if (ip < 2*kernpar) {
        if (k2f[ip] == ifree + ilamfree1) {
          ilamfree1++;
          ilamfree++;
        }
      } else {
        if (k2f[ip] == ifree + ilamfree1 + ilamfree2) {
          ilamfree2++;
          ilamfree++;
        }
      }
    }


  	free2kern = (int *)malloc((ifree + ilamfree) * sizeof(int));
  	kern2free = (int *)malloc(3 * kernpar * sizeof(int));

  	for (int ip = 0; ip != 3 * kernpar; ip++) {
  		kern2free[ip] = k2f[ip];
  	}
  	for (int iz = 0; iz != ifree + ilamfree; iz++) {
  		free2kern[iz] = f2k[iz];
  	}

  	// int iz = 0; for (int ip = 0; ip != 3 * kernpar; ip++) if (comp[ip]) { kern2free[ip] = iz; free2kern[iz] = ip; iz++; }
  	// else kern2free[ip] = -1;


  	extract_pfadinfo(pfad_index, path_info);




  		// Analysis by individuals

  	//double *g2 = 0; if (!(g2 = (double *)malloc(indi * sizeof(double)))) { printf("Allocation failure\n");	exit_status = -1; }
  	//double *likeli = 0; if (!(likeli = (double *)malloc(indi * sizeof(double)))) { printf("Allocation failure\n");	exit_status = -1; }


  	//by_individuals(daten, kerntree, beta, g2, likeli, rst);

  	//nnodes berechnen
  	int *nnodes = 0; nnodes = (int *)malloc(indi*kernpar * sizeof(int));
  	int nz, ntau;
  	make_nodes_by_ind(idaten, kerntree, nodes_per_par, nz, nnodes, ntau);

  	//NZ und NTAU Positions berechnen
  	int trialno = static_cast<int>(daten.size());
  	int *nz_position = 0; nz_position = (int *)malloc(trialno*nodemax * sizeof(int));
  	int *ntau_position = 0; ntau_position = (int *)malloc(2 * trialno*nodemax * sizeof(int));
  	make_positions(daten, nnodes, nz_position, ntau_position);

  	//nppr berechnen, factor definieren
  	nppr = (int *)malloc(indi*respno * sizeof(int));

  	for (int t = 0; t != indi * respno; t++) nppr[t] = 0;
  	for (int x = 0; x != static_cast<int>(daten.size()); x++)
  	{
  		int t = daten[x].person; int r = cat2resp[daten[x].category];
  		nppr[t*respno + r]++;
  	}





  	alphaoff = igroup * respno + 1 + (respno*(respno + 1)) / 2;
  	sigalphaoff = alphaoff + indi * respno;
  	restparsno = sigalphaoff + indi;

  	double *beta = 0;	beta = (double *)malloc(indi*ifree * sizeof(double));;
  	double *lambdas = 0; lambdas = (double *)malloc(indi*(ilamfree) * sizeof(double));
  	double *restpars = 0; restpars = (double *)malloc(restparsno * sizeof(double));

  	//compute individual fits for starting points

  	for (int ip = 0; ip != kernpar; ip++) if (!comp[ip]) {
  		consts[ip] = gsl_cdf_ugaussian_Pinv(consts[ip]);
  	}

  	if (generate_or_diagnose) tby_individuals(daten, kerntree, beta, lambdas, restpars, rst);

  	Rprintf("\nStart sampling from the posterior distribution:\n\n");

  	n_all_parameters = ifree*igroup     + ilamfree*igroup       + ((ifree+ilamfree)*(ifree+ilamfree+1))/2       + indi*ifree         + indi*ilamfree           + restparsno;
  	n_bridge_parameters = n_all_parameters + ifree + ilamfree + respno;

  	gsl_rng_memcpy(rst1, rst);
  	if (generate_or_diagnose) gibbs_times_new(daten, nnodes, nz, nz_position, beta, ntau, ntau_position,
  																						rst1, rst2, rst3, rst4, rst5, rst6, rst7, rst8, rst9, rst10, rst11, rst12, rst13, rst14, rst15, rst16,
  																						lambdas, restpars);


  	if (lambdas) free(lambdas);
  	if (restpars) free(restpars);
  	if (beta) free(beta);

  	Rprintf("\nCalculating some diagnostics. This might take some time.\n\n");

  	diagnosis(daten, idaten, kerntree, rst);
  	// char x; std::cin >> x;


  	if (cat2tree) free(cat2tree);
  	free(t2group);
  	//if (a) free(a);
  	if (ar) free(ar);
  	//if (b) free(b);
  	if (branch) free(branch);
  	if (nodes_per_par) free(nodes_per_par);
  	if (nodes_per_tree) free(nodes_per_tree);
  	if (tree_and_node2par) free(tree_and_node2par);
  	//if (g2) free(g2);
  	//if (likeli) free(likeli);
  	if (nnodes) free(nnodes);
  	if (idaten) free(idaten);

  	if (comp) free(comp);
  	if (nz_position) free(nz_position);
  	if (ntau_position) free(ntau_position);
  	if (drin) free(drin);
  	if (ndrin) free(ndrin);
  	if (nppr) free(nppr);

  	if (free2kern) free(free2kern);
  	if (kern2free) free(kern2free);
  	if (consts) free(consts);
  	if (pfad_index) free(pfad_index);
  	gsl_rng_free(rst);
  	gsl_rng_free(rst1);
  	gsl_rng_free(rst2);
  	gsl_rng_free(rst3);
  	gsl_rng_free(rst4);
  	gsl_rng_free(rst5);
  	gsl_rng_free(rst6);
  	gsl_rng_free(rst7);
  	gsl_rng_free(rst8);
  	gsl_rng_free(rst9);
  	gsl_rng_free(rst10);
  	gsl_rng_free(rst11);
  	gsl_rng_free(rst12);
  	gsl_rng_free(rst13);
  	gsl_rng_free(rst14);
  	gsl_rng_free(rst15);
  	gsl_rng_free(rst16);
  	return exit_status;
  }

}




namespace drtmpt {

  bool PROG_BAR_FLAG;


  //CPUs to be used in computing DIC
  int DIC_CPUs;// = NOTHREADS;

  //CPUs to be used if INITIALIZE = 1
  int INIT_CPUs;// = NOTHREADS;


  // Other global variables


  Node* trees[13] = {0,0,0,0,0,0,0,0,0,0,0,0,0};

  double* supersig = 0;
  int nhamil;
  int phase;

  int *cat2tree = 0;
  int *ar = 0;
  int *branch = 0;
  int *nodes_per_tree = 0;
  int *tree_and_node2par = 0;
  int* tree_and_node2map = 0;
  //bool *comp = 0;
  //int ifree[3], ifreeg, ifreemax;
  int ifreeg, ifreemax;
  int icomp[3], icompg;
  int* nnodes = 0;
  int* ndrin = 0, * drin = 0, * cdrin = 0, * ncdrin = 0, * pfadmax = 0;
  int *nppr = 0;
  int *tau_by_node=0;
  int* n_per_subj = 0;

  int iavwoff, irmuoff, ilamoff, isigoff;

  int ireps;// = IREP;

  int RMAX_reached;

  //int* kern2free = 0;
  int *free2comp = 0;
  // double *consts = 0;
  double *monitor = 0;

  int* map = 0;
  int* comb = 0;
  int no_patterns;
  std::vector<double> rtmins;

  int* mapmavw;
  int* mapavw;

  transform avwtrans[3];



  int ntau;

  int sample_size;

  double muplus;

  gsl_matrix* supsig;
  gsl_matrix* sigisqrt;


  //number of nodes per combination m and person nnodes(t,m);
  void make_nodes_by_ind(const std::vector<trial> & daten, int kerntree, int* nodes_per_tree, int* nnodes, int* n_per_subj) {
    ntau = 0;
    for (int i = 0; i != indi * no_patterns; i++) nnodes[i] = 0;
    for (int t = 0; t != indi; t++) n_per_subj[t] = 0;

    for (int x = 0; x != datenzahl; x++) {
      trial one = daten[x]; int t = one.person, itree = one.tree;
      n_per_subj[t]++;
      for (int n = 0; n != nodes_per_tree[itree]; n++) {
        int ia = dTREE_AND_NODE2PAR(itree, n, 0);
        int iv = dTREE_AND_NODE2PAR(itree, n, 1);
        int iw = dTREE_AND_NODE2PAR(itree, n, 2);
        int m = dMAP(ia, iv, iw);
        dNNODES(t, m)++;
        ntau += 2;
      }
    }
  }

  //map parameter combinations on parameters (comb) and reverse (map); assign tree, node and parameter type to index on ifree[type] scale (tree_and_node2par) and on combination no (tree_and_node2map)
  void make_map(int kerntree, int& no_patterns, int* tree_and_node2map) {
    if (!(map = (int*)malloc(ifree[0] * ifree[1] * ifree[2] * sizeof(int)))) { Rprintf("Allocation failure\n"); }
    if (!(comb = (int*)malloc(3 * kernpar * sizeof(int)))) { Rprintf("Allocation failure\n"); }

    for (int i = 0; i != ifree[0] * ifree[1] * ifree[2]; i++) map[i] = -1;
    for (int i = 0; i != 3 * kernpar; i++) comb[i] = -1;

    no_patterns = 0;
    for (int itree = 0; itree != kerntree; itree++)
      for (int n = 0; n != nodes_per_tree[itree]; n++)
      {
        bool neu = true;
        int ia = dTREE_AND_NODE2PAR(itree, n, 0);
        int iv = dTREE_AND_NODE2PAR(itree, n, 1);
        int iw = dTREE_AND_NODE2PAR(itree, n, 2);
        for (int m = 0; (m != no_patterns) && (neu); m++) {
          if ((ia == dCOMB(m, 0)) && (iv == dCOMB(m, 1)) && (iw == dCOMB(m, 2))) neu = false;
        }
        if (neu) {
          dMAP(ia, iv, iw) = no_patterns;
          dCOMB(no_patterns, 0) = ia; dCOMB(no_patterns, 1) = iv; dCOMB(no_patterns, 2) = iw;
          no_patterns++;
        }
      }
      for (int itree = 0; itree != kerntree; itree++)
        for (int n = 0; n != nodes_per_tree[itree]; n++)
        {
          int ia = dTREE_AND_NODE2PAR(itree, n, 0);
          int iv = dTREE_AND_NODE2PAR(itree, n, 1);
          int iw = dTREE_AND_NODE2PAR(itree, n, 2);

          dTREE_AND_NODE2MAP(itree, n) = dMAP(ia, iv, iw);
        }
  }

  //compute positions of tau in double* alltaus by tree, node, and parameter type (threshold, drift, start point)
  void make_positions(const std::vector<trial> & daten, int* tau_by_node) {
  #define dLOFFSET(T,I) loffset[T*no_patterns+I]
    int* loffset = 0;
    if (!(loffset = (int*)malloc(indi * no_patterns * sizeof(int)))) { Rprintf("Allocation failure\n");  }
    int* ltemp = 0; if (!(ltemp = (int*)malloc(indi * no_patterns * sizeof(int)))) { Rprintf("Allocation failure\n");  }


  #define dLTEMP(T,I) ltemp[T*no_patterns+I]


    for (int i = 0; i != indi * no_patterns; i++) loffset[i] = ltemp[i] = 0;
    int jj = 0;
    for (int im = 0; im != no_patterns; im++)
      for (int t = 0; t != indi; t++) {
        dLOFFSET(t, im) = jj;
        jj += 2* dNNODES(t, im);
      }

      for (int i = 0; i != 2 * nodemax * datenzahl; i++) tau_by_node[i] = -1;

    for (int i = 0; i != datenzahl; i++) {
      int itree = daten[i].tree, t = daten[i].person;
      for (int r = 0; r != nodes_per_tree[itree]; r++) {
        int ia = dTREE_AND_NODE2PAR(itree, r, 0);
        int iv = dTREE_AND_NODE2PAR(itree, r, 1);
        int iw = dTREE_AND_NODE2PAR(itree, r, 2);
        int im = dMAP(ia, iv, iw);
        dTAU_BY_NODE(i, r, 0) = dLTEMP(t, im) + dLOFFSET(t, im); dLTEMP(t, im)++;
        dTAU_BY_NODE(i, r, 1) = dLTEMP(t, im) + dLOFFSET(t, im); dLTEMP(t, im)++;
      }
    }
    for (int t = 0; t != indi; t++) for (int i = 0; i != no_patterns; i++) {
      //		if (dLTEMP(t, i) != 2 * NNODES(t, i)) std::cout << "L_PROBLEM" << setw(12) << t << setw(12) << i << setw(12) << dLTEMP(t,i) << setw(12) << NNODES(t,i) << std::endl;
    }
    if (ltemp) free(ltemp);
    if (loffset) free(loffset);
  }

  //frequency data by person and category
  void make_idaten(const std::vector<trial> & daten, int *idaten) {
    trial one;
    for (int i = 0; i != indi * kerncat; i++) idaten[i] = 0;
    for (int i = 0; i != datenzahl; i++) {
      one = daten[i];
      dIDATEN(one.person, one.category)++;
    }
  }

  //number of responses by person and response
  void compute_nppr(const std::vector<trial> & daten) {
    //nppr berechnen,
    if (!(nppr = (int*)malloc(indi * respno * sizeof(int)))) { Rprintf("Allocation failure\n"); }
    for (int t = 0; t != indi * respno; t++) nppr[t] = 0;
    for (int x = 0; x != datenzahl; x++)
    {
      int t = daten[x].person; 
      int r = cat2resp[daten[x].category];
      nppr[t * respno + r]++;
    }
  }

  //number of interior edges (n,o) on paths to category C (ncdrin) and list of them (cdrin); number of nodes in path no. k ending in c (ndrin) and list of them (drin)
  void make_drin_cdrin() {

    for (int i = 0; i != kerncat * 2 * nodemax * 2; i++) cdrin[i] = -1;
    for (int i = 0; i != kerncat; i++) ncdrin[i] = 0;

    for (int j = 0; j != kerncat; j++)
      for (int r = 0; r != nodes_per_tree[cat2tree[j]]; r++) {
        bool flag[2] = { false,false };
        for (int k = 0; k != branch[j]; k++) if (dAR(j, k, r) != 0) {
          int pm = (1 + dAR(j, k, r)) / 2;
          if (!(flag[pm])) {
            dCDRIN(j, dNCDRIN(j), 0) = r;
            dCDRIN(j, dNCDRIN(j), 1) = pm;
            dNCDRIN(j)++;
            flag[pm] = true; if (flag[1 - pm]) break;
          }
        }
      }


      for (int i = 0; i != kerncat * zweig * nodemax; i++) drin[i] = -1;
    for (int i = 0; i != kerncat * zweig; i++) ndrin[i] = 0;

    for (int j = 0; j != kerncat; j++)
      for (int k = 0; k != branch[j]; k++)
        for (int r = 0; r != nodes_per_tree[cat2tree[j]]; r++) if (dAR(j, k, r) != 0) {
          dDRIN(j, k, dNDRIN(j, k)) = r;
          dNDRIN(j, k)++;
        }
        for (int j = 0; j!= kerncat; j++) {
          int k = branch[j];
          pfadmax[j] = 0;
          for (int p = 0; p != k; p++) pfadmax[j] = std::max(dNDRIN(j, p), pfadmax[j]);
        }
  }

  //map diffusion-model parameters by type and index ip on ifree[type]-scales on position in hampar
  void make_parameter_maps(int* mapmavw, int* mapavw) {
    int jj = 0;
    for (int ig = 0; ig != igroup; ig++)
      for (int type = 0; type != 3; type++)
        for (int ip = 0; ip != ifree[type]; ip++)
          if (dCOMP(type, ip)) {
            dmapMAVW(ig, type, ip) = jj; jj++;
          }

          for (int t = 0; t != indi; t++)
            for (int type = 0; type != 3; type++)
              for (int ip = 0; ip != ifree[type]; ip++)
                if (dCOMP(type, ip)) {
                  dmapAVW(t, type, ip) = jj; jj++;
                }
  }

  //minimum response time associated with parameter combinations, nodes, and outcome
  void make_rtmins(const std::vector<trial> & daten, std::vector<double>& rtmins) {
    std::vector<double> tmincat(indi * kerncat, GSL_POSINF);

    //	for (int t = 0; t != indi; t++) for (int j = 0; j != kerncat; j++) tmincat.push_back(GSL_POSINF);

    for (int i = 0; i != indi * no_patterns * 2; i++) rtmins.push_back(GSL_POSINF);

    for (int x = 0; x != datenzahl; x++) {
      trial one = daten[x];
      int t = one.person;
      int j = one.category;
      tmincat[t * kerncat + j] = fmin(tmincat[t * kerncat + j], one.rt / 1000.0);
    }

    for (int t = 0; t != indi; t++) for (int j = 0; j != kerncat; j++) {
      int itree = cat2tree[j];
      for (int k = 0; k != branch[j]; k++) {
        for (int in = 0; in != dNDRIN(j, k); in++) {
          int n = dDRIN(j, k, in);
          int ia, iv, iw;
          ia = dTREE_AND_NODE2PAR(itree, n, 0);
          iv = dTREE_AND_NODE2PAR(itree, n, 1);
          iw = dTREE_AND_NODE2PAR(itree, n, 2);
          int m = dMAP(ia, iv, iw);
          int pm = (dAR(j, k, n) == 1) ? 1 : 0;
          rtmins[t * no_patterns * 2 + m * 2 + pm] = fmin(rtmins[t * no_patterns * 2 + m * 2 + pm], tmincat[t * kerncat + j]/10.0);
        }
      }
    }
  }



  // MAIN
  int
  main_d() {
    //main(int argc, char* argv[]) {

    // set some variables
    ireps = IREP;
    DIC_CPUs = MAXTHREADS;
    INIT_CPUs = MAXTHREADS;
    PROG_BAR_FLAG = true;
    nhamil = 0;
    phase = 1;
    RMAX_reached = 0;
    
    
    std::vector<trial> daten;
    int exit_status = 0;
    long int seed = static_cast<long>(std::time(nullptr)); seed = abs(seed * seed);
    gsl_rng* rst1;   rst1 = gsl_rng_alloc(T_rng); gsl_rng_set(rst1, seed);
    long int n = gsl_rng_max(rst1) / 2;
    seed = gsl_rng_uniform_int(rst1, n) + 1;
    gsl_rng* rst2;   rst2 = gsl_rng_alloc(T_rng); gsl_rng_set(rst2, seed);
    seed = gsl_rng_uniform_int(rst1, n) + 1;
    gsl_rng* rst3;   rst3 = gsl_rng_alloc(T_rng); gsl_rng_set(rst3, seed);
    seed = gsl_rng_uniform_int(rst1, n) + 1;
    gsl_rng* rst4;   rst4 = gsl_rng_alloc(T_rng); gsl_rng_set(rst4, seed);

    gsl_rng* rst;   rst = gsl_rng_alloc(T_rng);
    gsl_rng_memcpy(rst, rst1);

    int kerntree;
    lies(daten);
    datenzahl = static_cast<int>(daten.size());
    set_ns(daten, indi, kerntree, kerncat, igroup);
    if (!(cat2tree = (int*)malloc(kerncat * sizeof(int)))) { Rprintf("Allocation failure\n"); }
    set_cat2tree(daten, cat2tree);
    // nur f�r generate
    int* idaten = 0; if (!(idaten = (int*)malloc(indi * kerncat * sizeof(int)))) { Rprintf("Allocation failure\n"); }
    make_idaten(daten, idaten);
    
    // Model Design
    if (!(ar = (int*)malloc(kerncat * zweig * nodemax * sizeof(int)))) { Rprintf("Allocation failure\n"); }
    if (!(branch = (int*)malloc(kerncat * sizeof(int)))) { Rprintf("Allocation failure\n"); }
    if (!(nodes_per_tree = (int*)malloc(kerntree * sizeof(int)))) { Rprintf("Allocation failure\n"); }
    if (!(tree_and_node2par = (int*)malloc(kerntree * nodemax * 3 * sizeof(int)))) { Rprintf("Allocation failure\n"); }
    if (!(tree_and_node2map = (int*)malloc(kerntree * nodemax * sizeof(int)))) { Rprintf("Allocation failure\n"); }
    if (!(drin = (int*)malloc(kerncat * zweig * nodemax * sizeof(int)))) { Rprintf("Allocation failure\n"); }
    if (!(ndrin = (int*)malloc(kerncat * zweig * sizeof(int)))) { Rprintf("Allocation failure\n"); }
    if (!(cdrin = (int*)malloc(kerncat * 2 * (nodemax * 2) * sizeof(int)))) { Rprintf("Allocation failure\n"); }
    if (!(ncdrin = (int*)malloc(kerncat * sizeof(int)))) { Rprintf("Allocation failure\n"); }
    if (!(pfadmax = (int*)malloc(kerncat * sizeof(int)))) { Rprintf("Allocation failure\n"); }
    // if (!(kern2free = (int*)malloc(kernpar * 3 * sizeof(int)))) { Rprintf("Allocation failure\n"); }
    // Model specifications
    // Parameter: beta_comp yes/no
    // if (!(comp = (bool*)malloc(3 * kernpar * sizeof(bool)))) { Rprintf("Allocation failure\n"); }
    // if (!(consts = (double*)malloc(3 * kernpar * sizeof(double)))) { Rprintf("Allocation failure\n"); };
    
    model_design(kerntree, ar, branch, nodes_per_tree, tree_and_node2par);
    make_drin_cdrin();
    avwtrans[0] = prep_transform(1.0e-2, 1.0e2, 0.8, 0.2);
    avwtrans[1] = prep_transform(-1.0e2, 1.0e2, 0.0, 1.0);
    avwtrans[2] = prep_transform(1.0e-3, 0.999, 0.5, 0.1);
    make_map(kerntree, no_patterns, tree_and_node2map);
    compute_nppr(daten);
    nnodes = (int*)malloc(indi * no_patterns * sizeof(int));
    n_per_subj = (int*)malloc(indi * sizeof(int));
    make_nodes_by_ind(daten, kerntree, nodes_per_tree, nnodes, n_per_subj);
    //NTAU Positions berechnen
    if (!(tau_by_node = (int*)malloc(2 * datenzahl * nodemax * sizeof(int)))) { Rprintf("Allocation failure\n"); }
    make_positions(daten, tau_by_node);
    if (!(t2group = (int*)malloc(indi * sizeof(int)))) { Rprintf("Allocation failure\n"); }
    if (!(ng = (int*)calloc(igroup, sizeof(int)))) { Rprintf("Allocation failure\n"); }
    set_t2group(daten, t2group, ng);
    make_rtmins(daten, rtmins);
    mapmavw = (int*)calloc(igroup * ifreemax * 3, sizeof(int));
    mapavw = (int*)calloc(indi * ifreemax * 3, sizeof(int));
    make_parameter_maps(mapmavw, mapavw);
    iavwoff = igroup * icompg;
    irmuoff = (indi + igroup) * icompg;
    ilamoff = irmuoff + igroup * respno;
    isigoff = ilamoff + indi * respno;
    nhamil = (igroup + indi) * (icompg + respno) + indi;
    n_all_parameters = icompg * igroup + indi * icompg + (icompg * (icompg + 1)) / 2 + respno * igroup + (respno + 1) * indi + (respno * (respno + 1)) / 2 + 1;
    //                    ma,mv,mw                 a,v,w            sig                           rmu     lambdas+sig_t           gam                     omega
    supsig = gsl_matrix_alloc(n_all_parameters, n_all_parameters);
    sigisqrt = gsl_matrix_alloc(n_all_parameters, n_all_parameters);
    if (generate_or_diagnose) gibbs_times_new(daten, rst1, rst2, rst3, rst4);
    diagnosis(daten, idaten, kerntree, rst);

    if (cat2tree) free(cat2tree);
    if (ar) free(ar);
    if (branch) free(branch);
    if (nodes_per_tree) free(nodes_per_tree);
    if (tree_and_node2par) free(tree_and_node2par);
    if (tree_and_node2map) free(tree_and_node2map);

    if (idaten) free(idaten);
    //if (comp) free(comp);
    if (tau_by_node) free(tau_by_node);
    if (drin) free(drin);
    if (ndrin) free(ndrin);
    if (cdrin) free(cdrin);
    if (ncdrin) free(ncdrin);
    if (pfadmax) free(pfadmax);
    if (nnodes) free(nnodes);
    if (n_per_subj) free(n_per_subj);
    if (nppr) free(nppr);
    if (map) free(map);
    if (comb) free(comb);
    free(t2group);
    free(ng);

    //if (kern2free) free(kern2free);
    if (free2comp) free(free2comp);
    //if (consts) free(consts);
    gsl_rng_free(rst);
    gsl_rng_free(rst1);
    gsl_rng_free(rst2);
    gsl_rng_free(rst3);
    gsl_rng_free(rst4);

    gsl_matrix_free(supsig);
    gsl_matrix_free(sigisqrt);
    free(mapavw);
    free(mapmavw);

    return exit_status;
  }

}
