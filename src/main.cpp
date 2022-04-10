// authors: Christoph Klauer and Raphael Hartmann

#include "rts.h"
#include <ctime>
#include "main.h"


// Globale Variablen
int *cat2tree=0;
int kernpar;
int kerncat;
int indi;
int zweig;
int *ar=0;
int *branch=0;
int nodemax;
int *nodes_per_tree=0;
int *tree_and_node2par=0;
bool *comp=0;
int ifree, ilamfree;
int ipred=0;
int *ndrin=0, *drin=0;
// int n_all_parameters;
int *nppr=0;
// int n_bridge_parameters;

int RMAX_reached = 0;
bool BURNIN_flag = true;

int igroup;
int *t2group=0;
int ireps=IREP;
int *cat2resp=0;
int respno;
int alphaoff;
int sigalphaoff;
int restparsno;
int *free2kern=0;
int *kern2free=0;
double *consts=0;

int *pfad_index = 0;
vector<pfadinfo> path_info;


void set_ns(vector<trial> daten, int &indi, int &kerntree, int &kerncat, int &igroup, int &ntot)
{
	indi = 0; kerntree = 0; kerncat = 0; ntot = 0; igroup = 0;
	for (int i = 0; i != static_cast<int>(daten.size()); i++) {
		trial one = daten[i];
		if (one.person > indi) indi = one.person;
		igroup = std::max(igroup, one.group);
	}
	indi++; kerntree++; kerncat++; igroup++;
	ntot = static_cast<int>(daten.size());
	std::ifstream info(MODEL);
	info >> zweig;
	info >> kernpar;
	info >> nodemax;
	info >> kerntree;
	info >> kerncat;
	info.close();
}

void set_cat2tree(vector<trial> &daten, int *cat2tree)
{
	std::ifstream info(MODEL); int schrott;
	for (int j = 0; j != 5; j++) info >> schrott;
	for (int j = 0; j != kerncat; j++) { info >> cat2tree[j]; cat2tree[j]--; }
	int ntot = static_cast<int>(daten.size());
	for (int i = 0; i != ntot; ++i) {
		daten[i].tree = cat2tree[daten[i].category];
	}
	info.close();
}

void set_t2group(vector<trial> daten, int *t2group) {
	int ntot = static_cast<int>(daten.size());
	for (int i = 0; i != ntot; ++i) {
		trial one = daten[i];
		t2group[one.person] = one.group;
	}
}

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

void make_positions(vector<trial> daten, int *nnodes, int *nz_position, int *ntau_position) {
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
	//std::cout << "L_Test" << setw(12) << 19.2 << setw(12) << 3 << std::endl;
	//printf ("L_Test %12F%12d \n", 15, 3);
	// Rprintf ("L_TEST%12.2f%12d", 14.2, 3); Rprintf ("L_TEST%12.2f\n", 14.2);

	vector<trial> daten;
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
	if (DEBUG) Rprintf("%d\n", gsl_rng_size(rst));
	// std::cout << gsl_rng_name(rst) << std::endl;
	if (DEBUG) Rprintf("%s\n", gsl_rng_name(rst));
	// std::cout << gsl_rng_uniform(rst) << std::endl;
	if (DEBUG) Rprintf("%f\n", gsl_rng_uniform(rst));
	// char xx; std::cin >> xx;




	// Data
	int kerntree, ntot;
	lies(daten);

	//REPEAT:
	set_ns(daten, indi, kerntree, kerncat, igroup, ntot);

	cat2tree = (int *)malloc(kerncat * sizeof(int));
	set_cat2tree(daten, cat2tree);
	t2group = (int *)malloc(indi * sizeof(int));
	set_t2group(daten, t2group);
	int *idaten = 0; idaten = (int *)malloc(indi*kerncat * sizeof(int));
	make_idaten(daten, idaten);
	//	std::cout<< std::endl;
	//	for (int t=0;t!=indi;t++) {for (int j=0;j!=kerncat;j++) std::cout << setw(4) << IDATEN(t,j); std::cout<< std::endl;}

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
