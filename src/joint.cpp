#include "rts.h"

namespace drtmpt {

  //contributions to overall likelihood by diffusion-model parts; phase <=2
  double joint_likelihood(double* scale, gsl_vector* p, int* nips, gsl_vector* hampar, double* tavw, double* sig, double* sigi,  double* alltaus, double* dstore, double liknorm1) {

  	double h = 0;

  	gsl_vector_view t1 = gsl_vector_subvector(hampar, 0, igroup * icompg);
  	gsl_blas_ddot(&t1.vector, &t1.vector, &h);


  	h *= -PRIOR * 0.5;

  	double temp;
  	gsl_matrix_view SINV = gsl_matrix_view_array(sigi, icompg, icompg);
  	gsl_vector* zwischenv = gsl_vector_alloc(indi * icompg);
  	gsl_matrix_view zwischenm = gsl_matrix_view_vector(zwischenv, indi, icompg);
  	gsl_vector_view stackedv = gsl_vector_subvector(hampar, iavwoff, indi * icompg);
  	gsl_matrix_view stackedm = gsl_matrix_view_vector(&stackedv.vector, indi, icompg);
  	gsl_blas_dsymm(CblasRight, CblasLower, -0.5, &SINV.matrix, &stackedm.matrix, 0.0, &zwischenm.matrix);
  	gsl_blas_ddot(&stackedv.vector, zwischenv, &temp);
  	h += temp;

  /*	for (int t = 0; t != indi; t++) {
  		int jj = 0;
  		int itoff = t * 3 * ifreemax;
  		gsl_vector_view stacked = gsl_vector_subvector(hampar, iavwoff + t * icompg, icompg);
  		double temp;
  		gsl_vector_set_zero(zwischen);
  		gsl_blas_dsymv(CblasLower, -0.5, &SINV.matrix, &stacked.vector, 0.0, zwischen);
  		gsl_blas_ddot(&stacked.vector, zwischen, &temp);
  		h += temp;
  	}
  */

  	double xemp = liknorm1 / ntau;

  	int jj = 0;
  	for (int im = 0; im != no_patterns; im++) {
  		int iaa = dCOMB(im, 0);
  		int ivv = dCOMB(im, 1);
  		int iww = dCOMB(im, 2);
  		if ((dCOMP(0, iaa)) || (dCOMP(1, ivv)) || (dCOMP(2, iww)))
  			for (int t = 0; t != indi; t++) {
  				int itoff = t * 3 * ifreemax;
  				int ia = iaa + itoff;
  				int iv = ivv + itoff + ifreemax;
  				int iw = iww + itoff + 2 * ifreemax;

  				double aa = tavw[ia];
  				double va = tavw[iv];
  				double wa = tavw[iw];

  				int nntim = dNNODES(t, im);
  				for (int i = 0; i != nntim; i++) for (int pm = 0; pm != 2; pm++) {
  					h += dstore[jj] - xemp;
  					jj++;
  				}
  				for (int pm = 0; pm != 2; pm++) {
  					int n = dNIPS(t, pm, im);
  					if (n != 0) {
  						h -= n * logprob_upperbound(pm, aa, va, wa);
  					}
  				}
  			}
  		else for (int t = 0; t != indi; t++) {
  			int n = dNNODES(t, im);
  			jj += 2 * n; h -= 2 * n * xemp;
  		}
  	}

  	double current_proposed_k = 0.0; // current_k = 0.0;

  	gsl_vector_view p1 = gsl_vector_subvector(p, 0, iavwoff);
  	gsl_vector_view p2 = gsl_vector_view_array(scale, iavwoff);
  	gsl_vector* ptemp = gsl_vector_alloc(iavwoff);
  	gsl_vector_memcpy(ptemp, &p1.vector);
  	gsl_vector_mul(ptemp, &p2.vector);
  	gsl_blas_ddot(ptemp, ptemp, &temp);
  	current_proposed_k -= 0.5 * temp;
  	gsl_vector_free(ptemp);



  	gsl_matrix_view S = gsl_matrix_view_array(sig, icompg, icompg);
  	gsl_vector_view stackedv2 = gsl_vector_subvector(p, iavwoff, indi * icompg);
  	gsl_matrix_view stackedm2 = gsl_matrix_view_vector(&stackedv2.vector, indi, icompg);
  	gsl_blas_dsymm(CblasRight, CblasLower, -0.5, &S.matrix, &stackedm2.matrix, 0.0, &zwischenm.matrix);
  	gsl_blas_ddot(&stackedv2.vector, zwischenv, &temp);
  	current_proposed_k += temp;
  	gsl_vector_free(zwischenv);
  	return h + current_proposed_k;
  }

  //contribution to overall likelihood by motor-time parts and momentum variables
  double rjoint_likelihood(const std::vector<trial> & daten, double* rscale, double* sl, gsl_vector* p, double* rest, gsl_vector* hampar, double* loglambdas, double* gam, double* gami, double omega, double liknorm2) {


  	double temp = 0;

  	for (int x = 0; x != datenzahl; x++) {
  		int t = daten[x].person;
  		int r = cat2resp[daten[x].category], igr = t2group[t] * respno + r, itr = t * respno + r;
  		temp -= gsl_log1p(gsl_pow_2((rest[x] - gsl_vector_get(hampar, irmuoff + igr) - gsl_vector_get(hampar, ilamoff + itr)) / gsl_vector_get(hampar, isigoff + t)) / degf);
  	}

  	temp *= (degf + 1.0) / 2.0;

  	double xemp = liknorm2 / datenzahl;



  	for (int t = 0; t != indi; t++) {
  		int sdt = indi * respno + t;
  		double dsdt = gsl_vector_get(hampar, isigoff + t);
  		temp -= priordf * omega / (2 * gsl_pow_2(dsdt));

  		double lalogs = -loglambdas[t];
  		temp += (priordf + 1) * lalogs;
  		for (int r = 0; r != respno; r++) {
  			int igr = t2group[t] * respno + r, itr = t * respno + r;
  			temp -= nppr[itr] * log(gsl_cdf_tdist_P((gsl_vector_get(hampar, irmuoff + igr) + gsl_vector_get(hampar, ilamoff + itr)) / dsdt, degf));
  			temp += nppr[itr] * (lalogs - xemp);
  		}
  	}

  	gsl_matrix_view GINV = gsl_matrix_view_array(gami, respno, respno);
  	gsl_vector_view stacked = gsl_vector_subvector(hampar, ilamoff,indi* respno);
  	gsl_matrix_view stackedm = gsl_matrix_view_vector(&stacked.vector, indi, respno);
  	gsl_vector* zwischen = gsl_vector_alloc(indi*respno);
  	gsl_matrix_view zwischenm = gsl_matrix_view_vector(zwischen, indi, respno);

  	double help;
  //	gsl_vector_set_zero(zwischen);
  	gsl_blas_dsymm(CblasRight, CblasLower, -0.5, &GINV.matrix, &stackedm.matrix, 0.0, &zwischenm.matrix);
  	gsl_blas_ddot(&stacked.vector, zwischen, &help);
  	temp += help;



  	for (int ig = 0; ig != igroup; ig++) for (int r = 0; r != respno; r++) {
  		int igr = ig * respno + r;
  		temp -= 0.5 * gsl_pow_2((gsl_vector_get(hampar, irmuoff + igr) - mu_prior) / rsd);
  	}

  	double current_proposed_k = 0.0; // current_k = 0.0;

  	gsl_vector_view p1 = gsl_vector_subvector(p, irmuoff, igroup*respno);
  	gsl_vector_view p2 = gsl_vector_view_array(rscale, igroup*respno);
  	gsl_vector* ptemp = gsl_vector_alloc(igroup*respno);
  	gsl_vector_memcpy(ptemp, &p1.vector);
  	gsl_vector_mul(ptemp, &p2.vector);
  	double xtemp;
  	gsl_blas_ddot(ptemp, ptemp, &xtemp);
  	gsl_vector_free(ptemp);
  	current_proposed_k -= 0.5 * xtemp;

  	gsl_matrix_view G = gsl_matrix_view_array(gam, respno, respno);
  	gsl_vector_view stacked2 = gsl_vector_subvector(p, ilamoff,indi * respno);
  	gsl_matrix_view stackedm2 = gsl_matrix_view_vector(&stacked2.vector, indi, respno);
  	gsl_blas_dsymm(CblasRight, CblasLower, -0.5, &G.matrix, &stackedm2.matrix, 0.0, &zwischenm.matrix);
  	gsl_blas_ddot(&stacked2.vector, zwischen, &xtemp);
  	current_proposed_k += xtemp;

  	gsl_vector_free(zwischen);

  	gsl_vector_view p3 = gsl_vector_subvector(p, isigoff, indi);
  	gsl_vector_view p4 = gsl_vector_view_array(sl, indi);
  	gsl_vector* ptemp2 = gsl_vector_alloc(indi);
  	gsl_vector_memcpy(ptemp2, &p3.vector);
  	gsl_vector_mul(ptemp2, &p4.vector);
  	gsl_blas_ddot(ptemp2, ptemp2, &xtemp);
  	gsl_vector_free(ptemp2);
  	current_proposed_k -= 0.5 * xtemp;

  	return(current_proposed_k + temp);
  }

  //derivatives log-likelihood by motor-time related parameters
  void   dhudlam(const std::vector<trial> & daten, double* rest, gsl_vector* hampar, double* gami, double omega, gsl_vector* dhampar) {

  	double* dlam = (double*)calloc(indi * (respno + 1), sizeof(double));
  	double* dhrmu = (double*)calloc(igroup * respno, sizeof(double));



  	for (int x = 0; x != datenzahl; x++) {
  		int t = daten[x].person;
  		int r = cat2resp[daten[x].category];
  		int itr = t * respno + r, igr = t2group[t] * respno + r;
  		double xa = gsl_vector_get(hampar, irmuoff + igr) + gsl_vector_get(hampar, ilamoff + itr) - rest[x];
  		double xx = xa / (1.0 + gsl_pow_2(xa / gsl_vector_get(hampar, isigoff + t)) / degf);
  		dlam[itr] += xx;
  		dlam[indi * respno + t] -= xa * xx;
  	}

  	double dings = (degf + 1.0) / (1.0 * degf);
  	for (int t = 0; t != indi; t++) {
  		int it = indi * respno + t;
  		double sig = gsl_vector_get(hampar, isigoff + t);
  		dlam[it] *= dings;

  		dlam[it] -= priordf * omega;
  		dlam[it] /= gsl_pow_3(sig);
  		dlam[it] += (priordf + 1) / sig;
  		for (int r = 0; r != respno; r++) {
  			int itr = t * respno + r, igr = t2group[t] * respno + r;
  			dlam[itr] *= dings / gsl_pow_2(sig);
  			double mu = gsl_vector_get(hampar, irmuoff + igr) + gsl_vector_get(hampar, ilamoff + itr);

  			double factor = gsl_ran_tdist_pdf(mu / sig, degf) / gsl_cdf_tdist_P(mu / sig, degf);
  			dlam[itr] -= -nppr[itr] / sig * factor;


  			dlam[it] -= -nppr[itr] * (-mu / gsl_pow_2(sig)) * factor;
  			dlam[it] += nppr[itr] / sig;

  		}
  	}
  	for (int r = 0; r != respno; r++) {
  		for (int t = 0; t != indi; t++) dhrmu[t2group[t] * respno + r] += dlam[t * respno + r];
  		for (int ig = 0; ig != igroup; ig++) {
  			int igr = ig * respno + r;
  			dhrmu[igr] += (gsl_vector_get(hampar, irmuoff + igr) - mu_prior) / rsd / rsd;
  			//	dhrmu[igr] = -dhrmu[igr];
  		}
  	}



  	gsl_matrix_view gg = gsl_matrix_view_array(gami, respno, respno);
  	gsl_vector_view ll = gsl_vector_subvector(hampar, ilamoff, indi * respno);
  	gsl_matrix_view lll = gsl_matrix_view_vector(&ll.vector, indi, respno);

  	gsl_matrix_view dlll = gsl_matrix_view_array(dlam, indi, respno);

  	gsl_blas_dsymm(CblasRight, CblasLower, 1.0, &gg.matrix, &lll.matrix, 1.0, &dlll.matrix);

  	gsl_vector_view t1 = gsl_vector_subvector(dhampar, ilamoff, indi * (respno + 1));
  	gsl_vector_view t2 = gsl_vector_view_array(dlam, indi * (respno + 1));
  	gsl_vector_memcpy(&t1.vector, &t2.vector);
  	gsl_vector_view t3 = gsl_vector_subvector(dhampar, irmuoff, igroup * respno);
  	gsl_vector_view t4 = gsl_vector_view_array(dhrmu, igroup * respno);
  	gsl_vector_memcpy(&t3.vector, &t4.vector);

  	free(dlam); free(dhrmu);
  }


  //one leapfrog cycle
  void Leapfrog(int* nips, double* scale, gsl_vector* hampar, double* tavw, double* tlams, gsl_vector* dhampar, double* sig, double* sigi, const std::vector<trial> & daten, double* rscale, double* sl, double* loglambdas, double* gam, double* gami,
  	double* alltaus, double* dstore, double* rest, double omega, double eps, gsl_vector* p) {

  	gsl_blas_daxpy(-eps / 2, dhampar, p);


  	//Parameters

  	gsl_vector* ptemp = gsl_vector_alloc(nhamil);
  	gsl_vector_memcpy(ptemp, p);
  	gsl_vector_view p1 = gsl_vector_subvector(ptemp, 0, iavwoff);
  	gsl_vector_view p2 = gsl_vector_view_array(scale, iavwoff);

  	gsl_vector_mul(&p1.vector, &p2.vector);
  	gsl_vector_mul(&p1.vector, &p2.vector);

  	gsl_matrix_view S = gsl_matrix_view_array(sig, icompg, icompg);

  	gsl_vector_view pt = gsl_vector_subvector(p, iavwoff, indi * icompg);
  	gsl_matrix_view ptm = gsl_matrix_view_vector(&pt.vector, indi, icompg);
  	gsl_vector_view p3a = gsl_vector_subvector(ptemp, iavwoff, indi * icompg);
  	gsl_matrix_view p3b = gsl_matrix_view_vector(&p3a.vector, indi, icompg);
  	gsl_blas_dsymm(CblasRight, CblasLower, 1.0, &S.matrix, &ptm.matrix, 0.0, &p3b.matrix);


  	gsl_vector_view p4 = gsl_vector_subvector(ptemp, irmuoff, igroup * respno);
  	gsl_vector_view p5 = gsl_vector_view_array(rscale, igroup * respno);

  	gsl_vector_mul(&p4.vector, &p5.vector);
  	gsl_vector_mul(&p4.vector, &p5.vector);

  	gsl_vector_view p7 = gsl_vector_subvector(ptemp, isigoff, indi);
  	gsl_vector_view p8 = gsl_vector_view_array(sl, indi);

  	gsl_vector_mul(&p7.vector, &p8.vector);
  	gsl_vector_mul(&p7.vector, &p8.vector);


  	gsl_matrix_view G = gsl_matrix_view_array(gam, respno, respno);
  	gsl_vector_view ppt = gsl_vector_subvector(p, ilamoff, indi * respno);
  	gsl_matrix_view pptm = gsl_matrix_view_vector(&ppt.vector, indi, respno);
  	gsl_vector_view p6a = gsl_vector_subvector(ptemp, ilamoff, indi * respno);
  	gsl_matrix_view p6b = gsl_matrix_view_vector(&p6a.vector, indi, respno);
  	gsl_blas_dsymm(CblasRight, CblasLower, 1.0, &G.matrix, &pptm.matrix, 0.0, &p6b.matrix);

  	gsl_blas_daxpy(eps, ptemp, hampar);

  	for (int t = 0; t != indi; t++) {
  		int iofft = iavwoff + t * icompg, itfr = 3 * t * ifreemax;
  		for (int type = 0; type != 3; type++) {
  			int ift = ifree[type];
  			for (int ip = 0; ip != ift; ip++) if (dCOMP(type, ip)) {
  				int iavw = itfr + type * ifreemax + ip;
  				tavw[iavw] = logit(avwtrans[type], gsl_vector_get(hampar, dmapMAVW(t2group[t], type, ip)) + gsl_vector_get(hampar, dmapAVW(t, type, ip)));
  			}
  		}
  	}

  	for (int t = 0; t != indi; t++) for (int r = 0; r != respno; r++)
  		tlams[t * respno + r] = gsl_vector_get(hampar, irmuoff + t2group[t] * respno + r) + gsl_vector_get(hampar, ilamoff + t * respno + r);

  	for (int t = 0; t != indi; t++) {
  		int takt = isigoff + t;
  		if (gsl_vector_get(hampar, takt) < 0) {
  			gsl_vector_set(hampar, takt, -gsl_vector_get(hampar, takt));
  			gsl_vector_set(p, takt, -gsl_vector_get(p, takt));
  //			std::cout << "bounce back ";
  		}
  		loglambdas[t] = log(gsl_vector_get(hampar, takt));
  	}

  	gsl_vector_free(ptemp);

  	dhudwien(nips, hampar, tavw, sigi, alltaus, dstore, dhampar);
  	dhudlam(daten, rest, hampar, gami, omega, dhampar);

  	gsl_blas_daxpy(-eps / 2, dhampar, p);

  }

  //one Hamiltonian cycle
  double step0(int* nips, double* scale, double* sig, double* sigi, gsl_vector* dhampar, const std::vector<trial> & daten, double* rscale,  double* rest, double* gam, double* gami, double omega, double* sl, double* alltaus, struct Theta* theta,  gsl_vector* p, double u, int v, int j, double eps, int &n, int &s, double liknorm, double liknorm2) {
  	double* loglambdas = theta->loglambda;
  	double* tavw = theta->tavw;
  	double* tlams = theta->tlams;
  	gsl_vector* hampar = theta->hampar;
  	double* dstore = 0; if (!(dstore = (double*)malloc(ntau * sizeof(double)))) { Rprintf("Allocation failure2\n"); }
  	Leapfrog(nips, scale, hampar, tavw, tlams, dhampar, sig, sigi, daten, rscale, sl, loglambdas,  gam, gami, alltaus, dstore, rest, omega, v * eps, p);
  	double temp = joint_likelihood(scale, p, nips, hampar, tavw, sig, sigi, alltaus, dstore, liknorm) + rjoint_likelihood(daten, rscale, sl, p, rest, hampar, loglambdas, gam, gami, omega, liknorm2);
  	n = (u <= temp) ? 1 : 0;
  	s = (u - 1000.0 < temp) ? 1 : 0;
  	free(dstore);
  	return temp;
  }

  //inner product for U-turn criterion
  double inner_product(gsl_vector* p, gsl_vector* hamparp, gsl_vector* hamparm) {
  	gsl_vector* stacked = gsl_vector_alloc(nhamil);
  	gsl_vector_memcpy(stacked, hamparp);

  	gsl_vector_sub(stacked, hamparm);
  	double inner;

  	gsl_blas_ddot(stacked, p, &inner);
  	gsl_vector_free(stacked);
  	return inner;
  }

  //Build tree in NUTs algorithm
  void buildtree(int* nips, double* scale, double* sig, double* sigi, const std::vector<trial> & daten, double* rscale, double* rest, double* gam, double* gami,
  	double omega, double* sl, 	 double* alltaus, struct Theta* theta, struct Theta* thetadash, gsl_vector* dhampar, gsl_vector* p,
  	double u, int v, int j, double eps, gsl_rng* rst, double liknorm, double liknorm2, int& ndash, int& sdash, int& nadash, double& adash, bool adapt) {

  	std::stack<struct Node*> stack;
  	int n, s = 1;
  	Node* root = (j<13)? trees[j]: make_tree2(j);
  	std::vector<struct store> speicher;
  	int jp1 = j + 1;
  	for (int i = 0; i != jp1; i++) speicher.push_back(newstore());

  	do {
  		while (root != NULL) {
  			if (root->right != NULL) stack.push(root->right);
  			stack.push(root);
  			root = root->left;
  		}
  		root = stack.top(); stack.pop();
  		if ((root->right != NULL) && (!(stack.empty())) && (root->right == stack.top())) {
  			stack.pop();
  			stack.push(root); root = root->right;
  		}
  		else {
  			if ((root->status < 2) && (s == 1)) {
  				double temp = step0(nips, scale, sig, sigi, dhampar, daten, rscale, rest, gam, gami, omega, sl, alltaus, theta, p, u, v, j, eps, n, s, liknorm, liknorm2);
  				if (adapt) {
  					if (temp > 0) temp = 0.0;
  					speicher[root->index].a = temp;
  					speicher[root->index].na = 1;
  				}
  				speicher[root->index].n = n; speicher[root->index].s = s;
  				if (root->status == 0) {
  					thetacopy(speicher[root->index].theta, theta);
  					pcopy(speicher[root->index].p, p);
  				}
  				thetacopy(speicher[root->index].thetadash, theta);
  			};
  			if (root->status == 2) {
  				if (speicher[root->index].s == 1) {
  					int i = root->index, ipl = i + root->level;
  					int ndd = speicher[i + 1].n, nd = speicher[i].n;
  					//					if ((ndd + nd) == 0) std::cout << "und nu?";
  					if ((ndd > 0) && ((ndd + nd) * oneuni(rst) <= ndd)) {
  						thetacopy(speicher[i].thetadash, speicher[i + 1].thetadash);
  					}
  					speicher[i].s = speicher[i + 1].s;
  					if (speicher[i].s == 1)
  						if ((v * inner_product(speicher[i].p, theta->hampar, speicher[i].theta->hampar) < 0.0) ||
  							(v * inner_product(p, theta->hampar, speicher[i].theta->hampar) < 0.0))
  						{
  							speicher[i].s = 0; s = 0;
  						}
  					speicher[i].n += ndd;
  					if (adapt) {
  						speicher[i].a = logsum(speicher[i].a, speicher[i + 1].a);
  						speicher[i].na += speicher[i + 1].na;
  					}
  				}
  			}
  			root = NULL;
  		}
  	} while (!stack.empty());

  	thetacopy(thetadash, speicher[0].thetadash);
  	ndash = speicher[0].n; sdash = speicher[0].s;
  	if (adapt) { adash = speicher[0].a; nadash = speicher[0].na; }
  	for (int i = 0; i != jp1; i++) {
  		remove_Theta(speicher[i].theta);
  		remove_Theta(speicher[i].thetadash);
  		gsl_vector_free(speicher[i].p);
  	}
  }

  //the sampler; phase <= 2
  bool hnuts(double* scale, int* nips, gsl_vector* hampar, double* tavw, double* tlams, double* sig, double* sigi, gsl_matrix* Ltminusx,
  	const std::vector<trial> & daten, double* rscale, double* sl, double* rest, double * loglambdas,  double* gam, double* gami, gsl_matrix* Ltminusr,
  	double omega, double* alltaus, double &liknorm1, double &liknorm2, double& activeeps, double& epsm, double& Hobjective, int m, gsl_rng* rst) {


  	double lamtest = gsl_vector_get(hampar, 0);

  	bool adapt = (phase == 1) && (m <= PHASE1);
  	int xn = isigoff + indi;

  	double* dstore = 0; if (!(dstore = (double*)malloc(ntau * sizeof(double)))) { Rprintf("Allocation failure2\n"); }
  	gsl_vector* p = gsl_vector_alloc(nhamil);
  	gsl_vector* pp = gsl_vector_alloc(nhamil);
  	gsl_vector* pm = gsl_vector_alloc(nhamil);



  	gsl_vector* dhamparp = gsl_vector_alloc(nhamil);
  	gsl_vector* dhamparm = gsl_vector_alloc(nhamil);

  	dhudlam(daten, rest, hampar, gami, omega, dhamparp);
  	dhudwien(nips, hampar, tavw, sigi, alltaus, dstore, dhamparp);

  	gsl_vector_memcpy(dhamparm, dhamparp);


  	int icig = icompg * igroup;
  	for (int i = 0; i != icig; i++) {
  		gsl_vector_set(p,i, onenorm(rst) / scale[i]);
  	}

  	int iaicin = iavwoff + icompg * indi;
  	for (int i = iavwoff; i != iaicin; i++) gsl_vector_set(p, i,onenorm(rst));


  	gsl_vector_view ptemp = gsl_vector_subvector(p, iavwoff, indi * icompg);
  	gsl_matrix_view pptemp = gsl_matrix_view_vector(&ptemp.vector, indi, icompg);

  	gsl_blas_dtrmm(CblasRight, CblasLower, CblasTrans, CblasNonUnit, 1.0, Ltminusx, &pptemp.matrix);


  	liknorm1 += joint_likelihood(scale, p, nips,hampar, tavw, sig, sigi, alltaus, dstore, liknorm1);
  	free(dstore);

  	int ioff = icompg * (indi + igroup);

  	int irig = respno * igroup;
  	for (int i = 0; i != irig; i++) {
  		gsl_vector_set(p,ioff + i, onenorm(rst) / rscale[i]);
  	}

  	ioff += respno * igroup;
  	int iorin = ioff + respno * indi;
  	for (int i = ioff; i != iorin; i++) gsl_vector_set(p , i, onenorm(rst));

  	gsl_vector_view ptempn = gsl_vector_subvector(p, ioff, indi * respno);
  	gsl_matrix_view pptempn = gsl_matrix_view_vector(&ptempn.vector, indi, respno);

  	gsl_blas_dtrmm(CblasRight, CblasLower, CblasTrans, CblasNonUnit, 1.0, Ltminusr, &pptempn.matrix);

  	ioff += respno * indi;
  	for (int t = 0; t != indi; t++) {
  		int takt = ioff + t;
  		gsl_vector_set(p,takt, onenorm(rst) / sl[t]);
  	}

  	liknorm2 += rjoint_likelihood(daten, rscale, sl, p, rest, hampar, loglambdas, gam, gami, omega, liknorm2);

  	Theta* thetap = newTheta();
  	Theta* thetam = newTheta();
  	Theta* thetadash = newTheta();
  	Theta* thetac = (struct Theta*)malloc(sizeof(struct Theta));
  	thetac->hampar = hampar;  thetac->tavw = tavw; thetac->loglambda = loglambdas;
  	thetac->tlams = tlams;
  	thetacopy(thetap, thetac); thetacopy(thetam, thetac);
  	int n, s;
  	pcopy(pp, p); pcopy(pm, p);
  	n = 1;
  	s = 1;

  	int na; double a;

  	double u = log(oneuni(rst));
  	int j = 0;
  	while (s == 1) {
  		int ndash, sdash;
  		int v = (oneuni(rst) <= 0.5) ? -1 : 1;
  		if (v == -1) {
  			buildtree(nips, scale, sig, sigi, daten, rscale, rest, gam, gami, omega, sl, alltaus,
  				thetam, thetadash, dhamparm, pm, u, v, j, activeeps, rst, liknorm1, liknorm2, ndash, sdash, na, a, adapt);
  		}
  		else {
  			buildtree(nips, scale, sig, sigi, daten, rscale, rest, gam, gami, omega, sl, alltaus,
  				thetap, thetadash, dhamparp, pp, u, v, j, activeeps, rst, liknorm1, liknorm2, ndash, sdash, na, a, adapt);
  		}
  		if ((sdash == 1) && (n * oneuni(rst) <= ndash))
  		{
  			thetacopy(thetac, thetadash);
  		}
  		n += ndash;
  		j = j + 1;
  		if (j == maxtreedepth1_3) sdash = 0.0;
  		s = sdash;
  		if (s == 1)
  			if ((inner_product(pm, thetap->hampar, thetam->hampar) < 0.0) ||
  				(inner_product(pp, thetap->hampar, thetam->hampar) < 0.0))
  				s = 0;
  	}
  	if (adapt) {
  		double xtemp = 1.0 / (m + t0plus);
  		Hobjective = (1.0 - xtemp) * Hobjective + xtemp * (0.60 - exp(a) / na);
  		xtemp = muplus - sqrt(1.0 * m) / gamplus * Hobjective;
  		double xtemp2 = pow(m * 1.0, -kapplus);
  		epsm = xtemp2 * xtemp + (1 - xtemp2) * epsm;
  		activeeps = exp(xtemp);
  	}
  //	std::cout << setw(5) <<j << setw(5) << n;
  	free(thetac);
  	remove_Theta(thetam);
  	remove_Theta(thetap);
  	remove_Theta(thetadash);
  	gsl_vector_free(p); gsl_vector_free(pp); gsl_vector_free(pm);
  	gsl_vector_free(dhamparp); gsl_vector_free(dhamparm);

  	return(gsl_vector_get(hampar,0) != lamtest);
  }

}
