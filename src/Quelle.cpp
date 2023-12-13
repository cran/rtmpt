#include "rts.h"

namespace drtmpt {

  //diffusion-model contributions to likelihood (phase> 2)
  double joint_likelihood2(int* nips, gsl_vector* hampar, double* tavw,  double* alltaus, double* dstore, double liknorm1) {

  	double h = 0;

  	gsl_vector_view t1 = gsl_vector_subvector(hampar, 0, igroup * icompg);
  	gsl_blas_ddot(&t1.vector, &t1.vector, &h);
  	h *= -PRIOR * 0.5;

  	double temp;

  	gsl_vector_view stackedv = gsl_vector_subvector(hampar, iavwoff, indi * icompg);
  	gsl_blas_ddot(&stackedv.vector, &stackedv.vector, &temp);
  	h += -0.5 * temp;

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
  //Beschleunigung m�glich
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
  	return h ;
  }

  //motor-time densities contributions to model likelihood (phase > 2)
  double rjoint_likelihood2(const std::vector<trial> & daten, double* rest, gsl_vector* hampar, double* tlams, double* explambdas, double omega, double liknorm2) {

  	double temp = 0;

  	for (int x = 0; x != datenzahl; x++) {
  		int t = daten[x].person;
  		int r = cat2resp[daten[x].category], igr = t2group[t] * respno + r, itr = t * respno + r;
  		temp -= gsl_log1p(gsl_pow_2((rest[x] - tlams[itr]) / explambdas[t]) / degf);
  	}
  	temp *= (degf + 1.0) / 2.0;

  	double xemp = liknorm2 / datenzahl;



  	for (int t = 0; t != indi; t++) {
  		temp -= priordf * omega / (2 * gsl_pow_2(explambdas[t]));
  		// Jacoby-Faktor
  		temp += gsl_vector_get(hampar, isigoff + t);
  		// Jacoby end;
  		double lalogs = -gsl_vector_get(hampar, isigoff + t);
  		temp += (priordf + 1) * lalogs;
  		for (int r = 0; r != respno; r++) {
  			int igr = t2group[t] * respno + r, itr = t * respno + r;
  			temp -= nppr[itr] * log(gsl_cdf_tdist_P((tlams[itr]) / explambdas[t], degf));
  			temp += nppr[itr] * (lalogs - xemp);
  		}
  	}

  	double help;
  	gsl_vector_view stackedv = gsl_vector_subvector(hampar, ilamoff, indi * respno);
  	gsl_blas_ddot(&stackedv.vector, &stackedv.vector, &help);

  	temp += -0.5 * help;

  	for (int ig = 0; ig != igroup; ig++) for (int r = 0; r != respno; r++) {
  		int igr = ig * respno + r;
  		temp -= 0.5 * gsl_pow_2((gsl_vector_get(hampar, irmuoff + igr) - mu_prior) / rsd);
  	}
  	return(temp);
  }

  //contribution by momentum variables to likelihood
  double joint_likeli3(gsl_vector* p, double liknorm3) {

  	double current_proposed_k;

  	gsl_vector* zwischen = gsl_vector_calloc(n_all_parameters);
  	gsl_blas_dsymv(CblasLower, -0.5, supsig, p, 0.0, zwischen);
  	gsl_blas_ddot(p, zwischen, &current_proposed_k);


  	double help = current_proposed_k - liknorm3;
  	gsl_vector_free(zwischen);
  	return(help);
  }

  //for likeli4 and likeli5 see lkj.cpp

  //derivatives of joint likelihood by diffusion-model parameters; phase > 2
  void dhudwien2(int* nips, gsl_vector* hampar, double* tavw,  double* alltaus, double* dstore, gsl_vector* dhampar) {


  	gsl_vector_view t1 = gsl_vector_subvector(dhampar, 0, (indi + igroup) * icompg);
  	gsl_vector_set_zero(&t1.vector);

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

  				double tausum = 0.0;
  				double da = 0.0, dv = 0.0, dw = 0.0;
  				int nntim = dNNODES(t, im);
  				for (int i = 0; i != nntim; i++) for (int pm = 0; pm != 2; pm++) {
  					dstore[jj] = dwiener_d(alltaus[jj], aa, va, wa, accuracy * 1.2);

  					if (dCOMP(0, iaa)) da -= dadwiener_d(alltaus[jj], aa, va, wa, dstore[jj]);// gsl_vector_set(dhampar, dmapAVW(t, 0, iaa), gsl_vector_get(dhampar, dmapAVW(t, 0, iaa)) - dadwiener_d(alltaus[jj], aa, va, wa, d));
  					if (dCOMP(2, iww)) dw -= dwdwiener_d(alltaus[jj], aa, va, wa, dstore[jj]); // gsl_vector_set(dhampar, dmapAVW(t, 2, iww), gsl_vector_get(dhampar, dmapAVW(t, 2, iww)) - dwdwiener_d(alltaus[jj], aa, va, wa, d));

  					if (dCOMP(1, ivv)) tausum += fabs(alltaus[jj]);
  					jj++;
  				}
  				if (dCOMP(1, ivv)) {
  					double temp = -aa * (2 * wa - 1.0) * dNNODES(t, im) - va * tausum;
  					dv -= temp;
  					//					gsl_vector_set(dhampar, dmapAVW(t, 1, ivv), gsl_vector_get(dhampar, dmapAVW(t, 1, ivv)) - temp);
  				}
  				for (int pm = 0; pm != 2; pm++) {
  					int n = dNIPS(t, pm, im);
  					if (n != 0) {
  						double dav;
  						if ((dCOMP(0, iaa)) || (dCOMP(1, ivv))) {
  							dav = davlogprob_upperbound(pm, aa, va, wa);
  							if (dCOMP(0, iaa)) da += n * dalogprob_upperbound(pm, aa, va, wa, dav); // gsl_vector_set(dhampar, dmapAVW(t, 0, iaa), gsl_vector_get(dhampar, dmapAVW(t, 0, iaa)) + n * dalogprob_upperbound(pm, aa, va, wa, dav));
  							if (dCOMP(1, ivv)) dv += n * dvlogprob_upperbound(pm, aa, va, wa, dav); //gsl_vector_set(dhampar, dmapAVW(t, 1, ivv), gsl_vector_get(dhampar, dmapAVW(t, 1, ivv)) + n * dvlogprob_upperbound(pm, aa, va, wa, dav));
  						}
  						if (dCOMP(2, iww)) dw += n * dwlogprob_upperbound(pm, aa, va, wa); //gsl_vector_set(dhampar, dmapAVW(t, 2, iww), gsl_vector_get(dhampar, dmapAVW(t, 2, iww)) + n * dwlogprob_upperbound(pm, aa, va, wa));
  					}
  				}
  				int ima = dmapAVW(t, 0, iaa), imv = dmapAVW(t, 1, ivv), imw = dmapAVW(t, 2, iww);
  				if (dCOMP(0, iaa)) gsl_vector_set(dhampar, ima, gsl_vector_get(dhampar, ima) + da);
  				if (dCOMP(1, ivv)) gsl_vector_set(dhampar, imv, gsl_vector_get(dhampar, imv) + dv);
  				if (dCOMP(2, iww)) gsl_vector_set(dhampar, imw, gsl_vector_get(dhampar, imw) + dw);
  			}
  		else for (int t = 0; t != indi; t++) jj += 2 * dNNODES(t, im);
  	}

  	gsl_vector* temp = gsl_vector_alloc(indi * icompg); jj = 0;
  	for (int t = 0; t != indi; t++)
  		for (int type = 0; type != 3; type++) {
  			int ift = ifree[type];
  			for (int ip = 0; ip != ift; ip++) if (dCOMP(type, ip)) {
  				gsl_vector_set(temp, jj++, dlogit(avwtrans[type], invlogit(avwtrans[type], dTAVW(t, type, ip))));
  			}
  		}

  	gsl_vector_view t2 = gsl_vector_subvector(dhampar, iavwoff, indi * icompg);

  	gsl_vector_mul(&t2.vector, temp);

  	gsl_vector_free(temp);

  	gsl_vector_view tm2 = gsl_vector_subvector(hampar, 0, igroup * icompg);
  	gsl_vector_view dtm2 = gsl_vector_subvector(dhampar, 0, igroup * icompg);

  	for (int t = 0; t != indi; t++) {
  		gsl_vector_view dgr = gsl_vector_subvector(dhampar, t2group[t] * icompg, icompg);
  		gsl_vector_view dindi = gsl_vector_subvector(dhampar, iavwoff + t * icompg, icompg);
  		gsl_vector_add(&dgr.vector, &dindi.vector);
  	}

  	gsl_blas_daxpy(PRIOR, &tm2.vector, &dtm2.vector);
  	// what is missing from dhampar for person-specific deviations is: 1. Chain-rule for transformation using Cholesky-factor of Sigma 2. standard-normal prior for transformed person-specific deviations
  	// These contributions can be found in void dmvnlkjdy
  }


  //derivatives of model likelihood by motor-time related parameters
  void   dhudlam2(const std::vector<trial> & daten, double* rest, gsl_vector* hampar, double* tlams, double* explambda, double omega, gsl_vector* dhampar) {

  	double* dlam = (double*)calloc(indi * (respno + 1), sizeof(double));
  	double* dhrmu = (double*)calloc(igroup * respno, sizeof(double));

  	for (int x = 0; x != datenzahl; x++) {
  		int t = daten[x].person;
  		int r = cat2resp[daten[x].category];
  		int itr = t * respno + r, igr = t2group[t] * respno + r;
  		double xa = tlams[itr] - rest[x];
  		double xx = xa / (1.0 + gsl_pow_2(xa / explambda[t]) / degf);
  		dlam[itr] += xx;
  		dlam[indi * respno + t] -= xa * xx;
  	}

  	double dings = (degf + 1.0) / (1.0 * degf);
  	for (int t = 0; t != indi; t++) {
  		int it = indi * respno + t;
  		double sig = explambda[t];
  		dlam[it] *= dings;

  		dlam[it] -= priordf * omega;
  		dlam[it] /= gsl_pow_3(sig);
  		dlam[it] += (priordf + 1) / sig;
  		for (int r = 0; r != respno; r++) {
  			int itr = t * respno + r, igr = t2group[t] * respno + r;
  			dlam[itr] *= dings / gsl_pow_2(sig);
  			double mu = tlams[itr];

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
  			dhrmu[igr] += (gsl_vector_get(hampar, irmuoff + igr) - mu_prior) / gsl_pow_2(rsd);
  		}
  	}

  	gsl_vector_view t3 = gsl_vector_subvector(dhampar, ilamoff, indi * (respno + 1));
  	gsl_vector_view t4 = gsl_vector_view_array(dlam, indi * (respno + 1));
  	gsl_vector_memcpy(&t3.vector, &t4.vector);
  	gsl_vector_view t5 = gsl_vector_subvector(dhampar, irmuoff, igroup * respno);
  	gsl_vector_view t6 = gsl_vector_view_array(dhrmu, igroup * respno);
  	gsl_vector_memcpy(&t5.vector, &t6.vector);

  	// Korrektur wegen Jacoby-Faktor
  	//kann mit blas beschleunigt werden
  	for (int t = 0; t != indi; t++) gsl_vector_set(dhampar, isigoff + t, (-1.0) + gsl_vector_get(dhampar, isigoff + t) * explambda[t]);
  	free(dlam); free(dhrmu);
  	// what is missing from dhampar for person-specific deviations is: 1. Chain-rule for transformation using Cholesky-factor of Sigma 2. standard-normal prior for transformed person-specific deviations
  // These contributions can be found in void dmvnlkjdy
  }

  // dhudext see lkj.cpp

  //one leapfrog cycle
  void Leapfrog2(int* nips, gsl_vector* hampar, std::vector<double>& zt, std::vector<double>& zr, gsl_matrix* wt, gsl_matrix* wr, double* tavw, double* tlams, gsl_vector* dhampar, const std::vector<trial>& daten, double* explambdas,  double* alltaus, double* dstore, double* rest, double& omega, double eps, gsl_vector* p) {

  	gsl_blas_daxpy(-eps / 2, dhampar, p);

  	gsl_blas_dsymv(CblasLower, eps, supsig, p, 1.0, hampar);

  	make_tavwtlams(0, hampar, zt, wt, tavw);
  	make_tavwtlams(1, hampar, zr, wr, tlams);

  	for (int t = 0; t != indi; t++) explambdas[t] = exp(gsl_vector_get(hampar, isigoff + t));

  	omega = exp(gsl_vector_get(hampar, n_all_parameters - 1));

  	dhudwien2(nips, hampar, tavw, alltaus, dstore, dhampar);
  	dhudlam2(daten, rest, hampar, tlams, explambdas, omega, dhampar);
  	dhudext(hampar, explambdas, zt, zr, wt, wr, etat, etar, dhampar);

  	gsl_blas_daxpy(-eps / 2, dhampar, p);
  }



  double step0(int* nips, gsl_vector* dhampar, const std::vector<trial> & daten,  double* rest,
  	 double* alltaus, struct Theta* theta, gsl_vector* p, double u, int v, int j, double eps, int& n, int& s, double liknorm[6]) {

  	double* dstore = 0; if (!(dstore = (double*)malloc(ntau * sizeof(double)))) { Rprintf("Allocation failure2\n"); }
  	double* loglambdas = theta->loglambda;
  	double* tavw = theta->tavw;
  	double* tlams = theta->tlams;
  	gsl_vector* hampar = theta->hampar;
  	double omega;
  	std::vector<double> zt;
  	std::vector<double> zr;
  	gsl_matrix* wt = gsl_matrix_alloc(icompg, icompg);
  	gsl_matrix* wr = gsl_matrix_alloc(respno, respno);

  	Leapfrog2(nips, hampar, zt, zr, wt, wr, tavw, tlams, dhampar, daten,  loglambdas, alltaus, dstore, rest, omega, v * eps, p);

  	double temp = joint_likelihood2(nips, hampar, tavw, alltaus, dstore, liknorm[0]) +
  		rjoint_likelihood2(daten, rest, hampar, tlams, loglambdas, omega, liknorm[1]) +
  		joint_likeli3(p, liknorm[2]) +
  		joint_likeli4(0, hampar, zt, wt, etat, taut, liknorm[3]) +
  		joint_likeli4(1, hampar, zr, wr, etar, taur, liknorm[4]) +
  		joint_likeli5(hampar, loglambdas, liknorm[5]);
  	n = (u <= temp) ? 1 : 0;
  	s = (u - 1000.0 < temp) ? 1 : 0;
  	free(dstore);
  	gsl_matrix_free(wt);
  	gsl_matrix_free(wr);
  	return temp;
  }

  //one Hamiltonian cycle
  double inner_product2(gsl_vector* p, gsl_vector* hamparp, gsl_vector* hamparm) {
  	gsl_vector* stacked = gsl_vector_alloc(n_all_parameters);
  	gsl_vector_memcpy(stacked, hamparp);

  	gsl_vector_sub(stacked, hamparm);

  	double inner;
  	gsl_blas_ddot(stacked, p, &inner);
  	gsl_vector_free(stacked);
  	return inner;
  }

  //Build tree in NUTs algorithm
  void buildtree2(int* nips, const std::vector<trial> & daten, double* rest,
  	 double* alltaus, struct Theta* theta, struct Theta* thetadash, gsl_vector* dhampar, gsl_vector* p,
  	double u, int v, int j, double eps, gsl_rng* rst, double liknorm[6], int& ndash, int& sdash, int& nadash, double& adash, bool adapt) {

  	std::stack<struct Node*> stack;
  	int n, s = 1;
  	Node* root = (j < 13) ? trees[j] : make_tree2(j);
  	int jp1 = j + 1;
  	std::vector<struct store> speicher; for (int i = 0; i != jp1; i++) speicher.push_back(newstore());

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
  				double temp = step0(nips, dhampar, daten, rest, alltaus, theta, p, u, v, j, eps, n, s, liknorm);
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
  					int i = root->index;
  					int ndd = speicher[i + 1].n, nd = speicher[i].n;
  					if ((ndd > 0) && ((ndd + nd) * oneuni(rst) <= ndd)) {
  						thetacopy(speicher[i].thetadash, speicher[i + 1].thetadash);
  					}
  					speicher[i].s = speicher[i + 1].s;
  					if (speicher[i].s == 1) {
  						if ((v * inner_product2(speicher[i].p, theta->hampar, speicher[i].theta->hampar) < 0.0) ||
  							(v * inner_product2(p, theta->hampar, speicher[i].theta->hampar) < 0.0))
  						{
  							speicher[i].s = 0; s = 0;
  						}
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

  //the sampler; phase > 2
  bool hnuts2(int* nips, gsl_vector* hampar, double* tavw, double* tlams,  const std::vector<trial> & daten, double* rest, double* explambdas,
  	 double* alltaus, double liknorm[6], double& activeeps, double& epsm, double& Hobjective, int m, bool save, gsl_rng* rst) {

  	double lamtest = gsl_vector_get(hampar, 0);
  	int interval = std::max(PHASE2, 5 * n_all_parameters);
  //	interval = 100;
  	interval = (interval / ireps + 1) * ireps;
  	m = ((m-1) % interval) + 1;
  	bool adapt = ((m <= PHASE1) && (!(save)) && (phase == 3));

  	double* dstore = 0; if (!(dstore = (double*)malloc(ntau * sizeof(double)))) { Rprintf("Allocation failure2\n"); }
  	gsl_vector* p = gsl_vector_alloc(n_all_parameters);
  	gsl_vector* pp = gsl_vector_alloc(n_all_parameters);
  	gsl_vector* pm = gsl_vector_alloc(n_all_parameters);

  	gsl_vector* dhamparp = gsl_vector_calloc(n_all_parameters);
  	gsl_vector* dhamparm = gsl_vector_calloc(n_all_parameters);

  	gsl_matrix* wt = gsl_matrix_calloc(icompg, icompg);
  	gsl_matrix* wr = gsl_matrix_calloc(respno, respno);
  	std::vector<double> zt;
  	std::vector<double> zr;

  	make_tavwtlams(0, hampar, zt, wt, tavw);
  	make_tavwtlams(1, hampar, zr, wr, tlams);
  	double omega = exp(gsl_vector_get(hampar, n_all_parameters - 1));

  //	test0(nips,  dhamparp, daten, rest,
  //		 alltaus, hampar, pp, 0.001,  liknorm);

  	dhudwien2(nips, hampar, tavw, alltaus, dstore, dhamparp);
  	dhudlam2(daten, rest, hampar, tlams, explambdas, omega, dhamparp);
  	dhudext(hampar, explambdas, zt, zr, wt, wr, etat, etar, dhamparp);
  	gsl_vector_memcpy(dhamparm, dhamparp);

  	for (int i = 0; i != n_all_parameters; i++) gsl_vector_set(p, i, onenorm(rst));

  	gsl_blas_dtrmv(CblasLower, CblasTrans, CblasNonUnit, sigisqrt, p);

  	liknorm[0] += joint_likelihood2(nips, hampar, tavw, alltaus, dstore, liknorm[0]);
  	free(dstore);

  	liknorm[1] += rjoint_likelihood2(daten, rest, hampar, tlams, explambdas, omega, liknorm[1]);

  	liknorm[2] += joint_likeli3(p, liknorm[2]);

  	liknorm[3] += joint_likeli4(0, hampar, zt, wt, etat, taut,liknorm[3]);

  	liknorm[4] += joint_likeli4(1, hampar, zr, wr, etar, taur, liknorm[4]);

  	liknorm[5] += joint_likeli5(hampar, explambdas, liknorm[5]);

  	Theta* thetap = newTheta();
  	Theta* thetam = newTheta();
  	Theta* thetadash = newTheta();
  	Theta* thetac = (struct Theta*)malloc(sizeof(struct Theta));
  	thetac->hampar = hampar; thetac->tavw = tavw; thetac->loglambda = explambdas; thetac->tlams = tlams;
  	thetacopy(thetap, thetac); thetacopy(thetam, thetac);
  	int n, s;
  	pcopy(pp, p); pcopy(pm, p);
  	n = 1;
  	s = 1;

  	int na; double a;

  	double u = log(oneuni(rst));
  	double epstemp = activeeps;
  	int j = 0;
  	while (s == 1) {
  		int ndash, sdash;
  		int v = (oneuni(rst) <= 0.5) ? -1 : 1;
  		if (v == -1) {
  			buildtree2(nips,daten, rest, alltaus,
  				thetam, thetadash, dhamparm, pm, u, v, j, epstemp, rst, liknorm, ndash, sdash, na, a, adapt);
  		}
  		else {
  			buildtree2(nips, daten, rest, alltaus,
  				thetap, thetadash, dhamparp, pp, u, v, j, epstemp, rst, liknorm, ndash, sdash, na, a, adapt);
  		}
  		if ((sdash == 1) && (n * oneuni(rst) <= ndash))
  		{
  			thetacopy(thetac, thetadash);
  		}
  		n += ndash;
  		j = j + 1;
  		if (j == (phase<=3?maxtreedepth1_3:maxtreedepth4)) sdash = 0.0;
  		s = sdash;
  		if (s == 1)
  			if ((inner_product2(pm, thetap->hampar, thetam->hampar) < 0.0) ||
  				(inner_product2(pp, thetap->hampar, thetam->hampar) < 0.0))
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
  //	std::cout << setw(5) << j << setw(5) << n;
  	free(thetac);
  	remove_Theta(thetam);
  	remove_Theta(thetap);
  	remove_Theta(thetadash);
  	gsl_vector_free(p); gsl_vector_free(pp); gsl_vector_free(pm);
  	gsl_vector_free(dhamparp);
  	gsl_vector_free(dhamparm);
  	gsl_matrix_free(wt);
  	gsl_matrix_free(wr);

  	return(gsl_vector_get(hampar, 0) != lamtest);
  }

  //estimate variance-covariance matrix of posterior distribution of parameters
  void make_supersigs(int anz, double* parmonstore, gsl_matrix* supsig, gsl_matrix* sigisqrt) {
  	// pool
  	gsl_matrix_view ssig = gsl_matrix_view_array(supersig, NOTHREADS, n_all_parameters * n_all_parameters);
  	gsl_vector* ones = gsl_vector_alloc(NOTHREADS);
  	gsl_vector* temp = gsl_vector_alloc(n_all_parameters*n_all_parameters);
  	gsl_vector_set_all(ones, 1.0/(anz*NOTHREADS));
  	gsl_blas_dgemv(CblasTrans, 1.0, &ssig.matrix, ones, 0.0, temp);

  	gsl_vector* sums = gsl_vector_alloc(n_all_parameters);

  	gsl_vector_set_zero(sums);
  	gsl_vector_view parmst = gsl_vector_view_array(parmonstore, NOTHREADS * 2 * n_all_parameters);
  	for (int k = 0; k != NOTHREADS; k++) {
  		gsl_vector_view parmstk = gsl_vector_subvector(&parmst.vector, k * 2 * n_all_parameters, n_all_parameters);
  		gsl_vector_add(sums, &parmstk.vector);
  	}
  	gsl_vector* devs = gsl_vector_alloc(n_all_parameters);
  	gsl_matrix_view temps = gsl_matrix_view_vector(temp, n_all_parameters, n_all_parameters);
  	gsl_matrix_memcpy(supsig, &temps.matrix);

  	for (int k = 0; k != NOTHREADS; k++) {
  		gsl_vector_view parmstk = gsl_vector_subvector(&parmst.vector, k * 2 * n_all_parameters, n_all_parameters);
  		gsl_vector_memcpy(devs, &parmstk.vector);
  		gsl_blas_daxpy(-1.0 / NOTHREADS, sums, devs);
  		gsl_blas_dsyr(CblasLower, 1.0/NOTHREADS, devs, supsig);
  	}
  	// probably superfluous
  	for (int i = 0; i != n_all_parameters; i++)
  		for (int j = 0; j <= i; j++) {
  			double xtemp = gsl_matrix_get(supsig, i, j);
  	//		if ((i< icompg) && (i == j)) xtemp *= 1.05;
  			gsl_matrix_set(supsig, j, i, xtemp);

  		}
  	gsl_vector_free(ones);
  	gsl_vector_free(temp);
  	gsl_vector_free(sums);
  	gsl_vector_free(devs);

  	gsl_matrix_memcpy(sigisqrt, supsig);
  	gsl_linalg_cholesky_decomp1(sigisqrt);
  	gsl_linalg_tri_lower_invert(sigisqrt);
  }

}
