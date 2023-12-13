// authors: Christoph Klauer and Raphael Hartmann
#include "rts.h"

namespace ertmpt {
  
  bool restart;
  std::vector<trial> itdaten;
  
  void make_pij_for_one_trial_new_new(trial one, double *x, double *pij, double &pj) {
  
  #define DRIN(J,K,X) drin[J*zweig*nodemax+K*nodemax + X]
  #define NDRIN(J,K) ndrin[J*zweig+K]
  #define AR(I,J,R) ar[I*zweig*nodemax + J*nodemax + R]
  #define TREE_AND_NODE2PAR(ITREE,R) tree_and_node2par[ITREE*nodemax+R]
  
  	double d0ha;
  	int j = one.category, itree = one.tree;
  
  
  	pj = 0.0;
  	for (int k = 0; k != branch[j]; k++) {
  		//			 pij[k]=1.0;
  		for (int ir = 0; ir != NDRIN(j, k); ir++) {
  			int r = DRIN(j, k, ir); int ix = AR(j, k, r); int ip = TREE_AND_NODE2PAR(itree, r);
  			//			 d0ha=(ix>0)?x[ip]:(1.0-x[ip]);
  			if (comp[ip]) {
  				int iz = kern2free[ip];
  				d0ha = (ix > 0) ? lnnorm(x[iz]) : lnnorm(-x[iz]);
  			}
  			else {
  				d0ha = (ix > 0) ? lnnorm(consts[ip]) : lnnorm(-consts[ip]);
  			}
  			pij[k] = pij[k] + (d0ha);
  		}
  		pj = (k == 0) ? pij[0] : logsum(pj, pij[k]);
  		if (!(pj == pj) || !(std::isfinite(pj))) {
  			if (DEBUG) Rprintf("pj is %g\n", pij[k]);
  			pj = -sqrt(DBL_MAX);
  			restart = true;
  		}
  	}
  }
  /*
  void make_tij_for_one_trial_new_new(trial one, double *lambdas, double rmu, double rsig, double xsi, double *pij) {
  
  	int itree = one.tree; int j = one.category; double rt = one.rt / 1000.0;
  
  
  	for (int k = 0; k != branch[j]; k++) {
  		int pfadlength = NDRIN(j, k);
  
  		double *lams = 0; lams = (double *)malloc(pfadlength * sizeof(double));
  		int complength = 0;
  		for (int ir = 0; ir != pfadlength; ir++) {
  			int r = DRIN(j, k, ir); int ip = TREE_AND_NODE2PAR(itree, r); int pm = (AR(j, k, r) > 0) ? 1 : 0;
  			if (comp[ip + (1 + pm)*kernpar]) {
  				int iz = kern2free[(1 + pm)*kernpar + ip] - ifree; lams[complength] = lambdas[iz];
  				complength++;
  			}
  		}
  		if (complength == 0) { pij[k] = (-0.5*gsl_pow_2((rt - rmu) / rsig)) / sqr2pi / rsig - xsi; }
  		if (complength == 1) {
  			double lam = lams[0]; pij[k] = (logexgaussian(lam, rmu, rsig, rt) - xsi);
  
  		}
  		if (complength >= 2) {
  			double *loglams = 0; loglams = (double *)malloc(complength * sizeof(double));
  
  			double hypo_plus = 0, hypo_minus = 0; bool f_plus = true, f_minus = true;
  			for (int ir = 0; ir != complength; ir++) loglams[ir] = log(lams[ir]);
  
  			for (int ir = 0; ir != complength; ir++) {
  				int sign = 1; double factor = 0.0;
  				for (int ij = 0; ij != complength; ij++) if (ij != ir) {
  					if (lams[ij] - lams[ir] < 0) sign = -sign;
  					factor = factor + loglams[ij] - log(fabs(lams[ij] - lams[ir]));
  				}
  				double temp = factor + logexgaussian(lams[ir], rmu, rsig, rt);
  				if (sign > 0) { if (f_plus) { hypo_plus = temp; f_plus = false; } else hypo_plus = logsum(hypo_plus, temp); }
  				if (sign < 0) { if (f_minus) { hypo_minus = temp; f_minus = false; } else hypo_minus = logsum(hypo_minus, temp); }
  			}
  			double dif = hypo_minus - hypo_plus;
  			if (dif < 0) 	pij[k] = hypo_plus + gsl_log1p(-exp(hypo_minus - hypo_plus)) - xsi; else pij[k] = -1.0e10;
  
  
  			if (hypo_plus < hypo_minus) {
  				if (DEBUG) Rprintf("pij[k] < 0; pf = 2, pij[k]: %g hm: %g hp: %g rsig: %g rmu: %g\n", pij[k], hypo_minus, hypo_plus, rsig, rmu);
  				//std::cout << "pij[k] < 0; pf = 2, pij[k]: " << pij[k] << " " << "hm: " << hypo_minus << " hp: " << hypo_plus << " " << " rsig: " << rsig << " rmu: " << rmu << std::endl;
  				if (DEBUG) {for (int ir = 0; ir != pfadlength; ir++) Rprintf("%.4g ", loglams[ir]);}
  				if (DEBUG) Rprintf("%g\n", xsi);
  				pij[k] = -sqrt(DBL_MAX);
  				restart = true;
  
  			}
  			if (loglams) free(loglams);
  		}
  
  		free(lams);
  		if (DEBUG) {if ((pfadlength > 5) || (pfadlength == 0)) Rprintf("pfadlength\n");}
  	}
  
  }*/
  
  #define PFAD_INDEX(J,K) pfad_index[J*zweig+K]
  
  void make_tij_for_repetitions(trial one, double *lambdas, double rmu, double rsig, double xsi, double *pij) {
  
  	/*int itree = one.tree;*/ int j = one.category; double rt = one.rt / 1000.0;
  
  
  	for (int k = 0; k != branch[j]; k++) {
  		int pfadlength = NDRIN(j, k);
  
  		double *lams = 0; lams = (double *)malloc(pfadlength * sizeof(double));
  		int complength = 0;
  		/*if (PFAD_INDEX(j,k)==-1) {
  			for (int ir = 0; ir != pfadlength; ir++) {
  				int r = DRIN(j, k, ir); int ip = TREE_AND_NODE2PAR(itree, r); int pm = (AR(j, k, r) > 0) ? 1 : 0;
  				if (comp[ip + (1 + pm)*kernpar]) {
  					int iz = kern2free[(1 + pm)*kernpar + ip] - ifree; lams[complength] = lambdas[iz];
  					complength++;
  				}
  			}
  		}
  		else
  		{*/
  		int ipfad = PFAD_INDEX(j, k);
  		pfadinfo akt_pfad = path_info[ipfad];
  		complength = akt_pfad.a;
  		for (int ir = 0; ir != akt_pfad.a; ir++) {
  			int ip = akt_pfad.pfad_par[ir];
  			int pm = akt_pfad.pm[ir];
  			int iz = kern2free[(1 + pm)*kernpar + ip] - ifree; lams[ir] = lambdas[iz];
  		}
  		// }
  		if (complength == 0) { pij[k] = (-0.5*gsl_pow_2((rt - rmu) / rsig)) / sqr2pi / rsig - xsi; }
  		// if ((complength == 1) && (PFAD_INDEX(j,k) == -1)) {
  		// 	double lam = lams[0]; pij[k] = (logexgaussian(lam, rmu, rsig, rt) - xsi);
  		// }
  		if ((complength == 1) /*&& (PFAD_INDEX(j, k) > -1)*/) {
  			// printf("in ==1\n");
  			// int ipfad = PFAD_INDEX(j, k);
  			// pfadinfo akt_pfad = path_info[ipfad];
  			// if (complength != akt_pfad.a) Rprintf("complength != a");
  			double lam = lams[0];
  			double hplus, hminus;
  			loggammagaussian(akt_pfad.r[0] - 1, lam, rmu, rsig, rt, hplus, hminus);
  			double temp = logdiff(hplus,hminus) +log(lam) * akt_pfad.r[0];
  			if (temp == GSL_NEGINF) {
  				pij[k] = -sqrt(DBL_MAX);
  				restart = true;
  			}
  			else pij[k] = temp - xsi;
  			// printf("end ==1\n");
  		}
  /*		if ((complength >= 2) && (PFAD_INDEX(j, k) == -1)) {
  			double *loglams = 0; loglams = (double *)malloc(complength * sizeof(double));
  
  			double hypo_plus = 0, hypo_minus = 0; bool f_plus = true, f_minus = true;
  			for (int ir = 0; ir != complength; ir++) loglams[ir] = log(lams[ir]);
  
  			for (int ir = 0; ir != complength; ir++) {
  				int sign = 1; double factor = 0.0;
  				for (int ij = 0; ij != complength; ij++) if (ij != ir) {
  					if (lams[ij] - lams[ir] < 0) sign = -sign;
  					factor = factor + loglams[ij] - log(fabs(lams[ij] - lams[ir]));
  				}
  				double temp = factor + logexgaussian(lams[ir], rmu, rsig, rt);
  				if (sign > 0) { if (f_plus) { hypo_plus = temp; f_plus = false; } else hypo_plus = logsum(hypo_plus, temp); }
  				if (sign < 0) { if (f_minus) { hypo_minus = temp; f_minus = false; } else hypo_minus = logsum(hypo_minus, temp); }
  			}
  			double dif = hypo_minus - hypo_plus;
  			if (dif < 0) 	pij[k] = hypo_plus + gsl_log1p(-exp(hypo_minus - hypo_plus)) - xsi; else pij[k] = -1.0e10;
  
  
  			if (hypo_plus < hypo_minus) {
  				if (DEBUG) Rprintf("pij[k] < 0; pf = 2, pij[k]: %g hm: %g hp: %g rsig: %g rmu: %g\n", pij[k], hypo_minus, hypo_plus, rsig, rmu);
  				//std::cout << "pij[k] < 0; pf = 2, pij[k]: " << pij[k] << " " << "hm: " << hypo_minus << " hp: " << hypo_plus << " " << " rsig: " << rsig << " rmu: " << rmu << std::endl;
  				if (DEBUG) {for (int ir = 0; ir != pfadlength; ir++) Rprintf("%.4g ", loglams[ir]);}
  				if (DEBUG) Rprintf("%g\n", xsi);
  				pij[k] = -sqrt(DBL_MAX);
  				restart = true;
  
  			}
  			if (loglams) free(loglams);
  		}*/
  		if ((complength >= 2) /*&& (PFAD_INDEX(j, k) > -1)*/) {
  			// printf("in >= 2\n");
  			double *loglams = 0; loglams = (double *)malloc(complength * sizeof(double));
  			for (int ir = 0; ir != complength; ir++) loglams[ir] = log(lams[ir]);
  			// int ipfad = PFAD_INDEX(j, k);
  			// pfadinfo akt_pfad = path_info[ipfad];
  			double temp = logf_tij(akt_pfad.a, akt_pfad.r, lams,loglams, rmu, rsig, rt);
  			if (temp == GSL_NEGINF) {
  				pij[k] = -sqrt(DBL_MAX);
  				restart = true;
  			}
  			else pij[k] = (temp)-xsi;
  			if (loglams) free(loglams);
  			// printf("end >= 2\n");
  		}
  
  		free(lams);
  	}
  
  }
  
  void trans(int n, double *x, double *pars, bool inverse) {
  	if (!(inverse))
  		for (int i = 0; i != n; i++) x[i] = pars[1 + i] + (pars[1 + n + i] - pars[1 + i]) / (1 + exp(-x[i]));
  	else
  		for (int i = 0; i != n; i++) x[i] = -log((pars[1 + n + i] - x[i]) / (x[i] - pars[1 + i]));
  }
  
  
  double objfun(const gsl_vector * y, void * params)
  {
  	double *pars = (double *)params;
  	int n = static_cast<int>(trunc(pars[0]));
  	double *x = 0; x = (double *)malloc(n * sizeof(double));
  	for (int i = 0; i != n; i++) x[i] = gsl_vector_get(y, i);
  	trans(n, x, pars, false);
  	if (DEBUG) {if (!(x[1] == x[1])) Rprintf("x[1] is NaN\n");}
  
  
  	double *lambdas = 0; lambdas = (double *)malloc(ilamfree * sizeof(double));
  	double *x_for_all = 0; x_for_all = (double *)malloc(ifree * sizeof(double));
  	double *pij = 0; pij = (double *)malloc(zweig * sizeof(double));
  
  	int trialno = static_cast<int>(itdaten.size());
  
  	for (int i = 0; i != ilamfree; i++) {
  		if (DEBUG) if (!(x[ifree + i] == x[ifree + i])) {
  			Rprintf("lambdas[%d] is %g\n", ifree + i, x[ifree + i]);
  			//std::cout << "lambdas[" << ifree + i << "] is " << x[ifree + i] << std::endl;
  		}
  		lambdas[i] = x[ifree + i];
  	}
  	for (int i = 0; i != ifree; i++) x_for_all[i] = (x[i]);
  
  	double rmu = x[ifree + ilamfree], rsig = x[ifree + ilamfree + 1];
  	//   std::cout << rmu << " " << rsig << " " << x[12]; char xx; std::cin>> xx;
  	double xsi = lnnorm(rmu / rsig);
  	if (DEBUG) if (!(rmu == rmu)) Rprintf("rmu is %g\n", rmu);
  	double loglik = 0.0;
  
  	for (int x = 0; x != trialno; x++) {
  		trial one = itdaten[x];  one.person = 0;
  
  		make_tij_for_repetitions(one, lambdas, rmu, rsig, xsi, pij);
  		double p;
  		make_pij_for_one_trial_new_new(one, x_for_all, pij, p);
  		loglik += -2 * p;
  	}
  	if (!std::isfinite(loglik)) {
  		restart = true; loglik = -1.0e10;
  		if (DEBUG) printf("unfortunate\n");
  	}
  	else restart = false;
  	free(x);
  	free(pij);
  	free(lambdas);
  	free(x_for_all);
  
  	return(loglik);
  
  }
  
  
  
  void tby_individuals(std::vector<trial> daten, int kerntree, double *beta, double *lambdas, double *restpars, gsl_rng *rst) {
  #define BETA(T,I) beta[T*ifree+I]
  
  	double *pars = 0, *x = 0, *xsave = 0;
  	int n = ifree + ilamfree + 2;
  	x = (double *)malloc(n * sizeof(double));
  	xsave = (double *)malloc(n * sizeof(double));
  	pars = (double *)malloc((2 * n + 1) * sizeof(double));
  
  
  	double oldfit;
  
  
  	for (int ip = 0; ip != ifree; ip++) { pars[1 + ip] = -3; pars[1 + n + ip] = 3; }
  	for (int il = 0; il != ilamfree; il++) {
  		pars[1 + ifree + il] = 1.0;
  		pars[1 + n + ifree + il] = 200;
  	}
  	pars[1 + n - 2] = -1.0; pars[1 + n + n - 2] = 1.0;
  	pars[1 + n - 1] = 1.0e-3; pars[1 + n + n - 1] = 1.0;
  	pars[0] = n * 1.0;
  
  	int imax = 1000 * n;
  	double size;
  
  	if (!DEBUG) Rprintf("\nCalculating initial values:\n");
  
  	double progress = 0.0;
  	int ML_bar = 50;
  	if (!DEBUG) {
  		Rprintf("[");
  		for (int i = 0; i < ML_bar; i++) Rprintf(" ");
  		Rprintf("] 0%%");
  	}
  
  	for (int t = 0; t != indi; t++) {
  		restart = false;
  		itdaten.clear();
  		for (int ix = 0; ix != static_cast<int>(daten.size()); ix++) if (daten[ix].person == t) itdaten.push_back(daten[ix]);
  
  		progress = 1.0*(t+1)/indi;
  
  		double mean = 0.0, sx2 = 0.0; int ndat = static_cast<int>(itdaten.size());
  		for (int ix = 0; ix != ndat; ix++) {
  			double xd = itdaten[ix].rt / 1000.0;
  			mean += xd;
  			sx2 += gsl_pow_2(xd);
  		}
  		mean /= ndat;
  		sx2 /= ndat;
  		sx2 -= gsl_pow_2(mean); sx2 = sqrt(sx2);
  
  		oldfit = DBL_MAX;
  
  		for (int ix = 0; ix != 4; ix++) {
  		Vonvorn:
  			//	for (int i = 0; i < ifree; i++) x[i] = BETA(t, i) + 0.1*(oneuni(rst)-0.5);
  			for (int i = 0; i != ifree; i++) x[i] = oneuni(rst) - 0.5;
  			x[n - 2] = mean * 0.7 + 0.01*(oneuni(rst) - 0.5);
  			x[n - 1] = sx2 / 5;
  			for (int i = 0; i != ilamfree; i++) x[ifree + i] = pars[1 + ifree + i] + 10.0*oneuni(rst);
  
  			for (int i = 0; i != n; i++) {
  				x[i] = gsl_max(x[i], pars[1 + i] + 0.01*(1 + oneuni(rst)));
  				x[i] = gsl_min(x[i], pars[1 + n + i] - 0.01*(1 + oneuni(rst)));
  			}
  
  			trans(n, x, pars, true);
  
  			size_t iter = 0;
  			int status;
  			const gsl_multimin_fminimizer_type *T;
  			gsl_multimin_fminimizer *s;
  
  			gsl_vector *xx;
  			gsl_multimin_function my_func;
  
  			gsl_vector *step;
  			step = gsl_vector_alloc(n);
  			for (int i = 0; i != n; i++) gsl_vector_set(step, i, 0.1);
  
  			my_func.n = n;
  			my_func.f = objfun;
  			my_func.params = pars;
  
  			xx = gsl_vector_alloc(n);
  			for (int i = 0; i != n; i++) gsl_vector_set(xx, i, x[i]);
  
  			T = gsl_multimin_fminimizer_nmsimplex2;
  			s = gsl_multimin_fminimizer_alloc(T, n);
  
  
  			gsl_multimin_fminimizer_set(s, &my_func, xx, step);
  
  
  			do
  			{
  				iter++;
  				status = gsl_multimin_fminimizer_iterate(s);
  				if (status)
  					break;
  				size = gsl_multimin_fminimizer_size(s);
  				status = gsl_multimin_test_size(size, 1e-3);
  
  				if (restart)
  				{
  					gsl_vector_free(xx);
  					gsl_multimin_fminimizer_free(s);
  					gsl_vector_free(step);
  					if (DEBUG) Rprintf("Von vorn\n");
  					restart = false;
  					goto Vonvorn;
  				}
  
  				if (status == GSL_SUCCESS) {
  					if (DEBUG){
  						Rprintf("Minimum found at:\n");
  						Rprintf("%5d %10.5f\n", static_cast<int>(iter),
  							//			gsl_vector_get(s->x, 0),
  							//			gsl_vector_get(s->x, 1),
  							s->fval);
  					}
  				}
  
  				if (static_cast<int>(iter) == imax) {
  					if (DEBUG){
  						Rprintf("Minimum not found:\n");
  						Rprintf("%5d %10.5f\n", static_cast<int>(iter),
  							//			gsl_vector_get(s->x, 0),
  							//			gsl_vector_get(s->x, 1),
  							s->fval);
  					}
  				}
  			} while (status == GSL_CONTINUE && static_cast<int>(iter) < imax);
  
  			if (s->fval < oldfit) {
  				oldfit = s->fval;
  				for (int i = 0; i != n; i++) xsave[i] = gsl_vector_get(s->x, i);
  			}
  			gsl_vector_free(xx);
  			gsl_multimin_fminimizer_free(s);
  			gsl_vector_free(step);
  		}
  
  
  		trans(n, xsave, pars, false);
  
  		// int iz = 0;
  		for (int ip = 0; ip != ifree; ip++) { BETA(t, ip) = xsave[ip]; }
  		for (int ip = ifree; ip != ifree + ilamfree; ip++) lambdas[t*ilamfree + ip - ifree] = xsave[ip];
  		restpars[t + 3] = xsave[n - 2];
  		restpars[t + indi + 3] = gsl_pow_2(xsave[n - 1]);
  #pragma optimize("", off)
  		if (!DEBUG) {
  			int pos = static_cast<int>(ML_bar * progress);
  			Rprintf("\r[");
  			for (int i = 0; i < ML_bar; i++) {
  			  if (i < pos) {
  					Rprintf("=");
  				} else if (i == pos) {
  					Rprintf(">");
  				} else Rprintf(" ");
  			}
  			Rprintf("] %d%%", static_cast<int>(progress * 100.0));
  		}
  #pragma optimize("", on)
  
  	}
  	Rprintf("\n\n");
  
  	free(x);
  	free(xsave);
  	free(pars);
  }
  
}



namespace drtmpt {
  
  struct one_params {
    int n; //icompg + 2
    int m; //ifreeg + 2
    bool restart;
    std::vector<std::vector<double>> icdaten; //vector of response categories with entries given by vectors of reaction times observed for this category given a fixed person
  };
  
  
  void make_p_ind_cat(const std::vector<double> &rts, int j, double* x, std::vector<double>& ps) {
    int itree = cat2tree[j];
    std::vector < std::vector<double>> p; p.clear();
    double mu = x[ifreeg], sig = x[ifreeg + 1];
    
    int bj = branch[j];
    for (int k = 0; k != bj; k++) {
      int pfadlength = dNDRIN(j, k);
      double* a = (double*)malloc(pfadlength * sizeof(double));
      double* v = (double*)malloc(pfadlength * sizeof(double));
      double* w = (double*)malloc(pfadlength * sizeof(double));
      int* low_or_up = (int*)malloc(pfadlength * sizeof(int));
      for (int ir = 0; ir != pfadlength; ir++) {
        int r = dDRIN(j, k, ir);
        low_or_up[ir] = dAR(j, k, r);
        a[ir] = x[dTREE_AND_NODE2PAR(itree, r, 0)];
        v[ir] = x[ifree[0] + dTREE_AND_NODE2PAR(itree, r, 1)];
        w[ir] = x[ifree[0] + ifree[1] + dTREE_AND_NODE2PAR(itree, r, 2)];
      }
      
      std::vector<double> pbranch; pbranch.clear();
      convolution2(rts, pfadlength, low_or_up, a, v, w, mu, sig, pbranch);
      p.push_back(pbranch);
      free(a); free(v); free(w); free(low_or_up);
    }
    ps.clear();
    int rtss = rts.size();
    for (int x = 0; x != rtss; x++) {
      double temp = 0.0;
      for (int k = 0; k != bj; k++) temp += p[k][x];
      ps.push_back(log(temp));
    }
  }
  
  
  double objfun(const gsl_vector* y, void* params)
  {
    one_params* pars = (one_params*)params;
    int n = (pars->n), m = (pars->m);
    bool restart = (pars->restart);
    std::vector<std::vector<double>> icdaten = (pars->icdaten);
    double* x = 0; if (!(x = (double*)malloc(m * sizeof(double)))) { Rprintf("Allocation failure\n"); }
    //	for (int i = 0; i != n; i++) x[i] = gsl_vector_get(y, i);
    int jj = 0, iz = 0;
    for (int itype = 0; itype != 3; itype++) {
      int ift = ifree[itype];
      for (int ip = 0; ip != ift; ip++)
        if (dCOMP(itype, ip)) x[iz++] = logit(avwtrans[itype], gsl_vector_get(y, jj++));
        else x[iz++] = dCONSTS(itype, ip);
    }
    x[ifreeg] = gsl_vector_get(y, icompg);
    x[ifreeg + 1] = exp(gsl_vector_get(y, icompg + 1));
    
    if (!(x[1] == x[1])) Rprintf("x[1] is NaN\n");
    
    double mu = x[ifreeg], sig = x[ifreeg + 1];
    double xsi = log(gsl_cdf_tdist_P(mu / sig, degf) * sig);
    
    double loglik = 0.0;
    
    for (int j = 0; j != kerncat; j++) {
      std::vector<double> ps;
      if (icdaten[j].size() > 0) make_p_ind_cat(icdaten[j], j, x, ps);
      int icdjs = icdaten[j].size();
      for (int i = 0; i != icdjs; i++) loglik += -2 * ps[i];
      loglik += 2 * icdjs * xsi;
    }
    
    if (!isfinite(loglik)) { restart = true; loglik = -1.0e10; }
    else restart = false;
    pars->restart = restart;
    free(x);
    return(loglik);
  }
  
  
  
  void tby_individuals(const std::vector<trial>& daten, double* avw, double* lambdas, gsl_rng* rst) {
    
    double size;
    std::vector<double> temp; temp.clear();
    std::vector<std::vector<double>> itcdaten(indi * kerncat, temp);
    
    for (int x = 0; x != datenzahl; x++) {
      trial one = daten[x]; int t = one.person; int c = one.category;
      itcdaten[t * kerncat + c].push_back(one.rt / 1000.0);
    }
    
    int imax = 4000;
    
    double progress = 0.0;
    int ML_bar = 50;
    if (PROG_BAR_FLAG) {
      Rprintf("\nCalculating initial values:\n");
      
      Rprintf("[");
      for (int i = 0; i < ML_bar; i++) Rprintf(" ");
      Rprintf("] 0%%");
      
    }
    
    /* prepare threads */
    int NThreads = INIT_CPUs;
    int maxThreads = std::thread::hardware_concurrency();
    if (maxThreads == 0) Rprintf("Could not find out number of threads. Taking 2 threads.\n");
    int suppThreads = maxThreads == 0 ? 2 : maxThreads;
    int AmntOfThreads = suppThreads > NThreads ? NThreads : suppThreads;
    int NperThread = indi / AmntOfThreads;
    AmntOfThreads = AmntOfThreads > indi ? indi : AmntOfThreads;
    NperThread = AmntOfThreads == indi ? 1 : NperThread;
    std::vector<std::thread> threads(AmntOfThreads-1);
    
    /* starting threads while ... */
    for (int ithread = 0; ithread < AmntOfThreads-1; ithread++) {
      threads[ithread] = std::thread([&, ithread]() {
        for (int t = ithread*NperThread; t < (ithread+1)*NperThread; t++) {
          double oldfit = GSL_POSINF;
          
          gsl_vector* xx = gsl_vector_alloc(icompg + 2);
          gsl_vector* xsave = gsl_vector_alloc(icompg + 2);
          
          bool restart = false;
          std::vector<std::vector<double>> icdaten; icdaten.clear();
          int ndat = 0;
          for (int j = 0; j != kerncat; j++) {
            icdaten.push_back(itcdaten[t * kerncat + j]);
            ndat += static_cast<int>(itcdaten[t * kerncat + j].size());
          }
          
          int n = icompg + 2;
          int m = ifreeg + 2;
          one_params pars = { n, m, restart, icdaten };
          
          double mean = 0.0, sx2 = 0.0;
          for (int j = 0; j != kerncat; j++) {
            int icdjs = icdaten[j].size();
            for (int ix = 0; ix != icdjs; ix++) {
              double xd = icdaten[j][ix];
              mean += xd;
              sx2 += gsl_pow_2(xd);
            }
          }
          mean /= ndat;
          sx2 /= ndat;
          sx2 -= gsl_pow_2(mean); sx2 = sqrt(sx2);
          
          for (int ix = 0; ix != 1; ix++) {
            Vonvorn:
            gsl_vector* xx = gsl_vector_alloc(icompg + 2);
            for (int i = 0; i != icompg; i++) gsl_vector_set(xx, i, oneuni(rst) - 0.5);
            gsl_vector_set(xx, icompg, mean * 0.7 + 0.01 * (oneuni(rst) - 0.5));
            gsl_vector_set(xx, icompg + 1, log(sx2 / 5.0));
            
            size_t iter = 0;
            int status;
            const gsl_multimin_fminimizer_type* T;
            gsl_multimin_fminimizer* s;
            
            gsl_multimin_function my_func;
            
            gsl_vector* step;
            step = gsl_vector_alloc(icompg + 2);
            for (int i = 0; i != icompg + 2; i++) gsl_vector_set(step, i, 0.1);
            
            my_func.n = n;
            my_func.f = objfun;
            my_func.params = &pars;
            
            T = gsl_multimin_fminimizer_nmsimplex2;
            s = gsl_multimin_fminimizer_alloc(T, n);
            
            gsl_multimin_fminimizer_set(s, &my_func, xx, step);
            
            do
            {
              iter++;
              status = gsl_multimin_fminimizer_iterate(s);
              if (status)
                break;
              size = gsl_multimin_fminimizer_size(s);
              status = gsl_multimin_test_size(size, 1e-3);
              
              if (restart)
              {
                gsl_vector_free(xx);
                gsl_multimin_fminimizer_free(s);
                gsl_vector_free(step);
                //Rprintf("Likelihood not finite; trying again\n");
                restart = false;
                goto Vonvorn;
              }
              
              if (status == GSL_SUCCESS) {
                //Rprintf("Minimum found at:\n");
                //Rprintf("%5d %10.5f\n", int(iter), s->fval);
              }
              
              if (int(iter) == imax) {
                //Rprintf("Minimum not found:\n");
                //Rprintf("%5d %10.5f\n", int(iter), s->fval);
              }
            } while (status == GSL_CONTINUE && static_cast<int>(iter) < imax);
            
            if (s->fval < oldfit) {
              oldfit = s->fval;
              gsl_vector_memcpy(xsave, xx);
            }
            gsl_vector_free(xx);
            gsl_multimin_fminimizer_free(s);
            gsl_vector_free(step);
          }
          
          int iz = 0;
          for (int type = 0; type != 3; type++)
            for (int ip = 0; ip != ifree[type]; ip++)
              if (dCOMP(type, ip)) dAVW(t, type, ip) = gsl_vector_get(xsave,iz++);
              lambdas[t] = gsl_vector_get(xsave,icompg);
              lambdas[t + indi] = exp(gsl_vector_get(xsave,icompg + 1));
              gsl_vector_free(xsave);
        }
      });
    }
    
    /* ... the main thread also runs */
    for (int t = (AmntOfThreads-1)*NperThread; t < indi; t++) {
      double oldfit = GSL_POSINF;
      
      progress = 1.0*(t+1-(AmntOfThreads-1)*NperThread)/(indi-(AmntOfThreads-1)*NperThread);
      
      gsl_vector* xx = gsl_vector_alloc(icompg + 2);
      gsl_vector* xsave = gsl_vector_alloc(icompg + 2);
      
      bool restart = false;
      std::vector<std::vector<double>> icdaten; icdaten.clear();
      int ndat = 0;
      for (int j = 0; j != kerncat; j++) {
        icdaten.push_back(itcdaten[t * kerncat + j]);
        ndat += static_cast<int>(itcdaten[t * kerncat + j].size());
      }
      
      int n = icompg + 2;
      int m = ifreeg + 2;
      one_params pars = { n, m, restart, icdaten };
      
      double mean = 0.0, sx2 = 0.0;
      for (int j = 0; j != kerncat; j++) {
        int icdjs = icdaten[j].size();
        for (int ix = 0; ix != icdjs; ix++) {
          double xd = icdaten[j][ix];
          mean += xd;
          sx2 += gsl_pow_2(xd);
        }
      }
      mean /= ndat;
      sx2 /= ndat;
      sx2 -= gsl_pow_2(mean); sx2 = sqrt(sx2);
      
      for (int ix = 0; ix != 1; ix++) {
        Vonvorn:
        gsl_vector* xx = gsl_vector_alloc(icompg + 2);
        for (int i = 0; i != icompg; i++) gsl_vector_set(xx, i, oneuni(rst) - 0.5);
        gsl_vector_set(xx, icompg, mean * 0.7 + 0.01 * (oneuni(rst) - 0.5));
        gsl_vector_set(xx, icompg + 1, log(sx2 / 5.0));
        
        size_t iter = 0;
        int status;
        const gsl_multimin_fminimizer_type* T;
        gsl_multimin_fminimizer* s;
        
        gsl_multimin_function my_func;
        
        gsl_vector* step;
        step = gsl_vector_alloc(icompg + 2);
        for (int i = 0; i != icompg + 2; i++) gsl_vector_set(step, i, 0.1);
        
        my_func.n = n;
        my_func.f = objfun;
        my_func.params = &pars;
        
        T = gsl_multimin_fminimizer_nmsimplex2;
        s = gsl_multimin_fminimizer_alloc(T, n);
        
        gsl_multimin_fminimizer_set(s, &my_func, xx, step);
        
        do
        {
          iter++;
          status = gsl_multimin_fminimizer_iterate(s);
          if (status)
            break;
          size = gsl_multimin_fminimizer_size(s);
          status = gsl_multimin_test_size(size, 1e-3);
          
          if (restart)
          {
            gsl_vector_free(xx);
            gsl_multimin_fminimizer_free(s);
            gsl_vector_free(step);
            //Rprintf("Lielihood not finite; trying again\n");
            restart = false;
            goto Vonvorn;
          }
          
          if (status == GSL_SUCCESS) {
            //Rprintf("Minimum found at:\n");
            //Rprintf("%5d %10.5f\n", int(iter), s->fval);
          }
          
          if (int(iter) == imax) {
            //Rprintf("Minimum not found:\n");
            //Rprintf("%5d %10.5f\n", int(iter), s->fval);
          }
        } while (status == GSL_CONTINUE && static_cast<int>(iter) < imax);
        
        if (s->fval < oldfit) {
          oldfit = s->fval;
          gsl_vector_memcpy(xsave, xx);
        }
        gsl_vector_free(xx);
        gsl_multimin_fminimizer_free(s);
        gsl_vector_free(step);
      }
      
      int iz = 0;
      for (int type = 0; type != 3; type++)
        for (int ip = 0; ip != ifree[type]; ip++)
          if (dCOMP(type, ip)) dAVW(t, type, ip) = gsl_vector_get(xsave,iz++);
          lambdas[t] = gsl_vector_get(xsave,icompg);
          lambdas[t + indi] = exp(gsl_vector_get(xsave,icompg + 1));
          gsl_vector_free(xsave);
          
          //#pragma optimize("", off)
          if (PROG_BAR_FLAG) {
            int pos = static_cast<int>(ML_bar * progress);
            Rprintf("\r[");
            for (int i = 0; i < ML_bar; i++) {
              if (i < pos) {
                Rprintf("=");
              } else if (i == pos) {
                Rprintf(">");
              } else Rprintf(" ");
            }
            Rprintf("] %d%%", static_cast<int>(progress * 100.0));
          }
          //#pragma optimize("", on)
    }
    
    for (int ithread = 0; ithread < AmntOfThreads-1; ithread++) {
      threads[ithread].join();
    }
    Rprintf("\n");
    
  }

}
