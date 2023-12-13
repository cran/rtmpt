// authors: Christoph Klauer and Raphael Hartmann
#include "rts.h"
#include <atomic>

namespace ertmpt {
  
  std::atomic<int> atm_int(0);
  
  double  mlamb(int t, int pm, int ip, double *lambdas, double *lams)
  {
  	int ik = kernpar * (1 + pm) + ip;
  	if (comp[ik]) {
  		int iz = kern2free[ik] - ifree;
  		return exp(lambdas[t*ilamfree + iz] * lams[ifree + iz]);
  	}
  	else return(1.0e10);
  }
  
  void push(int ithread, int n_value_store, int n_all_parameters,double *factor, double *mu, double *lams, double *rhos, double *beta, double* lambdas,
  	double *restpars, double *slams, double *valuestore, double *parmon, double *parmonstore) {
  	int offset = ithread * n_value_store;
  	for (int i = 0; i != ifree * igroup; i++) valuestore[offset + i] = mu[i];
  	offset += ifree * igroup;
  	for (int i = 0; i != ifree + ilamfree; i++) valuestore[offset + i] = lams[i];
  	offset += ifree + ilamfree;
  	for (int i = 0; i != respno; i++) valuestore[offset + i] = slams[i];
  	offset += respno;
  	for (int i = 0; i != ilamfree * igroup; i++) valuestore[offset + i] = rhos[i];
  	offset += ilamfree * igroup;
  
  	for (int i = 0; i != indi * ifree; i++) valuestore[offset + i] = beta[i];
  	offset += indi * ifree;
  	for (int i = 0; i != indi * ilamfree; i++) valuestore[offset + i] = lambdas[i];
  	offset += indi * ilamfree;
  	for (int i = 0; i != restparsno; i++) valuestore[offset + i] = restpars[i];
  	offset += restparsno;
  	for (int i = 0; i != indi * respno; i++) valuestore[offset + i] = factor[i];
  	offset += indi * respno;
  
  	offset = ithread * 2 * n_all_parameters;
  	for (int i = 0; i != 2 * n_all_parameters; i++) parmonstore[offset + i] = parmon[i];
  	offset = 0;
  }
  
  void pop(int ithread, int n_value_store, int n_all_parameters,double *factor, double *mu, double *lams, double *rhos, double *beta, double* lambdas,
  	double *restpars, double *slams, double *valuestore, double *parmon, double *parmonstore) {
  
  	int offset = ithread * n_value_store;
  	for (int i = 0; i != ifree * igroup; i++) mu[i] = valuestore[offset + i];
  	offset += ifree * igroup;
  	for (int i = 0; i != ifree + ilamfree; i++) lams[i] = valuestore[offset + i];
  	offset += ifree + ilamfree;
  	for (int i = 0; i != respno; i++) slams[i] = valuestore[offset + i];
  	offset += respno;
  	for (int i = 0; i != ilamfree * igroup; i++) rhos[i] = valuestore[offset + i];
  	offset += ilamfree * igroup;
  
  	for (int i = 0; i != indi * ifree; i++) beta[i] = valuestore[offset + i];
  	offset += indi * ifree;
  	for (int i = 0; i != indi * ilamfree; i++) lambdas[i] = valuestore[offset + i];
  	offset += indi * ilamfree;
  	for (int i = 0; i != restparsno; i++) restpars[i] = valuestore[offset + i];
  	offset += restparsno;
  	for (int i = 0; i != indi * respno; i++) factor[i]= valuestore[offset + i] ;
  	offset += indi * respno;
  
  	offset = ithread * 2 * n_all_parameters;
  	for (int i = 0; i != 2 * n_all_parameters; i++) parmon[i] = parmonstore[offset + i];
  	offset = 0;
  }
  
  
  
  
  void make_rhos(int *nnodes, double *lambdas, double *lams, double *taus, double* rhos, gsl_rng *rst) {
  	//	NagError fail;	INIT_FAIL(fail);
  	double prior = pr_shape_exp_mu_beta;
  
  #define NNODES(I,J) nnodes[I*kernpar+J]
  #define RHOS(IG,IZ) rhos[(IG)*ilamfree+IZ]
  #define LAMBDAS(T,IZ) lambdas[T*ilamfree+IZ] //PM=0 negativ; PM=1 positiv
  
  	double *n = 0;	n = (double *)calloc(igroup * ilamfree, sizeof(double));
  	double *te = 0;	te = (double *)calloc(indi * ilamfree, sizeof(double));
  	double *p = 0;	p = (double *)calloc(igroup * ilamfree, sizeof(double));
  
  
  // Rprintf("hier 1\n");
  	int jj = 0;
  
  	for (int ip = 0; ip != kernpar; ip++) if ((comp[kernpar + ip]) || (comp[ip + 2 * kernpar])) {
      for (int t = 0; t != indi; t++) {
        int ig = t2group[t];
        for (int pm = 0; pm != 2; pm++) if (comp[ip + (1 + pm) * kernpar]) {
          int iz = kern2free[ip + (1 + pm) * kernpar] - ifree;
          n[ig * ilamfree + iz] += NNODES(t, ip) * 1.0;
        }
        for (int i = 0; i != NNODES(t, ip); i++) {
          for (int pm = 0; pm != 2; pm++) if (comp[ip + (1 + pm) * kernpar]) {
            int iz = kern2free[ip + (1 + pm) * kernpar] - ifree;
            te[t * ilamfree + iz] += taus[jj];
            jj++;
          }
        }
    	}
    }
  
    for (int iz = 0; iz != ilamfree; iz++) {
      int ip = free2kern[iz + ifree];
      int pm = (ip < 2 * kernpar) ? 0 : 1;
      ip -=  (1 + pm) * kernpar;
      for (int t = 0; t != indi; t++) {
        int ig = t2group[t];
        p[ig * ilamfree + iz] += mlamb(t, pm, ip, lambdas, lams) * te[t * ilamfree + iz];
      }
  
      for (int ig = 0; ig != igroup; ig++) {
        double x;
  			double a = n[ig * ilamfree + iz] + prior;
  			double b = p[ig * ilamfree + iz] + pr_rate_exp_mu_beta; //prior / 10.0;
        b = 1.0 / b;
        x = gsl_ran_gamma(rst, a, b);
        RHOS(ig, iz) = x;
      }
    }
  
  
  /*
  	for (int ip = 0; ip != kernpar; ip++) if ((comp[kernpar + ip]) || (comp[ip + 2 * kernpar])) {
  		//		int iz=index[ip];
  		for (int ig = 0; ig != igroup; ig++) { n[ig] = 0; p[ig] = 0.0; p[ig + igroup] = 0.0; }
  		for (int t = 0; t != indi; t++) {
  			int ig = t2group[t];
  			n[ig] += NNODES(t, ip)*1.0; te[ig] = 0.0; te[ig + igroup] = 0.0;
  			for (int i = 0; i != NNODES(t, ip); i++) {
  				for (int pm = 0; pm != 2; pm++) if (comp[ip + (1 + pm)*kernpar]) { te[pm*igroup + ig] += taus[jj];	jj++; }
  			}
  			for (int pm = 0; pm != 2; pm++) if (comp[ip + (1 + pm)*kernpar]) p[pm*igroup + ig] += mlamb(t, pm, ip, lambdas, lams)*te[pm*igroup + ig];
  		}
  		for (int ig = 0; ig != igroup; ig++)
  			for (int pm = 0; pm != 2; pm++) if (comp[ip + (1 + pm)*kernpar])
  			{
  				double x[1]; double a = n[ig] + prior; double b = p[pm*igroup + ig] + pr_rate_exp_mu_beta;//prior / 10.0;
  				b = 1.0 / b;
  				x[0] = gsl_ran_gamma(rst, a, b);
  				int iz = kern2free[(1 + pm)*kernpar + ip] - ifree;
  				RHOS(ig, iz) = x[0];
  			}
  	}
  */
  
    if (n) free(n);
  	if (p) free(p);
  	if (te) free(te);
  }
  
  void lambda_cond(double scale, double norm, double n, double alpha, double p, double *beta, double *sigi, double *lambdas, double *lams, int t, int iz, bool deriv, point &h) {
  #define BETA(T,I) beta[T*ifree+I]
  #define SIGI(I,J) sigi[(I)*(ifree+ilamfree)+J]
  	h.x = alpha; int ipar = ifree + iz;
  	alpha = alpha / scale;
  	if (!deriv) {
  		h.h = (alpha*lams[ipar])*n - p * exp(alpha*lams[ipar]) - norm;
  		int jj = 0;
  		for (int ix = 0; ix != ifree; ix++) { h.h -= (alpha)*SIGI(jj, ipar)*BETA(t, ix); jj++; }
  		for (int ix = 0; ix != ilamfree; ix++) if (jj != ipar) { h.h -= (alpha)*SIGI(ipar, jj)*LAMBDAS(t, ix); jj++; }
  		else { h.h -= 0.5*(alpha)*SIGI(ipar, jj)*(alpha); jj++; }
  	}
  	if (deriv) {
  		h.dh = (n - p * exp(alpha*lams[ipar]))*lams[ipar];
  		int jj = 0;
  		for (int ix = 0; ix != ifree; ix++) { h.dh -= SIGI(jj, ipar)*BETA(t, ix); jj++; }
  		for (int ix = 0; ix != ilamfree; ix++) if (jj != ipar) { h.dh -= SIGI(ipar, jj)*LAMBDAS(t, ix); jj++; }
  		else {
  			h.dh -= (alpha)*SIGI(ipar, jj); jj++;
  
  		}
  		h.dh /= scale;
  	}
  
  
  }
  
  
  
  void make_lambdas_new(int *nnodes, double *taus, double *beta, double *sigi, double *rhos, double *lambdas, double *lams, gsl_rng *rst) {
  
  	double* n = 0;	n = (double*)calloc(indi * ilamfree, sizeof(double));
  	double* p = 0;	p = (double*)calloc(indi * ilamfree, sizeof(double));
  
    int jj = 0;
  
    for (int ip = 0; ip != kernpar; ip++) if ((comp[kernpar + ip]) || (comp[ip + 2 * kernpar])) {
      for (int t = 0; t != indi; t++) {
        for (int pm = 0; pm != 2; pm++) if (comp[ip + (1 + pm) * kernpar]) {
          int iz = kern2free[(1 + pm) * kernpar + ip] - ifree;
          n[t * ilamfree + iz] += NNODES(t, ip);
        }
        for (int i = 0; i != NNODES(t, ip); i++) {
          for (int pm = 0; pm != 2; pm++) if (comp[ip + (1 + pm) * kernpar]) {
            int iz = kern2free[(1 + pm) * kernpar + ip] - ifree;
            p[t * ilamfree + iz] += taus[jj]; jj++;
          }
        }
      }
    }
    for (int iz = 0; iz != ilamfree; iz++) {
      for (int t = 0; t != indi; t++) {
        p[t * ilamfree + iz] *= RHOS(t2group[t], iz);
        // int ipar = iz + ifree;
        // double b = p[t * ilamfree + iz];
        double a = n[t * ilamfree + iz];
        double start = 0.0;
        //AGAIN:
        double step = 1, scale = (a > 0) ? sqrt(a) : 1.0; start = start * scale;
  			double totallow = -DBL_MAX;//GSL_NEGINF;
        double temp = ars(step, scale, totallow, a, p[t * ilamfree + iz], beta, sigi, lambdas, lams, t, iz, start, rst, lambda_cond);
        temp = temp / scale;
        LAMBDAS(t, iz) = temp;
      }
    }
  
  	if (n) free(n);
  	if (p) free(p);
  
  /*
  	for (int ip = 0; ip != kernpar; ip++) if ((comp[kernpar + ip]) || (comp[ip + 2 * kernpar])) {
  
  		for (int t = 0; t != indi; t++) {
  			double n = NNODES(t, ip)*1.0; double p[2] = { 0.0,0.0 };
  			for (int i = 0; i != NNODES(t, ip); i++) {
  				for (int pm = 0; pm != 2; pm++) if (comp[ip + (1 + pm)*kernpar]) { p[pm] += taus[jj]; jj++; }
  			}
  			for (int pm = 0; pm != 2; pm++) if (comp[ip + (1 + pm)*kernpar]) {
  				int iz = kern2free[(1 + pm)*kernpar + ip] - ifree;
  				p[pm] *= RHOS(t2group[t], iz);
  
  				int ipar = iz + ifree; double sm = SIGI(ipar, ipar);
  				double b = p[pm];
  				double a = n;
  				double start = 0.0;
  
  				double step = 1, scale = (n > 0) ? sqrt(n) : 1.0; start = start * scale; double totallow = -DBL_MAX;
  				double temp = ars(step, scale, totallow, n, p[pm], beta, sigi, lambdas, lams, t, iz, start, rst, lambda_cond);
  				temp = temp / scale;
  				LAMBDAS(t, iz) = temp;
  			}
  		}
  	}
  */
  }
  
  void lam2(double scale, double norm, double n, double alpha, double p, double *beta, double *sigi, double *lambdas, double *lams, int t, int iz, bool deriv, point &h) {
  // übergebe tau-Summen mal rho in Beta pm*indi + t
  
  	h.x = alpha; int ipar = iz;
  	alpha = 1 + alpha / scale;
  
  	if (!deriv) {
  	    double lax = LAMBDAS(0, iz);
  	    double temp = alpha * lax + log(beta[ipar]);
  	    for (int t = 1; t != indi; t++) {
  	        double lax = LAMBDAS(t, iz);
  	        temp = logsum(temp, alpha*lax + log(beta[t*ilamfree + ipar]));
  	    }
  
  	    if (alpha*n - norm > 0)
  	        h.h = elogdiff(log(alpha*n - norm), logsum(log(0.5*PRIOR*(alpha - 1.0)*(alpha - 1.0)), temp));
  	    else
  	        h.h = -rexp(logsum(log(norm - alpha * n), logsum(log(0.5*PRIOR*(alpha - 1.0)*(alpha - 1.0)), temp)));
  
  	}
  	if (deriv) {
  	    h.dh = log(fabs(n - PRIOR*(alpha - 1.0))); double temp_plus, temp_minus; bool fplus = true, fminus = true;
  	    for (int t = 0; t != indi; t++) {
  
  	        double lax = LAMBDAS(t, iz);
  	        if (lax > 0) {
  	            double temp = log(lax) + alpha * lax + log(beta[t*ilamfree + ipar]);
  	            if (fplus) { temp_plus = temp; fplus = false; }
  	            else temp_plus = logsum(temp_plus, temp);
  	        }
  	        else {
  	            double temp = log(-lax) + alpha * lax + log(beta[t*ilamfree + ipar]);
  	            if (fminus) { temp_minus = temp; fminus = false; }
  	            else temp_minus = logsum(temp_minus, temp);
  	        }
  	    }
  	    if (n - PRIOR*(alpha - 1.0) > 0)
  	        if (fminus) h.dh = (elogdiff(h.dh, temp_plus)); else if (fplus) h.dh = rexp(logsum(h.dh, temp_minus)); else { h.dh = logsum(h.dh, temp_minus); h.dh = elogdiff(h.dh, temp_plus); }
  	    else
  	        if (fminus) h.dh = -rexp(logsum(h.dh, temp_plus)); else if (fplus) h.dh = -elogdiff(h.dh, temp_minus); else { h.dh = logsum(h.dh, temp_plus); h.dh = -elogdiff(h.dh, temp_minus); }
  	    h.dh /= scale;
  	}
  }
  
  
  void make_lamb2(int *nnodes, double *taus, double *beta, double *sigi, double *rhos, double *lambdas, double *lams, gsl_rng *rst) {
  
  	double *b = 0;	b = (double *)calloc(ilamfree*indi , sizeof(double));
  	double *m = 0;	m = (double *)calloc(ilamfree , sizeof(double));
  	double *n = 0;	n = (double *)calloc(ilamfree , sizeof(double));
  	int jj = 0;
  
  	for (int ip = 0; ip != kernpar; ip++) if ((comp[kernpar + ip]) || (comp[ip + 2 * kernpar])) {
  		for (int t = 0; t != indi; t++) {
  			for (int pm = 0; pm != 2; pm++) if (comp[ip + (1 + pm) * kernpar]) {
  				int iz = kern2free[(1 + pm) * kernpar + ip] - ifree;
  				double lax = LAMBDAS(t, iz);
  				m[iz] += NNODES(t, ip);
  				n[iz] += NNODES(t, ip) * lax;
  			}
  			for (int i = 0; i != NNODES(t, ip); i++) {
  				for (int pm = 0; pm != 2; pm++) if (comp[ip + (1 + pm) * kernpar])
  				{
  					int iz = kern2free[(1 + pm) * kernpar + ip] - ifree;
  					b[t * ilamfree + iz] += taus[jj]; jj++;
  				}
  			}
  		}
  	}
  	for (int iz = 0; iz != ilamfree; iz++) {
  		for (int t = 0; t!= indi; t++) b[t * ilamfree + iz] *= (RHOS((t2group[t]), iz));
  		double start = 0.0; int xt = 0;
  		double step = 1, totallow;//=(0.1-1.0)*sqrt(abs(n[pm]));
  		totallow = -DBL_MAX;
  		double scale = sqrt(m[iz] / indi);
  		double temp = ars(step, scale, totallow, n[iz], n[iz], b, sigi, lambdas, lams, xt, iz, start, rst, lam2); //std::cout << temp<< std::endl;
  		lams[ifree + iz] = 1 + temp / scale;
  	}
  
  
  /*
  	for (int ip = 0; ip != kernpar; ip++) if ((comp[kernpar + ip]) || (comp[ip + 2 * kernpar])) {
  		double p[2] = { 0.0,0.0 }, n[2] = { 0.0,0.0 }, m[2] = { 0.0,0.0 };
  		for (int t = 0; t != indi; t++) {
  			for (int pm = 0; pm != 2; pm++) if (comp[ip + (1 + pm)*kernpar]) {
  				int iz = kern2free[(1 + pm)*kernpar + ip] - ifree;
  				double lax = LAMBDAS(t, iz);
  				m[pm] += NNODES(t, ip);
  				n[pm] += NNODES(t, ip)*lax;
  				b[t*(ilamfree)+iz] = 0.0;
  			}
  			for (int i = 0; i != NNODES(t, ip); i++) {
  				for (int pm = 0; pm != 2; pm++) if (comp[ip + (1 + pm)*kernpar])
  				{
  					int iz = kern2free[(1 + pm)*kernpar + ip] - ifree; b[t*(ilamfree)+iz] += taus[jj]; jj++;
  				}
  			}
  			for (int pm = 0; pm != 2; pm++)  if (comp[ip + (1 + pm)*kernpar])
  			{
  				int iz = kern2free[(1 + pm)*kernpar + ip] - ifree;
  				b[t*(ilamfree)+iz] *= (RHOS((t2group[t]), iz));
  			}
  		}
  		for (int pm = 0; pm != 2; pm++)  if (comp[ip + (1 + pm)*kernpar]) {
  
  			double start = 0.0; int xt = 0;
  			double step = 1, totallow;//=(0.1-1.0)*sqrt(abs(n[pm]));
  			totallow = -DBL_MAX;
  			double scale = sqrt(m[pm] / indi);
  			int iz = kern2free[(1 + pm)*kernpar + ip] - ifree;
  			double temp = ars(step, scale, totallow, n[pm], p[pm], b, sigi, lambdas, lams, xt, iz, start, rst, lam2); //std::cout << temp<< std::endl;
  			if (DEBUG) if (fabs(temp) > 1000) Rprintf("ars %g\n", 1 + temp / sqrt(fabs(m[pm])));//std::cout << "ars " << 1 + temp / sqrt(fabs(m[pm])) << std::endl;
  
  			lams[ifree + iz] = 1 + temp / scale;
  		}
  	}
  */
  
  	if (b) free(b);
  	if (m) free(m);
  	if (n) free(n);
  }
  
  void make_betas_new(double *mu, double *lams, double *beta, double *sigi, int *nnodes, double *z, double *lambdas, gsl_rng *rst) {
  	//	NagError fail;	INIT_FAIL(fail);
  
  	double* w = 0; w = (double*)malloc(ifree * sizeof(double));
  	double* hba = 0; hba = (double*)malloc(ifree * sizeof(double));
  	double* fig = 0; fig = (double*)malloc(indi * ifree * sizeof(double));
  	double* xfig = 0;	xfig = (double*)malloc(ifree * ifree * sizeof(double));
  	double* ba = 0;	ba = (double*)malloc(indi * ifree * sizeof(double));
  #define FIG(T,I) fig[T*ifree+I]
  #define XFIG(I,J) xfig[I*ifree+J]
  #define BA(T,I) ba[T*ifree+I]
  
  
  	for (int t = 0; t != indi; t++) for (int iz = 0; iz != ifree; iz++) { BA(t, iz) = 0.0; FIG(t, iz) = 0.0; }
  
  
  
  	int jj = -1;
  	for (int ip = 0; ip != kernpar; ip++) if (comp[ip]) {
  		int iz = kern2free[ip];
  		for (int t = 0; t != indi; t++) {
  			double figtiz = 0, batiz = 0;
  			double be = equation(t, ip, mu, lams, beta) - BETA(t, iz) * lams[iz];
  			figtiz += NNODES(t, ip);
  			for (int j = 0; j != NNODES(t, ip); j++) {
  				batiz += (z[++jj] - be);
  			}
  			FIG(t, iz) += figtiz * gsl_pow_2(lams[iz]);
  			BA(t, iz) += batiz * lams[iz];
  		}
  	}
  	for (int t = 0; t != indi; t++) {
  		for (int iz = 0; iz != ifree; iz++) {
  			w[iz] = BA(t, iz);
  			for (int jz = ifree; jz != ifree + ilamfree; jz++) w[iz] -= SIGI(iz, jz) * lambdas[t * ilamfree + jz - ifree];
  		}
  		for (int iz = 0; iz != ifree; iz++) {
  			for (int jz = 0; jz != ifree; jz++) if (iz != jz) XFIG(iz, jz) = SIGI(iz, jz); else XFIG(iz, iz) = FIG(t, iz) + SIGI(iz, iz);
  		}
  		bayesreg(ifree, w, xfig, hba, rst);
  		for (int iz = 0; iz != ifree; iz++)  BETA(t, iz) = hba[iz];
  	}
  
  
  	if (w) free(w);
  	if (hba) free(hba);
  	if (fig) free(fig);
  	if (xfig) free(xfig);
  	if (ba) free(ba);
  }
  
  
  
  
  
  
  void sample_sig(double *beta, double *lambdas, double *sig, double *sigi, gsl_rng *rst) {
  #define XY(I,J) xy[(I)*(ifree+ilamfree) + J]
  
  	double *xy = 0;	xy = (double *)malloc((indi + ifree + ilamfree + 1 + pr_df_add_inv_wish)*(ifree + ilamfree) * sizeof(double));
  	for (int i = 0; i != indi; i++) {
  		for (int iz = 0; iz != ifree; iz++) XY(i, iz) = BETA(i, iz);
  		for (int j = 0; j != ilamfree; j++) XY(i, ifree + j) = LAMBDAS(i, j);
  	}
  	invwis(indi, (ifree + ilamfree), xy, sig, sigi, pr_sf_scale_matrix_SIG, rst);
  	if (xy) free(xy);
  }
  
  #define TREE_AND_NODE2PAR(ITREE,R) tree_and_node2par[ITREE*nodemax+R]
  #define DRIN(J,K,X) drin[J*zweig*nodemax+K*nodemax + X]
  #define NDRIN(J,K) ndrin[J*zweig+K]
  #define AR(I,J,R) ar[I*zweig*nodemax + J*nodemax + R]
  
  
  
  void make_taus_one_trial_new_new(trial one, int itrial, int ipath, double *rhos, double *lambdas, double *lams, int* ntau_position, double *taus, double *rest, double *restpars, double* slams, gsl_rng *rst) {
  #define NTAU_POSITION(X,J,PM) ntau_position[2*X*nodemax+2*J+PM]
  
  
  	int t = one.person; int itree = one.tree; int j = one.category; double rt = one.rt / 1000.0; int resp = cat2resp[j];
  
  
  	// analyze path to find lambda_smallest
  
  
  
  	double la_min = -1.0; int pfadlength = NDRIN(j, ipath); int r_min = -1;
  	for (int x = 0; x != pfadlength; x++) {
  		int r = DRIN(j, ipath, x); double la_temp;
  		int ip = TREE_AND_NODE2PAR(itree, r);
  		int pm = (AR(j, ipath, r) > 0) ? 1 : 0;
  		if (comp[ip + (1 + pm)*kernpar]) {
  			int iz = kern2free[(1 + pm)*kernpar + ip] - ifree;
  			la_temp = mlamb(t, pm, ip, lambdas, lams)*RHOS(t2group[t], iz);
  			if ((la_min == -1) || (la_temp < la_min)) { la_min = la_temp; r_min = r; }
  		}
  	}
  	if (la_min == -1.0) pfadlength = 0;
  	if (pfadlength == 0) rest[itrial] = rt;
  	else {
  		double rest_old, upper;
  		double rmu = restpars[t2group[t] * respno + resp] + malpha(t, resp, restpars, slams) + la_min * restpars[t + sigalphaoff]; double rsig = sqrt(restpars[t + sigalphaoff]); double lower = -rmu; lower /= rsig;
  		double bound; if (pfadlength > 1) { bound = rt - rmu;   bound /= rsig; }
  
  		int ix = 0;
  	NEW: rest_old = rt;
  		if (pfadlength > 1) {
  			upper = rt;
  			for (int x = 0; x != pfadlength; x++) if (DRIN(j, ipath, x) != r_min) {
  				int r = DRIN(j, ipath, x); int ip = TREE_AND_NODE2PAR(itree, r);
  				int pm = (AR(j, ipath, r) > 0) ? 1 : 0;
  				if (comp[ip + (1 + pm)*kernpar]) {
  					int tau_pos = NTAU_POSITION(itrial, r, pm);
  					int iz = kern2free[(1 + pm)*kernpar + ip] - ifree;
  					double lambda = mlamb(t, pm, ip, lambdas, lams)*RHOS(t2group[t], iz);
  					if (lambda != la_min) taus[tau_pos] = truncexp(lambda - la_min, upper, rst);
  					else taus[tau_pos] = oneuni(rst)*upper;
  					rest_old -= taus[tau_pos];
  					if (rest_old < 0) goto NEW;
  				}
  			}
  		}
  
  		// truncated normal
  
  
  		upper = rest_old - rmu;  upper /= rsig;
  		if (pfadlength > 1) {
  			double cut = gsl_cdf_ugaussian_P(lower); cut = (gsl_cdf_ugaussian_P(upper) - cut) / (gsl_cdf_ugaussian_P(bound) - cut);
  			double u = oneuni(rst);
  			if (u > cut) { if (ix > 100) Rprintf("%13d", ix++); /*std::cout << std::setw(13) << ix++;*/ goto NEW; }
  		}
  		double zz = double_truncnorm(lower, upper, rst);
  		rest[itrial] = zz * rsig + rmu;
  		rest_old -= rest[itrial];
  
  
  
  		int pm = (AR(j, ipath, r_min) > 0) ? 1 : 0; int tau_pos = NTAU_POSITION(itrial, r_min, pm);  taus[tau_pos] = rest_old;
  	}
  
  	for (int r = 0; r != nodes_per_tree[itree]; r++) {
  		int ip = TREE_AND_NODE2PAR(itree, r);
  		if (AR(j, ipath, r) != 0) {
  			if ((AR(j, ipath, r) > 0) && (comp[ip + kernpar])) { int tau_pos = NTAU_POSITION(itrial, r, 0); int iz = kern2free[kernpar + ip] - ifree; double lambda = mlamb(t, 0, ip, lambdas, lams)*RHOS(t2group[t], iz);   taus[tau_pos] = oneexp(lambda, rst); }
  			if ((AR(j, ipath, r) < 0) && (comp[ip + 2 * kernpar])) { int tau_pos = NTAU_POSITION(itrial, r, 1); int iz = kern2free[2 * kernpar + ip] - ifree; double lambda = mlamb(t, 1, ip, lambdas, lams)*RHOS(t2group[t], iz); taus[tau_pos] = oneexp(lambda, rst); }
  		}
  		else for (int pm = 0; pm != 2; pm++)
  			if (comp[ip + (1 + pm)*kernpar]) {
  				int tau_pos = NTAU_POSITION(itrial, r, pm); int iz = kern2free[(1 + pm)*kernpar + ip] - ifree; double lambda = mlamb(t, pm, ip, lambdas, lams)*RHOS(t2group[t], iz);  taus[tau_pos] = oneexp(lambda, rst);
  			}
  	}
  }
  
  #define PFAD_INDEX(J,K) pfad_index[J*zweig+K]
  
  void make_tij_for_one_trial_new(trial one, double *rhos, double *lambdas, double *xlams, double *restpars, double* slams, double *pij) {
  
  	int t = one.person; /*int itree = one.tree;*/ int j = one.category; double rt = one.rt / 1000.0; int resp = cat2resp[j];
  
  	double rmu = restpars[t2group[t] * respno + resp] + malpha(t, resp, restpars, slams); double rsig = sqrt(restpars[t + sigalphaoff]);
  
  	for (int k = 0; k != branch[j]; k++) {
  		int pfadlength = NDRIN(j, k);
  		double *lams = 0; lams = (double *)malloc(pfadlength * sizeof(double));
  		int complength = 0;
  /*		if (PFAD_INDEX(j, k) == -1) {
  			for (int ir = 0; ir != pfadlength; ir++) {
  				int r = DRIN(j, k, ir); int ip = TREE_AND_NODE2PAR(itree, r); int pm = (AR(j, k, r) > 0) ? 1 : 0;
  				if (comp[ip + (1 + pm)*kernpar]) {
  					int iz = kern2free[(1 + pm)*kernpar + ip] - ifree; lams[complength] = mlamb(t, pm, ip, lambdas, xlams)*RHOS(t2group[t], iz);
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
  			int iz = kern2free[(1 + pm)*kernpar + ip] - ifree; lams[ir] = mlamb(t, pm, ip, lambdas, xlams)*RHOS(t2group[t], iz);
  		}
  		// }
  		if (complength == 0) { pij[k] = (-0.5*gsl_pow_2((rt - rmu) / rsig)) / sqr2pi / rsig; }
  		// if ((complength == 1) && (PFAD_INDEX(j, k) == -1)) {
  		// 	double lam = lams[0]; pij[k] = (logexgaussian(lam, rmu, rsig, rt));
  		// 	if (DEBUG) {if (pij[k] < 0.0) Rprintf("pij[k] < 0; pf = 1 %g\n", pij[k]);} //std::cout << "pij[k] < 0; pf =1" << pij[k] << std::endl;
  		// }
  		if (complength == 1 /*&& (PFAD_INDEX(j, k) > -1)*/) {
  			// int ipfad = PFAD_INDEX(j, k);
  			// pfadinfo akt_pfad = path_info[ipfad];
  			if (complength != akt_pfad.a) Rprintf("complength != a");
  			double lam = lams[0];
  			double hplus, hminus;
  			loggammagaussian(akt_pfad.r[0] - 1, lam, rmu, rsig, rt, hplus, hminus);
  			double temp = logdiff(hplus, hminus) + log(lam) * akt_pfad.r[0];
  			if (temp == GSL_NEGINF) {
  				pij[k] = temp;
  			}
  			else pij[k] = temp;
  		}
  		/*if ((complength >= 2) && (PFAD_INDEX(j, k) == -1)) {
  			//	    double *lams=0; if ( !(lams = NAG_ALLOC(pfadlength*sizeof(double)) )){printf("Allocation failure\n");exit_status = -1;}
  			double *loglams = 0; loglams = (double *)malloc(complength * sizeof(double));
  
  			double hypo_plus, hypo_minus; bool f_plus = true, f_minus = true;
  			for (int ir = 0; ir != complength; ir++) loglams[ir] = log(lams[ir]);
  
  			for (int ir = 0; ir != complength; ir++) {
  				int sign = 1; double factor = 0.0;
  				for (int ij = 0; ij != complength; ij++) if (ij != ir) { if (lams[ij] - lams[ir] < 0) sign = -sign; factor = factor + loglams[ij] - log(fabs(lams[ij] - lams[ir])); }
  				double temp = factor + logexgaussian(lams[ir], rmu, rsig, rt);
  				if (sign > 0) { if (f_plus) { hypo_plus = temp; f_plus = false; } else hypo_plus = logsum(hypo_plus, temp); }
  				if (sign < 0) { if (f_minus) { hypo_minus = temp; f_minus = false; } else hypo_minus = logsum(hypo_minus, temp); }
  			}
  
  			double dif = hypo_minus - hypo_plus;
  			if (dif < 0) 	pij[k] = hypo_plus + gsl_log1p(-exp(hypo_minus - hypo_plus)); else pij[k] = GSL_NEGINF;
  
  			if (dif>0) {
  				if (DEBUG) Rprintf("pij[k] < 0; pf = 2 %g %g %g  %g\n", pij[k], hypo_minus, hypo_plus, rsig);
  				//std::cout << "pij[k] < 0; pf =2" << pij[k] << " " << hypo_minus << hypo_plus
  				//	<< " " << " " << rsig << std::endl;
  				if (DEBUG) {for (int ir = 0; ir != complength; ir++) Rprintf("%15g", lams[ir]); Rprintf("\n");} //std::cout << std::setw(15) << lams[ir]; std::cout << std::endl;
  				pij[k] = -sqrt(DBL_MAX);
  			}
  
  			if (loglams) free(loglams);
  
  		}*/
  		if ((complength >= 2) /*&& (PFAD_INDEX(j, k) > -1)*/) {
  			double *loglams = 0; loglams = (double *)malloc(complength * sizeof(double));
  			for (int ir = 0; ir != complength; ir++) loglams[ir] = log(lams[ir]);
  			// int ipfad = PFAD_INDEX(j, k);
  			// pfadinfo akt_pfad = path_info[ipfad];
  			double temp = logf_tij(akt_pfad.a, akt_pfad.r, lams,loglams, rmu, rsig, rt);
  			if (temp == GSL_NEGINF) {
  				pij[k] = GSL_NEGINF;
  			}
  			else pij[k] = temp;
  			if (loglams) free(loglams);
  		}
  
  		if (lams) free(lams);
  	}
  
  }
  
  #define FACTOR(T,R) factor[T*respno+R]
  #define NPPR(T,R) nppr[(T)*respno + R]
  
  void gibbs_full_cycle(std::vector <trial> daten, double *factor, double *mu, double *lams, double *beta, double *sig, double *rhos, double *lambdas,
  	int ntau, int *ntau_position, double *taus, int nz, int *nz_position, int *nnodes, double *restpars, double *slams, bool xflag, gsl_rng *rst) {
  	//	NagError fail;	INIT_FAIL(fail);
  	double *x_for_all = 0; x_for_all = (double *)malloc(indi*kernpar * sizeof(double));
  	make_parameters_for_all(mu, lams, beta, x_for_all);
  
  	double *pij = 0; pij = (double *)malloc(zweig * sizeof(double));
  	double *z = 0; z = (double *)malloc(nz * sizeof(double));
  	double *rest = 0; rest = (double *)malloc(daten.size() * sizeof(double));
  	double *taui = 0; taui = (double *)malloc(respno*respno * sizeof(double));
  
  	int trialno = static_cast<int>(daten.size());
  
  	for (int x = 0; x != trialno; x++) {
  		double p;
  		trial one = daten[x];
  		make_tij_for_one_trial_new(one, rhos, lambdas, lams, restpars, slams, pij);
  		make_pij_for_one_trial(one, x_for_all, pij,p);
  		int ipath = make_path_for_one_trial(branch[one.category], pij,p, rst);
  		make_zs_one_trial(one, x, ipath, mu, lams, beta, nz_position, z, rst);
  		make_taus_one_trial_new_new(one, x, ipath, rhos, lambdas, lams, ntau_position, taus, rest, restpars, slams, rst);
  	}
  
  	make_rhos(nnodes, lambdas, lams, taus, rhos, rst);
  	make_rtau(restpars, taui, slams, rst);
  
  	if (xflag)
  		for (int t = 0; t != indi; t++)
  			for (int r = 0; r != respno; r++) {
  				double mu = malpha(t, r, restpars, slams) + restpars[t2group[t] * respno + r];
  				double rsig = sqrt(restpars[sigalphaoff + t]);
  				FACTOR(t, r) = lnnorm(mu / rsig)*NPPR(t, r);
  			}
  	make_rmu(daten, factor, rest, restpars, slams, rst);
  	make_slams(daten, factor, rest, restpars, slams, rst);
  	make_ralpha(daten, factor, rest, restpars, slams, taui, rst);
  	make_rsigalpha(daten, factor, rest, restpars, slams, xflag, rst);
  	make_rsig(daten, rest, restpars, rst);
  
  	make_mu(mu, lams, beta, nnodes, z, rst);
  	int npar = ifree + ilamfree;
  	double *sigi = 0; sigi = (double *)malloc(npar*npar * sizeof(double));
  	sample_sig(beta, lambdas, sig, sigi, rst);
  
  	make_lams(mu, lams, beta, nnodes, z, rst);
  	make_betas_new(mu, lams, beta, sigi, nnodes, z, lambdas, rst);
  	make_lamb2(nnodes, taus, beta, sigi, rhos, lambdas, lams, rst);
  	make_lambdas_new(nnodes, taus, beta, sigi, rhos, lambdas, lams, rst);
  
  	if (x_for_all) free(x_for_all);
  	if (pij) free(pij);
  	if (z) free(z);
  	if (sigi) free(sigi);
  	if (rest) free(rest);
  	if (taui) free(taui);
  }
  
  void on_screen3(int n_all_parameters, double *xwbr, double *parmon, double *beta, double rmax, int irun) {
  #define SIG(I,J) sig[I*(ifree+ilamfree)+J]
  #define XWBR(T,I) xwbr[(T-1)*n_all_parameters+I]
  #define PARMON(I,J) parmon[(I-1)*n_all_parameters+J]
  
  	double *sig = 0;	sig = (double *)malloc((ifree + ilamfree)*(ifree + ilamfree) * sizeof(double));;
  	int jz;
  	Rprintf("THETAS\nmean:"); //std::cout << "MUS" << std::endl;
  	jz = -1;
  	for (int ig = 0; ig != igroup; ig++) {
  		for (int jp = 0; jp != kernpar; jp++) if (comp[jp]) { jz = kern2free[jp] + ig * ifree; Rprintf("%15g", gsl_cdf_ugaussian_P(PARMON(1, jz))); /*std::cout << std::setw(15) << gsl_cdf_ugaussian_P(PARMON(1, jz));*/ }
  		else Rprintf("%15g", gsl_cdf_ugaussian_P(consts[jp])); //std::cout << std::setw(15) << gsl_cdf_ugaussian_P(consts[jp]);
  		Rprintf("\n");
  	}
  	Rprintf("Rhat:"); //std::cout << "R-statistic MUS" << std::endl;
  	jz = -1;
  	for (int ig = 0; ig != igroup; ig++) {
  		for (int jp = 0; jp != kernpar; jp++) if (comp[jp]) { jz = kern2free[jp] + ig * ifree; Rprintf("%15g", XWBR(3, jz)); /*std::cout << std::setw(15) << XWBR(3, jz);*/ }
  		else Rprintf("%15g", 0.0); //std::cout << std::setw(15) << 0.0;
  		Rprintf("\n");
  	}
  	Rprintf("--------\n");
  
  	Rprintf("LAMBDAS_MINUS\nmean:"); //std::cout << "RHOS MINUS" << std::endl;
  	jz = ifree * igroup;
  	for (int ig = 0; ig != igroup; ig++) {
  		for (int jp = 0; jp != kernpar; jp++) if (comp[jp + kernpar]) { jz = ifree * igroup + kern2free[kernpar + jp] - ifree + ig * ilamfree; Rprintf("%15g", PARMON(1, jz)); /*std::cout << std::setw(15) << PARMON(1, jz);*/ jz++; }
  		else Rprintf("%15g", 0.0); //std::cout << std::setw(15) << 0.0;
  		Rprintf("\n");
  	}
  	Rprintf("Rhat:"); //std::cout << "R-statistic RHOS MINUS" << std::endl;
  	jz = ifree * igroup;
  	for (int ig = 0; ig != igroup; ig++) {
  		for (int jp = 0; jp != kernpar; jp++) if (comp[kernpar + jp]) { jz = ifree * igroup + kern2free[kernpar + jp] - ifree + ig * ilamfree; Rprintf("%15g", XWBR(3, jz)); /*std::cout << std::setw(15) << XWBR(3, jz);*/ jz++; }
  		else Rprintf("%15g", 0.0); //std::cout << std::setw(15) << 0.0;
  		Rprintf("\n");
  	}
  	Rprintf("--------\n");
  
  	Rprintf("LAMBDAS_PLUS\nmean:"); //std::cout << "RHOS PLUS" << std::endl;
  	int js = jz;
  	for (int ig = 0; ig != igroup; ig++) {
  		for (int jp = 0; jp != kernpar; jp++) if (comp[jp + 2 * kernpar]) { jz = ifree * igroup + kern2free[2*kernpar + jp] -ifree + ig * ilamfree; Rprintf("%15g", PARMON(1, jz)); /*std::cout << std::setw(15) << PARMON(1, jz);*/ jz++; }
  		else Rprintf("%15g", 0.0); //std::cout << std::setw(15) << 0.0;
  		Rprintf("\n");
  	}
  	Rprintf("Rhat:"); //std::cout << "R-statistic RHOS PLUS" << std::endl;
  	jz = js;
  	for (int ig = 0; ig != igroup; ig++) {
  		for (int jp = 0; jp != kernpar; jp++) if (comp[jp + 2 * kernpar]) { jz = ifree * igroup + kern2free[2*kernpar + jp] -ifree + ig * ilamfree; Rprintf("%15g", XWBR(3, jz)); /*std::cout << std::setw(15) << XWBR(3, jz);*/ jz++; }
  		else Rprintf("%15g", 0.0); //std::cout << std::setw(15) << 0.0;
  		Rprintf("\n");
  	}
  	Rprintf("--------\n");
  
  
  	Rprintf("MU_GAMMAS, OMEGA^2\nmean:"); //std::cout << "RESTPARS" << std::endl;
  	for (int ig = 0; ig != igroup * respno + 1; ig++) 	Rprintf("%15g", PARMON(1, n_all_parameters - restparsno + ig)); //std::cout << std::setw(15) << PARMON(1, n_all_parameters - restparsno + ig);
  	Rprintf("\n"); //std::cout << std::endl;
  	Rprintf("Rhat:"); //std::cout << "R-statistic RESTPARS" << std::endl;
  	for (int ig = 0; ig != igroup * respno + 1; ig++) 	Rprintf("%15g", XWBR(3, n_all_parameters - restparsno + ig)); //std::cout << std::setw(15) << XWBR(3, n_all_parameters - restparsno + ig);
  	Rprintf("\n"); //std::cout << std::endl;
  	Rprintf("--------\n");
  
  	if ((rmax < RMAX) && !BURNIN_flag) RMAX_reached += 1;
  	else RMAX_reached = 0;
  	double pct_temp = (RMAX_reached>1) ? (100.0*ireps*(RMAX_reached-1)/(1.0*(THIN*SAMPLE_SIZE/NOTHREADS))) : 0.0;
  
  	Rprintf("max(Rhats): %g\n", rmax); //std::cout << "rmax " << rmax << std::endl;
  	if (!BURNIN_flag) Rprintf("Iterations: %d [sampling: %g%%]\n", (irun + 1)*ireps, pct_temp);
  	else Rprintf("Burnin: %d\n", BURNIN);
  	//std::cout << "Iterationen " << (irun + 1)*ireps << std::endl << std::endl;
  
  	if (RMAX_reached == 0 && !BURNIN_flag) Rprintf("Sampling starts when max(Rhats)<%g\n", RMAX);
  	if (RMAX_reached == 1) Rprintf("Sampling starts now.\n");
  
  	Rprintf("_____");
  	if (kernpar > igroup*respno) {
  		for (int i = 0; i < kernpar; i++) Rprintf("_______________");
  	} else {
  		for (int i = 0; i < (igroup*respno); i++) Rprintf("_______________");
  		Rprintf("_______________");
  	}
  	Rprintf("\n\n");
  	BURNIN_flag = false;
  
  	if (sig) free(sig);
  
  	R_CheckUserInterrupt();
  } // end on_screen3
  
  void belege_bridge(int ithread, int ix, int n_bridge_store, double *bridge_sample, double *mu, double *lams,
  	double *rhos, double *beta, double* lambdas, double *sig,
  	double *restpars, double *slams, double llik) {
  	int offset = (ithread * IREP + ix)*n_bridge_store;
  	for (int i = 0; i != ifree * igroup; i++) bridge_sample[offset + i] = mu[i];
  	offset += ifree * igroup;
  
  	for (int i = 0; i != ilamfree * igroup; i++) bridge_sample[offset + i] = rhos[i];
  	offset += ilamfree * igroup;
  
  	for (int i = 0; i != ifree + ilamfree; i++) bridge_sample[offset + i] = lams[i];
  	offset += ifree + ilamfree;
  
  	for (int t = 0; t != indi; t++) {
  		for (int j = 0; j != ifree; j++) bridge_sample[offset++] = BETA(t, j);
  		for (int j = 0; j != ilamfree; j++) bridge_sample[offset++] = LAMBDAS(t, j);
  	}
  	int n = ifree + ilamfree;
  	gsl_matrix *cx = gsl_matrix_alloc(n, n);
  
  
  	for (int j = 0; j != n; j++)
  		for (int i = j; i != n; i++) {
  			gsl_matrix_set(cx, i, j, SIG(i,j));
  			if (i != j) gsl_matrix_set(cx, j, i, SIG(i,j));
  		}
  	gsl_linalg_cholesky_decomp(cx);
  	for (int iz = 0; iz != ifree + ilamfree; iz++)
  		for (int jz = 0; jz <= iz; jz++) if (iz != jz) bridge_sample[offset++] = gsl_matrix_get(cx, iz, jz); else bridge_sample[offset++] = log(gsl_matrix_get(cx, iz, jz));
  	gsl_matrix_free(cx);
  
  	for (int i = 0; i != respno; i++) bridge_sample[offset + i] = slams[i];
  	offset += respno;
  
  
  	double *tau = (double *)malloc(respno*respno * sizeof(double));
  
  	double *temp = (double *)malloc(restparsno * sizeof(double));
  	for (int iz = 0; iz != restparsno; iz++) temp[iz] = restpars[iz];
  
  	n = respno;
  	gsl_matrix *cy = gsl_matrix_alloc(n, n);
  	int jj = 0;
  	for (int i = 0; i != n; i++)
  		for (int j = i; j !=n;j++) {
  			gsl_matrix_set(cy, i, j, restpars[igroup*respno+1+ jj]);
  			if (i != j) gsl_matrix_set(cy, j, i, restpars[igroup*respno + 1 + jj]);
  			jj++;
  		}
  	gsl_linalg_cholesky_decomp(cy);
  	jj = 0;
  	for (int iz = 0; iz != n; iz++)
  		for (int jz = 0; jz <= iz; jz++) {
  			if (iz != jz) temp[igroup*respno + 1 + jj] = gsl_matrix_get(cy, iz, jz); else
  				temp[igroup*respno + 1 + jj] = log(gsl_matrix_get(cy, iz, jz));
  			jj++;
  		}
  	gsl_matrix_free(cy); free(tau); free(temp);
  
  	for (int i = 0; i != restparsno; i++) bridge_sample[offset + i] = temp[i];
  	offset += restparsno;
  
  	bridge_sample[offset++] = llik;
  
  	if (offset - (ithread * IREP + ix)*n_bridge_store != n_bridge_store) Rprintf("Warnung: belege_bridge\n"); //std::cout << "Warnung: belege_bridge" << std::endl;
  }
  
  double loglik(std::vector<trial> daten, double *rhos, double *mu, double *beta, double *lambdas, double* lams, double *restpars, double *slams) {
  	double *x_for_all = 0; x_for_all = (double *)malloc(indi*kernpar * sizeof(double));
  	make_parameters_for_all(mu, lams, beta, x_for_all);
  
  	double *xsi = 0; xsi = (double *)malloc(indi*respno * sizeof(double));
  	double *pij = 0; pij = (double *)malloc(zweig * sizeof(double));
  
  
  	for (int t = 0; t != indi; t++) for (int r = 0; r != respno; r++) {
  		double rmu = restpars[t2group[t] * respno + r] + restpars[t*respno + r + alphaoff];
  		double rsig = sqrt(restpars[t + sigalphaoff]);
  		xsi[t*respno + r] = lnnorm(rmu / rsig);
  	}
  
  	double temp = 0;
  	int trialno = static_cast<int>(daten.size());
  	for (int x = 0; x != trialno; x++) {
  		double p;
  		trial one = daten[x];
  		int t = one.person; int r=cat2resp[one.category];
  		make_tij_for_one_trial_new(one, rhos, lambdas, lams, restpars, slams, pij);
  		make_pij_for_one_trial(one, x_for_all, pij, p);
  		temp += p-xsi[t*respno+r];
  	}
  	free(x_for_all);
  	free(xsi);
  	free(pij);
  	return temp;
  }
  
  void gibbs_and_monitor(std::vector<trial> daten,double* factor, double *mu, double *lams, double *beta, double *rhos, double *lambdas, int ntau, int *ntau_position, double *restpars, double *slams,
  	int *nnodes, int nz, int *nz_position, int offset, int n_all_parameters, double *parmon, gsl_rng *rst, int ithread,
  	bool save, double *sample, int n_bridge_store, double *bridge_sample) {
  	double *sig = 0;	sig = (double *)malloc((ifree + ilamfree)*(ifree + ilamfree) * sizeof(double));;
  	int *signs = 0; signs = (int *)malloc((ifree + ilamfree) * sizeof(int));
  	double *temp = 0; temp = (double *)malloc(n_all_parameters * sizeof(double));
  	double *taus = 0; taus = (double *)malloc(ntau * sizeof(double));
  
  
  
  	// n_all_parameters = ifree*igroup (mu) + 2*ilamfree*igroup (rhos) + (ifree+2*ilamfree)*(ifree+2*ilamfree+1)/2 (sig) + indi*ifree (betas) + indi*2*ilamfree (lambdas) + restparsno (restpars)
  	for (int i = 0; i != ireps; i++) {
  
  
  		bool xflag = (i == 0) && (offset == 0) ? true : false;
  		gibbs_full_cycle(daten,factor, mu, lams, beta, sig, rhos, lambdas, ntau, ntau_position, taus, nz, nz_position, nnodes, restpars, slams, xflag, rst);
  
  		int jj = -1;
  		for (int iz = 0; iz != ifree + ilamfree; iz++) signs[iz] = (lams[iz] >= 0) ? 1 : -1;
  		for (int j = 0; j != igroup; j++) for (int iz = 0; iz != ifree; iz++) temp[++jj] = mu[j*ifree + iz];
  
  
  		//	for (int iz=ifree;iz!=ifree+2*ilamfree;iz++) signs[iz]=(lams[iz-ifree+kernpar]>=0)?1:-1;
  
  		for (int j = 0; j != igroup; j++) for (int iz = 0; iz != ilamfree; iz++) temp[++jj] = rhos[j*ilamfree + iz];
  
  		for (int iz = 0; iz != ifree + ilamfree; iz++)
  			for (int jz = iz; jz != ifree + ilamfree; jz++) temp[++jj] = (iz != jz) ? SIG(iz, jz)*signs[iz] * signs[jz] / sqrt(SIG(iz, iz)*SIG(jz, jz)) :
  				sqrt(SIG(iz, iz))*signs[iz] * lams[iz];
  
  
  		for (int t = 0; t != indi; t++) 	for (int iz = 0; iz != ifree; iz++) temp[++jj] = lams[iz] * BETA(t, iz);
  
  		for (int t = 0; t != indi; t++) for (int iz = 0; iz != ilamfree; iz++) temp[++jj] = exp(lams[ifree + iz] * lambdas[t*ilamfree + iz]);
  
  		for (int it = 0; it != igroup * respno + 1; it++) temp[++jj] = restpars[it];
  		//		iz=-1;
  		int ioff = igroup * respno + 1;
  		for (int iz = 0; iz != respno; iz++) signs[iz] = (slams[iz] > 0) ? 1 : -1;
  		for (int ir = 0; ir != respno; ir++) for (int jr = ir; jr != respno; jr++) {
  			int iirr = ir * (respno - 1) - (ir*(ir - 1)) / 2 + ir, jjrr = jr * (respno - 1) - (jr*(jr - 1)) / 2 + jr, ijrr = ir * (respno - 1) - (ir*(ir - 1)) / 2 + jr;
  			temp[++jj] = (ir != jr) ? restpars[ioff + ijrr] * signs[ir] * signs[jr] / sqrt(restpars[ioff + iirr] * restpars[ioff + jjrr]) : sqrt(restpars[ioff + iirr])*signs[ir] * slams[ir];
  		}
  
  
  		for (int t = 0; t != indi; t++) for (int ir = 0; ir != respno; ir++) temp[++jj] = malpha(t, ir, restpars, slams);
  		for (int t = 0; t != indi; t++) temp[++jj] = restpars[sigalphaoff + t];
  
  
  		if ((save) && (i % THIN == 0)) {
  			int off = (ithread*IREP + i)*(n_all_parameters + 1);
  			for (int j = 0; j != n_all_parameters; j++) sample[off + j] = temp[j];
  			double llik = loglik(daten, rhos, mu, beta, lambdas, lams, restpars, slams);
  			sample[off + n_all_parameters] = llik;
  			if (for_bridge_flag) belege_bridge(ithread, i, n_bridge_store, bridge_sample, mu, lams, rhos, beta, lambdas, sig, restpars, slams, llik);
  		}
  		if ((i == 0) && (offset == 0)) for (int j = 0; j != 2 * n_all_parameters; j++) parmon[j] = 0.0;
  		double r = 1.0 / (i + offset + 1);
  #define PARMON(I,J) parmon[(I-1)*n_all_parameters+J]
  		for (int j = 0; j != n_all_parameters; j++) {
  			double dev = temp[j] - PARMON(1, j);
  			PARMON(2, j) += gsl_pow_2(dev)*(1.0 - r);
  			PARMON(1, j) += dev * r;
  		}
  	}
  	if (sig) free(sig);
  	if (signs) free(signs);
  	if (temp) free(temp);
  	if (taus) free(taus);
  
  }
  
  void initialize_new(std::vector <trial> daten, double *mu, double *lams, double* rhos, double *beta, double *lambdas, double *restpars, double *slams, gsl_rng *rst) {
  
  
  	for (int i = 0; i != ifree + ilamfree; i++) lams[i] = oneuni(rst) + 0.5;
  	for (int i = 0; i != respno; i++) slams[i] = oneuni(rst) + 0.5;
  
  	for (int ig = 0; ig != igroup; ig++)
  		for (int iz = 0; iz != ifree; iz++) {
  			mu[iz + ig * ifree] = 0.0; int ng = 0;
  			{
  				for (int t = 0; t != indi; t++) if (t2group[t] == ig) {
  					if (BETA(t, iz) < -4.0) BETA(t, iz) = -4.0 + 0.1*onenorm(rst); if (BETA(t, iz) > 4.0) BETA(t, iz) = 4.0 + 0.1*onenorm(rst);
  					mu[iz + ig * ifree] += BETA(t, iz); ng++;
  				}
  				mu[iz + ig * ifree] = (mu[iz + ig * ifree] + onenorm(rst)) / (ng + 1);
  			}
  		}
  
  	for (int t = 0; t != indi; t++) {
  
  		for (int iz = 0; iz != ifree; iz++) { BETA(t, iz) = (BETA(t, iz) - mu[iz + t2group[t] * ifree]) / lams[iz] + onenorm(rst); }
  	}
  
  	for (int ig = 0; ig != igroup; ig++)
  		for (int i = 0; i != ilamfree; i++) {
  			double mean = 0.0; int ng = 0;
  			for (int t = 0; t != indi; t++) if (t2group[t] == ig) {//if (LAMBDAS(t,pm,i) > 100.0) LAMBDAS(t,pm,i)=100;
  				mean += log(LAMBDAS(t, i)); ng++;
  			}
  			RHOS(ig, i) = exp((mean + onenorm(rst)) / (ng + 1));
  		}
  	for (int i = 0; i != ilamfree; i++)
  		for (int t = 0; t != indi; t++) {
  			double lax = LAMBDAS(t, i) / RHOS(t2group[t], i);
  			lax = log(lax) / lams[i + ifree] + onenorm(rst); LAMBDAS(t, i) = lax;
  		}
  
  	double *temp_rest = 0;  temp_rest = (double *)malloc(restparsno * sizeof(double));
  	for (int i = 0; i != restparsno; i++) { temp_rest[i] = restpars[i]; restpars[i] = 0.0; }
  
  
  	for (int ig = 0; ig != igroup; ig++) {
  		int ng = 0;
  		for (int t = 0; t != indi; t++) if (t2group[t] == ig) { restpars[ig*respno + 0] += temp_rest[t + 3]; ng++; }
  		double rsmu = restpars[ig*respno + 0];
  		for (int ir = 0; ir != respno; ir++) restpars[ig*respno + ir] = (rsmu + oneuni(rst)) / (ng + 1);
  	}
  	for (int t = 0; t != indi; t++) for (int r = 0; r != respno; r++) restpars[t*respno + r + alphaoff] = (temp_rest[t + 3] - restpars[t2group[t] * respno + r] + 0.01*onenorm(rst)) / slams[r];
  	for (int t = 0; t != indi; t++) if (temp_rest[t + indi + 3] < 0.4) restpars[t + sigalphaoff] = temp_rest[t + indi + 3] + 0.1*oneuni(rst); else restpars[t + sigalphaoff] = 0.5*oneuni(rst);
  	double gmu = 0; for (int t = 0; t != indi; t++) gmu += 1.0 / restpars[t + sigalphaoff]; restpars[1 + igroup * respno - 1] = 1.0 / ((gmu + 0.1*oneuni(rst)) / (indi + 1));
  	if (restpars[1 + igroup * respno - 1] > 0.05) restpars[1 + igroup * respno - 1] = 0.025 + 0.025*oneuni(rst);
  
  
  	if (DEBUG) {for (int it = 0; it != restparsno; it++) Rprintf("%g ", restpars[it]); /*std::cout << restpars[it] << " ";*/ Rprintf("\n");} //char xs; std::cin >> xs;
  
  
  
  	if (temp_rest) free(temp_rest);
  
  
  }
  
  
  void gibbs_times_new(std::vector<trial> daten, int *nnodes, int nz, int *nz_position, double *beta, int ntau, int *ntau_position,
  	gsl_rng *rst1, gsl_rng *rst2, gsl_rng *rst3, gsl_rng *rst4, gsl_rng *rst5, gsl_rng *rst6, gsl_rng *rst7, gsl_rng *rst8, gsl_rng *rst9, gsl_rng *rst10, gsl_rng *rst11, gsl_rng *rst12, gsl_rng *rst13, gsl_rng *rst14, gsl_rng *rst15, gsl_rng *rst16,
  	double* lambdas, double* restpars) {
  	bool do_burnin = (BURNIN > 0);
  	std::ofstream raus;
  	// std::ofstream raus_bridge;
  
  	double *lams = 0; lams = (double *)malloc((ifree + ilamfree) * sizeof(double));
  	double *slams = 0; slams = (double *)malloc(respno * sizeof(double));
  	double *rhos = 0; rhos = (double *)malloc(ilamfree*igroup * sizeof(double));
  	double *mu = 0; mu = (double *)malloc(ifree*igroup * sizeof(double));
  	double *factor = 0; factor = (double *)malloc(indi*respno * sizeof(double));
  
  
  	gsl_rng *xst;   xst = gsl_rng_alloc(T_rng);
  
  
  	// n_all_parameters = ifree * igroup + ilamfree * igroup + ((ifree + ilamfree)*(ifree + ilamfree + 1)) / 2 + indi * ifree + indi * ilamfree + restparsno;
  	int n_value_store = (indi + igroup + 1)*ifree + (indi + igroup + 1)*ilamfree + restparsno + respno + indi * respno;
  	int n_bridge_store = (n_value_store - indi * respno) + ((ifree + ilamfree + 1)*(ifree + ilamfree)) / 2 + 1;
  
  	double *parmon = 0; parmon = (double *)malloc(2 * n_all_parameters * sizeof(double));
  	for (int i = 0; i != 2 * n_all_parameters; i++) parmon[i] = 0.0;
  
  	double *valuestore = 0; valuestore = (double *)malloc(NOTHREADS*n_value_store * sizeof(double));
  	double *parmonstore = 0; parmonstore = (double *)malloc(NOTHREADS * 2 * n_all_parameters * sizeof(double));
  
  
  	int satemp = NOTHREADS * IREP*(n_all_parameters+1);
  	double *sample = 0; sample = (double *)malloc(satemp * sizeof(double));
  
  
  	// f�r bridge_sampler
  	double *bridge_sample = 0;
  	if (for_bridge_flag)
  		bridge_sample = (double *)malloc(NOTHREADS*IREP*n_bridge_store * sizeof(double));
  
  	//
  
  	// push f�r rstoreNOTHREADS
  	for (int ithread = 0; ithread != NOTHREADS; ithread++) {
  		push(ithread, n_value_store, n_all_parameters, factor, mu, lams, rhos, beta, lambdas, restpars, slams, valuestore, parmon, parmonstore);
  	}
  
  
  	for (int ithread = 0; ithread != NOTHREADS; ithread++) {
  		// pop f�r rst
  		switch (ithread + 1) {
  		case 1: gsl_rng_memcpy(xst, rst1); break;
  		case 2: gsl_rng_memcpy(xst, rst2); break;
  		case 3: gsl_rng_memcpy(xst, rst3); break;
  		case 4: gsl_rng_memcpy(xst, rst4); break;
  		case 5: gsl_rng_memcpy(xst, rst5); break;
  		case 6: gsl_rng_memcpy(xst, rst6); break;
  		case 7: gsl_rng_memcpy(xst, rst7); break;
  		case 8: gsl_rng_memcpy(xst, rst8); break;
  		case 9: gsl_rng_memcpy(xst, rst9); break;
  		case 10: gsl_rng_memcpy(xst, rst10); break;
  		case 11: gsl_rng_memcpy(xst, rst11); break;
  		case 12: gsl_rng_memcpy(xst, rst12); break;
  		case 13: gsl_rng_memcpy(xst, rst13); break;
  		case 14: gsl_rng_memcpy(xst, rst14); break;
  		case 15: gsl_rng_memcpy(xst, rst15); break;
  		case 16: gsl_rng_memcpy(xst, rst16); break;
  		}
  		pop(ithread, n_value_store, n_all_parameters, factor, mu, lams, rhos, beta, lambdas, restpars, slams, valuestore, parmon, parmonstore);
  		initialize_new(daten, mu, lams, rhos, beta, lambdas, restpars, slams, xst);
  		for (int i = 0; i != 2 * n_all_parameters; i++) parmon[i] = 0.0;
  		push(ithread, n_value_store, n_all_parameters, factor, mu, lams, rhos, beta, lambdas, restpars, slams, valuestore, parmon, parmonstore);
  		switch (ithread + 1) {
  			case 1: gsl_rng_memcpy(rst1, xst); break;
  			case 2: gsl_rng_memcpy(rst2, xst); break;
  			case 3: gsl_rng_memcpy(rst3, xst); break;
  			case 4: gsl_rng_memcpy(rst4, xst); break;
  			case 5: gsl_rng_memcpy(rst5, xst); break;
  			case 6: gsl_rng_memcpy(rst6, xst); break;
  			case 7: gsl_rng_memcpy(rst7, xst); break;
  			case 8: gsl_rng_memcpy(rst8, xst); break;
  			case 9: gsl_rng_memcpy(rst9, xst); break;
  			case 10: gsl_rng_memcpy(rst10, xst); break;
  			case 11: gsl_rng_memcpy(rst11, xst); break;
  			case 12: gsl_rng_memcpy(rst12, xst); break;
  			case 13: gsl_rng_memcpy(rst13, xst); break;
  			case 14: gsl_rng_memcpy(rst14, xst); break;
  			case 15: gsl_rng_memcpy(rst15, xst); break;
  			case 16: gsl_rng_memcpy(rst16, xst); break;
  		}
  	}
  	gsl_rng_free(xst);
  	free(lams); free(slams); free(rhos); free(mu); free(factor);
  	// main loop GIBBS
  RESTART:
  	double *xwbr = 0; xwbr = (double *)malloc(3 * n_all_parameters * sizeof(double));
  	for (int i = 0; i != 3 * n_all_parameters; i++) xwbr[i] = 0.0;
  	double rmax = 0.0;
  	int irun = -1;
  
  
  	bool save = false;
  	if (do_burnin) ireps = BURNIN; else ireps = IREP;
  	// double *complete_sample = 0;
  	// double *complete_bridge = 0;
  	int ioff = 0;
  WEITER: irun++;
  	int offset = irun * ireps;
  
  
  	// int maxThreads = std::thread::hardware_concurrency();
    // if (maxThreads == 0) Rprintf("Could not find out number of threads. Taking 2 threads.\n");
    // int suppThreads = maxThreads == 0 ? 2 : maxThreads;
    // int AmntOfThreads = suppThreads > NOTHREADS ? NOTHREADS : suppThreads;
    // int NperThread = N / AmntOfThreads;
    std::vector<std::thread> threads(NOTHREADS-1);
  
    /* starting threads while ... */
    for (int ithread = 0; ithread < NOTHREADS-1; ithread++) {
  		threads[ithread] = std::thread([&, ithread]() {
      // threads[ithread] = std::thread([ithread,rst1,rst2,rst3,rst4,rst5,rst6,rst7,rst8,rst9,rst10,rst11,rst12,rst13,rst14,rst15,rst16, save,sample,ntau,ntau_position,nz_position,n_value_store,daten,nnodes,nz,offset,valuestore,parmonstore,xwbr,rmax,free2kern,kern2free,comp,cat2tree, kernpar,kerncat,indi,zweig,branch,nodemax,ar,nodes_per_tree,tree_and_node2par,ilamfree, ifree,ipred,ndrin,drin,path_info,pfad_index,n_all_parameters,nppr,igroup,t2group,ireps,cat2resp,respno,alphaoff,sigalphaoff,restparsno,consts,bridge_sample,n_bridge_store]() {
  
  			//double rmax2 = 0;
  			double *mu = 0, *lams = 0, *slams = 0, *beta = 0, *rhos = 0, *lambdas = 0, *restpars = 0, *factor = 0, *parmon = 0;
  			
  			mu = (double *)malloc(ifree*igroup * sizeof(double));
  			lams = (double *)malloc((ifree + ilamfree) * sizeof(double));
  			slams = (double *)malloc(respno * sizeof(double));
  			beta = (double *)malloc(indi*ifree * sizeof(double));
  
  			parmon = (double *)malloc(2 * n_all_parameters * sizeof(double));
  			rhos = (double *)malloc(ilamfree*igroup * sizeof(double));
  			lambdas = (double *)malloc(indi*ilamfree * sizeof(double));
  			restpars = (double *)malloc(restparsno * sizeof(double));
  			factor = (double *)malloc(indi*respno * sizeof(double));
  
  			gsl_rng *rst;
  			rst = gsl_rng_alloc(T_rng);
  			switch (ithread + 1) {
  				case 1: gsl_rng_memcpy(rst, rst1); break;
  				case 2: gsl_rng_memcpy(rst, rst2); break;
  				case 3: gsl_rng_memcpy(rst, rst3); break;
  				case 4: gsl_rng_memcpy(rst, rst4); break;
  				case 5: gsl_rng_memcpy(rst, rst5); break;
  				case 6: gsl_rng_memcpy(rst, rst6); break;
  				case 7: gsl_rng_memcpy(rst, rst7); break;
  				case 8: gsl_rng_memcpy(rst, rst8); break;
  				case 9: gsl_rng_memcpy(rst, rst9); break;
  				case 10: gsl_rng_memcpy(rst, rst10); break;
  				case 11: gsl_rng_memcpy(rst, rst11); break;
  				case 12: gsl_rng_memcpy(rst, rst12); break;
  				case 13: gsl_rng_memcpy(rst, rst13); break;
  				case 14: gsl_rng_memcpy(rst, rst14); break;
  				case 15: gsl_rng_memcpy(rst, rst15); break;
  				case 16: gsl_rng_memcpy(rst, rst16); break;
  			}
  			pop(ithread, n_value_store, n_all_parameters, factor, mu, lams, rhos, beta, lambdas, restpars, slams, valuestore, parmon, parmonstore);
  			gibbs_and_monitor(daten,factor, mu, lams, beta, rhos, lambdas, ntau, ntau_position, restpars, slams, nnodes,
  				nz, nz_position, offset, n_all_parameters, parmon, rst, ithread, save, sample, n_bridge_store, bridge_sample);
  			push(ithread, n_value_store, n_all_parameters, factor, mu, lams, rhos, beta, lambdas, restpars, slams, valuestore, parmon, parmonstore);
  			switch (ithread + 1) {
  				case 1: gsl_rng_memcpy(rst1, rst); break;
  				case 2: gsl_rng_memcpy(rst2, rst); break;
  				case 3: gsl_rng_memcpy(rst3, rst); break;
  				case 4: gsl_rng_memcpy(rst4, rst); break;
  				case 5: gsl_rng_memcpy(rst5, rst); break;
  				case 6: gsl_rng_memcpy(rst6, rst); break;
  				case 7: gsl_rng_memcpy(rst7, rst); break;
  				case 8: gsl_rng_memcpy(rst8, rst); break;
  				case 9: gsl_rng_memcpy(rst9, rst); break;
  				case 10: gsl_rng_memcpy(rst10, rst); break;
  				case 11: gsl_rng_memcpy(rst11, rst); break;
  				case 12: gsl_rng_memcpy(rst12, rst); break;
  				case 13: gsl_rng_memcpy(rst13, rst); break;
  				case 14: gsl_rng_memcpy(rst14, rst); break;
  				case 15: gsl_rng_memcpy(rst15, rst); break;
  				case 16: gsl_rng_memcpy(rst16, rst); break;
  			}
  			int ido = 2;
  
  			while (ithread != atm_int);
  			if (ithread == 0) ido = 1;
  			int iter = offset + ireps;
  			r_statistic(ido, n_all_parameters, ithread, iter, parmon, xwbr, rmax);
  			atm_int++;
  
  			gsl_rng_free(rst);
  			free(mu); free(lams); free(slams); free(beta); free(parmon); free(rhos); free(lambdas); free(restpars); free(factor);
      });
    }
  
    /* the main thread also runs */
    {
    	double *mu = 0, *lams = 0, *slams = 0, *beta = 0, *rhos = 0, *lambdas = 0, *restpars = 0, *factor = 0;//, *parmon = 0;
    	mu = (double *)malloc(ifree*igroup * sizeof(double));
    	lams = (double *)malloc((ifree + ilamfree) * sizeof(double));
    	slams = (double *)malloc(respno * sizeof(double));
    	beta = (double *)malloc(indi*ifree * sizeof(double));
    
    	//parmon = (double *)malloc(2 * n_all_parameters * sizeof(double));
    	rhos = (double *)malloc(ilamfree*igroup * sizeof(double));
    	lambdas = (double *)malloc(indi*ilamfree * sizeof(double));
    	restpars = (double *)malloc(restparsno * sizeof(double));
    	factor = (double *)malloc(indi*respno * sizeof(double));
    
    	gsl_rng *rst;
    	rst = gsl_rng_alloc(T_rng);
    	switch (NOTHREADS) {
    		case 1: gsl_rng_memcpy(rst, rst1); break;
    		case 2: gsl_rng_memcpy(rst, rst2); break;
    		case 3: gsl_rng_memcpy(rst, rst3); break;
    		case 4: gsl_rng_memcpy(rst, rst4); break;
    		case 5: gsl_rng_memcpy(rst, rst5); break;
    		case 6: gsl_rng_memcpy(rst, rst6); break;
    		case 7: gsl_rng_memcpy(rst, rst7); break;
    		case 8: gsl_rng_memcpy(rst, rst8); break;
    		case 9: gsl_rng_memcpy(rst, rst9); break;
    		case 10: gsl_rng_memcpy(rst, rst10); break;
    		case 11: gsl_rng_memcpy(rst, rst11); break;
    		case 12: gsl_rng_memcpy(rst, rst12); break;
    		case 13: gsl_rng_memcpy(rst, rst13); break;
    		case 14: gsl_rng_memcpy(rst, rst14); break;
    		case 15: gsl_rng_memcpy(rst, rst15); break;
    		case 16: gsl_rng_memcpy(rst, rst16); break;
    	}
    	pop(NOTHREADS-1, n_value_store, n_all_parameters, factor, mu, lams, rhos, beta, lambdas, restpars, slams, valuestore, parmon, parmonstore);
    	gibbs_and_monitor(daten,factor, mu, lams, beta, rhos, lambdas, ntau, ntau_position, restpars, slams, nnodes,
    		nz, nz_position, offset, n_all_parameters, parmon, rst, NOTHREADS-1, save, sample, n_bridge_store, bridge_sample);
    	push(NOTHREADS-1, n_value_store, n_all_parameters, factor, mu, lams, rhos, beta, lambdas, restpars, slams, valuestore, parmon, parmonstore);
    	switch (NOTHREADS) {
    		case 1: gsl_rng_memcpy(rst1, rst); break;
    		case 2: gsl_rng_memcpy(rst2, rst); break;
    		case 3: gsl_rng_memcpy(rst3, rst); break;
    		case 4: gsl_rng_memcpy(rst4, rst); break;
    		case 5: gsl_rng_memcpy(rst5, rst); break;
    		case 6: gsl_rng_memcpy(rst6, rst); break;
    		case 7: gsl_rng_memcpy(rst7, rst); break;
    		case 8: gsl_rng_memcpy(rst8, rst); break;
    		case 9: gsl_rng_memcpy(rst9, rst); break;
    		case 10: gsl_rng_memcpy(rst10, rst); break;
    		case 11: gsl_rng_memcpy(rst11, rst); break;
    		case 12: gsl_rng_memcpy(rst12, rst); break;
    		case 13: gsl_rng_memcpy(rst13, rst); break;
    		case 14: gsl_rng_memcpy(rst14, rst); break;
    		case 15: gsl_rng_memcpy(rst15, rst); break;
    		case 16: gsl_rng_memcpy(rst16, rst); break;
    	}
    	int ido = 2;
    
    	while (NOTHREADS-1 != atm_int);
    	ido = 3;
    	int iter = offset + ireps;
    	r_statistic(ido, n_all_parameters, NOTHREADS-1, iter, parmon, xwbr, rmax);
    	atm_int.store(0);
    
    	gsl_rng_free(rst);
    	free(mu); free(lams); free(slams); free(beta); free(rhos); free(lambdas); free(restpars); free(factor);
    }
  
    /* join threads */
    for (int j = 0; j < NOTHREADS-1; j++) {
      threads[j].join();
    }
  
  
  
  // #ifdef _OPENMP
  // #pragma omp parallel for num_threads(NOTHREADS) ordered lastprivate(parmon) shared(rst1,rst2,rst3,rst4,rst5,rst6,rst7,rst8,rst9,rst10,rst11,rst12,rst13,rst14,rst15,rst16, save,sample,ntau,ntau_position,nz_position,n_value_store,daten,nnodes,nz,offset,valuestore,parmonstore,xwbr,rmax,free2kern,kern2free,comp,cat2tree, kernpar,kerncat,indi,zweig,branch,nodemax,ar,nodes_per_tree,tree_and_node2par,ilamfree, ifree,ipred,ndrin,drin,path_info,pfad_index,n_all_parameters,nppr,igroup,t2group,ireps,cat2resp,respno,alphaoff,sigalphaoff,restparsno,consts,bridge_sample,n_bridge_store)
  // #endif
  // 	for (int ithread = 0; ithread < NOTHREADS; ithread++) {
  // 		double *mu = 0, *lams = 0, *slams = 0, *beta = 0, *rhos = 0, *lambdas = 0, *restpars = 0, *factor = 0;
  // 		mu = (double *)malloc(ifree*igroup * sizeof(double));
  // 		lams = (double *)malloc((ifree + ilamfree) * sizeof(double));
  // 		slams = (double *)malloc(respno * sizeof(double));
  // 		beta = (double *)malloc(indi*ifree * sizeof(double));
  //
  // 		parmon = (double *)malloc(2 * n_all_parameters * sizeof(double));
  // 		rhos = (double *)malloc(ilamfree*igroup * sizeof(double));
  // 		lambdas = (double *)malloc(indi*ilamfree * sizeof(double));
  // 		restpars = (double *)malloc(restparsno * sizeof(double));
  // 		factor = (double *)malloc(indi*respno * sizeof(double));
  //
  //
  // 		gsl_rng *rst;
  // 		rst = gsl_rng_alloc(T_rng);
  // 		switch (ithread + 1) {
  // 			case 1: gsl_rng_memcpy(rst, rst1); break;
  // 			case 2: gsl_rng_memcpy(rst, rst2); break;
  // 			case 3: gsl_rng_memcpy(rst, rst3); break;
  // 			case 4: gsl_rng_memcpy(rst, rst4); break;
  // 			case 5: gsl_rng_memcpy(rst, rst5); break;
  // 			case 6: gsl_rng_memcpy(rst, rst6); break;
  // 			case 7: gsl_rng_memcpy(rst, rst7); break;
  // 			case 8: gsl_rng_memcpy(rst, rst8); break;
  // 			case 9: gsl_rng_memcpy(rst, rst9); break;
  // 			case 10: gsl_rng_memcpy(rst, rst10); break;
  // 			case 11: gsl_rng_memcpy(rst, rst11); break;
  // 			case 12: gsl_rng_memcpy(rst, rst12); break;
  // 			case 13: gsl_rng_memcpy(rst, rst13); break;
  // 			case 14: gsl_rng_memcpy(rst, rst14); break;
  // 			case 15: gsl_rng_memcpy(rst, rst15); break;
  // 			case 16: gsl_rng_memcpy(rst, rst16); break;
  // 		}
  // 		pop(ithread, n_value_store, n_all_parameters, factor, mu, lams, rhos, beta, lambdas, restpars, slams, valuestore, parmon, parmonstore);
  // 		gibbs_and_monitor(daten,factor, mu, lams, beta, rhos, lambdas, ntau, ntau_position, restpars, slams, nnodes,
  // 			nz, nz_position, offset, n_all_parameters, parmon, rst, ithread, save, sample, n_bridge_store, bridge_sample);
  // 		push(ithread, n_value_store, n_all_parameters, factor, mu, lams, rhos, beta, lambdas, restpars, slams, valuestore, parmon, parmonstore);
  // 		switch (ithread + 1) {
  // 			case 1: gsl_rng_memcpy(rst1, rst); break;
  // 			case 2: gsl_rng_memcpy(rst2, rst); break;
  // 			case 3: gsl_rng_memcpy(rst3, rst); break;
  // 			case 4: gsl_rng_memcpy(rst4, rst); break;
  // 			case 5: gsl_rng_memcpy(rst5, rst); break;
  // 			case 6: gsl_rng_memcpy(rst6, rst); break;
  // 			case 7: gsl_rng_memcpy(rst7, rst); break;
  // 			case 8: gsl_rng_memcpy(rst8, rst); break;
  // 			case 9: gsl_rng_memcpy(rst9, rst); break;
  // 			case 10: gsl_rng_memcpy(rst10, rst); break;
  // 			case 11: gsl_rng_memcpy(rst11, rst); break;
  // 			case 12: gsl_rng_memcpy(rst12, rst); break;
  // 			case 13: gsl_rng_memcpy(rst13, rst); break;
  // 			case 14: gsl_rng_memcpy(rst14, rst); break;
  // 			case 15: gsl_rng_memcpy(rst15, rst); break;
  // 			case 16: gsl_rng_memcpy(rst16, rst); break;
  // 		}
  // 		int ido = 2;
  //
  //
  // #ifdef _OPENMP
  // #pragma omp ordered
  // #endif
  // 		{
  // 			if (ithread == 0) ido = 1; if (ithread + 1 == NOTHREADS) ido = 3;
  // 			int iter = offset + ireps;
  // 			r_statistic(ido, n_all_parameters, ithread, iter, parmon, xwbr, rmax);
  // 		}
  // 		gsl_rng_free(rst);
  // 		free(mu); free(lams); free(slams); free(beta); free(rhos); free(lambdas); free(restpars); free(factor);
  // 	}
  
  
  	on_screen3(n_all_parameters, xwbr, parmon, beta, rmax, irun);
  
  
  
  	if (do_burnin) { do_burnin = false; goto RESTART; }
  	if ((save) && (rmax > RMAX)) { save = false; ioff = 0; if (complete_sample) free(complete_sample); if (complete_bridge) free(complete_bridge); goto WEITER; }
  	if (!(save) && (rmax > RMAX))  goto WEITER; else {
  		if (!(save)) {
  			complete_sample = (double *)malloc(SAMPLE_SIZE*(n_all_parameters+1) * sizeof(double));
  			save = true; raus.open(RAUS);
  			raus << std::setprecision(12);
  			raus << std::setw(5) << SAMPLE_SIZE << " " << n_all_parameters+1 << std::endl;
  			raus.close();
  			if (for_bridge_flag) {
  				// 	raus_bridge.open("raus_bridge");
  				// 	raus_bridge << std::setprecision(12);
  				// 	raus_bridge << std::setw(5) << SAMPLE_SIZE << " " << n_bridge_store << std::endl;
  				// 	raus_bridge.close();
  				complete_bridge = (double *)malloc(SAMPLE_SIZE*n_bridge_store * sizeof(double));
  			}
  			goto WEITER;
  		}
  
  	}
  #define SAMPLE(I,IP) sample[(I)*(n_all_parameters+1)+IP]
  #define COMPLETE_SAMPLE(S,I,P) complete_sample[(S)*(SAMPLE_SIZE/NOTHREADS)*(n_all_parameters+1) + (I)*(n_all_parameters+1)+P]
  #define BRIDGE_SAMPLE(I,IP) bridge_sample[(I)*n_bridge_store+IP]
  #define COMPLETE_BRIDGE(S,I,P) complete_bridge[(S)*(SAMPLE_SIZE/NOTHREADS)*n_bridge_store + (I)*n_bridge_store+P]
  
  	if (save) {
  		for (int is = 0; is != (NOTHREADS); is++)
  			for (int i = 0; i != IREP / THIN; i++) {
  				for (int j = 0; j != n_all_parameters+1; j++)
  					COMPLETE_SAMPLE(is, ioff*IREP / THIN + i, j) = SAMPLE(is*IREP + i * THIN, j);
  				if (for_bridge_flag)
  					for (int j = 0; j != n_bridge_store; j++)
  						COMPLETE_BRIDGE(is, ioff*IREP / THIN + i, j) = BRIDGE_SAMPLE(is*IREP + i * THIN, j);
  
  
  			}
  
  		ioff += 1;
  		if (ioff*NOTHREADS*IREP / THIN < SAMPLE_SIZE) goto WEITER;
  		raus.open(RAUS, std::ofstream::app);
  		raus << std::setprecision(12);
  		// if (for_bridge_flag) {
  		// 	raus_bridge.open("raus_bridge", std::ofstream::app);
  		// 	raus_bridge << std::setprecision(12);
  		// }
  		for (int is = 0; is != (NOTHREADS); is++)
  			for (int i = 0; i != (SAMPLE_SIZE / NOTHREADS); i++) {
  				for (int j = 0; j != n_all_parameters+1; j++)
  					raus << std::setw(20) << COMPLETE_SAMPLE(is, i, j);
  				raus << std::endl;
  				// if (for_bridge_flag) {
  				// 	for (int j = 0; j != n_bridge_store; j++)
  				// 		raus_bridge << std::setw(20) << COMPLETE_BRIDGE(is, i, j);
  				// 	raus_bridge << std::endl;
  				// }
  			}
  	}
  	raus.close();
  	// if (for_bridge_flag) raus_bridge.close();
  
  
  
  	//	if (lams) free(lams);
  	//	if (slams) free(slams);
  	//	if (mu) free(mu);
  	if (valuestore) free(valuestore);
  	if (parmonstore) free(parmonstore);
  	if (parmon) free(parmon);
  	if (xwbr) free(xwbr);
  	//	if (rhos) free(rhos);
  	//	if (factor) free(factor);
  	if (sample) free(sample);
  	if (bridge_sample) free(bridge_sample);
  	// if (complete_sample) free(complete_sample);
  	// if (complete_bridge) free(complete_bridge);
  }

}




namespace drtmpt {
  
  // sample path
  
  //omega - residual variance Gibbs sampler phase <=2
  void make_romega(gsl_vector* hampar, double* explambdas, double& omega, gsl_rng* rst) {
    //	double priordf = 2.0;
    //	double prioralpha = 0.0025;
    //	double priorbeta = 0.5;
    double s = 0.0;
    for (int t = 0; t != indi; t++) s += (phase >= 3) ? 1 / gsl_pow_2(explambdas[t]) : 1 / gsl_pow_2(gsl_vector_get(hampar, isigoff + t));
    
    s = s * priordf;
    double alpha = (indi * priordf) / 2.0 + prioralpha, beta = s / 2.0 + priorbeta;
    omega = gsl_ran_gamma(rst, alpha, 1.0 / beta);
  }
  
  //gibbs step individual deviations variance-covariance phase <=2
  void invwis(int cases, int nvar, double* xx, double* ssig, double* sigi, gsl_matrix* Ltminus, int nprior, double* xi, gsl_rng* rst) {
  // #define SSIG(I,J) ssig[I*nvar + J]
  // #define XX(T,J) xx[T*nvar + J]
  #define XB(T,J) xb[T*nvar+J]
    double* xb; if (!(xb = (double*)malloc(nvar * (cases + nvar + nprior) * sizeof(double)))) { Rprintf("Allocation failure\n"); }
    gsl_matrix* cx = gsl_matrix_alloc(nvar, nvar); // initializes to zero
    gsl_matrix_view XXX = gsl_matrix_view_array(xx, cases + nvar + nprior, nvar);
    
    
    gsl_matrix_view rows = gsl_matrix_submatrix(&XXX.matrix, 0, 0, cases, nvar);
    gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, &rows.matrix, 0.0, cx);
    
    
    
    gsl_vector_view diag = gsl_matrix_diagonal(cx);
    gsl_vector_view XI = gsl_vector_view_array(xi, nvar);
    gsl_blas_daxpy(2 * huang, &XI.vector, &diag.vector);
    
    gsl_linalg_cholesky_decomp1(cx);
    gsl_linalg_tri_lower_invert(cx);
    
    int nloop = nvar * (cases + nvar + nprior);
    for (int ih = 0; ih != nloop; ih++) xb[ih] = onenorm(rst);
    gsl_matrix_view XB = gsl_matrix_view_array(xb, nvar, cases + nvar + nprior);
    gsl_blas_dtrmm(CblasLeft, CblasLower, CblasTrans, CblasNonUnit, 1.0, cx, &XB.matrix);
    gsl_blas_dsyrk(CblasLower, CblasNoTrans, 1.0, &XB.matrix, 0.0, cx);
    
    gsl_matrix_view Xsigi = gsl_matrix_view_array(sigi, nvar, nvar);
    
    
    for (int j = 0; j != nvar; j++)
      for (int i = j; i != nvar; i++) {
        if (i != j) gsl_matrix_set(cx, j, i, gsl_matrix_get(cx, i, j));
      }
      gsl_matrix_memcpy(&Xsigi.matrix, cx);
    gsl_linalg_cholesky_decomp1(cx);
    
    
    if (phase < 3) {
      gsl_matrix_memcpy(Ltminus, cx);
    }
    
    gsl_linalg_cholesky_invert(cx);
    gsl_matrix_view tssig = gsl_matrix_view_array(ssig, nvar, nvar);
    gsl_matrix_memcpy(&tssig.matrix, cx);
    
    gsl_matrix_free(cx);
    free(xb);
    
  }
  
  //rgam - covariance motor times gibbs sampler phase<=2
  void make_rgam(gsl_vector* hampar, double* gam, double* gami, gsl_matrix* cr,  double* bi, gsl_rng* rst) {
    
    double* xr = 0;	if (!(xr = (double*)malloc((indi + respno + huang - 1) * (respno) * sizeof(double)))) { Rprintf("Allocation failure\n"); }
    gsl_vector_view t0 = gsl_vector_view_array(xr, (indi + respno + huang - 1) * respno);
    gsl_vector_view xxr = gsl_vector_subvector(&t0.vector, 0, indi * respno);
    gsl_vector_view t1 = gsl_vector_subvector(hampar, ilamoff, indi * respno);
    gsl_vector_memcpy(&xxr.vector, &t1.vector);
    
    invwis(indi, respno, xr, gam, gami, cr, huang - 1, bi, rst);
    double tausq = gsl_pow_2(taur);
    for (int i = 0; i != respno; i++) bi[i] = gsl_ran_gamma(rst, (respno + huang) * 0.5, 1.0 / (huang * dGAMI(i, i) + 1.0 / tausq));
    if (xr) free(xr);
  }
  
  //sig - covariance person random effects in a, v, and w, phase <=2
  void sample_sig(gsl_vector* hampar, double* sig, double* sigi, gsl_matrix* cx, double* ai, gsl_rng* rst) {
    
    
    double* xy = 0;	if (!(xy = (double*)malloc((indi + icompg + huang - 1) * icompg * sizeof(double)))) { Rprintf("Allocation failure\n"); }
    gsl_vector_view t0 = gsl_vector_view_array(xy, (indi + icompg + huang - 1) * icompg);
    gsl_vector_view xxy = gsl_vector_subvector(&t0.vector,0, indi * icompg);
    gsl_vector_view t1 = gsl_vector_subvector(hampar, igroup * icompg, indi * icompg);
    gsl_vector_memcpy(&xxy.vector, &t1.vector);
    
    invwis(indi, icompg, xy, sig, sigi, cx, huang - 1, ai, rst);
    double tausq = gsl_pow_2(taut);
    for (int i = 0; i != icompg; i++) ai[i] = gsl_ran_gamma(rst, (icompg + huang) * 0.5, 1.0 / (huang * dSIGI(i, i) + 1.0 / tausq));
    if (xy) free(xy);
  }
  
  double double_trunct(double lower, double upper, double plow, double help, gsl_rng* rst) {
    double result;
    if ((help > 0.01) || ((upper-lower) > 4.0)) result = gsl_cdf_tdist_Pinv(plow + oneuni(rst) * help, degf);
    else
    {
      double rho = 1.0, zz;
      double d = 0.0;
      if (lower * upper >= 0) {
        if (lower > 0) d = gsl_log1p(gsl_pow_2(lower) / degf);
        else
          if (upper < 0) d = gsl_log1p(gsl_pow_2(upper) / degf);
      }
      
      WEITER:   zz = oneuni(rst) * (upper - lower) + lower;
      if (lower * upper < 0) rho = exp(-0.5 * (degf + 1.0) * (gsl_log1p(gsl_pow_2(zz) / degf)));
      else
        rho = exp(0.5 * (degf + 1.0) * (d - gsl_log1p(gsl_pow_2(zz) / degf)));
      double u = oneuni(rst);
      if (u > rho) goto WEITER;
      result = zz;
    }
    
    return result;
  }
  
  void accept(int pfadlength, double* taua, double* tauc, double rest_a, double&  rest_c, double aa, double& ac) {
    gsl_vector_view ta = gsl_vector_view_array(taua, pfadlength);
    gsl_vector_view tc = gsl_vector_view_array(tauc, pfadlength);
    gsl_vector_memcpy(&tc.vector, &ta.vector);
    rest_c = rest_a;
    ac = aa;
    
  }
  
  void make_taus_met_hast(double rt, int pfadlength, int t, double* as, double* vs, double* ws, int* maps, int* low_or_up, double trmu, double tsig, double* taus, double& rest_old, ars_archiv& ars_store, gsl_rng* rst)
  {
    double* taun = (double*)malloc(pfadlength * sizeof(double));
    double lower = -trmu / tsig, upper = (rt - trmu) / tsig;
    double plow = gsl_cdf_tdist_P(lower, degf);
    double help = gsl_cdf_tdist_P(upper, degf) - plow;
    int x = pfadlength - 1;
    double ac = dwiener_d(taus[x] * low_or_up[x], as[x], vs[x], ws[x], accuracy);
    
    for (int i = 0; i != 100; i++) {
      
      //#pragma omp atomic
      //		MONITOR(0, 4)++;
      int icount = 0;
      NEW2: double tau_out = rt;
      icount++;
      if (icount > 1000000) goto END;
      double xt = double_trunct(lower, upper, plow, help, rst);
      xt = trmu + tsig * xt;
      tau_out -= xt;
      if (tau_out <= 0) goto NEW2;
      
      for (int x = 0; x != pfadlength - 1; x++)
      {
        taun[x] = make_rwiener(t, maps[x], (low_or_up[x] + 1) / 2, ars_store, rt, as[x], vs[x], ws[x], rst);
        tau_out -= taun[x];
        if (tau_out <= 0) {
          goto NEW2;
        }
      }
      int x = pfadlength - 1;
      taun[x] = tau_out;
      double aa = dwiener_d(tau_out * low_or_up[x], as[x], vs[x], ws[x], accuracy);
      
      double a = aa - ac;
      if ((aa > ac) || (log(oneuni(rst)) <= a)) {
        accept(pfadlength, taun, taus, xt, rest_old, aa, ac);
        //#pragma omp atomic
        //			MONITOR(1, 4)++;
      }
    }
    END:	free(taun);
  }
  
  // sample taus on path
  void make_taus_integrated(double rt, int pfadlength, int t, double* as, double* vs, double* ws, int* maps, int* low_or_up, double trmu, double tsig, double* taus, double& rest_old, ars_archiv& ars_store, gsl_rng* rst)
  {
    int icount=0;
    int acc = 0;
    
    double* tauc = (double*)malloc(pfadlength * sizeof(double));
    for (int x = 0; x != pfadlength; x++) tauc[x] = taus[x];
    double restc = rest_old;
    
    double cut = 0.0;
    if (rt < trmu) {
      cut = gsl_log1p(gsl_pow_2((rt - trmu) / tsig) / degf); // std::cout << "1";
    }
    else if (trmu < 0) {
      cut = gsl_log1p(gsl_pow_2((trmu) / tsig) / degf); //std::cout << "2";
    }
    
    double neuwc = -0.5 * (degf + 1.0) * (gsl_log1p(gsl_pow_2((restc - trmu) / tsig) / degf) - cut);
    
    //#pragma omp atomic
    //	MONITOR(0, 3)++;
    
    NEW2: rest_old = rt;
    
    for (int x = 0; x != pfadlength; x++)
    {
      taus[x] = make_rwiener(t, maps[x], (low_or_up[x] + 1) / 2, ars_store, rt, as[x], vs[x], ws[x], rst);
      rest_old -= taus[x];
      if (rest_old <= 0) {
        goto NEW2;
        
      }
    }
    
    double neuw_m = -0.5 * (degf + 1.0) * (gsl_log1p(gsl_pow_2((rest_old - trmu) / tsig) / degf) - cut);
    double u = log(oneuni(rst));
    if (u > neuw_m) {
      double a = neuw_m - neuwc;
      if ((neuw_m > neuwc) || (log(oneuni(rst)) <= a)) {
        accept(pfadlength, taus, tauc, rest_old, restc, neuw_m, neuwc);
        acc++;
      }
      icount++;
      if (((phase >= 3) || (icount <= 1000000)) && (icount <= 10000000)) goto NEW2; // || (phase3)) goto NEW;
      else {
        //			std::cout << "zu dumm";
        make_taus_met_hast(rt, pfadlength, t, as, vs, ws, maps, low_or_up, trmu, tsig, tauc, restc, ars_store, rst);
        for (int x = 0; x != pfadlength; x++) taus[x] = tauc[x];
        rest_old = restc;
        //#pragma omp atomic
        //			MONITOR(0, 2) += icount;
        //#pragma omp atomic
        //			MONITOR(1, 2) += acc;
        //#pragma omp atomic
        //			MONITOR(1, 3)++;
      }
    }
    for (int x = 0; x != pfadlength; x++) taus[x] *= low_or_up[x];
    free(tauc);
  }
  
  
  
  // sample taus not on path
  void make_taus_one_trial(trial one, int itrial, int ipath, double* tavw, int* tau_by_node, double* alltaus, ars_archiv& ars_store, gsl_rng* rst) {
    int itree = one.tree; int c = one.category, t = one.person;
    
    int npt = nodes_per_tree[itree];
    for (int n = 0; n != npt; n++) for (int pm = 0; pm != 2; pm++) {
      int low_or_up = (pm == 0) ? -1 : 1;
      if (dAR(c, ipath, n) != low_or_up)
      {
        int ia, iv, iw;
        ia = dTREE_AND_NODE2PAR(itree, n, 0);
        iv = dTREE_AND_NODE2PAR(itree, n, 1);
        iw = dTREE_AND_NODE2PAR(itree, n, 2);
        int m = dMAP(ia, iv, iw);
        double a = dTAVW(t, 0, ia);
        double v = dTAVW(t, 1, iv);
        double w = dTAVW(t, 2, iw);
        int tau_pos = dTAU_BY_NODE(itrial, n, pm);
        alltaus[tau_pos] = low_or_up * make_rwiener(one.person, m, pm, ars_store, GSL_POSINF, a, v, w, rst);
      }
    }
  }
  
  
  //number of nodes not on path
  void make_nips(const std::vector<trial> & daten, int* paths, int* nips) {
    
    for (int t = 0; t != indi; t++) for (int pm = 0; pm != 2; pm++) for (int m = 0; m != no_patterns; m++)	dNIPS(t, pm, m) = 0;
    
    for (int x = 0; x != datenzahl; x++) {
      int itree = daten[x].tree, t = daten[x].person, c = daten[x].category;
      int npt = nodes_per_tree[itree];
      for (int n = 0; n != npt; n++) {
        int im = dTREE_AND_NODE2MAP(itree, n);
        for (int pm = 0; pm != 2; pm++) {
          int low_or_up = (pm == 0) ? -1 : 1;
          if (dAR(c, paths[x], n) != low_or_up) dNIPS(t, pm, im)++;
        }
      }
    }
  }
  
  //update number of nodes not on path after sampling new path
  void update_nips(trial one, int path, int newpath, int* nips) {
    int t = one.person, c = one.category, itree = one.tree;
    int ndc = dNCDRIN(c);
    for (int in = 0; in != ndc; in++) {
      int n = dCDRIN(c, in, 0), pm = dCDRIN(c, in, 1);
      int low_or_up = (pm == 0) ? -1 : 1;
      if (dAR(c, path, n) != dAR(c, newpath, n)) {
        int m = dTREE_AND_NODE2MAP(itree, n);
        if (dAR(c, newpath, n) == low_or_up) dNIPS(t, pm, m)--;
        // neuer Pfad drin, daher alter nicht
        else if (dAR(c, path, n) == low_or_up) dNIPS(t, pm, m)++;
        // alter Pfad drin, daher neuer nicht
        
      }
    }
  }
  
  
  //path probabilities Carlin & Chibb 1995
  double fypgtau_and_path(int pfadlength, double* as, double* vs, double* ws, double trmu, double tsig, double* taus, double delta) {
    if (delta < 0) return  GSL_NEGINF;
    else
    {
      double  prob = 0.0;
      for (int in = 0; in != pfadlength; in++) {
        int pm = (taus[in] > 0) ? 1 : 0;
        prob += logprob_upperbound(pm, as[in], vs[in], ws[in]);
      }
      prob += log_tdist_pdf(trmu, tsig, degf, delta);
      return prob;
    }
  }
  
  //data augmentation: sample path and taus
  void make_path(trial one, int* nips, int itrial, int& path, gsl_vector* hampar, double* tavw, double* tlams, double* explambda, double* alltaus, double* rest, ars_archiv& ars_store, gsl_rng* rst) {
    int c = one.category, k = branch[c], tree = one.tree, t = one.person, r = cat2resp[c], new_path = -1;
    double rt = one.rt / 1000.0;
    
    double trmu = (phase >= 3)? tlams[t * respno + r]: gsl_vector_get(hampar, irmuoff + t2group[t] * respno + r) +
      gsl_vector_get(hampar, ilamoff + t * respno + r),
      tsig = (phase >= 3) ? explambda[t] : gsl_vector_get(hampar, isigoff + t);
    int pfadma = pfadmax[c];
    
    
    double* as = 0; if (!(as = (double*)malloc(pfadma * sizeof(double)))) { Rprintf("Allocation failure\n"); }
    double* vs = 0; if (!(vs = (double*)malloc(pfadma * sizeof(double)))) { Rprintf("Allocation failure\n"); }
    double* ws = 0; if (!(ws = (double*)malloc(pfadma * sizeof(double)))) { Rprintf("Allocation failure\n"); }
    int* low_or_up = 0; if (!(low_or_up = (int*)malloc(pfadma * sizeof(int)))) { Rprintf("Allocation failure\n"); }
    int* maps = 0; if (!(maps = (int*)malloc(pfadma * sizeof(int)))) { Rprintf("Allocation failure\n"); }
    double* taus = 0; if (!(taus = (double*)malloc(pfadma * sizeof(double)))) { Rprintf("Allocation failure\n"); }
    
    double* pj = 0; if (!(pj = (double*)malloc(k * sizeof(double)))) { Rprintf("Allocation failure\n"); }
    double* delta = 0; if (!(delta = (double*)malloc(k * sizeof(double)))) { Rprintf("Allocation failure\n"); }
    
    
    int akt_pflength = dNDRIN(c, path);
    
    // taus f�r path belegen:
    delta[path] = rt;
    for (int in = 0; in != akt_pflength; in++) {
      int n = dDRIN(c, path, in);
      low_or_up[in] = dAR(c, path, n);
      
      int ia = dTREE_AND_NODE2PAR(tree, n, 0); as[in] = dTAVW(t, 0, ia);
      int iv = dTREE_AND_NODE2PAR(tree, n, 1); vs[in] = dTAVW(t, 1, iv);
      int iw = dTREE_AND_NODE2PAR(tree, n, 2); ws[in] = dTAVW(t, 2, iw);
      maps[in] = dMAP(ia, iv, iw);
      int pm = (1 + dAR(c, path, n)) / 2;
      taus[in] = fabs(alltaus[dTAU_BY_NODE(itrial, n, pm)]);
      delta[path] -= taus[in];
      //			if (delta[path] < 0) std::cout << "aha";
    }
    make_taus_integrated(rt, akt_pflength, t, as, vs, ws, maps, low_or_up, trmu, tsig, taus, delta[path], ars_store, rst);
    for (int in = 0; in != akt_pflength; in++) {
      int n = dDRIN(c, path, in), pm = (1 + dAR(c, path, n)) / 2;
      alltaus[dTAU_BY_NODE(itrial, n, pm)] = taus[in];
    }
    make_taus_one_trial(one, itrial, path, tavw, tau_by_node, alltaus, ars_store, rst);
    
    if (k > 1) {
      
      pj[path] = fypgtau_and_path(akt_pflength, as, vs, ws, trmu, tsig, taus, delta[path]);
      
      double pjsum = GSL_NEGINF;
      // Pfadwahl multinomial probabilities
      for (int j = 0; j != k; j++) if (j != path) {
        delta[j] = rt;
        int ncj = dNDRIN(c, j);
        for (int in = 0; in != ncj; in++) {
          int n = dDRIN(c, j, in);
          int pm = (1 + dAR(c, j, n)) / 2;
          int ip = dTREE_AND_NODE2PAR(tree, n, 0); as[in] = dTAVW(t, 0, ip);
          ip = dTREE_AND_NODE2PAR(tree, n, 1); vs[in] = dTAVW(t, 1, ip);
          ip = dTREE_AND_NODE2PAR(tree, n, 2); ws[in] = dTAVW(t, 2, ip);
          int help = dTAU_BY_NODE(itrial, n, pm);
          taus[in] = alltaus[help];
          delta[j] -= fabs(taus[in]);
        }
        pj[j] = fypgtau_and_path(dNDRIN(c, j), as, vs, ws, trmu, tsig, taus, delta[j]);
      }
      for (int j = 0; j != k; j++) pjsum = logsum(pjsum, pj[j]);
      double u = pjsum + log(oneuni(rst));
      new_path = 0;
      double temp = pj[new_path];
      while (u > temp) {
        new_path++;
        temp = logsum(temp, pj[new_path]);
      }
  
      if (path != new_path) update_nips(one, path, new_path, nips);
    }
    else new_path = 0;
    path = new_path; rest[itrial] = delta[path];
    
    
    if (as) free(as);
    if (vs) free(vs);
    if (ws) free(ws);
    if (low_or_up) free(low_or_up);
    if (maps) free(maps);
    if (taus) free(taus);
    if (pj) free(pj);
    if (delta) free(delta);
  }

}
