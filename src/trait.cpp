// authors: Christoph Klauer
#include "rts.h"

extern "C" {
#include "recurse.h"
}


void invwis(int cases, int nvar, double *xx, double *ssig, double *sigi, double eps, gsl_rng *rst) {
#define SSIG(I,J) ssig[I*nvar + J]
#define XX(T,J) xx[T*nvar + J]
#define XB(T,J) xb[T*nvar+J]
	double *xb;

	gsl_matrix *cx = gsl_matrix_alloc(nvar, nvar);
	//gsl_vector *s = gsl_vector_alloc(nvar);


	xb = (double *)malloc(nvar*(cases + nvar + 1 + pr_df_add_inv_wish) * sizeof(double));


	for (int j = 0; j != nvar; j++)
		for (int i = j; i != nvar; i++) {
			SSIG(i, j) = 0.0;
			for (int t = 0; t != cases; t++) {
				SSIG(i, j) += XX(t, i)*XX(t, j);
				//only lower triangular is needed
			}
			if (i == j) SSIG(i, i) += eps;
			gsl_matrix_set(cx, i, j, SSIG(i, j));
			if (i != j) gsl_matrix_set(cx, j, i, SSIG(i, j));
		}

	gsl_linalg_cholesky_decomp(cx);
	gsl_linalg_tri_lower_invert_dings(cx);

	for (int ih = 0; ih != nvar * (cases + nvar + 1 + pr_df_add_inv_wish); ih++) xb[ih] = onenorm(rst);


	for (int t = 0; t != cases + nvar + 1 + pr_df_add_inv_wish; t++) {
		for (int j = 0; j != nvar; j++) {
			XX(t, j) = 0.0;
			for (int i = j; i != nvar; i++) XX(t, j) += gsl_matrix_get(cx, i, j) * XB(t, i);
			//	XX(t,j)*=gsl_vector_get(s, j);
		}
	}

	for (int j = 0; j != nvar; j++)
		for (int i = j; i != nvar; i++) {
			double temp = 0.0;
			for (int t = 0; t != cases + nvar + 1 + pr_df_add_inv_wish; t++) {
				temp += XX(t, i)*XX(t, j);
			}
			gsl_matrix_set(cx, i, j, temp);
			if (i != j) gsl_matrix_set(cx, j, i, temp);
			sigi[i*nvar + j] = temp;
			if (i != j) sigi[j*nvar + i] = temp;
		}

	gsl_linalg_cholesky_decomp(cx);
	gsl_linalg_cholesky_invert(cx);


	for (int j = 0; j != nvar; j++)
		for (int i = j; i != nvar; i++) {
			//	SSIG(i, j) = gsl_vector_get(s, i)*gsl_vector_get(s, j)*gsl_matrix_get(cx, i, j);
			SSIG(i, j) = gsl_matrix_get(cx, i, j);
			if (i != j) SSIG(j, i) = SSIG(i, j);
		}

	gsl_matrix_free(cx);
	//gsl_vector_free(s);
	free(xb);

}

double onenorm(gsl_rng *rst) {
	return gsl_ran_ugaussian(rst);
}


double truncnorm(double b, gsl_rng *rst) {
	double temp;
	if (b >= 0.0) {
		do  temp = onenorm(rst);  while (temp < -b);
		temp += b;
	}
	else temp = gsl_ran_ugaussian_tail(rst, -b) + b;
	return temp;
}

void bayesreg(int n, double *mean, double *sigma, double *out, gsl_rng *rst) {
#define NTIG(I,J) ntig[I*n+J]


	double *xb = 0, *hout = 0, *ntig = 0;
	xb = (double *)malloc(n * sizeof(double));
	hout = (double *)malloc(n * sizeof(double));
	ntig = (double *)malloc(n*n * sizeof(double));

	gsl_matrix *cx = gsl_matrix_alloc(n, n);
	//gsl_vector *s  = gsl_vector_alloc(n);

	for (int j = 0; j != n; j++)
		for (int i = j; i != n; i++) {
			gsl_matrix_set(cx, i, j, sigma[i*n + j]);
			if (i != j) gsl_matrix_set(cx, j, i, sigma[i*n + j]);
		}
	gsl_linalg_cholesky_decomp(cx);
	//	gsl_matrix_memcpy(cy, cx);
	gsl_linalg_tri_lower_invert_dings(cx);

	for (int ih = 0; ih != n; ih++) xb[ih] = onenorm(rst);

	for (int i = 0; i != n; i++) out[i] = hout[i] = 0.0;
	for (int j = 0; j != n; j++) {
		for (int i = j; i != n; i++) hout[j] += gsl_matrix_get(cx, i, j)*xb[i];
		//hout[j] *= gsl_vector_get(s, j);
	}

	// cx = Lt-1 * L-1 = cxt * cx

	for (int i = 0; i != n; i++)
		for (int k = i; k != n; k++) {
			NTIG(i, k) = 0.0;
			for (int j = k; j != n; j++)
				NTIG(i, k) += gsl_matrix_get(cx, j, i)*gsl_matrix_get(cx, j, k);
			//	NTIG(i, k) = gsl_vector_get(s, i)*gsl_vector_get(s, k)*NTIG(i, k);
			if (i != k) NTIG(k, i) = NTIG(i, k);
		}




	for (int i = 0; i != n; i++) {
		for (int j = 0; j != n; j++) out[i] += NTIG(i, j)*mean[j];
		out[i] += hout[i];
	}

	if (xb) free(xb);
	if (hout) free(hout);
	free(ntig);
	gsl_matrix_free(cx);
	// gsl_vector_free(s);
}


#define BETA(I,J) beta[I*ifree+J]


double equation(int t, int ip, double *mu, double *lams, double *beta) {

	double xmu;
	if (comp[ip]) {
		int iz = kern2free[ip];
		xmu = (igroup > 1) ? mu[iz + ifree * t2group[t]] : mu[iz];
		xmu += lams[iz] * BETA(t, iz);
	}
	else xmu = consts[ip];
	return xmu;
}


void make_pij_for_individual(double *x, double *pij, double *pj) {
	// berechnet  pj, und pij fuer Individuum t, entspricht altem estimate
#define PIJ(I,J) pij[I*zweig+J]
#define A(I,J,K) a[I*zweig*kernpar + J*kernpar + K]
#define B(I,J,K) b[I*zweig*kernpar + J*kernpar + K]
#define DRIN(J,K,X) drin[J*zweig*nodemax+K*nodemax + X]
#define NDRIN(J,K) ndrin[J*zweig+K]
#define AR(I,J,K) ar[I*zweig*nodemax + J*nodemax + K]
#define TREE_AND_NODE2PAR(I,J) tree_and_node2par[I*nodemax+J]

	double d0ha;

	for (int j = 0; j != kerncat; j++) {
		pj[j] = 0.0;
		for (int k = 0; k != branch[j]; k++) {
			PIJ(j, k) = 1.0;
			for (int xr = 0; xr != NDRIN(j, k); xr++) {
				int r = DRIN(j, k, xr);
				int ia = AR(j, k, r); int ip = TREE_AND_NODE2PAR(cat2tree[j], r);
				d0ha = (ia > 0) ? x[ip] : 1 - x[ip];
				PIJ(j, k) *= d0ha;
			}
			pj[j] += PIJ(j, k);
		}
	}
	for (int j = 0; j != kerncat; j++)
		if (pj[j] != 0.0) for (int k = 0; k != branch[j]; k++) PIJ(j, k) = PIJ(j, k) / pj[j];
		else for (int k = 0; k != branch[j]; k++) PIJ(j, k) = 1.0 / (1.0*branch[j]);
}


#define NNODES(T,IP)  nnodes[T*kernpar + IP]


void make_mu(double *mu, double *lams, double *beta, int *nnodes, double *z, gsl_rng *rst) {
	double *mean = 0;	mean = (double *)malloc((igroup) * sizeof(double));
	double *xtx = 0;	xtx = (double *)malloc((igroup) * sizeof(double));

	int jj = -1;
	for (int iz = 0; iz != ifree; iz++) {
		int ip = free2kern[iz];
		for (int i = 0; i != igroup; i++) { mean[i] = 0.0; xtx[i] = 0.0; }
		for (int t = 0; t != indi; t++) {
			int itg = t2group[t];
			xtx[itg] += NNODES(t, ip);
			double rest = lams[iz] * BETA(t, iz);
			for (int j = 0; j != NNODES(t, ip); j++) {
				jj++;
				double xh = z[jj] - rest;
				mean[itg] += xh;
			}
		}

		for (int ix = 0; ix != igroup; ix++) {
			xtx[ix] += PRIOR;
			mu[ix*ifree + iz] = mean[ix] / xtx[ix] + onenorm(rst) / sqrt(xtx[ix]);
		}
	}
	if (xtx) free(xtx);
	if (mean) free(mean);
}

void make_lams(double *mu, double *lams, double *beta, int *nnodes, double *z, gsl_rng *rst) {

	double *w = 0;	w = (double *)malloc(ifree * sizeof(double));
	double *u = 0;	u = (double *)malloc(ifree * sizeof(double));



	int jj = -1;
	for (int iz = 0; iz != ifree; iz++) {
		int ip = free2kern[iz];
		w[iz] = 0.0; u[iz] = PRIOR;
		for (int t = 0; t != indi; t++) {
			double uiz = 0, wiz = 0;
			double be = equation(t, ip, mu, lams, beta) - BETA(t, iz)*lams[iz];
			uiz += NNODES(t, ip);
			for (int j = 0; j != NNODES(t, ip); j++) {
				wiz += (z[++jj] - be);
			}
			u[iz] += uiz * gsl_pow_2(BETA(t, iz));
			w[iz] += wiz * BETA(t, iz);
		}
	}

	for (int iz = 0; iz != ifree; iz++) {
		if (u[iz] <= 0) u[iz] = DBL_MIN;
		lams[iz] = (PRIOR + w[iz]) / u[iz] + onenorm(rst) / sqrt(u[iz]);
	}

	if (w) free(w);
	if (u) free(u);
}
