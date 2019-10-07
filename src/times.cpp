// authors: Christoph Klauer and Raphael Hartmann
#include "rts.h"


#define NNODES(I,J) nnodes[I*kernpar+J]
#define AGAMMA(PM,IZ) agamma[PM*ilamfree+IZ]
#define BGAMMA(PM,IZ) bgamma[PM*ilamfree+IZ]


void make_idaten(vector<trial> daten, int *idaten) {
#define IDATEN(I,J) idaten[I*kerncat + J]
	trial one;
	for (int i = 0; i != indi * kerncat; ++i) idaten[i] = 0;
	for (int i = 0; i != int(daten.size()); ++i) {
		one = daten[i];
		IDATEN(one.person, one.category)++;
	}
}


double oneexp(double lambda, gsl_rng *rst) {
	double a, help[1];
	if (lambda > 0) {
		a = 1.0 / lambda;
		help[0] = gsl_ran_exponential(rst, a);
	}
	else {
		help[0] = DBL_MAX;
		if (DEBUG) Rprintf("oneexp");
	}
	return help[0];
}

double oneuni(gsl_rng *rst) {
	return gsl_rng_uniform_pos(rst);
}

double rexp(double x) {
	double result;
	if (x <= 700.0)
		result = exp(x);
	else
		result = exp(700.00);
	return result;
}
//end function rexp

double truncexp(double lambda, double upper, gsl_rng *rst) {
	double u; double help;
	if (fabs(lambda*upper) <= 1.0e-5) {
	NEW:    u = oneuni(rst); double z = upper * oneuni(rst);
		if ((lambda > 0) && (u >= exp(-lambda * z)))  goto NEW;
		if ((lambda < 0) && (u >= exp(lambda*(upper - z)))) goto NEW;
		help = z;
	}
	else {
		u = oneuni(rst);
		double temp = log(u) - lambda * upper;
		if (temp >= 700.0) help = temp / (-lambda);
		else help = -gsl_log1p(-u * (1.0 - exp(-lambda * upper))) / lambda;
	}
	return help;
}

void make_parameters_for_all(double *mu, double *lams, double *beta, double *x_for_all) {
	for (int t=0;t!=indi;t++) for (int i=0;i!=kernpar;i++)
		x_for_all[t*kernpar+i]=equation(t,i, mu,lams,beta);
}


double logsum(double xa, double xb) {
	double temp;
	if (xa == GSL_NEGINF) return xb;
	if (xb == GSL_NEGINF) return xa;
	if (xa > xb) temp = xa + gsl_log1p(exp(xb - xa));
	else if (xb > -DBL_MAX)
		temp = xb + gsl_log1p(exp(xa - xb));
	else temp = -DBL_MAX;
	return temp;
}


#define LAMBDAS(T,PM,IZ) lambdas[T*ilamfree+IZ] //PM=0 negativ; PM=1 positiv
#define TREE_AND_NODE2PAR(ITREE,R) tree_and_node2par[ITREE*nodemax+R]
#define DRIN(J,K,X) drin[J*zweig*nodemax+K*nodemax + X]
#define NDRIN(J,K) ndrin[J*zweig+K]
#define AR(I,J,R) ar[I*zweig*nodemax + J*nodemax + R]

double logexgaussian(double lam, double mu, double sd, double t)
{
	double temp;
	temp = log(lam) + lam * (mu + 0.5*lam*gsl_pow_2(sd) - t);
	double xsi = -mu - lam * gsl_pow_2(sd); xsi /= sd; double b = lnnorm(t / sd + xsi);
	double dif = lnnorm(xsi) - b;
	if (dif < 0) temp += b + gsl_log1p(-exp(dif)); else
	{//std::cout << "hoppla dif logexgaussian" << std::endl;
		temp = GSL_NEGINF;
	}
	return temp;
}



void make_pij_for_one_trial(trial one, double *x_for_all, double *pij, double &pj) {
	// berechnet  p


#define X_FOR_ALL(T,IP) x_for_all[T*kernpar+IP]


	double d0ha;
	int j = one.category, t = one.person, itree = one.tree;

	//		pj = 0.0;
	for (int k = 0; k != branch[j]; k++) {
		//			 pij[k]=1.0;
		for (int ir = 0; ir != NDRIN(j, k); ir++) {
			int r = DRIN(j, k, ir); int ix = AR(j, k, r); int ip = TREE_AND_NODE2PAR(itree, r);
			d0ha = (ix > 0) ? lnnorm(X_FOR_ALL(t, ip)) : lnnorm(-X_FOR_ALL(t, ip));
			pij[k] = pij[k] + d0ha;
		}
		pj = (k == 0) ? pij[0] : logsum(pj, pij[k]);
	}

	if (!(pj == pj) || !(std::isfinite(pj))) {
		if (branch[j] > 1) {
		if (DEBUG) Rprintf("pj is %g\n", pj);
		if (DEBUG) for (int k = 0; k != branch[j]; k++) Rprintf("%20g\n", pij[k]);
		// std::cout << "pj" << " is " << pj << std::endl;
		}
		pj = -sqrt(DBL_MAX);
		for (int k = 0; k != branch[j]; k++) pij[k] = log(1.0 / (1.0*branch[j])) - pj;
	}
}

int make_path_for_one_trial(int branchno, double *pij, double pj, gsl_rng *rst) {
	int exit_status = 0;
	int help = 0; double temp;
	if (branchno > 1) {
		double u = log(oneuni(rst)) + pj; temp = pij[help];
		while (u > temp) {
			help++;
			if (DEBUG) {if (help > branchno - 1) { Rprintf("Achtung non-multinomial"); exit_status = -1; }}
			temp =logsum(temp, pij[help]);
		}
	}
	return(help);

}

void make_zs_one_trial(trial one, int itrial, int ipath, double *mu, double *lams, double *beta, int *nz_position, double *z, gsl_rng *rst)
{

#define NZ_POSITION(X,J) nz_position[X*nodemax+J]

	double be; int t = one.person; int itree = one.tree; int j = one.category;
	for (int r = 0; r != nodes_per_tree[itree]; r++) {
		int ip = TREE_AND_NODE2PAR(itree, r);
		if (comp[ip]) {
			be = equation(t, ip, mu, lams, beta);  int z_pos = NZ_POSITION(itrial, r);
			if (AR(j, ipath, r) > 0) z[z_pos] = truncnorm(be, rst);
			if (AR(j, ipath, r) < 0) z[z_pos] = -truncnorm(-be, rst);
			if (AR(j, ipath, r) == 0)z[z_pos] = onenorm(rst) + be;
		}
	}
}

double double_truncnorm(double lower, double upper, gsl_rng *rst) {
	double result;
	int icount = 0;
	double plow = gsl_cdf_ugaussian_P(lower);
	double help = gsl_cdf_ugaussian_P(upper) - plow;
	if (help > 0.1)
	{
		double u = oneuni(rst); result = gsl_cdf_ugaussian_Pinv(plow + u * help);
	}
	else
	{
		double zz;
		if ((lower > 0.0) && (upper - lower > 0.4))
		{
			do { result = lower + truncnorm(-lower, rst); } while (result >= upper);
		}
		else
			if ((upper < 0.0) && (upper - lower) > 0.4)
			{
				do { result = upper - truncnorm(upper, rst); if (icount > 0) if (DEBUG) Rprintf("%d upper %g\n", icount++, upper); /*std::cout << icount++ << "upper " << upper << std::endl;*/ } while (result <= lower);
			}
			else
			{
				double rho;
			WEITER:   zz = oneuni(rst)*(upper - lower) + lower;
				if (lower*upper < 0) {
					rho = exp(-0.5*gsl_pow_2(zz));
				} else if (lower > 0) {
					rho = exp(0.5*(gsl_pow_2(lower) - gsl_pow_2(zz)));
				} else if (upper < 0) {
					rho = exp(0.5*(gsl_pow_2(upper) - gsl_pow_2(zz)));
				}
				double u = oneuni(rst);
				if (u > rho) goto WEITER;
				result = zz;
			}
	}
	return result;
}
