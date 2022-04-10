#include "rts.h"
#include <gsl/gsl_sf.h>

#define TREE_AND_NODE2PAR(ITREE,R) tree_and_node2par[ITREE*nodemax+R]
#define DRIN(J,K,X) drin[J*zweig*nodemax+K*nodemax + X]
#define NDRIN(J,K) ndrin[J*zweig+K]
#define AR(I,J,R) ar[I*zweig*nodemax + J*nodemax + R]
#define PFAD_INDEX(J,K) pfad_index[J*zweig+K]


int succ(int x, int k) {
	int temp;
	if (x == k - 1) temp = x + 2;
	else temp = x + 1;
	return temp;
}


void init_step(int k, int b, int *iz, int l_aktuell) {
	int x = -1; x = succ(x, k);
	iz[x] = l_aktuell;
	iz[k] = 0;
	for (int i = succ(x, k); i != b; i++) iz[i] = 0;
}

bool step(int k, int a, int *iz, int ll) {
	int x = -1; x = succ(x, k); int y = 0;
	while ((succ(x, k) < a) && ((iz[x]==0) || (iz[succ(x, k)] == ll))) {
		y += iz[x];
		x = succ(x, k);
	}
	bool temp = (succ(x, k) < a);
	if (temp) {
		y += iz[x];
		iz[succ(x, k)] += 1;
		init_step(k, succ(x, k), iz, y - 1);
	}
	return temp;

}

void logPhikl(int k, int a, vector<int> r, double *lams, double *loglams, int l, double &hypoplus, double &hypominus) {
	int *iz = 0; iz = (int *)malloc(a * sizeof(int));

	hypoplus = hypominus = GSL_NEGINF;
//	double zwi = 0;
	bool fplus = true, fminus = true;
	init_step(k, a, iz, l - 1);
	bool temp = true;
	while (temp)
	{
		if (temp) {
			double oneterm = 0.0;
			int sign = ((l - 1) % 2 == 0) ? 1 : -1;
			for (int j = 0; j != a; j++) if (j != k) {
				int iexp = r[j] + iz[j];
				if ((iexp % 2 ==1) && (lams[j] - lams[k] < 0)) sign = -sign;
					oneterm += gsl_sf_lnchoose(iexp - 1, iz[j]) - log(fabs(lams[j] - lams[k]))* iexp; //+ r[j]*loglams[j]

			}
			// oneterm += loglams[k]*r[k];
			if (sign > 0) {
				if (fplus) { hypoplus = oneterm; fplus = false; }
				else hypoplus = logsum(hypoplus, oneterm);
			}
			if (sign < 0) {
				if (fminus) { hypominus = oneterm; fminus = false; }
				else hypominus = logsum(hypominus, oneterm);
			}
			/*	for (int i = 0; i != a; i++) std::cout << std::setw(4) << iz[i];
						std::cout << std::endl;
						char x; std::cin >> x;
			 */
			// for (int i = 0; i != a; i++) Rprintf("%4d", iz[i]); Rprintf("\n");
		}
		temp = step(k, a, iz, l - 1);
	}

	// if ((l - 1) % 2 == 1) {
	// 	double temp = hypoplus;
	// 	hypoplus = hypominus;
	// 	hypominus = temp;
	// }

	free(iz);
//	hypoplus += gsl_sf_lnfact(l - 1);
//	hypominus += gsl_sf_lnfact(l - 1);
}

void loggammagaussian(int n, double lam, double mu, double sd, double t, double &hplus, double &hminus)
{
	hplus = hminus = GSL_NEGINF;
	if (n == 0) {
		hplus = (logexgaussian(lam, mu, sd, t)) - log(lam);
	}
	if (n == 1) {
		double x1 = mu + lam * sd * sd;
		double x2 = t - mu;
		double x3 = 0.5 * x1 * x1 / sd / sd;
		double x4 = 0.5 * (x1 - t) * (x1 - t) / sd / sd;
		double x6 = 0.5 * x2 * x2 / sd / sd;


		int sign_temp2 = (x1 - t) > 0 ? 1 : -1;

		double temp = x1 / sd, b, c;
		if (temp > 0) {
			c = lnnorm(-temp);
			b = lnnorm(-temp + t / sd);
		}
		else {
			b = lnnorm(temp);
			c = lnnorm(temp - t / sd);
		}
		double d = log(fabs(x1 - t));

		b = b + d + 0.5*M_LNPI + M_LN2 + x4 - x6;
		c = c + d + 0.5*M_LNPI + M_LN2 + x4 - x6;

		double lnsd = log(sd);

		double plus, minus;

		plus = 0.5*M_LN2 + lnsd - x6;
		minus = 0.5*M_LN2 + lnsd + x4 - x3 - x6;


		if (sign_temp2 == 1) {
			plus = logsum(plus, c);
			minus = logsum(minus, b);
		}
		else {
			plus = logsum(plus, b);
			minus = logsum(minus, c);
		}
		hplus = plus - M_LN2 - 0.5 * M_LNPI;
		hminus = minus - M_LN2 - 0.5 * M_LNPI;
	}

	if (n == 2) {
		double x1 = mu + lam * sd * sd;
		double x2 = t - mu;
		double x3 = 0.5 * x1 * x1 / sd / sd;
		double x4 = 0.5 * (x1 - t) * (x1 - t) / sd / sd;
		double x6 = 0.5 * x2 * x2 / sd / sd;

		int sign1 = (x1 - t) > 0 ? 1 : -1;
		int sign2 = (x1 - 2 * t) > 0 ? 1 : -1;

		double lnsd = log(sd);

		double plus, minus;
		plus = minus = GSL_NEGINF;
		double d = log(fabs(x1 - t));

		if (sign1==1) minus = -x6 + lnsd + d + 0.5* M_LN2;
		else plus = -x6 + lnsd + d + 0.5 * M_LN2;

		if (sign2==1) plus = logsum(plus,x4 - x3 - x6 + 0.5 * M_LN2 + lnsd + log(fabs(x1 - 2 * t)));
		else minus = logsum(minus, x4 - x3 - x6 + 0.5 * M_LN2 + lnsd + log(fabs(x1 - 2 * t)));

		double b = x4 - x6 + 0.5 * M_LNPI + log(sd * sd + (x1 - t) * (x1 - t)) + M_LN2;
		double temp = x1 / sd;
		if (temp > 0) {
			plus = logsum(plus, b + lnnorm(-temp + t / sd));
			minus = logsum(minus, b + lnnorm(-temp));
		}
		else
		{
			plus = logsum(plus, b + lnnorm(temp));
			minus = logsum(minus, b + lnnorm(temp-t/sd));
		}
		hplus = plus - 2*M_LN2 - 0.5 * M_LNPI;
		hminus = minus - 2*M_LN2 - 0.5 * M_LNPI;
	}

	if (n == 3) {
		double xx = 0.5 / sd / sd;
		double x1 = mu + lam * sd * sd;

		double x2 = (x1 * x1) * xx;
		// double x3 = (t * t) * xx;
		double x4 = t * lam;
		double x5 = (mu * mu) * xx;
		// double x6 = (2.0 * t * mu) * xx;

		double d1 = 3.0*t*t + 2.0*sd*sd - 3.0*t*x1 + x1*x1;
		double ln_abs_d1 = log(fabs(d1));

		double d2 = t-x1;
		double ln_abs_d2 = log(fabs(d2));

		int sign1 = d1 > 0 ? 1 : -1;
		int sign2 = d2 > 0 ? 1 : -1;

		double lnsd = log(sd);

		double plus, minus;
		plus = minus = GSL_NEGINF;

		if (sign1==1) plus = -x4-x5 + 0.5*M_LN2 + lnsd + ln_abs_d1;
		else minus = -x4-x5 + 0.5*M_LN2 + lnsd + ln_abs_d1;

		minus = logsum( minus, -(mu-t)*(mu-t)*xx + 0.5*M_LN2 + lnsd + log(2*sd*sd+(x1-t)*(x1-t)) );

		double b = x2-x4-x5 + 0.5*M_LNPI + ln_abs_d2 + log(3.0*sd*sd+(x1-t)*(x1-t)) + M_LN2;

		double temp = x1 / sd;
		if (temp > 0)
		{
			if (sign2==1) {
				plus = logsum(plus, b + lnnorm(-temp));
				minus = logsum(minus, b + lnnorm(-temp + t / sd));
			}
			else {
				plus = logsum(plus, b + lnnorm(-temp + t / sd));
				minus = logsum(minus, b + lnnorm(-temp));
			}

		}
		else
		{
			if (sign2==1) {
				plus = logsum(plus, b + lnnorm(temp - t / sd));
				minus = logsum(minus, b + lnnorm(temp));
			}
			else {
				plus = logsum(plus, b + lnnorm(temp));
				minus = logsum(minus, b + lnnorm(temp - t / sd));
			}
		}
// Achtung in Formel Mathematica f�hrendes "Minus" vergessen und hier nun durch Vertauschen von plus und minus ber�cksichtigt!
		hplus = minus - 2 * M_LN2 - log(3.0) - 0.5 * M_LNPI;
		hminus = plus - 2 * M_LN2 - log(3.0) - 0.5 * M_LNPI;
	}

		if (n == 4) {
			double xx = 0.5 / sd / sd;
			double x1 = mu + lam * sd * sd;

			double x2 = (x1 * x1) * xx;
			double x3 = (t * t) * xx;
			double x4 = t * lam;
			double x5 = (mu * mu) * xx;
			double x6 = (2.0 * t * mu) * xx;
			double sd2 = gsl_pow_2(sd), sd4 = gsl_pow_2(sd2);

			double d1 = -4.0 * gsl_pow_3(t) + 6.0 * gsl_pow_2(t) * x1 - 4.0 * t * (2.0 * sd2 + gsl_pow_2(x1)) + x1 * (5.0 * sd2 + gsl_pow_2(x1));
			double ln_abs_d1 = log(fabs(d1));

			double d2 = x1 - t;
			double ln_abs_d2 = log(fabs(d2));

			double t_mu = t - mu;
			double lam2 = gsl_pow_2(lam);
			double d3 = gsl_pow_4(t_mu) - 2.0 * gsl_pow_2(t_mu) * (-3.0 + 2.0 * lam * t_mu) * sd2 +
				3.0 * (1.0 + 2.0 * lam * t_mu * (-2.0 + lam * t_mu)) * sd4 +
				2.0 * lam2 * (3.0 - 2.0 * lam * t_mu) * sd4 * sd2 + gsl_pow_2(lam2 * sd4);
			double ln_abs_d3 = log(fabs(d3));

			int sign1 = d1 > 0 ? 1 : -1;
			int sign2 = d2 > 0 ? 1 : -1;
			int sign3 = d3 > 0 ? 1 : -1;

			double lnsd = log(sd);

			double plus, minus;
			plus = minus = GSL_NEGINF;

			if (sign1 == 1) plus = -x4 - x5 + 0.5 * M_LN2 + lnsd + ln_abs_d1;
			else minus = -x4 - x5 + 0.5 * M_LN2 + lnsd + ln_abs_d1;

			if (sign2 == 1) minus = logsum(minus, -x3 - x5 + x6 + 0.5 * M_LN2 + lnsd + ln_abs_d2 + log(5.0 * sd2 + (x1 - t) * (x1 - t)));
			else plus = logsum(plus, -x3 - x5 + x6 + 0.5 * M_LN2 + lnsd + ln_abs_d2 + log(5.0 * sd2 + (x1 - t) * (x1 - t)));

			double b = x2 - x4 - x5 + 0.5 * M_LNPI + ln_abs_d3 + M_LN2;

			double temp = x1 / sd;
			if (temp > 0)
			{
				if (sign3==1) {
					plus = logsum(plus, b + lnnorm(-temp + t / sd));
					minus = logsum(minus, b + lnnorm(-temp));
				}
				else {
					plus = logsum(plus, b + lnnorm(-temp));
					minus = logsum(minus, b + lnnorm(-temp + t / sd));
				}

			}
			else
			{
				if (sign3==1) {
					plus = logsum(plus, b + lnnorm(temp));
					minus = logsum(minus, b + lnnorm(temp - t / sd));
				}
				else {
					plus = logsum(plus, b + lnnorm(temp - t / sd));
					minus = logsum(minus, b + lnnorm(temp));
				}
			}
			hplus = plus - 4*M_LN2-log(3.0) - 0.5 * M_LNPI;
			hminus = minus - 4*M_LN2-log(3.0) - 0.5 * M_LNPI;
		}

}

/*
double logf_tij(int a, vector<int> r, double *lams, double *loglams, double mu, double sd, double t) {
//	double b = 0.0;
//	for (int i = 0; i != a; i++) b += log(lams[i])* r[i];
	double temp = 0.0;
	double hypoplus=GSL_NEGINF, hypominus=GSL_NEGINF;

	for (int k = 0; k != a; k++) for (int l = 0; l != r[k]; l++) {
		double plus1, minus1;
		double plus2, minus2;
		logPhikl(k, a, r, lams,loglams, l + 1, plus1, minus1);
		loggammagaussian(r[k] - 1 - l, lams[k], mu, sd, t, plus2, minus2);
		hypoplus = logsum(hypoplus, plus1 + plus2);// -gsl_sf_lnfact(l));
		hypoplus = logsum(hypoplus, minus1 + minus2); // -gsl_sf_lnfact(l));
		hypominus = logsum(hypominus, plus1 + minus2);// -gsl_sf_lnfact(l));
		hypominus = logsum(hypominus, minus1 + plus2);// -gsl_sf_lnfact(l));
	}

	if (hypoplus <= hypominus) {
		printf("logf_tij - Problem\n");
		for (int i = 0; i != a; i++) printf("%5d%12g%12g%12g\n", r[i], lams[i], mu, sd);
		printf("%g %g %g\n", hypoplus, hypominus, t);
//		std::cout << "logf_tij - Problem";
//		for (int i = 0; i != a; i++) std::cout << setw(5) << r[i] << setw(12) << lams[i] << setw(12) << mu << setw(12) << sd << std::endl;
//		std::cout << hypoplus << " " << hypominus << " " << t << std::endl;
//		char x; std::cin >> x;

		temp = GSL_NEGINF;
	}
	else temp = logdiff(hypoplus, hypominus);
//	temp = temp + b;
	return temp;
}
*/

bool trouble_shooter(int &a, vector<int> &r, double* lams, double* loglams) {
	// case 1 sehr gro�es lambda wird weggelassen!
	bool found_problem = false;
	int j = -1;
	for (int i = 0; i != a; i++) if (lams[i] > 1000.0) { j = i; found_problem = true; }
	if (found_problem) { r[j] = 0; }
	if (!found_problem) {
		int j1, j2;
		for (int i = 0; i != a; i++) for (int j = i + 1; j != a; j++) if (fabs(lams[i] - lams[j]) < 0.1) {
			j1 = i; j2 = j; found_problem = true;
		}
		if (found_problem) {
			r[j1] += r[j2]; r[j2] = 0;
			lams[j1] = (lams[j1] + lams[j2]) / 2.0;
			loglams[j1] = log(lams[j1]);
		}
	}
	if (found_problem) {
		int j = -1;
		vector<int> newr;
		for (int i = 0; i != a; i++) {
			if (r[i]!=0) {
				j++;
				newr.push_back(r[i]);
				lams[j] = lams[i];
				loglams[j] = loglams[i];
			}

		}
		a--;
		r = newr;
	}
	return found_problem;
}

double logf_tij(int a, vector<int> r, double *lams, double *loglams, double mu, double sd, double t) {
	double b = 0.0;
	for (int i = 0; i != a; i++) b += log(lams[i])* r[i];
	NEW: double temp = 0.0;
	vector<double> hypoplus, hypominus;

	double plus1, minus1;
	double plus2, minus2;
	for (int k = 0; k != a; k++) for (int l = 0; l != r[k]; l++) {

		logPhikl(k, a, r, lams,loglams, l + 1, plus1, minus1);

		loggammagaussian(r[k] - 1 - l, lams[k], mu, sd, t, plus2, minus2);

		if ((plus1!=GSL_NEGINF) && (plus2!=GSL_NEGINF)) hypoplus.push_back(plus1 + plus2);// -gsl_sf_lnfact(l));
		if ((minus1 != GSL_NEGINF) && (minus2 != GSL_NEGINF)) hypoplus.push_back(minus1 + minus2); // -gsl_sf_lnfact(l));
		if ((plus1 != GSL_NEGINF) && (minus2 != GSL_NEGINF)) hypominus.push_back(plus1 + minus2);// -gsl_sf_lnfact(l));
		if ((minus1 != GSL_NEGINF) && (plus2 != GSL_NEGINF)) hypominus.push_back(minus1 + plus2);// -gsl_sf_lnfact(l));
	}

	std::sort(hypoplus.begin(), hypoplus.end());
	std::sort(hypominus.begin(), hypominus.end());
	double hplus = GSL_NEGINF, hminus = GSL_NEGINF;
	for (int i = 0; i != static_cast<int>(hypoplus.size()); i++) hplus=logsum(hplus, hypoplus[i]);
	for (int i = 0; i != static_cast<int>(hypominus.size()); i++) hminus=logsum(hminus, hypominus[i]);

	if (hplus <= hminus) {
		// printf("logf_tij - Problem\n");
		// for (int i = 0; i != a; i++) printf("%5d%12g%12g%12g\n", r[i], lams[i], mu, sd);
		// printf("%g %g %g\n", hplus, hminus, t);
		//		char x; std::cin >> x;
		//		for (int i = 0; i != hypoplus.size(); i++) std::cout << setw(20) << hypoplus[i]; std::cout << std::endl;
		//		for (int i = 0; i != hypominus.size(); i++) std::cout << setw(20) << (hypominus[i]); std::cout << std::endl;
		//		std::cout << hypoplus - hypominus << std::endl;

		temp = GSL_NEGINF;
		if (trouble_shooter(a, r, lams, loglams)) {
			if (r.size() == 1) {
				loggammagaussian(r[0] - 1, lams[0], mu, sd, t, plus2, minus2);
				if (plus2 <= minus2) temp = GSL_NEGINF; else temp = logdiff(plus2, minus2);
			}
			else goto NEW;
		}
		// if (fabs(lams[0] - lams[1]) < -0.01) {
		// 	a = 1; r[0] = (r[1] > r[0]) ? r[1] : r[0];
		// 	loggammagaussian(2, lams[0], mu, sd, t, hplus, hminus);
		// 	temp = logdiff(hplus, hminus);
		// 	std::cout << "hallo" << setw(20) << temp << std::endl;
		// }
	}
	else temp = logdiff(hplus, hminus);
	temp = temp + b;
	if (DEBUG) if (!std::isfinite(temp)) Rprintf("f_tij\n");
	return temp;
}


/*
void test() {
	int a = 5;
//nt l = 3;
//nt k = 3;
	int *iz = 0; iz = (int *)malloc(a * sizeof(int));
	vector<int> r = { 2,1,2,2,1 };
	double *lams = 0; lams = (double *)malloc(a * sizeof(double));


	lams[0] = 1.0;
	lams[1] = 2.0;
	lams[2] = 3.0;
	lams[3] = 4.0;
	lams[4] = 5.0;

	double mu = 0.3;
	double sd = 0.4;
	double t = 0.4;

	const int N_d = 1000;
	double dens[N_d];
	double integ = 0;
	for (int i = 0; i != N_d; i++) {
		dens[i] = f_tij(a, r, lams, mu, sd, (i + 1)*1.0 / 10.0)/exp(lnnorm(mu/sd));
		integ += dens[i] * 1.0 / 10.0;
	}
	Rprintf("%g\n", integ);


	for (int i = 0; i != 100; i++) {
		double t = i * 1.0 / 100;	Rprintf("%5g%20g\n", t, gammagaussian(0, lams[4], 0.3, 0.04, t));
	}
	Rprintf("%g\n", f_tij(a, r, lams, mu, sd, t));
		//char x; std::cin >> x;
	// init_step(k, a, iz, r, l);
	// for (int i = 0; i != a; i++) std::cout << std::setw(4) << iz[i];
	// std::cout << std::endl;
	// char x; std::cin >> x;
	// bool temp = true;
	// while (temp)
	// {
	// 	temp = step(k, a, iz, r, l);
	// 	if (temp) for (int i = 0; i != a; i++) std::cout << std::setw(4) << iz[i];
	// 	std::cout << std::endl;
	// 	char x; std::cin >> x;
	// }

	free(iz);
	free(lams);
}
*/


void extract_pfadinfo(int *pfad_index,  vector<pfadinfo> &path_info) {

	int* counts = 0; counts = (int*)malloc(2 * kernpar * sizeof(int));
  path_info.clear();

	for (int c = 0; c != kerncat; c++) {
			int tree = cat2tree[c];
			for (int j = 0; j != branch[c]; j++) {
					int pl = NDRIN(c, j);

					for (int ip = 0; ip != 2 * kernpar; ip++) counts[ip] = 0;

					for (int in = 0; in != pl; in++) {
							int n = DRIN(c, j, in);
							int ip = TREE_AND_NODE2PAR(tree, n);
							int pm = (AR(c, j, n) > 0) ? 1 : 0;
							if (comp[ip + kernpar * (1 + pm)]) {
									int ipp = free2kern[kern2free[ip + (1 + pm) * kernpar]];
									counts[ipp - kernpar ]++;
							}
					}
					//            bool temp = false;
					//            for (int ip = 0; (ip != 2*kernpar) && !(temp); ip++) temp = (counts[ip] > 1);

					//            if (temp) {
					pfadinfo one;
					one.a = 0;
					for (int ip = 0; ip != kernpar; ip++)
							for (int pm = 0; pm != 2; pm++) {
									if (counts[ip + pm * kernpar] > 0) {
											one.r.push_back(counts[ip + pm * kernpar]);
											one.pfad_par.push_back(ip);
											one.pm.push_back(pm);
											one.a++;
									}
							}
					PFAD_INDEX(c, j) = path_info.size();
					path_info.push_back(one);
					//            }
					//            else PFAD_INDEX(c, j) = -1;
			}
	}

	if(counts) free(counts);



	// int *counts = 0; counts = (int *)malloc(2*kernpar * sizeof(int));
	// path_info.clear();
	// for (int c = 0; c != kerncat; c++) {
	// 	int tree = cat2tree[c];
	// 	for (int j = 0; j != branch[c]; j++) {
	// 		int pl = NDRIN(c, j);
	//
	// 		for (int ip = 0; ip != 2*kernpar; ip++) counts[ip] = 0;
	//
	// 		for (int in = 0; in != pl; in++) {
	// 			int n = DRIN(c, j, in);
	// 			int ip = TREE_AND_NODE2PAR(tree, n);
	// 			int pm = (AR(c, j, n) > 0) ? 1 : 0;
	// 			if (comp[ip + kernpar*(1+pm)]) counts[ip+pm*kernpar]++;
	// 		}
	// 		// if (temp) {
	// 		pfadinfo one;
	//
	// 		one.a = 0;
	// 		for (int ip = 0; ip != kernpar; ip++)
	// 			for (int pm = 0; pm != 2; pm++) {
	// 				if (counts[ip + pm * kernpar] > 0) {
	// 					one.r.push_back(counts[ip + pm * kernpar]);
	// 					one.pfad_par.push_back(ip);
	// 					one.pm.push_back(pm);
	// 					one.a++;
	// 				}
	// 			}
	// 		PFAD_INDEX(c, j) = path_info.size();
	// 		path_info.push_back(one);
	// 		// }
	// 		// else PFAD_INDEX(c, j) = -1;
	// 	}
	// }
	// free(counts);
}
