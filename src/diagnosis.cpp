//implements summary statistics on posterior distribution
//goodness-of-fit tests
//specific contrast
// authors: Christoph Klauer and Raphael Hartmann

#include "rts.h"

// namespace rtsNS {

std::ofstream tests_out;

#define SAMPLE(I,IP) sample[(I)*(n_all_parameters+1)+IP]
#define TREE_AND_NODE2PAR(ITREE,R) tree_and_node2par[ITREE*nodemax+R]
#define DRIN(J,K,X) drin[J*zweig*nodemax+K*nodemax + X]
#define NDRIN(J,K) ndrin[J*zweig+K]
#define AR(I,J,R) ar[I*zweig*nodemax + J*nodemax + R]
#define LAMBDAS(T,IZ) lambdas[T*ilamfree+IZ] //PM=0 negativ; PM=1 positiv
#define RHOS(IG,IZ) rhos[(IG)*ilamfree+IZ]

void lies(int n_all_parameters, double *sample) {

	std::ifstream rein(RAUS);
	int is, in;
	rein >> is >> in;
	if (is != SAMPLE_SIZE) Rprintf("HM\n"); if (in != (n_all_parameters+1)) Rprintf("HO\n");
	for (int i = 0; i != is; i++)
		for (int j = 0; j != in; j++) rein >> SAMPLE(i, j);
	rein.close();
}

void hdi(int length, double * parameter, double p, double iv[2]) {
	int inc = int(p*length) + 1;
	int nci = length - inc;
	int imin = -1;
	double tempold = parameter[length - 1] - parameter[0];
	for (int i = 0; i != nci; i++) {
		double temp = parameter[i + inc] - parameter[i];
		if (temp < tempold) {
			tempold = temp;
			imin = i;
		}
	}
	iv[0] = parameter[imin];
	iv[1] = parameter[imin + inc];

}

#define BETA(T,IP) beta[T*ifree+IP]

void belege_beta(double *sample, int is, double *beta) {
	for (int t = 0; t != indi; t++) for (int iz = 0; iz != ifree; iz++) BETA(t, iz) = SAMPLE(is, iz + ifree * t2group[t]) +
		SAMPLE(is, (iz + ilamfree * igroup + ((ifree + ilamfree)*(ifree + ilamfree + 1)) / 2 + (t + igroup)*ifree));
}


void belege_lambdas_rhos(double *sample, int is, double *rhos, double *lambdas) {
	for (int iz = 0; iz != ilamfree * igroup; iz++) rhos[iz] = SAMPLE(is, ifree*igroup + iz);
	for (int t = 0; t != indi; t++) for (int iz = 0; iz != ilamfree; iz++)
		lambdas[t*ilamfree + iz] = SAMPLE(is, ifree*igroup + ilamfree * igroup + ((ifree + ilamfree)*(ifree + ilamfree + 1)) / 2 + indi * ifree + t * ilamfree + iz);
}

void belege_nur_lambdas(double *sample, int is, double *lambdas) {
	for (int t = 0; t != indi; t++)
		for (int iz = 0; iz != ilamfree; iz++)
			lambdas[t*ilamfree + iz] = SAMPLE(is, ifree*igroup + iz + t2group[t] * ilamfree)*
			SAMPLE(is, ifree*igroup + ilamfree * igroup + ((ifree + ilamfree)*(ifree + ilamfree + 1)) / 2 + indi * ifree + t * ilamfree + iz);
}

void quantiles(vector<trial> daten, int n_all_parameters, double *sample) {
	double qv[5];
	double *temp = 0; temp = (double *)malloc(SAMPLE_SIZE * sizeof(double));
	//std::streamsize prec = cout.precision(); std::cout << std::setprecision(4);
	if (save_diagnose) tests_out << std::setprecision(4);
	Rprintf("theta per group [median, 95 and 99%% HDI]\n"); if (save_diagnose) tests_out << "MUs per group" << std::endl;
	for (int ig = 0; ig != igroup; ig++)
		for (int ip = 0; ip != kernpar; ip++) if (comp[ip]) {
			for (int j = 0; j != SAMPLE_SIZE; j++) temp[j] = gsl_cdf_ugaussian_P(SAMPLE(j, kern2free[ip] + ig * ifree));
			gsl_sort(temp, 1, SAMPLE_SIZE);
			qv[2] = gsl_stats_median_from_sorted_data(temp, 1, SAMPLE_SIZE);
			double iv[2]; hdi(SAMPLE_SIZE, temp, 0.95, iv); qv[1] = iv[0]; qv[3] = iv[1]; hdi(SAMPLE_SIZE, temp, 0.99, iv); qv[0] = iv[0]; qv[4] = iv[1];
			Rprintf("%3d", ip + 1); for (int iq = 0; iq != 5; iq++) Rprintf("%12.4g", qv[iq]); Rprintf("\n");
			if (save_diagnose) { tests_out << setw(3) << ip + 1; for (int iq = 0; iq != 5; iq++) tests_out << setw(12) << qv[iq]; tests_out << std::endl; }
		}

	Rprintf("tau- in ms per group [median, 95 and 99%% HDI]\n"); if (save_diagnose) tests_out << "rho minus per group" << std::endl;
	for (int ig = 0; ig != igroup; ig++)
		for (int ip = 0; ip != kernpar; ip++) if (comp[kernpar + ip]) {
			for (int j = 0; j != SAMPLE_SIZE; j++) temp[j] = 1000.0 / SAMPLE(j, igroup*ifree + ig * ilamfree + kern2free[kernpar + ip] - ifree);
			gsl_sort(temp, 1, SAMPLE_SIZE);
			qv[2] = gsl_stats_median_from_sorted_data(temp, 1, SAMPLE_SIZE);
			double iv[2]; hdi(SAMPLE_SIZE, temp, 0.95, iv); qv[1] = iv[0]; qv[3] = iv[1]; hdi(SAMPLE_SIZE, temp, 0.99, iv); qv[0] = iv[0]; qv[4] = iv[1];
			Rprintf("%3d", ip + 1); for (int iq = 0; iq != 5; iq++) Rprintf("%12.4g", qv[iq]); Rprintf("\n");
			if (save_diagnose) { tests_out << setw(3) << ip + 1; for (int iq = 0; iq != 5; iq++) tests_out << setw(12) << qv[iq]; tests_out << std::endl; }
		}
	Rprintf("tau+ in ms per group [median, 95 and 99%% HDI]\n"); if (save_diagnose) tests_out << "rho plus per group" << std::endl;
	for (int ig = 0; ig != igroup; ig++)
		for (int ip = 0; ip != kernpar; ip++) if (comp[ip + 2 * kernpar]) {
			for (int j = 0; j != SAMPLE_SIZE; j++) temp[j] = 1000.0 / SAMPLE(j, igroup*ifree + ig * ilamfree + kern2free[2 * kernpar + ip] - ifree);
			gsl_sort(temp, 1, SAMPLE_SIZE);
			qv[2] = gsl_stats_median_from_sorted_data(temp, 1, SAMPLE_SIZE);
			double iv[2]; hdi(SAMPLE_SIZE, temp, 0.95, iv); qv[1] = iv[0]; qv[3] = iv[1]; hdi(SAMPLE_SIZE, temp, 0.99, iv); qv[0] = iv[0]; qv[4] = iv[1];
			Rprintf("%3d", ip + 1); for (int iq = 0; iq != 5; iq++) Rprintf("%12.4g", qv[iq]); Rprintf("\n");
			if (save_diagnose) { tests_out << setw(3) << ip + 1; for (int iq = 0; iq != 5; iq++) tests_out << setw(12) << qv[iq]; tests_out << std::endl; }
		}

	Rprintf("SD and CORR of process params [median, 95 and 99%% HDI]\n"); if (save_diagnose) tests_out << "SIG" << std::endl;
	int iz = igroup * (ifree + ilamfree) - 1;
	for (int ix = 0; ix != ifree + ilamfree; ix++)
		for (int jz = ix; jz != ifree + ilamfree; jz++) {
			iz++;
			for (int j = 0; j != SAMPLE_SIZE; j++) temp[j] = SAMPLE(j, iz);
			gsl_sort(temp, 1, SAMPLE_SIZE);
			qv[2] = gsl_stats_median_from_sorted_data(temp, 1, SAMPLE_SIZE);
			double iv[2]; hdi(SAMPLE_SIZE, temp, 0.95, iv); qv[1] = iv[0]; qv[3] = iv[1]; hdi(SAMPLE_SIZE, temp, 0.99, iv); qv[0] = iv[0]; qv[4] = iv[1];
			Rprintf("%3d%3d", free2kern[ix]+1, free2kern[jz] + 1); for (int iq = 0; iq != 5; iq++) Rprintf("%12.4g", qv[iq]); Rprintf("\n");
			if (save_diagnose) { tests_out << setw(3) << free2kern[ix] + 1 << setw(3) << free2kern[jz] + 1; for (int iq = 0; iq != 5; iq++) tests_out << setw(12) << qv[iq]; tests_out << std::endl; }
		}

	//Rprintf("Restpars: Mean per group, Residual Variance in motor/encoding, Due to individual differences in motor/encoding\n");
	Rprintf("mu_gamma in ms per group [median, 95 and 99%% HDI]\n");
	if (save_diagnose) tests_out << "Restpars: Mean per group, Residual Variance in motor/encoding, Due to individual differences in motor/encoding" << std::endl;
	iz = ifree * igroup + ilamfree * igroup + ((ifree + ilamfree)*(ifree + ilamfree + 1)) / 2 + indi * ifree + indi * ilamfree;

	for (int ir = 0; ir != 1 + igroup * respno; ir++) {
		if (ir != igroup*respno) for (int j = 0; j != SAMPLE_SIZE; j++) temp[j] = 1000.0*SAMPLE(j, (iz + ir));
		if (ir == igroup*respno) for (int j = 0; j != SAMPLE_SIZE; j++) temp[j] = SAMPLE(j, (iz + ir));
		gsl_sort(temp, 1, SAMPLE_SIZE);
		qv[2] = gsl_stats_median_from_sorted_data(temp, 1, SAMPLE_SIZE);
		double iv[2]; hdi(SAMPLE_SIZE, temp, 0.95, iv); qv[1] = iv[0]; qv[3] = iv[1]; hdi(SAMPLE_SIZE, temp, 0.99, iv); qv[0] = iv[0]; qv[4] = iv[1];
		if (ir == igroup*respno) {
			Rprintf("omega^2 [median, 95 and 99%% HDI]\n");
		}
		for (int iq = 0; iq != 5; iq++) Rprintf("%12.4g", qv[iq]); Rprintf("\n");
		if (save_diagnose) { for (int iq = 0; iq != 5; iq++) tests_out << setw(12) << qv[iq]; tests_out << std::endl; }
	}

	Rprintf("SD and CORR of motor/encoding params [median, 95 and 99%% HDI]\n"); if (save_diagnose) tests_out << "RSIG" << std::endl;
	iz = n_all_parameters - restparsno + igroup * respno + 1 - 1;
	for (int ip = 0; ip != respno; ip++)
		for (int jp = ip; jp != respno; jp++) {
			iz++;
			for (int j = 0; j != SAMPLE_SIZE; j++) temp[j] = SAMPLE(j, iz);
			gsl_sort(temp, 1, SAMPLE_SIZE);
			qv[2] = gsl_stats_median_from_sorted_data(temp, 1, SAMPLE_SIZE);
			double iv[2]; hdi(SAMPLE_SIZE, temp, 0.95, iv); qv[1] = iv[0]; qv[3] = iv[1]; hdi(SAMPLE_SIZE, temp, 0.99, iv); qv[0] = iv[0]; qv[4] = iv[1];
			Rprintf("%3d%3d", ip + 1, jp + 1); for (int iq = 0; iq != 5; iq++) Rprintf("%12.4g", qv[iq]); Rprintf("\n");
			if (save_diagnose) { tests_out << setw(3) << ip + 1 << setw(3) << jp + 1; for (int iq = 0; iq != 5; iq++) tests_out << setw(12) << qv[iq]; tests_out << std::endl; }
		}

	// Daten zum Vergleich:
	double s = 0.0; int no_trials = int(daten.size());
	double *u = 0; u = (double *)malloc(indi * sizeof(double));
	int *nj = 0; nj = (int *)malloc(indi * sizeof(int));
	for (int t = 0; t != indi; t++) { u[t] = 0.0; nj[t] = 0; }

	for (int i = 0; i != no_trials; i++) { u[daten[i].person] += daten[i].rt / 1000.0; nj[daten[i].person]++; }
	for (int t = 0; t != indi; t++) { u[t] /= nj[t]; }

	for (int i = 0; i != no_trials; i++) s += gsl_pow_2(daten[i].rt / 1000.0 - u[daten[i].person]); s = s / (no_trials - 1);
	double grand = 0.0; for (int t = 0; t != indi; t++) grand += (u[t] * nj[t]) / no_trials;
	double salph = 0.0; for (int t = 0; t != indi; t++) salph += gsl_pow_2(u[t] - grand); salph = salph / (indi - 1);
	//Rprintf("Daten: Mean, Residual Variance, Due to Individual differences\n");
	Rprintf("RT: mean, variance, residual variance\n");
	Rprintf("%12.4g%12.4g%12.4g\n", grand, s, salph);
	if (save_diagnose) {
		tests_out << "Daten: Mean, Residual Variance, Due to Individual differences " << std::endl;
		tests_out << setw(12) << grand << setw(12) << s << setw(12) << salph << std::endl;
	}

	free(temp);
	free(u);
	free(nj);
}

void make_pij_for_one_trial_new(trial one, double *x_for_all, double *pij, double &pj) {
	// berechnet  p


#define X_FOR_ALL(T,IP) x_for_all[T*kernpar+IP]


	double d0ha;
	int j = one.category, t = one.person, itree = one.tree;


	// pj = 0.0;
	for (int k = 0; k != branch[j]; k++) {
//			 pij[k]=1.0;
		for (int ir = 0; ir != NDRIN(j, k); ir++) {
			int r = DRIN(j, k, ir); int ix = AR(j, k, r); int ip = TREE_AND_NODE2PAR(itree, r);
			d0ha = (ix > 0) ? lnnorm(X_FOR_ALL(t, ip)) : lnnorm(- X_FOR_ALL(t, ip));
			pij[k] = pij[k] + d0ha;
		}
		if (k == 0) pj = pij[k]; else pj = logsum(pj,pij[k]);
	}
	if (DEBUG) if (!(std::isfinite(pj))) {
		Rprintf("pj in diagnosis is not finite %d\n", branch[j]);
	}
}

#define PFAD_INDEX(J,K) pfad_index[J*zweig+K]

void make_tij_for_one_trial_new(trial one, double *rhos, double *lambdas, double *restpars, double *pij) {

	int t = one.person; /*int itree = one.tree;*/ int j = one.category; double rt = one.rt / 1000.0; int resp = cat2resp[j];

	double rmu = restpars[t2group[t] * respno + resp] + restpars[alphaoff + t * respno + resp]; double rsig = sqrt(restpars[t + sigalphaoff]);

	for (int k = 0; k != branch[j]; k++) {
		int pfadlength = NDRIN(j, k);
		double *lams = 0; lams = (double *)malloc(pfadlength * sizeof(double));
		int complength = 0;
/*		if (PFAD_INDEX(j, k) == -1){
			for (int ir = 0; ir != pfadlength; ir++) {
				int r = DRIN(j, k, ir); int ip = TREE_AND_NODE2PAR(itree, r); int pm = (AR(j, k, r) > 0) ? 1 : 0;
				if (comp[ip + (1 + pm)*kernpar]) {
					int iz = kern2free[(1 + pm)*kernpar + ip] - ifree; lams[complength] = LAMBDAS(t, iz)*RHOS(t2group[t], iz);
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
			int iz = kern2free[(1 + pm)*kernpar + ip] - ifree; lams[ir] = LAMBDAS(t,iz)*RHOS(t2group[t], iz);
		}
		// }
		if (complength == 0) { pij[k] = (-0.5*gsl_pow_2((rt - rmu) / rsig)) / sqr2pi / rsig; }
		// if ((complength == 1) && (PFAD_INDEX(j, k) == -1)) {
		// 	double lam = lams[0]; pij[k] = (logexgaussian(lam, rmu, rsig, rt));
		// 	if (DEBUG) {if (pij[k] < 0.0) Rprintf("pij[k] < 0; pf = 1 %.4g\n", pij[k]);}
		// }
		if ((complength == 1) /*&& (PFAD_INDEX(j, k) > -1)*/) {
			// int ipfad = PFAD_INDEX(j, k);
			// pfadinfo akt_pfad = path_info[ipfad];
			// if (complength != akt_pfad.a) Rprintf("complength != a");
			double lam = lams[0];
			double hplus, hminus;
			loggammagaussian(akt_pfad.r[0] - 1, lam, rmu, rsig, rt, hplus, hminus);
			double temp = logdiff(hplus, hminus) + log(lam) * akt_pfad.r[0];
			if (temp == GSL_NEGINF) {
				pij[k] = temp;
			}
			else pij[k] = (temp);
		}
		/*if ((complength >= 2) && (PFAD_INDEX(j, k) == -1)) {
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

			pij[k] = logdiff(hypo_plus, hypo_minus);
			if (hypo_plus<hypo_minus) {
				//std::cout << "pij[k] < 0; pf =2"<< pij[k] << " " << hypo_minus  << hypo_plus << " " << " " << rsig <<std::endl;
				for (int ir = 0; ir != complength; ir++) Rprintf("%15.4g", lams[ir]); Rprintf("\n");
				pij[k] = GSL_NEGINF;
			}

			free(loglams);
		}*/
		if ((complength >= 2) /*&& (PFAD_INDEX(j, k) > -1)*/) {
			double *loglams = 0; loglams = (double *)malloc(complength * sizeof(double));
			for (int ir = 0; ir != complength; ir++) loglams[ir] = log(lams[ir]);
			// int ipfad = PFAD_INDEX(j, k);
			// pfadinfo akt_pfad = path_info[ipfad];
			double temp = logf_tij(akt_pfad.a, akt_pfad.r, lams, loglams, rmu, rsig, rt);
			if (temp ==GSL_NEGINF) {
				pij[k] = temp;

			}
			else pij[k] = temp;
			if (loglams) free(loglams);
		}

		free(lams);
		//	 if ((pfadlength>5) || (pfadlength==0)) cout<< "pfadlength" << std::endl;
	}
}


void dic(int N_all_parameters, vector <trial> daten, double *beta, double *sample) {

	std::ofstream log_lik(LOGLIK);
	if (log_lik_flag) log_lik << std::setprecision(12);
	double dbar = 0.0, pd = 0.0, pv = 0.0;

	double *x_for_all = 0; x_for_all = (double *)malloc(indi*kernpar * sizeof(double));
	double *xbar = 0; xbar = (double *)malloc(n_all_parameters * sizeof(double));
	double *pij = 0; pij = (double *)malloc(zweig * sizeof(double));
	double *lambdas = 0; lambdas = (double *)malloc(ilamfree*indi * sizeof(double));
	double *rhos = 0; rhos = (double *)malloc(ilamfree*igroup * sizeof(double));
	double *restpars = 0; restpars = (double *)malloc(restparsno * sizeof(double));
	//		double *temp=0; if (!(temp=NAG_ALLOC(SAMPLE_SIZE*sizeof(double)))){printf("Allocation failure\n");exit_status = -1;}

	int trialno = int(daten.size());
	for (int i = 0; i != n_all_parameters; i++) xbar[i] = 0.0;
	for (int is = 0; is != SAMPLE_SIZE; is++) {

		for (int i = 0; i != n_all_parameters; i++) xbar[i] += SAMPLE(is, i) / (SAMPLE_SIZE);

		// belege x_for_all:
		belege_beta(sample, is, beta);
		for (int t = 0; t != indi; t++)  for (int ip = 0; ip != kernpar; ip++) x_for_all[ip + t * kernpar] = comp[ip] ? /*gsl_cdf_ugaussian_P*/(BETA(t, kern2free[ip])) : /*gsl_cdf_ugaussian_P*/(consts[ip]);
		// belege rhos, lambdas, restpars
		belege_lambdas_rhos(sample, is, rhos, lambdas);
		for (int ix = 0; ix != restparsno; ix++) restpars[ix] = SAMPLE(is, n_all_parameters - restparsno + ix);
		double persample = 0.0;
		for (int x = 0; x != trialno; x++) {
			trial one = daten[x];  int t = one.person; int r = cat2resp[one.category]; double rmu = restpars[t2group[t] * respno + r] + restpars[t*respno + r + alphaoff]; double rsig = sqrt(restpars[t + sigalphaoff]);
			double xsi = /*gsl_cdf_ugaussian_P*/lnnorm(rmu / rsig);
			make_tij_for_one_trial_new(one, rhos, lambdas, restpars, pij);
			double p;
			make_pij_for_one_trial_new(one, x_for_all, pij, p);
			p -= xsi;
			if (p==GSL_NEGINF)
				Rprintf("DIC loglik Problem\n");
			if (log_lik_flag) log_lik << setw(20) << p;
			persample += -2 * (p);
		}
		if (log_lik_flag) log_lik << std::endl;
		dbar += persample / (SAMPLE_SIZE);
		pv += gsl_pow_2(persample) / (SAMPLE_SIZE);
	}

	double xn = SAMPLE_SIZE * 1.0;
	pv = pv - gsl_pow_2(dbar); pv = xn / (xn - 1.0)*pv; pv = 0.5*pv;

	//belege x_for_all und so fort mit xbar
	for (int t = 0; t != indi; t++) for (int iz = 0; iz != ifree; iz++)  BETA(t, iz) = xbar[iz + t2group[t] * ifree] +
		xbar[(iz + ilamfree * igroup + ((ifree + ilamfree)*(ifree + ilamfree + 1)) / 2 + (t + igroup)*ifree)];
	for (int t = 0; t != indi; t++) 	for (int ip = 0; ip != kernpar; ip++) if (comp[ip])
		x_for_all[ip + t * kernpar] = /*gsl_cdf_ugaussian_P*/(BETA(t, kern2free[ip])); else
		x_for_all[ip + t * kernpar] = /*gsl_cdf_ugaussian_P*/(consts[ip]);
	// belege rhos, lambdas, restpars mit xbar
	for (int iz = 0; iz != ilamfree * igroup; iz++) rhos[iz] = xbar[ifree*igroup + iz];
	for (int t = 0; t != indi; t++) for (int iz = 0; iz != ilamfree; iz++)
		lambdas[t*ilamfree + iz] = xbar[ifree*igroup + ilamfree * igroup + ((ifree + ilamfree)*(ifree + ilamfree + 1)) / 2 + indi * ifree + t * ilamfree + iz];
	for (int ix = 0; ix != restparsno; ix++) restpars[ix] = xbar[n_all_parameters - restparsno + ix];
	//
	for (int x = 0; x != trialno; x++) {
		trial one = daten[x];  int t = one.person;  int r = cat2resp[one.category]; double rmu = restpars[t2group[t] * respno + r] + restpars[t*respno + r + alphaoff]; double rsig = sqrt(restpars[t + sigalphaoff]);
		double xsi = /*gsl_cdf_ugaussian_P*/lnnorm(rmu / rsig);
		make_tij_for_one_trial_new(one, rhos, lambdas, restpars, pij);
		double p;
		make_pij_for_one_trial_new(one, x_for_all, pij, p);
		p -= xsi;
		if (p==GSL_NEGINF) Rprintf("DIC loglik Problem in pd\n");
		pd += -2 * p;
	}


	pd = dbar - pd;

	//std::cout << std::setprecision(8);
	if (save_diagnose) tests_out << std::setprecision(8);

	Rprintf("DIC1, DIC2, pd, pv:\n");
	Rprintf("%15.8g%15.8g%15.8g%15.8g\n", pd + dbar, pv + dbar, pd, pv);
	if (save_diagnose) {
		tests_out << "DIC1, DIC2, pd, pv: " << std::endl;
		tests_out << setw(15) << pd + dbar << setw(15) << pv + dbar << setw(15) << pd << std::endl << setw(15) << pv << std::endl;
	}
	log_lik.close();
	free(pij);
	free(xbar);
	free(x_for_all);
	free(rhos);
	free(lambdas);
	free(restpars);
}



void test(double *t1, double *t2, string what) {
	double p = 0.0; double out[2] = { 0.0,0.0 };
	double r, count;
	for (int i = 0; i != SAMPLE_SIZE; i++) {
		r = 1.0 / (i + 1);
		count = (t1[i] < t2[i]) ? 1.0 : 0.0;
		out[0] += (t1[i] - out[0])*r; out[1] += (t2[i] - out[1])*r; p += (count - p)*r;
	}
	Rprintf("\n");
	Rprintf("%s\n", what.c_str());
	//std::cout << std::setprecision(4);
	Rprintf("%12.4g%12.4g%12.4g\n", out[0], out[1], p);
	if (save_diagnose) {
		tests_out << std::endl;
		tests_out << what << std::endl;
		tests_out << std::setprecision(4);
		tests_out << setw(12) << out[0] << setw(12) << out[1] << setw(12) << p << std::endl;
	}
	double iv[2];
	for (int is = 0; is != SAMPLE_SIZE; is++) t1[is] = (t1[is] - t2[is]);
	gsl_sort(t1, 1, SAMPLE_SIZE);
	hdi(SAMPLE_SIZE, t1, .95, iv);
	Rprintf("95%% HDI\n"); if (save_diagnose) tests_out << "95% HDI" << std::endl;
	for (int iq = 0; iq != 2; iq++) Rprintf("%12.4g", iv[iq]); Rprintf("\n");
	if (save_diagnose) {
		for (int iq = 0; iq != 2; iq++) tests_out << setw(12) << iv[iq]; tests_out << std::endl;
	}
}

void aggregate(int n_all_parameters, int kerntree, int *idaten, vector<trial> daten, int *nks, int *jks, int* tree2cat, double *beta, double *sample, gsl_rng *rst) {

#define IDATEN(T,I) idaten[T*kerncat + I]

#define NKS(T,IT) nks[T*kerntree+IT]
#define DRIN(J,K,X) drin[J*zweig*nodemax+K*nodemax + X]
#define NDRIN(J,K) ndrin[J*zweig+K]

//#define AR(I,J,K) ar[I*zweig*nodemax + J*nodemax + K]
//#define TREE_AND_NODE2PAR(I,J) tree_and_node2par[I*nodemax+J]
#define PIJ(I,J) pij[I*zweig+J]
#define TREE2CAT(IT,J) tree2cat[IT*kerncat+J]

	//	Nag_ModeRNG mode = Nag_GenerateWithoutReference;Nag_OrderType order = Nag_RowMajor;NagError fail;	INIT_FAIL(fail);

	double *t1 = 0; t1 = (double *)malloc(SAMPLE_SIZE * sizeof(double));
	double *t2 = 0; t2 = (double *)malloc(SAMPLE_SIZE * sizeof(double));

	int *obs = 0; obs = (int *)malloc(kerncat * sizeof(int));
	double *expe = 0; expe = (double *)malloc(kerncat * sizeof(double));
	int *rep = 0; rep = (int *)malloc(kerncat * sizeof(int));

	int *sobs = 0; sobs = (int *)malloc(kerncat*igroup * sizeof(int));
	double *sexp = 0; sexp = (double *)malloc(kerncat*igroup * sizeof(double));
	int *srep = 0; srep = (int *)malloc(kerncat*igroup * sizeof(int));

	double *tobs = 0; tobs = (double *)malloc(kerncat * sizeof(double));
	double *texp = 0; texp = (double *)malloc(kerncat * sizeof(double));
	double *trep = 0; trep = (double *)malloc(kerncat * sizeof(double));

	double *stobs = 0; stobs = (double *)malloc(kerncat*igroup * sizeof(double));
	double *stexp = 0; stexp = (double *)malloc(kerncat*igroup * sizeof(double));
	double *strep = 0; strep = (double *)malloc(kerncat*igroup * sizeof(double));

	double *pij = 0; pij = (double *)malloc(zweig*kerncat * sizeof(double));
	double *onepij = 0; onepij = (double *)malloc(zweig * sizeof(double));
	double *x = 0; x = (double *)malloc(kernpar * sizeof(double));
	double *lambdas = 0; lambdas = (double *)malloc(ilamfree*indi * sizeof(double));
	double *tdaten = 0; tdaten = (double *)malloc(indi*kerncat * sizeof(double));

	int *nobs = 0; nobs = (int *)malloc(kerncat*igroup * sizeof(int));
	int *nrep = 0; nrep = (int *)malloc(kerncat*igroup * sizeof(int));

	double *d = 0; d = (double *)malloc(kerncat * sizeof(double));
	double *x1 = 0; x1 = (double *)malloc(kerncat*igroup*SAMPLE_SIZE * sizeof(double));
	double *x2 = 0; x2 = (double *)malloc(kerncat*igroup*SAMPLE_SIZE * sizeof(double));
	unsigned int *drep = 0; drep = (unsigned int *)malloc(kerncat * sizeof(unsigned int));
	int *ng = 0; ng = (int *)malloc(igroup * sizeof(int));

#define X1(IS,IG,J) x1[IS*igroup*kerncat+ IG*kerncat+J]
#define X2(IS,IG,J) x2[IS*igroup*kerncat+IG*kerncat+J]
#define SOBS(IG,J) sobs[IG*kerncat+J]
#define SEXP(IG,J) sexp[IG*kerncat+J]
#define SREP(IG,J) srep[IG*kerncat+J]
#define STOBS(IG,J) stobs[IG*kerncat+J]
#define STEXP(IG,J) stexp[IG*kerncat+J]
#define STREP(IG,J) strep[IG*kerncat+J]
#define NOBS(IG,J) nobs[IG*kerncat+J]
#define NREP(IG,J) nrep[IG*kerncat+J]


	for (int is = 0; is != SAMPLE_SIZE; is++) {
		for (int j = 0; j != kerncat * igroup; j++) { sexp[j] = 0.0; sobs[j] = srep[j] = 0; x1[j + is * kerncat*igroup] = 0.0; x2[j + is * kerncat*igroup] = 0.0; }
		// compute exp pro Person und rep pro Person; aggregate, compute chi-square
		belege_beta(sample, is, beta);
		for (int t = 0; t != indi; t++) {
			for (int ip = 0; ip != kernpar; ip++) x[ip] = comp[ip] ? gsl_cdf_ugaussian_P(BETA(t, kern2free[ip])) : gsl_cdf_ugaussian_P(consts[ip]);
			make_pij_for_individual(x, pij, expe); for (int j = 0; j != kerncat; j++) drep[j] = 0;
			for (int it = 0; it != kerntree; it++) {
				for (int j = 0; j != jks[it]; j++) d[j] = expe[TREE2CAT(it, j)];
				gsl_ran_multinomial(rst, jks[it], NKS(t, it), d, drep);
				for (int j = 0; j != jks[it]; j++) rep[TREE2CAT(it, j)] = drep[j];
			}
			int ig = t2group[t];
			for (int j = 0; j != kerncat; j++) { expe[j] = NKS(t, cat2tree[j])*expe[j]; obs[j] = IDATEN(t, j); }
			for (int j = 0; j != kerncat; j++) { SEXP(ig, j) += expe[j]; SOBS(ig, j) += obs[j]; SREP(ig, j) += rep[j]; X1(is, ig, j) += (obs[j] * 1.0) / NKS(t, cat2tree[j]); X2(is, ig, j) += (rep[j] * 1.0) / NKS(t, cat2tree[j]); }
			//			for (int it=0;it!=kerntree;it++){int no=0,ne=0; for (int j=0;j!=jks[it];j++) {no+=obs[TREE2CAT(it,j)];ne+=rep[TREE2CAT(it,j)];}; if (no!=ne) std::cout<< "hello" << std::endl;}

		}
		t1[is] = 0.0; t2[is] = 0.0; for (int ig = 0; ig != igroup; ig++) for (int j = 0; j != kerncat; j++) { t1[is] += gsl_pow_2(SOBS(ig, j) - SEXP(ig, j)) / SEXP(ig, j); t2[is] += gsl_pow_2(SREP(ig, j) - SEXP(ig, j)) / SEXP(ig, j); }
	}
	test(t1, t2, "Posterior predictive check: frequencies");

	for (int ig = 0; ig != igroup; ig++) { ng[ig] = 0; for (int t = 0; t != indi; t++) ng[ig] += (t2group[t] == ig); }

	double qv[5];
	// int *correct = 0; correct = (int *)malloc(kerncat * sizeof(int));
	// for (int j = 0; j != kerncat; j++) correct[j] = ((j % 6 == 1) || (j % 6 == 2) || (j % 6 == 4)) ? 1 : 0;
	//std::cout << std::setprecision(4);
	for (int ig = 0; ig != igroup; ig++)
		for (int j = 0; j != kerncat; j++) {
			for (int is = 0; is != SAMPLE_SIZE; is++) t1[is] = X1(is, ig, j) / ng[ig];
			gsl_sort(t1, 1, SAMPLE_SIZE);
			qv[2] = gsl_stats_median_from_sorted_data(t1, 1, SAMPLE_SIZE);
			Rprintf("%3d", j); Rprintf("%12.4g", qv[2]);
			for (int is = 0; is != SAMPLE_SIZE; is++) t2[is] = X2(is, ig, j) / ng[ig];
			gsl_sort(t2, 1, SAMPLE_SIZE);
			qv[2] = gsl_stats_median_from_sorted_data(t2, 1, SAMPLE_SIZE);
			double iv[2]; hdi(SAMPLE_SIZE, t2, 0.95, iv); qv[1] = iv[0]; qv[3] = iv[1]; hdi(SAMPLE_SIZE, t2, 0.99, iv); qv[0] = iv[0]; qv[4] = iv[1];
			Rprintf("%12.4g%12.4g%12.4g\n", qv[1], qv[2], qv[3]);
		}


	// Das Ganze f�r die Zeiten

	for (int j = 0; j != kerncat * indi; j++) tdaten[j] = 0.0;
#define TDATEN(T,J) tdaten[T*kerncat+J]
	for (int it = 0; it != int(daten.size()); it++) {
		trial one = daten[it]; TDATEN(one.person, one.category) += one.rt / 1000.0;
	}
	for (int t = 0; t != indi; t++) for (int j = 0; j != kerncat; j++) if (IDATEN(t, j) > 0) TDATEN(t, j) /= IDATEN(t, j);



	for (int is = 0; is != SAMPLE_SIZE; is++) {
		for (int j = 0; j != kerncat * igroup; j++) { stobs[j] = stexp[j] = strep[j] = 0.0; nobs[j] = nrep[j] = 0; }
		// compute exp pro Person und rep pro Person; aggregate, compute chi-square
		belege_beta(sample, is, beta);
		belege_nur_lambdas(sample, is, lambdas);

		for (int j = 0; j != kerncat; j++) nobs[j] = nrep[j] = 0;

		for (int t = 0; t != indi; t++) {
			for (int ip = 0; ip != kernpar; ip++) x[ip] = comp[ip] ? gsl_cdf_ugaussian_P(BETA(t, kern2free[ip])) : gsl_cdf_ugaussian_P(consts[ip]);
			make_pij_for_individual(x, pij, expe);

			for (int j = 0; j != kerncat; j++) {
				int r = cat2resp[j];
				double mu = SAMPLE(is, (n_all_parameters - restparsno + t2group[t] * respno + r)) + SAMPLE(is, (n_all_parameters - restparsno + alphaoff + t * respno + r)), sig = sqrt(SAMPLE(is, (n_all_parameters - indi + t))); double z = -mu / sig;
				z = 1.0 / sqrt(2.0*M_PI)*exp(-0.5*gsl_pow_2(z)) / (1.0 - gsl_cdf_ugaussian_P(z));
				texp[j] = 0.0;
				for (int k = 0; k != branch[j]; k++)
					for (int xr = 0; xr != NDRIN(j, k); xr++) {
						int r = DRIN(j, k, xr); int ip = TREE_AND_NODE2PAR(cat2tree[j], r);  int pm = (AR(j, k, r) > 0) ? 1 : 0;
						if (comp[(1 + pm)*kernpar + ip]) {
							double lambda = LAMBDAS(t, (kern2free[(1 + pm)*kernpar + ip] - ifree)); texp[j] += PIJ(j, k)*1.0 / lambda;
						}
					}
				texp[j] += mu + sig * z;

				// compute exp und variance
			}
			for (int j = 0; j != kerncat; j++) drep[j] = 0;
			for (int it = 0; it != kerntree; it++) {
				for (int j = 0; j != jks[it]; j++) d[j] = expe[TREE2CAT(it, j)];
				gsl_ran_multinomial(rst, jks[it], NKS(t, it), d, drep);
				//				nag_rand_gen_multinomial(order, mode, 1,NKS(t,it),jks[it],d, r, lr,rst,drep,jks[it],&fail); if (fail.code != NE_NOERROR){printf("Error from nag_rand_gen_multinomial (g05tgc).\n%s\n",fail.message);exit_status = 1;}
				for (int j = 0; j != jks[it]; j++) rep[TREE2CAT(it, j)] = drep[j];
			}
			for (int j = 0; j != kerncat; j++) {
				int r = cat2resp[j];
				double mu = SAMPLE(is, (n_all_parameters - restparsno + t2group[t] * respno + r)) + SAMPLE(is, (n_all_parameters - restparsno + alphaoff + t * respno + r)), sig = sqrt(SAMPLE(is, (n_all_parameters - indi + t)));
				trep[j] = 0.0;
				tobs[j] = TDATEN(t, j);
				for (int ir = 0; ir != rep[j]; ir++) {
					double temp = 0.0;
					for (int k = 0; k != branch[j]; k++) onepij[k] = log(PIJ(j, k));
					int ipath = make_path_for_one_trial(branch[j], onepij, 0.0, rst);
					for (int xr = 0; xr != NDRIN(j, ipath); xr++) {
						int r = DRIN(j, ipath, xr); int ip = TREE_AND_NODE2PAR(cat2tree[j], r);  int pm = (AR(j, ipath, r) > 0) ? 1 : 0;
						if (comp[(1 + pm)*kernpar + ip]) {
							double lambda = LAMBDAS(t, (kern2free[(1 + pm)*kernpar + ip] - ifree));
							temp += oneexp(lambda, rst);
						}
					}
					temp += truncnorm(mu / sig, rst)*sig;
					trep[j] += temp;
				}
				if (rep[j] > 0) trep[j] /= rep[j]; else trep[j] = 0.0;
			}
			int ig = t2group[t];
			for (int j = 0; j != kerncat; j++) {
				STEXP(ig, j) += texp[j]; if (tobs[j] > 0) { STOBS(ig, j) += tobs[j]; NOBS(ig, j)++; } if (trep[j] > 0) { STREP(ig, j) += trep[j]; NREP(ig, j)++; }
			}
		}

		for (int ig = 0; ig != igroup; ig++) for (int j = 0; j != kerncat; j++) { STEXP(ig, j) /= ng[ig]; STOBS(ig, j) /= NOBS(ig, j); if (NREP(ig, j) > 0) STREP(ig, j) /= NREP(ig, j); X1(is, ig, j) = STOBS(ig, j); X2(is, ig, j) = STREP(ig, j); }
		t1[is] = 0.0; t2[is] = 0.0;

		for (int j = 0; j != kerncat * igroup; j++) {
			t1[is] += gsl_pow_2(stobs[j] - stexp[j]) / stexp[j]; t2[is] += gsl_pow_2(strep[j] - stexp[j]) / stexp[j];
		}
	}
	test(t1, t2, "Posterior predictive check: latencies");



	//std::cout << std::setprecision(4);
	for (int ig = 0; ig != igroup; ig++)
		for (int j = 0; j != kerncat; j++) {
			for (int is = 0; is != SAMPLE_SIZE; is++) t1[is] = X1(is, ig, j);
			gsl_sort(t1, 1, SAMPLE_SIZE);
			qv[2] = gsl_stats_median_from_sorted_data(t1, 1, SAMPLE_SIZE);
			Rprintf("%3d", j); Rprintf("%12.4g", qv[2]);
			for (int is = 0; is != SAMPLE_SIZE; is++) t2[is] = X2(is, ig, j);
			gsl_sort(t2, 1, SAMPLE_SIZE);
			qv[2] = gsl_stats_median_from_sorted_data(t2, 1, SAMPLE_SIZE);
			double iv[2]; hdi(SAMPLE_SIZE, t2, 0.95, iv); qv[1] = iv[0]; qv[3] = iv[1]; hdi(SAMPLE_SIZE, t2, 0.99, iv); qv[0] = iv[0]; qv[4] = iv[1];
			Rprintf("%12.4g%12.4g%12.4g\n", qv[1], qv[2], qv[3]);
		}
	// if (correct) free(correct);

	free(t1);
	free(t2);
	free(obs);
	free(expe);
	free(rep);
	free(sobs);
	free(sexp);
	free(srep);
	free(tobs);
	free(texp);
	free(trep);
	free(stobs);
	free(stexp);
	free(strep);
	free(pij);
	free(onepij);
	free(x);
	free(lambdas);
	free(tdaten);
	free(nobs);
	free(nrep);
	free(d);
	free(drep);
	free(x1);
	free(x2);
	free(ng);
}




void correlation(double *sample, int *index1, int *index2)
{
#define SIGMA(I,J) sigma[I*(ilamfree+ifree)+J]
	double qv[5];
	double *temp = 0; temp = (double *)malloc(SAMPLE_SIZE * sizeof(double));
	double *sigma = 0; sigma = (double *)malloc((ilamfree + ifree)*(ilamfree + ifree) * sizeof(double));
	for (int is = 0; is != SAMPLE_SIZE; is++) {
		int iz = (igroup)*(ifree + ilamfree) - 1;
		for (int ip = 0; ip != ilamfree + ifree; ip++)
			for (int jp = ip; jp != ilamfree + ifree; jp++) {
				iz++;
				SIGMA(ip, jp) = SAMPLE(is, iz); SIGMA(jp, ip) = SIGMA(ip, jp);
			}
		for (int ip = 0; ip != ilamfree + ifree; ip++) for (int jp = 0; jp != ilamfree + ifree; jp++) if (ip != jp)
			SIGMA(ip, jp) = sqrt(SIGMA(ip, ip))*SIGMA(ip, jp)*sqrt(SIGMA(jp, jp));
		double sig1 = 0.0, sig2 = 0.0, cov = 0.0;
		for (int ip = 0; ip != ilamfree + ifree; ip++) for (int jp = 0; jp != ilamfree + ifree; jp++) {
			sig1 += index1[ip] * index1[jp] > 0 ? SIGMA(ip, jp) : 0;
			sig2 += index2[ip] * index2[jp] > 0 ? SIGMA(ip, jp) : 0;
			cov += index1[ip] * index2[jp] > 0 ? SIGMA(ip, jp) : 0;
		}
		temp[is] = cov / sqrt(sig1*sig2);
	}
	gsl_sort(temp, 1, SAMPLE_SIZE);
	qv[2] = gsl_stats_median_from_sorted_data(temp, 1, SAMPLE_SIZE);
	double iv[2]; hdi(SAMPLE_SIZE, temp, 0.95, iv); qv[1] = iv[0]; qv[3] = iv[1]; hdi(SAMPLE_SIZE, temp, 0.99, iv); qv[0] = iv[0]; qv[4] = iv[1];
	Rprintf("Corr"); for (int iq = 0; iq != 5; iq++) Rprintf("%12.4g", qv[iq]); Rprintf("\n");
	if (save_diagnose) { tests_out << "Corr "; for (int iq = 0; iq != 5; iq++) tests_out << setw(12) << qv[iq]; tests_out << std::endl; }

	free(sigma);
	free(temp);
}

void groupwise(double *sample) {


	double *t1 = 0; t1 = (double *)malloc(SAMPLE_SIZE * sizeof(double));
	double *t2 = 0; t2 = (double *)malloc(SAMPLE_SIZE * sizeof(double));

	for (int ip = 0; ip != ifree; ip++)
	{
		for (int is = 0; is != SAMPLE_SIZE; is++) {
			//				par.push_back(gsl_cdf_ugaussian_P(SAMPLE(is,ip))-gsl_cdf_ugaussian_P(SAMPLE(is,ip+ifree)));
			t1[is] = gsl_cdf_ugaussian_P(SAMPLE(is, ip));
			t2[is] = gsl_cdf_ugaussian_P(SAMPLE(is, ip + ifree));
		}
		test(t1, t2, "group-tests mu");
	}

	for (int ip = 0; ip != ilamfree; ip++)
	{

		for (int is = 0; is != SAMPLE_SIZE; is++) {
			//				par.push_back(1000.0/SAMPLE(is,ip+ifree*igroup)-1000.0/SAMPLE(is,ip+ifree*igroup+2*ilamfree));
			t1[is] = 1000.0 / SAMPLE(is, ip + ifree * igroup);
			t2[is] = 1000.0 / SAMPLE(is, ip + ifree * igroup + ilamfree);
		}
		test(t1, t2, "group-tests pho");
	}
	int    iz = ifree * igroup + ilamfree * igroup + ((ifree + ilamfree)*(ifree + ilamfree + 1)) / 2 + indi * ifree + indi * ilamfree;
	for (int ip = 0; ip != respno; ip++)
	{

		for (int is = 0; is != SAMPLE_SIZE; is++) {
			//				par.push_back(SAMPLE(is,iz+ip)-SAMPLE(is,ip+iz+ respno));
			t1[is] = SAMPLE(is, iz + ip);
			t2[is] = SAMPLE(is, iz + ip + respno);
		}
		test(t1, t2, "group-tests residual");
	}
	iz = ifree * igroup + ilamfree * igroup + ((ifree + ilamfree)*(ifree + ilamfree + 1)) / 2 + indi * ifree + indi * ilamfree;
//--------
	// if (respno > 1) {
	// 	for (int ig = 0; ig != igroup; ig++)
	// 	{
	//
	// 		for (int is = 0; is != SAMPLE_SIZE; is++) {
	// 			//				par.push_back(SAMPLE(is,iz+ip)-SAMPLE(is,ip+iz+ respno));
	// 			t1[is] = SAMPLE(is, iz + 0);
	// 			t2[is] = SAMPLE(is, iz + 1 + ig * respno);
	// 		}
	// 		test(t1, t2, "within-group residuals");
	// 	}
	// }

	free(t1);
	free(t2);

}

void diagnosis(vector<trial> daten, int *idaten, int kerntree, gsl_rng *rst) {
	int *nks = 0; nks = (int *)malloc(indi*kerntree * sizeof(int));
	int *jks = 0; jks = (int *)malloc(kerntree * sizeof(int));
	int *tree2cat = 0; tree2cat = (int *)malloc(kerntree*kerncat * sizeof(int));
	double *beta = 0;	beta = (double *)malloc(indi*ifree * sizeof(double));

	n_all_parameters = ifree * igroup + ilamfree * igroup + ((ifree + ilamfree)*(ifree + ilamfree + 1)) / 2 + indi * ifree + indi * ilamfree + restparsno;

	double *sample = 0;
	sample = (double *)malloc(SAMPLE_SIZE*(n_all_parameters+1) * sizeof(double));
	lies(n_all_parameters, sample);
	if (save_diagnose) tests_out.open(diagn_tests);
	quantiles(daten, n_all_parameters, sample);
	// make nks

	for (int t = 0; t != indi; t++) for (int it = 0; it != kerntree; it++) NKS(t, it) = 0;
	for (int t = 0; t != indi; t++) for (int j = 0; j != kerncat; j++) NKS(t, cat2tree[j]) += IDATEN(t, j);
	for (int it = 0; it != kerntree; it++) jks[it] = 0;
	for (int j = 0; j != kerncat; j++) { TREE2CAT(cat2tree[j], jks[cat2tree[j]]) = j; jks[cat2tree[j]]++; }
	/*
	int *index1 = 0; index1 = (int *)malloc((ilamfree + ifree) * sizeof(int));
	int *index2 = 0; index2 = (int *)malloc((ilamfree + ifree) * sizeof(int));

	for (int i = 0; i != ilamfree + ifree; i++) index1[i] = index2[i] = 0;
	for (int i = 0; i != ilamfree / 2; i++) index1[ifree + i] = index2[ifree + ilamfree / 2 + i] = 1;
	correlation(sample, index1, index2);
	for (int i = 0; i != ilamfree + ifree; i++) index1[i] = index2[i] = 0;
	for (int i = 0; i != 3; i++) index1[i] = index2[ifree + i] = index2[ifree + ilamfree / 2 + i] = 1;
	correlation(sample, index1, index2);
	*/
	dic(n_all_parameters, daten, beta, sample);
	aggregate(n_all_parameters, kerntree, idaten, daten, nks, jks, tree2cat, beta, sample, rst);
	//    SSE_by_individual(n_all_parameters,kerntree, idaten, daten,nks,jks,tree2cat, beta,index,sample,rst);

	// if (igroup > 1) groupwise(sample);


	if (save_diagnose) tests_out.close();

	free(nks);
	// free(index1);
	// free(index2);
	free(jks);
	free(beta);

	free(tree2cat);
	free(sample);
}

// }
