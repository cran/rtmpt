#include "rts.h"

namespace drtmpt {

  // reshullfe hampar: group mu_theta, deviations_theta, group_mu_r, deviations_r, sig_t in response,
  //partial correlations theta, sig_theta, partial correlatons mu_r, gamma_r, log(omega^2)


  //from model parameters to partial correlations
  void from_y_to_z(int flag, gsl_vector* hampar, std::vector<double> &z) {
  	int k; int offset = nhamil;
  	z.clear();
  	if (flag == 0) {
  		k = icompg;
  	}
  	else {
  		k = respno;
  		offset += (icompg * (icompg - 1)) / 2 + icompg;
  	}
  	int kc2 = (k * (k - 1)) / 2;
  	for (int i = 0; i != kc2; i++)
  		z.push_back(tanh(gsl_vector_get(hampar,offset + i)));
  }

  //from partial correlations to Cholesky factor of correlation matrix
  void from_z_to_w(int flag, const std::vector<double> &z, gsl_matrix* w) {
  	int k;
  	if (flag == 0) {
  		k = icompg;
  	}
  	else {
  		k = respno;
  	}
  	//lower triangular
  	gsl_matrix_set(w, 0, 0, 1.0);
  	int jj = 0;
  	for (int i = 1; i < k; ++i) {
  		gsl_matrix_set(w, i, 0, z[jj++]);
  		double sum_sqs = gsl_pow_2(gsl_matrix_get(w, i, 0));
  		for (int j = 1; j < i; ++j) {
  			gsl_matrix_set(w, i, j, z[jj++] * sqrt(1.0 - sum_sqs));
  			sum_sqs += gsl_pow_2(gsl_matrix_get(w, i, j));
  		}
  		gsl_matrix_set(w, i, i, sqrt(1.0 - sum_sqs));
  	}
  }

  //from correlation matrix and standard deviations to variance-covariance matrix (flag = 0: diffusion parameters; flag = 0: motor-time parameters)
  void from_w_to_sig_sigi(int flag, gsl_vector* hampar, gsl_matrix* w, double* sig) {
  	int k; int offset = nhamil + (icompg * (icompg - 1)) / 2;
  	// extrahiere standardabweichungen
  	if (flag == 0) {
  		k = icompg;
  	}
  	else {
  		k = respno;
  		offset += icompg + (respno * (respno - 1)) / 2;
  	}
  	std::vector<double> wsigs;
  	for (int i = 0; i != k; i++) wsigs.push_back(exp(gsl_vector_get(hampar, offset + i)));
  	for (int i = 0; i != k; i++)
  		for (int j = 0; j <= i; j++)
  			gsl_matrix_set(w, i, j, gsl_matrix_get(w, i, j) * wsigs[i]);

  	gsl_matrix_view t1 = gsl_matrix_view_array(sig, k, k);
  	gsl_matrix_transpose_memcpy(&t1.matrix, w);
  	gsl_blas_dtrmm(CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, 1.0, w, &t1.matrix);
  }

  //putting the above together: from model parameters to variance-covariance matrices (flag = 0: diffusion parameters; flag = 0: motor-time parameters)
  void make_sigs(int flag, gsl_vector* hampar, double* sig) {
  	int k;
  	if (flag == 0) {
  		k = icompg;
  	}
  	else {
  		k = respno;
  	}
  	std::vector<double> z;
  	from_y_to_z(flag, hampar, z);
  	gsl_matrix* w = gsl_matrix_calloc(k, k);
  	from_z_to_w(flag, z, w);
  	from_w_to_sig_sigi(flag, hampar, w, sig);
  	gsl_matrix_free(w);
  }

  //make diffusion-model parameters from model parameters
  void make_tavwtlams(int flag, gsl_vector* hampar, std::vector<double>& z, gsl_matrix* w, double* tpars) {
  	int k, ioff; int offset = nhamil + (icompg * (icompg - 1)) / 2;
  	int offset2 = nhamil;
  	if (flag == 0)
  	{
  		k = icompg;
  		ioff = iavwoff;
  	}
  	else {
  		k = respno;
  		ioff = ilamoff;
  		offset += icompg + (respno * (respno - 1)) / 2;
  		offset2 += (icompg * (icompg - 1)) / 2 + icompg;
  	}

  	from_y_to_z(flag, hampar, z);


  	from_z_to_w(flag, z, w);

  	gsl_vector* wsigs = gsl_vector_alloc(k);
  	for (int i = 0; i != k; i++) gsl_vector_set(wsigs, i, exp(gsl_vector_get(hampar, offset + i)));

  	gsl_vector_view thetv = gsl_vector_subvector(hampar, ioff, indi * k);
  	gsl_matrix_view thetm = gsl_matrix_view_vector(&thetv.vector, indi, k);
  	gsl_matrix* thetmtr = gsl_matrix_alloc(k, indi);
  	gsl_matrix_transpose_memcpy(thetmtr, &thetm.matrix);
  	gsl_matrix* sd = gsl_matrix_calloc(k, k);
  	gsl_vector_view sddiag = gsl_matrix_diagonal(sd);
  	gsl_vector_memcpy(&sddiag.vector, wsigs);

  	gsl_blas_dtrmm(CblasRight, CblasLower, CblasNoTrans, CblasNonUnit, 1.0, w, sd);

  	gsl_blas_dtrmm(CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, 1.0, sd, thetmtr);

  	// pro Gruppe, so dass nur einmal pro Gruppe mu extrahiert wird?
  	for (int t = 0; t != indi; t++) {
  		gsl_vector_view mu = gsl_vector_subvector(hampar, (flag == 0) ? t2group[t] * k : irmuoff + t2group[t] * k, k);
  		gsl_vector_view thetind = gsl_matrix_column(thetmtr, t);
  		gsl_vector_add(&thetind.vector, &mu.vector);
  		if (flag == 0) {
  			int jj = 0;
  			for (int type = 0; type != 3; type++)
  				for (int ip = 0; ip != ifree[type]; ip++) {
  					if (dCOMP(type, ip))
  						tpars[t * 3 * ifreemax + type * ifreemax + ip] =
  							logit(avwtrans[type], gsl_vector_get(&thetind.vector, jj++));
  				}
  		}
  		else {
  			gsl_vector_view tparsv = gsl_vector_view_array(tpars, indi * k);
  			gsl_vector_view tparsind = gsl_vector_subvector(&tparsv.vector, t * k, k);
  			gsl_vector_memcpy(&tparsind.vector, &thetind.vector);
  		}
  	}
  	gsl_matrix_free(sd);
  	gsl_vector_free(wsigs);
  	gsl_matrix_free(thetmtr);
  }

  //from variance-covariance matrix to Cholesky factor correlation matrix
  void from_sig_to_w(int flag, gsl_vector* hampar, gsl_matrix* w, double* sig) {
  	int k; int offset = nhamil + (icompg * (icompg - 1)) / 2; int ioff = iavwoff;
  	// extrahiere standardabweichungen
  	if (flag == 0) {
  		k = icompg;
  	}
  	else {
  		k = respno;
  		offset += icompg + (respno * (respno - 1)) / 2;
  		ioff = ilamoff;
  	}
  	// make_sig
  	std::vector<double> wsigs;
  	for (int i = 0; i != k; i++) {
  		wsigs.push_back(sqrt(sig[i*k + i]));
  		gsl_vector_set(hampar, offset + i, log(wsigs[i]));
  	}

  	gsl_matrix_view t1 = gsl_matrix_view_array(sig, k, k);
  	gsl_matrix_memcpy(w, &t1.matrix);
  	gsl_linalg_cholesky_decomp1(w);

  	gsl_matrix* wi = gsl_matrix_alloc(k, k);
  	gsl_matrix_memcpy(wi, w);
  	gsl_linalg_tri_lower_invert(wi);

  	gsl_vector_view stackedv = gsl_vector_subvector(hampar, ioff, indi * k);
  	gsl_matrix_view stackedm = gsl_matrix_view_vector(&stackedv.vector, indi, k);

  	gsl_blas_dtrmm(CblasRight, CblasLower, CblasTrans, CblasNonUnit, 1.0, wi, &stackedm.matrix);

  	for (int i = 0; i != k; i++)
  		for (int j = 0; j <= i; j++)
  			gsl_matrix_set(w, i, j, gsl_matrix_get(w, i, j) / wsigs[i]);

  	gsl_matrix_free(wi);
  }

  //from Cholesky factor correlation matrix to partial correlations
  void from_w_to_z(int flag, std::vector<double>& z, gsl_matrix* w) {
  	int k;
  	if (flag == 0) {
  		k = icompg;
  	}
  	else {
  		k = respno;
  	}
  	//lower triangular
  //	int jj = 0;
  	for (int i = 1; i < k; ++i) {
  		z.push_back(gsl_matrix_get(w, i, 0));
  		double sum_sqs = gsl_pow_2(gsl_matrix_get(w, i, 0));
  		for (int j = 1; j < i; ++j) {
  			double temp = sqrt(1.0 - sum_sqs);
  			if (temp > 0) z.push_back(gsl_matrix_get(w, i, j) / temp); else z.push_back(0.0);
  			sum_sqs += gsl_pow_2(gsl_matrix_get(w, i, j));
  		}
  	}
  }

  //from partical correlations to model parameters
  void from_z_to_y(int flag, gsl_vector* hampar, const std::vector<double>& z) {
  	int k; int offset = nhamil;
  	if (flag == 0) {
  		k = icompg;
  	}
  	else {
  		k = respno;
  		offset += (icompg * (icompg - 1)) / 2 + icompg;
  	}
  	int kc2 = (k * (k - 1)) / 2;
  	for (int i = 0; i != kc2; i++) {
  		double temp = z[i];
  		gsl_vector_set(hampar, offset + i, 0.5 * log((1 + temp) / (1 - temp)));
  	}
  }

  //putting the above together: from variance-covariance matrices to model parameters
  void make_hampar_from_sig(int flag, double* sig, gsl_vector* hampar) {
  	int k;
  	if (flag == 0) k = icompg;
  	else {
  		k = respno;
  	}
  	gsl_matrix* w = gsl_matrix_calloc(k, k);
  	from_sig_to_w(flag, hampar, w, sig);
  	std::vector<double> z;
  	from_w_to_z(flag, z, w);
  	from_z_to_y(flag, hampar, z);

  	gsl_matrix_free(w);
  }

  //contributions to overall model likelihood from LKJ-related parameters
  double joint_likeli4(int flag, gsl_vector* hampar, const std::vector<double> &z, gsl_matrix* w, double eta,
  	double taut, double liknorm4) {
  	//Stan reference manual & Stan Functions reference
  	// lkj: det(corr)^(eta-1); jacoby transfrom (1-z^2_ij)^0.5*(K-i-1); dy/dz 1-z^2_ij; MVN -0.5 * log(det(sigma));
  	int k; int offset = nhamil + (icompg * (icompg - 1)) / 2;
  	// extrahiere standardabweichungen
  	if (flag == 0) {
  		k = icompg;
  	}
  	else {
  		k = respno;
  		offset += icompg + (respno * (respno - 1)) / 2;
  	}
  	std::vector<double> wsigs;
  	for (int i = 0; i != k; i++) wsigs.push_back(exp(gsl_vector_get(hampar, offset + i)));

  	std::vector<double> logw; double logwsum = 0.0;
  	for (int i = 1; i != k; i++) {
  		double temp = log(gsl_matrix_get(w, i, i));
  		logw.push_back(temp);
  		logwsum += temp;
  	}

  	double temp = 0.0;
  	for (int i = 1; i != k; i++) temp -= (i + 1) * logw[i - 1];
  	temp += (k + 2 * eta - 2) * logwsum;
  //		(k - (i + 1) + 2 * eta - 2) * log(gsl_matrix_get(w, i, i));
  // -2 - indi statt -2 wegen mvn - 0.5 log(detsig), hier r�ckg�ngig wegen Parametrisierung als Sigma * Standard Normal

  	int jj = 0;
  	for (int i = 1; i != k; i++) {
  		temp += log1p(-gsl_pow_2(z[jj++]));
  		for (int j = 1; j != i; j++) {
  			//			temp += log1p(-gsl_pow_2(z[jj++]));
  			temp += log1p(-gsl_pow_2(z[jj]));
  			//			double zwi = 0.0;
  			//			for (int l = 0; l != j; l++) zwi += gsl_pow_2(gsl_matrix_get(w, i, l));
  			//			temp += 0.5 * log1p(-zwi);
  			if (z[jj] != 0.0) temp += log(gsl_matrix_get(w, i, j) / z[jj]);
  			else {
  				double zwi = 0.0;
  				for (int l = 0; l != j; l++) zwi += gsl_pow_2(gsl_matrix_get(w, i, l));
  				temp += 0.5 * log1p(-zwi);

  			}
  			jj++;
  		}
  	}

  	for (int i = 0; i != k; i++) {
  		double x = gsl_vector_get(hampar, offset + i);
  		double sig = wsigs[i]; double sigtau = gsl_pow_2(sig / taut);
  		temp += x -  log1p(sigtau); //half_cauchy on sig-diagonals * jacoby factor
  	}
  	temp -= liknorm4;
  	return(temp);
  }

  //contributions to overall model likelihood from Omega^2
  double joint_likeli5(gsl_vector* hampar, double* explambdas, double liknorm6) {
  	double x = gsl_vector_get(hampar, n_all_parameters-1);
  //	double priordf = 2.0;

  	double alpha = (indi * priordf) / 2.0 + prioralpha, beta = priorbeta;

  	double temp;
  	temp = (alpha)*x - beta * exp(x);
  	temp -= liknorm6;
  	return temp;
  }

  //derivatives cholesky factor by partial correlations
  void dwdz(int flag, const std::vector<double> &z, gsl_matrix* w, std::vector<gsl_matrix*>& dwz) {
  	int k; if (flag == 0) k = icompg; else k = respno;
  	int jj = 0;

  	for (int i = 1; i != k; i++) {
  		double sum_sqs = gsl_pow_2(gsl_matrix_get(w, i, 0));
  		gsl_matrix_set(dwz[0], i, 0, 1.0); jj++;
  		for (int j = 1; j <= i; j++) {
  			for (int jd = 0; jd != j; jd++) {
  				double temp = 0;
  				for (int l = jd; l <= j - 1; l++) {
  					temp += gsl_matrix_get(w, i, l) * gsl_matrix_get(dwz[jd], i, l);
  				}
  				if (j < i) temp *= (sum_sqs < 1.0)? z[jj] / sqrt(1.0 - sum_sqs): 0.0;
  				else temp /= sqrt(1.0 - sum_sqs);
  				gsl_matrix_set(dwz[jd], i, j, -temp);
  			}
  			if (j < i) {
  				gsl_matrix_set(dwz[j], i, j, sqrt(1.0 - sum_sqs));
  				sum_sqs += gsl_pow_2(gsl_matrix_get(w, i, j));
  				jj++;
  			}
  		}
  	}
  }

  //derivatives of diffusion models and lkj prior by involved parameters; Jacoby factors diffusion; lkj contributions (flag = 0: diffusion-model part; flag = 1: motor-times part)
  void dmvnlkjdy(int flag, const std::vector<double>& z, double eta, gsl_vector* hampar, const std::vector<gsl_matrix*>& dwz, gsl_matrix* w, gsl_vector* dhampar) {
  	int k, ioff; int offset = nhamil + (icompg * (icompg - 1)) / 2;
  	int offset2 = nhamil;

  	if (flag == 0)
  	{
  		k = icompg;
  		ioff = iavwoff;
  	}
  	else {
  		k = respno;
  		ioff = ilamoff;
  		offset += icompg + (respno * (respno - 1)) / 2;
  		offset2 += (icompg * (icompg - 1)) / 2 + icompg;
  	}

  	//d diffusion dw und dsigs

  	gsl_matrix* temp = gsl_matrix_calloc(k, k);
  	gsl_vector* wsigs = gsl_vector_alloc(k);
  	for (int i = 0; i != k; i++) gsl_vector_set(wsigs, i, exp(gsl_vector_get(hampar, offset + i)));

  	gsl_matrix* sd = gsl_matrix_calloc(k, k);
  	gsl_vector_view sddiag = gsl_matrix_diagonal(sd);
  	gsl_vector_memcpy(&sddiag.vector, wsigs);

  	gsl_vector_view tempdiag = gsl_matrix_diagonal(temp);
  	gsl_vector_memcpy(&tempdiag.vector, wsigs);

  	gsl_matrix* sum0 = gsl_matrix_calloc(k, k);
  	gsl_vector* sum = gsl_vector_calloc(k);
  	gsl_vector* zwt = gsl_vector_alloc(k);

  	for (int t = 0; t != indi; t++) {
  		gsl_vector_view dindi = gsl_vector_subvector(dhampar, ioff + t * k, k);
  		gsl_vector_view thetindi = gsl_vector_subvector(hampar, ioff + t * k, k);
  		gsl_blas_dger(1.0, &dindi.vector, &thetindi.vector, sum0);

  		gsl_vector_memcpy(zwt, &thetindi.vector);
  		gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, w, zwt);
  		gsl_vector_mul(zwt, &dindi.vector);
  		gsl_vector_add(sum, zwt);
  	}

  	gsl_blas_dtrmm(CblasRight, CblasLower, CblasNoTrans, CblasNonUnit, 1.0, sum0, temp);
  	gsl_matrix_free(sum0);

  	gsl_vector_free(zwt);
  	// chain rule
  	gsl_vector_mul(sum, wsigs);
  	gsl_vector_view dsigs = gsl_vector_subvector(dhampar, offset, k);
  	gsl_vector_add(&dsigs.vector, sum);
  	gsl_vector_free(sum);


  	//d diffusion dthet

  	gsl_vector_view dthetv = gsl_vector_subvector(dhampar, ioff, indi * k);
  	gsl_matrix_view dthetm = gsl_matrix_view_vector(&dthetv.vector, indi, k);



  	gsl_blas_dtrmm(CblasRight, CblasLower, CblasNoTrans, CblasNonUnit, 1.0, w, sd);

  	gsl_blas_dtrmm(CblasRight, CblasLower, CblasNoTrans, CblasNonUnit, 1.0, sd, &dthetm.matrix);
  	gsl_matrix_free(sd);

  	// standard normal
  	gsl_vector_view thetv = gsl_vector_subvector(hampar, ioff, indi * k);
  	gsl_vector_add(&dthetv.vector, &thetv.vector);


  	// d/dz log(LKJ)
  	int jj = 0;
  	for (int s = 1; s != k; s++) for (int t = 0; t != s; t++) {
  		double zij = 1 - gsl_pow_2(z[jj]);
  		double x = 0;
  		for (int l = t; l <= s; l++)  x += gsl_matrix_get(temp, s, l) * gsl_matrix_get(dwz[t], s, l);

  		x -= (k - (s + 1) + 2 * eta - 2) * gsl_matrix_get(dwz[t], s, s) / gsl_matrix_get(w, s, s);

  		for (int l = t + 1; l < s; l++) x -= (gsl_matrix_get(w, s, l) != 0.0 || gsl_matrix_get(dwz[t], s, l) != 0.0) ? gsl_matrix_get(dwz[t], s, l) / gsl_matrix_get(w, s, l) : 0.0;

  		x *= zij; x += 2 * z[jj];
  //		if (isnan(x)) std::cout << "dmvndlkj";
  		gsl_vector_set(dhampar, offset2 + jj, x);
  		jj++;
  	}
  	gsl_vector_free(wsigs);
  	gsl_matrix_free(temp);
  }

  //derivatives of half-Cauchy priors by (log) standard deviations
  void dhudsigs(int flag, gsl_vector* hampar, double tau, gsl_vector* dhampar) {
  	int	k, offset = nhamil + (icompg * (icompg - 1))/2;
  	if (flag == 0)
  	{
  		k = icompg;
  	}
  	else {
  		k = respno;
  		offset += icompg + (respno * (respno-1))/2;
  	}

  //cauchy-prior

  	for (int i = 0; i != k; i++) {
  		double x = gsl_vector_get(hampar, offset + i);
  		double sig = gsl_pow_2(exp(x) / tau);
  		gsl_vector_set(dhampar, offset + i, -1.0  + 2 * sig / (1 + sig));
  // + indi wegen d II /d lambda hier r�ckg�ngig
  	}
  }

  //putting the above derivatives together
  void dhudext(gsl_vector* hampar, double* explambda, const std::vector<double>& zt, const std::vector<double>& zr, gsl_matrix* wt, gsl_matrix* wr,
               double etat, double etar, gsl_vector* dhampar) {

    // d / dsigs
    for (int flag = 0; flag != 2; flag++) {
      int k;
      if (flag == 0) k = icompg; else k = respno;
      // d / dsigs
      if (flag == 0) dhudsigs(flag, hampar, taut, dhampar);
      else dhudsigs(flag, hampar, taur, dhampar);

      std::vector<gsl_matrix*> dwz;
      dwz.clear();
      for (int l = 0; l != k; l++) dwz.push_back(gsl_matrix_calloc(k, k));
      dwdz(flag, (flag == 0) ? zt : zr, (flag == 0) ? wt : wr, dwz);


      dmvnlkjdy(flag, (flag == 0) ? zt : zr, (flag == 0) ? etat : etar, hampar, dwz, (flag == 0) ? wt : wr, dhampar);


      for (int l = 0; l != k; l++) gsl_matrix_free(dwz[l]);
    }

    // d/log(domega-square)

    double x = gsl_vector_get(hampar, n_all_parameters - 1);
    //	double priordf = 2.0;
    double s = 0.0;
    for (int t = 0; t != indi; t++) s += (phase >= 3) ? 1 / gsl_pow_2(explambda[t]) : 1 / gsl_pow_2(gsl_vector_get(hampar, isigoff + t));

    s = s * priordf; //	s=s*adf;
    double alpha = (indi * priordf) / 2.0 + prioralpha, beta = s / 2.0 + priorbeta;
    gsl_vector_set(dhampar, n_all_parameters - 1, -alpha + beta * exp(x));
  }

}
