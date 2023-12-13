#include "rts.h"

namespace drtmpt {

  //log density normal
  double lognormal(double x) {
  	return	-0.5*x*x - M_LN_SQRT_PI - 0.5*M_LN2;
  }

  //log Mill's ratio
  double logMill(double x) {
  	double m;
  	if (x > 1.0e5) return -log(x);
  	m = lnnorm(-x) - lognormal(x);
  	return m;
  }

  //number of components needed in Wiener density small time series; from Gondan et al. "evenfaster.R"
  double ks(double t, double w, double eps)
  {
  	double K1 = (sqrt(2 * t) + w) / 2;
  	double u_eps = fmin(-1, M_LN2 + M_LNPI + 2 * log(t) + 2 * (eps)); // # Safe bound so that

  	double	arg = -t * (u_eps - sqrt(-2 * u_eps - 2.0)); //# sqrt(x) with x > 0
  	double 	K2 = (arg > 0) ? 0.5 * (sqrt(arg) - w) : K1;
  	return ceil(fmax(K1, K2));
  }

  //log Wiener density; small time series
  double logfsw(double t, double w, int K)
  {
  	if (w == 0) return GSL_NEGINF;
  	double	fplus = GSL_NEGINF, fminus = GSL_NEGINF, twot = 2 * t;
  	if (K > 0)
  		for (int k = K; k >= 1; k--) {
  			double temp1 = w + 2 * k, temp2 = w - 2 * k;

  			fplus = logsum(log(temp1) - gsl_pow_2(temp1) / twot, fplus);
  			fminus = logsum(log(-temp2) - gsl_pow_2(temp2) / twot, fminus);
  		}
  	fplus = logsum(log(w) - gsl_pow_2(w) / twot, fplus);
  	return  -0.5 * M_LN2 - M_LN_SQRT_PI - 1.5 * log(t) + logdiff(fplus, fminus);
  }


  ///number of components needed in Wiener density large time series; nach Rwiener.R, source file dwiener.c
  double kl(double q, double v, double w, double err) {
  	// calculate number of terms needed for large t
  	double K1 = 1.0 / (M_PI * sqrt(q)), K2=0.0;
  	double temp = -2 * (log(M_PI * q) + err);
  	if (temp>=0) K2 = sqrt(temp/(gsl_pow_2(M_PI) * q));
  	return ceil(fmax(K1,K2));
  }

  //log Wiener density; large time series
  double logfl(double q, double v, double w, int K) {

  	if (w == 0) return GSL_NEGINF;
  	double fplus = GSL_NEGINF, fminus = GSL_NEGINF;

  	double halfq = q / 2;
  	for (int k = K; k >= 1; k--) {
  		double temp = k * M_PI;
  		double check = sin(temp * w);
  		if (check > 0) fplus = logsum(log(k) - gsl_pow_2(temp) * halfq + log(check), fplus);
  		else fminus = logsum(log(k) - gsl_pow_2(temp) * halfq + log(-check), fminus);
  	}
  	return	logdiff(fplus, fminus) + M_LNPI;
  }

  //log Wiener density
  double dwiener_d(double q, double a, double vn, double wn, double eps)
  {
  	double kll, kss, ans, v, w;
  	double errziel = eps;
  	double err = eps * 1.1;
  	int zahl = 0;

  	if (q >= 0) {
  		w = 1.0 - wn;
  		v = -vn;
  	}
  	else {
  		q = fabs(q);
  		w = wn;
  		v = vn;
  	}

  	double q_asq = q / gsl_pow_2(a);

  	double lg1 = (-v * a * w - (gsl_pow_2(v)) * q / 2) - 2 * log(a);
  NEW:	ans = 0;
  	double es = (err - lg1);
  	kss = ks(q_asq, w, es);
  	double el = es;
  	kll = kl(q_asq, v, w, el);


  	if (2 * kss < kll)  // if small t is better
  		ans = lg1 + logfsw(q_asq, w, static_cast<int>(kss));
  	else
  		// if large t is better...
  		ans = lg1 + logfl(q_asq, v, w, static_cast<int>(kll));

  	zahl++; // MONITOR(0, 5)++;
  	if (zahl == 10) {
  //		std::cout << "Zahl = 10 " << setw(20) << q << setw(20) << a << setw(20) << vn << setw(20) << wn << setw(20) << eps << std::endl;
  		return ans;
  	}
  	if (err - ans > errziel) {
  		//		MONITOR(1, 5)++;
  		err = (isfinite(ans)) ? errziel * (1 + zahl * 0.1) + ans : 2 * err;
  		goto NEW;
  	}
  	return ans;
  }



  // number of components for Wiener cumulative distribution function, small time series; from Gondan et al. ever faster
  double Ks(double t, double v, double a, double w, double eps)
  {
  	double K1 = 0.5 * (fabs(v) / a * t - w);
  	double arg = fmax(0, fmin(1, exp(v*a*w + gsl_pow_2(v)*t / 2 + (eps)) / 2));
  	double K2 = (arg==0)?GSL_POSINF:(arg==1)?GSL_NEGINF: -sqrt(t) / 2 / a * gsl_cdf_ugaussian_Pinv(arg);
  	return ceil(fmax(K1, K1 + K2));
  }

  // log Wiener cumulative distribution function, small time series
  double logFs(double t, double v, double a, double w, int K)
  {
  	double	fplus = GSL_NEGINF, fminus = GSL_NEGINF;
  	double sqt = sqrt(t), temp = -v * a * w - gsl_pow_2(v) * t / 2;
  	double vt = v * t;
  	for (int k = K - 1; k >= 0; k--)
  	{
  		double rj = a * (2 * k + w);
  		double dj = lognormal(rj / sqt);
  		double pos1 = dj + logMill((rj - vt) / sqt);
  		double pos2 = dj + logMill((rj + vt) / sqt);
  		fplus = logsum(logsum(pos1, pos2), fplus);
  		rj = a * (2.0 * k + 2.0 - w);
  		dj = lognormal(rj / sqt);
  		double neg1 = dj + logMill((rj - vt) / sqt);
  		double neg2 = dj + logMill((rj + vt) / sqt);
  		fminus = logsum(logsum(neg1, neg2), fminus);
  	}
  	return logdiff(fplus, fminus) + temp;
  }

  // number of components for Wiener cumulative distribution function, large time series; following Burton et al.
  double Kl(double t, double v, double a, double w, double err) {
  	double api = a / M_PI, vsq = gsl_pow_2(v);
  	double sqrtL1 = sqrt(1 / t) * api;
  	double sqrtL2 = sqrt(fmax(1.0, -2 / t * gsl_pow_2(api) * (err+log(M_PI*t / 2 * (vsq + gsl_pow_2(M_PI / a))) + v * a * w + vsq * t / 2)));
  	return ceil(fmax(sqrtL1, sqrtL2));
  }

  // log Wiener cumulative distribution function, large time series
  double logFl(double q, double v, double a, double w, int K)
  {
  	double fplus = GSL_NEGINF, fminus = GSL_NEGINF;
  	double la = log(a), lv = log(fabs(v));

  	double F = GSL_NEGINF;

  	for (int k = K; k >= 1; k--) {
  		double temp0 = log(k * 1.0), temp1 = k * M_PI, temp2 = temp1 * w;
  		double check = sin(temp2);
  		if (check > 0) {
  			double temp = temp0 - logsum(2 * lv, 2 * (temp0 + M_LNPI - la)) - 0.5 * gsl_pow_2(temp1 / a) * q + log(check);
  			fplus = logsum(temp, fplus);
  		}
  		else if (check < 0)
  		{
  			double temp = temp0 - logsum(2 * lv, 2 * (temp0 + M_LNPI - la)) - 0.5 * gsl_pow_2(temp1 / a) * q + log(-check);
  			fminus = logsum(temp, fminus);
  		}
  	}
  	F = logdiff(fplus, fminus);
  	// Note: logFl is  really log(prob_upperbound + 2*pi/a/a*exp(F))
  	return F - v * a * w - 0.5 * gsl_pow_2(v) * q;
  }

  //// log Wiener cumulative distribution function
  double logFjoint_lower(double q, double a, double v, double w) {
  	if (q == 0) return(GSL_NEGINF);
  	double Kll, Kss, ans;
  	double err = accuracy;
  	if (gsl_isinf(q)) return logprob_upperbound(0, a, v, w);

  	double es = err;
  	Kss = Ks(q, v, a, w, es);
  	double lg = M_LN2 + M_LNPI - 2.0 * log(a);
  	double el = es;

  	Kll = Kl(q, v, a, w, el);

    	if (3 * Kss < Kll)
  	{
  		ans = logFs(q, v, a, w, static_cast<int>(Kss));
  	}
    	else
  	{
  		ans = logdiff(logprob_upperbound(0, a, v, w), lg + logFl(q, v, a, w, static_cast<int>(Kll)));
  	}
  	return ans;
  }

  //log of probability of reaching bound pm in Wiener diffusion
  double logprob_upperbound(int pm, double a, double v, double w)
  {
  	if (pm == 1) { v = -v; w = 1.0 - w; }
  	if (fabs(v) == 0.0) return log1p(-w);


  	double prob;
  	double e = (-2.0 * v * a * (1.0 - w));
  	if (e < 0) {
  		double tt;
  		tt = log1pem1(e) - logdiff(2 * v * a * w, e);
  		prob = tt;
  	}
  	else {
  		double tt;
  		tt = log1pem1(-e) - log1pem1(2 * v * a);
  		prob = tt;
  	}
  	return prob;
  }

  //mean first passage times at lower bound
  double lower_bound_time(double a, double vn, double wn) {
  	double temp, amw = a * (1 - wn);
  	if (abs(vn) < 1e-5) {
  		temp = (gsl_pow_2(a) - gsl_pow_2(amw)) / 3.0;

  	}
  	else {
  		temp = a  / tanh(a*vn) - amw  /tanh(vn*amw);
  		temp /= vn;
  	}
  	return temp;
  }

  //mean first passage times at bound pm
  double exp_mean(int pm, double a, double v, double w) {
  	if (pm == 1) { v = -v; w = 1 - w; }
  	return lower_bound_time(a, v, w);

  }

  //Phi-function from Grasman et al., JMP, 2009
  double phi(double x, double y, double v) {
  	double vx = 2 * v * x, vy = 2 * v * y;
  	if (vy > vx) return (exp(vx) * (gsl_expm1(vy - vx)));
  	else return (-exp(vy) * (gsl_expm1(vx - vy)));
  }


  //Variance of first passage times at lower bound
  double lower_bound_var(double a, double vn, double wn) {
  	double z = a * wn, phiza = phi(z, a, vn), vn3 = gsl_pow_3(vn);
  	double temp = (-2 * a * phi(0, z, vn) * (2 * vn * a * phi(z, 2 * a, vn) + phi(0, a, vn) * phiza)) * exp(2 * vn * a) / (vn3 * gsl_pow_2(phi(0, a, vn) * phiza));
  	temp += (4 * vn * z * (2 * a - z) * exp(2 * vn * (z + a)) + z * phi(2 * z, 2 * a, vn)) / vn3 / gsl_pow_2(phiza);
  	if (temp <= 0) {
  //		std::cout << "! " << setw(20) << a << setw(20) << vn << setw(20) << wn << setw(20) << temp << std::endl;
  		temp = 0.1;
  	}
  	return temp;
  }

}
