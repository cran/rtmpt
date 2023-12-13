#include "rts.h"

namespace drtmpt {

  // compute needed number of components to sum for small time series for derivatives by a and t
  double dtks(double t, double w, double eps)
  {
  	// t = t/a^2; eps = log(eps * a^2)
  	double K1 = (sqrt(3 * t) + w) / 2;
  	double u_eps = fmin(-1, (log(8.0 / 27.0) + M_LNPI + 4.0*log(t) + 2.0*eps)/3); // # Safe bound so that
  	double	arg = -3 * t * (u_eps - sqrt(-2 * u_eps - 2)); //# sqrt(x) with x > 0
  	double 	K2 = (arg > 0) ? 0.5 * (sqrt(arg) - w) : K1;
  	return ceil(fmax(K1, K2));
  }

  //derivative of log(f_wiener) by t; small time series
  void logdtfsw(double t, double w, int K, double& erg, int& newsign)
  {
  	double	fplus = GSL_NEGINF, fminus = GSL_NEGINF, twot=2*t;
  	{
  		if (K > 0)
  			for (int k = K; k >= 1; k--) {
  				double temp1 = w + 2 * k, temp2 = w - 2 * k;
  				fplus = logsum(3 * log(temp1) - gsl_pow_2(temp1) / twot, fplus);
  				fminus = logsum(3 * log(-temp2) - gsl_pow_2(temp2) / twot, fminus);
  			}
  	}
  	fplus = logsum(3 * log(w) - gsl_pow_2(w) / twot, fplus);
  	erg = (fplus > fminus) ? logdiff(fplus, fminus) : logdiff(fminus, fplus);
  	newsign = (fplus > fminus) ? 1 : -1;
  }

  // compute needed number of components to sum for large time series for derivatives by a and t
  double dtkl(double q, double v, double a, double err) {
  	// t = t/a^2 eps = log(eps * a^2)

  	double K1 = sqrt(3 / q)/M_PI;

  	double u_eps = fmin(-1, err + log(0.6) + M_LNPI + 2.0 * log(q) ); // # Safe bound so that

  	double	arg = -2.0/M_PISQ/q*(u_eps - sqrt(-2 * u_eps - 2)); //# sqrt(x) with x > 0
  	double 	kl = (arg > 0) ? sqrt(arg) : K1;
  	return ceil(fmax(kl,K1));
  }


  //derivative of log(f_wiener) by t; large time series
  void logdtfl(double q, double w, int K, double& erg, int& newsign) {
  	double fplus = GSL_NEGINF, fminus = GSL_NEGINF;

  	double  halfq = q / 2;
  	for (int k = K; k >= 1; k--) {
  		double temp = M_PI * k, zwi = sin(temp * w);
  		if (zwi > 0) {
  			fplus = logsum(3 * log(k) - gsl_pow_2(temp) * halfq + log(zwi), fplus);

  		}
  		if (zwi < 0) {
  			fminus = logsum(3 * log(k) - gsl_pow_2(temp) * halfq + log(-zwi), fminus);

  		}
  	}
  	erg = (fplus > fminus) ? logdiff(fplus, fminus) : logdiff(fminus, fplus);
  	newsign = (fplus > fminus) ? 1 : -1;
  }

  //derivative of log(f_wiener) by a
  double dadwiener_d(double q, double a, double vn, double wn, double d)
  {
  	double kll, kss, ans, v, w;
  	double errziel = accuracy;
  	double err = errziel*1.2;
  	double la = log(a), lq = log(fabs(q));


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

  	double ans0 =  - v * w;
  	double lg1 = -v * a*w - gsl_pow_2(v)*q / 2 - 2 * la;
  	double factor = lg1 - 3 * la;


  NEW:	double es = (err - lg1 + d);
  	es = es + la; //wegen q_asq in daks
  	es = es -  M_LN2 +  2*la -  lq; // wegen Verwendung von dtks anstelle von daks

  	kss = dtks(q_asq, w, es);

  	double el = (err - lg1 + d);
  	el = el + la; 	// wegen q_asq in dakl
  	el = el - M_LN2 + 2*la - lq; // wegen verwendung von dtkl anstelle von dakl
  	kll = dtkl(q_asq, v, a, el);

  	// if small t is better
  	if (2 * kss < kll) {
  		double erg; int newsign;
  		logdtfsw(q_asq, w, static_cast<int>(kss), erg, newsign);
  		ans = ans0 + 1 / a - newsign * exp(-0.5*M_LN2 - M_LN_SQRT_PI - 2.5*lq + 4 * la + lg1 + erg - d);
  	}
  	// if large t is better...
  	else {
  		double erg; int newsign;
  		logdtfl(q_asq, w, static_cast<int>(kll), erg, newsign);
  		ans = ans0 - 2.0 / a + newsign * exp(lq + factor + 3 * M_LNPI + erg - d);
  	}

  	double temp = log(fabs(ans)) + d;
  	if (temp < d) {
  		double check = temp - d;
  		if (err - check > errziel) {
  			err = errziel*1.2 + check;
  			goto NEW;
  		}
  	}
  	double check = temp + M_LN2 - d;

  //	MONITOR(0, 7)++;
  	if (err + check > errziel) {
  		err = errziel*1.2 - check;
  		d = dwiener_d(-q, a, v, w, err);
  		goto NEW;
  	}
  	return ans;
  }

  //derivative of log(f_wiener) by t
  double dtdwiener_d(double q, double a, double v, double w, double &d)
  {
  	double kll, kss, ans;
  	double errziel = accuracy;
  	double err = errziel*1.3;


  	double q_asq = q / gsl_pow_2(a);

  	double ans0 = -gsl_pow_2(v) / 2.0;
  	double la = 2 * log(a);
  	double lg1 = -v * a * w - gsl_pow_2(v) * q / 2 - la;
  	double factor = lg1 - la;




  NEW:	double es = err - lg1 + d;
  	es = es + la; //wegen q_asq in dtk1
  	kss = dtks(q_asq, w, es);
  	double el = (err - lg1 + d) ;
  	el = el +la; 	// wegen q_asq in dtkl
  	kll = dtkl(q_asq, v, a, el);

  	// if small t is better
  	if (2*kss < kll)
  	{
  		double erg; int newsign;
  		logdtfsw(q_asq, w, static_cast<int>(kss), erg, newsign);
  		ans = ans0 - 1.5 / q + newsign * exp(factor - 1.5 * M_LN2 - M_LN_SQRT_PI - 3.5 * log(q_asq) + erg - d);
  	}
  	// if large t is better...
  	else
  	{
  		double erg; int newsign;
  		logdtfl(q_asq, w, static_cast<int>(kll), erg, newsign);
  		ans = ans0 - newsign * exp(factor + 3 * M_LNPI - M_LN2 + erg - d);

  	}

  	// error |f' - f'True| / f'
  	double temp = log(fabs(ans)) +d;
  	if (temp < d) {
  		double check = temp - d;
  		if (err - check > errziel) {
  			err = errziel*1.3 + check;
  			goto NEW;
  		}
  	}
  	// error |f'/f - f'True/fTrue|
  	double check = temp + M_LN2 - d;

  //	MONITOR(0, 6)++;

  	if (err  + check > errziel) {
  //		MONITOR(1, 6)++;
  		err = errziel*1.3 - check;
  		d = dwiener_d(-q, a, v, w, err);
  		goto NEW;
  	}
  	return ans;
  }

  // compute needed number of components to sum for small time series for derivatives by w
  double dwks(double t, double w, double eps)
  {
  	double K1 = (sqrt(3 * t) + w) / 2;
  	double u_eps = fmin(-1, 2 * (eps)+M_LN2 + M_LNPI + 2.0 * log(t)); // # Safe bound so that
  	double	arg = -t * (u_eps - sqrt(-2 * u_eps - 2)); //# sqrt(x) with x > 0
  	double 	K2 = (arg > 0) ? (sqrt(arg) + w) / 2 : K1;
  	return ceil(fmax(K1, K2));
  }

  //derivative of log(f_wiener) by w; small time series
  void logdwfsw(double t, double w, int K, double& erg, int& sign)
  {
  	double	fplus = GSL_NEGINF, fminus = GSL_NEGINF, twot = 2 * t;
  	for (int k = K; k >= 1; k--) {
  		double temp1 = gsl_pow_2(w + 2 * k), temp2 = gsl_pow_2(w - 2 * k), temp3 = temp1 - t, temp4 = temp2 - t;
  		if (temp3 > 0) fplus = logsum(log(temp3) - temp1 / twot, fplus);
  		else if (temp3 < 0) fminus = logsum(log(-(temp3)) - temp1 / twot, fminus);
  		if (temp4 > 0) fplus = logsum(log(temp4) - temp2 / twot, fplus);
  		else if (temp4 < 0) fminus = logsum(log(-(temp4)) - temp2 / twot, fminus);
  	}
  	double temp = gsl_pow_2(w), temp1 = temp - t;
  	if (temp1 > 0) fplus = logsum(log(temp1) - temp / twot, fplus);
  	else if (temp1 < 0) fminus = logsum(log(-(temp1)) - temp / twot, fminus);
  	erg = (fplus > fminus) ? logdiff(fplus, fminus) : logdiff(fminus, fplus);
  	sign = (fplus < fminus) ? -1 : 1;
  	// fehlt 1/sqrt(2*pi*t^5)
  }


  // compute needed number of components to sum for large time series for derivatives by w
  double dwkl(double q, double v, double err)
  {

  	double K1 = sqrt(2 / q) / M_PI;
  	double u_eps = fmin(-1, log(4.0 / 9.0) + 2 * M_LNPI + 3.0 * log(q) + 2.0 * (err)); // # Safe bound so that
  	double	arg = -(u_eps - sqrt(-2 * u_eps - 2)); //# sqrt(x) with x > 0
  	double 	K2 = (arg > 0) ? 1.0 / M_PI * sqrt(arg / q) : K1;
  	return ceil(fmax(K1, K2));
  }

  //derivative of log(f_wiener) by w; large time series
  void logdwfl(double q, double v, double w, int K, double& erg, int& sign) {

  	double fplus = GSL_NEGINF, fminus = GSL_NEGINF;

  	double  halfq = q / 2.0;

  	for (int k = K; k >= 1; k--) {
  		double temp = M_PI * k;
  		double x = cos(temp * w);
  		if (x > 0)
  			fplus = logsum(2 * log(k) - gsl_pow_2(temp) * halfq + log(x), fplus);
  		else if (x < 0)
  			fminus = logsum(2 * log(k) - gsl_pow_2(temp) * halfq + log(-x), fminus);
  	}
  	erg = (fplus > fminus) ? logdiff(fplus, fminus) : logdiff(fminus, fplus);
  	sign = (fplus < fminus) ? -1 : 1;
  	// fehlt pi^2
  }


  //derivative of log(f_wiener) by w
  double dwdwiener_d(double q, double a, double vn, double wn, double d)
  {
  	double kll, kss, ans, v, w;
  	double errziel = accuracy;
  	double err = errziel * 1.2;

  	int sign = 1;

  	if (q >= 0) {
  		w = 1 - wn;
  		v = -vn;
  		sign = -1;
  	}
  	else {
  		q = fabs(q);
  		w = wn;
  		v = vn;
  	}


  	double q_asq = q / gsl_pow_2(a);
  	double ans0 = -v * a;
  	double lg1 = (-v * a * w - gsl_pow_2(v) * (q) / 2) - 2 * log(a);
  	double ls = -lg1 + d;
  	double ll = -lg1 + d;
  NEW:	double es = (err);
  	kss = dwks(q_asq, w, es + ls);

  	kll = dwkl(q_asq, v, es + ll);

  		 // if small t is better
  	if (2 * kss < kll) {

  		double erg; int newsign;
  		logdwfsw(q_asq, w, static_cast<int>(kss), erg, newsign);
  		ans = ans0 - newsign * exp(erg - ls - 2.5 * log(q_asq) - .5 * M_LN2 - .5 * M_LNPI);
  	}
  	// if large t is better...
  	else {

  		double erg; int newsign;
  		logdwfl(q_asq, v, w, static_cast<int>(kll), erg, newsign);
  		ans = ans0 + newsign * exp(erg - ll + 2.0 * M_LNPI);
  	}
  //	if (!gsl_finite(ans))
  //		std::cout << "w " << a << setw(20) << vn << setw(20) << wn << setw(20) << q << std::endl;

  	// relative error in d/dw f (not d/dw log(f))
  	double temp = log(fabs(ans)) + d;
  	if (temp < d) {
  		double check = temp - d;
  		if (err - check > errziel) {
  			err = errziel * 1.2 + check;
  			goto NEW;
  		}
  	}
  	double check = temp + M_LN2 - d;

  	//  	MONITOR(0, 8)++;

  	if (err + check > errziel) {
  		err = errziel * 1.2 - check;
  		d = dwiener_d(-q, a, v, w, err);
  		goto NEW;
  	}
  	return ans * sign;
  }

  //helper function derivative of logprob_upperbound by a and v
  double davlogprob_upperbound(int pm, double a, double v, double w) {
  	const double em1 = log1p(-1.1 * 1.0e-8);
  	double tt;
  	if (pm == 1) {
  		w = 1.0 - w;
  		v = -v;
  	}

  	if (fabs(v) == 0.0) return(-w);

  	if (v < 0) {


  		double emw = (2.0 * v * a * (1.0 - w)), ew = 2 * a * v * w, e = 2 * a * v;

  		if (((emw >= em1) || (ew >= em1)) || (e >= em1)) return(-w);
  		tt = M_LN2 + emw - log1pem1(emw);
  		double temp = log1pem1(ew) - log1pem1(e);
  		double lw = log(w);
  		if (lw > temp) {
  			tt += logdiff(lw, temp);
  			tt = exp(tt);
  		}
  		else {
  			tt += logdiff(temp, lw);
  			tt = -exp(tt);
  		}


  	}
  	else {
  		double emw = (-2.0 * v * a * (1.0 - w)), e = (-2 * a * v);

  		if ((emw >= em1) || (e >= em1)) return(-w);
  		tt = M_LN2 - log1pem1(emw);
  		double temp = logdiff(emw, e) - log1pem1(e);
  		if (log(w) > temp) {
  			tt += logdiff(log(w), temp);
  			tt = -exp(tt);
  		}
  		else {
  			tt += logdiff(temp, log(w));
  			tt = exp(tt);
  		}
  	}
  	if (gsl_finite(tt)) return(tt); else
  	{
  //		std::cout << "davlogprob " << setw(20) << a << setw(20) << v << setw(20) << w << std::endl;
  		return(GSL_NEGINF);
  	}

  	return tt;
  }

  //derivative by logprob_upperbound by a
  double dalogprob_upperbound(int pm, double a, double v, double w, double dav) {

  	if (fabs(v) == 0.0) return 0.0;
  	double tt;
  	tt = dav * v;
  	tt = (pm == 1) ? -tt : tt;
  	if (gsl_finite(tt)) return(tt); else
  	{
  //		std::cout << "dalogprob " << setw(20) << a << setw(20) << v << setw(20) << w << std::endl;
  		return(GSL_NEGINF);
  	}
  	return tt;
  }

  //derivative by logprob_upperbound by v
  double dvlogprob_upperbound(int pm, double a, double v, double w, double dav) {

  	//	if (fabs(v) == 0.0) return 0.0;
  	int sign = 1; if (pm == 1) sign = -1;
  	double tt = dav * a * sign;
  	if (gsl_finite(tt)) return(tt); else
  	{
  //		std::cout << "dvlogprob " << setw(20) << a << setw(20) << v << setw(20) << w << std::endl;
  		return(GSL_NEGINF);
  	}
  	return tt;
  }

  //derivative by logprob_upperbound by w
  double dwlogprob_upperbound(int pm, double a, double v, double w) {

  	double tt = 1.0;
  	if (pm == 1) {
  		w = 1.0 - w;
  		v = -v;
  		tt = -1.0;
  	}

  	if (fabs(v) == 0.0) return -tt / (1.0 - w);

  	if (v < 0) {
  		double e = (2.0 * v * a * (1.0 - w));
  		double temp = M_LN2 + e + log(fabs(v)) + log(a) - log1pem1(e);
  		tt *= -exp(temp);
  	}
  	else
  	{
  		double e = -(2.0 * v * a * (1.0 - w));
  		double temp = M_LN2 + log(fabs(v)) + log(a) - log1pem1(e);
  		tt *= -exp(temp);
  	}
  	return tt;
  }

}
