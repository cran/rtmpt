//implements adaptive rejection sampling
// authors: Christoph Klauer

#include "rts.h"
#include <mutex>

std::mutex mtx_R_CUI;

// namespace rtsNS {

	struct piece {
		double z;
		double slope;
		double absc;
		double center;
	};




	void generate_intervals(double totallow, vector<point> h, vector<piece> &lower, vector<piece> &upper) {
		int k = static_cast<int>(h.size());

		lower.clear(); upper.clear(); piece low, up;
		for (int j = 0; j != k; j++) {
			double z;
			if (j == 0) z = totallow;
			else z = (h[j].h - h[j - 1].h - h[j].x*h[j].dh + h[j - 1].x*h[j - 1].dh) / (h[j - 1].dh - h[j].dh);
			up.z = z;
			up.slope = h[j].dh; up.absc = h[j].h; up.center = h[j].x;
			upper.push_back(up);
			if (j == 0) low.z = totallow;
			else low.z = h[j - 1].x;
			lower.push_back(low);
		}
		low.z = h[k - 1].x; lower.push_back(low);
	}

	void update_intervals(double totallow, point new_point, vector<point> &h, vector<piece> &lower, vector<piece> &upper) {
		vector<point> temp; temp.clear();
		double x = new_point.x;
		int i = 0; int k = static_cast<int>(h.size());
		while ((i != k) && (x > h[i].x))  i++;
		h.insert(h.begin() + i, new_point);
		piece low;
		int j = i + 1;
		low.z = h[i].x;
		lower.insert(lower.begin() + j, low);
		j = i;
		piece up;
		double z;
		if (j == 0) z = totallow;
		else z = (h[j].h - h[j - 1].h - h[j].x*h[j].dh + h[j - 1].x*h[j - 1].dh) / (h[j - 1].dh - h[j].dh);
		up.z = z;
		up.slope = h[j].dh; up.absc = h[j].h; up.center = h[j].x;
		if (i < k) upper[i] = up; else upper.push_back(up);
		if (i < k) {
			j = i + 1;
			z = (h[j].h - h[j - 1].h - h[j].x*h[j].dh + h[j - 1].x*h[j - 1].dh) / (h[j - 1].dh - h[j].dh);
			up.z = z;
			up.slope = h[j].dh; up.absc = h[j].h; up.center = h[j].x;
			upper.insert(upper.begin() + j, up);
		}
	}

	double fun_upper(double x, vector<piece> upper) {
		int i = 1; int k = static_cast<int>(upper.size());
		while ((i != k) && (x >= upper[i].z)) i++;
		i = i - 1;
		double t = upper[i].absc + upper[i].slope * (x - upper[i].center);
		return t;
	}

	double fun_lower(double x, vector<point> h, vector<piece> lower) {
		int i = 1; int k = static_cast<int>(lower.size());
		while ((i != k) && (x >= lower[i].z)) i++;
		i = i - 1; double t;
		if ((i == 0) || (i == k - 1)) t = -DBL_MAX;
		else t = ((h[i].x - x)*h[i - 1].h + (x - h[i - 1].x)*h[i].h) / (h[i].x - h[i - 1].x);

		return t;
	}

	double elogdiff(double xa, double xb) {
		double result;
		if (xa > xb) result = rexp(xa + gsl_log1p(-exp(xb - xa)));
		else
			if (xb > xa) result = -rexp(xb + gsl_log1p(-exp(xa - xb)));
			else result = 0.0;
		return result;
	}

	double logdiff(double xa, double xb) {
		double result;
		if (xa > xb) result = (xa + gsl_log1p(-exp(xb - xa)));
		else
			if (xb > xa) result = (xb + gsl_log1p(-exp(xa - xb)));
			else result = -DBL_MAX;
		return result;
	}

	double inverse_distribution(double xstar, vector<piece> upper, bool &flag) {
		double sum = 0, t; vector<double> s;
		int k = static_cast<int>(upper.size());
		for (int i = 0; i != k; i++) {
			if (i == 0) t = fun_upper(upper[i + 1].z, upper);
			else
				if (i < k - 1) {
					double sl = upper[i].slope;
					t = upper[i].absc - upper[i].center *sl + logdiff(upper[i + 1].z * sl,
						upper[i].z*sl);
				}
				else t = (fun_upper(upper[i].z, upper));
			t -= log(fabs(upper[i].slope));
			//		if (t<0) std::cout<< "Problem";
			if (i == 0) sum = t;
			else sum = logsum(sum, t);
			s.push_back(sum);
		}

		int j = 0;
		double temp = log(xstar) + sum;
		while (temp > s[j]) j++;
		if (j > k - 1) { Rprintf("Wie das?"); }

		double sl = upper[j].slope;
		double help = log(fabs(sl)); int sign = sl > 0 ? 1 : -1;

		if (j > 0) temp = logdiff(temp, s[j - 1]);
		help = help + temp - upper[j].absc + upper[j].center*sl;
		if (sign == 1) temp = logsum(help, upper[j].z*sl);
		else temp = logdiff(upper[j].z*sl, help);
		t = temp / sl;

		if (t < upper[j].z) {
			/*		std::cout << std::endl;
					std::cout<<"nanu " << j << " " << k-1 << " " << t << " " << upper[j].z << " " << upper[j+1].z << " " << s[j-1]
						<< " " << upper[j].slope << " " << upper[j].absc << " " << temp << " "
						<< fun_upper(upper[j].z,upper) << " " << fun_upper(upper[j+1].z,upper)  << std::endl;
					std::cout<<std::endl;// else if ((j+1<k) && (t>upper[j+1].z)) std::cout << "nanu2";
			*/	//char x; std::cin >> x;
		t = upper[j].z;
		flag = true;
		}

		return t;
	}


	double ars(double step, double &scale, double totallow, double n, double xp, double *beta, double *sigi, double *lambdas, double *lams, int tt, int iz, double start, gsl_rng *rst, void gamma_prior(double scale, double norm, double n, double alpha, double p, double *beta, double *sigi, double *lambdas, double *lams, int t, int iz, bool deriv, point &h)) {
		// intialize h's
		double norm = 0.0;
	NEW:
		bool flag = false;
		vector<point> h; vector<piece> lower, upper; h.clear(); lower.clear(); upper.clear();
		double w, t, xstar;
		point one;
		one.x = start;
		gamma_prior(scale, norm, n, one.x, xp, beta, sigi, lambdas, lams, tt, iz, false, one); norm = one.h + norm;

		gamma_prior(scale, norm, n, one.x, xp, beta, sigi, lambdas, lams, tt, iz, true, one);

		double ub, lb;
		point high, low;
		int sign = one.dh > 0 ? 1 : -1;
		// make ldh >2.0 <5.0
		for (int i = 0; i != 2; i++) {
			int cnt_while = 0;
			double dh = sign * one.dh;
			if ((dh <= 2.0) || (dh >= 5.0))
			{
				if (dh <= 2.0) {
					lb = one.x;
					while (dh <= 2.0) {
						one.x -= sign * step;
						gamma_prior(scale, norm, n, one.x, xp, beta, sigi, lambdas, lams, tt, iz, true, one);
						dh = sign * one.dh;
						cnt_while++; if (cnt_while % 1024 == 0) {
						  std::lock_guard<std::mutex> guard(mtx_R_CUI);
						  R_CheckUserInterrupt();
						}
					}
					ub = one.x;
				}
				else {
					if (dh >= 5.0) {
						ub = one.x;
						while (dh >= 5.0) {
							one.x += sign * step;
							gamma_prior(scale, norm, n, one.x, xp, beta, sigi, lambdas, lams, tt, iz, true, one);
							dh = sign * one.dh;
							cnt_while++; if (cnt_while % 1024 == 0) {
						  std::lock_guard<std::mutex> guard(mtx_R_CUI);
						  R_CheckUserInterrupt();
						}
						}
						lb = one.x;
					}
				}
				while ((dh <= 2.0) || (dh >= 5.0)) {
					one.x = (lb + ub) / 2.0;
					gamma_prior(scale, norm, n, one.x, xp, beta, sigi, lambdas, lams, tt, iz, true, one);
					dh = sign * one.dh;
					if (dh <= 2.0) { lb = one.x; }
					if (dh >= 5.0) { ub = one.x; }
					cnt_while++; if (cnt_while % 1024 == 0) {
						  std::lock_guard<std::mutex> guard(mtx_R_CUI);
						  R_CheckUserInterrupt();
						}
				}
			}
			if (sign == 1) low.x = one.x; else high.x = one.x;
			sign = -sign;
		}
		gamma_prior(scale, norm, n, (low.x + high.x) / 2, xp, beta, sigi, lambdas, lams, tt, iz, false, one);
		norm = one.h + norm;
		gamma_prior(scale, norm, n, low.x, xp, beta, sigi, lambdas, lams, tt, iz, true, one); gamma_prior(scale, norm, n, low.x, xp, beta, sigi, lambdas, lams, tt, iz, false, one); h.push_back(one);
		gamma_prior(scale, norm, n, high.x, xp, beta, sigi, lambdas, lams, tt, iz, true, one); gamma_prior(scale, norm, n, high.x, xp, beta, sigi, lambdas, lams, tt, iz, false, one); h.push_back(one);


		generate_intervals(totallow, h, lower, upper);
	WEITER:	xstar = oneuni(rst);
		xstar = inverse_distribution(xstar, upper, flag);
		if (flag) {
			scale = scale * 10; start /= 10; Rprintf("NEW0 in ars"); goto NEW;
		}
		w = log(oneuni(rst)); t = fun_upper(xstar, upper); double s = fun_lower(xstar, h, lower);
		if (w <= (s - t))  goto STOP;
		one.x = xstar; gamma_prior(scale, norm, n, xstar, xp, beta, sigi, lambdas, lams, tt, iz, false, one);
		if (w <= (one.h - t))  goto STOP;

		gamma_prior(scale, norm, n, xstar, xp, beta, sigi, lambdas, lams, tt, iz, true, one);

		update_intervals(totallow, one, h, lower, upper);
		goto WEITER;
	STOP:

		return xstar;

	}

// }
