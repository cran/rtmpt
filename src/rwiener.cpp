#include "rts.h"

namespace drtmpt {
  
  //upper hull in adaptive rejection sampling
  double fun_upper(int k, double x, const std::vector<piece> &upper) {
    int i = 1;
    //	int k = static_cast<int>(upper.size());
    while ((i != k) && (x >= upper[i].z)) i++;
    i = i - 1;
    double t = upper[i].absc + upper[i].slope * (x - upper[i].center);
    return t;
  }
  
  //generate approximations in adaptive rejection sampling
  void generate_intervals(int& k, double totallow, const std::vector<point> &h, std::vector<piece>& lower, std::vector<piece>& upper, std::vector<double>& s) {
    k = static_cast<int>(h.size());
    
    lower.clear(); upper.clear(); piece low, up;
    for (int j = 0; j != k; j++) {
      double z;
      if (j == 0) z = totallow;
      else z = (h[j].h - h[j - 1].h - h[j].x * h[j].dh + h[j - 1].x * h[j - 1].dh) / (h[j - 1].dh - h[j].dh);
      
      up.z = z;
      up.slope = h[j].dh; up.absc = h[j].h; up.center = h[j].x;
      upper.push_back(up);
      if (j == 0) low.z = totallow;
      else low.z = h[j - 1].x;
      lower.push_back(low);
    }
    low.z = h[k - 1].x; lower.push_back(low);
    
    double sum = GSL_NEGINF, t; s.clear();
    for (int i = 0; i != k; i++) {
      if (i == 0) t = fun_upper(k, upper[i + 1].z, upper);
      else
        if (i < k - 1) {
          double sl = upper[i].slope;
          t = upper[i].absc - upper[i].center * sl + ( (sl > 0)? logdiff(upper[i + 1].z * sl,
                                                        upper[i].z * sl) : logdiff(upper[i].z * sl, upper[i + 1].z * sl));
          
        }
        else {
          t = (fun_upper(k, upper[i].z, upper));
        }
        t -= log(fabs(upper[i].slope));
        
        sum = logsum(sum, t);
        s.push_back(sum);
    }
  }
  
  //update approximations in adaptive rejection sampling if newpoint (= xstar) is added
  bool update_intervals(int k, double totallow, point new_point, std::vector<point>& h, std::vector<piece>& lower, std::vector<piece>& upper, std::vector<double>& s) {
    double x = new_point.x;
    bool flag = false;
    int i = 0;
    k = static_cast<int>(h.size());
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
    else z = (h[j].h - h[j - 1].h - h[j].x * h[j].dh + h[j - 1].x * h[j - 1].dh) / (h[j - 1].dh - h[j].dh);
    up.z = z;
    
    up.slope = h[j].dh; up.absc = h[j].h; up.center = h[j].x;
    if (i < k) upper[i] = up; else upper.push_back(up);
    
    
    if (i < k) {
      j = i + 1;
      z = (h[j].h - h[j - 1].h - h[j].x * h[j].dh + h[j - 1].x * h[j - 1].dh) / (h[j - 1].dh - h[j].dh);
      up.z = z;
      up.slope = h[j].dh; up.absc = h[j].h; up.center = h[j].x;
      upper.insert(upper.begin() + j, up);
    }
    
    k = k + 1;
    
    double sum = 0, t;
    std::vector<double> sold = s;
    
    {
      if (i > 1) sum = sold[i - 2];
      int iup = (i + 1 == k) ? i + 1 : i + 2;
      int ilow = (i == 0) ? 0 : i - 1;
      for (int j = ilow; j != iup; j++) {
        if (j == 0) t = fun_upper(k, upper[j + 1].z, upper);
        else
          if (j < k - 1) {
            double sl = upper[j].slope;
            t = upper[j].absc - upper[j].center * sl + ((sl > 0)? logdiff(upper[j + 1].z * sl,
                                                         upper[j].z * sl) : logdiff(upper[j].z * sl, upper[j + 1].z * sl));
          }
          else {
            t = (fun_upper(k, upper[j].z, upper));
          }
          t -= log(fabs(upper[j].slope));
          if (j == 0) sum = t;
          else sum = logsum(sum, t);
          if (j != i) s[j] = sum; else s.insert(s.begin() + j, sum);
      }
    }
    if (i + 1 < k) {
      if (sold[i] < s[i + 1]) {
        //			std::cout << "Denkfehler " << setw(20) << sold[i] << setw(20) << s[i + 1] << std::endl;
        flag = true;
      }
      double temp = logdiff(sold[i], s[i + 1]);
      for (int j = i + 2; j < k; j++) {
        s[j] = logdiff(sold[j - 1], temp);
      }
    }
    return flag;
  }
  
  
  //lower hull in adpative rejection sampling
  double fun_lower(int k, double x, const std::vector<point> &h, const std::vector<piece> &lower) {
    int i = 1;
    k = k + 1;
    while ((i != k) && (x >= lower[i].z)) i++;
    i = i - 1; double t;
    if ((i == 0) || (i == k - 1)) t = GSL_NEGINF;
    else t = ((h[i].x - x) * h[i - 1].h + (x - h[i - 1].x) * h[i].h) / (h[i].x - h[i - 1].x);
    
    return t;
  }
  
  //invert linear approximation of density in adaptive rejection sampling
  double inverse_distribution(int k, double xstar, const std::vector<piece> &upper, std::vector<double> s, double bound, bool& flag) {
    double sum = 0, t;
    
    if (bound == GSL_POSINF) sum = s[k - 1];
    else {
      if (bound <= upper[k - 1].z) {
        //			std::cout << "Problem in inverse" << setw(20) << k << setw(20) << bound << setw(20) << upper[k-1].z;
        flag = true;
      }
      double sl = upper[k - 1].slope;
      t = upper[k - 1].absc - upper[k - 1].center * sl + ((sl > 0) ? logdiff(bound * sl, upper[k - 1].z * sl) : logdiff(upper[k - 1].z * sl, bound * sl));
      t -= log(fabs(sl));
      s[k - 1] = logsum(t, s[k - 2]);
      sum = s[k - 1];
    }
    int j = 0;
    double temp = log(xstar) + sum;
    while (temp > s[j]) j++;
    //	if (j > k - 1) { std::cout << "Wie das?"; char x; std::cin >> x; }
    
    
    double sl = upper[j].slope;
    double help = log(fabs(sl)); int sign = sl > 0 ? 1 : -1;
    if (isnan(sl)) {
      flag = true;
    }
    
    if (j > 0) temp = logdiff(temp, s[j - 1]);
    help = help + temp - upper[j].absc + upper[j].center * sl;
    if (sign == 1) temp = logsum(help, upper[j].z * sl);
    else temp = logdiff(upper[j].z * sl, help);
    t = temp / sl;
    
    if (t < upper[j].z) {
      //		std::cout << std::endl;
      //		std::cout << "nanu " << j << " " << k - 1 << " " << t << " " << upper[j].z << " " << upper[j + 1].z << " " << s[j - 1]
      //			<< " " << upper[j].slope << " " << upper[j].absc << " " << temp << " "
      //			<< fun_upper(k, upper[j].z, upper) << " " << fun_upper(k, upper[j + 1].z, upper) << std::endl;
      //		std::cout << std::endl;
      t = upper[j].z;
      flag = true;
    }
    return t;
  }
  
  //transformed Wiener density for adaptive rejection sampling
  void wiener_comp(double start, double scale, double norm, double alpha, double a, double v, double w, point& h) {
    h.x = alpha;
    double t = exp(alpha * scale + start);
    
    h.h = dwiener_d(-t, a, v, w, accuracy*1.2);
    h.dh = dtdwiener_d(t, a, v, w, h.h);
    h.h += start + h.x * scale + log(scale) - norm;
    h.dh = (1.0 + t * h.dh) * scale;
  }
  
  //define order
  bool compare(point h1, point h2) {
    return (h1.x < h2.x);
  }
  
  //init adaptive rejection sampling
  void initialize_ars(int t, double* tavw, ars_archiv& ars_store) {
    
    for (int m = 0; m != no_patterns; m++) {
      double a = dTAVW(t, 0, dCOMB(m, 0));
      double v = dTAVW(t, 1, dCOMB(m, 1));
      double w = dTAVW(t, 2, dCOMB(m, 2));
      
      for (int pm = 0; pm != 2; pm++) {
        double norm = 0.0, step = 1.0;
        if (pm == 1) {
          v = -v;
          w = 1 - w;
        }
        
        double t0 = exp_mean(0, a, v, w);
        double start = log(t0);
        
        double scale;
        if (fabs(v) > 0.01)
          scale = sqrt(lower_bound_var(a, v, w)) / 1.0;
        else {
          int sign = (v > 0) ? 1 : -1;
          scale = sqrt(lower_bound_var(a, sign * 0.01, w)) / 1.0;
        }
        
        //			scale = fmax(scale, 0.05);
        scale = log(t0 + scale) - start;
        
        point begun, one;
        double ub, lb;
        point high, low;
        begun.x = 0.0;
        wiener_comp(start, scale, norm, begun.x, a, v, w, begun);
        //			if (fabs(begun.dh) > 10) { std::cout << "begun " << begun.dh << setw(20) << a << setw(20) << v << setw(20) << w << std::endl; }
        norm = begun.h; begun.h = 0.0;
        int sign = begun.dh > 0 ? 1 : -1;
        // make ldh >2.0 <5.0
        std::vector<point> h;
        for (int i = 0; i != 2; i++) {
          one = begun;
          double dh = sign * one.dh;
          if ((dh <= 2.0) || (dh >= 5.0))
          {
            if (dh <= 2.0) {
              while (dh <= 2.0) {
                double dold = dh;
                lb = one.x;
                one.x -= sign * step;
                
                wiener_comp(start, scale, norm, one.x, a, v, w, one);
                if ((abs(one.dh) > 2.0) && (abs(one.dh) < 5)) h.push_back(one);
                dh = sign * one.dh;
                if (dold > dh) {
                  //							std::cout << "convex2?";
                  //char xx; std::cin >> xx;
                }
              }
              ub = one.x;
            }
            else {
              while (dh >= 5.0) {
                double dold = dh;
                ub = one.x;
                one.x += sign * step;
                wiener_comp(start, scale, norm, one.x, a, v, w, one);
                if ((abs(one.dh) > 2.0) && (abs(one.dh) < 5)) h.push_back(one);
                dh = sign * one.dh;
                if (dh > dold) {
                  //							std::cout << "convex2?";
                }
              }
              lb = one.x;
            }
          }
          while ((dh <= 2.0) || (dh >= 5.0)) {
            one.x = (lb + ub) / 2.0;
            wiener_comp(start, scale, norm, one.x, a, v, w, one);
            if (abs(one.dh) < 5) h.push_back(one);
            dh = sign * one.dh;
            if (dh <= 2.0) { lb = one.x; }
            if (dh >= 5.0) { ub = one.x; }
          }
          if (sign == 1) low = one; else high = one;
          sign = -sign;
        }
        
        wiener_comp(start, scale, norm, (low.x + high.x) / 2, a, v, w, one);
        
        h.push_back(one); h.push_back(low); h.push_back(high);
        
        double bound = (log(rtmins[t * no_patterns * 2 + m * 2 + pm]) - start) / scale;
        if (high.x > bound) {
          one.x = bound - step;
          wiener_comp(start, scale, norm, one.x, a, v, w, one);
          h.push_back(one);
        }
        if (low.x > bound) {
          one.x = bound - 2 * step;
          wiener_comp(start, scale, norm, one.x, a, v, w, one);
          h.push_back(one);
        }
        
        sort(h.begin(), h.end(), compare);
        
        std::vector<point> xemp; xemp.clear();
        xemp.push_back(h[0]);
        int hs = h.size();
        for (int i = 1; i != hs; i++) if (h[i].x > h[i - 1].x) xemp.push_back(h[i]);
        h = xemp;
        
        int k = 0;
        std::vector<piece> lower, upper; std::vector<double> s;
        generate_intervals(k, GSL_NEGINF, h, lower, upper, s);
        
        ars_store.startstore.push_back(start);
        ars_store.scalestore.push_back(scale);
        ars_store.normstore.push_back(norm);
        ars_store.hstore.push_back(h);
        ars_store.lowerstore.push_back(lower);
        ars_store.upperstore.push_back(upper);
        ars_store.sstore.push_back(s);
      }
    }
  }
  
  //adaptive rejection sampling
  double arst(int t, int m, int pm, ars_archiv& ars_store, double scale, double totallow, double start, double bound, double a, double v, double w, gsl_rng* rst,
              void generic2(double start, double scale, double norm, double alpha, double a, double v, double w, point& h)) {
    int index = t * no_patterns * 2 + 2 * m + pm;
    double norm = ars_store.normstore[index]; bool flag;
    
    flag = false;
    std::vector<point> h; std::vector<piece> lower, upper; std::vector<double> s;
    double ww, tt, xstar;
    point one;
    
    h = ars_store.hstore[index];
    lower = ars_store.lowerstore[index];
    upper = ars_store.upperstore[index];
    s = ars_store.sstore[index];
    
    int k = static_cast<int>(h.size());
    bool update = false;
    
    if (bound < GSL_POSINF)
    {
      int l = 0;
      while ((l != k) && (h[l].x < bound)) l++;
      k = l;
    }
    
    double ss;
    
    WEITER:
      //	MONITOR(0, 2)++;
      xstar = oneuni(rst);
    
    xstar = inverse_distribution(k, xstar, upper, s, bound, flag);
    if (flag) {
      xstar = GSL_NEGINF;
      goto END;
    }
    ww = log(oneuni(rst)); tt = fun_upper(k, xstar, upper);  ss = fun_lower(k, xstar, h, lower);
    if (ww <= (ss - tt))  goto STOP;
    one.x = xstar; generic2(start, scale, norm, xstar, a, v, w, one);
    flag = update_intervals(k, totallow, one, h, lower, upper, s);
    if (flag) {
      xstar = GSL_NEGINF;
      goto END;
    }
    k = k + 1;
    update = true;
    if (ww <= (one.h - tt))  goto STOP;
    //	MONITOR(1, 2)++;
    goto WEITER;
    STOP:
      if (update) {
        ars_store.hstore[index] = h;
        ars_store.lowerstore[index] = lower;
        ars_store.upperstore[index] = upper;
        ars_store.sstore[index] = s;
      }
      END:	return xstar;
      
  }
  
  //generate random Wiener fisrt passage time at bound pm based on inverse cdf method
  double rwiener_diag(int pm, double bound, double a, double v, double w, gsl_rng* rst) {
    double eps = 1e-5;
    double pmid = 0;
    double qmin = 0;
    double qmax = bound;
    double q = gsl_isinf(bound) ? 1.0 : bound / 2;
    double p = log(oneuni(rst));
    
    double qold;
    
    if (pm == 1) {
      v = -v;
      w = 1.0 - w;
    }
    double total = (gsl_isinf(bound)) ? logprob_upperbound(0, a, v, w) : logFjoint_lower(bound, a, v, w);
    
    do {
      qold = q;
      pmid = logFjoint_lower(q, a, v, w) - total;
      if (p <= pmid) {
        qmax = q;
        q = qmin + (qmax - qmin) / 2.0;
      }
      else {
        qmin = q;
        q = gsl_isinf(qmax) ? q * 2 : qmin + (qmax - qmin) / 2.0;
      }
      
    } while (fabs(q - qold) > eps);
    return(q);
  }
  
  //generate wiener variate using adaptive rejection sampling with fallback based on inverting cumulative distribution function
  double make_rwiener(int t, int m, int pm, ars_archiv& ars_store, double bound, double a, double v, double w, gsl_rng* rst) {
    double temp;
    int index = t * no_patterns * 2 + 2 * m + pm;
    double start = ars_store.startstore[index];
    double  scale = ars_store.scalestore[index];
    double bound2 = (bound == GSL_POSINF) ? bound : (log(bound) - start) / scale;
    if (pm == 1) { v = -v; w = 1 - w; }
    temp = arst(t, m, pm, ars_store, scale, GSL_NEGINF, start, bound2, a, v, w, rst, wiener_comp);
    if (temp != GSL_NEGINF) temp = exp(start + temp * scale);
    else {
      //		std::cout << "arts hat nicht geklappt";
      temp = rwiener_diag(0, bound, a, v, w, rst);
    }
    return temp;
  }

}
