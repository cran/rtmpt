#include "rts.h"
#include <atomic>

namespace drtmpt {

  std::atomic<int> curr_order(0);
  //one cycle of Gibbs-Hamiltonian algorithm
  void gibbs_full_cycle(bool& change, ars_archiv& ars_store, const std::vector<trial> & daten, int* nips, gsl_vector* hampar, double* tavw, double* tlams, double* ai, double* sig,
                        double* loglambdas, double* bi, double* alltaus, double* rest, double* gam, double& omega, int* paths,
                        double liknorm[6], double& activeeps, double& epsm, double& Hobjective, int m, bool save, gsl_rng* rst) {


    double* sigi = 0; if (!(sigi = (double*)malloc(icompg * icompg * sizeof(double)))) { Rprintf("Allocation failure2\n"); }
    double* gami = 0; if (!(gami = (double*)malloc(respno * respno * sizeof(double)))) { Rprintf("Allocation failure2\n"); }


    //#pragma omp atomic
    //	MONITOR(0, 0) += 1;

    if (change) {
      ars_store.hstore.clear(); ars_store.lowerstore.clear(); ars_store.startstore.clear(); ars_store.upperstore.clear(); ars_store.scalestore.clear(); ars_store.normstore.clear(); ars_store.sstore.clear();
      for (int t = 0; t != indi; t++) initialize_ars(t, tavw, ars_store);
    }

    for (int x = 0; x != datenzahl; x++) {
      int old_path = paths[x];
      make_path(daten[x], nips, x, paths[x], hampar, tavw, tlams, loglambdas, alltaus, rest, ars_store, rst);
      //#pragma omp atomic
      //		MONITOR(0, 1)++;
      //		if (paths[x] != old_path) MONITOR(1, 1)++;
    }
    if (phase <= 2) {
      gsl_matrix* Ltminusx = gsl_matrix_alloc(icompg, icompg);
      gsl_matrix* Ltminusr = gsl_matrix_alloc(respno, respno);
      sample_sig(hampar, sig, sigi, Ltminusx, ai, rst);
      make_rgam(hampar, gam, gami, Ltminusr, bi, rst);
      make_romega(hampar, loglambdas, omega, rst);

      double* scale = 0; if (!(scale = (double*)malloc(icompg * igroup * sizeof(double)))) { Rprintf("Allocation failure2\n"); }
      for (int ig = 0; ig != igroup; ig++) for (int i = 0; i != icompg; i++) scale[ig * icompg + i] = sqrt(dSIG(i, i) / ng[ig]);
      double* rscale = 0; if (!(rscale = (double*)malloc(respno * igroup * sizeof(double)))) { Rprintf("Allocation failure2\n"); }
      for (int ig = 0; ig != igroup; ig++) for (int i = 0; i != respno; i++) rscale[ig * respno + i] = sqrt(dGAM(i, i) / ng[ig]);
      double* sl = 0; if (!(sl = (double*)malloc(indi * sizeof(double)))) { Rprintf("Allocation failure2\n"); }
      for (int t = 0; t != indi; t++) {
        double nt = n_per_subj[t];
        nt = nt / (nt - 2.0);
        sl[t] = sqrt(nt * omega);
      }
      change = hnuts(scale, nips, hampar, tavw, tlams, sig, sigi, Ltminusx, daten, rscale, sl, rest, loglambdas, gam, gami, Ltminusr, omega, alltaus, liknorm[0], liknorm[1], activeeps, epsm, Hobjective, m, rst);
      free(scale); free(rscale); free(sl);
      gsl_matrix_free(Ltminusx); gsl_matrix_free(Ltminusr);
    }
    else
      change = hnuts2(nips, hampar, tavw, tlams, daten, rest, loglambdas, alltaus, liknorm, activeeps, epsm, Hobjective, m, save, rst);

    // if (!change)
      //#pragma omp atomic
      //		MONITOR(1, 0) += 1;


      if (sigi) free(sigi);
      if (gami) free(gami);

  }


  //ireps cycles; save into sample if phase 4; update for R statstic; update for posterior variance/covariance structure in phases 2 and 3
  void gibbs_and_monitor(const std::vector<trial> & daten, int* nips, gsl_vector* hampar, double* tavw, double* tlams, double* ai, double* loglambdas, double* bi, double* alltaus, double* rest, int* paths,  double liknorm[6],
                         double& activeeps, double& epsm, double& Hobjective,
                         int offset, int n_all_parameters, double* parmon, gsl_rng* rst, int ithread,
                         bool save, double* sample) {
    double* sig = 0;	if (!(sig = (double*)malloc(icompg * icompg * sizeof(double)))) { Rprintf("Allocation failure\n"); }
    double* gam = 0;	if (!(gam = (double*)malloc(respno * respno * sizeof(double)))) { Rprintf("Allocation failure\n"); }
    double* temp = 0; if (!(temp = (double*)malloc(n_all_parameters * sizeof(double)))) { Rprintf("Allocation failure\n"); }
    ars_archiv ars_store;
    double omega;

    bool change = true;

    //reset parmon if phase change
    if (offset == 0) {
      gsl_vector_view tx = gsl_vector_view_array(parmon, 2 * n_all_parameters);
      gsl_vector_set_zero(&tx.vector);
    }
    
    //compute ireps cycles
    for (int i = 0; i != ireps; i++) {
      gibbs_full_cycle(change, ars_store, daten, nips, hampar, tavw, tlams, ai, sig, loglambdas, bi, alltaus, rest, gam, omega, paths, liknorm, activeeps, epsm, Hobjective, offset + i + 1, save, rst);
      
      //save parameters in temp in the phase 3 and 4 parameterization
      gsl_vector_view t1 = gsl_vector_view_array(temp, n_all_parameters);
      
      if (phase >= 3) {
        gsl_vector_memcpy(&t1.vector, hampar);
      } else {// do we need this in phase == 1 ?
        gsl_vector_view t2 = gsl_vector_subvector(&t1.vector, 0, nhamil - indi);
        gsl_vector_view t3 = gsl_vector_subvector(hampar, 0, nhamil - indi);
        gsl_vector_memcpy(&t2.vector, &t3.vector);
        // in these procedures we also move from theta'' to theta''' (see paper)
        make_hampar_from_sig(0, sig, &t1.vector);
        make_hampar_from_sig(1, gam, &t1.vector);
        gsl_vector_set(&t1.vector, n_all_parameters - 1, log(omega));
        

        gsl_vector_view t4 = gsl_vector_subvector(&t1.vector, nhamil - indi, indi);
        gsl_vector_view t5 = gsl_vector_view_array(loglambdas, indi);
        gsl_vector_memcpy(&t4.vector, &t5.vector);
      }
      
      // n_all_parameters = icompg*igroup + indi*icompg + (icompg*(icompg+1))/2 + respno*igroup + (respno+1)*indi + (respno*(respno+1))/2 +1
      //                    ma,mv,mw  a,v,w            sig                   rmu     lambdas+sig_t          gam            omega

      // save in double* sample if seriously sampling
      if ((save) && (i % THIN == 0)) {
        // at this point phase must equal 3 or larger
        int off = (ithread * IREP + i) * (n_all_parameters);
        gsl_vector_view ts = gsl_vector_view_array(sample, NOTHREADS * IREP * (n_all_parameters));
        gsl_vector_view ts2 = gsl_vector_subvector(&ts.vector, off, n_all_parameters);

        gsl_vector* hamp_temp = gsl_vector_alloc(n_all_parameters);
        gsl_vector_memcpy(hamp_temp, hampar);

        for (int flag = 0; flag != 2; flag++) {
          int k, ioff, offsig;
          if (flag == 0) {
            k = icompg;
            ioff = iavwoff;
            offsig = nhamil + (icompg * (icompg - 1)) / 2;
          }
          else {
            k = respno;
            ioff = ilamoff;
            offsig = nhamil + (icompg * (icompg + 1)) / 2 + (respno * (respno - 1)) / 2;
          }
          gsl_vector* wsigs = gsl_vector_alloc(k);
          for (int i = 0; i != k; i++) gsl_vector_set(wsigs, i, exp(gsl_vector_get(hampar, offsig + i)));

          gsl_vector_view thetv = gsl_vector_subvector(hamp_temp, ioff, indi * k);
          gsl_matrix_view thetm = gsl_matrix_view_vector(&thetv.vector, indi, k);

          gsl_matrix* sd = gsl_matrix_calloc(k, k);
          gsl_vector_view sddiag = gsl_matrix_diagonal(sd);
          gsl_vector_memcpy(&sddiag.vector, wsigs);
          std::vector<double> z;
          from_y_to_z(flag, hamp_temp, z);
          gsl_matrix* w = gsl_matrix_alloc(k, k);
          from_z_to_w(flag, z, w);
          gsl_blas_dtrmm(CblasRight, CblasLower, CblasNoTrans, CblasNonUnit, 1.0, w, sd);
          gsl_blas_dtrmm(CblasRight, CblasLower, CblasTrans, CblasNonUnit, 1.0, sd, &thetm.matrix);
          gsl_matrix_view t1 = gsl_matrix_view_array((flag == 0) ? sig : gam, k, k);
          gsl_matrix_transpose_memcpy(&t1.matrix, sd);
          gsl_blas_dtrmm(CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, 1.0, sd, &t1.matrix);
          gsl_matrix_free(w);
          gsl_matrix_free(sd);
          gsl_vector_free(wsigs);
        }
        int jj = nhamil;
        for (int ip = 0; ip != icompg; ip++) for (int jp = ip; jp != icompg; jp++) gsl_vector_set(hamp_temp, jj++, dSIG(ip, jp));
        for (int ip = 0; ip != respno; ip++) for (int jp = ip; jp != respno; jp++) gsl_vector_set(hamp_temp, jj++, dGAM(ip, jp));
        gsl_vector_set(hamp_temp, jj, exp(gsl_vector_get(hampar, jj)));
        gsl_vector_memcpy(&ts2.vector, hamp_temp);
        gsl_vector_free(hamp_temp);
      }
      
      //update r statistics
      int m = i + offset + 1;
      double r = 1.0 / m;


      gsl_vector_view vparm = gsl_vector_view_array(parmon, 2 * n_all_parameters);
      gsl_vector_view vparm1 = gsl_vector_subvector(&vparm.vector, 0, n_all_parameters);
      gsl_vector_view vparm2 = gsl_vector_subvector(&vparm.vector, n_all_parameters, n_all_parameters);
      gsl_vector* dev1 = gsl_vector_alloc(n_all_parameters);
      

      gsl_vector_memcpy(dev1, &t1.vector);
      gsl_blas_daxpy(-1.0, &vparm1.vector, dev1);
      gsl_vector* help = gsl_vector_alloc(n_all_parameters);
      gsl_vector_memcpy(help, dev1);
      gsl_vector_mul(help, help);
      gsl_blas_daxpy(1.0 - r, help, &vparm2.vector);
      gsl_blas_daxpy(r, dev1, &vparm1.vector);
      gsl_vector_free(help);
      
      //update posterior variance/covariance matrix in phase 2 and 3
      if ((phase >= 2) && (phase < 4)) {

        gsl_vector* dev2 = gsl_vector_alloc(n_all_parameters);
        gsl_vector_memcpy(dev2, &t1.vector);
        gsl_blas_daxpy(-1.0, &vparm1.vector, dev2);

        gsl_vector_view ss = gsl_vector_view_array(supersig, NOTHREADS * n_all_parameters * n_all_parameters);
        gsl_vector_view ssv = gsl_vector_subvector(&ss.vector, ithread * n_all_parameters * n_all_parameters, n_all_parameters * n_all_parameters);
        gsl_matrix_view ssm = gsl_matrix_view_vector(&ssv.vector, n_all_parameters, n_all_parameters);
        gsl_blas_dger(1.0, dev1, dev2, &ssm.matrix);
        gsl_vector_free(dev2);
      }
      gsl_vector_free(dev1);
    }
    if (sig) free(sig);
    if (gam) free(gam);
    if (temp) free(temp);
  }

  //helper frunction: reparameterize in going from phase 2 to 3
  void transit_from2_to3(int n_all_parameters, double* parmonstore, int n_value_store, double* valuestore, gsl_rng* rst1) {
    double liknorm[6]; double activeeps, epsm, Hobjective;
    for (int ithread = 0; ithread != NOTHREADS; ithread++) {
      double* tavw = (double*)malloc(ifreemax * 3 * indi * sizeof(double)); double* loglambdas = (double*)malloc(indi * sizeof(double));
      double* tlams = (double*)malloc(indi * respno * sizeof(double));
      gsl_vector* hampar = gsl_vector_alloc((phase <= 2) ? nhamil : n_all_parameters);
      int* paths = (int*)malloc(datenzahl * sizeof(int)); int* nips = (int*)malloc(no_patterns * 2 * indi * sizeof(int));
      double* ai = (double*)malloc(icompg * sizeof(double)); 	double* bi = (double*)malloc(respno * sizeof(double));
      double* parmon = (double*)malloc(2 * n_all_parameters * sizeof(double));
      double* alltaus = (double*)malloc(ntau * sizeof(double));
      double* rest = (double*)malloc(datenzahl * sizeof(double));

      pop(ithread, n_value_store, n_all_parameters, hampar, tavw, tlams, ai, loglambdas, bi, alltaus, rest, datenzahl, paths, nips, liknorm, activeeps, epsm, Hobjective, valuestore, parmon, parmonstore);
      Hobjective = 0.0; epsm = 0.0;
      for (int ii = 0; ii != 6; ii++) {liknorm[ii] = 0.0;} muplus = log(10 * activeeps);
      double omega;

      gsl_vector_view t1 = gsl_vector_view_array(loglambdas, indi);
      gsl_vector_view t2 = gsl_vector_subvector(hampar, isigoff, indi);
      gsl_vector_swap(&t1.vector, &t2.vector);
      make_romega(hampar, loglambdas, omega, rst1);

      double* sig = 0; if (!(sig = (double*)malloc(icompg * icompg * sizeof(double)))) { Rprintf("Allocation failure2\n"); }
      double* gam = 0; if (!(gam = (double*)malloc(respno * respno * sizeof(double)))) { Rprintf("Allocation failure2\n"); }
      double* sigi = 0; if (!(sigi = (double*)malloc(icompg * icompg * sizeof(double)))) { Rprintf("Allocation failure2\n"); }
      double* gami = 0; if (!(gami = (double*)malloc(respno * respno * sizeof(double)))) { Rprintf("Allocation failure2\n"); }
      gsl_matrix* Ltminusx = gsl_matrix_alloc(icompg, icompg);
      gsl_matrix* Ltminusr = gsl_matrix_alloc(respno, respno);
      sample_sig(hampar, sig, sigi, Ltminusx, ai, rst1);
      make_rgam(hampar, gam, gami, Ltminusr, bi, rst1);
      make_hampar_from_sig(0, sig, hampar);
      make_hampar_from_sig(1, gam, hampar);

      gsl_vector_set(hampar, n_all_parameters - 1, log(omega));

      push(ithread, n_value_store, n_all_parameters, hampar, tavw, tlams, ai, loglambdas, bi, alltaus, rest, datenzahl, paths, nips, liknorm, activeeps, epsm, Hobjective, valuestore, parmon, parmonstore);
      gsl_vector_free(hampar); free(tavw); free(tlams);  free(paths); free(nips); free(loglambdas);   free(ai); free(bi);  free(alltaus); free(rest);
      free(parmon);
      free(sig); free(sigi); free(gam); free(gami); gsl_matrix_free(Ltminusx); gsl_matrix_free(Ltminusr);
    }
  }

  //the sampler
  void gibbs_times_new(const std::vector<trial> & daten,  gsl_rng* rst1, gsl_rng* rst2, gsl_rng* rst3, gsl_rng* rst4) {

    int irun, ioff, factor;
    //settings for adaptive choice of stepsize
    double xeps = 0.095;
    muplus = log(10 * xeps);
    double epsm = 0.0; double activeeps = xeps; double Hobjective = 0.0;
    double liknorm[6] = { 0 * 6 };
    monitor = (double*)calloc(2 * 2 * 10, sizeof(double));

    //minimum number of samples for posterior variance/covariance estimate in phase 2 and for online updating in phase 3
    int interval = std::max(PHASE2, 5 * n_all_parameters);
    //	interval = 100;
    interval = (interval / ireps + 1) * ireps;

    int n_value_store = n_all_parameters + icompg * indi + icompg + indi +       respno + respno * indi + datenzahl + no_patterns * 2 * indi + ntau + datenzahl + 6 + 3;
    //                                        tavw         ai       loglambdas		  bi   	 tlams  paths         nips		alltaus, rest,	liknorm,  activeeps,epsm,Hobjective
    double* valuestore = 0; if (!(valuestore = (double*)malloc(NOTHREADS * n_value_store * sizeof(double)))) { Rprintf("Allocation failure\n"); }
    double* parmonstore = 0; if (!(parmonstore = (double*)malloc(NOTHREADS * 2 * n_all_parameters * sizeof(double)))) { Rprintf("Allocation failure\n"); }

    int satemp = NOTHREADS * IREP * (n_all_parameters);
    double* sample = 0; if (!(sample = (double*)malloc(satemp * sizeof(double)))) { Rprintf("Allocation failure\n"); }

    //initialize parameters 0 = randomly 1 = based on individual max. lik. estimates
    if (!goon) initialize(INITIALIZE, daten, xeps, parmonstore, n_value_store, valuestore, rst1, rst2, rst3, rst4);
    // main loop GIBBS
    double* xwbr = 0; if (!(xwbr = (double*)calloc(3 * n_all_parameters, sizeof(double)))) { Rprintf("Allocation failure\n"); }
    bool save = false;
    // double* complete_sample = 0;
    int sample_size2;
    preptrees(trees);

    RESTART:
      gsl_vector_view t0 = gsl_vector_view_array(xwbr, 3 * n_all_parameters);
    gsl_vector_set_zero(&t0.vector);

    gsl_vector_view t1 = gsl_vector_view_array(parmonstore, 2 * n_all_parameters * NOTHREADS);
    gsl_vector_set_zero(&t1.vector);
    double rmax = 0.0; int imax = -1; double epshelp;
    irun = -1;


    sample_size = SAMPLE_SIZE;
    ioff = 0;
    factor = 1;
    //ioff: how many packages have been saved in units of (ireps * NOTHREADS)/thin
    //factor: how many have been saved at the end of phase 4 in units of SAMPLE_SIZE
    if (goon) {
      phase = 4;
      save = true;
      sample_size = ADDITION;
      ioff = 0; factor = 1;
      std::ifstream raustemp;
      raustemp.open(RAUS);
      raustemp >> sample_size2 >> n_all_parameters;
      //		raustemp.close();
      std::ofstream rausFiles;
      for (int is = 0; is != NOTHREADS; is++) {
        rausFiles.open(std::string(TMPDIR) + "temp" + std::to_string(is));
        rausFiles << std::setprecision(12);
        for (int j = 0; j != sample_size2 / NOTHREADS; j++) {
          for (int i = 0; i != n_all_parameters; i++) {
            double x;
            raustemp >> x;
            rausFiles << std::setw(20) << x;
          }
          rausFiles << std::endl;
        }
        rausFiles.close();
        rausFiles.clear();
      }
      raustemp.close();
      raustemp.clear();
      std::string tempPath = std::string(TMPDIR) + "temp";
      std::rename(RAUS, tempPath.c_str());
      pop_continue(n_value_store, irun, valuestore, parmonstore, rst1, rst2, rst3, rst4);
      if (!(complete_sample = (double*)malloc(sample_size * (n_all_parameters) * sizeof(double)))) { Rprintf("Allocation failure\n"); }
    }
    // reicht nicht iresp = IREP generell, wenn Phase1 und Phase2 Vielfaches davon?
    if (phase == 1) ireps = std::min(IREP, PHASE1); else if (phase == 2) ireps = std::min(IREP, PHASE2); else ireps = IREP;

    WEITER: irun++;
    int offset = irun * ireps;
    //set up NOTHREADS chains and run chains
    std::vector<std::thread> threads(NOTHREADS-1);
    
    /* starting threads while ... */
    for (int ithread = 0; ithread < NOTHREADS-1; ithread++) {
      threads[ithread] = std::thread([&, ithread]() {
        // #pragma omp parallel for ordered shared(phase, offset, ireps, supsig, sigisqrt, supersig, epshelp, rst1,rst2,rst3,rst4,save,sample,ntau,n_value_store,daten ,valuestore,parmonstore,xwbr,imax,rmax,monitor)
        // 	for (int ithread = 0; ithread < NOTHREADS; ithread++) {
        double* tavw = 0;  double* tlams = 0; double* loglambdas = 0; gsl_vector* hampar = gsl_vector_alloc((phase <= 2) ? nhamil : n_all_parameters);
        int* paths = 0; int* nips = 0; double* ai = 0; double* bi = 0; double* parmon = 0;
        double* alltaus = 0; double* rest = 0;
        gsl_rng* rst;
        rst = gsl_rng_alloc(T_rng);
        parmon = (double*)malloc(2 * n_all_parameters * sizeof(double));

        tavw = (double*)malloc(ifreemax * 3 * indi * sizeof(double));
        tlams = (double*)malloc(respno * indi * sizeof(double));
        loglambdas = (double*)malloc(indi * sizeof(double));
        paths = (int*)malloc(datenzahl * sizeof(int));
        nips = (int*)malloc(no_patterns * 2 * indi * sizeof(int));
        ai = (double*)malloc(icompg * sizeof(double));
        bi = (double*)malloc(respno * sizeof(double));
        alltaus = (double*)malloc(ntau * sizeof(double));
        rest = (double*)malloc(datenzahl * sizeof(double));

        double liknorm[6];
        double epsm, activeeps, Hobjective;

        //run sampler ireps times
        switch (ithread + 1) {
        case 1: gsl_rng_memcpy(rst, rst1); break;
        case 2: gsl_rng_memcpy(rst, rst2); break;
        case 3: gsl_rng_memcpy(rst, rst3); break;
        case 4: gsl_rng_memcpy(rst, rst4); break;
        }
        pop(ithread, n_value_store, n_all_parameters, hampar, tavw, tlams, ai, loglambdas, bi, alltaus, rest, datenzahl, paths, nips, liknorm, activeeps, epsm, Hobjective, valuestore, parmon, parmonstore);
        gibbs_and_monitor(daten, nips, hampar, tavw, tlams, ai, loglambdas, bi, alltaus, rest, paths, liknorm, activeeps, epsm, Hobjective, offset, n_all_parameters, parmon, rst, ithread, save, sample);
        push(ithread, n_value_store, n_all_parameters, hampar, tavw, tlams, ai, loglambdas, bi, alltaus, rest, datenzahl, paths, nips, liknorm, activeeps, epsm, Hobjective, valuestore, parmon, parmonstore);
        switch (ithread + 1) {
        case 1: gsl_rng_memcpy(rst1, rst); break;
        case 2: gsl_rng_memcpy(rst2, rst); break;
        case 3: gsl_rng_memcpy(rst3, rst); break;
        case 4: gsl_rng_memcpy(rst4, rst); break;
        }
        
        int ido = 2;
        //		std::cout << setw(5) << ithread << setw(20) << activeeps << std::endl;
        //r statitstics
        while (curr_order.load() != ithread) {}
        // #pragma omp ordered
        {
          if (ithread == 0) ido = 1; //if (ithread + 1 == NOTHREADS) ido = 3;
          r_statistic(ido, n_all_parameters, ithread, offset + ireps, parmon, xwbr, rmax, imax);

          //prepare adapt stepsize in phase 1 and phase 3 (each time posterior variance/covariance matrix has been updated)
          if (((((offset + ireps) % interval) == PHASE1) && (phase % 2 == 1)) && (!(save))) {
            if (ithread == 0) epshelp = 0.0;
            epshelp += epsm;
            valuestore[(ithread + 1) * n_value_store - 3] = exp(epsm);
          }
        }
        curr_order++;
        gsl_rng_free(rst);
        gsl_vector_free(hampar); free(tavw); free(tlams);  free(paths); free(nips); free(loglambdas);   free(ai); free(bi); free(parmon);
        free(alltaus); free(rest);
        
      });
    }

    /* ... the main thread also runs */
    {
      double* tavw = 0;  double* tlams = 0; double* loglambdas = 0; gsl_vector* hampar = gsl_vector_alloc((phase <= 2) ? nhamil : n_all_parameters);
      int* paths = 0; int* nips = 0; double* ai = 0; double* bi = 0; double* parmon = 0;
      double* alltaus = 0; double* rest = 0;
      gsl_rng* rst;
      rst = gsl_rng_alloc(T_rng);
      parmon = (double*)malloc(2 * n_all_parameters * sizeof(double));
      tavw = (double*)malloc(ifreemax * 3 * indi * sizeof(double));
      tlams = (double*)malloc(respno * indi * sizeof(double));
      loglambdas = (double*)malloc(indi * sizeof(double));
      paths = (int*)malloc(datenzahl * sizeof(int));
      nips = (int*)malloc(no_patterns * 2 * indi * sizeof(int));
      ai = (double*)malloc(icompg * sizeof(double));
      bi = (double*)malloc(respno * sizeof(double));
      alltaus = (double*)malloc(ntau * sizeof(double));
      rest = (double*)malloc(datenzahl * sizeof(double));

      // double liknorm[6];
      // double epsm, activeeps, Hobjective;
      //run sampler ireps times
      switch (NOTHREADS) {
      case 1: gsl_rng_memcpy(rst, rst1); break;
      case 2: gsl_rng_memcpy(rst, rst2); break;
      case 3: gsl_rng_memcpy(rst, rst3); break;
      case 4: gsl_rng_memcpy(rst, rst4); break;
      }
      pop(NOTHREADS-1, n_value_store, n_all_parameters, hampar, tavw, tlams, ai, loglambdas, bi, alltaus, rest, datenzahl, paths, nips, liknorm, activeeps, epsm, Hobjective, valuestore, parmon, parmonstore);
      gibbs_and_monitor(daten, nips, hampar, tavw, tlams, ai, loglambdas, bi, alltaus, rest, paths, liknorm, activeeps, epsm, Hobjective, offset, n_all_parameters, parmon, rst, NOTHREADS-1, save, sample);
      push(NOTHREADS-1, n_value_store, n_all_parameters, hampar, tavw, tlams, ai, loglambdas, bi, alltaus, rest, datenzahl, paths, nips, liknorm, activeeps, epsm, Hobjective, valuestore, parmon, parmonstore);
      switch (NOTHREADS) {
      case 1: gsl_rng_memcpy(rst1, rst); break;
      case 2: gsl_rng_memcpy(rst2, rst); break;
      case 3: gsl_rng_memcpy(rst3, rst); break;
      case 4: gsl_rng_memcpy(rst4, rst); break;
      }

      int ido = 2;
      //		std::cout << setw(5) << ithread << setw(20) << activeeps << std::endl;
      //r statitstics
      while (curr_order.load() != NOTHREADS-1) {}
      // #pragma omp ordered
      {
        ido = 3;
        // if (ithread == 0) ido = 1; if (ithread + 1 == NOTHREADS) ido = 3;
        r_statistic(ido, n_all_parameters, NOTHREADS-1, offset + ireps, parmon, xwbr, rmax, imax);

        //prepare adapt stepsize in phase 1 and phase 3 (each time posterior variance/covariance matrix has been updated)
        if (((((offset + ireps) % interval) == PHASE1) && (phase % 2 == 1)) && (!(save))) {
          // if (NOTHREADS-1 == 0) epshelp = 0.0;
          epshelp += epsm;
          valuestore[(NOTHREADS) * n_value_store - 3] = exp(epsm);
        }
      }
      curr_order.store(0);
      gsl_rng_free(rst);
      gsl_vector_free(hampar); free(tavw); free(tlams);  free(paths); free(nips); free(loglambdas);   free(ai); free(bi); //free(parmon);
      free(alltaus); free(rest);

    }

    /* join threads */
    for (int ithread = 0; ithread < NOTHREADS-1; ithread++) {
      threads[ithread].join();
    }
    
    // #pragma omp parallel for ordered shared(phase, offset, ireps, supsig, sigisqrt, supersig, epshelp, rst1,rst2,rst3,rst4,save,sample,ntau,n_value_store,daten ,valuestore,parmonstore,xwbr,imax,rmax,monitor)
    // 	for (int ithread = 0; ithread < NOTHREADS; ithread++) {
    // 		double* tavw = 0;  double* tlams = 0; double* loglambdas = 0; gsl_vector* hampar = gsl_vector_alloc((phase <= 2) ? nhamil : n_all_parameters);
    // 		int* paths = 0; int* nips = 0; double* ai = 0; double* bi = 0; double* parmon = 0;
    // 		double* alltaus = 0; double* rest = 0;
    // 		gsl_rng* rst;
    // 		rst = gsl_rng_alloc(T_rng);
    // 		parmon = (double*)malloc(2 * n_all_parameters * sizeof(double));
    //
    // 		if (!(tavw = (double*)malloc(ifreemax * 3 * indi * sizeof(double)))) { Rprintf("Allocation failure\n"); }
    // 		if (!(tlams = (double*)malloc(respno * indi * sizeof(double)))) { Rprintf("Allocation failure\n"); }
    // 		if (!(loglambdas = (double*)malloc(indi * sizeof(double)))) { Rprintf("Allocation failure\n"); }
    // 		if (!(paths = (int*)malloc(datenzahl * sizeof(int)))) { Rprintf("Allocation failure2\n"); }
    // 		if (!(nips = (int*)malloc(no_patterns * 2 * indi * sizeof(int)))) { Rprintf("Allocation failure\n"); }
    // 		if (!(ai = (double*)malloc(icompg * sizeof(double)))) { Rprintf("Allocation failure\n"); }
    // 		if (!(bi = (double*)malloc(respno * sizeof(double)))) { Rprintf("Allocation failure\n"); }
    // 		if (!(alltaus = (double*)malloc(ntau * sizeof(double)))) { Rprintf("Allocation failure\n"); }
    // 		if (!(rest = (double*)malloc(datenzahl * sizeof(double)))) { Rprintf("Allocation failure\n"); }
    //
    // 		double liknorm[6];
    // 		double epsm, activeeps, Hobjective;
    //
    // 		//run sampler ireps times
    // 		switch (ithread + 1) {
    // 		case 1: gsl_rng_memcpy(rst, rst1); break;
    // 		case 2: gsl_rng_memcpy(rst, rst2); break;
    // 		case 3: gsl_rng_memcpy(rst, rst3); break;
    // 		case 4: gsl_rng_memcpy(rst, rst4); break;
    // 		}
    // 		pop(ithread, n_value_store, n_all_parameters, hampar, tavw, tlams, ai, loglambdas, bi, alltaus, rest, datenzahl, paths, nips, liknorm, activeeps, epsm, Hobjective, valuestore, parmon, parmonstore);
    // 		gibbs_and_monitor(daten, nips, hampar, tavw, tlams, ai, loglambdas, bi, alltaus, rest, paths, liknorm, activeeps, epsm, Hobjective, offset, n_all_parameters, parmon, rst, ithread, save, sample);
    // 		push(ithread, n_value_store, n_all_parameters, hampar, tavw, tlams, ai, loglambdas, bi, alltaus, rest, datenzahl, paths, nips, liknorm, activeeps, epsm, Hobjective, valuestore, parmon, parmonstore);
    // 		switch (ithread + 1) {
    // 		case 1: gsl_rng_memcpy(rst1, rst); break;
    // 		case 2: gsl_rng_memcpy(rst2, rst); break;
    // 		case 3: gsl_rng_memcpy(rst3, rst); break;
    // 		case 4: gsl_rng_memcpy(rst4, rst); break;
    // 		}
    //
    // 		int ido = 2;
    // //		std::cout << setw(5) << ithread << setw(20) << activeeps << std::endl;
    // 		//r statitstics
    // #pragma omp ordered
    // 		{
    // 			if (ithread == 0) ido = 1; if (ithread + 1 == NOTHREADS) ido = 3;
    // 			r_statistic(ido, n_all_parameters, ithread, offset + ireps, parmon, xwbr, rmax, imax);
    //
    // 			//prepare adapt stepsize in phase 1 and phase 3 (each time posterior variance/covariance matrix has been updated)
    // 			if (((((offset + ireps) % interval) == PHASE1) && (phase % 2 == 1)) && (!(save))) {
    // 				if (ithread == 0) epshelp = 0.0;
    // 				epshelp += epsm;
    // 				valuestore[(ithread + 1) * n_value_store - 3] = exp(epsm);
    // 			}
    // 		}
    // 		gsl_rng_free(rst);
    // 		gsl_vector_free(hampar); free(tavw); free(tlams);  free(paths); free(nips); free(loglambdas);   free(ai); free(bi); free(parmon);
    // 		free(alltaus); free(rest);
    // 	}
    double* parmon = (double*)malloc(2 * n_all_parameters * sizeof(double));
    gsl_vector_view t3 = gsl_vector_view_array(parmonstore, 2 * n_all_parameters * NOTHREADS);
    gsl_vector_view t4 = gsl_vector_subvector(&t3.vector, 0, 2 * n_all_parameters);
    gsl_vector_view t5 = gsl_vector_view_array(parmon, 2 * n_all_parameters);
    gsl_vector_memcpy(&t5.vector, &t4.vector);
    
    R_CheckUserInterrupt();
    //show interim results
    on_screen3(n_all_parameters, xwbr, parmon, consts, rmax, imax, irun);
    free(parmon);
    R_CheckUserInterrupt();
    
    //from phase 1 to phase 2
    if (((phase == 1) && (offset + ireps >= PHASE1)) && (rmax <= 50.0)) {
      supersig = (double*)calloc(NOTHREADS * n_all_parameters * n_all_parameters, sizeof(double));
      phase = 2;
      goto RESTART;
    }
    
    //from phase 2 to phase 3
    if (((phase == 2) && (offset + ireps >= interval)) && (rmax <= 50.0)) {
      make_supersigs(offset + ireps, parmonstore, supsig, sigisqrt);
      gsl_vector_view ty = gsl_vector_view_array(supersig, NOTHREADS * n_all_parameters * n_all_parameters);
      gsl_vector_set_zero(&ty.vector);
      phase = 3;
      transit_from2_to3(n_all_parameters, parmonstore, n_value_store, valuestore, rst1);
      //		char x; std::cin >> x;
      goto RESTART;
    }
    
    //in phase 3,  recalibrate
    if (((phase == 3) && ((offset + ireps) % interval) == 0)) {//&& (!save)) {
      make_supersigs(offset + ireps, parmonstore, supsig, sigisqrt);

      for (int ii = 0; ii != NOTHREADS; ii++) {
        valuestore[(ii + 1) * n_value_store - 2] = 0.0;
        valuestore[(ii + 1) * n_value_store - 1] = 0.0;
      }
      activeeps = valuestore[n_value_store - 3];
      muplus = log(10 * activeeps);
    }



    //start sampling if phase 3, new supersig has been computed, and rmax <= RMAX! i.e. go to phase 4
    if ((!save) && (rmax <= RMAX) && (offset + ireps >= interval) && ((offset + ireps)%interval >= PHASE1) && (phase == 3)) {
      gsl_vector_view ty = gsl_vector_view_array(supersig, NOTHREADS * n_all_parameters * n_all_parameters);
      gsl_vector_set_zero(&ty.vector);
      save = true;
      phase = 4;
      std::remove(RAUS);
      if (!(complete_sample = (double*)malloc(sample_size * (n_all_parameters) * sizeof(double)))) { Rprintf("Allocation failure\n"); }
      goto RESTART;
    }

    if (!save) goto WEITER;

    //save package (NOTHREAD*ireps)/THIN
  #define COMPLETE_SAMPLE(S,I,P) complete_sample[(S)*(sample_size/NOTHREADS)*(n_all_parameters) + (I)*(n_all_parameters)+P]
    if (save) {
      int irt = ireps / THIN;
      for (int is = 0; is != (NOTHREADS); is++)
        for (int i = 0; i != irt; i++) {
          for (int j = 0; j != n_all_parameters; j++)
            COMPLETE_SAMPLE(is, ioff * IREP / THIN + i, j) = dSAMPLE(is * IREP + i * THIN, j);
        }

        ioff += 1;
      if (ioff * NOTHREADS * ireps / THIN < sample_size) goto WEITER;
    }


    // logistics for saving in one file ordered by threads
    std::ofstream rausFiles;

    // zur Erinnerung:	if (goon) sample_size = ADDITION;
    int sano = sample_size / NOTHREADS;

    for (int is = 0; is != NOTHREADS; is++) {
      rausFiles.open(std::string(TMPDIR) + "raus" + std::to_string(is), std::ios_base::app);
      rausFiles << std::setprecision(12);
      for (int i = 0; i != sano; i++) {
        for (int j = 0; j != n_all_parameters; j++)
          rausFiles << std::setw(20) << COMPLETE_SAMPLE(is, i, j);
        rausFiles << std::endl;
      }
      rausFiles.close();
      rausFiles.clear();
    }
    if ((rmax > RMAX) && !goon) { ioff = 0; factor += 1;  goto WEITER; }

    //for diagnosis.cpp
    // Reihenfolge umkehren : Zuerst alt, dann fortsetzung, wegen Autokorrelationen
    if (goon) {
      std::ifstream inputFiles;
      for (int is = 0; is != NOTHREADS; is++) {
        inputFiles.open(std::string(TMPDIR) + "temp" + std::to_string(is));
        rausFiles.open(std::string(TMPDIR) + "raus" + std::to_string(is), std::ios_base::app);
        rausFiles << std::setprecision(12);
        rausFiles << inputFiles.rdbuf();
        inputFiles.close();
        inputFiles.clear();
        rausFiles.close();
        rausFiles.clear();
        std::string t = std::string(TMPDIR) + "temp" + std::to_string(is);
        std::remove(t.c_str());
      }
      std::string tempPath = std::string(TMPDIR) + "temp";
      std::remove(tempPath.c_str());
    }

    if (goon) sample_size = sample_size2 + ADDITION; else sample_size = factor * SAMPLE_SIZE;

    std::ofstream outputFile(RAUS);
    outputFile << std::setw(5) << sample_size << " " << n_all_parameters << std::endl;

    std::ifstream inputFiles;
    for (int is = 0; is != NOTHREADS; is++) {
      inputFiles.open(std::string(TMPDIR) + "raus" + std::to_string(is));
      outputFile << inputFiles.rdbuf();
      inputFiles.close();
      inputFiles.clear();
      std::string s = std::string(TMPDIR) + "raus" + std::to_string(is);
      std::remove(s.c_str());
    }
    outputFile.close();

    push_continue(n_value_store, irun, valuestore, parmonstore, rst1, rst2, rst3, rst4);

    removetrees(trees);

    if (monitor) free(monitor);
    if (valuestore) free(valuestore);
    if (parmonstore) free(parmonstore);
    if (xwbr) free(xwbr);
    if (sample) free(sample);
    //if (complete_sample) free(complete_sample);
    if (supersig) free(supersig);
    R_CheckUserInterrupt();
  }

}
