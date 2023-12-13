#include "rts.h"

namespace drtmpt {

  //create new theta
  struct Theta* newTheta() {
  	struct Theta* theta = (struct Theta*)malloc(sizeof(struct Theta));
  	int n = ifreemax * 3 * indi;
  	theta->hampar = gsl_vector_alloc((phase <= 2) ? nhamil : n_all_parameters);
  	theta->tavw = (double*)malloc(n * sizeof(double));
  	theta->loglambda = (double*)malloc(indi * sizeof(double));
  	theta->tlams = (double*)malloc(indi * respno * sizeof(double));
  	return theta;
  }

  //delete theta and free allocated memory
  void remove_Theta(Theta*& theta) {
  	gsl_vector_free(theta->hampar);
  	if (theta->loglambda) free(theta->loglambda);
  	if (theta->tavw) free(theta->tavw);
  	if (theta->tlams) free(theta->tlams);
  	if (theta) free(theta);
  }

  //copy theta
  void thetacopy(Theta*& thetadest, Theta* thetasrc) {
  	if (thetadest == NULL) thetadest = newTheta();
  	gsl_vector_view t1 = gsl_vector_view_array(thetadest->loglambda, indi);
  	gsl_vector_view t2 = gsl_vector_view_array(thetasrc->loglambda, indi);
  	gsl_vector_memcpy(&t1.vector, &t2.vector);
  	int n = 3 * ifreemax * indi;
  	gsl_vector_view t3 = gsl_vector_view_array(thetadest->tavw, n);
  	gsl_vector_view t4 = gsl_vector_view_array(thetasrc->tavw, n);
  	gsl_vector_memcpy(&t3.vector, &t4.vector);
  	n = indi * respno;
  		gsl_vector_view t5 = gsl_vector_view_array(thetadest->tlams, n);
  		gsl_vector_view t6 = gsl_vector_view_array(thetasrc->tlams, n);
  		gsl_vector_memcpy(&t5.vector, &t6.vector);
  	gsl_vector_memcpy(thetadest->hampar, thetasrc->hampar);
  }

  //create new node
  struct Node* newNode() {
  	struct Node* node = (struct Node*)malloc(sizeof(struct Node));
  	node->status = 0;
  	node->level = 0;
  	node->index = 0;
  	node->left = NULL;
  	node->right = NULL;
  	return node;
  }

  //destroy node
  void removenode(struct Node*& node) {
  	if (node == NULL) // std::cout << "Node gleich Null"
  		;
  	else {
  		free(node);
  		node = NULL;
  	}
  }

  //copy vector of momentum variables
  void pcopy(gsl_vector* pdest, gsl_vector* psrc) {
  	gsl_vector_memcpy(pdest, psrc);
  }

  //helper function in creating tree for NUTs sampler
  void count_increment(std::vector<bool>& index, int& count) {
  	if (index.size() == 0) {
  		count = 0;
  		index.push_back(false);
  		return;
  	}

  	int first_zero = 0;
  	while ((first_zero < static_cast<int>(index.size())) && index[first_zero]) first_zero++;
  	if (first_zero == static_cast<int>(index.size())) {
  		index.flip();
  		index.push_back(true);
  		count = 1;
  	}
  	else {
  		int fzp1 = first_zero + 1;
  		for (int i = 0; i != fzp1; i++) index[i] = !index[i];
  		count -= first_zero - 1;
  	}
  }

  //create tree of depth j
  struct Node* make_tree2(int j) {
  	std::vector<Node*> current, lower;
  	int levelsize = static_cast<int>(pow(2, j));
  	int level = 0;
  	std::vector<bool> index; index.clear();  int count;
  	for (int i = 0; i != levelsize; i++) {
  		count_increment(index, count);
  		struct Node* help = newNode();
  		if (GSL_IS_EVEN(i)) help->status = 0; else help->status = 1;
  		help->index = count; help->level = level;
  		current.push_back(help);
  	}
  	while (current.size() > 1) {
  		lower = current;
  		index.clear();
  		current.clear();
  		level++;
  		levelsize = static_cast<int>(lower.size()) / 2;
  		for (int i = 0; i != levelsize; i++) {
  			count_increment(index, count);
  			struct Node* help = newNode();
  			help->status = 2;
  			help->left = lower[2 * i];
  			help->right = lower[2 * i + 1];
  			help->index = count; help->level = level;
  			current.push_back(help);
  		}
  	}
  	return current[0];
  }


  //set up trees of depth up to 12
  void preptrees(Node* trees[13]) {
  	for (int j = 0; j != 13; j++) trees[j] = make_tree2(j);
  }

  //traverse tree
  void postorderTraversal(struct Node* &node) {
  	if (node == NULL)
  		return;
  	postorderTraversal(node->left);
  	postorderTraversal(node->right);
  	removenode(node);
  }

  //destroy tree structures
  void removetrees(Node* trees[13]) {
  	for (int j = 0; j != 13; j++) postorderTraversal(trees[j]);
  }

  //create storage for NUTs sampler
  struct store newstore() {
  	store speicher;
  	speicher.n = 0;
  	speicher.s = 0;
  	speicher.na = 0;
  	speicher.a = 0.0;
  	speicher.p = gsl_vector_alloc((phase <= 2) ? nhamil : n_all_parameters);
  	//	speicher.p = 0; if (!(speicher.p = (gsl_vector*)malloc(((icompg + respno) * (igroup + indi) + indi) * sizeof(double)))) { printf("Allocation failure\n"); }
  	speicher.theta = newTheta();
  	speicher.thetadash = newTheta();
  	return(speicher);
  }

  //derivatives of model likelihood by diffusion-model parameters
  void dhudwien(int* nips, gsl_vector* hampar, double* tavw, double* sigi,  double* alltaus, double* dstore, gsl_vector* dhampar) {


  	gsl_vector_view t1 = gsl_vector_subvector(dhampar, 0, (indi + igroup) * icompg);
  	gsl_vector_set_zero(&t1.vector);

  	int jj = 0;
  	for (int im = 0; im != no_patterns; im++) {
  		int iaa = dCOMB(im, 0);
  		int ivv = dCOMB(im, 1);
  		int iww = dCOMB(im, 2);
  		if ((dCOMP(0, iaa)) || (dCOMP(1, ivv)) || (dCOMP(2, iww)))
  			for (int t = 0; t != indi; t++) {
  				int itoff = t * 3 * ifreemax;
  				int ia = iaa + itoff;
  				int iv = ivv + itoff + ifreemax;
  				int iw = iww + itoff + 2 * ifreemax;

  				double aa = tavw[ia];
  				double va = tavw[iv];
  				double wa = tavw[iw];

  				double tausum = 0.0;
  				double da = 0.0, dv = 0.0, dw = 0.0;
  				int nntim = dNNODES(t, im);
  				for (int i = 0; i != nntim; i++) for (int pm = 0; pm != 2; pm++) {
  					dstore[jj] = dwiener_d(alltaus[jj], aa, va, wa, accuracy * 1.2);

  					if (dCOMP(0, iaa)) da -= dadwiener_d(alltaus[jj], aa, va, wa, dstore[jj]);
  					if (dCOMP(2, iww)) dw -= dwdwiener_d(alltaus[jj], aa, va, wa, dstore[jj]);

  					if (dCOMP(1, ivv)) tausum += fabs(alltaus[jj]);
  					jj++;
  				}
  				if (dCOMP(1, ivv)) {
  					double temp = -aa * (2 * wa - 1.0) * dNNODES(t, im) - va * tausum;
  					dv -= temp;
  				}
  				for (int pm = 0; pm != 2; pm++) {
  					int n = dNIPS(t, pm, im);
  					if (n != 0) {
  						double dav;
  						if ((dCOMP(0, iaa)) || (dCOMP(1, ivv))) {
  							dav = davlogprob_upperbound(pm, aa, va, wa);
  							if (dCOMP(0, iaa)) da += n * dalogprob_upperbound(pm, aa, va, wa, dav);
  							if (dCOMP(1, ivv)) dv += n * dvlogprob_upperbound(pm, aa, va, wa, dav);
  						}
  						if (dCOMP(2, iww)) dw += n * dwlogprob_upperbound(pm, aa, va, wa);
  					}
  				}
  				int ima = dmapAVW(t, 0, iaa), imv = dmapAVW(t, 1, ivv), imw = dmapAVW(t, 2, iww);
  				if (dCOMP(0, iaa)) gsl_vector_set(dhampar, ima, gsl_vector_get(dhampar, ima) + da);
  				if (dCOMP(1, ivv)) gsl_vector_set(dhampar, imv, gsl_vector_get(dhampar, imv) + dv);
  				if (dCOMP(2, iww)) gsl_vector_set(dhampar, imw, gsl_vector_get(dhampar, imw) + dw);
  			}
  		else for (int t = 0; t != indi; t++) jj += 2 * dNNODES(t, im);
  	}

  	gsl_vector* temp = gsl_vector_alloc(indi * icompg); jj = 0;
  	for (int t = 0; t!= indi; t++)
  		for (int type = 0; type != 3; type++) {
  			int ift = ifree[type];
  			for (int ip = 0; ip != ift; ip++) if (dCOMP(type, ip)) {
  				gsl_vector_set(temp, jj++, dlogit(avwtrans[type], gsl_vector_get(hampar, dmapMAVW(t2group[t], type, ip)) + gsl_vector_get(hampar, dmapAVW(t, type, ip))));
  			}
  		}

  	gsl_vector_view t2 = gsl_vector_subvector(dhampar, iavwoff, indi * icompg);

  	gsl_vector_mul(&t2.vector, temp);

  	gsl_vector_free(temp);

  	gsl_vector_view tm2 = gsl_vector_subvector(hampar, 0, igroup * icompg);
  	gsl_vector_view dtm2 = gsl_vector_subvector(dhampar, 0, igroup * icompg);

  	for (int t = 0; t != indi; t++) {
  		gsl_vector_view dgr = gsl_vector_subvector(dhampar, t2group[t] * icompg, icompg);
  		gsl_vector_view dindi = gsl_vector_subvector(dhampar, iavwoff + t * icompg, icompg);
  		gsl_vector_add(&dgr.vector, &dindi.vector);
  	}

  	gsl_blas_daxpy(PRIOR, &tm2.vector, &dtm2.vector);

  	gsl_matrix_view SINV = gsl_matrix_view_array(sigi, icompg, icompg);

  	gsl_vector_view t3 = gsl_vector_subvector(hampar, iavwoff, indi * icompg);
  	gsl_matrix_view xavwavw = gsl_matrix_view_vector(&t3.vector, indi, icompg);

  	gsl_matrix_view dstacked = gsl_matrix_view_vector(&t2.vector, indi, icompg);

  	gsl_blas_dsymm(CblasRight, CblasLower, 1.0, &SINV.matrix, &xavwavw.matrix, 1.0, &dstacked.matrix);
  }

}
