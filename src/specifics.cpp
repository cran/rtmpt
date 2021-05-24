// authors: Christoph Klauer and Raphael Hartmann

#include "rts.h"
// namespace rtsNS {

	void lieszeile(std::ifstream& rein) {
		char x = ' '; while (x != '\n') rein.get(x);
	}



	//Tukey for outlier treatment

	// void tukey(vector<trial>& daten) {
	// 	int indi, kerntree, kerncat, ntot;
	// 	set_ns(daten, indi, kerntree, kerncat, igroup, ntot);
	//
	// 	vector<double> o; double ff[3]; double lower, upper; vector<trial> temp; temp.clear();
	//
	// 	for (int t = 0; t != indi; t++) {
	// 		o.clear();
	// 		for (int i = 0; i != int(daten.size()); i++) if (daten[i].person == t) o.push_back(daten[i].rt*1.0);
	//
	// 		int n = int(o.size());
	// 		sort(o.begin(), o.end());
	//
	// 		for (int q = 0; q != 3; ++q) {
	// 			double p = (q + 1)*0.25;
	// 			int j = int(p*(n + 1));
	// 			double f = p * (n + 1) - j;
	// 			ff[q] = (j >= 1) ? ((j < n) ? (1 - f)*o[j - 1] + f * o[j] : o[j - 1]) : o[0];
	// 		}
	// 		upper = ff[2] + 3.0*(ff[2] - ff[0]); lower = ff[0] - 3.0*(ff[2] - ff[0]);
	// 		for (int i = 0; i != int(daten.size()); i++) if ((daten[i].person == t) && (daten[i].rt >= lower) && (daten[i].rt <= upper)) temp.push_back(daten[i]);
	// 	}
	// 	Rprintf("%7d%7d\n", daten.size(), temp.size());
	// 	//char x; std::cin >> x;
	// 	daten = temp;
	// }

	//read in data and set category-to-response mapping

	void lies(vector<trial> &daten)
	{
		std::ifstream rein(DATA);
		lieszeile(rein);
		trial one;


		while (!rein.eof()) {
			rein >> one.person >> one.group >> one.tree >> one.category >> one.rt; one.item = 0;
			if (!rein.eof()) daten.push_back(one);
		}
		rein.close();

		int indi, kerntree, kerncat, ntot;
		set_ns(daten, indi, kerntree, kerncat, igroup, ntot);

		bool format = false;
		for (int ix = 0; ix != ntot; ix++) {
			trial one = daten[ix];
			if ((one.tree > 0) && (one.category == 0)) {
				format = true;
				break;
			}
		}
		if (format) {
			int indi, kerntree, kerncat, ntot;
			set_ns(daten, indi, kerntree, kerncat, igroup, ntot);
			int *cat2tree = 0, *tree2cat = 0;
			cat2tree = (int *)calloc(kerncat, sizeof(int));
			tree2cat = (int *)calloc(kerntree, sizeof(int));
			std::ifstream info(MODEL); int schrott;
			for (int j = 0; j != 5; j++) info >> schrott;
			for (int j = 0; j != kerncat; j++) { info >> cat2tree[j]; cat2tree[j]--; }
			info.close();
			//		int start = -1;
			for (int it = 0; it != kerntree; it++) {
				for (int j = 0; j != kerncat; j++) if (cat2tree[j] == it) tree2cat[it] = std::max(tree2cat[it], j);
				//			int temp = tree2cat[it];
				//			tree2cat[it] -= start;
				//			start = temp;
				tree2cat[it]++;
			}
			for (int ix = 0; ix != ntot; ix++) {
				daten[ix].category = (daten[ix].tree == 0) ? daten[ix].category : tree2cat[daten[ix].tree - 1] + daten[ix].category;
			}
		}


		vector<int>  o;
		int n = static_cast<int>(daten.size());
		for (int i = 0; i != static_cast<int>(daten.size()); i++) o.push_back(daten[i].rt);
		sort(o.begin(),o.end());
		if (DEBUG) Rprintf("%6d%6d%6d%6d%6d\n", o[0], o[n/4], o[n/2], o[(3*n)/4], o[n-1]);

		if (DEBUG) {for (int i=0;i!=20;i++) Rprintf("%6d", o[i]); Rprintf("\n");}
		if (DEBUG) {for (int i=0;i!=20;i++) Rprintf("%6d", o[n-1-i]); Rprintf("\n");}

		kerncat=nKERN;
		cat2resp = (int *)calloc(kerncat, sizeof(int));
		for (int i=0;i!=kerncat;i++) cat2resp[i]= CatToResp[i]; respno = nRESP;
	  // char y; std::cin >> y;
	}


	void model_design(int kerntree, int *ar, int *branch, int *nodes_per_par, int *nodes_per_tree, int *tree_and_node2par) {
	// #define A(I,J,K) a[I*zweig*kernpar + J*kernpar + K]
	// #define B(I,J,K) b[I*zweig*kernpar + J*kernpar + K]
	#define AR(I,J,K) ar[(I)*zweig*nodemax + J*nodemax + K]
	#define NODES_PER_PAR(I,J) nodes_per_par[I*kernpar + J]
	#define TREE_AND_NODE2PAR(I,J) tree_and_node2par[(I)*nodemax+(J)]
	#define DRIN(J,K,X) drin[J*zweig*nodemax+K*nodemax + X]
	#define NDRIN(J,K) ndrin[J*zweig+K]
		// zweig,kernpar,nodemax sind definiert;
/// ALT begin
		// int schrott;
		// std::ifstream info(MODEL); for (int j = 0; j != 5 + kerncat; j++) info >> schrott;
		// for (int j = 0; j != kerncat; j++) info >> branch[j];
		// //for (int j = 0; j != kerncat; j++) std::cout << setw(3) << branch[j];
		//
		// for (int it = 0; it != kerntree; it++) for (int in = 0; in != nodemax; in++) { info >> TREE_AND_NODE2PAR(it, in); TREE_AND_NODE2PAR(it, in)--; }
		// for (int it = 0; it != kerntree; it++) info >> nodes_per_tree[it];
		//
		// for (int i = 0; i != kerncat * zweig*nodemax; i++) ar[i] = 0;
		// for (int in = 0; in != nodemax; in++) for (int ip = 0; ip != zweig; ip++) for (int j = 0; j != kerncat; j++) info >> AR(j, ip, in);
		//
		// info.close();
		//
		// for (int i = 0; i != kerntree; i++) for (int j = 0; j != kernpar; j++) NODES_PER_PAR(i, j) = 0;
		// for (int i = 0; i != kerntree; i++) for (int r = 0; r != nodes_per_tree[i]; r++) {
		// 	NODES_PER_PAR(i, TREE_AND_NODE2PAR(i, r))++;
		// }
		//
		// for (int i = 0; i != kerncat * zweig*nodemax; i++) drin[i] = 0;
		// for (int i = 0; i != kerncat * zweig; i++) ndrin[i] = 0;
		// for (int j = 0; j != kerncat; j++)
		// 	for (int k = 0; k != branch[j]; k++)
		// 		for (int r = 0; r != nodes_per_tree[cat2tree[j]]; r++) if (AR(j, k, r) != 0) {
		// 			DRIN(j, k, NDRIN(j, k)) = r;
		// 			NDRIN(j, k)++;
		// 		}
/// ALT end
		bool auto_or_eqns = false;
		if (!auto_or_eqns) {
			int schrott;
			std::ifstream info(MODEL); for (int j=0;j!= 5+kerncat;j++) info>> schrott;
			for (int j=0;j!=kerncat;j++) info>> branch[j];
			if (DEBUG) {for (int j=0;j!=kerncat;j++) Rprintf("%3d", branch[j]); Rprintf("\n");}

			for (int it=0;it!=kerntree;it++) for (int in=0;in!=nodemax;in++) { info >> TREE_AND_NODE2PAR(it,in); TREE_AND_NODE2PAR(it,in)--; }
			for (int it=0;it!=kerntree;it++) info>> nodes_per_tree[it];

			for (int i=0;i!=kerncat*zweig*nodemax;i++) ar[i]=0;
			for (int in=0;in!=nodemax;in++) for (int ip=0;ip!=zweig;ip++) for (int j=0;j!=kerncat;j++) info>>AR(j,ip,in);

			info.close();
		}

		for (int i=0;i!=kerntree;i++) for (int j=0;j!=kernpar;j++) NODES_PER_PAR(i,j)=0;
		for (int i=0;i!=kerntree;i++) for (int r=0;r!=nodes_per_tree[i];r++) {
			NODES_PER_PAR(i,TREE_AND_NODE2PAR(i,r))++;
		}



		for (int i=0;i!=kerncat*zweig*nodemax;i++) drin[i]=0;
		for (int i=0;i!=kerncat*zweig;i++) ndrin[i]=0;
		for (int j=0;j!=kerncat;j++)
			for (int k=0;k!=branch[j];k++)
				for (int r=0;r!=nodes_per_tree[cat2tree[j]];r++) if (AR(j,k,r)!=0) {
					DRIN(j,k,NDRIN(j,k))=r;
					NDRIN(j,k)++;
				}
		// Konstanten SETZEN
		for (int i=0; i!=kernpar; i++) {
			if ((ConstProb[i]<=0.0) || (ConstProb[i]>=1.0)) comp[i]=true;
			else comp[i] = false;
		}
		//    comp[3]=comp[4]=false;

		ifree = 0;
		for (int i=0;i!=kernpar;i++) if (comp[i]) ifree++;
		// GEGEBENENFALLS konstanten setzen
		for (int i=0;i!=kernpar;i++) if (!comp[i]) consts[i]=ConstProb[i]; else	consts[i]=-1.0;


		// ZEITKRITISCHE PROZESSE SETZEN
		for (int i=0;i!=kernpar;i++) comp[kernpar+i] = CompMinus[i]==0 ? false : true;
		for (int i=0;i!=kernpar;i++) comp[2*kernpar+i] = CompPlus[i]==0 ? false : true;
		//comp[9]=comp[17]=false;
		ilamfree = 0;
		for (int i = 0; i != 2 * kernpar; i++) if (comp[kernpar + i]) ilamfree++;

	}

// }
