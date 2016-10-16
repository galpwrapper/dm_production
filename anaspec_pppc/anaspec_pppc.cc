#include"anaspec_pppc.h"
#include"mydebug.h"
#include<iostream>
#include<sstream>
#include<fstream>
#include<cmath>


using std::istringstream;
using std::ifstream;
using std::log10;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;

pppc::pppc_branch branch_map[anaspec::branch_choice_size] = {
  pppc::ele, pppc::mu, pppc::tau,
  pppc::W_boson,
  pppc::quark, //The quark channel of PPPC denotes a light quark
  pppc::charm, pppc::bottom, pppc::top,
  pppc::gluon,
  pppc::four_e, pppc::four_mu, pppc::four_tau,
  pppc::pppc_branch_size, pppc::pppc_branch_size // for four_pi or four_pi0 channel, set all the result to zero
};

const vector <vector <pppc::pppc_product> > product_map = { {pppc::positrons}, {pppc::antiprotons}, {pppc::gammas}, {pppc::neutrinos_e, pppc::neutrinos_mu, pppc::neutrinos_tau}, {pppc::antideuterons} };
const vector <vector <double> > product_norm = { {1}, {1}, {1}, {2.0/3, 2.0/3, 2.0/3}, {1} };

anaspec_pppc::anaspec_pppc(): anaspec(), loaded(product_choice_size, false) { }
anaspec_pppc::anaspec_pppc(anaspec::product_choice product_): anaspec(product_), loaded(product_choice_size, false) {
}
anaspec_pppc::anaspec_pppc(double *branch_, anaspec::product_choice product_): anaspec(branch_, product_), loaded(product_choice_size, false) {

}
anaspec_pppc::anaspec_pppc(const std::vector <double> &branch_, anaspec::product_choice product_): anaspec(branch_, product_), loaded(product_choice_size, false) {

}

double anaspec_pppc::ask(double E, double mdm, bool cumulate) {
  printDebugMsg("Routine", ">>ask: E, m_dm, cumulate = %f, %f, %d", E, mdm, cumulate);
  double result = 0;
  if (!loaded[product]) load(product);
  for(unsigned i = 0; i < product_map[product].size(); i++)
    result += pppc_ask(E, mdm, product_map[product][i], cumulate) * product_norm[product][i];
  printDebugMsg("Routine", "<<ask: %f", result);
  return result;
}
double anaspec_pppc::pppc_ask(double E, double mdm, pppc::pppc_product prod, bool cumulate) {
  printDebugMsg("Routine", ">>pppc_ask: E, m_dm, cumulate = %f, %f, %d, %d", E, mdm, prod, cumulate);
  double result = 0, x = E/mdm, realx;
  for (unsigned i = 0; i < branch_choice_size; i++)
    if(branch[i] && branch_map[i] < pppc::pppc_branch_size) {
      if(cumulate && x > low_x) {
        realx = x > 1 ? 1 : x;
        result += branch[i] * cumutable[prod][branch_map[i]].lnask(realx, mdm);
      } else if (x > low_x && x <= 1)
        result += branch[i] * table[prod][branch_map[i]].lnask(x, mdm) / mdm;
    }
  printDebugMsg("Routine", "<<pppc_ask: %f", result);
  return result;
}

int anaspec_pppc::load(pppc::pppc_product prod) {
#ifdef DATDIR
  string filename = DATDIR;
#else
  string filename = ".";
#endif
  filename += "/AtProduction_" + pppc::product_file[prod] +  ".dat";
  printDebugMsg("Routine", ">>load: filename = %s", filename.c_str());
  cerr << "loading PPPC file " << filename << endl;
  ifstream dats(filename.c_str());

  string line;
  Table2D tab[pppc::pppc_branch_size], cumtab[pppc::pppc_branch_size];
  vector <double> xaxis, mass;

  getline(dats, line);

  int mind = -1, xind = -1;
  double Nx[pppc::pppc_branch_size],
         mold = 0, mnew, xnew, val; // xold = 0
  printDebugMsg("Routine", "Go to getlines");
  while (getline(dats, line)) {
    istringstream iss(line);
    iss >> mnew; iss >> xnew;

    if (mnew != mold) {
      xind = -1;
      mind++;
      mass.push_back(mnew);
      for (unsigned i = 0; i < pppc::pppc_branch_size; i++) Nx[i] = 0;
    }
    xind++;
    if (mind <= 0) xaxis.push_back(pow(10, xnew));

    for (unsigned i_branch = 0; i_branch < pppc::pppc_branch_size; i_branch++) {
      iss >> val;
      double dndx = val / (xaxis[xind] * log(10));

      if (xind != 0) Nx[i_branch] += (xaxis[xind] - xaxis[xind - 1]) * dndx;

      tab[i_branch].insval(xaxis[xind], mass[mind], dndx);
      cumtab[i_branch].insval(xaxis[xind], mass[mind], Nx[i_branch]);
    }

    mold = mnew;// xold = xnew;
  }

  add_threshold_mass(tab);

  for (unsigned i_branch = 0; i_branch < pppc::pppc_branch_size; i_branch++) {
    printDebugMsg("Routine", "Init branch %d %d %d", i_branch, prod, anaspec::deuterons);
    table[prod][i_branch].tabling(tab[i_branch]);
    cumutable[prod][i_branch].tabling(cumtab[i_branch]);
  }

  low_x = (tab[0].xaxis.begin())->first / 2;

  printDebugMsg("Routine", "<<load");
  return 0;
}

int anaspec_pppc::load(anaspec::product_choice prod) {
  if (loaded[prod]) return 0;
  loaded[prod] = true;
  for (unsigned i = 0; i < product_map[prod].size(); i++)
    load(product_map[prod][i]);
  return 1;
}

int anaspec_pppc::add_threshold_mass(Table2D *tab) const {
  printDebugMsg("Routine", ">>add_threshold_mass");
  const vector <double> threshold_mass = { 0, 0, 0, 0, 0, 0, 0, 0, 0, //e, mu, tau
    0, 0, 0, //quark, charm, bottom; We set those mass smaller than 5GeV to be zero, and do nothing to these branch
    173.07, 80.4, 80.4, 80.4, 91.2, 91.2, 91.2, //t, W, Z
    0, 0, 126, 0, 0, 0, 0, 0, 0 };//gluon, gamma, higgs, nu, four_leptons

  for(unsigned i_branch = 0; i_branch < pppc::pppc_branch_size; i_branch++) {
    if (threshold_mass[i_branch] == 0) continue;

    double nextmass = tab[i_branch].yaxis.upper_bound(threshold_mass[i_branch])->first;
    for(Table2D::LineIter i = tab[i_branch].xaxis.begin(); i != tab[i_branch].xaxis.end(); i++) {
      tab[i_branch].insval(i->first, threshold_mass[i_branch], tab[i_branch].value[i->first][nextmass]);
    }
  }

  printDebugMsg("Routine", ">>add_threshold_mass");
  return 0;
}
