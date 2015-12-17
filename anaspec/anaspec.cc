#include"anaspec.h"
#include"mydebug.h"
#include<iostream>
#include<sstream>
#include<fstream>
#include<cmath>

using std::string;
using std::istringstream;
using std::ifstream;
using std::log10;
using std::vector;
using std::cout;
using std::endl;

const string anaspec::prodname[product_choice_size] = { "posi", "pbar", "gamma", "nu", "deuterons" },
    anaspec::braname[branch_choice_size] = { "e", "mu", "tau", "ww", "uu", "cc", "bb", "tt", "gg", "foure", "fourmu", "fourtau" };

#define NAME(TYPE) ENUMNAMECLS(TYPE, anaspec)
#define X(a) #a,
NAME(product_choice) = {
#include "enumdef/product_choice.def"
};
NAME(branch_choice) = {
#include "enumdef/branch_choice.def"
};
#undef X
#undef NAME

#define X(TYPE) ENUMANDSTR(anaspec::TYPE, anaspec::TYPE##_name)
X(product_choice) X(branch_choice)
#undef X

const bool anaspec::newformat[product_choice_size][branch_choice_size] = { {1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0},
                                                                        {1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0},
                                                                        {1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0},
                                                                        {1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0},
                                                                        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0} },
      anaspec::exist[product_choice_size][branch_choice_size] = { {1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0},
                                                               {0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0},
                                                               {1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0},
                                                               {1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0},
                                                               {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0} },
      anaspec::withcum[product_choice_size][branch_choice_size] = { {1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0},
                                                                 {0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0},
                                                                 {1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0},
                                                                 {1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0},
                                                                 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0} };

#define BRANCHING\
  for (int i = 0; i < branch_choice_size; i++) branch[i] = branch_[i]
#define INILOADED memset(loaded, false, branch_choice_size * product_choice_size * sizeof(bool))

anaspec::anaspec() {
  INILOADED;
}
anaspec::anaspec(anaspec::product_choice product_): product(product_) {
  INILOADED;
}
anaspec::anaspec(double *branch_, anaspec::product_choice product_): product(product_) {
  BRANCHING;
  INILOADED;
}
anaspec::anaspec(const vector <double> &branch_, anaspec::product_choice product_): product(product_) {
  BRANCHING;
  INILOADED;
}

int anaspec::rebranch(double *branch_) {
  BRANCHING;
  return 0;
}
int anaspec::rebranch(const vector <double> &branch_) {
  BRANCHING;
  return 0;
}
int anaspec::choose(anaspec::product_choice product_) {
  product = product_;
  return 0;
}

double anaspec::ask(double E, double mdm, bool cumulate) {
  printDebugMsg("Routine", ">>ask: E, m_dm, cumulate = %f, %f, %d", E, mdm, cumulate);
  double sum = 0,
         realx;
  double upbound = newformat[product][0] ? 0.999999 : 0.9772;

  for (int i = 0; i < branch_choice_size; i++) {
    printDebugMsg("Routine", "calculating branch %d with sum %f", i, sum);

    if (exist[product][i] && branch[i]) {
      if (!loaded[product][i]) load(product, (anaspec::branch_choice)i);
      if (cumulate && withcum[product][i] && E / mdm > 1e-5) {
        realx = (E / mdm) > 1 ? 1 : (E / mdm);
        sum += branch[i] * cumutable[product][i].lnask(realx, mdm);
      } else if(E / mdm > 1e-5 && E / mdm < upbound)
        sum += branch[i] * table[product][i].lnask(E / mdm, mdm) / mdm;
    }

  }

  printDebugMsg("Routine", "<<ask: %f", sum);
  return sum;
}

int anaspec::load(anaspec::product_choice prod, anaspec::branch_choice branch) {
  load(prod, branch, false);
  if(withcum[prod][branch]) load(prod, branch, true);
  loaded[prod][branch] = true;
  return 0;
}

int anaspec::load(anaspec::product_choice prod, anaspec::branch_choice branch, bool cumulate) {
#ifdef DATDIR
  string filename = DATDIR;
#else
  string filename = ".";
#endif
  filename += "/" + prodname[prod] + "_" + braname[branch] + (cumulate ? "_cum" : "" ) + ".dat";
  printDebugMsg("Routine", ">>load: filename = %s", filename.c_str());
  ifstream dats(filename.c_str());

  const double mfactor = 2;
  Table2D tab;
  string line;
  double x, y, val;
  vector <double> yaxis;

  if(newformat[prod][branch]) {
    getline(dats, line);
    istringstream iss(line);
    while(iss >> y) yaxis.push_back(y);
  }

  while (getline(dats, line)) {
    istringstream iss(line);
    iss >> x;
    x *= log(10);

    if(newformat[prod][branch])
      for(unsigned i = 0; (iss >> val) && i < yaxis.size(); i++) tab.insval(x, yaxis[i] * log(10), val * log(10));
    else
      for (y = 1; iss >> val; y += log10(mfactor)) tab.insval(x, y * log(10), val * log(10)); //times log(10) to change the base from 10 to e
  }

  cumulate ? cumutable[prod][branch].lntabling(tab) : table[prod][branch].lntabling(tab);
  printDebugMsg("Routine", "<<load");
  return 0;
}
