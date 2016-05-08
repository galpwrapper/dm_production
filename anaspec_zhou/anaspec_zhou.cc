#ifndef NO_ROOT

#include <iostream>
#include <cstring>
#include "TFile.h"
#include "TGraph.h"

#include "w_interp.h"
#include "anaspec_zhou.h"
#define ROOTNAME "DM_Spectrum.root"

using std::cout;
using std::endl;

namespace zhou {
  zhou_branch branch_map[] = {
    e, mu, tau,
    W, q,
    (zhou_branch)-1, bottom, top,
    (zhou_branch)-1,
    four_e, four_mu, four_tau
  };
  zhou_product product_map[] = {
    electron,
    antiproton,
    (zhou_product)-1,
    (zhou_product)-1,
    (zhou_product)-1
  };
}

const bool anaspec_zhou::exist[product_choice_size][branch_choice_size] = {{1,1,1,1,0,0,1,0,0,1,1,1},
                                                                           {0,1,0,1,1,0,1,1,0,0,0,0},
                                                                           {0,0,0,0,0,0,0,0,0,0,0,0},
                                                                           {0,0,0,0,0,0,0,0,0,0,0,0},
                                                                           {0,0,0,0,0,0,0,0,0,0,0,0}
};

const double anaspec_zhou::eMass = 0.5109990615e-3;
const double anaspec_zhou::pMass = 939.e-3;

anaspec_zhou::anaspec_zhou(): anaspec() {
  memset(tg, NULL, sizeof(tg));
  memset(loaded, false, sizeof(loaded));
}

double anaspec_zhou::ask(double E, double mdm, bool cumulate) {
  if (cumulate) {
    cout<<"cumulate was not usefull now"<<endl;
    return 0;
  }
  double result = 0;
  for (int i=0; i<branch_choice_size; i++) {
    if (branch[i] && exist[product][i]) {
      if (!loaded[product][i]) {
        if (load(zhou::product_map[product], zhou::branch_map[i])) {
          loaded[product][i] = true;
        } else {
          cout<<"[ERROR]: anaspec_zhou::ask"<<endl;
          cout<<"         load spectrum fail"<<endl;
          return 0;
        }
      }
      result += zhou_ask(E, mdm, zhou::product_map[product],zhou::branch_map[i]);
    }else if (branch[i] && !exist[product][i]) {
      cout<<"[ERROR]: anaspec_zhou::ask"<<endl;
      cout<<"         product or branch not exist"<<endl;
    }
  }
  return result;
}

bool anaspec_zhou::load(zhou::zhou_product prod, zhou::zhou_branch bran) {
#ifdef DATDIR
  string filename = DATDIR;
#else
  string filename = ".";
#endif
  filename = filename + "/" + ROOTNAME;
  string tgname = "G_" + zhou::branch_name[bran] + "_" + zhou::product_name[prod];
  TFile f(filename.c_str());
  f.GetObject(tgname.c_str(), tg[prod][bran]);
  f.Close();
  if (tg[prod][bran])
    return true;
  else
    return false;
}

double anaspec_zhou::zhou_ask(double E, double mdm, zhou::zhou_product prod, zhou::zhou_branch bran) {
  int n = tg[prod][bran]->GetN();
  double mass = 0;
  if (prod == zhou::electron) {
    mass = eMass;
  } else if (prod == zhou::antiproton) {
    mass = pMass;
  }
  double xtemp = (E + mass) / mdm,
         result = 0,
         * x = tg[prod][bran]->GetX(),
         * dndx = tg[prod][bran]->GetY();
  if (xtemp > 1) {
    return 0;
  }
  if (WMATH::interp_1D_loglog(x, dndx, n, &xtemp, &result, 1)) {
    return result / mdm;
  } else {
    cout<<"[ERROR]: anaspec_zhou::zhou_ask"<<endl;
    cout<<"         interp_1D_loglog fail"<<endl;
    return 0;
  }
}

#endif
