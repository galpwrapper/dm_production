#ifndef _ANASPEC_zhou_H
#define _ANASPEC_zhou_H

#ifndef NO_ROOT

#include "anaspec.h"
#include "TGraph.h"

/*********************************************************************
This class is used to calculate generation spectrum from the
annihilation of darkmatter Eg and mdm are both in [GeV], the result
returned by function ask is in [GeV^-1].
With the result of zhou
The cumulated mode will return the value of \int_0^{E_0} dN/dE dE
*********************************************************************/

using std::string;

namespace zhou{
  enum zhou_branch {e, mu, tau, q, bottom, top, W, Z, H, ZH,
                    four_e, four_mu, four_tau, zhou_branch_size};
  enum zhou_product {electron, antiproton, zhou_product_size};

  const string branch_name[zhou_branch_size] = {"2e", "2mu", "2tau", "2q", "2b", "2t", "2w", "2z", "2h", "zh",
                                                "4e", "4mu", "4tau"};
  const string product_name[zhou_product_size] = {"electron", "antiproton"};
}

class anaspec_zhou : public anaspec {
  public:
    anaspec_zhou();
    double ask(double E, double mdm, bool cumulate = false);

  private:
    TGraph *tg[zhou::zhou_product_size][zhou::zhou_branch_size];
    static const bool exist[product_choice_size][branch_choice_size];
    bool loaded[product_choice_size][branch_choice_size];
    static const double eMass; //Electron mass in Gev
    static const double pMass; //Proton mass in Gev

    bool load(zhou::zhou_product prod, zhou::zhou_branch bran);
    double zhou_ask(double E, double mdm, zhou::zhou_product prod, zhou::zhou_branch bran);
};

#endif // NO_ROOT

#endif
