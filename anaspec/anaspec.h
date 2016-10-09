#ifndef _ANASPEC_H
#define _ANASPEC_H
#include<string>
#include<vector>
#include"Interp2D.h"
#include"enum_utils.h"

/*********************************************************************
This class is used to calculate generation spectrum from the
annihilation of darkmatter Eg and mdm are both in [GeV], the result
returned by function ask is in [GeV^-1].
Modified from anaspec.f(yuan)    linsj20140308
The cumulated mode will return the value of \int_0^{E_0} dN/dE dE
*********************************************************************/
class anaspec {
public:
#define X(a) a,
  enum product_choice {
#include "enumdef/product_choice.def"
  };
  enum branch_choice {
#include "enumdef/branch_choice.def"
  };
#undef X

#define X(TYPE) static ENUMNAME(TYPE)
  X(product_choice); X(branch_choice);
#undef X

  anaspec();
  anaspec(anaspec::product_choice product_);
  anaspec(double *branch_, anaspec::product_choice product_);
  anaspec(const std::vector <double> &branch_, anaspec::product_choice product_);

  static unsigned getenum(const std::string &enum_name);
  int rebranch(double *branch_);
  int rebranch(const std::vector <double> &branch_);
  int choose(anaspec::product_choice product_);
  virtual double ask(double E, double mdm, bool cumulate = false);

protected:
  anaspec::product_choice product;
  double branch[branch_choice_size];

private:
  static const bool exist[product_choice_size][branch_choice_size],
               withcum[product_choice_size][branch_choice_size];

  Interp2D table[product_choice_size][branch_choice_size],
           cumutable[product_choice_size][branch_choice_size];

  bool loaded[product_choice_size][branch_choice_size];

  int load(anaspec::product_choice prod, anaspec::branch_choice bra);
  int load(anaspec::product_choice prod, anaspec::branch_choice bra, bool cumutable);
};

#define X(TYPE) ENUMDECLEAR(anaspec::TYPE)
X(product_choice) X(branch_choice)
#undef X

#endif // for #ifndef _ANASPEC_H
