#ifndef TRANSITMASK_H
#define TRANSITMASK_H
#include <iostream>
#include <math.h>
using namespace std;
#include "maskarray.h"
class TransitMask :  public MaskArray{
  public:
    TransitMask(const double *cmap, int nbin, int ntr,double high,double low);
    TransitMask(const TransitMask &TM);
    ~TransitMask();
    void TransitModel(const double *Ma);
    void OutputTransit(double *M, double & high, double &low) const;
    void PrintTransit() const;
    void Chisquare(const double *modle, double & err) const;
    TransitMask &operator = (const TransitMask &Ma);
  private:
    double *model_;
    double high_;
    double low_;
    double* transit_;
};
#endif
