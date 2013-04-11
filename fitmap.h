#ifndef FITMAP_H
#define FITMAP_H
#include<math.h>
#include<iostream>
#include "transitmask.h"
using namespace std;
class Fitmap{
  /* nx is the number of bins, ny is the number of transits
   *
   *
   * */
  public:
    Fitmap(double *cmap, int nx, int ny, bool pyflag=true);
    ~Fitmap();
    void FitSingleMask();
    void FitSingleFlat(int qint=10);
    void FitTTV(int nshift=0, int qint=10);
    void StdOutput();
    double Error();
    TransitMask *cmap_;
    //double *cmap_;
  private:
    void FitTran_(int qint, TransitMask &newtran, int &indmin,double &err);
    void SingleMask_(int q,int index,int indey);
    void SingleFlat_(int q,int indey);
    void AllTranMask_(int q);
    //void Shuffle_(int ns, int indey, const double* oldmap, const double* oldmask, double* newmap,double* newmask) const;
    //void CopyMask_(const double* oldmap, const double* oldmask, double* newmap,double* newmask) const;
    //void Chisquare_(const double *arrx, const double *arry, double & err);
    double *mask_,*maskoned_,*model_,*dipmask_;
    int nbin_;
    int ntr_;
    //int nonzeros_;
    double high_,low_;
};
#endif
