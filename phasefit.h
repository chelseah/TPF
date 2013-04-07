#ifndef PHASEFIT_H
#define PHASEFIT_H
#include <iostream>
using namespace std;
class PhaseLc{
  /* Phase fold a lightcurve use given information.
   * And create a color map of transit (nbin(intran)\timesntran), the width 
   * of intransit is given by user.
   * can ask for finer bins in the intran part.
   * */
  public:
    PhaseLc(double *time, int lt, double *mag, int lm, double period, double epoch, int bin = 200);
    ~PhaseLc();
    void Foldts();
    void TransitColor(double q=0.04, int bin=0);
    double *time_;
    double *mag_;
    double *phase_;
    double *magbin_;
  private:
    double period_;
    double epoch_;
    int lbin_;
    int lt_;
    int ntran_;
};
#endif
