#ifndef BLSSIMPLE_H
#define BLSSIMPLE_H
#include <iostream>
#include <math.h>
using namespace std;
class BlsSpec{
  /* Create BLSspectrum for transits with width and epoch and depth 
   * information. 
   * */
  public:
    BlsSpec(double *time, int lt, double *mag, int lm, double pmax, double pmin, int fn, double epoch, double qvar, double depth,int bin = 200);

    ~BlsSpec();
    void GenSpec();
    //void BlsAnal(double q=0.04, int bin=0);
    //void OutputSpec(double *sr, int lbin);
    void StandOutput(int precision=7);

  private:
    void Foldts_();
    double CalSr_();
    int lt_;
    int lbin_;
    int fn_;
    double *time_;
    double *mag_;
    double *spec_;
    double *freqs_;
    double *magbin_;
    double *count_;
    double fmin_;
    double fmax_;
    double epoch_;
    double qvar_;
    double dip_;
    double period_;
};
#endif
