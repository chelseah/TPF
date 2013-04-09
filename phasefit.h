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
    /*\brief A constructor which take the information from stric periodic 
     * search and initialize the class.
     *  
     *  Input: 
     *    time (double *): time series of the lightcurve.
     *    lt (int): the length of the time series. 
     *    double *time, int lt directly corresponding to python 1d numpy array.
     *    mag (double *): flux series of the lightcurve (in magnitude is 
     *      recommended).
     *    lm (int):  the length of the magnitude series. (need to be the same 
     *      with the time series).
     *    double *mag, int lm works the same as double *time, int lt when 
     *      face with python.
     *    period (double): the period of the lightcurve from stric periodic 
     *      search.
     *    epoch(double): the transit center of the 0th transit. 
     *    bin(int) (option) : the total number of bins used in the phase 
     *      space (from 0 to 1, with 0.5 corresponding to the center of transit)
     * */
    PhaseLc(double *time, int lt, double *mag, int lm, double period, double epoch, int bin = 200);

    ~PhaseLc();
    void Foldts();
    /* \brief Fold the lightcurve in phase space, the result is stored in 
     * the class instead of return to the user.
     * */
    void TransitColor(double q=0.04, int bin=0);
    /* \brief Fold the lightcurve in phase ntran space, and cut out the 
     * required width with the required bin size. The result is stored in 
     * the class instead of return to the user.
     *
     *  Input: 
     *    q (double [0,1]) (optional, default=0.04): the width in the 
     *      phase space need to zoom in.
     *    bin (int) (optional,default is the number of bins used when 
     *      initialized the class): the number of bins for the all phase 
     *      (from 0 to 1)to use in the zoomed in region. The actual number 
     *      of bins used in the zoomed in region in 2q*bin.
     * */
    void OutputPhase(double *magbin, int lbin);
    /* \brief Output the phase data to an array phase. 
     *
     *  Output:
     *    phase (double *): array containing the phase information.
     *    lp (int): the length of phase array.
     *    User is responsible for the create and free of phase. 
     *    phase and lp can be directly linked to a python 1d array.
     * */
    void OutputColor(double *color, int nbin,int ntran);
    /* \brief Output the color data to an array phase. 
     *
     *  Output:
     *    color (double *): array [nbin*ntran] containing the color
     *    information.
     *    nbin (int): the x-dim of the color array.
     *    ntran (int): the y-dim of the color array.
     *    User is responsible for the create and free of color. 
     *    phase and nbin,ntran can be directly linked to a python 2d array.
     * */

    void StandOutputPhase(int precision=7);
    /* \brief Output the phase data to standard IO
     *  
     *    Input: 
     *      precision (int) (optional, default=7): the number of precision 
     *      to use when output.
     *    Format:
     *      #Phase magbin
     *      #[1]    [2]  
     * */

    void StandOutputColor(int precision=7);
     /* \brief Output the color data to standard IO
     *  
     *    Input: 
     *      precision (int) (optional, default=7): the number of precision 
     *      to use when output. 
     * */

  private:
    double *time_;
    double *mag_;
    double *phase_;
    double *color_;
    double *magbin_;
    double period_;
    double epoch_;
    int lbin_;
    int lt_;
    int ntran_;
    int intran_;
};
#endif
