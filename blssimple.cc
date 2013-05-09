#include "blssimple.h"

BlsSpec::BlsSpec(double *time, int lt, double *mag, int lm, double pmax, double pmin, int fn, double epoch, double qvar, double depth, int bin)
  :lt_(lt),fmin_(1./pmax),fmax_(1./pmin),fn_(fn),epoch_(epoch),qvar_(qvar),dip_(depth),lbin_(bin){
  /*allocate meomery for time and mag*/
  time_ = new double [lt];
  mag_ = new double [lm];
  freqs_ = new double [fn];
  spec_ = new double [fn];
  period_ = -1;
  double df = (fmax_-fmin_)/(fn_-1);
  double avmag = 0;
  for (int i = 0; i < lt; i++){
    time_[i] = time[i];
    mag_[i] = mag[i];
    avmag+=mag_[i]/lt_;
  }
  for (int i=0; i<lt;i++){
    mag_[i]-=avmag;
  }
  for (int i=0; i<fn_; i++){
    freqs_[i] = fmin_+df*i;
    spec_[i] = 0;
  }
  magbin_ = new double [lbin_];
  count_ = new double [lbin_];
  for (int i=0; i<lbin_; i++){
    magbin_[i]=0;
    count_[i]=0;
  }
}
BlsSpec::~BlsSpec(){
  /*free space*/
  delete [] time_;
  delete [] mag_;
  delete [] freqs_;
  delete [] spec_;
  delete [] magbin_;
  delete [] count_;
}

void BlsSpec::GenSpec(){
  double sr;
  for (int i=0; i<fn_; i++){
   // int i=78573;
    cout << freqs_[i]<< endl;
    period_ = 1.0/freqs_[i];
    Foldts_();
    
    sr = CalSr_();
    spec_[i] = sr;
  }
}

double BlsSpec::CalSr_(){
/*fit it with flat OOTV and known depth, think about other ways later*/
  int nintran = (int)lbin_*qvar_; 
  int qi = (int)round(lbin_/2.-nintran/2.);
  int qe = qi+nintran-1;
  int count=0;
  double mout=0,min,err=0,s=0,r=0;
  //the bls statistic
  for (int i=qi; i<=qe; i++){
    s+=magbin_[i];
    r+=count_[i];
  }
  //cout << s << " " << r << endl;
  return s*s/(r*(1.0-r));

  //for(int i = 0; i<qi;i++){
  //  mout+=magbin_[i];
  //  if(magbin_[i]!=0){
  //    count+=1;
  //  }
  //}
  //for(int i = qe+1; i<lbin_;i++){
  //  mout+=magbin_[i];
  //  if(magbin_[i]!=0){
  //    count+=1;
  //  }
  //}
  //mout/=count;
  //min = mout+dip_;
  ////cout << mout << " " <<min << " " << qi << " " << qe << endl;
  //for(int i = 0; i<lbin_; i++){
  //  if(qi<=i && i<=qe){
  //    err+=fabs(magbin_[i]-min);
  //  } else {
  //    err+=fabs(magbin_[i]-mout);
  //  }
  //}
  //return (1.0/err*(mout*lbin_));
}

void BlsSpec::Foldts_(){
  double f = 1.0/period_,ph;
  /*temperary stor the number of data points in each bin for later on 
    compute the mean*/
  //double *count=new double [lbin_]; 
  int k;
  for (k=0; k<lbin_; k++){
    magbin_[k]=0;
    count_[k]=0;
  }
  /*standard phase folding with transit at phase=0.5*/
  for (int i=0; i< lt_; i++){
    ph = (time_[i]-epoch_-0.5*period_)*f-(int)((time_[i]-epoch_-0.5*period_)*f);
    if(ph<0) ph+=1;
    k = (int)(lbin_*ph);
    //magbin_[k]+=mag_[i];
    magbin_[k]+=mag_[i]/lt_; //bls-z
    //count[k]+=1;
    count_[k]+=1.0/lt_;  //bls-y
  }
//  for (k=0; k<lbin_; k++){
//    /*do not compute the mean for the gap*/
//    if(count[k]!=0) {
//      magbin_[k]/=count[k];
//     // cout << k << " " << magbin_[k] << endl;
//    }
//  }
//  delete [] count;
  return;
}

void BlsSpec::StandOutput(int precision){
  cout.precision(precision);
  cout << "#freq sr" <<endl;
  cout << "#[1] [2]" << endl;
  for (int i=0;i<fn_; i++){
    cout << freqs_[i] << " " << spec_[i]  << endl; 
  }
  return;
}

