#include "phasefit.h"
PhaseLc::PhaseLc(double* time, int lt, double* mag, int lm, double period, double epoch, int bin)
  :period_(period),epoch_(epoch),lbin_(bin),lt_(lt){
  /*allocate meomery for time and mag*/
  time_ = new double [lt];
  mag_ = new double [lm];
  color_= NULL; /*Do not allocate space for color_ is the color map function 
                  is not called*/
  intran_= 0;   /*The number of bins in transit, need to be figured out 
                  by later settings*/
  double tmax=0; /*In case time is not sorted, might have trouble with 
                   tmin as well, potential bug.*/
  for (int i = 0; i < lt; i++){
    if(time[i]>tmax){
      tmax=time[i];
    }
    time_[i] = time[i];
    mag_[i] = mag[i];
  }
  phase_ = new double [lbin_]; 
  magbin_ = new double [lbin_];
  for (int i=0; i<lbin_; i++){
    phase_[i]=1.0/(lbin_-1)*i; //right now this seems useless, but later on
                              //we will want onevenly binned phases
    magbin_[i]=0;
  }
  ntran_ = (tmax-epoch-0.5*period)/period+1; /*compute the number of the 
                                               transits*/
//  cout << ntran_<< " " << lbin_ << " " <<(int)((time_[lt-1]-epoch-0.5*period)/period_)<< " " << time_[lt-1] << endl;

  }

PhaseLc::~PhaseLc(){
  /*free space*/
  if(color_){ /*since color_ is not necessay non zeros, so 
               need to justify here.*/
    delete [] color_;
  }
  delete [] time_;
  delete [] mag_;
  delete [] phase_;
  delete [] magbin_;
}

void PhaseLc::Foldts(){
  double f = 1.0/period_,ph;
  /*temperary stor the number of data points in each bin for later on 
    compute the mean*/
  double *count=new double [lbin_]; 
  int k;
  for (k=0; k<lbin_; k++){
    count[k]=0;
  }
  /*standard phase folding with transit at phase=0.5*/
  for (int i=0; i< lt_; i++){
    ph = (time_[i]-epoch_-0.5*period_)*f-(int)((time_[i]-epoch_-0.5*period_)*f);
    if(ph<0) ph+=1;
    k = (int)(lbin_*ph);
    magbin_[k]+=mag_[i];
    count[k]+=1;
  }
  for (k=0; k<lbin_; k++){
    /*do not compute the mean for the gap*/
    if(count[k]!=0) {
      magbin_[k]/=count[k];
    }
  }
  delete [] count;
  return;
}

void PhaseLc::TransitColor(double q, int bin){
  double f = 1.0/period_,ph;
  if(bin==0) bin=lbin_;
  /*compute the boundary for the region to zoom in(or out?)*/
  int is = (int)(bin*(0.5-q)); 
  int ie = (int)(bin*(0.5+q));
  int j,k,indextran, tempindex=0;
  intran_=ie-is+1; /*compute the width.*/
  if (color_){
    delete [] color_; /*free the old space to avoid leaking.*/
  }
  color_=new double [intran_*ntran_];
  /*temporary space for storing number counts.*/
  double * count = new double [intran_*ntran_];
  for (j=0;j<ntran_;j++){
    for (k=0; k<intran_; k++){
      color_[k+j*intran_]=0;
      count[k+j*intran_]=0;
    }
  }

  for (int i=0; i< lt_; i++){
    /*the index of the transit*/
    indextran=(int)((time_[i]-epoch_-0.5*period_)*f);
    ph = (time_[i]-epoch_-0.5*period_)*f-indextran;
    if(ph<0) ph+=1;
    k = (int)(bin*ph);
    if((is<=k) && (k<=ie)){
      color_[(k-is)+indextran*intran_]+=mag_[i];    
      count[(k-is)+indextran*intran_]+=1;   
      if(((k-is)+indextran*intran_)>tempindex){
        tempindex=(k-is)+indextran*intran_; /*is become 0*/
      }
    }
  }

  for (j=0;j<ntran_;j++){
    for (k=0; k<intran_;k++){
      if(count[k+j*intran_]!=0){
        color_[k+j*intran_] = color_[k+j*intran_]/count[k+j*intran_]; 
      } 
    }
  }
  delete [] count;
  return;
}

void PhaseLc::OutputPhase(double *magbin, int lbin){
  for (int i=0; i<lbin_; i++){
    magbin[i] = magbin_[i];
  }
  return;
}
void PhaseLc::OutputColor(double *color, int nbin, int ntran){
  for (int j=0; j<ntran_; j++){
    for (int i=0; i<intran_; i++){
      color[i+j*intran_] = color_[i+j*intran_];
    }
  }
  return;
}
void PhaseLc::StandOutputPhase(int precision){
  cout.precision(precision);
  cout << "#PHASE MAG" <<endl;
  cout << "#[1] [2]";
  for (int i=0;i<lbin_; i++){
    cout << phase_[i] << " " << magbin_[i]  << endl; 
  }
  return;
}

void PhaseLc::StandOutputColor(int precision){
  cout.precision(precision);
  for (int j=0; j<ntran_; j++){
    for (int i=0; i<intran_; i++){
      cout << color_[i+j*intran_] << " ";
    }
      cout << "\n"; 
  }
  return;
}
