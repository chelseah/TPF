#include "phasefit.h"
PhaseLc::PhaseLc(double* time, int lt, double* mag, int lm, double period, double epoch, int bin)
  :period_(period),epoch_(epoch),lbin_(bin),lt_(lt){
  time_ = new double [lt];
  mag_ = new double [lm];
  for (int i = 0; i < lt; i++){
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
  ntran_ = (time[lt-1]-epoch-0.5*period)/period+1;
  //cout << ntran_<< " " << lbin_ << " " <<(time_[lt-1]-epoch-0.5*period)/period_<< endl;
}

PhaseLc::~PhaseLc(){
  delete [] time_;
  delete [] mag_;
  delete [] phase_;
  delete [] magbin_;
}

void PhaseLc::Foldts(){
  double f = 1.0/period_,ph;
  double *count=new double [lbin_];
  int k;
  for (k=0; k<lbin_; k++){
    count[k]=0;
  }
  for (int i=0; i< lt_; i++){
    ph = (time_[i]-epoch_-0.5*period_)*f-(int)((time_[i]-epoch_-0.5*period_)*f);
    if(ph<0) ph+=1;
    k = (int)(lbin_*ph);
    magbin_[k]+=mag_[i];
    count[k]+=1;
  }
  for (k=0; k<lbin_; k++){
    if(count[k]!=0) {
      magbin_[k]/=count[k];
    //  cout<< phase_[k] << " " << magbin_[k]  << endl;
    }
  }
  delete [] count;
  return;
}

void PhaseLc::TransitColor(double q, int bin){
  double f = 1.0/period_,ph;
  if(bin==0) bin=lbin_;
  int is = (int)(bin*(0.5-q));
  int ie = (int)(bin*(0.5+q));
  int j,k,indextran, intran=ie-is+1;
  //cout<< is << " " << ie << " " << intran << " " << ntran_<<endl;
  double * color=new double [intran*ntran_];
  double * count = new double [intran*ntran_];
  for (j=0;j<ntran_;j++){
    for (k=0; k<intran; k++){
      color[k+j*intran]=0;
      count[k+j*intran]=0;
    }
  }

  for (int i=0; i< lt_; i++){
    indextran=(int)((time_[i]-epoch_-0.5*period_)*f);
    ph = (time_[i]-epoch_-0.5*period_)*f-indextran;
    if(ph<0) ph+=1;
    k = (int)(bin*ph);
    if((is<=k) && (k<=ie)){
      color[(k-is)+indextran*intran]+=mag_[i];    
      count[(k-is)+indextran*intran]+=1;   
    //  cout<< count[(k-is)+indextran*intran] << " " << k-is+indextran*intran << " " << indextran << endl;
    }
  }
  for (j=0;j<ntran_;j++){
    for (k=0; k<intran;k++){
      if(count[k+j*intran]!=0){
        color[k+j*intran] = color[k+j*intran]/count[k+j*intran]; 
      } 
        cout<<color[k+j*intran]<< " ";
    }
        cout<<" "<<endl;
  }
  delete [] color;
  delete [] count;
  return;
}
