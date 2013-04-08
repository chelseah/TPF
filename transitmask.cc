#include "transitmask.h"
TransitMask :: TransitMask(const double *cmap, int nbin, int ntr,double high,double low) : MaskArray(cmap, nbin, ntr),high_(high),low_(low) {
  transit_ = new double [nx_*ny_];
  for (int j=0; j<ny_; j++){
    for (int i=0; i<nx_; i++){
      transit_[i+j*nx_]=0;
    }
  }
}
TransitMask :: TransitMask(const TransitMask &TM):MaskArray(TM){
  transit_ = new double [nx_*ny_];
  TM.OutputTransit(transit_,high_,low_);
}
TransitMask ::~TransitMask(){
  delete [] transit_;
}
void TransitMask::TransitModel(const double *Ma){
  double temphigh=0;
  int count = 0;
  for (int j=0; j<ny_; j++){
    for (int i=0; i<nx_; i++){
      temphigh += value_[i+j*nx_]*Ma[i+j*nx_];
      count += mask_[i+j*nx_]*Ma[i+j*nx_];
    }
  }
  high_ = temphigh/count;
  low_ = (nonzeros_*mean_-high_*count)/(nonzeros_-count);

  for (int j=0; j<ny_; j++){
    for (int i=0; i<nx_; i++){
      transit_[i+j*nx_] = (Ma[i+j*nx_]*high_+(1-Ma[i+j*nx_])*low_)*mask_[i+j*nx_]; 
    }
  }
}
void TransitMask::OutputTransit(double *M, double &high, double &low) const{
  for (int j=0; j<ny_; j++){
    for (int i=0; i<nx_; i++){
      M[i+j*nx_] = transit_[i+j*nx_]; 
    }
  }
  high = high_;
  low = low_;
}
void TransitMask::PrintTransit() const{
  for (int j=0;j<ny_;j++){
    for (int i=0;i<nx_;i++){
      cout << transit_[i+j*nx_] << " "; 
    }
    cout << " " << endl; 
  }
}
void TransitMask::Chisquare(const double* model, double &err) const{
  err=0;
  for (int j=0;j<ny_;j++){
    for (int i=0;i<nx_;i++){
      //cout << value_[i+j*nx_] << " " << model[i+j*nx_] << endl;
      err+=pow((value_[i+j*nx_]-model[i+j*nx_]),2.0);
    }
  }
  err/=(nx_*ny_);
  err = sqrt(err);
  return;
}
TransitMask &TransitMask::operator = (const TransitMask &Ma){
  if(this==&Ma){
    return *this;
  }
  Ma.Dimen(nx_,ny_);
  double* value = new double [nx_*ny_];
  double* mask = new double [nx_*ny_];
  Ma.View(value);
  Ma.Mask(mask);
  Copy(value,mask,nx_,ny_);
  high_=0;
  low_=0;
  delete [] value;
  delete [] mask;
}
