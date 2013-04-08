#include "maskarray.h"
MaskArray::MaskArray(const double *M, int nx, int ny)
  : nx_(nx),ny_(ny){
    value_=new double [nx_*ny_];
    mask_=new double [nx_*ny_];
    onedmeanx_=new double [ny_];
    onedmeany_=new double [nx_];
    for (int i=0;i<nx_;i++){
      onedmeany_[i] = 0;
    }
    for (int i=0;i<ny_;i++){
      onedmeanx_[i] = 0;
    }
    mean_ = 0;
    nonzeros_ =0;
    for (int j=0;j<ny_;j++){
      for (int i=0;i<nx_;i++){
        value_[i+j*nx_] = M[i+j*nx_]; //hack for python input,need to change later
        if(value_[i+j*nx_] == 0){
          mask_[i+j*nx_] = 0;
        } else {
          mask_[i+j*nx_] = 1;
          nonzeros_+=1;
          mean_+=value_[i+j*nx_];
        }
      //   cout << mask_[i+j*nx_] << " ";
      }
     // cout << "\n";
    }
    if(nonzeros_==0){
      cout<<"no valid values in the map"<<endl;
      exit(0);
    }
    mean_/=nonzeros_;
    return;
}
MaskArray::MaskArray(const MaskArray &Ma){
  double mean;
  Ma.Dimen(nx_,ny_);
  if(value_){
    delete [] value_;
    value_=new double [nx_*ny_];
  }
  if(mask_){
    delete [] mask_;
    mask_=new double [nx_*ny_];
  }
  if(onedmeanx_){
    delete [] onedmeanx_;
    onedmeanx_=new double [ny_];
  }
  if(onedmeany_){
    delete [] onedmeany_;
    onedmeany_=new double [nx_];
  }
  for (int i=0;i<nx_;i++){
      onedmeany_[i] = 0;
  }
  for (int i=0;i<ny_;i++){
    onedmeanx_[i] = 0;
  }

  Ma.Mean(mean);
  mean_ = mean;
  Ma.View(value_);
  Ma.Mask(mask_);
  return;
}
MaskArray::~MaskArray(){ 
  delete [] value_;
  delete [] mask_;
  delete [] onedmeanx_;
  delete [] onedmeany_;
}
void MaskArray::Copy(double *value, double *mask, int nx, int ny){
  if(nx!=nx_ || ny!=ny_){
    cout << "dimension not match!" << endl;
    exit(1);
  }
  nonzeros_=0;
  for (int j=0;j<ny_;j++){
    for (int i=0;i<nx_;i++){
      value_[i+j*nx_] = value[i+j*nx_]; //hack for python input,need to change later
      mask_[i+j*nx_] = mask[i+j*nx_];
      nonzeros_+=mask_[i+j*nx_];
      mean_+=value_[i+j*nx_]*mask_[i+j*nx_];
    }
  }
  mean_/=nonzeros_;
}
void MaskArray::Dimen(int nx, int ny) const{
  nx = nx_;
  ny = ny_;
}
void MaskArray::View(double *M) const{
    for (int j=0;j<ny_;j++){
      for (int i=0;i<nx_;i++){
        M[i+j*nx_] = value_[i+j*nx_]; //hack for python input,need to change later
      }
    }
}
void MaskArray::Mask(double *M) const{
    for (int j=0;j<ny_;j++){
      for (int i=0;i<nx_;i++){
        M[i+j*nx_] = mask_[i+j*nx_]; //hack for python input,need to change later
      }
    }
}

void MaskArray::MaskView(double *M) const{
    for (int j=0;j<ny_;j++){
      for (int i=0;i<nx_;i++){
        M[i+j*nx_] = mask_[i+j*nx_]*M[i+j*nx_]; //hack for python input,need to change later
      }
    }
}
void MaskArray::Mean(double &mean) const{
  mean = mean_;
  return;
}
void MaskArray::OneDMean(int dim, double *onemean){
  double* countx_=new double [ny_];
  double* county_=new double [nx_];
  
  for (int j=0;j<ny_;j++){
    for (int i=0;i<nx_;i++){
      onedmeanx_[j] += value_[i+j*nx_];
      countx_[j] += mask_[i+j*nx_];
      onedmeany_[i] += value_[i+j*nx_];
      county_[i] += mask_[i+j*nx_];
    }
  }
  for (int j=0;j<ny_;j++){
    if(countx_[j]!=0){
      onedmeanx_[j]/=countx_[j];
      if(dim==0){
        onemean[j]=onedmeanx_[j];
      }
    }
  }
  for (int i=0;i<nx_;i++){
    if(county_[i]!=0){
      onedmeany_[i]/=county_[i];
      if(dim==1){
        onemean[i]=onedmeany_[i];
      }
    }
  }
  delete [] countx_;
  delete [] county_;
}
void MaskArray::OneDMask(int dim, double *onemean){
  double* countx_=new double [ny_];
  double* county_=new double [nx_];
  
  for (int j=0;j<ny_;j++){
    for (int i=0;i<nx_;i++){
      countx_[j] += mask_[i+j*nx_];
      county_[i] += mask_[i+j*nx_];
    }
  }
  if(dim==0){
    for (int j=0;j<ny_;j++){
      onemean[j]=countx_[j];
    }
  }
 
  if(dim==1){
    for (int i=0;i<nx_;i++){
      onemean[i]=county_[i];
    }
  }
  delete [] countx_;
  delete [] county_;
}

void MaskArray::Shuffle(int step, int indey,MaskArray *Ma) const{
  //do nothing
  int newx;
  
  double* newvalue=new double [nx_*ny_];
  double* newmask=new double [nx_*ny_];
  for (int j=0;j<ny_;j++){
    for(int i=0;i<nx_;i++){
      if(j==indey){
        if((i+step)<0){
          newx = nx_+i+step;
        } else {
          if((i+step)>=nx_){
            newx = (i+step)-nx_;
          } else {
          newx = i+step;
          }
        }
        newvalue[i+j*nx_] = value_[newx+j*nx_];
        newmask[i+j*nx_] = mask_[newx+j*nx_];
      }else{     
        newvalue[i+j*nx_] = value_[i+j*nx_];
        newmask[i+j*nx_] = mask_[i+j*nx_];
      }
    }
  }
  Ma->Copy(newvalue,newmask,nx_,ny_);
  delete [] newvalue;
  delete [] newmask;
  return;
}
void MaskArray::StandOutput() const{ 
  for (int j=0;j<ny_;j++){
    for (int i=0;i<nx_;i++){
      cout << value_[i+j*nx_] << " "; 
    }
    cout << " " << endl; 
  }

}
