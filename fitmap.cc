#include"fitmap.h"
Fitmap::Fitmap(double *cmap, int nx, int ny)
  :nbin_(nx),ntr_(ny){
    nonzeros_=0; 
    cmap_=new double [nx*ny];
    mask_=new double [nx*ny];
    maskoned_=new double [nx];
    model_ = new double [nx*ny];    
    for (int i=0;i<nbin_;i++){
      maskoned_[i] = 0;
    }
    mean_ = 0;
    max_ = cmap_[0];
    for (int j=0;j<ntr_;j++){
      for (int i=0;i<nbin_;i++){
        model_[i+j*nbin_]=0;
        cmap_[i+j*nbin_] = cmap[j+i*ntr_]; //hack for python input,need to change later
        if(cmap_[i+j*nbin_] == 0){
          mask_[i+j*nbin_] = 0;
        } else {
          mask_[i+j*nbin_] = 1;
          maskoned_[i]+=1;
          nonzeros_+=1;
          mean_+=cmap_[i+j*nbin_];
          if(cmap_[i+j*nbin_]>max_) max_=cmap_[i+j*nbin_];
        }
      }
    }
    if(nonzeros_==0){
      cout<<"no valid values in the map"<<endl;
      exit(0);
    }
    mean_/=nonzeros_;
    high_= mean_;
    low_ = mean_;
}
Fitmap::~Fitmap(){
  delete [] cmap_;
  delete [] mask_;
  delete [] maskoned_;
  delete [] model_;
}
double Fitmap::Error(){
  double err=0;  
  Chisquare_(model_,cmap_,err);
  return err;
}
void Fitmap::FitSingleMask(){
  // assume that the outlier part is only high without low
  /* Step1 Let's decide which transit contributes to the signal first*/
  int index=nbin_,indey=ntr_,imin=0,imax=0;
  for (int i=0; i<nbin_/2; i++){
    if(maskoned_[(int)(nbin_/2)-i]==1){
      index=(int)(nbin_/2)-i;
      break;
    } else {
      if(maskoned_[(int)(nbin_/2)+i]==1){
        index=(int)(nbin_/2)+i;
        break;
      }
    }
  }
  if (index==nbin_ or index==0){
    cout<< "not be able to fit with a signle Mask model"<<endl;
    high_=0;
    low_=0;
    return;
  } 
  for (int j=0;j<ntr_;j++){
    if(mask_[index+j*nbin_]==1){
      indey=j;
      break;
    }
  }
  high_=0;
  imin = index;
  for (int i=index;i>0;i--){      
    if(mask_[i+indey*nbin_]==0){
      imin=i+1;
      break;
    }
  }
  imax = index+1;
  for (int i=index+1;i<nbin_;i++){
    if(mask_[i+indey*nbin_]==0){
      imax=i-1;
      break;
    }
  }
  index=(int)(nbin_/2);
  /* Step2 Now find the best model by varying q*/
  int qint = max((imax-index),(index-imin)),q, indmin=0; 
  double err=0,temperr=0,temphigh=0,templow=0;
  for (q=0; q<=qint; q++){
    int count=0;
    temphigh=0;
    for (int i=index-q;i<=index+q;i++){
      temphigh+=cmap_[i+indey*nbin_];
      count+=mask_[i+indey*nbin_];
    }
    temphigh/=count;
    templow = (nonzeros_*mean_-high_*count)/(nonzeros_-count);
    SingleMask_(q,index,indey,temphigh,templow);
    if(q==0){
      Chisquare_(model_,cmap_,err);
      temperr=err;
      high_=temphigh;
      low_=templow;
    } else {
      Chisquare_(model_,cmap_,temperr);
      if(temperr<err){
        err=temperr;
        indmin = q;
        high_=temphigh;
        low_=templow;
      }
    }
  }
  SingleMask_(indmin,index,indey);
  return;
}

void Fitmap::FitSingleFlat(int qint){
  //the model is one outlier plus flats in other transits
  int indey=0,indmin=0,index=(int)(nbin_/2),indeymin=0;
  double temphigh,templow;
  double err=0,temperr;
  double* meanoned = new double [ntr_];
  temphigh = 0;
  for (int q=0;q<=qint;q++){
    for (int j=0; j< ntr_;j++){
      meanoned[j]=0;
      int count=0;
      for (int i=index-q; i<=index+q; i++){
        meanoned[j]+=mask_[i+j*nbin_]*cmap_[i+j*nbin_];
        count+=mask_[i+j*nbin_];
      }
      if(count>0){
        meanoned[j]/=count;
      }
      if(meanoned[j]>temphigh){
        indey=j;
        temphigh=meanoned[j];
      }
    }
    templow = (nonzeros_*mean_-temphigh*(2*q+1))/(nonzeros_-(2*q+1));
    SingleFlat_(q,indey,temphigh,templow);
    if(q==0){
      Chisquare_(model_,cmap_,err);
      temperr=err;
      high_=temphigh;
      indeymin= indey;
      low_=templow;
    } else {
      Chisquare_(model_,cmap_,temperr);
      if(temperr<err){
        err=temperr;
        indmin = q;
        indeymin= indey;
        high_=temphigh;
        low_=templow;
      }
    }
 //   cout << q << " " << err << " " << indey <<endl;
  }
  SingleFlat_(indmin,indeymin);
  delete [] meanoned;
}
void Fitmap::FitTran_(int qint, const double *cmap, const double* mask, int &indmin,double &err){
  int q=0,index=(int)(nbin_/2)+1;
  double temphigh,templow,temperr;
  err =0;
  for (q=0; q<=qint; q++){
      int count=0;
      temphigh=0;
      for(int j=0;j<ntr_;j++){
        for (int i=index-q;i<=index+q;i++){
          temphigh+=cmap[i+j*nbin_];
          count+=mask[i+j*nbin_];
        }
      }
      temphigh/=count;
      templow = (nonzeros_*mean_-temphigh*count)/(nonzeros_-count);
      TransitMask_(q,temphigh,templow);
      if(q==0){
        Chisquare_(model_,cmap,err);
        temperr=err;
        high_=temphigh;
        low_=templow;
      } else {
        Chisquare_(model_,cmap,temperr);
        if(temperr<err){
          err=temperr;
          indmin = q;
          high_=temphigh;
          low_=templow;
        }
      }
      //cout << q << " " << temperr << endl;
  }
  TransitMask_(indmin);
}
void Fitmap::FitTTV(int nshift, int qint){
  //nshift=0 is the same as not allow TTV
  int q=0,indmin=0,nsmin=0;
  double err=0,temperr=0;
  FitTran_(qint,cmap_,mask_,indmin,err);
  if(nshift>0){
    double *newmap = new double [nbin_*ntr_];
    double *newmask = new double [nbin_*ntr_];
    double *tempmap = new double [nbin_*ntr_];
    double *tempmask = new double [nbin_*ntr_];
    CopyMask_(cmap_,mask_,newmap,newmask);
    CopyMask_(cmap_,mask_,tempmap,tempmask);
    for (int j=1;j<ntr_;j++){
      for (int ns=1;ns<=nshift; ns++){
        Shuffle_(ns,j,newmap,newmask,tempmap,tempmask);
        //if(j==4 && ns ==3){
        //for (int k=0;k<ntr_;k++){
        //  for (int i=0;i<nbin_;i++){
        //   cout << tempmask[i+k*nbin_] << " "; 
        //  }
        //  cout << "\n" ;
        //}
        //}
        FitTran_(qint,tempmap,tempmask,q,temperr); 
        if(j==4 && ns ==3){
        for (int k=0;k<ntr_;k++){
          for (int i=0;i<nbin_;i++){
           cout << model_[i+k*nbin_] << " "; 
          }
          cout << "\n" ;
        }
        }
        if(temperr<err){
          err=temperr;
          indmin = q;
          nsmin = ns;
        }
        Shuffle_(-ns,j,newmap,newmask,tempmap,tempmask);
        FitTran_(qint,tempmap,tempmask,q,temperr);
        if(temperr<err){
          err=temperr;
          indmin = q;
          nsmin = -ns;
        }
        //cout << temperr << " " << q << " " << ns << " " << j << " " << high_ << " " << low_<< endl;
      }
      Shuffle_(nsmin,j,newmap,newmask,tempmap,tempmask);
      CopyMask_(tempmap,tempmask,newmap,newmask);
    }
    FitTran_(qint,newmap,newmask,indmin,temperr);
    //cout<< indmin << " " << nsmin << endl; 
    delete [] newmap;
    delete [] newmask;
    delete [] tempmap;
    delete [] tempmask;
  }
  TransitMask_(indmin);
  return;
}

void Fitmap::Shuffle_(int ns, int indey,const double* oldmap,const double* oldmask,double *newmap ,double *newmask ) const {
  //do nothing
  int newx;
  for (int j=0;j<ntr_;j++){
    for(int i=0;i<nbin_;i++){
      if(j==indey){
        if((i+ns)<0){
          newx = nbin_+i+ns;
        } else {
          if((i+ns)>=nbin_){
            newx = (i+ns)-nbin_;
          } else {
          newx = i+ns;
          }
        }
        newmap[i+j*nbin_] = oldmap[newx+j*nbin_];
        newmask[i+j*nbin_] = oldmask[newx+j*nbin_];
      }else{     
        newmap[i+j*nbin_] = oldmap[i+j*nbin_];
        newmask[i+j*nbin_] = oldmask[i+j*nbin_];
      }
    }
  }
  return;
}
void Fitmap::CopyMask_(const double* oldmap,const double* oldmask,double *newmap ,double *newmask) const{
    for (int j=0;j<ntr_;j++){
      for(int i=0;i<nbin_;i++){
        newmap[i+j*nbin_] = oldmap[i+j*nbin_];
        newmask[i+j*nbin_] = oldmask[i+j*nbin_];
      }
    }
}
void Fitmap::StdOutput(){
 
  for (int j=0;j<ntr_;j++){
    for (int i=0;i<nbin_;i++){
      cout << model_[i+j*nbin_] << " "; 
    }
    cout << " " << endl; 
  }
}
void Fitmap::SingleMask_(int q,int index,int indey, double high, double low){
  if(high==0){
    high=high_;
  }
  if(low==0){
    low=low_;
  }
  for (int j=0;j<ntr_;j++){
         for (int i=0;i<nbin_;i++){
           if(j==indey){
             if(((index-q)<=i) && (i<=(index+q))){
               model_[i+j*nbin_]=high*mask_[i+j*nbin_];
             } else {
               model_[i+j*nbin_]=mask_[i+j*nbin_]*low;
             }
           } else {
               model_[i+j*nbin_]=mask_[i+j*nbin_]*low;
           }
         }
      }
}

void Fitmap::SingleFlat_(int q, int indey, double high, double low){
  if(high==0){
    high=high_;
  }
  if(low==0){
    low=low_;
  }
  int index=(int)(nbin_/2),imin=index-q,imax=index+q;
  for (int j=0;j<ntr_;j++){
    for (int i=0;i<nbin_;i++){
      if(j==indey){
        if((imin<=i) && (i<=imax)){
          model_[i+j*nbin_]=high*mask_[i+j*nbin_];
        } else {
          model_[i+j*nbin_]=mask_[i+j*nbin_]*low;
        }
      } else {
        model_[i+j*nbin_]=mask_[i+j*nbin_]*low;
      }
    }
  }
}

void Fitmap::TransitMask_(int q, double high, double low){
  if(high==0){
    high=high_;
  }
  if(low==0){
    low=low_;
  }
  int index=(int)(nbin_/2);
  for (int j=0;j<ntr_;j++){
     for (int i=0;i<nbin_;i++){
        if(((index-q)<=i) && (i<=(index+q))){
           model_[i+j*nbin_]=high*mask_[i+j*nbin_];
        } else {
           model_[i+j*nbin_]=mask_[i+j*nbin_]*low;
        }
    }
  }
}
void Fitmap::Chisquare_(const double* arrx, const double* arry, double &err){
  err=0;
  for (int j=0;j<ntr_;j++){
    for (int i=0;i<nbin_;i++){
      err+=pow((arrx[i+j*nbin_]-arry[i+j*nbin_]),2.0);
    }
  }
  err/=(ntr_*nbin_);
  err = sqrt(err);
}
