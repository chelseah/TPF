#include"fitmap.h"
Fitmap::Fitmap(double *cmap, int nx, int ny)
  :nbin_(nx),ntr_(ny){
    double *pycmap = new double [nx*ny];//hack python
    for (int j=0;j<ntr_;j++){
      for (int i=0; i<nbin_;i++){
        pycmap[i+j*nbin_] = cmap[j+i*ntr_];
      }
    }
    cmap_=new TransitMask(pycmap, nx, ny,0.0,0.0);
    model_ = new double [nx*ny];   
    mask_ = new double [nx*ny];  
    dipmask_ = new double [nx*ny];  
    maskoned_ = new double [nx];
    cmap_->OneDMask(1,maskoned_);
    cmap_->Mask(mask_);

    delete [] pycmap;
}
Fitmap::~Fitmap(){
  delete cmap_;
  delete [] model_;
  delete [] maskoned_;
  delete [] mask_;
}
double Fitmap::Error(){
  double err=0;  
  cmap_->Chisquare(model_,err);
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
    SingleMask_(q,index,indey);
    cmap_->TransitModel(dipmask_);
    cmap_->OutputTransit(model_,temphigh, templow);
     if(q==0){
      cmap_->Chisquare(model_,err);
      temperr=err;
      high_=temphigh;
      low_=templow;
    } else {
      cmap_->Chisquare(model_,temperr);
      if(temperr<err){
        err=temperr;
        indmin = q;
        high_=temphigh;
        low_=templow;
      }
    }
  }
  SingleMask_(indmin,index,indey);
  cmap_->TransitModel(dipmask_);
  cmap_->OutputTransit(model_,temphigh, templow);
  return;
}

void Fitmap::FitSingleFlat(int qint){
  //the model is one outlier plus flats in other transits
  int indey=0,indmin=0,indeymin=0;
  double temphigh,templow,meanoned;
  double err=0,temperr;
  temphigh = 0;
  for (int q=0;q<=qint;q++){
    for (int j=0; j< ntr_;j++){
      SingleFlat_(q,j);
      cmap_->TransitModel(dipmask_);
      cmap_->OutputTransit(model_,meanoned, templow);
      if(meanoned>temphigh){
        indey=j;
        temphigh=meanoned;
      }
    }
    SingleFlat_(q,indey);
    cmap_->TransitModel(dipmask_);
    cmap_->OutputTransit(model_,temphigh, templow);
    if(q==0){
      cmap_->Chisquare(model_,err);
      temperr=err;
      high_=temphigh;
      indeymin= indey;
      low_=templow;
    } else {
      cmap_->Chisquare(model_,temperr);
      if(temperr<err){
        err=temperr;
        indmin = q;
        indeymin= indey;
        high_=temphigh;
        low_=templow;
      }
    }
  }
  SingleFlat_(indmin,indeymin);
  cmap_->TransitModel(dipmask_);
  cmap_->OutputTransit(model_,temphigh, templow);
}
void Fitmap::FitTran_(int qint, TransitMask &newtran, int &indmin,double &err){
  int q=0;
  double temphigh,templow,temperr;
  err =0;
  for (q=0; q<=qint; q++){
      AllTranMask_(q);
      newtran.TransitModel(dipmask_);
      newtran.OutputTransit(model_,temphigh, templow);
      if(q==0){
        newtran.Chisquare(model_,err);
        temperr=err;
        high_=temphigh;
        low_=templow;
      } else {
        newtran.Chisquare(model_,temperr);
        if(temperr<err){
          err=temperr;
          indmin = q;
          high_=temphigh;
          low_=templow;
        }
      }
  }
  AllTranMask_(indmin);
  newtran.TransitModel(dipmask_);
  newtran.OutputTransit(model_,temphigh, templow);
  return;
}
void Fitmap::FitTTV(int nshift, int qint){
  //nshift=0 is the same as not allow TTV
  int q=0,indmin=0,nsmin=0;
  double err=0,temperr=0;
  FitTran_(qint,*cmap_,indmin,err);
  if(nshift>0){
    TransitMask newtran = TransitMask(*cmap_);
    TransitMask temptran = TransitMask(*cmap_);
    double newmean;
    newtran.Mean(newmean);
    temptran.Mean(newmean);
    for (int j=1;j<ntr_;j++){
      for (int ns=1;ns<=nshift; ns++){
        temptran.Mean(newmean);
        FitTran_(qint,temptran,q,temperr); 
        newtran.Shuffle(ns,j,&temptran);
        if(temperr<err){
          err=temperr;
          indmin = q;
          nsmin = ns;
        }
        newtran.Shuffle(-ns,j,&temptran);
        FitTran_(qint,temptran,q,temperr);
        if(temperr<err){
          err=temperr;
          indmin = q;
          nsmin = -ns;
        }
      }
        newtran.Shuffle(nsmin,j,&temptran);
        newtran = temptran;
    }
    FitTran_(qint,newtran,indmin,temperr);

  }
  return;
}

void Fitmap::StdOutput(){
 
  for (int j=0;j<ntr_;j++){
    for (int i=0;i<nbin_;i++){
      cout << model_[i+j*nbin_] << " "; 
    }
    cout << " " << endl; 
  }
}
void Fitmap::SingleMask_(int q,int index,int indey){
  for (int j=0;j<ntr_;j++){
    for (int i=0;i<nbin_;i++){
      if(j==indey){
        if(((index-q)<=i) && (i<=(index+q))){
          dipmask_[i+j*nbin_]=1;
         } else {
          dipmask_[i+j*nbin_]=0;
         }
      } else {
         dipmask_[i+j*nbin_]=0;
      }
    }
  }
}

void Fitmap::SingleFlat_(int q, int indey){
  int index=(int)(nbin_/2),imin=index-q,imax=index+q;
  for (int j=0;j<ntr_;j++){
    for (int i=0;i<nbin_;i++){
      if(j==indey){
        if((imin<=i) && (i<=imax)){
          dipmask_[i+j*nbin_]=1;
        } else {
          dipmask_[i+j*nbin_]=0;
        }
      } else {
          dipmask_[i+j*nbin_]=0;
      }
    }
  }
}

void Fitmap::AllTranMask_(int q){
  int index=(int)(nbin_/2);
  for (int j=0;j<ntr_;j++){
     for (int i=0;i<nbin_;i++){
        if(((index-q)<=i) && (i<=(index+q))){
          dipmask_[i+j*nbin_] = 1;
        } else {
          dipmask_[i+j*nbin_] = 0;
        }
    }
  }
}

