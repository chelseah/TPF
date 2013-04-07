#ifndef MASKARRAY_H
#define MASKARRAY_H
#include "iostream"
class MaskArray {
  public:
    MaskArray(double *M, int nx, int ny);
    MaskArray(const MaskArray *Ma);
    ~MaskArray();
    double Mean();
    void Dimen(int nx, int ny) const;
    void View(double *M) const; //return a view of value;
    void Mask(double *M) const; //return a view of mask;
    void MaskView(double *M) const; //return a masked riew of M;
    void Shuffle(int step, int indey, MaskArray *Ma)const;
    void OneDMean(int dim, double *onemean);
    void Copy(double *value, double *mask);
    void StandOutput() const;
  private:
    double *value_;
    double *mask_;
    double mean_;
    double *onedmeanx_;
    double *onedmeany_;
    int nx_;
    int ny_;
}
#endif
