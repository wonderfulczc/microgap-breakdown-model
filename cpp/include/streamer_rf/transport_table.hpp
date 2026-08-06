#pragma once
#include <algorithm>
#include <stdexcept>
#include <vector>
namespace streamer_rf {
class TransportTable {
 public:
  TransportTable(std::vector<double> field,std::vector<double> value):x_(std::move(field)),y_(std::move(value)){
    if(x_.size()<2||x_.size()!=y_.size()) throw std::invalid_argument("invalid transport table");
  }
  double interpolate(double x) const {
    if(x<x_.front()||x>x_.back()) throw std::out_of_range("E/N outside transport table");
    auto k=std::upper_bound(x_.begin(),x_.end(),x)-x_.begin(); if(k==0)return y_[0]; if(k==x_.size())return y_.back();
    const double t=(x-x_[k-1])/(x_[k]-x_[k-1]); return y_[k-1]+t*(y_[k]-y_[k-1]);
  }
 private: std::vector<double>x_,y_;
};
}
