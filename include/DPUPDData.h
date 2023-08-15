#ifndef DPUPDDATA_H
#define DPUPDDATA_H

#include <boost/multi_array.hpp>
#include <vector>
#include <DPConst.h>
#include <DPDefs.h>

namespace DP {

// Object used to store contents of a UPD file
class UPDData { // TODO: move to an appropriate location
public:
  enum upd_const_t {
    UPD_MALE = 1,
    UPD_FEMALE = 2,
    UPD_YEAR_START = 1970,
    UPD_YEAR_FINAL = 2049
  };
  
  static const int upd_years_basepop = 4;
  static const int upd_years_overall = 80; // This must equal (upd_year_final - upd_year_start + 1)
   
  // Constructors and destructors
  UPDData();
  ~UPDData();
  
  // read the contents of a UPD file
  // return 0 on success, <0 on failure
  int read(const std::string &upd_filename);
  
  // Accessors. In accessors, year t=0 denotes 1970
  inline double basepop(const int t, const int s, const int a) const {return _basepop[t][s][a];}
  
  inline double lx(const int t, const int s, const int a) const {return _lx[t][s][a];}
  inline double ex(const int t, const int s, const int a) const {return _ex[t][s][a];}
  inline double Sx(const int t, const int s, const int a) const {return _Sx[t][s][a];}
  
  inline double tfr(const int t) const {return _tfr[t];}
  inline double srb(const int t) const {return _srb[t];}
  
  // Access using age 15 <= a < 50
  inline double pasfrs(const int t, const int a) const {return _pasfrs[t][a - DP::AGE_BIRTH_MIN];}
  
  inline double migration(const int t, const int s, const int a) {return _migration[t][s][a];}
  
private:
  int read_basepop(std::ifstream &fin);
  int read_lfts(std::ifstream &fin);
  int read_series(std::ifstream &fin, time_series_t &series, const std::string &end_tag);
  int read_pasfrs(std::ifstream &fin);
  int read_migration(std::ifstream &fin);
  
  static const char _delim = ',';
  
  // Note: _lx, _ex, and _Sx are indexed by 0..81, with 81 denoting ages >80,
  // whereas other age-indexed arrays just use 0..80 with 80 denoting ages >=80
  year_sex_age_t _basepop;
  year_sex_age_t _lx;
  year_sex_age_t _ex;
  year_sex_age_t _Sx;
  time_series_t  _tfr;
  time_series_t  _srb;
  year_age_t     _pasfrs;
  year_sex_age_t _migration;
};

} // end namespace DP

#include <DPUPDData_impl.h>

#endif // DPUPDDATA_H