#ifndef DPUPDDATA_IMPL_H
#define DPUPDDATA_IMPL_H

#include <boost/multi_array.hpp>
#include <fstream>
#include <GBUtil.H>

namespace DP {

UPDData::UPDData()
  : _basepop(year_sex_age_t(boost::extents[upd_years_basepop][DP::N_SEX][DP::N_AGE])),
    _lx(year_sex_age_t(boost::extents[upd_years_overall][DP::N_SEX][DP::N_AGE + 1])),
    _ex(year_sex_age_t(boost::extents[upd_years_overall][DP::N_SEX][DP::N_AGE + 1])),
    _Sx(year_sex_age_t(boost::extents[upd_years_overall][DP::N_SEX][DP::N_AGE + 1])),
    _tfr(time_series_t(upd_years_overall)),
    _srb(time_series_t(upd_years_overall)),
    _pasfrs(year_age_t(boost::extents[upd_years_overall][DP::N_AGE_BIRTH])),
    _migration(year_sex_age_t(boost::extents[upd_years_overall][DP::N_SEX][DP::N_AGE])) {
}

UPDData::~UPDData() {}

int UPDData::read(const std::string &upd_filename) {
  int errcode(0);
  std::string line_buffer;
  
  // Read the UPD file contents into memory, one string per line
  std::ifstream fin(upd_filename.c_str());
  if (fin.good()) {
    while (std::getline(fin, line_buffer)) {
      if (GB::contains(line_buffer, "<basepop>"))        errcode = read_basepop(fin);
      else if (GB::contains(line_buffer, "<lfts>"))      errcode = read_lfts(fin);
      else if (GB::contains(line_buffer, "<tfr>"))       errcode = read_series(fin, _tfr, "</tfr>");
      else if (GB::contains(line_buffer, "<srb>"))       errcode = read_series(fin, _srb, "</srb>");
      else if (GB::contains(line_buffer, "<pasfrs>"))    errcode = read_pasfrs(fin);
      else if (GB::contains(line_buffer, "<migration>")) errcode = read_migration(fin);
      
      if (errcode < 0) return errcode;
    }
  }
  
  return 0;
}

int UPDData::read_basepop(std::ifstream &fin) {
  const std::string end_tag("</basepop>");
  int t, s, a;
  std::string line_buffer;
  std::vector<std::string> tokens;
  
  // read the column names
  std::getline(fin, line_buffer);
  while (std::getline(fin, line_buffer) && !GB::contains(line_buffer, end_tag)) {
    tokens = GB::split(line_buffer, _delim);
    
    t = std::atoi(tokens[0].c_str());
    s = std::atoi(tokens[1].c_str());
    a = std::atoi(tokens[2].c_str());
    
    t = (t - UPD_YEAR_START) / 5; // map from {1970, 1975, 1980, 1985} to {0, 1, 2, 3}
    s = (s == UPD_MALE) ? DP::MALE : DP::FEMALE;
    
    _basepop[t][s][a] = std::atof(tokens[3].c_str());
  }
  
  return(GB::contains(line_buffer, end_tag) ? 0 : -1);
}

int UPDData::read_lfts(std::ifstream &fin) {
  const std::string end_tag("</lfts>");
  int t, s, a;
  std::string line_buffer;
  std::vector<std::string> tokens;
  
  std::getline(fin, line_buffer); // read the column names
  while (std::getline(fin, line_buffer) && !GB::contains(line_buffer, end_tag)) {
    tokens = GB::split(line_buffer, _delim);
    
    t = std::atoi(tokens[0].c_str());
    s = std::atoi(tokens[1].c_str());
    a = std::atoi(tokens[2].c_str());
    
    t = t - UPD_YEAR_START;
    s = (s == UPD_MALE) ? DP::MALE : DP::FEMALE;
    
    _lx[t][s][a] = std::atof(tokens[3].c_str());
    _ex[t][s][a] = std::atof(tokens[4].c_str());
    _Sx[t][s][a] = std::atof(tokens[5].c_str());
  }
  
  return(GB::contains(line_buffer, end_tag) ? 0 : -1);
}

int UPDData::read_series(std::ifstream &fin, time_series_t &series, const std::string &end_tag) {
  int t;
  std::string line_buffer;
  std::vector<std::string> tokens;
  
  std::getline(fin, line_buffer); // read the column names
  while (std::getline(fin, line_buffer) && !GB::contains(line_buffer, end_tag)) {
    tokens = GB::split(line_buffer, _delim);
    
    t = std::atoi(tokens[0].c_str());
    t = t - UPD_YEAR_START;
    
    series[t] = std::atof(tokens[1].c_str());
  }
  
  return(GB::contains(line_buffer, end_tag) ? 0 : -1);
}

int UPDData::read_pasfrs(std::ifstream &fin) {
  const std::string end_tag("</pasfrs>");
  int t, a;
  std::string line_buffer;
  std::vector<std::string> tokens;
  
  std::getline(fin, line_buffer); // read the column names
  while (std::getline(fin, line_buffer) && !GB::contains(line_buffer, end_tag)) {
    tokens = GB::split(line_buffer, _delim);
    
    t = std::atoi(tokens[0].c_str());
    a = std::atoi(tokens[1].c_str());
    
    t = t - UPD_YEAR_START;
    a = a - DP::AGE_BIRTH_MIN;
    
    _pasfrs[t][a] = std::atof(tokens[2].c_str());
  }
  
  return(GB::contains(line_buffer, end_tag) ? 0 : -1);
}

int UPDData::read_migration(std::ifstream &fin) {
  const std::string end_tag("</migration>");
  int t, s, a;
  std::string line_buffer;
  std::vector<std::string> tokens;
  
  std::getline(fin, line_buffer); // read the column names
  while (std::getline(fin, line_buffer) && !GB::contains(line_buffer, end_tag)) {
    tokens = GB::split(line_buffer, _delim);
    
    t = std::atoi(tokens[0].c_str());
    s = std::atoi(tokens[1].c_str());
    a = std::atoi(tokens[2].c_str());
    
    t = t - UPD_YEAR_START;
    s = (s == UPD_MALE) ? DP::MALE : DP::FEMALE;
    
    _migration[t][s][a] = std::atof(tokens[3].c_str());
  }
  
  return(GB::contains(line_buffer, end_tag) ? 0 : -1);
}

} // end namespace DP

#endif // DPUPDDATA_IMPL_H