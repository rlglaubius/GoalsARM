#ifndef POPULATION_H
#define POPULATION_H

#include <vector>
#include <boost/multi_array.hpp>
#include <DPConst.h>

namespace DP {

// Class for population structures
class Population {
public:
  // +-+ Nested types +-+
  typedef boost::multi_array_ref<double, 4> adult_neg_t; // HIV-negative adults, stratified by year, sex, age, risk
  typedef boost::multi_array_ref<double, 6> adult_hiv_t; // HIV-positive adults, stratified by year, sex, age, risk, CD4, and care status
  typedef boost::multi_array_ref<double, 3> child_neg_t; // HIV-negative children, stratified by year, sex and age
  typedef boost::multi_array_ref<double, 5> child_hiv_t; // HIV-positive children, stratified by year, sex, age, CD4, and care status

  // +-+ Methods +-+
  // Constructors
  Population(const int year_min, const int year_max);
  ~Population();

  /// Share memory for storing population sizes
  /// @param ptr_adult_neg HIV-negative adults by year, sex, age, and risk
  /// @param ptr_adult_hiv HIV-positive adults by year, sex, age, risk, CD4, and care status
  /// @param ptr_child_neg HIV-negative children by year, sex, age
  /// @param ptr_child_hiv HIV-positive children by year, sex, age, CD4, and care status
  void share_storage(
      double* ptr_adult_neg,
      double* ptr_adult_hiv,
      double* ptr_child_neg,
      double* ptr_child_hiv);
  
  // Accessors
  int year_first() const;
  int year_final() const;
  int num_years() const;

  // Population accessors: "get" methods
  inline double adult_neg(int t, int s, int a, int r) const { return (*_adult_neg)[t][s][a][r]; }
  inline double adult_hiv(int t, int s, int a, int r, int h, int d) const { return (*_adult_hiv)[t][s][a][r][h][d]; }
  inline double child_neg(int t, int s, int a) const { return (*_child_neg)[t][s][a]; }
  inline double child_hiv(int t, int s, int a, int h, int d) const { return (*_child_hiv)[t][s][a][h][d]; }

  // Population accessors: "set" methods
  inline double& adult_neg(int t, int s, int a, int r) { return (*_adult_neg)[t][s][a][r]; }
  inline double& adult_hiv(int t, int s, int a, int r, int h, int d) { return (*_adult_hiv)[t][s][a][r][h][d]; }
  inline double& child_neg(int t, int s, int a) { return (*_child_neg)[t][s][a]; }
  inline double& child_hiv(int t, int s, int a, int h, int d) { return (*_child_hiv)[t][s][a][h][d]; }
  
  // Convenience functions
  void initialize(double value); // Set all compartment sizes = value

private:
  int _year_first;
  int _year_final;
  int _n_year;

  // Population state variables
  adult_neg_t* _adult_neg; // HIV-negative adults
  adult_hiv_t* _adult_hiv; // HIV-positive adults
  child_neg_t* _child_neg; // HIV-negative children
  child_hiv_t* _child_hiv; // HIV-positive children
};

} // END namespace DP

#endif // POPULATION_H