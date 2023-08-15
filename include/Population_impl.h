#ifndef POPULATION_IMPL_H
#define POPULATION_IMPL_H

namespace DP {

    Population::Population(const int year_min, const int year_max)
        : _adult_neg(NULL), _adult_hiv(NULL), _child_neg(NULL), _child_hiv(NULL) {
      _year_first = year_min;
      _year_final = year_max;
      _n_year = year_max - year_min + 1;
    }

    Population::~Population() {
      if (_adult_neg) { delete _adult_neg; }
      if (_adult_hiv) { delete _adult_hiv; }
      if (_child_neg) { delete _child_neg; }
      if (_child_hiv) { delete _child_hiv; }
    }

    void Population::share_storage(
        double* ptr_adult_neg,
        double* ptr_adult_hiv,
        double* ptr_child_neg,
        double* ptr_child_hiv) {
        const int n_year(year_final() - year_first() + 1);
        _adult_neg = new adult_neg_t(ptr_adult_neg, boost::extents[n_year][N_SEX_MC][N_AGE_ADULT][N_POP]);
        _adult_hiv = new adult_hiv_t(ptr_adult_hiv, boost::extents[n_year][N_SEX_MC][N_AGE_ADULT][N_POP][N_HIV_ADULT][N_DTX]);
        _child_neg = new child_neg_t(ptr_child_neg, boost::extents[n_year][N_SEX_MC][N_AGE_CHILD]);
        _child_hiv = new child_hiv_t(ptr_child_hiv, boost::extents[n_year][N_SEX_MC][N_AGE_CHILD][N_HIV_CHILD][N_DTX]);
    }

    int Population::year_first() const {
      return _year_first;
    }

    int Population::year_final() const {
      return _year_final;
    }

    int Population::num_years() const {
      return _n_year;
    }

    void Population::initialize(double value) {
      std::fill_n(_child_neg->data(), _child_neg->num_elements(), value);
      std::fill_n(_child_hiv->data(), _child_hiv->num_elements(), value);
      std::fill_n(_adult_neg->data(), _adult_neg->num_elements(), value);
      std::fill_n(_adult_hiv->data(), _adult_hiv->num_elements(), value);
    }

} // END namespace DP

#endif // POPULATION_IMPL_H
