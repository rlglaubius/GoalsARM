#ifndef DPDEFS_H
#define DPDEFS_H

#include <boost/multi_array.hpp>

namespace DP {
	// Many typedefs specified here map to the same underlying boost::multi_array. We
	// maintain distinct typedefs (e.g., year_age_t, year_sex_t) for readability when
	// instantiating these array types. To avoid a huge number of typedefs, our convention
	// is to define new types using the ordering
	//
	// 1. time (usually year)
	// 2. sex
	// 3. age
	// 4. behavioral risk group
	// 5. HIV disease stage
	// 6. HIV diagnosis and treatment state
	//
	// For example, if you plan to define a type to store indicators
	// by time, risk, and HIV disease stage, you would define this as
	// typedef boost::multi_array<double, 3> year_risk_hiv_t;
	// To optimize performance, define nested loops in this same ordering,
	// e.g., for (t in years) {for (r in risk) {for (h in hiv_stages) {}}}.
	// Doing so optimizes performance by minimizing cache misses.

	typedef boost::multi_array<double, 3> year_sex_age_t;
	typedef boost::multi_array<double, 3> year_sex_pop_t;
	typedef boost::multi_array<double, 4> year_sex_age_pop_t;
	typedef boost::multi_array<double, 2> year_age_t;
	typedef boost::multi_array<double, 2> year_sex_t;
	typedef boost::multi_array<double, 4> year_sex_hiv_dtx_t;
	typedef boost::multi_array<double, 4> year_sex_age_pop_t;
	typedef boost::multi_array<double, 5> year_sex_age_hiv_dtx_t;
	typedef boost::multi_array<double, 2> year_bond_t;
	typedef std::vector<double> time_series_t;
	typedef std::vector<int>    time_series_int_t;

	typedef boost::multi_array<double, 2> sex_age_t;
	typedef boost::multi_array<double, 2> sex_hiv_t;
	typedef boost::multi_array<double, 3> sex_age_hiv_t;
	typedef boost::multi_array<double, 2> sex_pop_t;

	typedef boost::multi_array<double, 5> mixing_matrix_t;

	// boost::multi_array instances manage their own memory. multi_array_ref
	// instances provide a multi_array interface, but use memory allocated
	// elsewhere. We can use these *_ref_t types to use data in preallocated
	// memory without copying those data
	typedef boost::multi_array_ref<double, 2> sex_pop_ref_t;
	typedef boost::multi_array_ref<double, 2> year_age_ref_t;
	typedef boost::multi_array_ref<double, 2> year_sex_ref_t;
	typedef boost::multi_array_ref<double, 2> year_dtx_ref_t;
	typedef boost::multi_array_ref<double, 4> year_sex_age_pop_ref_t;

	typedef boost::multi_array_ref<double, 1> time_series_ref_t;
	typedef boost::multi_array_ref<int, 1> time_series_int_ref_t;

	// alias for a two-dimensional array with rows corresponding
	// to CD4 categories and columns for (sex,age) pairs
	typedef boost::multi_array_ref<double, 2> cd4_sex_age_ref_t;

	// More general typedefs
	typedef boost::multi_array_ref<double, 4> array4d_ref_t;

} // END namespace DP

#endif // DPDEFS_H