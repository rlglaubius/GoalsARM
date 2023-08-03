#include <algorithm>
#include <iostream>
#include <DPProjection.H>

enum result_t {TEST_SUCCESS = 0, TEST_FAILURE = 1};

result_t test_births() {
	const int year_first(1970), year_final(1971), num_years(2);
	const double tfr(7.33), srb(101.4);
	const double pasfrs[DP::N_AGE_BIRTH] = {
		0.0271, 0.0271, 0.0271, 0.0271, 0.0271,
		0.0393, 0.0393, 0.0393, 0.0393, 0.0393,
		0.0404, 0.0404, 0.0404, 0.0404, 0.0404,
		0.0356, 0.0356, 0.0356, 0.0356, 0.0356,
		0.0283, 0.0283, 0.0283, 0.0283, 0.0283,
		0.0177, 0.0177, 0.0177, 0.0177, 0.0177,
		0.0116, 0.0116, 0.0116, 0.0116, 0.0116};

	const double pop_t0[DP::N_AGE_BIRTH] = {
		46937, 45325, 43840, 42144, 40549,
		42038, 42973, 42724, 43281, 42085,
		40528, 38887, 37321, 35867, 34529,
		33343, 32357, 31368, 30192, 28948,
		27840, 26854, 25888, 24968, 24092,
		23256, 22488, 21721, 20890, 20059,
		19272, 18473, 17824, 17438, 17226};

	const double pop_t1[DP::N_AGE_BIRTH] = {
		48652, 46632, 45019, 43533, 41835,
		40239, 41700, 42613, 42355, 42898,
		41703, 40154, 38526, 36970, 35528,
		34200, 33022, 32043, 31060, 29891,
		28657, 27558, 26578, 25617, 24702,
		23829, 23001, 22238, 21475, 20654,
		19829, 19047, 18254, 17607, 17216};

	int a, b;
	double births;

	boost::multi_array<double, 3> child_neg(boost::extents[num_years][DP::N_SEX_MC][DP::N_AGE_CHILD]);
	boost::multi_array<double, 5> child_hiv(boost::extents[num_years][DP::N_SEX_MC][DP::N_AGE_CHILD][DP::N_HIV][DP::N_DTX]);
	boost::multi_array<double, 4> adult_neg(boost::extents[num_years][DP::N_SEX_MC][DP::N_AGE_ADULT][DP::N_POP]);
	boost::multi_array<double, 6> adult_hiv(boost::extents[num_years][DP::N_SEX_MC][DP::N_AGE_ADULT][DP::N_POP][DP::N_HIV][DP::N_DTX]);

	std::fill_n(child_neg.data(), child_neg.num_elements(), 0.0);
	std::fill_n(child_hiv.data(), child_hiv.num_elements(), 0.0);
	std::fill_n(adult_neg.data(), adult_neg.num_elements(), 0.0);
	std::fill_n(adult_hiv.data(), adult_hiv.num_elements(), 0.0);

	DP::Projection proj(year_first, year_final);
	proj.pop.share_storage(adult_neg.data(), adult_hiv.data(), child_neg.data(), child_hiv.data());

	proj.dat.tfr(year_final - year_first, tfr);
	proj.dat.srb(year_final - year_first, srb);
	for (b = 0; b < DP::N_AGE_BIRTH; ++b) {
		a = b + DP::AGE_BIRTH_MIN;
		proj.dat.pasfrs(year_final - year_first, a, pasfrs[b]);
		proj.pop.adult_neg(year_first - year_first, DP::FEMALE, b, DP::POP_NOSEX) = pop_t0[b];
		proj.pop.adult_neg(year_final - year_first, DP::FEMALE, b, DP::POP_NOSEX) = pop_t1[b];
	}

	births = proj.calc_births(year_final - year_first);

	return fabs(births - 251855) < 0.5 ? TEST_SUCCESS : TEST_FAILURE;
}

result_t test_births_exposed() {
	const int year_first(1970), year_final(1971);
	DP::Projection proj(year_first, year_final);


	// set up population
	// set up PASFR
	// set up TFR
	// set up FRR
	// calc
	// check against known good value
	// return 0 on match, -1 on failure
	return TEST_SUCCESS;
}

struct TestRecord {
	std::string name;
	result_t (*func)(void);
};

static TestRecord TestRegistry[] = {
	{"test_births", test_births},
	{"", NULL}
};

int main(int argc, char **argv) {
	result_t result;

	for (int k(0); TestRegistry[k].name != ""; ++k) {
		result = TestRegistry[k].func();
		std::cout << TestRegistry[k].name << '\t' << (result == TEST_SUCCESS ? "Pass" : "Fail") << std::endl;
	}

	return 0;
}
