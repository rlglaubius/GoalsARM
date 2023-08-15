#include <algorithm>
#include <iostream>
#include <GoalsARM.H>

enum result_t {TEST_SUCCESS = 0, TEST_FAILURE = 1};

void setup_projection(DP::Projection& proj) {
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

	const double prev[DP::N_AGE_BIRTH] = {
		0.037, 0.037, 0.037, 0.037, 0.037,
		0.132, 0.132, 0.132, 0.132, 0.132,
		0.155, 0.155, 0.155, 0.155, 0.155,
		0.181, 0.181, 0.181, 0.181, 0.181,
		0.170, 0.170, 0.170, 0.170, 0.170,
		0.179, 0.179, 0.179, 0.179, 0.179,
		0.133, 0.133, 0.133, 0.133, 0.133};

	const double frr_age[DP::N_AGE_BIRTH] = {
		1.24, 1.24, 1.24, 1.24, 1.24,
		0.86, 0.86, 0.86, 0.86, 0.86,
		0.83, 0.83, 0.83, 0.83, 0.83,
		0.82, 0.82, 0.82, 0.82, 0.82,
		0.60, 0.60, 0.60, 0.60, 0.60,
		0.60, 0.60, 0.60, 0.60, 0.60,
		0.60, 0.60, 0.60, 0.60, 0.60};

	const double frr_art[DP::N_AGE_BIRTH] = {
		1.09, 1.09, 1.09, 1.09, 1.09,
		1.01, 1.01, 1.01, 1.01, 1.01,
		0.92, 0.92, 0.92, 0.92, 0.92,
		0.81, 0.81, 0.81, 0.81, 0.81,
		0.64, 0.64, 0.64, 0.64, 0.64,
		0.64, 0.64, 0.64, 0.64, 0.64,
		0.64, 0.64, 0.64, 0.64, 0.64};

	const double frr_cd4[DP::N_HIV_ADULT] = {1.00, 1.00, 0.96, 0.85, 0.61, 0.38, 0.30};
	const double frr_loc(0.96);
	const double cd4[DP::N_HIV_ADULT] = {0.0000, 0.5735, 0.2186, 0.1592, 0.0420, 0.0058, 0.0009};
	const double dtx[DP::N_DTX] = {0.40, 0.10, 0.10, 0.05, 0.05, 0.30};
	const int yidx_first(0), yidx_final(1);
	int a, b, h, d;

	proj.dat.tfr(yidx_final, tfr);
	proj.dat.srb(yidx_final, srb);
	for (b = 0; b < DP::N_AGE_BIRTH; ++b) {
		a = b + DP::AGE_BIRTH_MIN;
		proj.dat.pasfrs(yidx_final, a, pasfrs[b]);

		proj.dat.frr_age_no_art(yidx_final, b, frr_loc * frr_age[b]);
		proj.dat.frr_age_on_art(b, frr_loc * frr_art[b]);

		proj.pop.adult_neg(yidx_first, DP::FEMALE, b, DP::POP_NOSEX) = pop_t0[b] * (1.0 - prev[b]);
		proj.pop.adult_neg(yidx_final, DP::FEMALE, b, DP::POP_NOSEX) = pop_t1[b] * (1.0 - prev[b]);
		for (h = DP::HIV_ADULT_MIN; h <= DP::HIV_ADULT_MAX; ++h) {
			for (d = DP::DTX_MIN; d <= DP::DTX_MAX; ++d) {
				proj.pop.adult_hiv(yidx_first, DP::FEMALE, b, DP::POP_NOSEX, h, d) = pop_t0[b] * prev[b] * cd4[h] * dtx[d];
				proj.pop.adult_hiv(yidx_final, DP::FEMALE, b, DP::POP_NOSEX, h, d) = pop_t1[b] * prev[b] * cd4[h] * dtx[d];
			}
		}
	}

	for (h = DP::HIV_ADULT_MIN; h <= DP::HIV_ADULT_MAX; ++h)
		proj.dat.frr_cd4_no_art(h, frr_cd4[h]);
}

result_t test_births() {
	const int year_first(1970), year_final(1971), num_years(2);
	const double target_births(251855), tolerance(0.5);
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
	setup_projection(proj);

	births = proj.calc_births(year_final - year_first);
	return fabs(births - target_births) < tolerance ? TEST_SUCCESS : TEST_FAILURE;
}

result_t test_births_hiv_exposed() {
	const int year_first(1970), year_final(1971), num_years(2);
	const double target_births(26895), tolerance(0.5);
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
	setup_projection(proj);

	births = proj.calc_births_hiv_exposed(year_final - year_first);
	return fabs(births - target_births) < tolerance ? TEST_SUCCESS : TEST_FAILURE;
}

struct TestRecord {
	std::string name;
	result_t (*func)(void);
};

static TestRecord TestRegistry[] = {
	{"test_births",             test_births},
	{"test_births_hiv_exposed", test_births_hiv_exposed},
	// Add new tests here
	{"", NULL} // Sentinel value, do not remove
};

int main(int argc, char **argv) {
	result_t result;

	for (int k(0); TestRegistry[k].name != ""; ++k) {
		result = TestRegistry[k].func();
		std::cout << (result == TEST_SUCCESS ? "Pass" : "Fail") << " : " << TestRegistry[k].name << std::endl;
	}

	return 0;
}
