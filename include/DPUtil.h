#ifndef DPUTIL_H
#define DPUTIL_H

#include <boost/math/distributions/lognormal.hpp>
#include <DPData.h>

namespace DP {

	// +=+ interface +=+

	// +-+ ModelData input helpers +-+
	
	// The functions below are used to initialize ModelData inputs using
	// transformed input parameters. These are most helpful when the values that
	// make sense to the user may not be the most computationally-efficient way
	// to parameterize the model. For example, the user may have survey data on
	// the age at sexual debut, while the model calculations are designed to use
	// annual rates of sexual debut 

	// Initialize debut rates from the median age at sexual debut
	// @param dat ModelData instance to initialize
	// @param sex Sex of the population
	// @param age Median age at sexual debut
	template<typename popsize_t>
	void set_median_age_debut(ModelData<popsize_t>& dat, const int sex, const double age);

	// Initialize marriage/cohabitation rates from the median age at first union
	// @param dat ModelData instance to initialize
	// @param sex Sex of the population
	// @param age Median age at first marriage/cohabitation
	template<typename popsize_t>
	void set_median_age_union(ModelData<popsize_t>& dat, const int sex, const double age);

	// Initialize union dissolution rates from the mean duration of unions
	// @param dat ModelData instance to initialize
	// @param duration Union duration in years
	template<typename popsize_t>
	void set_mean_union_duration(ModelData<popsize_t>& dat, const double duration);

	// Initialize key population turnover rates from mean durations
	// @param dat ModelData instance to initialize
	// @param sex Sex of the population
	// @param pop Key population
	// @param duration Mean time in population in years
	template<typename popsize_t>
	void set_mean_keypop_duration(const int sex, const pop_t pop, const double duration);

	// Initialize the input age distribution of a key population
	// @param dat Model data instance to initialize
	// @param sex The key population sex
	// @param pop The key population
	// @param loc Lognormal distribution location parameter
	// @param shp Lognormal distribution shape parameter
	template<typename popsize_t>
	void set_keypop_age(ModelData<popsize_t>& dat, const int sex, const pop_t pop, const double loc, const double shp);

	// Initialize transmission probabilities per sex act
	// @param f2m female-to-male transmission probability, as % (0 to 100)
	// @param m2f male-to-female odds ratio relative to f2m
	// @param m2m male-to-male odds ratio relative to f2m
	// @param primary odds ratio of transmission during untreated primary infection
	// @param chronic odds ratio of transmission during untreated chronic infection
	// @param symptom odds ratio of transmission during untreated symptomatic infection
	// @param art_supp odds ratio of transmission on ART with suppressed viral load
	// @param art_fail odds ratio of transmission on ART with unsuppressed viral load
	template<typename popsize_t>
	void set_transmission(ModelData<popsize_t>& dat,
						  const double pct_f2m,
						  const double or_m2f,
						  const double or_m2m,
						  const double primary,
						  const double chronic,
						  const double symptom,
						  const double art_supp,
						  const double art_fail);

	// Initialize adult disease progression and HIV-related mortality rates off ART
	// 
	// Inputs are 2-dimensional arrays where rows correspond to CD4 categories
	// HIV_PRIMARY..HIV_000_050 and columns correspond to sex-age pairs. The
	// first four columns are for males aged 15-24, 25-34, 35-44, 45+, the next
	// four are for females of those ages.
	//
	// @param dat Model data instance to initialize
	// @param dist CD4 distribution after primary infection - no row for primary infection
	// @param prog HIV progression rates with untreated HIV - no row for CD4<50
	// @param mort HIV-related mortality rates with untreated HIV
	template<typename popsize_t>
	void set_adult_prog_from_10yr(ModelData<popsize_t>& dat, cd4_sex_age_ref_t& dist, cd4_sex_age_ref_t& prog, cd4_sex_age_ref_t& mort);

	// Initialize adult HIV-related mortality rates on ART
	// 
	// Inputs art1, art2, and art3 are 2-dimensional arrays where rows
	// correspond to CD4 categories HIV_PRIMARY..HIV_000_050 and columns
	// correspond to sex-age pairs. The first four columns are for males aged
	// 15-24, 25-34, 35-44, 45+, the next four are for females of those ages.
	//
	// @param dat Model data instance to initialize
	// @param art1 HIV-related mortality rates when on ART for [0,6) months
	// @param art2 HIV-related mortality rates when on ART for [6,12) months
	// @param art3 HIV-related mortality rates when on ART for [12,\infty) months
	// @param mrr  HIV-related mortality rate ratios over time on ART for [0,12) or [12,\infty) months
	template<typename popsize_t>
	void set_adult_art_mort_from_10yr(ModelData<popsize_t>& dat, cd4_sex_age_ref_t& art1, cd4_sex_age_ref_t& art2, cd4_sex_age_ref_t& art3, year_dtx_ref_t& mrr);

	/// Initialize adult ART eligibility by CD4 count thresholds
	/// @param dat Model data instance to initialize
	/// @param cd4 The highest CD4 cell count eligible for ART by year
	template<typename popsize_t>
	void set_adult_art_eligibility_from_cd4(ModelData<popsize_t>& dat, time_series_int_ref_t& cd4);

	// Initialize numbers of CLHIV aging in using Spectrum outputs
	//
	// @param dat Model data instance to initialize
	// @param clhiv 2d array with one row per year and 84 columns. Columns
	// correspond to (sex, CD4, state) tuples. There are seven states:
	// 0 Perinatal infection
	// 1 HIV acquired during breastfeeding at [0,6) months after delivery
	// 2 HIV acquired during breastfeeding at [7,12) months after delivery
	// 3 HIV acquired during breastfeeding at [12,\infty) months after delivery
	// 4 On ART [0,6) months
	// 5 On ART [7,12) months
	// 6 On ART [12,\infty) months
	// Transmission timing information is lost in the model at age 15. We assume
	// all CLHIV who are not on ART at age 15 are previously treated.
	template<typename popsize_t>
	void set_clhiv_agein(ModelData<popsize_t>& dat, boost::multi_array_ref<double, 2>& clhiv);

	// +=+ implementation +=+

	template<typename popsize_t>
	void set_median_age_debut(ModelData<popsize_t>& dat, const int sex, const double age) {
		// rate = log(2) / age; prop = 1.0 - exp(-rate); thus prop = 1.0 - 2^(-1/age)
		dat.debut_prop(sex, 1.0 - std::exp2(-1.0 / (age - DP::AGE_ADULT_MIN)));
	}

	template<typename popsize_t>
	void set_median_age_union(ModelData<popsize_t>& dat, const int sex, const double age) {
		dat.union_prop(sex, 1.0 - std::exp2(-1.0 / (age - DP::AGE_ADULT_MIN)));
	}

	template<typename popsize_t>
	void set_mean_union_duration(ModelData<popsize_t>& dat, const double duration) {
		dat.split_prop(1.0 - std::exp(-1.0 / duration));
	}

	template<typename popsize_t>
	void set_mean_keypop_duration(ModelData<popsize_t>& dat, const int sex, const pop_t pop, const double duration) {
		dat.keypop_exit_prop(sex, pop, 1.0 - std::exp(-1.0 / duration));
	}

	template<typename popsize_t>
	void set_keypop_age(ModelData<popsize_t>& dat, const int sex, const pop_t pop, const double loc, const double shp) {
		boost::math::lognormal dist(loc, shp);
		double denom(cdf(dist, DP::N_AGE_ADULT - 1)); // exclude 80+
		for (int age(0); age < DP::N_AGE_ADULT - 1; ++age)
			dat.keypop_age_dist(sex, age, pop, (cdf(dist, age+1) - cdf(dist, age)) / denom);
		dat.keypop_age_dist(sex, DP::AGE_ADULT_MAX, 0.0);
	}

	template<typename popsize_t>
	void set_transmission(ModelData<popsize_t>& dat,
						  const double transmit_f2m,
						  const double or_m2f,
						  const double or_m2m,
						  const double primary,
						  const double chronic,
						  const double symptom,
						  const double or_art_supp,
						  const double or_art_fail,
							const double or_sti_hiv_pos,
							const double or_sti_hiv_neg) {
		double base, mult, prob;
		double ratio_sex[DP::N_SEX][DP::N_SEX], ratio_hiv[DP::N_HIV_ADULT], ratio_vl[DP::N_VL];

		base = transmit_f2m; // female-to-male transmission probability

		ratio_sex[DP::FEMALE][DP::FEMALE] = 0.0;    // female-to-female transmission is not modeled
		ratio_sex[DP::FEMALE][DP::MALE  ] = 1.0;    // HIV+ female to HIV- male is the baseline group for odds ratios by sex
		ratio_sex[DP::MALE  ][DP::FEMALE] = or_m2f; // HIV+ male to HIV- female
		ratio_sex[DP::MALE  ][DP::MALE  ] = or_m2m; // male-to-male

		ratio_hiv[DP::STAGE_PRIMARY] = primary;
		ratio_hiv[DP::STAGE_CHRONIC] = chronic;
		ratio_hiv[DP::STAGE_SYMPTOM] = symptom;

		ratio_vl[DP::VL_OFF_ART] = 1.0;
		ratio_vl[DP::VL_SUCCESS] = or_art_supp;
		ratio_vl[DP::VL_FAILURE] = or_art_fail;

		for (int s_neg(0); s_neg < DP::N_SEX; ++s_neg)
			for (int s_pos(0); s_pos < DP::N_SEX; ++s_pos)
				for (int h(0); h < DP::N_STAGE; ++h)
					for (int v(0); v < DP::N_VL; ++v) {
						mult = ratio_sex[s_pos][s_neg] * ratio_hiv[h] * ratio_vl[v];
						prob = base * mult / (1.0 - base + base * mult);
						dat.hiv_risk_per_act(s_neg, s_pos, h, v, prob);
					}

		dat.effect_sti_hivpos(or_sti_hiv_pos);
		dat.effect_sti_hivneg(or_sti_hiv_neg);
	}

	template<typename popsize_t>
	void set_adult_prog_from_10yr(ModelData<popsize_t>& dat, cd4_sex_age_ref_t& dist, cd4_sex_age_ref_t& prog, cd4_sex_age_ref_t& mort) {
		const int n_age_group = 4;
		int row, col;

#ifndef SPECTRUM_CD4
		for (int h(DP::HIV_GEQ_500); h <= DP::HIV_000_050; ++h) {
			row = h - DP::HIV_GEQ_500;
			for (int a(DP::AGE_ADULT_MIN); a <= DP::AGE_ADULT_MAX; ++a) {
				col = std::min((a - DP::AGE_ADULT_MIN) / 10, n_age_group - 1);
				dat.hiv_dist(DP::MALE,   a, h, dist[row][col]);
				dat.hiv_dist(DP::FEMALE, a, h, dist[row][col + n_age_group]);
			}
		}
#else
		// When Spectrum CD4 categories are used, Spectrum CD4 categories are remapped
		// onto Goals categories. Since Goals only has inputs for 6 of 7 categories,
		// we split the input for Goals 200-350 across compartments to represent
		// Spectrum 200-250 and 250-350.
		for (int a(DP::AGE_ADULT_MIN); a <= DP::AGE_ADULT_MAX; ++a) {
			col = std::min((a - DP::AGE_ADULT_MIN) / 10, n_age_group - 1);
			dat.hiv_dist(DP::MALE, a, DP::HIV_PRIMARY, dist[0][col]);	       // Spectrum CD4>500 mapped to Goals Primary
			dat.hiv_dist(DP::MALE, a, DP::HIV_GEQ_500, dist[1][col]);        // Spectrum 350-500 mapped to Goals CD4>500
			dat.hiv_dist(DP::MALE, a, DP::HIV_350_500, dist[2][col] * 0.74); // Spectrum 250-350 mapped to Goals 350-500
			dat.hiv_dist(DP::MALE, a, DP::HIV_200_350, dist[2][col] * 0.26); // Spectrum 200-250 mapped to Goals 200-350
			dat.hiv_dist(DP::MALE, a, DP::HIV_100_200, dist[3][col]);        // Spectrum 100-200 mapped to Goals 100-200
			dat.hiv_dist(DP::MALE, a, DP::HIV_050_100, dist[4][col]);        // Spectrum  50-100 mapped to Goals  50-100
			dat.hiv_dist(DP::MALE, a, DP::HIV_000_050, dist[5][col]);        // Spectrum   0-50  mapped to Goals   0-50

			col += n_age_group;
			dat.hiv_dist(DP::FEMALE, a, DP::HIV_PRIMARY, dist[0][col]);
			dat.hiv_dist(DP::FEMALE, a, DP::HIV_GEQ_500, dist[1][col]);
			dat.hiv_dist(DP::FEMALE, a, DP::HIV_350_500, dist[2][col] * 0.74);
			dat.hiv_dist(DP::FEMALE, a, DP::HIV_200_350, dist[2][col] * 0.26);
			dat.hiv_dist(DP::FEMALE, a, DP::HIV_100_200, dist[3][col]);
			dat.hiv_dist(DP::FEMALE, a, DP::HIV_050_100, dist[4][col]);
			dat.hiv_dist(DP::FEMALE, a, DP::HIV_000_050, dist[5][col]);
		}
#endif

		for (int h(DP::HIV_PRIMARY); h <= DP::HIV_050_100; ++h) {
			for (int a(DP::AGE_ADULT_MIN); a <= DP::AGE_ADULT_MAX; ++a) {
				col = std::min((a - DP::AGE_ADULT_MIN) / 10, n_age_group - 1);
				dat.hiv_prog(DP::MALE,   a, h, prog[h][col]);
				dat.hiv_prog(DP::FEMALE, a, h, prog[h][col + n_age_group]);
			}
		}

		for (int h(DP::HIV_MIN); h <= DP::HIV_MAX; ++h) {
			for (int a(DP::AGE_ADULT_MIN); a <= DP::AGE_ADULT_MAX; ++a) {
				col = std::min((a - DP::AGE_ADULT_MIN) / 10, n_age_group - 1);
				dat.hiv_mort(DP::MALE,   a, h, mort[h][col]);
				dat.hiv_mort(DP::FEMALE, a, h, mort[h][col + n_age_group]);
			}
		}

	}

	template<typename popsize_t>
	void set_adult_art_mort_from_10yr(ModelData<popsize_t>& dat, cd4_sex_age_ref_t& art1, cd4_sex_age_ref_t& art2, cd4_sex_age_ref_t& art3, year_dtx_ref_t& mrr) {
		const int n_age_group = 4;
		int col;
		for (int t(0); t < dat.num_years(); ++t) {
			for (int a(0); a < DP::N_AGE_ADULT; ++a) {
				col = std::min(a / 10, n_age_group - 1);
				for (int h(DP::HIV_MIN); h <= DP::HIV_MAX; ++h) {
					dat.art_mort_adult(t, DP::MALE, a, h, DP::DTX_ART1, mrr[t][0] * art1[h][col]);
					dat.art_mort_adult(t, DP::MALE, a, h, DP::DTX_ART2, mrr[t][0] * art2[h][col]);
					dat.art_mort_adult(t, DP::MALE, a, h, DP::DTX_ART3, mrr[t][1] * art3[h][col]);

					dat.art_mort_adult(t, DP::FEMALE, a, h, DP::DTX_ART1, mrr[t][0] * art1[h][col + n_age_group]);
					dat.art_mort_adult(t, DP::FEMALE, a, h, DP::DTX_ART2, mrr[t][0] * art2[h][col + n_age_group]);
					dat.art_mort_adult(t, DP::FEMALE, a, h, DP::DTX_ART3, mrr[t][1] * art3[h][col + n_age_group]);
				}
			}
		}
	}

	template<typename popsize_t>
	void set_adult_art_eligibility_from_cd4(ModelData<popsize_t>& dat, time_series_int_ref_t& cd4) {
		int h;
		for (int t(0); t < dat.num_years(); ++t) {
			h = DP::HIV_ADULT_MIN;
			while ((h < DP::HIV_ADULT_MAX) && (DP::CD4_ADULT_LOWER[h] >= cd4[t]))
				++h;
			dat.art_first_eligible_stage_adult(t,h);
		}
	}

	template<typename popsize_t>
	void set_clhiv_agein(ModelData<popsize_t>& dat, boost::multi_array_ref<double, 2>& clhiv) {
		const int n_state(7), n_acq(4);
		double buffer[n_state], off_art;
		int row;

		// s_model = sex_map[s_input] transforms from input to Goals ARM sex ordering.
		int s_input, s_model;
		const int sex_map[] = { DP::MALE, DP::FEMALE };

		for (int t(0); t < dat.num_years(); ++t) {
			for (s_input = 0; s_input < 2; ++s_input) {
				s_model = sex_map[s_input];
				for (int h(0); h < DP::N_HIV_CHILD_PED; ++h) {
					// Cache inputs by HIV acquisition time or ART state
					for (int d(0); d < n_state; ++d) {
						row = s_input * (DP::N_HIV_CHILD_PED * n_state) + h * n_state + d;
						buffer[d] = clhiv[t][row];
					}

					// Sum children off ART by acquisition mode then assign them to previous ART
					off_art = 0.0;
					for (int d(0); d < n_acq; ++d)
						off_art += buffer[d];
					dat.clhiv_agein(t, s_model, h, DP::DTX_PREV_TX, off_art);
					dat.clhiv_agein(t, s_model, h, DP::DTX_UNAWARE, 0.0);
					dat.clhiv_agein(t, s_model, h, DP::DTX_AWARE,   0.0);

					// Transfer children on ART directly. We need to map input row indices (4..6) to
					// ART indices (DTX_ART1..DTX_ART3) to do this correctly
					for (int d(DP::DTX_ART1); d <= DP::DTX_ART3; ++d)
						dat.clhiv_agein(t, s_model, h, d, buffer[d - DP::DTX_ART1 + n_acq]);
				}
			}
		}
	}

} // END namespace DP

#endif // DPUTIL_H
