#include <algorithm>
#include <iostream>
#include <numeric>
#include <DPProjection.h>

namespace DP {

	Projection::Projection(const int year_start, const int year_final) 
	: pop(year_start, year_final),
		dth(year_start, year_final),
		dat(year_start, year_final),
		_last_valid_time(-1) {
		_year_first = year_start;
		_year_final = year_final;
		_num_years = year_final - year_start + 1;
	}

	Projection::~Projection() {}

	void Projection::initialize(const std::string &upd_filename) {
		dat.initialize(upd_filename);
	}

	void Projection::project(const int year_end) {
		const int time_end(std::min(year_end - year_first(), num_years()));
		const int time_bgn(std::max(_last_valid_time, 0) + 1);

		if (_last_valid_time < 0) {
			init_baseyear_population();
			calc_births_baseyear();
			calc_deaths_baseyear();
			calc_popsize(0);
		}

		for (int t(time_bgn); t <= time_end; ++t) {
			project_one_year(t);
			calc_popsize(t);
		}

		_last_valid_time = time_end;
	}

	void Projection::calc_popsize(const int t) {
		int a, b, d, h, r, s;
		popsize_t count;

		for (a = DP::AGE_CHILD_MIN; a <= DP::AGE_CHILD_MAX; ++a) {
			s = DP::FEMALE;
			count = pop.child_neg(t,s,a);
			for (h = DP::HIV_CHILD_MIN; h <= DP::HIV_CHILD_MAX; ++h)
				for (d = DP::DTX_MIN; d <= DP::DTX_MAX; ++d)
					count += pop.child_hiv(t,s,a,h,d);
			dat.popsize(t,s,a,count);

			s = DP::MALE;
			count = pop.child_neg(t, DP::MALE_U, a) + pop.child_neg(t, DP::MALE_C, a);
			for(h = DP::HIV_CHILD_MIN; h <= DP::HIV_CHILD_MAX; ++h)
				for(d = DP::DTX_MIN; d <= DP::DTX_MAX; ++d)
					count += pop.child_hiv(t, DP::MALE_U, a, h, d) + pop.child_hiv(t, DP::MALE_C, a, h, d);
			dat.popsize(t,s,a,count);
		}

		for (a = DP::AGE_ADULT_MIN; a <= DP::AGE_ADULT_MAX; ++a) {
			b = a - DP::AGE_ADULT_MIN;

			s = DP::FEMALE;
			count = 0.0;
			for (r = DP::POP_MIN; r <= DP::POP_MAX; ++r)
				count += pop.adult_neg(t, s, b, r);
			for (r = DP::POP_MIN; r <= DP::POP_MAX; ++r)
				for (h = DP::HIV_ADULT_MIN; h <= DP::HIV_ADULT_MAX; ++h)
					for (d = DP::DTX_MIN; d <= DP::DTX_MAX; ++d)
						count += pop.adult_hiv(t, s, b, r , h, d);
			dat.popsize(t,s,a,count);

			s = DP::MALE;
			count = 0.0;
			for (r = DP::POP_MIN; r <= DP::POP_MAX; ++r)
				count += pop.adult_neg(t, DP::MALE_U, b, r) + pop.adult_neg(t, DP::MALE_C, b, r);
			for (r = DP::POP_MIN; r <= DP::POP_MAX; ++r)
				for (h = DP::HIV_ADULT_MIN; h <= DP::HIV_ADULT_MAX; ++h)
					for (d = DP::DTX_MIN; d <= DP::DTX_MAX; ++d)
						count += pop.adult_hiv(t, DP::MALE_U, b, r, h, d) + pop.adult_hiv(t, DP::MALE_C, b, r, h, d);
			dat.popsize(t,s,a,count);
		}
	}

	void Projection::init_baseyear_population() {
		const int t(0), r(DP::POP_NOSEX);
		int a, s;

		pop.initialize(0.0); // Initialize population sizes to zero
		dth.initialize(0.0); // Initialize deaths to zero

		for (s = DP::SEX_MIN; s <= DP::SEX_MAX; ++s) {
			for (a = 0; a < DP::N_AGE_CHILD; ++a) {
				pop.child_neg(t, s, a) = dat.basepop(s,a);
			}
		}

		for (s = DP::SEX_MIN; s <= DP::SEX_MAX; ++s) {
			for (a = 0; a < DP::N_AGE_ADULT; ++a) {
				pop.adult_neg(t, s, a, r) = dat.basepop(s,a + DP::AGE_ADULT_MIN);
			}
		}

		init_baseyear_risk();
		init_baseyear_male_circumcision();
	}

	/// @pre pop.adult_neg(0, s, a, DP::POP_NOSEX) stores the whole population for sex s and age a. Males
	/// are not yet disaggregated by circumcision status.
	void Projection::init_baseyear_risk() {
		const int t(0);
		int a, k, r, s;
		double size_fert, size_curr, size_prev, scale;
		double p_naive, n_nosex, n_never, n_union, n_split, p_union;
		double kp_need, kp_have, kp_pool;
		double size_active[DP::N_SEX];
		double kp_stay_enter[DP::N_SEX][DP::N_POP];
		double kp_turn_enter[DP::N_SEX][DP::N_AGE_ADULT], kp_turn_leave[DP::N_SEX][DP::N_AGE_ADULT];
		double p_stay_enter[DP::N_SEX], p_turn_enter;
		double n_total[DP::N_SEX][DP::N_AGE_ADULT];

		int nt(0), ns(0);
		int kp_turn_sex[DP::N_SEX * DP::N_POP], kp_stay_sex[DP::N_SEX * DP::N_POP];
		int kp_turn_pop[DP::N_SEX * DP::N_POP], kp_stay_pop[DP::N_SEX * DP::N_POP];

		// Enumerate key populations with and without turnover ("turn" and "stay", respectively)
		for (s = DP::SEX_MIN; s <= DP::SEX_MAX; ++s)
			for (r = DP::POP_KEY_MIN; r < N_POP_SEX[s]; ++r)
				if (dat.keypop_stay(s, r)) {
					kp_stay_sex[ns] = s;
					kp_stay_pop[ns] = r;
					++ns;
				} else {
					kp_turn_sex[nt] = s;
					kp_turn_pop[nt] = r;
					++nt;
				}

		for (s = DP::SEX_MIN; s <= DP::SEX_MAX; ++s)
			for (a = 0; a < DP::N_AGE_ADULT; ++a)
				n_total[s][a] = pop.adult_neg(t, s, a, DP::POP_NOSEX);

		size_fert = 0.0;
		for (s = DP::SEX_MIN; s <= DP::SEX_MAX; ++s)
			for (a = 0; a < DP::N_AGE_BIRTH; ++a)
				size_fert += n_total[s][a];

		for (s = DP::SEX_MIN; s <= DP::SEX_MAX; ++s) {
			size_active[s] = 0.0;
			p_naive = 1.0;
			for (a = 0; a < DP::N_AGE_BIRTH; ++a) {
				p_naive *= 1.0 - dat.debut_prop(s);
				size_active[s] += n_total[s][a] * (1.0 - p_naive);
			}
		}

		for (s = DP::SEX_MIN; s <= DP::SEX_MAX; ++s) {
			for (a = 0; a < DP::N_AGE_ADULT; ++a) {
				kp_turn_enter[s][a] = 0.0;
				kp_turn_leave[s][a] = 0.0;
			}
		}

		p_stay_enter[DP::FEMALE] = 0.0;
		p_stay_enter[DP::MALE  ] = 0.0;
		for (k = 0; k < ns; ++k) { // key populations without turnover
			s = kp_stay_sex[k];
			r = kp_stay_pop[k];
			kp_stay_enter[s][r] = dat.keypop_size(s, r) * size_fert / size_active[s];
			p_stay_enter[s] += kp_stay_enter[s][r];
		}

		for (k = 0; k < nt; ++k) { // key populations with turnover
			s = kp_turn_sex[k];
			r = kp_turn_pop[k];

			scale = 0.0;
			for (a = 0; a < DP::N_AGE_BIRTH; ++a)
				scale += dat.keypop_age_dist(s, a, r);

			scale = dat.keypop_size(s, r) * size_fert / scale;
			for (a = 0; a < DP::N_AGE_ADULT; ++a)
				pop.adult_neg(t, s, a, r) = scale * dat.keypop_age_dist(s, a, r);

			kp_turn_enter[s][0] += pop.adult_neg(t, s, 0, r);
			kp_turn_leave[s][0] += pop.adult_neg(t, s, 0, r) * dat.keypop_exit_prop(s, r);
			size_curr = n_total[s][0];
			for (a = 1; a < DP::N_AGE_ADULT; ++a) { // calculate numbers entering or exiting by age
				size_prev = size_curr;
				size_curr = n_total[s][a];
				kp_need = pop.adult_neg(t, s, a,   r);
				kp_have = pop.adult_neg(t, s, a-1, r) * (1.0 - dat.keypop_exit_prop(s, r)) * size_curr / size_prev;
				if (kp_need > kp_have)
					kp_turn_enter[s][a] += kp_need - kp_have;
				else
					kp_turn_leave[s][a] += (kp_have - kp_need) * size_prev / size_curr;
				kp_turn_leave[s][a] += pop.adult_neg(t, s, a-1, r) * dat.keypop_exit_prop(s,r);
			}
		}

		for (s = DP::SEX_MIN; s <= DP::SEX_MAX; ++s) {
			a = 0;
			size_curr = n_total[s][a];
			p_turn_enter = kp_turn_enter[s][a] / size_curr;

			pop.adult_neg(t, s, a, DP::POP_NOSEX) = size_curr * (1.0 - dat.debut_prop(s) - p_turn_enter);
			pop.adult_neg(t, s, a, DP::POP_NEVER) = size_curr * dat.debut_prop(s) * (1.0 - dat.union_prop(s) - p_stay_enter[s]);
			pop.adult_neg(t, s, a, DP::POP_UNION) = size_curr * dat.debut_prop(s) * dat.union_prop(s);
			pop.adult_neg(t, s, a, DP::POP_SPLIT) = 0.0;

			for (a = 1; a < DP::N_AGE_ADULT; ++a) {
				size_prev = size_curr;
				size_curr = n_total[s][a];

				// Cache population sizes, adjusted for relative birth cohort sizes. The
				// adjustment size_curr / size_prev approximates survival and net migration
				// to the population aged a-1 in year t-1.
				n_nosex = pop.adult_neg(t, s, a - 1, DP::POP_NOSEX) * size_curr / size_prev;
				n_never = pop.adult_neg(t, s, a - 1, DP::POP_NEVER) * size_curr / size_prev;
				n_union = pop.adult_neg(t, s, a - 1, DP::POP_UNION) * size_curr / size_prev;
				n_split = pop.adult_neg(t, s, a - 1, DP::POP_SPLIT) * size_curr / size_prev;
				p_union = n_union / (n_union + n_split);

				kp_pool = n_nosex + n_never + n_union + n_split;
				p_turn_enter = kp_turn_enter[s][a] / kp_pool;

				// GP
				pop.adult_neg(t, s, a, DP::POP_NOSEX) = n_nosex * (1.0 - dat.debut_prop(s) - p_turn_enter);
				pop.adult_neg(t, s, a, DP::POP_NEVER) = n_nosex * dat.debut_prop(s) * (1.0 - dat.union_prop(s) - p_stay_enter[s])
					                                    + n_never * (1.0 - dat.union_prop(s) - p_turn_enter);
				pop.adult_neg(t, s, a, DP::POP_UNION) = n_nosex * dat.debut_prop(s) * dat.union_prop(s)
					                                    + n_never * dat.union_prop(s)
					                                    + n_union * (1.0 - dat.split_prop() - p_turn_enter)
					                                    + n_split * dat.union_prop(s)
					                                    + (size_curr / size_prev) * kp_turn_leave[s][a] * p_union;
				pop.adult_neg(t, s, a, DP::POP_SPLIT) = n_union * dat.split_prop()
					                                    + n_split * (1.0 - dat.union_prop(s) - p_turn_enter)
					                                    + (size_curr / size_prev) * kp_turn_leave[s][a] * (1.0 - p_union);
			}
		}

		// Key populations without turnover
		for (k = 0; k < ns; ++k) {
			s = kp_stay_sex[k];
			r = kp_stay_pop[k];
			size_curr = n_total[s][0];
			pop.adult_neg(t, s, 0, r) = size_curr * dat.debut_prop(s) * kp_stay_enter[s][r];
			for (a = 1; a < DP::N_AGE_ADULT; ++a) {
				size_prev = size_curr;
				size_curr = n_total[s][a];
				pop.adult_neg(t, s, a, r) = pop.adult_neg(t, s, a - 1, r) + pop.adult_neg(t, s, a-1, DP::POP_NOSEX) * dat.debut_prop(s) * kp_stay_enter[s][r];
				pop.adult_neg(t, s, a, r) *= size_curr / size_prev; // adjust for relative birth cohort sizes
			}
		}
	}

	void Projection::init_baseyear_male_circumcision() {
		const int t(0);
		double prop[DP::N_AGE];
		double n;
		int a, b, r;

		// Calculate the cumulative proportion circumcised in each age group if past
		// male circumcision uptake was constant at first-year input levels
		prop[DP::AGE_MIN] = dat.uptake_male_circumcision(t, DP::AGE_MIN);
		for (a = DP::AGE_MIN + 1; a <= DP::AGE_MAX; ++a)
			prop[a] = prop[a-1] + (1.0 - prop[a-1]) * dat.uptake_male_circumcision(t, a);

		// Children
		for (a = DP::AGE_CHILD_MIN; a <= DP::AGE_CHILD_MAX; ++a) {
			n = prop[a] * pop.child_neg(t, DP::MALE_U, a);
			pop.child_neg(t, DP::MALE_U, a) -= n;
			pop.child_neg(t, DP::MALE_C, a) += n;
		}

		// Adults
		for (a = DP::AGE_ADULT_MIN; a <= DP::AGE_ADULT_MAX; ++a) {
			b = a - DP::AGE_ADULT_MIN;
			for (r = DP::POP_MIN; r <= DP::POP_MAX; ++r) {
				n = prop[a] * pop.adult_neg(t, DP::MALE_U, b, r);
				pop.adult_neg(t, DP::MALE_U, b, r) -= n;
				pop.adult_neg(t, DP::MALE_C, b, r) += n;
			}
		}
	}

	void Projection::calc_births_baseyear() {
		const int t(0), s(DP::FEMALE);

		const double perc_m(dat.srb(t) / (dat.srb(t) + 100.0));
		const double perc_f(1.0 - perc_m);
		popsize_t female[DP::N_AGE];
		popsize_t denom(0.0), births(0.0);
		int a, b, d, h, r;

		for (a = 14; a <= DP::AGE_BIRTH_MAX; ++a)
			female[a] = 0.0;

		// Count reproductive-age females, including females who will be of reproductive age at next birthday
		a = 14;
		female[a] += pop.child_neg(t, s, a);

		for (h = DP::HIV_CHILD_MIN; h <= DP::HIV_CHILD_MAX; ++h)
			for (d = DP::DTX_MIN; d <= DP::DTX_MAX; ++d)
				female[a] += pop.child_hiv(t, s, a, h, d);

		for (a = DP::AGE_BIRTH_MIN; a <= DP::AGE_BIRTH_MAX; ++a) {
			b = a - DP::AGE_BIRTH_MIN;
			for (r = DP::POP_MIN; r <= DP::POP_MAX; ++r)
				female[a] += pop.adult_neg(t, s, b, r);
		}

		for (a = DP::AGE_BIRTH_MIN; a <= DP::AGE_BIRTH_MAX; ++a) {
			b = a - DP::AGE_BIRTH_MIN;
			for (r = DP::POP_MIN; r <= DP::POP_MAX; ++r)
				for (h = DP::HIV_ADULT_MIN; h <= DP::HIV_ADULT_MAX; ++h)
					for (d = DP::DTX_MIN; d <= DP::DTX_MAX; ++d)
						female[a] += pop.adult_hiv(t, s, b, r, h, d);
		}

		for (a = DP::AGE_BIRTH_MIN; a <= DP::AGE_BIRTH_MAX; ++a)
			denom += dat.pasfrs(t,a);

		for (a = DP::AGE_BIRTH_MIN; a <= DP::AGE_BIRTH_MAX; ++a)
			births += 0.5 * (female[a] + female[a-1] * dat.Sx(t,s,a-1)) * dat.pasfrs(t,a) * dat.tfr(t) / denom;

		dat.births(t, DP::MALE,   births * perc_m);
		dat.births(t, DP::FEMALE, births * perc_f);
	}

	void Projection::calc_deaths_baseyear() {
		const int t(0);
		double mort, deaths;
		int a, s;

		for (s = DP::SEX_MIN; s <= DP::SEX_MAX; ++s) {
			// age 0
			a = DP::AGE_MIN;
			mort = 1.0 - dat.Sx(t,s,a);
			deaths = dat.births(t,s) * mort;
			dat.deaths(t,s,a,deaths);

			// ages 1-79
			for (a = DP::AGE_MIN + 1; a < DP::AGE_MAX; ++a) {
				mort = 1.0 - dat.Sx(t,s,a);
				deaths = dat.basepop(s,a-1) * mort;
				dat.deaths(t,s,a,deaths);
			}

			// ages 80+
			a = DP::AGE_MAX;
			mort = 1.0 - dat.Sx(t,s,a);
			deaths = (dat.basepop(s,a-1) + dat.basepop(s,a)) * mort;
			dat.deaths(t,s,a,deaths);
		}
	}

	void Projection::project_one_year(const int t) {
		// We sequence risk calculations before HIV calculations so that we capture HIV
		// risk among those who debut sexually at age 15
		advance_one_year_demography(t);
		advance_one_year_risk(t);
		advance_one_year_male_circumcision(t);
		advance_one_year_hiv(t);
		insert_clhiv_agein(t);
		insert_endyear_migrants(t);
	}

	void Projection::advance_one_year_demography(const int t) {
		double buff[DP::N_HIV_CHILD];
		double surv, mort;
		int a, b, d, h, r, s, u;

		// Projection
		for (u = DP::SEX_MC_MIN; u <= DP::SEX_MC_MAX; ++u) {
			s = sex[u];

			// ages 1-14
			for (a = DP::AGE_CHILD_MIN + 1; a <= DP::AGE_CHILD_MAX; ++a) {
				surv = dat.Sx(t, s, a);
				mort = 1.0 - surv;

				pop.child_neg(t, u, a) = pop.child_neg(t - 1, u, a - 1) * surv;
				dth.child_neg(t, u, a) = pop.child_neg(t - 1, u, a - 1) * mort;

				for (h = DP::HIV_CHILD_MIN; h <= DP::HIV_CHILD_MAX; ++h) {
					for (d = DP::DTX_MIN; d <= DP::DTX_MAX; ++d) {
						pop.child_hiv(t, u, a, h, d) = pop.child_hiv(t - 1, u, a - 1, h, d) * surv;
						dth.child_hiv(t, u, a, h, d) = pop.child_hiv(t - 1, u, a - 1, h, d) * mort;
					}
				}
			}

			// age 15
			a = DP::AGE_ADULT_MIN;
			b = a - DP::AGE_ADULT_MIN;
			r = DP::POP_NOSEX; // Children are assumed sexually inactive before age 15
			surv = dat.Sx(t, s, a);
			mort = 1.0 - surv;

			pop.adult_neg(t, u, b, r) = pop.child_neg(t - 1, u, a - 1) * surv;
			dth.adult_neg(t, u, b, r) = pop.child_neg(t - 1, u, a - 1) * mort;
			for (h = DP::HIV_ADULT_MIN; h <= DP::HIV_ADULT_MAX; ++h) {
				for (d = DP::DTX_MIN; d <= DP::DTX_MAX; ++d) {
					pop.adult_hiv(t, u, b, r, h, d) = pop.child_hiv(t - 1, u, a - 1, h, d) * surv;
					dth.adult_hiv(t, u, b, r, h, d) = pop.child_hiv(t - 1, u, a - 1, h, d) * mort;
				}
			}

			// ages 16-79
			for (a = DP::AGE_ADULT_MIN + 1; a < DP::AGE_ADULT_MAX; ++a) {
				b = a - DP::AGE_ADULT_MIN;
				surv = dat.Sx(t, s, a);
				mort = 1.0 - surv;

				for (r = DP::POP_MIN; r <= DP::POP_MAX; ++r) {
					pop.adult_neg(t, u, b, r) = pop.adult_neg(t - 1, u, b - 1, r) * surv;
					dth.adult_neg(t, u, b, r) = pop.adult_neg(t - 1, u, b - 1, r) * mort;
				}
				for (r = DP::POP_MIN; r <= DP::POP_MAX; ++r)
					for (h = DP::HIV_ADULT_MIN; h <= DP::HIV_ADULT_MAX; ++h)
						for (d = DP::DTX_MIN; d <= DP::DTX_MAX; ++d) {
							pop.adult_hiv(t, u, b, r, h, d) = pop.adult_hiv(t - 1, u, b - 1, r, h, d) * surv;
							dth.adult_hiv(t, u, b, r, h, d) = pop.adult_hiv(t - 1, u, b - 1, r, h, d) * mort;
						}
			}

			// ages 80+
			a = DP::AGE_ADULT_MAX;
			b = a - DP::AGE_ADULT_MIN;
			const double surv_79(dat.Sx(t, s, a));
			const double surv_80(dat.Sx(t, s, a + 1));
			const double mort_79(1.0 - surv_79);
			const double mort_80(1.0 - surv_80);
			for (r = DP::POP_MIN; r <= DP::POP_MAX; ++r) {
				pop.adult_neg(t, u, b, r) = pop.adult_neg(t - 1, u, b - 1, r) * surv_79 + pop.adult_neg(t - 1, u, b, r) * surv_80;
				dth.adult_neg(t, u, b, r) = pop.adult_neg(t - 1, u, b - 1, r) * mort_79 + pop.adult_neg(t - 1, u, b, r) * mort_80;
			}
			for (r = DP::POP_MIN; r <= DP::POP_MAX; ++r)
				for (h = DP::HIV_ADULT_MIN; h <= DP::HIV_ADULT_MAX; ++h)
					for (d = DP::DTX_MIN; d <= DP::DTX_MAX; ++d) {
						pop.adult_hiv(t, u, b, r, h, d) = pop.adult_hiv(t - 1, u, b - 1, r, h, d) * surv_79 + pop.adult_hiv(t - 1, u, b, r, h, d) * surv_80;
						dth.adult_hiv(t, u, b, r, h, d) = pop.adult_hiv(t - 1, u, b - 1, r, h, d) * mort_79 + pop.adult_hiv(t - 1, u, b, r, h, d) * mort_80;
					}

			// redistribute 5 year-olds from CD4 percentages to numbers
			a = 5;
			for (d = DP::DTX_MIN; d <= DP::DTX_MAX; ++d) {
				for (h = DP::HIV_CHILD_MIN; h <= DP::HIV_CHILD_MAX; ++h)
					buff[h] = pop.child_hiv(t, u, a, h, d);
				for (h = DP::HIV_CHILD_MIN; h <= DP::HIV_CHILD_MAX; ++h) {
					pop.child_hiv(t, u, a, h, d) = 0.0;
					for (int i(DP::HIV_CHILD_MIN); i <= DP::HIV_CHILD_MAX; ++i)
						pop.child_hiv(t, u, a, h, d) += CD4_MAP_AGE_5[i][h] * buff[i];
				}
			}

			// redistribute 15 year-olds from child to adult CD4 categories
			a = DP::AGE_ADULT_MIN;
			b = a - DP::AGE_ADULT_MIN;
			for (r = DP::POP_MIN; r <= DP::POP_MAX; ++r) {
				for (d = DP::DTX_MIN; d <= DP::DTX_MAX; ++d) {
					for (h = DP::HIV_ADULT_MIN; h <= DP::HIV_ADULT_MAX; ++h)
						buff[h] = pop.adult_hiv(t, u, b, r, h, d);
					pop.adult_hiv(t, u, b, r, DP::HIV_PRIMARY, d) = 0.00;
					pop.adult_hiv(t, u, b, r, DP::HIV_GEQ_500, d) = buff[DP::HIV_PED_GEQ_1000] + buff[DP::HIV_PED_750_1000] + buff[DP::HIV_PED_500_750];
					pop.adult_hiv(t, u, b, r, DP::HIV_350_500, d) = buff[DP::HIV_PED_350_500];
					pop.adult_hiv(t, u, b, r, DP::HIV_200_350, d) = buff[DP::HIV_PED_200_350];
					pop.adult_hiv(t, u, b, r, DP::HIV_100_200, d) = 0.35 * buff[DP::HIV_PED_LT_200];
					pop.adult_hiv(t, u, b, r, DP::HIV_050_100, d) = 0.21 * buff[DP::HIV_PED_LT_200];
					pop.adult_hiv(t, u, b, r, DP::HIV_050_100, d) = 0.44 * buff[DP::HIV_PED_LT_200];
				}
			}

		}

		calc_births(t);
	}
	
	void Projection::advance_one_year_risk(const int t) {
		const double eps = std::numeric_limits<double>::epsilon(); // padding term to avoid divide-by-zero
		const sex_t umin[DP::N_SEX] = {DP::FEMALE, DP::MALE_U};
		const sex_t umax[DP::N_SEX] = {DP::FEMALE, DP::MALE_C};

		int nt(0), ns(0);
		int kp_turn_sex[DP::N_SEX * DP::N_POP], kp_stay_sex[DP::N_SEX * DP::N_POP];
		int kp_turn_pop[DP::N_SEX * DP::N_POP], kp_stay_pop[DP::N_SEX * DP::N_POP];

		int a, d, h, k, r, s, u;
		boost::multi_array<double, 2> n_total(boost::extents[DP::N_SEX][DP::N_AGE_ADULT]);
		boost::multi_array<double, 3> n_group(boost::extents[DP::N_SEX][DP::N_AGE_ADULT][DP::N_POP]);
		boost::multi_array<double, 3> dneg(boost::extents[DP::N_SEX_MC][DP::N_AGE_ADULT][DP::N_POP]);
		boost::multi_array<double, 5> dhiv(boost::extents[DP::N_SEX_MC][DP::N_AGE_ADULT][DP::N_POP][DP::N_HIV_ADULT][DP::N_DTX]);
		double n_nosex, n_never, n_union, n_split, n_active, size_pool, prop_fert, size_fert(0.0);
		double p_union, p_naive, p_enter, p_leave;
		double kp_need, kp_have;
		double kp_pool[DP::N_SEX][DP::N_AGE_ADULT];
		double p_turn_enter[DP::N_SEX][DP::N_AGE_ADULT]; // % recruitment to key populations with turnover by age
		double p_stay_enter_pop[DP::N_SEX * DP::N_POP];  // % recruitment to key populations without turnover
		double p_stay_enter[DP::N_SEX];

		std::fill_n(dneg.data(), dneg.num_elements(), 0.0);
		std::fill_n(dhiv.data(), dhiv.num_elements(), 0.0);
		std::fill_n(n_group.data(), n_group.num_elements(), eps);
		std::fill_n(n_total.data(), n_total.num_elements(), 0.0);

		// Enumerate key populations with and without turnover ("turn" and "stay", respectively)
		for (s = DP::SEX_MIN; s <= DP::SEX_MAX; ++s)
			for (r = DP::POP_KEY_MIN; r < N_POP_SEX[s]; ++r)
				if (dat.keypop_stay(s, r)) {
					kp_stay_sex[ns] = s;
					kp_stay_pop[ns] = r;
					++ns;
				} else {
					kp_turn_sex[nt] = s;
					kp_turn_pop[nt] = r;
					++nt;
				}

		// Cache population sizes by sex, age, and behavioral risk groups
		for (u = DP::SEX_MC_MIN; u <= DP::SEX_MC_MAX; ++u) {
			s = sex[u];
			for (a = 0; a < DP::N_AGE_ADULT; ++a) {
				for (r = DP::POP_MIN; r < DP::N_POP_SEX[s]; ++r) {
					n_group[s][a][r] += pop.adult_neg(t, u, a, r);
					for (h = DP::HIV_ADULT_MIN; h <= DP::HIV_ADULT_MAX; ++h)
						for (d = DP::DTX_MIN; d <= DP::DTX_MAX; ++d)
							n_group[s][a][r] += pop.adult_hiv(t, u, a, r, h, d);
				}
			}
		}

		for (s = DP::SEX_MIN; s <= DP::SEX_MAX; ++s) {
			for (a = 0; a < DP::N_AGE_ADULT; ++a)
				for (r = DP::POP_MIN; r < DP::N_POP_SEX[s]; ++r)
					n_total[s][a] += n_group[s][a][r];
		}
		
		// Calculate the size of the 15-49 population
		for (s = DP::SEX_MIN; s <= DP::SEX_MAX; ++s) {
			for (a = 0; a < DP::N_AGE_BIRTH; ++a)
				size_fert += n_total[s][a];
		}

		// Calculate the pool of people that new key population members
		// are recruited from
		for (s = DP::SEX_MIN; s <= DP::SEX_MAX; ++s) {
			for (a = 0; a < DP::N_AGE_ADULT; ++a)
				kp_pool[s][a] = n_group[s][a][DP::POP_NOSEX]
				              + n_group[s][a][DP::POP_NEVER]
				              + n_group[s][a][DP::POP_UNION]
				              + n_group[s][a][DP::POP_SPLIT];
		}

		// Initialize cache used to store recruitment to key populations with turnover
		for (s = DP::SEX_MIN; s <= DP::SEX_MAX; ++s) {
			for (a = 0; a < DP::N_AGE_ADULT; ++a)
				p_turn_enter[s][a] = 0.0;
		}

		// key populations with turnover
		for (k = 0; k < nt; ++k) {
			s = kp_turn_sex[k];
			r = kp_turn_pop[k];

			// Since the population size is specified as a % of the 15-49 population,
			// we need the % of key population members are 15-49.
			prop_fert = 0.0;
			for (a = 0; a < DP::N_AGE_BIRTH; ++a)
				prop_fert += dat.keypop_age_dist(s, a, r);

			for (a = 0; a < DP::N_AGE_ADULT; ++a) {			
				kp_need = size_fert * dat.keypop_size(s, r) * dat.keypop_age_dist(s, a, r) / prop_fert; // number we should have
				kp_have = n_group[s][a][r] * (1.0 - dat.keypop_exit_prop(s, r));                        // number we would have after turnover if no new recruits or extra exits

				if (kp_need > kp_have) { // need more recruits
					p_enter = (kp_need - kp_have) / kp_pool[s][a];
					p_leave = dat.keypop_exit_prop(s, r);
				} else {                 // need more exits
					p_enter = 0.0;
					p_leave = dat.keypop_exit_prop(s, r) + (kp_have - kp_need) / n_group[s][a][r];
				}
				p_turn_enter[s][a] += p_enter;

				// Calculate the changes in key population size and contribution of turnover to general population
				for (u = umin[s]; u <= umax[s]; ++u) {
					n_split = pop.adult_neg(t, u, a, DP::POP_SPLIT);
					n_union = pop.adult_neg(t, u, a, DP::POP_UNION);
					p_union = n_union / (eps + n_union + n_split);
					size_pool = pop.adult_neg(t, u, a, DP::POP_NOSEX)
						        + pop.adult_neg(t, u, a, DP::POP_NEVER)
						        + pop.adult_neg(t, u, a, DP::POP_UNION)
						        + pop.adult_neg(t, u, a, DP::POP_SPLIT);
					dneg[u][a][r] += p_enter * size_pool - p_leave * pop.adult_neg(t, u, a, r);
					dneg[u][a][DP::POP_UNION] += p_leave * pop.adult_neg(t, u, a, r) * p_union;
					dneg[u][a][DP::POP_SPLIT] += p_leave * pop.adult_neg(t, u, a, r) * (1.0 - p_union);
					for (h = DP::HIV_ADULT_MIN; h <= DP::HIV_ADULT_MAX; ++h) {
						for (d = DP::DTX_MIN; d <= DP::DTX_MAX; ++d) {
							n_split = pop.adult_hiv(t, u, a, DP::POP_SPLIT, h, d);
							n_union = pop.adult_hiv(t, u, a, DP::POP_UNION, h, d);
							p_union = n_union / (eps + n_union + n_split);
							size_pool = pop.adult_hiv(t, u, a, DP::POP_NOSEX, h, d)
								        + pop.adult_hiv(t, u, a, DP::POP_NEVER, h, d)
								        + pop.adult_hiv(t, u, a, DP::POP_UNION, h, d)
								        + pop.adult_hiv(t, u, a, DP::POP_SPLIT, h, d);
							dhiv[u][a][r][h][d] += p_enter * size_pool - p_leave * pop.adult_hiv(t, u, a, r, h, d);
							dhiv[u][a][DP::POP_UNION][h][d] += p_leave * pop.adult_hiv(t, u, a, r, h, d) * p_union;
							dhiv[u][a][DP::POP_SPLIT][h][d] += p_leave * pop.adult_hiv(t, u, a, r, h, d) * (1.0 - p_union);
						}
					}
				}
			}
		}

		// Calculate the proportion of people who enter key populations without turnover
		// at sexual debut. If we set this to p=dat.keypop_size(s,r), then p*100% of
		// sexually active people would be in the key population, but since not everyone
		// is sexually active, that population size would be too small. The calculation
		// below accounts for the % of the 15-49 population that is sexually active.
		p_stay_enter[DP::FEMALE] = 0.0;
		p_stay_enter[DP::MALE  ] = 0.0;
		for (k = 0; k < ns; ++k) {
			s = kp_stay_sex[k];
			r = kp_stay_pop[k];
			
			p_naive = 1.0;
			n_active = 0.0;
			for (a = 0; a < DP::N_AGE_BIRTH; ++a) {
				p_naive *= 1.0 - dat.debut_prop(s);
				n_active += n_total[s][a] * (1.0 - p_naive);
			}
			p_stay_enter_pop[k] = dat.keypop_size(s, r) * size_fert / n_active;
			p_stay_enter[s] += p_stay_enter_pop[k];
		}

		for (u = DP::SEX_MC_MIN; u <= DP::SEX_MC_MAX; ++u) {
			s = sex[u];
			for (a = 0; a < DP::N_AGE_ADULT; ++a) {
				n_nosex = pop.adult_neg(t, u, a, DP::POP_NOSEX);
				n_never = pop.adult_neg(t, u, a, DP::POP_NEVER);
				n_union = pop.adult_neg(t, u, a, DP::POP_UNION);
				n_split = pop.adult_neg(t, u, a, DP::POP_SPLIT);
				dneg[u][a][DP::POP_NOSEX] -= n_nosex * (dat.debut_prop(s) + p_turn_enter[s][a]);
				dneg[u][a][DP::POP_NEVER] += n_nosex * dat.debut_prop(s) * (1.0 - dat.union_prop(s) - p_stay_enter[s]);
				dneg[u][a][DP::POP_NEVER] -= n_never * (dat.union_prop(s) + p_turn_enter[s][a]);
				dneg[u][a][DP::POP_UNION] += n_nosex * dat.debut_prop(s) * dat.union_prop(s);
				dneg[u][a][DP::POP_UNION] += (n_never + n_split) * dat.union_prop(s);
				dneg[u][a][DP::POP_UNION] -= n_union * (dat.split_prop() + p_turn_enter[s][a]);
				dneg[u][a][DP::POP_SPLIT] += n_union * dat.split_prop();
				dneg[u][a][DP::POP_SPLIT] -= n_split * (dat.union_prop(s) + p_turn_enter[s][a]);
				for (h = DP::HIV_ADULT_MIN; h <= DP::HIV_ADULT_MAX; ++h) {
					for (d = DP::DTX_MIN; d <= DP::DTX_MAX; ++d) {
						n_nosex = pop.adult_hiv(t, u, a, DP::POP_NOSEX, h, d);
						n_never = pop.adult_hiv(t, u, a, DP::POP_NEVER, h, d);
						n_union = pop.adult_hiv(t, u, a, DP::POP_UNION, h, d);
						n_split = pop.adult_hiv(t, u, a, DP::POP_SPLIT, h, d);
						dhiv[u][a][DP::POP_NOSEX][h][d] -= n_nosex * (dat.debut_prop(s) + p_turn_enter[s][a]);
						dhiv[u][a][DP::POP_NEVER][h][d] += n_nosex * dat.debut_prop(s) * (1.0 - dat.union_prop(s) - p_stay_enter[s]);
						dhiv[u][a][DP::POP_NEVER][h][d] -= n_never * (dat.union_prop(s) + p_turn_enter[s][a]);
						dhiv[u][a][DP::POP_UNION][h][d] += n_nosex * dat.debut_prop(s) * dat.union_prop(s);
						dhiv[u][a][DP::POP_UNION][h][d] += (n_never + n_split) * dat.union_prop(s);
						dhiv[u][a][DP::POP_UNION][h][d] -= n_union * (dat.split_prop() + p_turn_enter[s][a]);
						dhiv[u][a][DP::POP_SPLIT][h][d] += n_union * dat.split_prop();
						dhiv[u][a][DP::POP_SPLIT][h][d] -= n_split * (dat.union_prop(s) + p_turn_enter[s][a]);
					}
				}
			}
		}

		// key populations without turnover
		for (k = 0; k < ns; ++k) {
			s = kp_stay_sex[k];
			r = kp_stay_pop[k];
			for (u = umin[s]; u <= umax[s]; ++u) {
				for (a = 0; a < DP::N_AGE_ADULT; ++a) {
					dneg[u][a][r] += pop.adult_neg(t, u, a, DP::POP_NOSEX) * dat.debut_prop(s) * p_stay_enter_pop[k];
					for (h = DP::HIV_ADULT_MIN; h <= DP::HIV_ADULT_MAX; ++h) {
						for (d = DP::DTX_MIN; d <= DP::DTX_MAX; ++d) {
							dhiv[u][a][r][h][d] += pop.adult_hiv(t, u, a, DP::POP_NOSEX, h, d) * dat.debut_prop(s) * p_stay_enter_pop[k];
						}
					}
				}
			}
		}

		// Finalize the changes.
		for (u = DP::SEX_MC_MIN; u <= DP::SEX_MC_MAX; ++u) {
			s = sex[u];
			for (a = 0; a < DP::N_AGE_ADULT; ++a)
				for (r = DP::POP_MIN; r < DP::N_POP_SEX[s]; ++r) {
					pop.adult_neg(t, u, a, r) += dneg[u][a][r];
					for (h = DP::HIV_ADULT_MIN; h <= DP::HIV_ADULT_MAX; ++h)
						for (d = DP::DTX_MIN; d <= DP::DTX_MAX; ++d)
							pop.adult_hiv(t, u, a, r, h, d) += dhiv[u][a][r][h][d];
				}
		}

	}

	void Projection::advance_one_year_male_circumcision(const int t) {
		int a, b, r, h, d;
		double puptake, nuptake;

		// Children
		for (a = DP::AGE_CHILD_MIN; a <= DP::AGE_CHILD_MAX; ++a) {
			puptake = dat.uptake_male_circumcision(t, a);
			nuptake = pop.child_neg(t, DP::MALE_U, a) * puptake;
			pop.child_neg(t, DP::MALE_U, a) -= nuptake;
			pop.child_neg(t, DP::MALE_C, a) += nuptake;
			for (h = DP::HIV_CHILD_MIN; h <= DP::HIV_CHILD_MAX; ++h) {
				for (d = DP::DTX_MIN; d <= DP::DTX_MAX; ++d) {
					nuptake = pop.child_hiv(t, DP::MALE_U, a, h, d) * puptake;
					pop.child_hiv(t, DP::MALE_U, a, h, d) -= nuptake;
					pop.child_hiv(t, DP::MALE_C, a, h, d) += nuptake;
				}
			}
		}

		// Adults
		for (a = DP::AGE_ADULT_MIN; a <= DP::AGE_ADULT_MAX; ++a) {
			b = a - DP::AGE_ADULT_MIN;
			puptake = dat.uptake_male_circumcision(t, a);
			for (r = DP::POP_MIN; r <= DP::POP_MAX; ++r) {
				nuptake = pop.adult_neg(t, DP::MALE_U, b, r) * puptake;
				pop.adult_neg(t, DP::MALE_U, b, r) -= nuptake;
				pop.adult_neg(t, DP::MALE_C, b, r) += nuptake;
				for (h = DP::HIV_ADULT_MIN; h <= DP::HIV_ADULT_MAX; ++h) {
					for (d = DP::DTX_MIN; d <= DP::DTX_MAX; ++d) {
						nuptake = pop.adult_hiv(t, DP::MALE_U, b, r, h, d) * puptake;
						pop.adult_hiv(t, DP::MALE_U, b, r, h, d) -= nuptake;
						pop.adult_hiv(t, DP::MALE_C, b, r, h, d) += nuptake;
					}
				}
			}
		}

	}

	void Projection::advance_one_year_hiv(const int time) {
		advance_one_year_hiv_adult(time);
		advance_one_year_hiv_child(time);
	}

	void Projection::advance_one_year_hiv_adult(const int time) {
		if (!dat.direct_incidence() && time == dat.seed_time()) {
			seed_epidemic(time, dat.seed_prevalence());
		}

		for (int step(0); step < DP::HIV_TIME_STEPS; ++step) {
			advance_one_step_hiv_adult(time, step);
			if (dat.direct_incidence()) {
				insert_adult_infections(time, step);
			} else {
				if (time >= dat.seed_time())
					calc_adult_infections(time, step);
			}
		}
	}

	void Projection::advance_one_year_hiv_child(const int t) {
		double births_exposed(calc_births_hiv(t));
		dat.births_hiv_exposed(t, births_exposed);
		// TODO: pediatric HIV infection and progression calculations
	}

	void Projection::advance_one_step_hiv_adult(const int t, const int step) {
		const double eps = std::numeric_limits<double>::epsilon(); // padding to avoid divide-by-zero

		int s, u, a, b, r, h, d;
		double influx[DP::N_HIV_ADULT][DP::N_DTX], efflux[DP::N_HIV_ADULT][DP::N_DTX];
		double art_exit[DP::N_HIV_ADULT];
		double prog_primary, off_art;
		double num_art[DP::N_SEX], num_off[DP::N_SEX];
		double art_mort_scale[DP::N_SEX][DP::N_AGE_ADULT][DP::N_HIV_ADULT];
		sex_hiv_t uptake_rate(boost::extents[DP::N_SEX][DP::N_HIV_ADULT]);

		calc_adult_art_uptake(t, step, uptake_rate);

		// Calculate scale factors for adjusting HIV-related mortality off ART based
		// on ART coverage. We sum over behavioral risk groups, male circumcision
		// status, and different off-ART compartment sizes so that these reductions
		// align with Spectrum.
		for (b = 0; b < DP::N_AGE_ADULT; ++b) {
			for (h = DP::HIV_ADULT_MIN; h <= DP::HIV_ADULT_MAX; ++h) {
				num_art[DP::MALE] = num_art[DP::FEMALE] = 0.0;
				num_off[DP::MALE] = num_off[DP::FEMALE] = 0.0;
				for (r = DP::POP_MIN; r <= DP::POP_MAX; ++r) {
					for (d = DP::DTX_ART_MIN; d <= DP::DTX_ART_MAX; ++d) {
						num_art[DP::MALE  ] += pop.adult_hiv(t, DP::MALE_U, b, r, h, d) + pop.adult_hiv(t, DP::MALE_C, b, r, h, d);
						num_art[DP::FEMALE] += pop.adult_hiv(t, DP::FEMALE, b, r, h, d);
					}
					for (d = DP::DTX_OFF_MIN; d <= DP::DTX_OFF_MAX; ++d) {
						num_off[DP::MALE  ] += pop.adult_hiv(t, DP::MALE_U, b, r, h, d) + pop.adult_hiv(t, DP::MALE_C, b, r, h, d);
						num_off[DP::FEMALE] += pop.adult_hiv(t, DP::FEMALE, b, r, h, d);
					}
				}
				for (s = DP::SEX_MIN; s <= DP::SEX_MAX; ++s)
					art_mort_scale[s][b][h] = 1.0 - num_art[s] / (num_art[s] + num_off[s] + eps);
			}			
		}

		for (u = DP::SEX_MC_MIN; u <= DP::SEX_MC_MAX; ++u) {
			s = sex[u];
			for (a = DP::AGE_ADULT_MIN; a <= DP::AGE_ADULT_MAX; ++a) {
				b = a - DP::AGE_ADULT_MIN;
				for (r = DP::POP_MIN; r <= DP::POP_MAX; ++r) {

					// Buffer patients who interrupt ART. ART_EXIT_STAGE maps
					// from baseline HIV stage h at ART initiation and ART duration d
					// to HIV stage after interruption.
					for (h = DP::HIV_MIN; h <= DP::HIV_MAX; ++h)
						art_exit[h] = 0.0;
					for (h = DP::HIV_MIN; h <= DP::HIV_MAX; ++h) {
						for (d = DP::DTX_ART_MIN; d <= DP::DTX_ART_MAX; ++d)
							art_exit[ART_EXIT_STAGE[h][d]] += pop.adult_hiv(t, u, b, r, h, d) * dat.art_exit_adult(t, s);
					}

					// Calculate flows
					// Not on ART
					// TODO: implement reductions in off-ART mortality proportional to ART coverage as in AIM
					for (d = DP::DTX_UNAWARE; d <= DP::DTX_PREV_TX; ++d) {
#ifndef SPECTRUM_CD4
						prog_primary = pop.adult_hiv(t, u, b, r, HIV_PRIMARY, d) * dat.hiv_prog(s, a, HIV_PRIMARY);

						influx[HIV_PRIMARY][d] = 0.0;
						influx[HIV_GEQ_500][d] = prog_primary * dat.hiv_dist(s, a, HIV_GEQ_500);
						influx[HIV_350_500][d] = prog_primary * dat.hiv_dist(s, a, HIV_350_500) + pop.adult_hiv(t, u, b, r, HIV_GEQ_500, d) * dat.hiv_prog(s, a, HIV_GEQ_500);
						influx[HIV_200_350][d] = prog_primary * dat.hiv_dist(s, a, HIV_200_350) + pop.adult_hiv(t, u, b, r, HIV_350_500, d) * dat.hiv_prog(s, a, HIV_350_500);
						influx[HIV_100_200][d] = prog_primary * dat.hiv_dist(s, a, HIV_100_200) + pop.adult_hiv(t, u, b, r, HIV_200_350, d) * dat.hiv_prog(s, a, HIV_200_350);
						influx[HIV_050_100][d] = prog_primary * dat.hiv_dist(s, a, HIV_050_100) + pop.adult_hiv(t, u, b, r, HIV_100_200, d) * dat.hiv_prog(s, a, HIV_100_200);
						influx[HIV_000_050][d] = prog_primary * dat.hiv_dist(s, a, HIV_000_050) + pop.adult_hiv(t, u, b, r, HIV_050_100, d) * dat.hiv_prog(s, a, HIV_050_100);

						efflux[HIV_PRIMARY][d] = pop.adult_hiv(t, u, b, r, HIV_PRIMARY, d) * (dat.hiv_prog(s, a, HIV_PRIMARY) + art_mort_scale[s][b][HIV_PRIMARY] * dat.hiv_mort(s, a, HIV_PRIMARY));
						efflux[HIV_GEQ_500][d] = pop.adult_hiv(t, u, b, r, HIV_GEQ_500, d) * (dat.hiv_prog(s, a, HIV_GEQ_500) + art_mort_scale[s][b][HIV_GEQ_500] * dat.hiv_mort(s, a, HIV_GEQ_500));
						efflux[HIV_350_500][d] = pop.adult_hiv(t, u, b, r, HIV_350_500, d) * (dat.hiv_prog(s, a, HIV_350_500) + art_mort_scale[s][b][HIV_350_500] * dat.hiv_mort(s, a, HIV_350_500));
						efflux[HIV_200_350][d] = pop.adult_hiv(t, u, b, r, HIV_200_350, d) * (dat.hiv_prog(s, a, HIV_200_350) + art_mort_scale[s][b][HIV_200_350] * dat.hiv_mort(s, a, HIV_200_350));
						efflux[HIV_100_200][d] = pop.adult_hiv(t, u, b, r, HIV_100_200, d) * (dat.hiv_prog(s, a, HIV_100_200) + art_mort_scale[s][b][HIV_100_200] * dat.hiv_mort(s, a, HIV_100_200));
						efflux[HIV_050_100][d] = pop.adult_hiv(t, u, b, r, HIV_050_100, d) * (dat.hiv_prog(s, a, HIV_050_100) + art_mort_scale[s][b][HIV_050_100] * dat.hiv_mort(s, a, HIV_050_100));
						efflux[HIV_000_050][d] = pop.adult_hiv(t, u, b, r, HIV_000_050, d) * dat.hiv_mort(s, a, HIV_000_050);
#else
						influx[DP::HIV_ADULT_MIN][d] = 0.0;
						for (h = DP::HIV_ADULT_MIN + 1; h <= DP::HIV_ADULT_MAX; ++h)
							influx[h][d] = pop.adult_hiv(t, u, b, r, h - 1, d) * dat.hiv_prog(s, a, h - 1);

						for (h = DP::HIV_ADULT_MIN; h < DP::HIV_ADULT_MAX; ++h)
							efflux[h][d] = pop.adult_hiv(t, u, b, r, h, d) * (dat.hiv_prog(s, a, h) + art_mort_scale[s][b][h] * dat.hiv_mort(s, a, h));
						efflux[DP::HIV_ADULT_MAX][d] = pop.adult_hiv(t, u, b, r, DP::HIV_ADULT_MAX, d) * art_mort_scale[s][b][h] * dat.hiv_mort(s, a, DP::HIV_ADULT_MAX);
#endif

						for (h = DP::HIV_ADULT_MIN; h <= DP::HIV_ADULT_MAX; ++h)
							dth.adult_hiv(t, u, b, r, h, d) += DP::HIV_STEP_SIZE * pop.adult_hiv(t, u, b, r, h, d) * art_mort_scale[s][b][h] * dat.hiv_mort(s, a, h);
					}

					// Adjust for ART interruption
					for (h = DP::HIV_MIN; h <= DP::HIV_MAX; ++h) {
						influx[h][DP::DTX_PREV_TX] += art_exit[h];
					}

					// Adjust for ART uptake
					for (h = DP::HIV_MIN; h <= DP::HIV_MAX; ++h) {
						for (d = DP::DTX_OFF_MIN; d <= DP::DTX_OFF_MAX; ++d) {
							efflux[h][d] += pop.adult_hiv(t, u, b, r, h, d) * uptake_rate[s][h];
						}
					}

					// On ART
					for (h = DP::HIV_MIN; h <= DP::HIV_MAX; ++h) {
						off_art = pop.adult_hiv(t, u, b, r, h, DP::DTX_UNAWARE)
							      + pop.adult_hiv(t, u, b, r, h, DP::DTX_AWARE)
							      + pop.adult_hiv(t, u, b, r, h, DP::DTX_PREV_TX);

						influx[h][DP::DTX_ART1] = off_art * uptake_rate[s][h];
						influx[h][DP::DTX_ART2] = pop.adult_hiv(t, u, b, r, h, DP::DTX_ART1) * dat.art_flow(DP::DTX_ART1);
						influx[h][DP::DTX_ART3] = pop.adult_hiv(t, u, b, r, h, DP::DTX_ART2) * dat.art_flow(DP::DTX_ART2);

						efflux[h][DP::DTX_ART1] = pop.adult_hiv(t, u, b, r, h, DP::DTX_ART1) * (dat.art_exit_adult(t, s) + dat.art_mort_adult(t, s, b, h, DP::DTX_ART1) + dat.art_flow(DP::DTX_ART1));
						efflux[h][DP::DTX_ART2] = pop.adult_hiv(t, u, b, r, h, DP::DTX_ART2) * (dat.art_exit_adult(t, s) + dat.art_mort_adult(t, s, b, h, DP::DTX_ART2) + dat.art_flow(DP::DTX_ART2));
						efflux[h][DP::DTX_ART3] = pop.adult_hiv(t, u, b, r, h, DP::DTX_ART3) * (dat.art_exit_adult(t, s) + dat.art_mort_adult(t, s, b, h, DP::DTX_ART3));

						for (d = DP::DTX_ART_MIN; d <= DP::DTX_ART_MAX; ++d)
							dth.adult_hiv(t, u, b, r, h, d) += DP::HIV_STEP_SIZE * pop.adult_hiv(t, u, b, r, h, d) * dat.art_mort_adult(t, s, b, h, d);
					}

					// Implement flows
					for (h = DP::HIV_ADULT_MIN; h <= DP::HIV_ADULT_MAX; ++h)
						for (d = DP::DTX_MIN; d <= DP::DTX_MAX; ++d)
							pop.adult_hiv(t, u, b, r, h, d) += DP::HIV_STEP_SIZE * (influx[h][d] - efflux[h][d]);
				}
			}
		}
	}

	void Projection::calc_adult_art_uptake(const int t, const int step, sex_hiv_t& uptake_rate) {
		std::fill_n(uptake_rate.data(), uptake_rate.num_elements(), 0.0);

		// Short-circuit uptake calculations if no one is on ART
		if (dat.art_num_adult( t, DP::FEMALE) == 0.0 && dat.art_num_adult( t, DP::MALE) == 0.0 &&
			  dat.art_prop_adult(t, DP::FEMALE) == 0.0 && dat.art_prop_adult(t, DP::MALE) == 0.0) {
			return;
		}

		const double eps = std::numeric_limits<double>::epsilon(); // padding to avoid divide-by-zero
		const int elig_first = dat.art_first_eligible_stage_adult(t);
		const double wgt_mort = dat.art_mort_weight();
		const double wgt_elig = 1.0 - dat.art_mort_weight();

		double eligible[DP::N_SEX] = {0.0, 0.0}; // eligible, on or off ART
		double elig_off[DP::N_SEX] = {0.0, 0.0}; // eligible, off ART
		double retained[DP::N_SEX] = {0.0, 0.0}; // anticipated retention at next timestep
		double uptake[DP::N_SEX] = {0.0, 0.0};
		double art_input[2][DP::N_SEX]; // art_input[k] is the ART input k = 0 or 1 years ago
		double elig_cd4[DP::N_SEX][DP::N_HIV_ADULT]; // number eligible off ART by sex and CD4 count
		double mort_cd4[DP::N_SEX][DP::N_HIV_ADULT]; // HIV-related deaths off ART by sex and CD4 count
		double init_cd4[DP::N_SEX][DP::N_HIV_ADULT]; // ART uptake by CD4 count
		double prop_mort[DP::N_HIV_ADULT], prop_elig[DP::N_HIV_ADULT];
		double norm_mort, norm_elig;
		double loss_rate, target, remaining;
		int a, b, d, h, k, r, s, u;

		for (s = DP::SEX_MIN; s <= DP::SEX_MAX; ++s) {
			for (h = DP::HIV_ADULT_MIN; h <= DP::HIV_ADULT_MAX; ++h) {
				elig_cd4[s][h] = 0.0;
				mort_cd4[s][h] = 0.0;
				init_cd4[s][h] = 0.0;
			}
		}

		// Count the number eligible and project retention
		for (u = DP::SEX_MC_MIN; u <= DP::SEX_MC_MAX; ++u) {
			s = sex[u];
			for (a = DP::AGE_ADULT_MIN; a <= DP::AGE_ADULT_MAX; ++a) {
				b = a - DP::AGE_ADULT_MIN;
				for (r = DP::POP_MIN; r <= DP::POP_MAX; ++r) {
					for (h = elig_first; h <= DP::HIV_ADULT_MAX; ++h) {
						for (d = DP::DTX_OFF_MIN; d <= DP::DTX_OFF_MAX; ++d) {
							elig_cd4[s][h] += pop.adult_hiv(t, u, b, r, h, d);
							mort_cd4[s][h] += pop.adult_hiv(t, u, b, r, h, d) * dat.hiv_mort(s, a, h);
						}
						for (d = DP::DTX_ART_MIN; d <= DP::DTX_ART_MAX; ++d) {
							loss_rate = (dat.art_exit_adult(t,s) + dat.art_mort_adult(t,s,b,h,d)) * DP::HIV_STEP_SIZE;
							retained[s] += pop.adult_hiv(t, u, b, r, h, d) * (1.0 - loss_rate);
							eligible[s] += pop.adult_hiv(t, u, b, r, h, d);
						}
					}
				}
			}
		}

		for (s = DP::SEX_MIN; s <= DP::SEX_MAX; ++s) {
			for (h = elig_first; h <= DP::HIV_ADULT_MAX; ++h)
				elig_off[s] += elig_cd4[s][h];
			eligible[s] += elig_off[s];
		}

		// We linearly interpolate target numbers on ART by time step. We use input
		// proportions when available and input numbers otherwise.
		k = t > 0 ? t - 1 : t;
		for (s = DP::SEX_MIN; s <= DP::SEX_MAX; ++s) {
			art_input[1][s] = (dat.art_prop_adult(k, s) > 0.0 ? dat.art_prop_adult(k, s) * eligible[s] : dat.art_num_adult(k, s));
			art_input[0][s] = (dat.art_prop_adult(t, s) > 0.0 ? dat.art_prop_adult(t, s) * eligible[s] : dat.art_num_adult(t, s));
		}

		// Calculate ART uptake as the difference between a target number and expected
		// number retained on ART. Targets depend on whether inputs are percentages or
		// numbers. For percentages, the target is a proportion of the distance to that
		// input coverage level. For absolute numbers, we interpolate the target from
		// consecutive Dec. 31 inputs.
		for (s = DP::SEX_MIN; s <= DP::SEX_MAX; ++s) {
			// This mechanism is based on Spectrum, where steps are indexed 1..10. Goals
			// ARM indexes steps by 0..9, so we add 1 to match.
			if (dat.art_prop_adult(t,s) > 0.0) { // Input in percentages
				target = retained[s] + (art_input[0][s] - retained[s]) * DP::HIV_STEP_SIZE * (step + 1);
			} else {                             // Input in absolute numbers
				target = art_input[1][s] + (art_input[0][s] - art_input[1][s]) * DP::HIV_STEP_SIZE * (step + 1);
			}

			// uptake[s] = target - retained[s] ideally, but cannot be negative and
			// cannot exceed the number of untreated eligible people
			uptake[s] = std::clamp(target - retained[s], 0.0, elig_off[s]);
		}

		// Calculate uptake by CD4 category
		for (s = DP::SEX_MIN; s <= DP::SEX_MAX; ++s) {
			norm_elig = std::accumulate(elig_cd4[s], elig_cd4[s] + DP::N_HIV_ADULT, eps);
			for (h = DP::HIV_ADULT_MIN; h <= DP::HIV_ADULT_MAX; ++h)
				prop_elig[h] = elig_cd4[s][h] / norm_elig;

			// Uptake based on expected mortality could suggest more people initiate
			// ART than are present in a compartment if the mortality rate exceeds 1/year.
			// The incremental calculation here ensures that if the mortality rate for
			// a given category exceeds 1, the resulting excess initiations roll over to
			// higher CD4 categories.
			remaining = uptake[s];
			for (h = DP::HIV_ADULT_MAX; h >= DP::HIV_ADULT_MIN; --h) {
				norm_mort = std::accumulate(mort_cd4[s], mort_cd4[s] + h + 1, 0.0);
				init_cd4[s][h] = norm_mort > 0 ? std::min(remaining * mort_cd4[s][h] / norm_mort, elig_cd4[s][h]) : 0.0;
				remaining -= init_cd4[s][h];
			}

			if (uptake[s] > 0.0) {
				for (h = DP::HIV_ADULT_MIN; h <= DP::HIV_ADULT_MAX; ++h) {
					prop_mort[h] = init_cd4[s][h] / uptake[s];
				}
			} else {
				for (h = DP::HIV_ADULT_MIN; h <= DP::HIV_ADULT_MAX; ++h) {
					prop_mort[h] = 0.0;
				}
			}

			for (h = DP::HIV_ADULT_MIN; h <= DP::HIV_ADULT_MAX; ++h) {
				init_cd4[s][h] = (wgt_mort * prop_mort[h] + wgt_elig * prop_elig[h]) * uptake[s];
			}
		}

		// Return an annualized uptake rate
		for (s = DP::SEX_MIN; s <= DP::SEX_MAX; ++s)
				for (h = elig_first; h <= DP::HIV_ADULT_MAX; ++h)
					uptake_rate[s][h] = (elig_cd4[s][h] > 0.0 ? HIV_TIME_STEPS * init_cd4[s][h] / elig_cd4[s][h] : 0.0);
	}

	void Projection::calc_adult_infections(const int t, const int step) {
		// TODO: This is quite slow. Can we approximate this well by doing calculations by age groups?
		// TODO: add code to precalculate quantities that are constant throughout a year (transmission probabilities, maybe mixing and balancing are close enough?)
		// TODO: needle-based transmission

		const double eps = std::numeric_limits<double>::epsilon(); // padding to avoid divide-by-zero
		const double e_condom = dat.effect_condom();

		int si, bi, ri, sj, bj, rj; // i refers to HIV- partner, j to the HIV+
		int ai, ui, uj, cj, hj, dj, vj, qij;

		double prop_transmit, vmmc_mult, new_hiv;
		double bal_numer, bal_denom;
		double bal_mix, per_act, acts_with, acts_wout, num_art, force_group;
		double canmix_numer[DP::N_SEX][DP::N_POP], prefer_numer[DP::N_SEX][DP::N_POP];
		double force[DP::N_SEX][DP::N_AGE_ADULT][DP::N_POP];
		double popsize[DP::N_SEX][DP::N_AGE_ADULT][DP::N_POP];
		double prev[DP::N_SEX][DP::N_AGE_ADULT][DP::N_POP][DP::N_STAGE][DP::N_VL];

		// We calculate transmission in heterosexual marital or cohabiting "unions" separately from "other"
		// partnerships that include same sex, casual, or commercial sexual partnerships. We make this distinction to
		// ensure unions and other partnerships are both balanced
		double mix_other, bal_other, canmix_denom, prefer_denom, assort;
		double supply_other[DP::N_SEX][DP::N_AGE_ADULT][DP::N_POP];
		double supply_pop_other[DP::N_SEX][DP::N_POP];
		double mix_pop_other[DP::N_SEX][DP::N_POP][DP::N_SEX][DP::N_POP];
		double force_other[DP::N_SEX][DP::N_AGE_ADULT][DP::N_POP];

		double mix_union, bal_union, union_denom;
		double prop_union[DP::N_SEX][DP::N_POP];
		double supply_union[DP::N_SEX][DP::N_AGE_ADULT][DP::N_POP];
		double supply_pop_union[DP::N_SEX][DP::N_POP];
		double mix_pop_union[DP::N_SEX][DP::N_POP][DP::N_SEX][DP::N_POP];
		double force_union[DP::N_SEX][DP::N_AGE_ADULT][DP::N_POP];

		// Transmission probability per partnership ptransmit[si][sj][qij][hj][dj]
		// si  HIV- partner sex
		// sj  HIV+ partner sex
		// qij Partnership type
		// hj  HIV+ partner infection stage
		// vj  HIV+ partner viral load
		// We assume transmission risk is independent of age after adjusting for the strata above
		double ptransmit[DP::N_SEX][DP::N_SEX][DP::N_BOND][DP::N_STAGE][DP::N_VL];

		if (step == 0) { // initialization at first step of year
			for (ui = 0; ui < DP::N_SEX_MC; ++ui) {
				for (bi = 0; bi < DP::N_AGE_ADULT; ++bi) {
					ai = bi + DP::AGE_ADULT_MIN;
					for (ri = 0; ri < DP::N_POP; ++ri)
						dat.new_hiv_infections(t, ui, ai, ri, 0.0);
				}
			}
		}

		for (si = DP::SEX_MIN; si <= DP::SEX_MAX; ++si) {
			prop_union[si][DP::POP_NEVER] = 0.0;
			prop_union[si][DP::POP_UNION] = 1.0;
			prop_union[si][DP::POP_SPLIT] = 0.0;
			prop_union[si][DP::POP_PWID ] = dat.keypop_married(si, DP::POP_PWID );
			prop_union[si][DP::POP_BOTH ] = dat.keypop_married(si, DP::POP_BOTH );
		}
		prop_union[DP::MALE][DP::POP_MSM] = dat.keypop_married(DP::MALE, DP::POP_MSM);
		prop_union[DP::MALE][DP::POP_TGW] = dat.keypop_married(DP::MALE, DP::POP_TGW);

		// calculate population sizes
		for (si = DP::SEX_MIN; si <= DP::SEX_MAX; ++si) {
			for (bi = 0; bi < DP::N_AGE_ADULT; ++bi)
				for (ri = DP::POP_MIN; ri < DP::N_POP_SEX[si]; ++ri)
					popsize[si][bi][ri] = 0.0;
		}

		for (uj = DP::SEX_MC_MIN; uj <= DP::SEX_MC_MAX; ++uj) {
			sj = sex[uj];
			for (bj = 0; bj < DP::N_AGE_ADULT; ++bj)
				for (rj = DP::POP_NEVER; rj < DP::N_POP_SEX[sj]; ++rj) {
					popsize[sj][bj][rj] += pop.adult_neg(t, uj, bj, rj);
					for (cj = DP::HIV_ADULT_MIN; cj <= DP::HIV_ADULT_MAX; ++cj)
						for (dj = DP::DTX_MIN; dj <= DP::DTX_MAX; ++dj)
							popsize[sj][bj][rj] += pop.adult_hiv(t, uj, bj, rj, cj, dj);
				}
		}

		// calculate partnership supply
		for (si = DP::SEX_MIN; si <= DP::SEX_MAX; ++si) {
			for (bi = 0; bi < DP::N_AGE_ADULT; ++bi)
				for (ri = DP::POP_NEVER; ri < DP::N_POP_SEX[si]; ++ri) {
					supply_other[si][bi][ri] = popsize[si][bi][ri] * dat.partner_rate(t, si, bi, ri);
					supply_union[si][bi][ri] = popsize[si][bi][ri] * prop_union[si][ri];
				}
		}

		for (si = DP::SEX_MIN; si <= DP::SEX_MAX; ++si) {
			for (ri = DP::POP_NEVER; ri < DP::N_POP_SEX[si]; ++ri) {
				supply_pop_other[si][ri] = 0.0;
				supply_pop_union[si][ri] = 0.0;
				for (bi = 0; bi < DP::N_AGE_ADULT; ++bi) {
					supply_pop_other[si][ri] += supply_other[si][bi][ri];
					supply_pop_union[si][ri] += supply_union[si][bi][ri];
				}
			}
		}

		// calculate mixing coefficient factors by behavioral risk group
		for (si = DP::SEX_MIN; si <= DP::SEX_MAX; ++si) {
			for (ri = DP::POP_NEVER; ri < DP::N_POP_SEX[si]; ++ri) {
				assort = dat.partner_assortativity(si, ri);
				canmix_denom = eps;
				prefer_denom = eps;
				union_denom = eps;
				for (sj = DP::SEX_MIN; sj <= DP::SEX_MAX; ++sj) {
					for (rj = DP::POP_NEVER; rj < DP::N_POP_SEX[sj]; ++rj) {
						canmix_numer[sj][rj] = supply_pop_other[sj][rj] * (dat.mix_structure(si, ri, sj, rj) > 0); // groups can mix
						prefer_numer[sj][rj] = supply_pop_other[sj][rj] * (dat.mix_structure(si, ri, sj, rj) > 1); // groups prefer to mix
						canmix_denom += canmix_numer[sj][rj];
						prefer_denom += prefer_numer[sj][rj];

						mix_pop_union[si][ri][sj][rj] = supply_pop_union[sj][rj] * (si != sj);
						union_denom += mix_pop_union[si][ri][sj][rj];
					}
				}
				for (sj = DP::SEX_MIN; sj <= DP::SEX_MAX; ++sj) {
					for (rj = DP::POP_NEVER; rj < DP::N_POP_SEX[sj]; ++rj) {
						mix_pop_other[si][ri][sj][rj] = (1.0 - assort) * canmix_numer[sj][rj] / canmix_denom + assort * prefer_numer[sj][rj] / prefer_denom;
						mix_pop_union[si][ri][sj][rj] /= union_denom;
					}
				}
			}
		}

		for (sj = 0; sj < DP::N_SEX; ++sj) {
			for (bj = 0; bj < DP::N_AGE_ADULT; ++bj)
				for (rj = 0; rj < DP::N_POP_SEX[sj]; ++rj) // TODO: start from rj = DP::POP_NEVER instead of DP::POP_NOSEX
					for (hj = 0; hj < DP::N_STAGE; ++hj)
						for (vj = 0; vj < DP::N_VL; ++vj)
							prev[sj][bj][rj][hj][vj] = 0.0;
		}

		// pre-calculate PLHIV among potential sex partners
		for (uj = 0; uj < DP::N_SEX_MC; ++uj) {
			sj = sex[uj];
			for (bj = 0; bj < DP::N_AGE_ADULT; ++bj) {
				for (rj = DP::POP_NEVER; rj < DP::N_POP_SEX[sj]; ++rj) {
					for (cj = 0; cj < DP::N_HIV_ADULT; ++cj) {
						hj = stage[cj];
						num_art = pop.adult_hiv(t, uj, bj, rj, cj, DP::DTX_ART1) + pop.adult_hiv(t, uj, bj, rj, cj, DP::DTX_ART2) + pop.adult_hiv(t, uj, bj, rj, cj, DP::DTX_ART3);
						prev[sj][bj][rj][hj][DP::VL_OFF_ART] += pop.adult_hiv(t, uj, bj, rj, cj, DP::DTX_UNAWARE);
						prev[sj][bj][rj][hj][DP::VL_OFF_ART] += pop.adult_hiv(t, uj, bj, rj, cj, DP::DTX_AWARE  );
						prev[sj][bj][rj][hj][DP::VL_OFF_ART] += pop.adult_hiv(t, uj, bj, rj, cj, DP::DTX_PREV_TX);
						prev[sj][bj][rj][hj][DP::VL_FAILURE] += num_art * (1.0 - dat.art_suppressed_adult(t, sj, bj));
						prev[sj][bj][rj][hj][DP::VL_SUCCESS] += num_art * dat.art_suppressed_adult(t, sj, bj);
					}
				}
			}
		}

		// Divide by population sizes to convert PLHIV to prevalence in potential sex partners
		for (sj = 0; sj < DP::N_SEX; ++sj) {
			for (bj = 0; bj < DP::N_AGE_ADULT; ++bj)
				for (rj = DP::POP_NEVER; rj < DP::N_POP_SEX[sj]; ++rj)
					for (hj = 0; hj < DP::N_STAGE; ++hj)
						for (vj = 0; vj < DP::N_VL; ++vj)
							prev[sj][bj][rj][hj][vj] /= popsize[sj][bj][rj];
		}

		// calculate transmission probabilities
		// TODO: doi:10.1002/14651858.CD003255 estimated condoms reduced incidence 80%, so could apply
		// condom effect to incidence rate instead of per-act probabilities. Results should be approximately
		// the same, but computation may be more efficient (fewer pow calls)
		// RG 2022-01-25: I tried optimizing by unrolling sex loops so that we could easily skip
		// the female-to-female transmission = 0 case. This did not improve performance, and may
		// have actually slowed down calculation.
		for (qij = 0; qij < DP::N_BOND; ++qij) { // TODO: calculate once per year instead of once per step
			acts_with = dat.sex_acts(qij) * dat.condom_freq(t, qij);
			acts_wout = dat.sex_acts(qij) - acts_with;
			for (si = 0; si < DP::N_SEX; ++si)
				for (sj = 0; sj < DP::N_SEX; ++sj)
					for (hj = 0; hj < DP::N_STAGE; ++hj)
						for (vj = 0; vj < DP::N_VL; ++vj) {
							per_act = dat.hiv_risk_per_act(si, sj, hj, vj);
							ptransmit[si][sj][qij][hj][vj] = 1.0 - pow(1.0 - per_act, acts_wout) * pow(1.0 - per_act * e_condom, acts_with);
						}
		}

		for (si = 0; si < DP::N_SEX; ++si) {
			for (bi = 0; bi < DP::N_AGE_ADULT; ++bi)
				for (ri = DP::POP_NEVER; ri < DP::N_POP_SEX[si]; ++ri) {
					force_other[si][bi][ri] = 0.0;
					force_union[si][bi][ri] = 0.0;
				}
		}

		// The loop below is hideously expensive. Before putting ANYTHING
		// in this loop, ask yourself whether it could be precalculated
		// outside this loop. If not, put the calculation at the highest
		// level of this loop possible.
		for (si = DP::SEX_MIN; si <= DP::SEX_MAX; ++si) {
			for (bi = 0; bi < DP::N_AGE_ADULT; ++bi) {
				for (ri = DP::POP_NEVER; ri < DP::N_POP_SEX[si]; ++ri) {
					for (sj = DP::SEX_MIN; sj <= DP::SEX_MAX; ++sj) {
						// we do not model female-to-female transmission, so we only
						// execute the inner loop when at least one partner is male
						if (si == DP::MALE || sj == DP::MALE) {
							for (bj = 0; bj < DP::N_AGE_ADULT; ++bj) {
								for (rj = DP::POP_NEVER; rj < DP::N_POP_SEX[sj]; ++rj) {
									// non-marital, non-cohabiting partnerships
									bal_denom = supply_other[si][bi][ri] * dat.partner_preference_age(si, bi, sj, bj) * mix_pop_other[si][ri][sj][rj];
									bal_numer = supply_other[sj][bj][rj] * dat.partner_preference_age(sj, bj, si, bi) * mix_pop_other[sj][rj][si][ri];
									bal_other = (bal_denom > 0.0 ? sqrt(bal_numer / bal_denom) : 0.0);
									mix_other = dat.partner_preference_age(si, bi, sj, bj) * mix_pop_other[si][ri][sj][rj];
									if (mix_other > 0.0 && bal_other > 0.0 && popsize[sj][bj][rj] > 0.0) {
										qij = DP::BOND_TYPE[si][ri][sj][rj];
										bal_mix = mix_other * bal_other;

										force_group = 0.0;
										for (hj = 0; hj < DP::N_STAGE; ++hj)
											for (vj = 0; vj < DP::N_VL; ++vj)
												force_group += ptransmit[si][sj][qij][hj][vj] * prev[sj][bj][rj][hj][vj]; // TODO: calc independent of age bi. Cache?
										force_other[si][bi][ri] += bal_mix * force_group;
									}

									// marital or cohabiting partnerships
									bal_denom = supply_union[si][bi][ri] * dat.partner_preference_age(si, bi, sj, bj) * mix_pop_union[si][ri][sj][rj];
									bal_numer = supply_union[sj][bj][rj] * dat.partner_preference_age(sj, bj, si, bi) * mix_pop_union[sj][rj][si][ri];
									bal_union = (bal_denom > 0.0 ? sqrt(bal_numer / bal_denom) : 0.0);
									mix_union = dat.partner_preference_age(si, bi, sj, bj) * mix_pop_union[si][ri][sj][rj];
									if (mix_union > 0.0 && bal_union > 0.0 && popsize[sj][bj][rj] > 0.0) {
										qij = DP::BOND_UNION;
										bal_mix = mix_union * bal_union;

										force_group = 0.0;
										for (hj = 0; hj < DP::N_STAGE; ++hj)
											for (vj = 0; vj < DP::N_VL; ++vj)
												force_group += ptransmit[si][sj][qij][hj][vj] * prev[sj][bj][rj][hj][vj]; // TODO: calc independent of age bi. Cache?
										force_union[si][bi][ri] += bal_mix * force_group;
									}
								}
							}
						}
					}
					force[si][bi][ri] = dat.partner_rate(t, si, bi, ri) * force_other[si][bi][ri] + prop_union[si][ri] * force_union[si][bi][ri];
				}
			}
		}

		double force_pwid[DP::N_SEX];
		force_pwid[DP::FEMALE] = dat.pwid_needle_sharing(t) * dat.pwid_infection_force(t, DP::FEMALE);
		force_pwid[DP::MALE  ] = dat.pwid_needle_sharing(t) * dat.pwid_infection_force(t, DP::MALE);

		// update the population and record new infections
		for (ui = 0; ui < DP::N_SEX_MC; ++ui) {
			si = sex[ui];
			vmmc_mult = (ui == DP::MALE_C) ? 1.0 - dat.effect_vmmc() : 1.0;
			for (bi = 0; bi < DP::N_AGE_ADULT; ++bi) {
				ai = bi + DP::AGE_ADULT_MIN;
				for (ri = DP::POP_NEVER; ri < DP::N_POP_SEX[si]; ++ri) {
					prop_transmit = 1.0 - exp(-HIV_STEP_SIZE * (force[si][bi][ri] * vmmc_mult + force_pwid[si] * (ri == DP::POP_PWID)));
					new_hiv = prop_transmit * pop.adult_neg(t, ui, bi, ri);
					pop.adult_neg(t, ui, bi, ri) -= new_hiv;
					pop.adult_hiv(t, ui, bi, ri, DP::HIV_PRIMARY, DP::DTX_UNAWARE) += new_hiv;
					dat.new_hiv_infections(t, ui, ai, ri, dat.new_hiv_infections(t, ui, ai, ri) + new_hiv);
				}
			}
		}
	}

	void Projection::insert_adult_infections(const int t, const int step) {
		// TODO: This was originally implemented when Spectrum calculated infections
		// once per year instead of once per timestep. Calculations that refer
		// to year t-1 could be done once annually instead of once per timestep.
		const double eps = 1e-8 / 3.0; // padding term used to avoid divide-by-zero issues without resorting to conditional logic
		const double irr_sex(dat.irr_sex(t));
		int s, u, a, b, r;
		double X[DP::N_SEX];
		double neg_age[DP::N_SEX][DP::N_AGE_ADULT], new_hiv_age[DP::N_SEX][DP::N_AGE_ADULT], new_hiv_pop[DP::N_POP];
		double new_hiv_sex[DP::N_SEX], new_hiv;
		double new_hiv_all[DP::N_SEX_MC][DP::N_AGE_ADULT][DP::N_POP];
		double scale, denom, wnum;

		// Add up HIV-negative X individuals by sex in the reproductive age range
		X[DP::FEMALE] = X[DP::MALE] = 0.0;
		for (a = DP::AGE_BIRTH_MIN; a <= DP::AGE_BIRTH_MAX; ++a) {
			b = a - DP::AGE_BIRTH_MIN;
			for (r = DP::POP_MIN; r <= DP::POP_MAX; ++r) {
				X[DP::FEMALE] += pop.adult_neg(t - 1, DP::FEMALE, b, r);
				X[DP::MALE  ] += pop.adult_neg(t - 1, DP::MALE_U, b, r) + pop.adult_neg(t-1, DP::MALE_C, b, r);
			}
		}

		new_hiv = DP::HIV_STEP_SIZE * dat.incidence(t) * (X[DP::FEMALE] + X[DP::MALE]);
		new_hiv_sex[DP::MALE  ] = new_hiv / (irr_sex * X[DP::FEMALE] + X[DP::MALE]) * X[DP::MALE  ];
		new_hiv_sex[DP::FEMALE] = new_hiv / (irr_sex * X[DP::FEMALE] + X[DP::MALE]) * X[DP::FEMALE] * irr_sex;

		for (b = 0; b < DP::N_AGE_ADULT; ++b) {
			neg_age[DP::MALE][b] = neg_age[DP::FEMALE][b] = 0.0;
			for (r = DP::POP_MIN; r <= DP::POP_MAX; ++r) {
				neg_age[DP::MALE  ][b] += pop.adult_neg(t, DP::MALE_U, b, r) + pop.adult_neg(t, DP::MALE_C, b, r);
				neg_age[DP::FEMALE][b] += pop.adult_neg(t, DP::FEMALE, b, r);
			}
		}

		for (s = DP::SEX_MIN; s <= DP::SEX_MAX; ++s) {
			denom = 0.0;
			for (a = DP::AGE_BIRTH_MIN; a <= DP::AGE_BIRTH_MAX; ++a) {
				b = a - DP::AGE_BIRTH_MIN;
				denom += neg_age[s][b] * dat.irr_age(t, s, a);
			}

			scale = new_hiv_sex[s] / denom;
			for (a = DP::AGE_ADULT_MIN; a <= DP::AGE_ADULT_MAX; ++a) {
				b = a - DP::AGE_ADULT_MIN;
				new_hiv_age[s][b] = scale * neg_age[s][b] * dat.irr_age(t, s, a);
			}
		}

		// Calculate new female infections
		for (b = 0; b < DP::N_AGE_ADULT; ++b) {
			a = b + DP::AGE_ADULT_MIN;
			for (r = DP::POP_MIN; r <= DP::POP_MAX; ++r)
				new_hiv_pop[r] = dat.irr_pop(t, DP::FEMALE, r) * pop.adult_neg(t, DP::FEMALE, b, r);
			denom = std::accumulate(new_hiv_pop, new_hiv_pop + DP::N_POP, 0.0);
			scale = new_hiv_age[DP::FEMALE][b] / denom;
			for (r = DP::POP_MIN; r <= DP::POP_MAX; ++r)
				new_hiv_all[DP::FEMALE][b][r] = scale * new_hiv_pop[r];
		}

		// Calculate new male infections
		for (b = 0; b < DP::N_AGE_ADULT; ++b) {
			a = b + DP::AGE_ADULT_MIN;
			for (r = DP::POP_MIN; r <= DP::POP_MAX; ++r)
				new_hiv_pop[r] = dat.irr_pop(t, DP::MALE, r) * (pop.adult_neg(t, DP::MALE_U, b, r) + pop.adult_neg(t, DP::MALE_C, b, r));
			denom = std::accumulate(new_hiv_pop, new_hiv_pop + DP::N_POP, 0.0);
			scale = new_hiv_age[DP::MALE][b] / denom;
			for (r = DP::POP_MIN; r <= DP::POP_MAX; ++r) {
				wnum = pop.adult_neg(t, DP::MALE_U, b, r) + (1.0 - dat.effect_vmmc()) * pop.adult_neg(t, DP::MALE_C, b, r) + eps;
				new_hiv_all[DP::MALE_U][b][r] = pop.adult_neg(t, DP::MALE_U, b, r) * scale * new_hiv_pop[r] / wnum;
				new_hiv_all[DP::MALE_C][b][r] = pop.adult_neg(t, DP::MALE_C, b, r) * scale * new_hiv_pop[r] / wnum * (1.0 - dat.effect_vmmc());
			}
		}

		// Distribute new infections
#ifndef SPECTRUM_CD4
		for (u = DP::SEX_MC_MIN; u <= DP::SEX_MC_MAX; ++u) {
			s = sex[u];
			for (b = 0; b < DP::N_AGE_ADULT; ++b) {
				a = b + DP::AGE_ADULT_MIN;
				for (r = DP::POP_MIN; r <= DP::POP_MAX; ++r) {
					dat.new_hiv_infections(t, u, a, r, dat.new_hiv_infections(t, u, a, r) + new_hiv_all[u][b][r]);
					pop.adult_neg(t, u, b, r) -= new_hiv_all[u][b][r];
					pop.adult_hiv(t, u, b, r, DP::HIV_PRIMARY, DP::DTX_UNAWARE) += new_hiv_all[u][b][r];
				}
			}
		}
#else
		for (u = DP::SEX_MC_MIN; u <= DP::SEX_MC_MAX; ++u) {
			s = sex[u];
			for (b = 0; b < DP::N_AGE_ADULT; ++b) {
				a = b + DP::AGE_ADULT_MIN;
				for (r = DP::POP_MIN; r <= DP::POP_MAX; ++r) {
					dat.new_hiv_infections(t, u, a, r, dat.new_hiv_infections(t, u, a, r) + new_hiv_all[u][b][r]);
					pop.adult_neg(t, u, b, r) -= new_hiv_all[u][b][r];
					for (int h = DP::HIV_ADULT_MIN; h <= DP::HIV_ADULT_MAX; ++h)
						pop.adult_hiv(t, u, b, r, h, DP::DTX_UNAWARE) += dat.hiv_dist(s, a, h) * new_hiv_all[u][b][r];
				}
			}
		}
#endif
	}

	void Projection::insert_clhiv_agein(const int t) {
		const double eps = std::numeric_limits<double>::epsilon(); // padding to avoid divide-by-zero
		const int a = 14;
		int s, h, d;
		const double numer = pop.child_neg(t, DP::MALE_C, a);
		const double denom = pop.child_neg(t, DP::MALE_C, a) + pop.child_neg(t, DP::MALE_U, a);
		const double p_circ = numer / (denom + eps);

		for (h = DP::HIV_CHILD_MIN; h <= DP::HIV_CHILD_MAX; ++h) {
			for (d = DP::DTX_MIN; d <= DP::DTX_MAX; ++d) {
				s = DP::FEMALE;
				pop.child_hiv(t, s, a, h, d) += dat.clhiv_agein(t, s, h, d);
				pop.child_neg(t, s, a) -= dat.clhiv_agein(t, s, h, d);

				s = DP::MALE;
				pop.child_hiv(t, DP::MALE_U, a, h, d) += dat.clhiv_agein(t, s, h, d) * (1.0 - p_circ);
				pop.child_hiv(t, DP::MALE_C, a, h, d) += dat.clhiv_agein(t,s,h,d) * p_circ;
				pop.child_neg(t, DP::MALE_U, a) -= dat.clhiv_agein(t, s, h, d) * (1.0 - p_circ);
				pop.child_neg(t, DP::MALE_C, a) -= dat.clhiv_agein(t, s, h, d) * p_circ;
			}
		}
	}

	void Projection::insert_endyear_migrants(const int t) {
		double migr;
		int a, b, d, h, r, s, u;

		calc_popsize(t);

		// Children (ages 0-14)
		for (u = SEX_MC_MIN; u <= DP::SEX_MC_MAX; ++u) {
			s = sex[u];
			for (a = DP::AGE_CHILD_MIN; a <= DP::AGE_CHILD_MAX; ++a) {
				migr = dat.migration(t, s, a) / dat.popsize(t, s, a);
				pop.child_neg(t, u, a) *= (1.0 + migr);
				for (h = DP::HIV_CHILD_MIN; h <= DP::HIV_CHILD_MAX; ++h)
					for (d = DP::DTX_MIN; d <= DP::DTX_MAX; ++d)
						pop.child_hiv(t, u, a, h, d) *= (1.0 + migr);
			}
		}

		// Adults (ages 15+)
		for (u = SEX_MC_MIN; u <= DP::SEX_MC_MAX; ++u) {
			s = sex[u];
			for (a = DP::AGE_ADULT_MIN; a <= DP::AGE_ADULT_MAX; ++a) {
				b = a - DP::AGE_ADULT_MIN;
				migr = dat.migration(t, s, a) / dat.popsize(t, s, a);
				for (r = DP::POP_MIN; r <= DP::POP_MAX; ++r) {
					pop.adult_neg(t, u, b, r) *= (1.0 + migr);
					for (h = DP::HIV_ADULT_MIN; h <= DP::HIV_ADULT_MAX; ++h)
						for (d = DP::DTX_MIN; d <= DP::DTX_MAX; ++d)
							pop.adult_hiv(t, u, b, r, h, d) *= (1.0 + migr);
				}
			}
		}
	}

	void Projection::seed_epidemic(const int t, const double prev) {
		double cases;
		for (int u(0); u < DP::N_SEX_MC; ++u) {
			for (int b(0); b < DP::N_AGE_ADULT; ++b)
				for (int r(0); r < DP::N_POP; ++r) {
					cases = prev * pop.adult_neg(t, u, b, r);
					pop.adult_neg(t, u, b, r) -= cases;
					pop.adult_hiv(t, u, b, r, DP::HIV_PRIMARY, DP::DTX_UNAWARE) += cases;
				}
		}
	}

	double Projection::calc_births_hiv(const int t) {
		int u;
		double discount, births, prev, scale, afrr, dfrr, hfrr;
		popsize_t nneg[2][DP::N_AGE_BIRTH], nhiv[2][DP::N_AGE_BIRTH], nart[2][DP::N_AGE_BIRTH];
		popsize_t ncd4[2][DP::N_AGE_BIRTH][DP::N_HIV_ADULT];
		double prob_mtct[DP::N_HIV_ADULT] = {0.15, 0.15, 0.27, 0.27, 0.37, 0.37, 0.37}; // TODO: hard-coded for testing, need to move these into the input workbook
		double size_abs[DP::N_HIV_ADULT], size_wgt[DP::N_HIV_ADULT];

		// Count reproductive-age women by age, HIV, and ART status in the current and previous year
		for (int i(0); i < 2; ++i) {
			u = t - 1 + i;
			for (int b(0); b < DP::N_AGE_BIRTH; ++b) {
				nneg[i][b] = nhiv[i][b] = nart[i][b] = 0.0;
				for (int h(0); h < DP::N_HIV_ADULT; ++h)
					ncd4[i][b][h] = 0.0;

				for (int r(0); r < DP::N_POP; ++r) {
					nneg[i][b] += pop.adult_neg(u, DP::FEMALE, b, r);
					for (int h(0); h < DP::N_HIV_ADULT; ++h) {
						for (int d(DP::DTX_UNAWARE); d <= DP::DTX_ART1; ++d)
							ncd4[i][b][h] += pop.adult_hiv(u, DP::FEMALE, b, r, h, d);
						for (int d(DP::DTX_ART2); d <= DP::DTX_ART3; ++d)
							nart[i][b] += pop.adult_hiv(u, DP::FEMALE, b, r, h, d);
					}
				}

				for (int h(0); h < DP::N_HIV_ADULT; ++h)
					nhiv[i][b] += ncd4[i][b][h];
			}
		}

		births = 0.0;
		for (int b(0); b < DP::N_AGE_BIRTH; ++b) {
			int a = b + DP::AGE_BIRTH_MIN;
			prev = (nhiv[1][b] + nart[1][b]) / (nhiv[1][b] + nart[1][b] + nneg[1][b]);
			afrr = dat.frr_age_no_art(t, b);
			dfrr = dat.frr_age_on_art(b);
			discount = (nart[1][b] + nart[0][b]) * dfrr;
			for (int h(0); h < DP::N_HIV_ADULT; ++h) {
				hfrr = dat.frr_cd4_no_art(h);
				discount += (ncd4[1][b][h] + ncd4[0][b][h]) * afrr * hfrr;
			}
			discount /= (nhiv[1][b] + nart[1][b] + nhiv[0][b] + nart[1][b]);
			scale = discount / (discount * prev + 1.0 - prev);

			births += scale * 0.5 * (nhiv[1][b] + nart[1][b] + nhiv[0][b] + nart[0][b]) * dat.tfr(t) * dat.pasfrs(t, a);
		}

		double hiv_perinatal_abs(0), hiv_perinatal_wgt(0), denom_abs, denom_wgt;
		std::fill(size_abs, size_abs + DP::N_HIV_ADULT, 0.0);
		std::fill(size_wgt, size_wgt + DP::N_HIV_ADULT, 0.0);
		for (int b(0); b < DP::N_AGE_BIRTH; ++b) {
			int a = b + DP::AGE_BIRTH_MIN;
			afrr = dat.frr_age_no_art(t, b);
			for (int h(0); h < DP::N_HIV_ADULT; ++h) {
				hfrr = dat.frr_cd4_no_art(h);
				size_abs[h] += ncd4[1][b][h];
				size_wgt[h] += ncd4[1][b][h] * afrr * hfrr * dat.pasfrs(t, a);
			}
		}
		denom_abs = std::accumulate(size_abs, size_abs + DP::N_HIV_ADULT, 0.0);
		denom_wgt = std::accumulate(size_wgt, size_wgt + DP::N_HIV_ADULT, 0.0);

		for (int h(0); h < DP::N_HIV_ADULT; ++h) {
			hiv_perinatal_abs += prob_mtct[h] * size_abs[h];
			hiv_perinatal_wgt += prob_mtct[h] * size_wgt[h];
		}
		hiv_perinatal_abs *= births / denom_abs;
		hiv_perinatal_wgt *= births / denom_wgt;

		//fprintf(stderr, "%s[%d] t=%d %12.4f %12.4f %12.4f\n", __FILE__, __LINE__, t + 1970 - 1, births, hiv_perinatal_abs, hiv_perinatal_wgt);

		return births;
	}

	void Projection::calc_births(const int t) {
		const int s(DP::FEMALE);
		const double perc_m(dat.srb(t) / (dat.srb(t) + 100.0));
		const double perc_f(1.0 - perc_m);
		popsize_t female[DP::N_AGE];
		popsize_t denom(0.0), births(0.0);

		int a, b, d, h, r;
		double surv, mort;

		for (a = DP::AGE_BIRTH_MIN; a <= DP::AGE_BIRTH_MAX; ++a)
			denom += dat.pasfrs(t,a);

		for (a = DP::AGE_BIRTH_MIN; a <= DP::AGE_BIRTH_MAX; ++a) {
			b = a - DP::AGE_BIRTH_MIN;
			female[a] = 0.0;
			for (r = DP::POP_MIN; r <= DP::POP_MAX; ++r)
				female[a] += 0.5 * (pop.adult_neg(t, s, b, r) + pop.adult_neg(t - 1, s, b, r));
			for (r = DP::POP_MIN; r <= DP::POP_MAX; ++r)
				for (h = DP::HIV_ADULT_MIN; h <= DP::HIV_ADULT_MAX; ++h)
					for (d = DP::DTX_MIN; d <= DP::DTX_MAX; ++d)
						female[a] += 0.5 * (pop.adult_hiv(t, s, b, r, h, d) + pop.adult_hiv(t - 1, s, b, r, h, d));

			births += female[a] * dat.tfr(t) * dat.pasfrs(t,a) / denom;
		}

		dat.births(t, DP::MALE,   births * perc_m);
		dat.births(t, DP::FEMALE, births * perc_f);

		// Note: since DP::MALE = DP::MALE_U, all newborn males are uncircumcised
		for (int u = DP::SEX_MIN; u <= DP::SEX_MAX; ++u) {
			surv = dat.Sx(t,u,0);
			mort = 1.0 - surv;
			pop.child_neg(t, u, 0) = dat.births(t,u) * surv;
			dth.child_neg(t, u, 0) = dat.births(t,u) * mort;
		}
	}

	void Projection::calc_deaths(const int t) {
		double mort, deaths;
		int a, s;

		for (s = DP::SEX_MIN; s <= DP::SEX_MAX; ++s) {
			// age 0
			a = DP::AGE_MIN;
			mort = 1.0 - dat.Sx(t,s,a);
			deaths = dat.births(t,s) * mort;
			dat.deaths(t,s,a,deaths);

			// ages 1-79
			for (a = DP::AGE_MIN+1; a <= DP::AGE_MAX-1; ++a) {
				mort = 1.0 - dat.Sx(t,s,a);
				deaths = dat.popsize(t-1,s,a-1) * mort;
				dat.deaths(t,s,a,deaths);
			}

			// ages 80+
			a = DP::AGE_MAX;
			const double mort_80(1.0 - dat.Sx(t,s,a));
			const double mort_81(1.0 - dat.Sx(t,s,a+1));
			deaths = dat.popsize(t-1,s,a-1) * mort_80 + dat.popsize(t-1,s,a) * mort_81;
			dat.deaths(t,s,a,deaths);
		}

	}

} // END namespace DP
