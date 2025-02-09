#ifndef DPDATA_H
#define DPDATA_H

#include <boost/multi_array.hpp>
#include <string>
#include <vector>
#include <GBDemogInterp.h>
#include <DPDefs.h>
#include <DPUPDData.h>

namespace DP {

	/// Template class for populations
	/// @param popsize_t is a real-valued type (e.g. double or single) used to store population sizes
	template<typename popsize_t>
	class ModelData {
	public:
		ModelData(const int year_start, const int year_final);
		~ModelData();

		// Initialize should only initialize model inputs, not model output
		void initialize(const std::string &upd_filename);

		// +=+ Batch initialization +=+
		// Initialize PASFRs by single age from inputs by five-year age group
		// pasfrs5y is a two-dimensional array, rows=years and columns=ages
		// pasfrs5y.shape()[0] must be 7 (ages 15-19, 20-24, 25-29, 30-34, 35-39, 40-44, 45-49)
		// pasfrs5y.shape()[1] must be at least the number of years in the projection
		void init_pasfrs_from_5yr(const year_age_ref_t pasfrs5y);

		// Initialize net migration by sex from inputs aggregated by five-year age group
		// netmigr5y is a two-dimensional array, rows=years and columns=ages
		// netmigr5y.shape()[0] must be 17 (ages 0-4, 5-9, ..., 75-79, 80+)
		// netmigr5y.shape()[1] must be at least the number of years in the projection
		void init_migr_from_5yr(const int sex, const year_age_ref_t netmigr5y);

		// Initialize incidence rate ratios by age from inputs aggregated by five-year age group
		// airr5y is a two-dimensional array, rows=years and columns=ages
		// airr5y.shape()[0] must be 17 (ages 0-4, 5-9, ..., 75-79, 80+)
		// airr5y.shape()[1] must be at least the number of years in the projection
		void init_age_irr_from_5yr(const int sex, const year_age_ref_t airr5y);

		// Convenience functions to access projection the time frame
		inline const int year_first() const { return _year_first; }
		inline const int year_final() const { return _year_final; }
		inline const int num_years() const { return _num_years; }

		// +=+ Memory transfer +=+
		// Goals ARM core does not manage its own memory for larger arrays. The functions below
		// transfer memory to the ModelData object for specific uses.

		/// Transfer an array for storing output numbers of births
		/// @param ptr_births pointer to an array of births by year and sex
		void share_births(double* ptr_births);

		/// Transfer an array for storing output numbers of births to mothers living with HIV
		/// @param ptr_births_exposed pointer to an array of births by year
		void share_births_exposed(double* ptr_births_exposed);

		/// Transfer an array for storing output numbers of new HIV infections
		/// @param ptr_newhiv pointer to an array of new infections by year, sex, age, and behavioral risk group
		void share_new_infections(double* ptr_newhiv);
		
		/// Transfer an array of partner rates
		/// @param partner_rate pointer to an array of partner rates by year, sex, age, and population
		void share_partner_rate(double* ptr_partner_rate);

		/// Transfer age mixing preference coefficients
		/// @param ptr_mix pointer to an array of mixing terms by sex1, age1, sex2, age2
		void share_age_mixing(double* ptr_mix);

		/// Transfer behavioral risk assortativity parameter values
		/// @param ptr_assort pointer to an array of assortativity values by sex and behavioral risk group
		void share_pop_assortativity(double* ptr_assort);

		/// Transfer memory for storing the risk of HIV transmission from needle sharing
		/// @param ptr_pwid_infection_force pointer to the force of infection by year and sex among PWID who share needles
		/// @param ptr_needle_sharing pointer to the proportion of PWID who share needles by year
		void share_pwid_risk(double* ptr_pwid_infection_force, double* ptr_needle_sharing);

		// +=+ Accessors +=+
		// In accessors, time t=0 denotes year year_start
		inline double basepop(const int s, const int a) const {return _basepop[s][a];}

		inline double migration(const int t, const int s, const int a) {return _migration[t][s][a];}
		inline void migration(const int t, const int s, const int a, const double value) {_migration[t][s][a] = value;}
 
		inline double lx(const int t, const int s, const int a) const {return _lx[t][s][a];}
		inline double ex(const int t, const int s, const int a) const {return _ex[t][s][a];}
		inline double Sx(const int t, const int s, const int a) const {return _Sx[t][s][a];}

		inline double tfr(const int t) const {return _tfr[t];}
		inline void tfr(const int t, const double value) {_tfr[t] = value;}

		inline double srb(const int t) const {return _srb[t];}
		inline void srb(const int t, const double value) {_srb[t] = value;}

		inline double uptake_male_circumcision(const int t, const int a) const {return _uptake_male_circumcision[t][a];}
		inline void uptake_male_circumcision(const int t, const int a, const double value) {_uptake_male_circumcision[t][a] = value;}

		inline double debut_prop(const int s) const {return _debut_prop[s];}
		inline void debut_prop(const int s, const double value) {_debut_prop[s] = value;}

		inline double union_prop(const int s) const {return _union_prop[s];}
		inline void union_prop(const int s, const double value) {_union_prop[s] = value;}

		inline double split_prop() const {return _split_prop;}
		inline void split_prop(const double value) {_split_prop = value;}

		inline double prop_debut_in_union(const int s) const {return _prop_debut_in_union[s];}
		inline void prop_debut_in_union(const int s, const double value) {_prop_debut_in_union[s] = value;}

		// r must be DP::POP_MSM, DP::POP_TGW, DP::POP_FSW, DP::POP_KEY, or DP::POP_PWID
		inline double keypop_exit_prop(const int s, const int r) const {return _keypop_exit_prop[s][r - DP::POP_KEY_MIN];}
		inline void keypop_exit_prop(const int s, const int r, const double value) {_keypop_exit_prop[s][r - DP::POP_KEY_MIN] = value;}

		// r must be DP::POP_MSM, DP::POP_TGW, DP::POP_FSW, DP::POP_KEY, or DP::POP_PWID
		inline double keypop_size(const int s, const int r) const {return _keypop_size[s][r - DP::POP_KEY_MIN];}
		inline void keypop_size(const int s, const int r, const double value) {_keypop_size[s][r - DP::POP_KEY_MIN] = value;}

		// r must be DP::POP_MSM, DP::POP_TGW, DP::POP_FSW, DP::POP_KEY, or DP::POP_PWID
		inline bool keypop_stay(const int s, const int r) const {return _keypop_stay[s][r - DP::POP_KEY_MIN];}
		inline void keypop_stay(const int s, const int r, const bool value) {_keypop_stay[s][r - DP::POP_KEY_MIN] = value;}

		// r must be DP::POP_MSM, DP::POP_TGW, DP::POP_FSW, DP::POP_KEY, or DP::POP_PWID
		inline double keypop_age_dist(const int s, const int a, const int r) const {return _keypop_age_dist[s][a][r - DP::POP_KEY_MIN];}
		inline void keypop_age_dist(const int s, const int a, const int r, const double value) {_keypop_age_dist[s][a][r - DP::POP_KEY_MIN] = value;}

		// r must be DP::POP_MSM, DP::POP_TGW, DP::POP_FSW, DP::POP_KEY, or DP::POP_PWID
		inline double keypop_married(const int s, const int r) const {return _keypop_married[s][r - DP::POP_KEY_MIN];}
		inline void keypop_married(const int s, const int r, const double value) {_keypop_married[s][r - DP::POP_KEY_MIN] = value;}

		// Access using age 15 <= a < 50
		inline double pasfrs(const int t, const int a) const {return _pasfrs[t][a - DP::AGE_BIRTH_MIN];}
		inline void pasfrs(const int t, const int a, const double value) {_pasfrs[t][a - DP::AGE_BIRTH_MIN] = value;}

		inline double mtct_rate(const int mtct_timing, const int mtct_regimen, const int mtct_cd4) const {return _mtct_rate[mtct_timing][mtct_regimen][mtct_cd4];}
		inline void mtct_rate(const int mtct_timing, const int mtct_regimen, const int mtct_cd4, const double value) {_mtct_rate[mtct_timing][mtct_regimen][mtct_cd4] = value;}

		// age should be of type mtct_t. We subtract 1 so that MTCT_MOS_00_02 is index 0
		inline double breastfeeding(const int t, const int arv, const int age) const {return _breastfeeding[t][arv][age-1];}
		inline void breastfeeding(const int t, const int arv, const int age, const double value) {_breastfeeding[t][arv][age-1] = value;}

		inline double births(const int t, const int s) const {return (*_births)[t][s];}
		inline void births(const int t, const int s, const popsize_t value) {(*_births)[t][s] = value;}

		inline double births_hiv_exposed(const int t) const {return (*_births_exposed)[t];}
		inline void births_hiv_exposed(const int t, const popsize_t value) {(*_births_exposed)[t] = value;}

		inline double deaths(const int t, const int s, const int a) const {return _deaths[t][s][a];}
		inline void deaths(const int t, const int s, const int a, const popsize_t value) {_deaths[t][s][a] = value;}
		inline const year_sex_age_t& deaths() const {return _deaths;}

		inline double popsize(const int t, const int s, const int a) const {return _popsize[t][s][a];}
		inline void popsize(const int t, const int s, const int a, const popsize_t value) {_popsize[t][s][a] = value;}

		inline double new_hiv_infections(const int t, const int s, const int a, const int r) const {return (*_new_hiv_infections)[t][s][a][r];}
		inline void new_hiv_infections(const int t, const int s, const int a, const int r, const double value) {(*_new_hiv_infections)[t][s][a][r] = value;}

		inline bool direct_incidence() const {return _direct_incidence;}
		inline void direct_incidence(const bool value) {_direct_incidence = value;}

		inline double incidence(const int t) const {return _incidence[t];}
		inline void incidence(const int t, const double value) {_incidence[t] = value;}

		inline double irr_sex(const int t) const {return _irr_sex[t];}
		inline void irr_sex(const int t, const double value) {_irr_sex[t] = value;}

		inline double irr_age(const int t, const int s, const int a) const {return _irr_age[t][s][a];}
		inline void irr_age(const int t, const int s, const int a, const double value) {_irr_age[t][s][a] = value;}

		inline double irr_pop(const int t, const int s, const int r) const {return _irr_pop[t][s][r];}
		inline void irr_pop(const int t, const int s, const int r, const double value) {_irr_pop[t][s][r] = value;}

		// This should be relative to the epidemic start year (e.g., seed_time=5 for a projection
		// that starts in 1970 and epidemic that starts in 1975)
		inline int seed_time() const {return _seed_time;}
		inline void seed_time(const int time) {_seed_time = time;}

		inline double seed_prevalence() const {return _seed_prev;}
		inline void seed_prevalence(const double prev) {_seed_prev = prev;}

		inline double partner_rate(const int t, const int s, const int a, const int r) const {return (*_partner_rate)[t][s][a][r];}
		inline void partner_rate(const int t, const int s, const int a, const int r, const double value) {(*_partner_rate)[t][s][a][r] = value;}

		inline double partner_preference_age(const int s1, const int a1, const int s2, const int a2) const {return (*_partner_preference_age)[s1][a1][s2][a2];}
		inline void partner_preference_age(const int s1, const int a1, const int s2, const int a2, const double value) {(*_partner_preference_age)[s1][a1][s2][a2] = value;}

		inline double partner_assortativity(const int s, const int r) const {return (*_partner_assortativity)[s][r];}
		inline void partner_assortativity(const int s, const int r, const double value) {(*_partner_assortativity)[s][r] = value;}

		inline int mix_structure(const int s1, const int r1, const int s2, const int r2) const {return _mix_structure[s1][r1][s2][r2];}
		inline void mix_structure(const int s1, const int r1, const int s2, const int r2, const int value) {_mix_structure[s1][r1][s2][r2] = value;}

		inline double sex_acts(const int bond) const {return _sex_acts[bond];}
		inline void sex_acts(const int bond, const double value) {_sex_acts[bond] = value;}

		inline double condom_freq(const int t, const int bond) const {return _condom_freq[t][bond];}
		inline void condom_freq(const int t, const int bond, const double value) {_condom_freq[t][bond] = value;}

		inline double sti_prev(const int t, const int s, const int a, const int r) const {return _sti_prev[t][s][a][r];}
		inline void sti_prev(const int t, const int s, const int a, const int r, const double value) {_sti_prev[t][s][a][r] = value;}

		inline double pwid_infection_force(const int t, const int s) const {return (*_pwid_infection_force)[t][s];}
		inline void pwid_infection_force(const int t, const int s, const double value) {(*_pwid_infection_force)[t][s] = value;}

		inline double pwid_needle_sharing(const int t) const {return (*_pwid_needle_sharing)[t];}
		inline void pwid_needle_sharing(const int t, const double value) {(*_pwid_needle_sharing)[t] = value;}

		inline double hiv_dist(const int s, const int a, const int h) const {return _hiv_dist[s][a][h];}
		inline void hiv_dist(const int s, const int a, const int h, const double value) {_hiv_dist[s][a][h] = value;}

		inline double hiv_prog(const int s, const int a, const int h) const {return _hiv_prog[s][a][h];}
		inline void hiv_prog(const int s, const int a, const int h, const double value) {_hiv_prog[s][a][h] = value;}

		inline double hiv_mort(const int s, const int a, const int h) const {return _hiv_mort[s][a][h];}
		inline void hiv_mort(const int s, const int a, const int h, const double value) {_hiv_mort[s][a][h] = value;}

		// transmission risk per sex act, indexed by HIV- partner sex s_neg and HIV+ partner sex s_pos, HIV stage h, and viral load status v
		// h is stage_t (primary/chronic/symptomatic stages), not hiv_t (CD4 stages)
		inline double hiv_risk_per_act(const int s_neg, const int s_pos, const int h, const int v) const {return _hiv_transmit[s_neg][s_pos][h][v];}
		inline void hiv_risk_per_act(const int s_neg, const int s_pos, const int h, const int v, const double value) {_hiv_transmit[s_neg][s_pos][h][v] = value;}

		inline double art_mort_adult(const int t, const int s, const int a, const int h, const int d) const {return _art_mort_adult[t][s][a][h][d];}
		inline void art_mort_adult(const int t, const int s, const int a, const int h, const int d, const double value)  {_art_mort_adult[t][s][a][h][d] = value;}

		inline double art_num_adult(const int t, const int s) const {return _art_num_adult[t][s];}
		inline void art_num_adult(const int t, const int s, const double value) {_art_num_adult[t][s] = value;}

		inline double art_prop_adult(const int t, const int s) const {return _art_prop_adult[t][s];}
		inline void art_prop_adult(const int t, const int s, const double value) {_art_prop_adult[t][s] = value;}

		inline double art_exit_adult(const int t, const int s) const {return _art_exit_adult[t][s];}
		inline void art_exit_adult(const int t, const int s, const double value) {_art_exit_adult[t][s] = value;}

		inline double art_suppressed_adult(const int t, const int s, const int a) const {return _art_suppressed_adult[t][s][a];}
		inline void art_suppressed_adult(const int t, const int s, const int a, const double value) {_art_suppressed_adult[t][s][a] = value;}

		inline int art_first_eligible_stage_adult(const int t) const {return _art_first_eligible_stage_adult[t];}
		inline void art_first_eligible_stage_adult(const int t, const int h) {_art_first_eligible_stage_adult[t] = h;}

		inline double art_mort_weight() const {return _art_mort_weight;}
		inline void art_mort_weight(const double value) {_art_mort_weight = value;}

		inline double art_flow(const int d) const {return _art_flow[d-DP::DTX_ART_MIN];}
		inline void art_flow(const int d, const double value) {_art_flow[d-DP::DTX_ART_MIN] = value;}

		inline double frr_age_no_art(const int t, const int a) const {return _frr_age_no_art[t][a];}
		inline void frr_age_no_art(const int t, const int a, const double value) {_frr_age_no_art[t][a] = value;}

		inline double frr_age_on_art(const int a) const {return _frr_age_on_art[a];}
		inline void frr_age_on_art(const int a, const double value) {_frr_age_on_art[a] = value;}

		inline double frr_cd4_no_art(const int h) const {return _frr_cd4_no_art[h];}
		inline void frr_cd4_no_art(const int h, const double value) {_frr_cd4_no_art[h] = value;}

		// Access the number of HIV+ pregnant women receiving PMTCT. Rgimen should be one of MTCT_RX_SDNVP..MTCT_RX_ART_LATE
		inline double pmtct_num(const int t, const int regimen) const {return _pmtct_num[t][regimen];}
		inline void pmtct_num(const int t, const int regimen, const double value) {_pmtct_num[t][regimen] = value;}

		// Access the proportion of HIV+ pregnant women receiving PMTCT. Regimen should be one of MTCT_RX_SDNVP..MTCT_RX_ART_LATE
		inline double pmtct_prop(const int t, const int regimen) const {return _pmtct_prop[t][regimen];}
		inline void pmtct_prop(const int t, const int regimen, const double value) {_pmtct_prop[t][regimen] = value;}

		// Access the proportion of HIV+ women retained on ART at delivery among those who started ART before their current pregnancy
		inline double pmtct_retained_art_before(const int t) const {return _pmtct_retained_art_before[t];}
		inline void pmtct_retained_art_before(const int t, const double value) {_pmtct_retained_art_before[t] = value;}

		// Access the proportion of HIV+ women retained on ART at delivery among those who started ART during their current pregnancy
		inline double pmtct_retained_art_during(const int t) const {return _pmtct_retained_art_during[t];}
		inline void pmtct_retained_art_during(const int t, const double value) {_pmtct_retained_art_during[t] = value;}

		// Access the proportion of HIV+ women retained on PMTCT regimens from month-to-month during breastfeeding
		inline double pmtct_retained_postnatal(const int t, const int regimen, const int month) const {return _pmtct_retained_postnatal[t][regimen][month];}
		inline void pmtct_retained_postnatal(const int t, const int regimen, const int month, const double value) {_pmtct_retained_postnatal[t][regimen][month] = value;}

		inline double clhiv_agein(const int t, const int s, const int h, const int d) const {return _clhiv_agein[t][s][h][d];}
		inline void clhiv_agein(const int t, const int s, const int h, const int d, const double value) {_clhiv_agein[t][s][h][d] = value;}

		inline double effect_sti_hivpos() const {return _effect_sti_hivpos;}
		inline void effect_sti_hivpos(const double value) {_effect_sti_hivpos = value;}

		inline double effect_sti_hivneg() const {return _effect_sti_hivneg;}
		inline void effect_sti_hivneg(const double value) {_effect_sti_hivneg = value;}

		inline double effect_vmmc() const {return _effect_vmmc;}
		inline void effect_vmmc(const double value) {_effect_vmmc = value;}

		inline double effect_condom() const {return _effect_condom;}
		inline void effect_condom(const double value) {_effect_condom = value;}

	private:
		// Model inputs
		int _year_first;
		int _year_final;
		int _num_years;

		sex_age_t      _basepop;
		year_sex_age_t _lx;
		year_sex_age_t _ex;
		year_sex_age_t _Sx;
		time_series_t  _tfr;
		time_series_t  _srb;
		year_age_t     _pasfrs;
		year_sex_age_t _migration;

		year_age_t     _uptake_male_circumcision;

		// Model inputs - behavioral risk group sizes and dynamics
		double _debut_prop[DP::N_SEX];          // proportion who debut sexually per year
		double _union_prop[DP::N_SEX];          // proportion who marry or start cohabitating per year
		double _split_prop;                     // proportion of marriages/cohabitating partnerships that end per year
		double _prop_debut_in_union[DP::N_SEX]; // proportion who first debut in a marital or cohabiting partnership
		double _keypop_exit_prop[DP::N_SEX][DP::N_POP_KEY];                 // proportion of key pops lost to turnover each year
		double _keypop_size[DP::N_SEX][DP::N_POP_KEY];                      // proportion of 15-49 sex s who are in population r
		bool   _keypop_stay[DP::N_SEX][DP::N_POP_KEY];                      // indicates whether key population membership is lifelong (true) or not (false)
		double _keypop_age_dist[DP::N_SEX][DP::N_AGE_ADULT][DP::N_POP_KEY]; // proportion of sex s & pop r who are age a
		double _keypop_married[DP::N_SEX][DP::N_POP_KEY];                   // proportion of key populations who are married

		// Model inputs - MTCT parameters
		double _mtct_rate[DP::N_MTCT][DP::N_MTCT_RX][DP::N_MTCT_CD4]; // MTCT rates by perinatal/breastfeeding timing, PMTCT regimen, and maternal CD4 count
		array3d_t _breastfeeding; // breastfeeding inputs by year, ARV status, and time since delivery

		bool _direct_incidence; // toggle for direct vs. mechanistic HIV incidence calculation

		// Model inputs - direct incidence. These values are undefined if direct incidence input *is not* used
		time_series_t  _incidence; // direct HIV incidence input
		time_series_t  _irr_sex;   // incidence rate ratios by sex, female-to-male
		year_sex_age_t _irr_age;   // incidence rate ratios by age
		year_sex_pop_t _irr_pop;   // incidence rate ratios by behavioral risk group

		// Model inputs - epidemic start date (ignored if direct incidence is enabled)
		int _seed_time;    // epidemic start date (number of years since _year_first)
		double _seed_prev; // initial HIV prevalence

		// Model inputs - sexual behavior. These values are undefined if direct incidence *is* used
		year_sex_age_pop_ref_t* _partner_rate;          // partner change rates by year, sex, age, and behavioral risk group
		sex_pop_ref_t*          _partner_assortativity; // assortativity weights by sex and behavioral risk. High weights mean more partnerships come from preferred groups.

		// Model inputs - mixing preferences by age
		// _partner_preference_age[s][a][z][b] is the proportion of partnerships (s,a) prefers from (z,b)
		array4d_ref_t* _partner_preference_age;

		// Model inputs - mixing structure by behavioral risk groups
		// _mix_structure[s1][r1][s2][r2] = value means:
		// value=0 : group (s1,r1) does not get partners from (s2,r2)
		// value=1 : group (s1,r1) may get partners from (s2,r2)
		// value=2 : group (s1,r1) prefers partners from (s2,r2)
		int _mix_structure[DP::N_SEX][DP::N_POP][DP::N_SEX][DP::N_POP];

		// Model inputs - sexual behavior within partnerships
		double _sex_acts[DP::N_BOND]; // sex acts per partner per year
		year_bond_t _condom_freq;     // condom use at last sex

		// Model inputs - STI symptom prevalence
		year_sex_age_pop_t _sti_prev;
		double _effect_sti_hivpos; // effect of STIs on HIV transmission if HIV+ partner STI symptomatic
		double _effect_sti_hivneg; // effect of STIs on HIV transmission if HIV- partner STI symptomatic

		// Model inputs - HIV risk from unsafe injecting practices
		year_sex_ref_t*    _pwid_infection_force; // Force of infection acting on PWID who share needles
		time_series_ref_t* _pwid_needle_sharing;  // Proportion of PWID who share needles

		// Model inputs - HIV natural history
		sex_age_hiv_t  _hiv_dist;  // distribution of disease stages at infection
		sex_age_hiv_t  _hiv_prog;  // disease progression rates with untreated HIV
		sex_age_hiv_t  _hiv_mort;  // mortality from untreated HIV

		// Model inputs - HIV transmission probability per sex act
		double _hiv_transmit[DP::N_SEX][DP::N_SEX][DP::N_STAGE][DP::N_VL];

		// Model inputs - adult ART
		year_sex_age_hiv_dtx_t _art_mort_adult;  // mortality rates on ART
		year_sex_t             _art_num_adult;   // numbers on ART
		year_sex_t             _art_prop_adult;  // percentage on ART
		year_sex_t             _art_exit_adult;  // ART interruption rates
		year_sex_age_t         _art_suppressed_adult; // proportion on ART who are virally suppressed
		double                 _art_flow[DP::N_ART]; // flow rates between ART duration categories
		double                 _art_mort_weight; // weight placed on expected mortality when allocating ART

		time_series_int_t      _art_first_eligible_stage_adult; // index of the earliest HIV stage that is eligible for ART by CD4 count threshold

		// Model inputs - HIV-related fertility effects
		year_age_t _frr_age_no_art;                  // fertility rate ratios for women by age when not on ART
		double     _frr_age_on_art[DP::N_AGE_BIRTH]; // fertility rate ratios for women by age when on ART
		double     _frr_cd4_no_art[DP::N_HIV_ADULT]; // fertility rate ratios for women by CD4 count when not on ART

		// Model inputs - PMTCT coverage and retention
		array2d_t _pmtct_num;                      // Number of HIV+ pregnant women receiving PMTCT, by regimen
		array2d_t _pmtct_prop;                     // Number of HIV+ pregnant women receiving PMTCT, by regimen
		time_series_t _pmtct_retained_art_before;  // Proportion of HIV+ pregnant women retained on ART at delivery among those who initiated before current pregnancy
		time_series_t _pmtct_retained_art_during;  // Proportion of HIV+ pregnant women retained on ART at delivery among those who initiated during current pregnancy
		array3d_t _pmtct_retained_postnatal;       // Proportion of HIV+ breastfeeding women retained on PMTCT from month-to-month, by year, regimen, and time since delivery (MTCT_MOS_xx_yy)

		// Model inputs - direct input of 14-year-old CLHIV. Undefined if direct CLHIV input is not used
		year_sex_hiv_dtx_t _clhiv_agein;

		// Model inputs - intervention effects (ART prevention effects are with transmission parameters)
		double _effect_vmmc;
		double _effect_condom;

		// Model output
		year_sex_ref_t*    _births;
		year_sex_age_t     _deaths;
		year_sex_age_t     _popsize;

		time_series_ref_t* _births_exposed;

		year_sex_age_pop_ref_t* _new_hiv_infections;
	};

	template<typename popsize_t>
	ModelData<popsize_t>::ModelData(const int year_start, const int year_final)
		: _year_first(year_start),
			_year_final(year_final),
			_num_years(year_final - year_start + 1),

			_basepop(sex_age_t(boost::extents[DP::N_SEX][DP::N_AGE])),
			_lx(year_sex_age_t(boost::extents[year_final - year_start + 1][DP::N_SEX][DP::N_AGE])),
			_ex(year_sex_age_t(boost::extents[year_final - year_start + 1][DP::N_SEX][DP::N_AGE])),
			_Sx(year_sex_age_t(boost::extents[year_final - year_start + 1][DP::N_SEX][DP::N_AGE + 1])),
			_tfr(time_series_t(year_final - year_start + 1)),
			_srb(time_series_t(year_final - year_start + 1)),
			_pasfrs(year_age_t(boost::extents[year_final - year_start + 1][DP::N_AGE_BIRTH])),
			_migration(year_sex_age_t(boost::extents[year_final - year_start + 1][DP::N_SEX][DP::N_AGE])),

			_uptake_male_circumcision(year_age_t(boost::extents[year_final - year_start + 1][DP::N_AGE])),

			_breastfeeding(array3d_t(boost::extents[year_final - year_start + 1][DP::N_BF_ARV][DP::N_MTCT_MOS-1])),
	
			_incidence(time_series_t(year_final - year_start + 1)),
			_irr_sex(time_series_t(year_final - year_start + 1)),
			_irr_age(year_sex_age_t(boost::extents[year_final - year_start + 1][DP::N_SEX][DP::N_AGE])),
			_irr_pop(year_sex_pop_t(boost::extents[year_final - year_start + 1][DP::N_SEX][DP::N_POP])),

			_partner_rate(NULL),
			_partner_preference_age(NULL),
			_partner_assortativity(NULL),

			_condom_freq(year_bond_t(boost::extents[year_final - year_start + 1][DP::N_BOND])),

			_sti_prev(year_sex_age_pop_t(boost::extents[year_final - year_start + 1][DP::N_SEX][DP::N_AGE_ADULT][DP::N_POP])),

			_hiv_dist(sex_age_hiv_t(boost::extents[DP::N_SEX][DP::N_AGE][DP::N_HIV])),
			_hiv_prog(sex_age_hiv_t(boost::extents[DP::N_SEX][DP::N_AGE][DP::N_HIV])),
			_hiv_mort(sex_age_hiv_t(boost::extents[DP::N_SEX][DP::N_AGE][DP::N_HIV])),

			_art_mort_adult(year_sex_age_hiv_dtx_t(boost::extents[year_final - year_start + 1][DP::N_SEX][DP::N_AGE_ADULT][DP::N_HIV][DP::N_DTX])),
			_art_num_adult( year_sex_t(boost::extents[year_final - year_start + 1][DP::N_SEX])),
			_art_prop_adult(year_sex_t(boost::extents[year_final - year_start + 1][DP::N_SEX])),
			_art_exit_adult(year_sex_t(boost::extents[year_final - year_start + 1][DP::N_SEX])),
			_art_suppressed_adult(year_sex_age_t(boost::extents[year_final - year_start + 1][DP::N_SEX][DP::N_AGE_ADULT])),
			_art_first_eligible_stage_adult(time_series_int_t(year_final - year_start + 1)),

			_frr_age_no_art(year_age_t(boost::extents[year_final - year_start + 1][DP::N_AGE_BIRTH])),

			_pmtct_num(array2d_t(boost::extents[year_final - year_start + 1][DP::N_MTCT_ARV_RX])),
			_pmtct_prop(array2d_t(boost::extents[year_final - year_start + 1][DP::N_MTCT_ARV_RX])),
			_pmtct_retained_art_before(time_series_t(year_final - year_start + 1)),
			_pmtct_retained_art_during(time_series_t(year_final - year_start + 1)),
			_pmtct_retained_postnatal(array3d_t(boost::extents[year_final - year_start + 1][DP::N_MTCT_ARV_RX][DP::N_MTCT_MOS])),

			_clhiv_agein(year_sex_hiv_dtx_t(boost::extents[year_final - year_start + 1][DP::N_SEX][DP::N_HIV][DP::N_DTX])),

			_births(NULL),
			_births_exposed(NULL),
			_deaths(year_sex_age_t(boost::extents[year_final - year_start + 1][DP::N_SEX][DP::N_AGE])),
			_popsize(year_sex_age_t(boost::extents[year_final - year_start + 1][DP::N_SEX][DP::N_AGE])),
			_new_hiv_infections(NULL)
	{
		art_flow(DP::DTX_ART1, 2.0); // 6 months in 1st ART state [0,6)  months
		art_flow(DP::DTX_ART2, 2.0); // 6 months in 2nd ART state [6,12) months
		art_flow(DP::DTX_ART3, 0.0); // absorbing state [12,\infty) months
	}

	template<typename popsize_t>
	ModelData<popsize_t>::~ModelData() {
		if (_births) {delete _births;}
		if (_births_exposed) {delete _births_exposed;}
		if (_partner_rate) {delete _partner_rate;}
		if (_partner_preference_age) {delete _partner_preference_age;}
		if (_partner_assortativity) {delete _partner_assortativity;}
		if (_new_hiv_infections) {delete _new_hiv_infections;}
	}

	template<typename popsize_t>
	void ModelData<popsize_t>::initialize(const std::string &upd_filename) {
		const int year_upd[4] = {1970, 1975, 1980, 1985}; // years when basepop is defined in UPD files
		int time_upd, time_dat, s, a, i;
		double wgt1, wgt2;
		UPDData upd;

		upd.read(upd_filename);

		// Calculate baseyear population by linearly interpolating between surrounding upd base populations
		if (_year_first >= year_upd[0] && _year_first < year_upd[1])
			i = 0;
		else if (_year_first >= year_upd[1] && _year_first < year_upd[2])
			i = 1;
		else
			i = 2;

		wgt2 = (_year_first - year_upd[i]) / static_cast<double>(year_upd[i+1] - year_upd[i]);
		wgt1 = 1.0 - wgt2;

		for (s = DP::SEX_MIN; s <= DP::SEX_MAX; ++s) {
			for (a = DP::AGE_MIN; a <= DP::AGE_MAX; ++a) {
				_basepop[s][a] = wgt1 * upd.basepop(i,s,a) + wgt2 * upd.basepop(i+1,s,a);
			}
		}

		for (time_dat = 0; time_dat < _num_years; ++time_dat) {
			time_upd = time_dat + _year_first - UPDData::UPD_YEAR_START;

			_tfr[time_dat] = upd.tfr(time_upd);
			_srb[time_dat] = upd.srb(time_upd);

			for (a = DP::AGE_BIRTH_MIN; a <= DP::AGE_BIRTH_MAX; ++a) {
				_pasfrs[time_dat][a - DP::AGE_BIRTH_MIN] = upd.pasfrs(time_upd,a);
			}

			for (s = DP::SEX_MIN; s <= DP::SEX_MAX; ++s) {
				for (a = DP::AGE_MIN; a <= DP::AGE_MAX; ++a) {
					_lx[time_dat][s][a] = upd.lx(time_upd,s,a);
					_ex[time_dat][s][a] = upd.ex(time_upd,s,a);
					

					_migration[time_dat][s][a] = upd.migration(time_upd,s,a);
				}
			}

			for (s = DP::SEX_MIN; s <= DP::SEX_MAX; ++s) {
				for (a = DP::AGE_MIN; a <= DP::AGE_MAX + 1; ++a)
					_Sx[time_dat][s][a] = upd.Sx(time_upd,s,a);
			}
		}
	}

	template<typename popsize_t>
	void ModelData<popsize_t>::init_pasfrs_from_5yr(const year_age_ref_t pasfrs5y) {
		size_t a5y;
		for (size_t t(0); t < _num_years; ++t)
			for (size_t a1y(DP::AGE_BIRTH_MIN); a1y <= DP::AGE_BIRTH_MAX; ++a1y) {
				a5y = (a1y - DP::AGE_BIRTH_MIN) / 5;
				pasfrs(t, a1y, 0.2 * pasfrs5y[t][a5y]);
			}
	}

	template<typename popsize_t>
	void ModelData<popsize_t>::init_migr_from_5yr(const int sex, const year_age_ref_t netmigr5y) {
		const size_t n_group(17);
		double csum(0.0);

		// Calculate first panel coefficients based on survival rates as in AIM
		double coeff[5][5][5];
		for (int r(0); r < 5; ++r)
			for (int c(0); c < 5; ++c)
				coeff[0][r][c] = 0.0;
		for (int p(1); p < 5; ++p) {
			for (int r(0); r < 5; ++r)
				for (int c(0); c < 5; ++c)
					coeff[p][r][c] = GB::COEFF_BEERS_ORDINARY[p][r][c];
		}
		coeff[0][0][0] = Sx(0, sex, 0) * 2000;
		for (int r(1); r < 5; ++r) coeff[0][r][0] = coeff[0][r-1][0] * Sx(0, sex, r);
		for (int r(0); r < 5; ++r) csum += coeff[0][r][0];
		for (int r(0); r < 5; ++r) coeff[0][r][0] /= csum;

		// TODO: reduce the amount of copying done
		double buff5y[n_group - 1];
		double buff1y[DP::N_AGE];
		for (size_t t(0); t < _num_years; ++t) {
			buff1y[DP::N_AGE - 1] = netmigr5y[t][n_group - 1]; // Copy 80+ net migrants directly
			for (size_t a(0); a < n_group - 1; ++a) buff5y[a] = netmigr5y[t][a];
			GB::demog_interp(buff5y, n_group-1, buff1y, DP::N_AGE-1, coeff);
			for (size_t a(0); a < DP::N_AGE; ++a) migration(t, sex, a, buff1y[a]);
		}
	}

	template<typename popsize_t>
	void ModelData<popsize_t>::init_age_irr_from_5yr(const int sex, const year_age_ref_t airr5y) {
		const size_t n_group(17);

		// TODO: reduce the amount of copying done
		double buff5y[n_group - 1];
		double buff1y[DP::N_AGE];
		for (size_t t(0); t < _num_years; ++t) {
			buff1y[DP::N_AGE - 1] = airr5y[t][n_group - 1]; // Copy 80+ IRR directly
			for (size_t a(0); a < n_group - 1; ++a) buff5y[a] = airr5y[t][a];
			GB::demog_interp(buff5y, n_group-1, buff1y, DP::N_AGE-1, GB::COEFF_BEERS_ORDINARY);
			for (size_t a(0); a < DP::N_AGE; ++a) irr_age(t, sex, a, std::max(buff1y[a], 0.0));
		}
	}

	template<typename popsize_t>
	void ModelData<popsize_t>::share_births(double* ptr_births) {
		_births = new year_sex_ref_t(ptr_births, boost::extents[year_final() - year_first() + 1][DP::N_SEX]);
	}

	template<typename popsize_t>
	void ModelData<popsize_t>::share_births_exposed(double* ptr_births_exposed) {
		_births_exposed = new time_series_ref_t(ptr_births_exposed, boost::extents[year_final() - year_first() + 1]);
	}

	template<typename popsize_t>
	void ModelData<popsize_t>::share_new_infections(double* ptr_newhiv) {
		_new_hiv_infections = new year_sex_age_pop_ref_t(ptr_newhiv, boost::extents[year_final() - year_first() + 1][DP::N_SEX_MC][DP::N_AGE][DP::N_POP]);
	}

	template<typename popsize_t>
	void ModelData<popsize_t>::share_partner_rate(double* ptr_partner_rate) {
		_partner_rate = new year_sex_age_pop_ref_t(ptr_partner_rate, boost::extents[year_final() - year_first() + 1][DP::N_SEX][DP::N_AGE_ADULT][DP::N_POP]);
	}

	template<typename popsize_t>
	void ModelData<popsize_t>::share_age_mixing(double* ptr_mix) {
		_partner_preference_age = new array4d_ref_t(ptr_mix, boost::extents[DP::N_SEX][DP::N_AGE_ADULT][DP::N_SEX][DP::N_AGE_ADULT]);
	}

	template<typename popsize_t>
	void ModelData<popsize_t>::share_pop_assortativity(double* ptr_assort) {
		_partner_assortativity = new sex_pop_ref_t(ptr_assort, boost::extents[DP::N_SEX][DP::N_POP]);
	}

	template<typename popsize_t>
	void ModelData<popsize_t>::share_pwid_risk(double* ptr_pwid_infection_force, double* ptr_needle_sharing) {
		_pwid_infection_force = new year_sex_ref_t(ptr_pwid_infection_force, boost::extents[year_final() - year_first() + 1][DP::N_SEX]);
		_pwid_needle_sharing  = new time_series_ref_t(ptr_needle_sharing, boost::extents[year_final() - year_first() + 1]);
	}
}

#endif // DPDATA_H