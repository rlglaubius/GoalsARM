#ifndef DPPROJECTION_H
#define DPPROJECTION_H

#include <DPConst.h>
#include <DPData.h>
#include <Population.h>

namespace DP {

	class Projection {
	public:
		typedef double popsize_t; // This sets the data type used to store populations.

		Population pop;
		Population dth;
		ModelData<popsize_t> dat;

		Projection(const int year_start, const int year_final);
		~Projection();

		void initialize(const std::string &upd_filename);

		// Perform model projection from year_start to year_end (year_end must not exceed year_final)
		// This is intended (e.g.) to enable faster model fitting when fitting data run out before year_final
		void project(const int year_end);

		// Accessors
		inline const int year_first() const {return _year_first;}
		inline const int year_final() const {return _year_final;}
		inline const int num_years() const {return _num_years;}

		double calc_births(const int time);

	private:
		// Projection is not default-constructible - start and final years must be specified
		Projection();

		void calc_popsize(const int time);

		void init_baseyear_population();
		void init_baseyear_risk();
		void init_baseyear_male_circumcision();
		void calc_births_baseyear();
		void calc_deaths_baseyear();

		void project_one_year(const int time);

		void advance_one_year_demography(const int time);
		void advance_one_year_risk(const int time);
		void advance_one_year_male_circumcision(const int time);
		void advance_one_year_hiv(const int time);
		void advance_one_year_hiv_adult(const int time);
		void advance_one_year_hiv_child(const int time);
		void advance_one_step_hiv_adult(const int time, const int step);
		void calc_adult_art_uptake(const int time, const int step, sex_hiv_t& rate);
		void calc_adult_infections(const int time, const int step);
		void insert_adult_infections(const int time, const int step);
		void insert_clhiv_agein(const int time);
		void insert_endyear_migrants(const int time);
		void seed_epidemic(const int time, const double prev);

		double calc_births_hiv(const int time);

		void calc_deaths(const int time);

		int _year_first;
		int _year_final;
		int _num_years;

		// if _last_valid_time < 0,  then the entire population projection in pop is invalid
		// if _last_valid_time >= 0, then the population projection in pop is valid through _last_valid_year and invalid afterwards
		int _last_valid_time;
	};

} // END namespace DP

#endif // DPPROJECTION_H
