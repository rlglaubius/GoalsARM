#ifndef GB_DEMOG_INTERP_H
#define GB_DEMOG_INTERP_H

namespace GB {

	const int NUM_PANEL = 5;
	const int COEFF_DIM = 5;

	const double COEFF_SPRAGUE[NUM_PANEL][COEFF_DIM][COEFF_DIM] = {
		{{ 0.3616, -0.2768,  0.1488, -0.0336,  0.0000},
		 { 0.2640, -0.0960,  0.0400, -0.0080,  0.0000},
		 { 0.1840,  0.0400, -0.0320,  0.0080,  0.0000},
		 { 0.1200,  0.1360, -0.0720,  0.0160,  0.0000},
		 { 0.0704,  0.1968, -0.0848,  0.0176,  0.0000}},

		{{ 0.0336,  0.2272, -0.0752,  0.0144,  0.0000},
		 { 0.0080,  0.2320, -0.0480,  0.0080,  0.0000},
		 {-0.0080,  0.2160, -0.0080,  0.0000,  0.0000},
		 {-0.0160,  0.1840,  0.0400, -0.0080,  0.0000},
		 {-0.0176,  0.1408,  0.0912, -0.0144,  0.0000}},

		{{-0.0128,  0.0848,  0.1504, -0.0240,  0.0016},
		 {-0.0016,  0.0144,  0.2224, -0.0416,  0.0064},
		 { 0.0064, -0.0336,  0.2544, -0.0336,  0.0064},
		 { 0.0064, -0.0416,  0.2224,  0.0144, -0.0016},
		 { 0.0016, -0.0240,  0.1504,  0.0848, -0.0128}},

		{{ 0.0000, -0.0144,  0.0912,  0.1408, -0.0176},
		 { 0.0000, -0.0080,  0.0400,  0.1840, -0.0160},
		 { 0.0000,  0.0000, -0.0080,  0.2160, -0.0080},
		 { 0.0000,  0.0080, -0.0480,  0.2320,  0.0080},
		 { 0.0000,  0.0144, -0.0752,  0.2272,  0.0336}},

		{{ 0.0000,  0.0176, -0.0848,  0.1968,  0.0704},
		 { 0.0000,  0.0160, -0.0720,  0.1360,  0.1200},
		 { 0.0000,  0.0080, -0.0320,  0.0400,  0.1840},
		 { 0.0000, -0.0080,  0.0400, -0.0960,  0.2640},
		 { 0.0000, -0.0336,  0.1488, -0.2768,  0.3616}}
	};

	const double COEFF_BEERS_ORDINARY[NUM_PANEL][COEFF_DIM][COEFF_DIM] = {
		{{ 0.3333, -0.1636, -0.0210,  0.0796, -0.0283},
		 { 0.2595, -0.0780,  0.0130,  0.0100, -0.0045},
		 { 0.1924,  0.0064,  0.0184, -0.0256,  0.0084},
		 { 0.1329,  0.0844,  0.0054, -0.0356,  0.0129},
		 { 0.0819,  0.1508, -0.0158, -0.0284,  0.0115}},
		
		{{ 0.0404,  0.2000, -0.0344, -0.0128,  0.0068},
		 { 0.0093,  0.2268, -0.0402,  0.0028,  0.0013},
		 {-0.0108,  0.2272, -0.0248,  0.0112, -0.0028},
		 {-0.0198,  0.1992,  0.0172,  0.0072, -0.0038},
		 {-0.0191,  0.1468,  0.0822, -0.0084, -0.0015}},
		
		{{-0.0117,  0.0804, 0.1570, -0.0284,  0.0027},
		 {-0.0020,  0.0160, 0.2200, -0.0400,  0.0060},
		 { 0.0050, -0.0280, 0.2460, -0.0280,  0.0050},
		 { 0.0060, -0.0400, 0.2200,  0.0160, -0.0020},
		 { 0.0027, -0.0284, 0.1570,  0.0804, -0.0117}},

		{{-0.0015, -0.0084,  0.0822, 0.1468, -0.0191},
		 {-0.0038,  0.0072,  0.0172, 0.1992, -0.0198},
		 {-0.0028,  0.0112, -0.0248, 0.2272, -0.0108},
		 { 0.0013,  0.0028, -0.0402, 0.2268,  0.0093},
		 { 0.0068, -0.0128, -0.0344, 0.2000,  0.0404}},

		{{ 0.0115, -0.0284, -0.0158,  0.1508, 0.0819},
		 { 0.0129, -0.0356,  0.0054,  0.0844, 0.1329},
		 { 0.0084, -0.0256,  0.0184,  0.0064, 0.1924},
		 {-0.0045,  0.0100,  0.0130, -0.0780, 0.2595},
		 {-0.0283,  0.0796, -0.0210, -0.1636, 0.3333}}
	};


	// Use osculatory interpolation to disaggregate values by five-year age
	// groups to single ages
	// @param x5 pointer to array of values by five-year age group
	// @param n5 length of x5
	// @param x1 return value: pointer to array of values by single year. 
	// @param n1 length of x1. n1 must be at least 5*n5
	// @param coeff Interpolation coefficients
	// @return EXIT_SUCCESS on success, EXIT_FAILURE on failure. See Details.
	// @details
	// Interpolation coefficients are available for Sprague (COEFF_SPRAGUE)
	// or Beers ordinary (COEF_BEERS_ORDINARY) methods.
	// 
	// Interpolation fails if x1 is too small to store the disaggregated values, or if n5 < 5. x1 is undefined when EXIT_FAILURE is returned.
	// 
	// Only the first n5 * 5 elements of x1 are initialized on success.
	template<typename RealType>
	int demog_interp(const RealType* const x5, const size_t n5, RealType* x1, const size_t n1, const double coeff[][COEFF_DIM][COEFF_DIM]);

	template<typename RealType>
	void demog_interp_split_group(const RealType* const x5, RealType* const x1, const double panel[][COEFF_DIM]);

	template<typename RealType>
	int demog_interp(const RealType* const x5, const size_t n5, RealType* x1, const size_t n1, const double coeff[][COEFF_DIM][COEFF_DIM]) {
		if ((5 * n5 > n1) || (n5 < NUM_PANEL))
			return EXIT_FAILURE;

		demog_interp_split_group(x5, x1, coeff[0]);
		demog_interp_split_group(x5, x1 + 5, coeff[1]);
		for (int k(2); k < n5 - 2; ++k) demog_interp_split_group(x5 + k - 2, x1 + COEFF_DIM * k, coeff[2]);
		demog_interp_split_group(x5 + n5 - 5, x1 + COEFF_DIM * (n5 - 2), coeff[3]);
		demog_interp_split_group(x5 + n5 - 5, x1 + COEFF_DIM * (n5 - 1), coeff[4]);

		return EXIT_SUCCESS;
	}

	template<typename RealType>
	void demog_interp_split_group(const RealType* const x5, RealType* const x1, const double panel[][COEFF_DIM]) {
		x1[0] = panel[0][0] * x5[0] + panel[0][1] * x5[1] + panel[0][2] * x5[2] + panel[0][3] * x5[3] + panel[0][4] * x5[4];
		x1[1] = panel[1][0] * x5[0] + panel[1][1] * x5[1] + panel[1][2] * x5[2] + panel[1][3] * x5[3] + panel[1][4] * x5[4];
		x1[2] = panel[2][0] * x5[0] + panel[2][1] * x5[1] + panel[2][2] * x5[2] + panel[2][3] * x5[3] + panel[2][4] * x5[4];
		x1[3] = panel[3][0] * x5[0] + panel[3][1] * x5[1] + panel[3][2] * x5[2] + panel[3][3] * x5[3] + panel[3][4] * x5[4];
		x1[4] = panel[4][0] * x5[0] + panel[4][1] * x5[1] + panel[4][2] * x5[2] + panel[4][3] * x5[3] + panel[4][4] * x5[4];
	}

} // end namespace GB

#endif // GB_DEMOG_INTERP_H
