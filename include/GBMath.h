#ifndef GBMATH_H
#define GBMATH_H

#include <cmath>
#include <boost/math/policies/policy.hpp>
#include <boost/math/distributions/detail/common_error_handling.hpp>

namespace GB {

    // Fisk log-logistic distribution with shift parameter. This implementation is
    // patterned after boost's lognormal distribution. This is a partial implementation,
    // and should be expanded if other parts of the standard boost interface for
    // probability distributions are needed.
    template <class RealType = double, class Policy = boost::math::policies::policy< > >
    class fisk_distribution {
    public:
        typedef RealType value_type;
        typedef Policy policy_type;

        fisk_distribution(RealType shape, RealType scale, RealType shift = 0)
            : _shape(shape), _scale(scale), _shift(shift) {}

        RealType shape() const {return _shape;}
        RealType scale() const {return _scale;}
        RealType shift() const {return _shift;}
    private:
        RealType _shape;
        RealType _scale;
        RealType _shift;
    };

    typedef fisk_distribution<> fisk;

    template<class RealType, class Policy>
    inline bool check_shape(const char* function, RealType const& shape, RealType* result, const Policy& pol) {
        if ((shape <= 0) || !(boost::math::isfinite)(shape)) {
            *result = boost::math::policies::raise_domain_error<RealType>(
                function,
                "Shape parameter is %1%, but must be > 0 !", shape, pol);
            return false;
        }
        return true;
    }

    template<class RealType, class Policy>
    inline RealType pdf(const fisk_distribution<RealType, Policy>& dist, const RealType& x) {
        const RealType shape(dist.shape()), scale(dist.scale()), shift(dist.shift());

        static const char* function = "GB::pdf(const fisk_distribution<%1%>&, %1%)";
        RealType result = 0;
        if (0 == boost::math::detail::check_scale(function, dist.scale(), &result, Policy()))
            return result;
        if (0 == check_shape(function, dist.shape(), &result, Policy()))
            return result;

        if (x <= dist.shift())
            return 0.0;

        double z = (x - shift) / scale;
        double numer = (shape / scale) * std::pow(z, shape - 1);
        double sqrtd = 1.0 + std::pow(z, shape);
        return numer / (sqrtd * sqrtd);
    }

    template<class RealType, class Policy>
    inline RealType cdf(const fisk_distribution<RealType, Policy>& dist, const RealType& x) {
        static const char* function = "GB::cdf(const fisk_distribution<%1%>&, %1%)";

        RealType result = 0;
        if (0 == boost::math::detail::check_scale(function, dist.scale(), &result, Policy()))
            return result;
        if (0 == check_shape(function, dist.shape(), &result, Policy()))
            return result;

        if (x <= dist.shift())
            return 0.0;

        return(1.0 / (1.0 + std::pow((x - dist.shift()) / dist.scale(), -dist.shape())));
    }

} // end namespace GB

#endif // GBMATH_H
