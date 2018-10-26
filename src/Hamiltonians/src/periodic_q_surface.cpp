//
// Created by Panagiotis Zestanakis on 25/10/18.
//

#include <iostream>
#include "periodic_q_surface.hpp"
#include "boost/math/constants/constants.hpp"
namespace Integrators
{
    namespace Geometry
    {
        using boost::math::double_constants::pi;


        template<typename T>
        inline bool  different_sign(T d1, T d2)
        {
            return (d1 < 0 && d2 >= 0) || (d1 > 0 && d2 <= 0);
        };

        template <typename T>
        inline bool not_too_far(T d1, T d2, T threshold = pi )
        {
            return std::abs(d1 - d2) < threshold;
        };

        template<typename T>
        inline bool  cross_zero_and_not_branch_cut(T d1, T d2)
        {
          return different_sign(d1,d2)&&not_too_far(d1,d2);
        }

        bool PeriodicQSurfaceCrossObserver::operator() (const State2& next_point) const noexcept
        {
          const double next_distance = periodicQDistanceCalculator_.distance(next_point.q());
          const bool crossed_line = cross_zero_and_not_branch_cut(distance_,next_distance);
          distance_ = next_distance;
          return crossed_line;
        }

        double PeriodicQSurfaceCrossObserver::distance () const noexcept
        {
          return distance_;
        }
    }
}