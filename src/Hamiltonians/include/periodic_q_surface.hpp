//
// Created by Panagiotis Zestanakis on 25/10/18.
//

#ifndef HAMILTONIANS_PERIODIC_Q_SURFACE_HPP
#define HAMILTONIANS_PERIODIC_Q_SURFACE_HPP

#include "details/periodic_q_distance.hpp"
#include "State.hpp"

namespace Integrators
{
    namespace Geometry
    {

        class PeriodicQSurfaceCrossObserver
        {
         private:
          Integrators::Internals::PeriodicQDistance periodicQDistanceCalculator_;
          mutable double distance_=0;
         public:
          explicit PeriodicQSurfaceCrossObserver(const State2& s_init): periodicQDistanceCalculator_{s_init.q()}
          {};
          bool operator()(const State2& next_point) const noexcept ;
          double distance() const noexcept ;

        };


    }
}
#endif //HAMILTONIANS_PERIODIC_Q_SURFACE_HPP
