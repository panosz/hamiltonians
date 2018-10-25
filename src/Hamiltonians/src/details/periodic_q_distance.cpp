//
// Created by Panagiotis Zestanakis on 25/10/18.
//

#include <algorithm>
#include "myUtilities/wrap.hpp"
#include "details/periodic_q_distance.hpp"

namespace Integrators
{
    namespace Internals
    {
        PeriodicQDistance::PeriodicQDistance (double q0)
            : q0_(q0)
        { }
        double PeriodicQDistance::distance (double q) const noexcept
        {
          return PanosUtilities::wrap_minus_pi_pi(q - q0_);
        }
    }


}