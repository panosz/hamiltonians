//
// Created by Panagiotis Zestanakis on 23/10/18.
//



#ifndef HAMILTONIANS_PERIODIC_Q_SURFACE_HPP
#define HAMILTONIANS_PERIODIC_Q_SURFACE_HPP

#include "myUtilities.hpp"

namespace Integrators
{

    class PeriodicQSurface
    {
      double q0_=0;
     public:
      explicit PeriodicQSurface (double q0)
          : q0_(q0)
      { };

      double distance(double q)
      {
        return PanosUtilities::wrap_minus_pi_pi(q-q0_);
      }

    };
}

#endif //HAMILTONIANS_PERIODIC_Q_SURFACE_HPP
