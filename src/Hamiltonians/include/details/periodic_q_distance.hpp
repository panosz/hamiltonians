//
// Created by Panagiotis Zestanakis on 23/10/18.
//



#ifndef HAMILTONIANS_PERIODIC_Q_DISTANCE_HPP
#define HAMILTONIANS_PERIODIC_Q_DISTANCE_HPP

namespace Integrators
{
    namespace Internals
    {

        class PeriodicQDistance {
          double q0_ = 0;
         public:
          explicit PeriodicQDistance (double q0);

          double distance (double q) const noexcept;

        };
    }

}

#endif //HAMILTONIANS_PERIODIC_Q_DISTANCE_HPP
