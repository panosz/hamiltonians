//
// Created by Panagiotis Zestanakis on 03/10/18.
//

#include "observer.hpp"
#include <utility>
#include <boost/math/constants/constants.hpp>

namespace Integrators
{
    namespace Observer
    {

        void PushBackObserver::operator() (Geometry::State2_Action s, double t)
        {
          if (!(count++ % every_))
            {
              s_.push_back(s);
              t_.push_back(t);
              }

        }
        PushBackObserver::PushBackObserver (std::vector<Geometry::State2_Action>& s, std::vector<double>& t, size_t every)
            :s_{s},t_{t},
            every_{every > 0 ? every : 1}
        {

        }

        bool crossZeroPositiveDirectionPredicate (double current_value, double previous_value)
        {
          return (current_value >= 0 && previous_value< 0);
        }
    }
}