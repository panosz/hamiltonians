//
// Created by Panagiotis Zestanakis on 03/10/18.
//

#include "observer.hpp"
#include <boost/math/constants/constants.hpp>

namespace Integrators
{
    namespace Observer
    {

        void PushBackObserver::operator() (Geometry::State2_Extended s)
        {
          if (!(count++ % every_))
            {
              s_.push_back(s);

              }

        }
        PushBackObserver::PushBackObserver (std::vector<Geometry::State2_Extended>& s, size_t every)
            :s_{s},
            every_{every > 0 ? every : 1}
        {

        }

        bool crossZeroPositiveDirectionPredicate (double current_value, double previous_value)
        {
          return (current_value >= 0 && previous_value< 0);
        }
    }
}