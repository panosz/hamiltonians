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
        PushBackObserver::PushBackObserver ( size_t every)
            :
            every_{every > 0 ? every : 1}
        {

        }
        std::vector<PushBackObserver::value_type> PushBackObserver::observations () const noexcept
        {
          return s_;
        }

        bool crossZeroPositiveDirectionPredicate (double current_value, double previous_value)
        {
          return (current_value >= 0 && previous_value< 0);
        }
    }
}