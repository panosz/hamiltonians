//
// Created by Panagiotis Zestanakis on 03/10/18.
//

#include "observer.hpp"
#include <utility>

namespace Integrators
{
    namespace Observer
    {
        PushBackObserver::PushBackObserver (std::vector<Geometry::State2>& s_out, std::vector<double>& t_out, size_t every)
            :
            sout_{s_out}, tout_{t_out}, every_{every > 0 ? every : 1}
        { }


        void PushBackObserver::operator() (Geometry::State2 s, double t) const
        {
          if (!(count++ % every_))
            {
              sout_.push_back(s);
              tout_.push_back(t);
            }

        }
    }
}