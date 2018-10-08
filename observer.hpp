//
// Created by Panagiotis Zestanakis on 03/10/18.
//

#ifndef HAMILTONIANS_OBSERVER_HPP
#define HAMILTONIANS_OBSERVER_HPP

#include "State2.hpp"
#include <vector>
#include <iostream>

namespace Integrators
{
    namespace Observer
    {

        class NoOpFunctor {
         public:
          void operator() (Geometry::State2&, double&)
          {
          };

        };

        class PushBackObserver {
         private:
          std::vector<Geometry::State2>& sout_;
          std::vector<double>& tout_;
          const size_t every_{};
          mutable size_t count = 0;
         public:
          PushBackObserver () noexcept = delete;
          PushBackObserver (const PushBackObserver&) noexcept = default;
          PushBackObserver (PushBackObserver&&) noexcept = default;

          PushBackObserver (std::vector<Geometry::State2>& sout, std::vector<double>& tout, size_t every);
          void operator() (Geometry::State2 s, double t) const;

        };

        template<typename ActionFunctor, typename SurfaceFunctor>
        class CrossZeroObserver {
         private:

          ActionFunctor actionFunctor_;
          SurfaceFunctor surfaceFunctor_;

          double previous_value_{0};
          double calculate_value (const Geometry::State2& s) const
          {
            return surfaceFunctor_(s);
          }

          bool event (double curr_value) const
          {
            return (curr_value >= 0 && previous_value_ < 0);
          }

          void action (const Geometry::State2& s, double t) const
          {
            actionFunctor_(s, t);
          }
         public:
          CrossZeroObserver () = delete;

          CrossZeroObserver (ActionFunctor af, SurfaceFunctor sf)
              :
              actionFunctor_{std::move(af)}, surfaceFunctor_{std::move(sf)}
          { };

          void operator() (const Geometry::State2& s, double t)
          {
            const auto curr_value = calculate_value(s);
            if (event(curr_value))
              action(s, t);
            previous_value_ = curr_value;
          }

        };

        template<typename ActionFunctor, typename SurfaceFunctor>
        auto makeCrossZeroObserver (ActionFunctor af, SurfaceFunctor sf)
        {
          return CrossZeroObserver<ActionFunctor, SurfaceFunctor>(af, sf);
        };

    }
}

#endif //HAMILTONIANS_OBSERVER_HPP
