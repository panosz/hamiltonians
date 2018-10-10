//
// Created by Panagiotis Zestanakis on 03/10/18.
//

#ifndef HAMILTONIANS_OBSERVER_HPP
#define HAMILTONIANS_OBSERVER_HPP

#include "State.hpp"
#include "line.hpp"
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
          std::vector<Geometry::State2>& s_;
          std::vector<double>& t_;
          size_t every_{};
          mutable size_t count = 0;
         public:
          PushBackObserver ()  = delete;
          PushBackObserver (const PushBackObserver&) = default;
          PushBackObserver (PushBackObserver&&) noexcept = default;

          PushBackObserver (std::vector<Geometry::State2>& s, std::vector<double>& t,size_t every);
          void operator() (Geometry::State2 s, double t);
          const std::vector<Geometry::State2>& s () const
          {
            return s_;
          };

          const std::vector<double>& t() const
          {
            return t_;
          }


        };

        bool crossZeroPositiveDirectionPredicate (double current_value, double previous_value);

        template<typename ActionFunctor, typename SurfaceFunctor>
        class CrossSurfaceObserver {
         private:

          ActionFunctor actionFunctor_;
          SurfaceFunctor surfaceFunctor_;
          PushBackObserver pushBackObserver_;

          double previous_distance_from_surface_{0};
          double distance_from_surface (const Geometry::State2& s) const
          {
            return surfaceFunctor_(s);
          }

          bool event (double distance_from_surface) const
          {
            return crossZeroPositiveDirectionPredicate(distance_from_surface, previous_distance_from_surface_);
          }

          void action (const Geometry::State2& s, double t, double distance)
          {
            const auto [s_out, t_out] = actionFunctor_(s, t, distance);
            pushBackObserver_(s_out,t_out);
          }
         public:
          CrossSurfaceObserver () = delete;

          CrossSurfaceObserver (ActionFunctor af, SurfaceFunctor sf, PushBackObserver pbo)
              : actionFunctor_{std::move(af)},
                surfaceFunctor_{std::move(sf)},
                pushBackObserver_{pbo}
          { };

          bool operator() (const Geometry::State2& s, double t)
          {
            bool event_occured = false;

            const auto distance = distance_from_surface(s);
            if (event(distance))
              {
                action(s, t, distance);
                event_occured = true;
              }
            previous_distance_from_surface_ = distance;
            return event_occured;
          }

          bool operator() (const std::pair<Geometry::State2,double>& s_t)
          {
            const auto& [s,t] = s_t;
            return operator()(s,t);
          }

        };

        template<typename ActionFunctor, typename SurfaceFunctor>
        auto makeCrossSurfaceObserver (ActionFunctor af, SurfaceFunctor sf,
                                       std::vector<Geometry::State2>& s_out,
                                       std::vector<double>& t_out, size_t every =0)
        {
          PushBackObserver pbo(s_out,t_out,every);
          return CrossSurfaceObserver<ActionFunctor, SurfaceFunctor>(af, sf,pbo);
        };

        template<typename ActionFunctor>
        auto makeCrossLineObserver (ActionFunctor af, const Geometry::Line& line, std::vector<Geometry::State2>& s_out,
                                    std::vector<double>& t_out, size_t every =0)
        {
          return makeCrossSurfaceObserver(af, line,s_out,t_out, every);
        }

    }
}

#endif //HAMILTONIANS_OBSERVER_HPP
