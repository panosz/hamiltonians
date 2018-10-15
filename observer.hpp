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
          void operator() (Geometry::State2_Action&, double&)
          {
          };

        };

        class PushBackObserver {
         private:
          std::vector<Geometry::State2_Action>& s_;
          std::vector<double>& t_;
          size_t every_{};
          mutable size_t count = 0;
         public:
          PushBackObserver ()  = delete;
          PushBackObserver (const PushBackObserver&) = default;
          PushBackObserver (PushBackObserver&&) noexcept = default;

          PushBackObserver (std::vector<Geometry::State2_Action>& s, std::vector<double>& t,size_t every);
          void operator() (Geometry::State2_Action s, double t);


        };

        bool crossZeroPositiveDirectionPredicate (double current_value, double previous_value);

        template<typename StepOnFunctor, typename DistanceFunctor, typename FilterObservationPredicate>
        class CrossSurfaceObserver {
         private:

          StepOnFunctor stepOnFunctor_;
          DistanceFunctor distanceFunctor_;
          PushBackObserver pushBackObserver_;
          FilterObservationPredicate validCrossingPredicate_;


          double previous_distance_from_surface_{0};


          bool crossing_detected (double distance_from_surface) const
          {
            return crossZeroPositiveDirectionPredicate(distance_from_surface, previous_distance_from_surface_);
          }

          /// \brief carries out the necessary operations after a crossing has been detected
          /// \param s the current position
          /// \param t the current time
          /// \param distance the distance from the surface
          /// \return true, if the crossing has been accepted and forwarded to the pushBackObserver
          bool after_crossing_action (const Geometry::State2_Action& s, double t, double distance)
          {

            const auto [s_out, t_out] = stepOnFunctor_(s, t, distance);
            if ( validCrossingPredicate_(s_out) )
              {
                pushBackObserver_(s_out, t_out);
                return true;
              }
            return false;
          }
         public:
          CrossSurfaceObserver () = delete;

          CrossSurfaceObserver (StepOnFunctor af, DistanceFunctor sf, PushBackObserver pbo, FilterObservationPredicate fop)
              : stepOnFunctor_{std::move(af)},
                distanceFunctor_{std::move(sf)},
                pushBackObserver_{pbo},
                validCrossingPredicate_{fop}
          { };

          bool operator() (const Geometry::State2_Action& s, double t)
          {
            bool crossing_accepted = false;
            const auto distance = distanceFunctor_(Geometry::State2{s});

            if (crossing_detected(distance))
                crossing_accepted = after_crossing_action(s, t, distance);


            previous_distance_from_surface_ = distance;
            return crossing_accepted;
          }

          bool operator() (const std::pair<Geometry::State2_Action,double>& s_t)
          {
            const auto& [s,t] = s_t;
            return operator()(s,t);
          }

        };

        template<typename StepOnFunctor, typename SurfaceFunctor, typename ValidCrossingPredicate>
        auto makeCrossSurfaceObserver (StepOnFunctor stepOnFunctor, SurfaceFunctor sf, ValidCrossingPredicate vcp,
                                       std::vector<Geometry::State2_Action>& s_out,
                                       std::vector<double>& t_out, size_t every =0)
        {
          PushBackObserver pbo(s_out,t_out,every);
          return CrossSurfaceObserver<StepOnFunctor, SurfaceFunctor, ValidCrossingPredicate>(stepOnFunctor, sf,pbo, vcp);
        }

        template<typename StepOnFunctor>
        auto makeCrossLineObserver (StepOnFunctor stepOnFunctor, const Geometry::Line& line, std::vector<Geometry::State2_Action>& s_out,
                                    std::vector<double>& t_out, size_t every =0)
        {
          return makeCrossSurfaceObserver(stepOnFunctor, line,[](auto &){return true;},s_out,t_out, every);
        }

    }
}

#endif //HAMILTONIANS_OBSERVER_HPP
