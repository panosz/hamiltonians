//
// Created by Panagiotis Zestanakis on 03/10/18.
//

#ifndef HAMILTONIANS_OBSERVER_HPP
#define HAMILTONIANS_OBSERVER_HPP

#include "State.hpp"
#include "line.hpp"
#include <vector>
#include <iostream>
#include <boost/range/algorithm/find_if.hpp>

namespace Integrators
{
    namespace Observer
    {



        class PushBackObserver {
         public:
          using value_type = Geometry::State2_Extended;
         private:

          std::vector<value_type> s_{};
          size_t every_ = 1;
          mutable size_t count = 0;
         public:
          PushBackObserver () =default;
          PushBackObserver (const PushBackObserver&) = default;
          PushBackObserver (PushBackObserver&&) noexcept = default;

          explicit PushBackObserver (size_t every);
          void operator() (value_type s);
          std::vector<value_type> observations() const noexcept;

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

            const auto s_out_extended = stepOnFunctor_(s, t, distance);
            if (validCrossingPredicate_(s_out_extended))
              {
                pushBackObserver_(s_out_extended);
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

          bool operator() (const std::pair<Geometry::State2_Action, double>& s_t)
          {
            const auto&[s, t] = s_t;
            return operator()(s, t);
          }

          auto observations() const noexcept
          {
            return pushBackObserver_.observations();
          }
        };

        template<typename StepOnFunctor, typename SurfaceFunctor, typename ValidCrossingPredicate>
        auto makeCrossSurfaceObserver (StepOnFunctor stepOnFunctor, SurfaceFunctor sf, ValidCrossingPredicate vcp,
                                        size_t every = 0)
        {
          PushBackObserver pbo(every);
          return CrossSurfaceObserver<StepOnFunctor, SurfaceFunctor, ValidCrossingPredicate>(stepOnFunctor, sf, pbo, vcp);
        }

        template<typename StepOnFunctor>
        auto
        makeCrossLineObserver (StepOnFunctor stepOnFunctor, const Geometry::Line& line,
                               size_t every = 0)
        {
          return makeCrossSurfaceObserver(stepOnFunctor, line, [] (auto&)
          { return true; }, every);
        }

        /// \brief Aplies the observer on the integration_range
        /// \tparam Observer a type defining a bool operator() (const IntegrationRange::value_type & s_t)
        /// \tparam IntegrationRange a boost range type
        /// \param observer
        /// \param integration_range

        template<typename Observer, typename IntegrationRange>
        void observe (Observer& observer, const IntegrationRange& integration_range)
        {
          for (const auto& s_t: integration_range)
            observer(s_t);
        }



        /// \brief Aplies the observer on the integration_range until the observer returns true
        /// \tparam Observer a type defining a bool operator() (const IntegrationRange::value_type & s_t)
        /// \tparam IntegrationRange a boost range type
        /// \param observer
        /// \param integration_range
        template<typename Observer, typename IntegrationRange>
        void observe_if (Observer& observer, const IntegrationRange& integration_range)
        {

          boost::range::find_if(integration_range, std::ref(observer));
        }


    }
}

#endif //HAMILTONIANS_OBSERVER_HPP
