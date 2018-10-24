//
// Created by Panagiotis Zestanakis on 08/10/18.
//

#ifndef HAMILTONIANS_INTEGRATION_HPP
#define HAMILTONIANS_INTEGRATION_HPP

#include "State.hpp"
#include "line.hpp"
#include "dynamic_system.hpp"
#include "observer.hpp"

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/iterator/times_time_iterator.hpp>

namespace Integrators
{
    template< typename StateType>
    using ErrorStepperType = boost::numeric::odeint::runge_kutta_cash_karp54<StateType, double, StateType,
        double, boost::numeric::odeint::vector_space_algebra>;

    template<typename StateType>
    using ControlledStepperType = boost::numeric::odeint::controlled_runge_kutta<ErrorStepperType< StateType> >;



    struct IntegrationOptions {
        double abs_err = 1.0e-16;
        double rel_err = 1.0e-14;
        double integration_time = 0;
        double initial_time_step = 1e-5;
        double distance_threshold = 1.0e-13;

        IntegrationOptions () = default;

        void set_abs_err (double error)
        {
          abs_err = error;
        }
        void set_rel_err (double error)
        {
          rel_err = error;
        }
        void set_integration_time (double i_time)
        {
          integration_time = i_time;
        }
        void set_initial_time_step (double dt_init)
        {
          initial_time_step = dt_init;
        }

        void set_distance_threshold (double ds)
        {
          distance_threshold = ds;
        }

    };

    template<typename DS>
    Geometry::State2_Extended step_back (const DS& system,
                                                          const Geometry::State2 direction,
                                                          const Geometry::State2_Action& s,
                                                          double t,
                                                          double distance)
    {
      using ErrorStepperType_Extended = ErrorStepperType<Geometry::State2_Extended >;

      auto state_extended = Geometry::State2_Extended{s};

      state_extended.t() = t;

      auto df = [&system, &direction] (const Geometry::State2_Extended& s_extended, Geometry::State2_Extended& dsdt_extended, double /*time*/)
      {

          dsdt_extended = system.dynamic_system_along_direction(direction,
                                                                Geometry::State2_Action{s_extended});
      };

      ErrorStepperType_Extended().do_step(
          df,
          state_extended,
          t,
          -distance);

      return state_extended;
    }

    template<typename DS>
    inline auto
    make_dynamic_system_integration_range (DS system, // not const &, see comment below
                                           Geometry::State2_Action& s_start,
                                           const IntegrationOptions& options)
    {

      const auto controlled_stepper = make_controlled(options.abs_err,
                                                      options.rel_err,
                                                      ErrorStepperType<Geometry::State2_Action >());


      //system should be passed by value to the closure, because integration_functor is coppied into the output range
      //and reference may dangle
      auto integration_functor = [sys = std::move(system)] (const Geometry::State2_Action& s, Geometry::State2_Action& dsdt, double /*t*/)
      {
          dsdt = sys.dynamic_system_Action(s);
      };

      return boost::make_iterator_range(
          make_adaptive_time_range(controlled_stepper,
                                   integration_functor,
                                   s_start, 0.0,
                                   options.integration_time,
                                   options.initial_time_step));

    }

    template<typename DS>
    inline auto
    make_interval_range (DS system,   // not const &, see comment below
                         Geometry::State2_Action& s_start,
                         const std::vector<double>& times,
                         const IntegrationOptions& options)
    {

      const auto controlled_stepper = make_controlled(options.abs_err,
                                                      options.rel_err,
                                                      ErrorStepperType<Geometry::State2_Action>());


      //system should be passed by value to the closure, because integration_functor is coppied into the output range
      //and reference may dangle

      auto integration_functor = [sys = std::move(system)]
          (const Geometry::State2_Action& s, Geometry::State2_Action& dsdt, double /*t*/)
      {
          dsdt = sys.dynamic_system_Action(s);
      };

      return boost::make_iterator_range(
          boost::numeric::odeint::make_times_time_range(controlled_stepper,
                                                        integration_functor,
                                                        s_start, times.begin(), times.end(), options.initial_time_step));

    }

    template<typename DS>
    inline auto
    make_interval_range (DS system,   // not const &, see comment below
                         Geometry::State2& s_start,
                         const std::vector<double>& times,
                         const IntegrationOptions& options)
    {

      const auto controlled_stepper = make_controlled(options.abs_err,
                                                      options.rel_err,
                                                      ErrorStepperType<Geometry::State2>());


      //system should be passed by value to the closure, because integration_functor is coppied into the output range
      //and reference may dangle

      auto integration_functor = [sys = std::move(system)]
          (const Geometry::State2& s, Geometry::State2& dsdt, double /*t*/)
      {
          dsdt = sys.dynamic_system(s);
      };

      return boost::make_iterator_range(
          boost::numeric::odeint::make_times_time_range(controlled_stepper,
                                                        integration_functor,
                                                        s_start, times.begin(), times.end(), options.initial_time_step));

    }


    template<typename DS>
    Geometry::Line make_init_cross_line (const DS& system, Geometry::State2 s_start)
    {
      const auto start_direction = system.dynamic_system(s_start);

      return Integrators::Geometry::Line(s_start, start_direction);
    }

    template<typename DS, typename FP>
    inline auto make_project_on_line_observer (DS system, // not const &. may dangle
                                               const Geometry::Line& line,
                                               FP filteringPredicate)
    {

      auto action_functor =
          [sys = std::move(system), direction = line.perpendicular_vector()]
              (Geometry::State2_Action s, double t, double current_distance)
          {
              return step_back(sys, direction, s, t, current_distance);
          };



      return Observer::makeProjectOnSurfaceObserver(action_functor, Geometry::LineCrossObserver(line), filteringPredicate);
    }


    class StateNear {
      double distance_;
      bool are_near (const Geometry::State2& s1, const Geometry::State2& s2) const noexcept
      {
        return magnitude(s1 - s2) <= distance_;
      }
     public:
      explicit StateNear (double distance)
          : distance_(distance)
      { }

      bool operator() (const Geometry::State2& s1, const Geometry::State2& s2) const noexcept
      {
        return are_near(s1, s2);
      }

      bool operator() (const Geometry::State2_Extended& s1, const Geometry::State2& s2) const noexcept
      {
        return are_near(Geometry::State2{s1}, s2);
      }

      bool operator() (const Geometry::State2& s1, const Geometry::State2_Extended& s2) const noexcept
      {
        return are_near(s1, Geometry::State2{s2});
      }

    };

    template<typename Ham>
    std::vector<Geometry::State2_Extended>
    calculate_crossings (const Ham& hamiltonian,
                         const Geometry::State2& s_start,
                         const Geometry::Line& cross_line,
                         const IntegrationOptions& options)
    {

      Geometry::State2_Action s_start_Action{s_start};

      const auto system = Dynamics::DynamicSystem{hamiltonian};
      const auto integration_range = Integrators::make_dynamic_system_integration_range(system, s_start_Action, options);

      auto observer = Integrators::make_project_on_line_observer(system, cross_line, [] (auto&)
      { return true; });
      observe(observer, integration_range);

      return observer.observations();
    }

    template<typename Ham>
    Geometry::State2_Extended
    calculate_first_crossing (const Ham& hamiltonian,
                              const Geometry::State2& s_start,
                              const Geometry::Line& cross_line,
                              const IntegrationOptions& options)
    {

      const auto system = Dynamics::DynamicSystem{hamiltonian};
      Geometry::State2_Action s_start_Action{s_start};

      const auto integration_range = Integrators::make_dynamic_system_integration_range(system, s_start_Action, options);

      auto observer = Integrators::make_project_on_line_observer(system, cross_line, [] (auto&)
      { return true; });

      observe_if(observer, integration_range);

      const auto observations = observer.observations();

      if (observations.empty())
        throw std::runtime_error("orbit never came back");

      return observations.front();
      }

    template<typename Ham>
    Geometry::State2_Extended
    come_back_home (const Ham& hamiltonian,
                    const Geometry::State2& s_start,
                    const IntegrationOptions& options)
    {
      const auto system = Dynamics::DynamicSystem{hamiltonian};

      const auto cross_line = make_init_cross_line(system, s_start);

      const auto is_back_predicate =
          [&s_home = s_start,
           & distance_threshold = options.distance_threshold] (auto& s)
          {
              return StateNear(distance_threshold)(s_home, s);
          };

      auto back_home_observer = Integrators::make_project_on_line_observer(system,
                                                                           cross_line,
                                                                           is_back_predicate);

      Geometry::State2_Action s_start_Action{s_start};

      const auto integration_range = Integrators::make_dynamic_system_integration_range(system, s_start_Action, options);

      observe_if(back_home_observer, integration_range);

      const auto observations = back_home_observer.observations();

      if (observations.empty())
        throw std::runtime_error("orbit never came back");

      return observations.front();

    }
}
#endif //HAMILTONIANS_INTEGRATION_HPP
