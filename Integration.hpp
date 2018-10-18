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
    using ErrorStepperType = boost::numeric::odeint::runge_kutta_cash_karp54<Geometry::State2_Action, double, Geometry::State2_Action,
        double, boost::numeric::odeint::vector_space_algebra>;

    using ErrorStepperType_Extended = boost::numeric::odeint::runge_kutta_cash_karp54<Geometry::State2_Extended, double, Geometry::State2_Extended,
        double, boost::numeric::odeint::vector_space_algebra>;

    using ControlledStepperType = boost::numeric::odeint::controlled_runge_kutta<ErrorStepperType>;

    struct IntegrationOptions {
        double abs_err = 1.0e-16;
        double rel_err = 1.0e-14;
        double integration_time = 0;
        double initial_time_step = 1e-5;

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

    };

    template<typename DS>
    std::pair<Geometry::State2_Action, double> step_back (const DS& system,
                                                          const Geometry::State2 direction,
                                                          const Geometry::State2_Action& s,
                                                          double t,
                                                          double distance)
    {

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

      const Geometry::State2_Action s_out{state_extended};
      const double t_out = state_extended.t();
      return std::make_pair(s_out, t_out);
    }

    template<typename DS>
    inline auto
    make_integration_range (DS system, // not const &, see comment below
                            Geometry::State2_Action& s_start,
                            const IntegrationOptions& options)
    {

      const auto controlled_stepper = make_controlled(options.abs_err, options.rel_err, ErrorStepperType());


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

      const auto controlled_stepper = make_controlled(options.abs_err, options.rel_err, ErrorStepperType());


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
    Geometry::Line make_init_cross_line (const DS& system, Geometry::State2 s_start)
    {
      const auto start_direction = system.dynamic_system(s_start);

      return Integrators::Geometry::Line(s_start, start_direction);
    }

    template<typename DS>
    inline auto make_init_surface_observer (DS system, // not const &. may dangle
                                            const Geometry::Line & cross_line,
                                            std::vector<Geometry::State2_Action>& s_out,
                                            std::vector<double>& t_out)
    {

      auto action_functor =
          [sys = std::move(system), direction = cross_line.perpendicular_vector()]
              (Geometry::State2_Action s, double t, double current_distance)
          {
              return step_back(sys, direction, s, t, current_distance);
          };

      const auto keep_all = [] (auto&)
      { return true; };

      return Observer::makeCrossSurfaceObserver(action_functor, cross_line, keep_all, s_out, t_out);
    }

}
#endif //HAMILTONIANS_INTEGRATION_HPP
