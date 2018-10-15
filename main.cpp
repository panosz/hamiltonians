#include <iostream>
#include "line.hpp"
#include "State.hpp"
#include "Hamiltonian.hpp"
#include "dynamic_system.hpp"
#include "observer.hpp"
#include "Integration.hpp"
#include "PoincareSurface.hpp"

#include <boost/math/constants/constants.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/iterator/times_time_iterator.hpp>
#include <boost/range.hpp>
#include <boost/range/algorithm/find_if.hpp>
#include "myUtilities.hpp"

#include <armadillo>

using namespace Integrators;
using namespace Integrators::Geometry;
using namespace boost::numeric::odeint;
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

template<typename Ham>
inline auto
make_integration_range (const Ham& hamiltonian,
                        State2_Action& s_start,
                        const IntegrationOptions& options)
{

  const auto controlled_stepper = make_controlled(options.abs_err, options.rel_err, ErrorStepperType());

  auto integration_functor = [&hamiltonian] (const Geometry::State2_Action& s, Geometry::State2_Action& dsdt, double /*t*/)
  {
      dsdt = Integrators::Dynamics::dynamic_system_Action(hamiltonian, s);
  };

  return boost::make_iterator_range(
      make_adaptive_time_range(controlled_stepper,
                               integration_functor,
                               s_start, 0.0,
                               options.integration_time,
                               options.initial_time_step));

}

template<typename Ham>
inline auto
make_interval_range (const Ham& hamiltonian,
                     State2_Action& s_start,
                     const std::vector<double>& times,
                     const IntegrationOptions& options)
{

  const auto controlled_stepper = make_controlled(options.abs_err, options.rel_err, ErrorStepperType());

  auto integration_functor = [&hamiltonian] (const Geometry::State2_Action& s, Geometry::State2_Action& dsdt, double /*t*/)
  {
      dsdt = Integrators::Dynamics::dynamic_system_Action(hamiltonian, s);
  };

  return boost::make_iterator_range(
      make_times_time_range(controlled_stepper,
                            integration_functor,
                            s_start, times.begin(), times.end(), options.initial_time_step));

}

template<typename PoincareSurfaceType>
inline auto make_init_surface_observer (const PoincareSurfaceType& poincareSurface,
                                        std::vector<Geometry::State2_Action>& s_out,
                                        std::vector<double>& t_out)
{
  auto action_functor = [&poincareSurface] (State2_Action s, double t, double current_distance)
  {
      return poincareSurface.step_back(s, t, current_distance);
  };

  const auto keep_all = [](auto &){return true;};

  return Integrators::Observer::makeCrossSurfaceObserver(action_functor, poincareSurface.cross_line(),keep_all, s_out, t_out);
}

template<typename Observer, typename IntegrationRange>
void observe_first (Observer& observer, const IntegrationRange& integration_range)
{

  boost::range::find_if(integration_range, observer);
}

template<typename Observer, typename IntegrationRange>
void observe (Observer& observer, const IntegrationRange& integration_range)
{
  for (const auto& s_t: integration_range)
    observer(s_t);
}

template<typename Ham>
std::pair<std::vector<State2_Action>, std::vector<double> >
calculate_crossings (const Ham& hamiltonian, Geometry::State2 s_start, const IntegrationOptions& options)
{

  const auto poincareSurface = make_PoincareSurface(hamiltonian, s_start, ErrorStepperType_Extended());

  const auto integration_range = make_integration_range(hamiltonian, s_start, options);

  std::vector<Geometry::State2_Action> s_out;
  std::vector<double> t_out;
  auto observer = make_init_surface_observer(poincareSurface, s_out, t_out);

  observe(observer, integration_range);

  return std::make_pair(s_out, t_out);
}

template<typename Ham>
std::pair<std::vector<State2_Action>, std::vector<double> >
calculate_first_crossing (const Ham& hamiltonian, Geometry::State2 s_start, const IntegrationOptions& options)
{

  const auto poincareSurface = make_PoincareSurface(hamiltonian, s_start,ErrorStepperType_Extended());


  Geometry::State2_Action s_start_Action{s_start,0};

   auto integration_range = make_integration_range(hamiltonian, s_start_Action, options);

  std::vector<Geometry::State2_Action> s_out;
  std::vector<double> t_out;
  auto observer = make_init_surface_observer(poincareSurface, s_out, t_out);

  observe_first(observer, integration_range);

  return std::make_pair(s_out, t_out);
}

int main ()
{

  const auto a = 1;
  const auto b = 1;

  const auto position = Integrators::Geometry::State2{1, -1};
  const auto perp_v = Integrators::Geometry::State2{a, b};
  const auto line = Integrators::Geometry::Line{position, perp_v};
  std::cout << line << std::endl;

  const auto pos = State2{0, 0};
  std::cout << "line value at " << pos << " is " << line(pos) << std::endl;

  {//integrate with adaptive step using make_controlled
    const State2 s_start{0.472035, 7.86664};

    const auto ham = Integrators::Hamiltonian::DuffingHamiltonian();

    IntegrationOptions options;

    options.set_integration_time(100.0);


//    const auto poincareSurface = make_PoincareSurface(ham, s_start,ErrorStepperType_Extended());


    Geometry::State2_Action s_start_Action{s_start,0};

    const auto[s_out, t_out] = calculate_first_crossing(ham, s_start, options);
//
    std::cout << "adaptive integration output:\n";

    for (size_t i = 0; i < t_out.size(); ++i)
      std::cout << t_out[i] << " " << s_out[i] << '\n';


    State2_Action s2 = s_start;

    auto times = PanosUtilities::linspace(0.0, t_out[0], 100);

    std::cout << "times:\n";
    for (const auto& t:times)
      std::cout << t << '\n';

    auto orbit_range = make_interval_range(ham, s2, times, options);

    std::cout << "orbit range:\n";
    for (const auto& s_t : orbit_range)
      {
        std::cout << s_t.first << " , " << s_t.second << '\n';
      }



    return 0;
  }

}