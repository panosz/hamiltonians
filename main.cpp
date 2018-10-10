#include <iostream>
#include "line.hpp"
#include "State.hpp"
#include "Hamiltonian.hpp"
#include "dynamic_system.hpp"
#include "observer.hpp"
#include "Integration.hpp"
#include <boost/math/constants/constants.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/range.hpp>
#include <boost/range/algorithm/find_if.hpp>

using namespace Integrators;
using namespace Integrators::Geometry;
using namespace boost::numeric::odeint;
using ErrorStepperType = boost::numeric::odeint::runge_kutta_cash_karp54<Geometry::State2, double, Geometry::State2,
    double, boost::numeric::odeint::vector_space_algebra>;

using ErrorStepperType_Extended = boost::numeric::odeint::runge_kutta_cash_karp54<Geometry::State3, double, Geometry::State3,
    double, boost::numeric::odeint::vector_space_algebra>;

using ControlledStepperType = boost::numeric::odeint::controlled_runge_kutta<ErrorStepperType>;

template<typename Ham>
class ClosedLineIntegratorStuff {
 private:
  Ham ham_;
  Geometry::State2 s_start_;
  Geometry::State2 direction_;
  Geometry::Line cross_line_;

 public:
  inline void
  directional_functor (const Geometry::State3& s_extended, Geometry::State3& dsdt_extended, double /*t*/) const
  {

    Geometry::State2 s{s_extended};
    dsdt_extended = Integrators::Dynamics::dynamic_system_along_direction(ham_, direction_, s);

  };

  inline void non_directional_functor (const Geometry::State2& s, Geometry::State2& dsdt, double /*t*/) const
  {
    dsdt = Integrators::Dynamics::dynamic_system(ham_, s);
  }

  std::pair<State2, double> step_back (State2 s, double t, double current_distance) const
  {

    auto state_extended = State3{s, t};

    auto df = [this] (const Geometry::State3& s_extended, Geometry::State3& dsdt_extended, double t)
    {
        directional_functor(s_extended, dsdt_extended, t);
    };

    ErrorStepperType_Extended().do_step(
        df,
        state_extended,
        t,
        -current_distance);

    s = State2{state_extended};
    t = state_extended.t();

    return std::make_pair(s, t);
  };

  Geometry::Line cross_line () const
  { return cross_line_; };

  ClosedLineIntegratorStuff (Ham ham, Geometry::State2 s_start, Geometry::State2 direction, Geometry::Line cross_line)
      : ham_{ham}, s_start_{s_start}, direction_{direction}, cross_line_{cross_line}
  {
  };
};

template<typename Ham>
inline auto make_ClosedLineIntegratorStuff (Ham hamiltonian, Geometry::State2 s_start)
{
  const auto start_direction = Integrators::Dynamics::dynamic_system(hamiltonian, s_start);
  const auto start_line = Integrators::Geometry::Line(s_start, start_direction);
  return ClosedLineIntegratorStuff<Ham>(hamiltonian, s_start, start_direction, start_line);

}

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
make_integration_range (const ClosedLineIntegratorStuff<Ham>& integrator_stuff, State2& s_start, const IntegrationOptions& options)
{

  const auto controlled_stepper = make_controlled(options.abs_err, options.rel_err, ErrorStepperType());

  auto integration_functor = [&integrator_stuff] (const Geometry::State2& s, Geometry::State2& dsdt, double t)
  { integrator_stuff.non_directional_functor(s, dsdt, t); };

  return boost::make_iterator_range(
      make_adaptive_time_range(controlled_stepper,
                               integration_functor,
                               s_start, 0.0,
                               options.integration_time,
                               options.initial_time_step));

}

template<typename Ham>
inline auto make_init_surface_observer (const ClosedLineIntegratorStuff<Ham>& integrator_stuff,
                                        std::vector<Geometry::State2>& s_out,
                                        std::vector<double>& t_out)
{
  auto action_functor = [&integrator_stuff] (State2 s, double t, double current_distance)
  {
      return integrator_stuff.step_back(s, t, current_distance);
  };

  return Integrators::Observer::makeCrossSurfaceObserver(action_functor, integrator_stuff.cross_line(), s_out, t_out);
}

template<typename Observer, typename IntegrationRange>
void observe_first (Observer& observer, const IntegrationRange& integration_range)
{

  boost::range::find_if(integration_range, observer);
};

template<typename Observer, typename IntegrationRange>
void observe (Observer& observer, const IntegrationRange& integration_range)
{
  for (const auto& s_t: integration_range)
    observer(s_t);
};

template<typename Ham>
std::pair<std::vector<State2>, std::vector<double> >
calculate_crossings (const Ham& hamiltonian, Geometry::State2 s_start, const IntegrationOptions& options)
{

  const auto integrator_stuff = make_ClosedLineIntegratorStuff(hamiltonian, s_start);

  const auto integration_range = make_integration_range(integrator_stuff, s_start, options);

  std::vector<Geometry::State2> s_out;
  std::vector<double > t_out;
  auto observer = make_init_surface_observer(integrator_stuff,s_out, t_out);

  observe(observer,integration_range);

  return std::make_pair(s_out,t_out);
};

template<typename Ham>
std::pair<std::vector<State2>, std::vector<double> >
calculate_first_crossing (const Ham& hamiltonian, Geometry::State2 s_start, const IntegrationOptions& options)
{

  const auto integrator_stuff = make_ClosedLineIntegratorStuff(hamiltonian, s_start);

  const auto integration_range = make_integration_range(integrator_stuff, s_start, options);

  std::vector<Geometry::State2> s_out;
  std::vector<double > t_out;
  auto observer = make_init_surface_observer(integrator_stuff,s_out, t_out);

  observe_first(observer,integration_range);

  return std::make_pair(s_out,t_out);
};



int main ()
{

  const auto a = 1;
  const auto b = 1;
  const auto c = -2;

  const auto position = Integrators::Geometry::State2{1, -1};
  const auto perp_v = Integrators::Geometry::State2{a, b};
  const auto line = Integrators::Geometry::Line{position, perp_v};
  std::cout << line << std::endl;

  const auto pos = State2{0, 0};
  std::cout << "line value at " << pos << " is " << line(pos) << std::endl;

  {//integrate with adaptive step using make_controlled
    const State2 s_start{1, 1};

    const auto ham = Integrators::Hamiltonian::HarmonicOscillator();

    IntegrationOptions options;

    options.set_integration_time(100.0);

    const auto[s_out, t_out] = calculate_first_crossing(ham,s_start,options);

    std::cout << "adaptive integration output:\n";

    for (int i = 0; i < t_out.size(); ++i)
      std::cout << t_out[i] / boost::math::double_constants::two_pi << " " << s_out[i] << '\n';
    }

  return 0;
}