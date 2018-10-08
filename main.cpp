#include <iostream>
#include "line.hpp"
#include "State2.hpp"
#include "Hamiltonian.hpp"
#include "dynamic_system.hpp"
#include "observer.hpp"
#include <boost/math/constants/constants.hpp>
#include <boost/numeric/odeint.hpp>

using namespace Integrators::Geometry;
using namespace boost::numeric::odeint;


using ErrorStepperType = boost::numeric::odeint::runge_kutta_cash_karp54<State2, double, State2,
    double, boost::numeric::odeint::vector_space_algebra>;
using ControlledStepperType = controlled_runge_kutta<ErrorStepperType>;

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
    State2 s_start{1, 1};
    const double abs_err = 1.0e-16;
    const double rel_err = 1.0e-14;
    auto controlled_stepper = make_controlled(abs_err, rel_err, ErrorStepperType());
    std::vector<double> t_out;
    std::vector<State2> s_out;

    const auto system = Integrators::Dynamics::system<Integrators::Hamiltonian::HarmonicOscillator>;
    Integrators::Geometry::State2 start_direction{};
    system(s_start, start_direction, 0);

    const auto start_line = Integrators::Geometry::Line{s_start, start_direction};
    const auto systemAlongDirection = Integrators::Dynamics::SystemAlongDirection<Integrators::Hamiltonian::HarmonicOscillator>
        (start_direction);

    const auto surface_functor = [&start_line] (const State2& s)
    { return start_line(s); };

    const auto push_back_observer = Integrators::Observer::PushBackObserver(s_out, t_out, 0);

    const auto step_back_functor = [&systemAlongDirection,&push_back_observer, &surface_functor] (State2 s, double t)
    {

        ErrorStepperType().do_step(
            systemAlongDirection,
            s,
            t,
            -surface_functor(s));

        push_back_observer(s, t);
    };

    auto observer = Integrators::Observer::CrossZeroObserver(step_back_functor, surface_functor);

    integrate_adaptive(controlled_stepper,
                       system, s_start, 0.0, 100.0, 0.01,
                       observer);
    //Integrators::Observer::PushBackObserver(s_out,t_out,10) );


    std::cout << "adaptive integration output:\n";
    for (int i = 0; i < t_out.size(); ++i)
      std::cout << t_out[i]/boost::math::double_constants::two_pi << " " << s_out[i] << '\n';

    std::cout << "final state = " << s_start << '\n';
  }

  return 0;
}