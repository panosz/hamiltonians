#include <iostream>
#include "line.hpp"
#include "Hamiltonian.hpp"
#include "dynamic_system.hpp"
#include "observer.hpp"
#include "Integration.hpp"
#include "PoincareSurface.hpp"

#include <boost/range/algorithm/find_if.hpp>
#include "myUtilities.hpp"

using namespace Integrators;
using namespace Integrators::Geometry;
using namespace boost::numeric::odeint;








template<typename Ham>
std::pair<std::vector<State2_Action>, std::vector<double> >
calculate_crossings (const Ham& hamiltonian, Geometry::State2 s_start, const IntegrationOptions& options)
{

  const auto system = Dynamics::DynamicSystem{hamiltonian};
  const auto cross_line = Integrators::make_init_cross_line(system,s_start);
  std::vector<Geometry::State2_Action> s_out;
  std::vector<double> t_out;
  auto observer = Integrators::make_init_surface_observer(system,cross_line, s_out, t_out);


  Geometry::State2_Action s_start_Action{s_start};
  s_start_Action.J() = 0;

  const auto integration_range = Integrators::make_integration_range(system, s_start_Action, options);



  observe(observer, integration_range);

  return std::make_pair(s_out, t_out);
}

template<typename Ham>
std::pair<std::vector<State2_Action>, std::vector<double> >
calculate_first_crossing (const Ham& hamiltonian, const Geometry::State2& s_start, const IntegrationOptions& options)
{

  const auto system = Dynamics::DynamicSystem{hamiltonian};
  const auto cross_line = Integrators::make_init_cross_line(system,s_start);
  std::vector<Geometry::State2_Action> s_out;
  std::vector<double> t_out;
  auto observer = Integrators::make_init_surface_observer(system,cross_line, s_out, t_out);


  Geometry::State2_Action s_start_Action{s_start};
  s_start_Action.J() = 0;

  const auto integration_range = Integrators::make_integration_range(system, s_start_Action, options);

  observe_if(observer, integration_range);

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


    Geometry::State2_Action s_start_Action{s_start};
    s_start_Action.J() = 0;

    const auto[s_out, t_out] = calculate_first_crossing(ham, s_start, options);
//
    std::cout << "adaptive integration output:\n";

    for (size_t i = 0; i < t_out.size(); ++i)
      std::cout << t_out[i] << " " << s_out[i] << '\n';

    State2_Action s2{s_start};
    s2.J() = 0;

    auto times = PanosUtilities::linspace(0.0, t_out[0], 100);

    std::cout << "times:\n";
    for (const auto& t:times)
      std::cout << t << '\n';


    auto orbit_range = make_interval_range(Dynamics::DynamicSystem(ham), s2, times, options);

    std::cout << "orbit range:\n";
    for (const auto& s_t : orbit_range)
      {
        std::cout << s_t.first << " , " << s_t.second << '\n';
      }

    return 0;
  }

}