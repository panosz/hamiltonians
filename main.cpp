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
std::vector<State2_Extended >
calculate_crossings (const Ham& hamiltonian, Geometry::State2 s_start, const IntegrationOptions& options)
{

  std::vector<Geometry::State2_Extended > s_out;

  const auto system = Dynamics::DynamicSystem{hamiltonian};
  const auto cross_line = Integrators::make_init_cross_line(system,s_start);
  auto observer = Integrators::make_init_surface_observer(system,cross_line, s_out);


  Geometry::State2_Action s_start_Action{s_start};
  s_start_Action.J() = 0;

  const auto integration_range = Integrators::make_integration_range(system, s_start_Action, options);



  observe(observer, integration_range);

  return s_out;
}

template<typename Ham>
std::vector<State2_Extended>
calculate_first_crossing (const Ham& hamiltonian, const Geometry::State2& s_start, const IntegrationOptions& options)
{
  std::vector<Geometry::State2_Extended> s_out{};

  const auto system = Dynamics::DynamicSystem{hamiltonian};
  const auto cross_line = Integrators::make_init_cross_line(system,s_start);
  auto observer = Integrators::make_init_surface_observer(system,cross_line, s_out);


  Geometry::State2_Action s_start_Action{s_start};
  s_start_Action.J() = 0;

  const auto integration_range = Integrators::make_integration_range(system, s_start_Action, options);

  observe_if(observer, integration_range);

  return s_out;
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


    Geometry::State2_Action s_start_Action{s_start};
    s_start_Action.J() = 0;

    const auto s_out_extended = calculate_first_crossing(ham, s_start, options);
//
    std::cout << "adaptive integration output:\n";

    for (const auto & sp : s_out_extended )
      std::cout << sp<< '\n';

    State2_Action s2{s_start};
    s2.J() = 0;

    auto times = PanosUtilities::linspace(0.0, s_out_extended[0].t(), 100);

    std::cout << "times:\n";
    for (const auto& t:times)
      std::cout << t << '\n';


    auto orbit_range = make_interval_range(Dynamics::DynamicSystem(ham), s2, times, options);

    std::cout << "orbit range:\n";
    for (const auto& s_t : orbit_range)
      {
        std::cout <<  s_t.second << "  "<< s_t.first  << '\n';
      }

    return 0;
  }

}