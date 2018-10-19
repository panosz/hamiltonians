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

class StateNear {
  double distance_;
  bool are_near (const State2& s1, const State2& s2) const noexcept
  {
    return magnitude(s1 - s2) <= distance_;
  }
 public:
  explicit StateNear (double distance)
      : distance_(distance)
  { }

  bool operator() (const State2& s1, const State2& s2) const noexcept
  {
    return are_near(s1, s2);
  }

  bool operator() (const State2_Extended& s1, const State2& s2) const noexcept
  {
    return are_near(State2{s1}, s2);
  }

  bool operator() (const State2& s1, const State2_Extended& s2) const noexcept
  {
    return are_near(s1, State2{s2});
  }

};

//template<typename Ham>
//inline auto
//make_hamiltonian_integration_range (const Ham& hamiltonian, Geometry::State2& s_start, const IntegrationOptions& options)
//{
//
//};

template<typename Ham>
std::vector<State2_Extended>
calculate_crossings (const Ham& hamiltonian, Geometry::State2 s_start, const Geometry::Line& cross_line, const IntegrationOptions& options)
{

  std::vector<Geometry::State2_Extended> s_out;

  Geometry::State2_Action s_start_Action{s_start};

  const auto system = Dynamics::DynamicSystem{hamiltonian};
  const auto integration_range = Integrators::make_dynamic_system_integration_range(system, s_start_Action, options);

  auto observer = Integrators::make_cross_line_observer(system, cross_line, s_out);
  observe(observer, integration_range);

  return s_out;
}

template<typename Ham>
std::vector<State2_Extended>
calculate_first_crossing (const Ham& hamiltonian,
                          const Geometry::State2& s_start,
                          const Geometry::Line& cross_line,
                          const IntegrationOptions& options)
{
  std::vector<Geometry::State2_Extended> s_out{};

  const auto system = Dynamics::DynamicSystem{hamiltonian};
  Geometry::State2_Action s_start_Action{s_start};

  const auto integration_range = Integrators::make_dynamic_system_integration_range(system, s_start_Action, options);

  auto observer = Integrators::make_cross_line_observer(system, cross_line, s_out);
  observe_if(observer, integration_range);

  return s_out;
}


int main ()
{


  //integrate with adaptive step using make_controlled
  const State2 s_start{0.472035, 7.86664};

  const Line cross_line{s_start, {-1, 0}};

  const auto hamiltonian = Integrators::Hamiltonian::DuffingHamiltonian();
  const auto system = Dynamics::DynamicSystem{hamiltonian};

  IntegrationOptions options;

  options.set_integration_time(100.0);

//  Geometry::State2_Action s_start_Action{s_start};
//  const auto integration_range = Integrators::make_dynamic_system_integration_range(system, s_start_Action, options);

  const auto s_out_extended = calculate_first_crossing(hamiltonian, s_start, cross_line, options);
//
  std::cout << "adaptive integration output:\n";

  for (const auto& sp : s_out_extended)
    std::cout << sp << '\n';

  State2_Action s2{s_start};
  s2.J() = 0;

  auto times = PanosUtilities::linspace(0.0, s_out_extended[0].t(), 100);

  auto orbit_range = make_interval_range(Dynamics::DynamicSystem(hamiltonian), s2, times, options);

  std::cout << "orbit range:\n";
  for (const auto& s_t : orbit_range)
    {
      std::cout << s_t.second << "  " << s_t.first << '\n';
    }

  return 0;

}