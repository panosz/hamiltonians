#include <iostream>
#include <stdexcept>

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



int main ()
{


  //integrate with adaptive step using make_controlled
  const State2 s_start{0.472035, 7.86664};

  const Line cross_line{s_start, {-1, 0}};

  const auto hamiltonian = Integrators::Hamiltonian::DuffingHamiltonian();

  IntegrationOptions options;

  options.set_integration_time(100.0);

  Geometry::State2_Action s_start_Action{s_start};


  const auto s_out_extended = come_back_home(hamiltonian, s_start, options);
//
  std::cout << "adaptive integration output:\n";
  std::cout << s_out_extended << '\n';

  State2_Action s2{s_start};
  s2.J() = 0;

  auto times = PanosUtilities::linspace(0.0, s_out_extended.t(), 100);

  auto orbit_range = make_interval_range(Dynamics::DynamicSystem(hamiltonian), s2, times, options);

  std::cout << "orbit range:\n";
  for (const auto& s_t : orbit_range)
    {
      std::cout << s_t.second << "  " << s_t.first << '\n';
    }

  return 0;

}