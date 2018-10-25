//
// Created by Panagiotis Zestanakis on 23/10/18.
//

#include <iostream>
#include <boost/range/combine.hpp>
#include <boost/math/constants/constants.hpp>

#include "line.hpp"
#include "Hamiltonian.hpp"
#include "dynamic_system.hpp"
#include "observer.hpp"
#include "Integration.hpp"
#include "action_angle.hpp"

#include "myUtilities/linspace.hpp"

using namespace Integrators;
using namespace Integrators::Geometry;
using namespace boost::numeric::odeint;

int main ()
{

  const State2 s_start{0.1, 0.8};

  const auto hamiltonian = Hamiltonian::PendulumHamiltonian{1, 1};

  auto times = PanosUtilities::linspace(0.0, 10, 100);

  IntegrationOptions options;

  options.set_integration_time(100.0);

  const auto analytical_action = hamiltonian.analytical_action(s_start);

  State2 s1= s_start;


  auto orbit_range = make_interval_range(Dynamics::DynamicSystem(hamiltonian), s1, times, options);

  for (const auto& s : orbit_range)
    std::cout << s.first << '\n';



  const auto actionAngleOrbit = calculate_action_angle_on_closed_orbit(hamiltonian,s_start,options,60);

  std::cout << "orbit range:\n";
  std::cout << "analytical_action = "<< analytical_action<<'\n';
  std::cout << "action = "<< (actionAngleOrbit.action_two_pi()*boost::math::double_constants::one_div_two_pi) <<'\n';

  for ( const auto s_t : boost::range::combine(actionAngleOrbit.positions(),actionAngleOrbit.theta()))
    {
      double theta;
      State2 s;
      boost::tie(s,theta) = s_t;
      std::cout << theta/actionAngleOrbit.omega() << "  " << s << '\n';
    }

  return 0;

}