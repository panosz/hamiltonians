#include <iostream>
#include <boost/range/combine.hpp>


#include "line.hpp"
#include "Hamiltonian.hpp"
#include "dynamic_system.hpp"
#include "observer.hpp"
#include "Integration.hpp"
#include "action_angle.hpp"



using namespace Integrators;
using namespace Integrators::Geometry;
using namespace boost::numeric::odeint;



int main ()
{

  const State2 s_start{0.472035, 7.86664};

  const Line cross_line{s_start, {-1, 0}};

  const auto hamiltonian = Integrators::Hamiltonian::DuffingHamiltonian();

  IntegrationOptions options;

  options.set_integration_time(100.0);

  const auto integration_end_point_without_closeness_filtering =
      calculate_first_crossing(hamiltonian, s_start, cross_line,options);

  std::cout <<" integration end point without closeness filtering : "
              << integration_end_point_without_closeness_filtering<<'\n';

  const auto actionAngleOrbit = calculate_action_angle_on_closed_orbit(hamiltonian,s_start,options,60);

  std::cout << "orbit range:\n";
  std::cout << "action = "<< actionAngleOrbit.action_two_pi()<<'\n';

  for ( const auto s_t : boost::range::combine(actionAngleOrbit.positions(),actionAngleOrbit.theta()))
    {
      double theta;
      State2 s;
      boost::tie(s,theta) = s_t;
      std::cout << theta/actionAngleOrbit.omega() << "  " << s << '\n';
    }

  return 0;

}