//
// Created by Panagiotis Zestanakis on 25/10/18.
//

#include <iostream>
#include <boost/range/combine.hpp>
#include <boost/math/constants/constants.hpp>

#include "periodic_q_surface.hpp"
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
  const State2 s_start{0, 0.5};

  const auto periodicQsurfObserver = PeriodicQSurfaceCrossObserver{s_start};

  const auto hamiltonian = Integrators::Hamiltonian::FreeParticle();

  const auto vec = PanosUtilities::linspace(0,30,200);


  IntegrationOptions options;

  TimeInterval t_interval{0,100,0.01};


  const auto crossings =
      calculate_crossings(hamiltonian, s_start, periodicQsurfObserver, t_interval, options);

  std::cout << " crossings : \n";
  for (const auto& c:crossings)
    std::cout << c << '\n';

  return 0;
}