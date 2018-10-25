//
// Created by Panagiotis Zestanakis on 25/10/18.
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
  const State2 s_start{0,0.5};

  const State2 s_cross{1,0};

  const Line cross_line{s_cross, {1, 0}};

  const auto hamiltonian = Integrators::Hamiltonian::FreeParticle();

  IntegrationOptions options;

  options.set_integration_time(100.0);

  const auto integration_end_point_without_closeness_filtering =
      calculate_first_crossing(hamiltonian, s_start, cross_line,options);

  std::cout <<" integration end point without closeness filtering : "
            << integration_end_point_without_closeness_filtering<<'\n';

  return 0;
}