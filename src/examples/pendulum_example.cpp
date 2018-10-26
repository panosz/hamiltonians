//
// Created by Panagiotis Zestanakis on 23/10/18.
//

#include <iostream>
#include <boost/range/adaptor/transformed.hpp>
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

struct ActionResult
{
    double analytical_action=0;
    double numerical_action = 0;
    double energy=0;
    double omega=0;
};

ActionResult calculate_action(const State2& s_init)
{
  ActionResult result{};

  const auto hamiltonian = Hamiltonian::PendulumHamiltonian{1, 1};

  result.analytical_action= hamiltonian.analytical_action(s_init);

  result.energy = hamiltonian.value(s_init);


  IntegrationOptions options;
  options.set_distance_threshold(1e-9);

  const auto integrationTime = Integrators::TimeInterval{0.0,1000};


  State2 s1= s_init;
  try
    {
      const auto actionAngleOrbit = calculate_action_angle_on_closed_orbit(hamiltonian, s1, integrationTime, options, 2);

      result.numerical_action = actionAngleOrbit.action_two_pi() * boost::math::double_constants::one_div_two_pi;
      result.omega = actionAngleOrbit.omega();
    }
  catch (std::exception & e)
    {
      std::cout<<e.what()<<'\n';
      std::cout<<"At orbit starting from point "<< s_init<<'\n';
      result.numerical_action=std::nan("");
      result.omega=std::nan("");

    }
  return result;

}

int main ()
{

  const auto x_init_vals = PanosUtilities::linspace(0.01,boost::math::double_constants::pi,30);

  std::vector<ActionResult> action_result_vector;

  boost::push_back(action_result_vector, x_init_vals | boost::adaptors::transformed( [](auto x){return calculate_action(State2{x,0});}));

  std::cout<<"engergy\tanalytical_action\tnumerical_action\tomega\n";
  for (const auto& a_r: action_result_vector)
    std::cout << a_r.energy<<'\t'
        << a_r.analytical_action<<'\t'
        << a_r.numerical_action<<'\t'
        << a_r.omega<<'\n';


  auto times = PanosUtilities::linspace(0.0, 10, 1000);

  auto s_escaping_init = Geometry::State2{1.62979, 0};
  auto orbit_range = make_interval_range(Dynamics::DynamicSystem(Hamiltonian::PendulumHamiltonian{1, 1}),
                                         s_escaping_init, times, Integrators::IntegrationOptions{});

  std::ofstream output("not_retunring.txt");

  output<< "not returning orbit\n";
  for (const auto & s: orbit_range)
    output<<s.first<<'\n';
//  for ( const auto s_t : boost::range::combine(actionAngleOrbit.positions(),actionAngleOrbit.theta()))
//    {
//      double theta;
//      State2 s;
//      boost::tie(s,theta) = s_t;
//      std::cout << theta/actionAngleOrbit.omega() << "  " << s << '\n';
//    }

  return 0;

}