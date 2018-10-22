#include <iostream>
#include <stdexcept>

#include "line.hpp"
#include "Hamiltonian.hpp"
#include "dynamic_system.hpp"
#include "observer.hpp"
#include "Integration.hpp"
#include "PoincareSurface.hpp"

#include <boost/math/constants/constants.hpp>

#include <boost/range/combine.hpp>
#include <boost/range/algorithm/find_if.hpp>
#include <boost/range/algorithm_ext/push_back.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include "myUtilities.hpp"

using namespace Integrators;
using namespace Integrators::Geometry;
using namespace boost::numeric::odeint;

class ActionAngleOrbit {
  double action_two_pi_;
  double omega_;
  std::vector<double> theta_;
  std::vector<Geometry::State2> positions_;

 public:
  ActionAngleOrbit (double action_two_pi, double omega, const std::vector<double>& theta, const std::vector<State2>& positions)
      : action_two_pi_(action_two_pi), omega_(omega), theta_(theta), positions_(positions)
  {
  }
  double action_two_pi () const
  {
    return action_two_pi_;
  }
  double omega () const
  {
    return omega_;
  }
  const std::vector<double>& theta () const
  {
    return theta_;
  }
  const std::vector<State2>& positions () const
  {
    return positions_;
  }
};

template<typename Ham>
ActionAngleOrbit calculate_action_angle_on_closed_orbit(Ham hamiltonian,
                                                        const State2& s_start,
                                                       IntegrationOptions& options,
                                                       size_t number_of_angles =100)
{
  Geometry::State2_Action s_start_Action{s_start};

  const auto s_out_extended = come_back_home(hamiltonian, s_start, options);

  const auto action = s_out_extended.J();
  const auto omega = boost::math::double_constants::two_pi / s_out_extended.t();


  State2 s2 = s_start;

  auto times = PanosUtilities::linspace(0.0, s_out_extended.t(), number_of_angles);

  auto orbit_range = make_interval_range(Dynamics::DynamicSystem(hamiltonian), s2, times, options);

  std::vector<State2> positions{};
  boost::push_back(positions, orbit_range|boost::adaptors::transformed([](const auto & p){return p.first;}));

  const auto theta = PanosUtilities::linspace(0.0, boost::math::double_constants::two_pi,number_of_angles);

  return ActionAngleOrbit{action, omega, theta, std::move(positions)};

}

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