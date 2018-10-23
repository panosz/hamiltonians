//
// Created by Panagiotis Zestanakis on 23/10/18.
//

#ifndef HAMILTONIANS_ACTION_ANGLE_HPP
#define HAMILTONIANS_ACTION_ANGLE_HPP

#include <boost/math/constants/constants.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm_ext/push_back.hpp>



#include "State.hpp"
#include "Integration.hpp"
#include "myUtilities.hpp"


namespace Integrators
{

    class ActionAngleOrbit {
      double action_two_pi_;
      double omega_;
      std::vector<double> theta_;
      std::vector<Geometry::State2> positions_;

     public:
      ActionAngleOrbit (double action_two_pi,
                        double omega,
                        const std::vector<double>& theta,
                        const std::vector<Geometry::State2>& positions);
      double action_two_pi () const;
      double omega () const;
      const std::vector<double>& theta () const;
      const std::vector<Geometry::State2>& positions () const
      {
        return positions_;
      }
    };

    template<typename Ham>
    ActionAngleOrbit calculate_action_angle_on_closed_orbit (Ham hamiltonian,
                                                             const Geometry::State2& s_start,
                                                             IntegrationOptions& options,
                                                             size_t number_of_angles = 100)
    {
      Geometry::State2_Action s_start_Action{s_start};

      const auto s_out_extended = come_back_home(hamiltonian, s_start, options);

      const auto action = s_out_extended.J();
      const auto omega = boost::math::double_constants::two_pi / s_out_extended.t();

      Geometry::State2 s2 = s_start;

      auto times = PanosUtilities::linspace(0.0, s_out_extended.t(), number_of_angles);

      auto orbit_range = make_interval_range(Dynamics::DynamicSystem(hamiltonian), s2, times, options);

      std::vector<Geometry::State2> positions{};
      boost::push_back(positions, orbit_range | boost::adaptors::transformed([] (const auto& p)
                                                                             { return p.first; }));

      const auto theta = PanosUtilities::linspace(0.0, boost::math::double_constants::two_pi, number_of_angles);

      return ActionAngleOrbit{action, omega, theta, std::move(positions)};

    }

}
#endif //HAMILTONIANS_ACTION_ANGLE_HPP
