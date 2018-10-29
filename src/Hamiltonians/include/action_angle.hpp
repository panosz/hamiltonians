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
#include "myUtilities/myUtilities.hpp"

namespace Integrators
{

    struct AnglesPositions {
        std::vector<double> theta{};
        std::vector<Geometry::State2> positions{};
    };

    class ActionAngleOrbit {
      double action_two_pi_{};
      double omega_{};
      std::vector<double> theta_{};
      std::vector<Geometry::State2> positions_{};

     public:
      ActionAngleOrbit () = default;
      ActionAngleOrbit (double action_two_pi,
                        double omega,
                        const std::vector<double>& theta,
                        const std::vector<Geometry::State2>& positions);

      ActionAngleOrbit (double action_two_pi, double omega, AnglesPositions anglesPositions);

      ActionAngleOrbit (double action_two_pi, double omega, AnglesPositions&& anglesPositions);


      double action_two_pi () const;
      double omega () const;
      const std::vector<double>& theta () const;
      const std::vector<Geometry::State2>& positions () const
      {
        return positions_;
      }
    };

    template<typename Ham>
    AnglesPositions map_positions_to_angles_along_orbit (Ham hamiltonian,
                                                         Geometry::State2 s_start,
                                                         double orbit_completion_time,
                                                         const IntegrationOptions& options,
                                                         size_t number_of_angles)
    {
      auto times = PanosUtilities::linspace(0.0, orbit_completion_time, number_of_angles);

      auto orbit_range = make_interval_range(Dynamics::DynamicSystem(hamiltonian), s_start, times, options);

      std::vector<Geometry::State2> positions{};
      boost::push_back(positions, orbit_range | boost::adaptors::transformed([] (const auto& p)
                                                                             { return p.first; }));

      return AnglesPositions{PanosUtilities::linspace(0.0, boost::math::double_constants::two_pi, number_of_angles),
                             positions};
    }

    template<typename Ham>
    ActionAngleOrbit calculate_action_angle_on_closed_orbit (Ham hamiltonian,
                                                             const Geometry::State2& s_start,
                                                             const TimeInterval& integrationTime,
                                                             const IntegrationOptions& options,
                                                             size_t number_of_angles = 100)
    {

      const auto s_out_extended = come_back_home_closed_orbit(hamiltonian, s_start, integrationTime, options);

      const auto action = s_out_extended.J();
      const auto period = s_out_extended.t();
      const auto omega = boost::math::double_constants::two_pi / period;

      const auto anglesPositions = map_positions_to_angles_along_orbit(hamiltonian,
                                                                       s_start,
                                                                       period,
                                                                       options,
                                                                       number_of_angles);


      return ActionAngleOrbit{action, omega, anglesPositions};

    }

    template<typename Ham>
    ActionAngleOrbit calculate_action_angle_on_periodic_orbit (Ham hamiltonian,
                                                               const Geometry::State2& s_start,
                                                               const TimeInterval& integrationTime,
                                                               const IntegrationOptions& options,
                                                               size_t number_of_angles = 100)
    {

      const auto s_out_extended = come_back_home_periodic_orbit(hamiltonian, s_start, integrationTime, options);

      const auto action = s_out_extended.J();
      const auto period = s_out_extended.t();
      const auto omega = boost::math::double_constants::two_pi / period;

      const auto anglesPositions = map_positions_to_angles_along_orbit(hamiltonian,
                                                                       s_start,
                                                                       period,
                                                                       options,
                                                                       number_of_angles);


      return ActionAngleOrbit{action, omega, anglesPositions};

    }

}
#endif //HAMILTONIANS_ACTION_ANGLE_HPP
