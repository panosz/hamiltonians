//
// Created by Panagiotis Zestanakis on 23/10/18.
//

#include "action_angle.hpp"
namespace Integrators
{

    ActionAngleOrbit::ActionAngleOrbit (double action_two_pi, double omega, const std::vector<double>& theta, const std::vector<Geometry::State2>& positions)
        : action_two_pi_(action_two_pi), omega_(omega), theta_(theta), positions_(positions)
    {
    }
    double ActionAngleOrbit::action_two_pi () const
    {
      return action_two_pi_;
    }
    double ActionAngleOrbit::omega () const
    {
      return omega_;
    }
    const std::vector<double>& ActionAngleOrbit::theta () const
    {
      return theta_;
    }
}