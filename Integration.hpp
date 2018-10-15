//
// Created by Panagiotis Zestanakis on 08/10/18.
//

#ifndef HAMILTONIANS_INTEGRATION_HPP
#define HAMILTONIANS_INTEGRATION_HPP

#include "line.hpp"
#include "Hamiltonian.hpp"
#include "dynamic_system.hpp"
#include <boost/numeric/odeint.hpp>


namespace Integrators
{



//    template <typename Ham>
//    Geometry::Line make_perpendicular_line_to_motion(const Geometry::State2& initial_position)
//    {
//      auto initial_derivative = Geometry::State2{};
//      Dynamics::dynamic_system<Ham>(initial_position, initial_derivative, 0/*t*/);
//
//      return Geometry::Line(initial_position,initial_derivative);
//    }

//    template <typename Ham>
//    class CrossLinePussBackFunctor
//    {
//
//     public:
//      CrossLinePussBackFunctor ()
//      { }
//    };

}
#endif //HAMILTONIANS_INTEGRATION_HPP
