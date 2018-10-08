//
// Created by Panagiotis Zestanakis on 03/10/18.
//

#include "Hamiltonian.hpp"
namespace Integrators
{
    namespace Hamiltonian
    {

        double HarmonicOscillator::value (const Geometry::State2& s) const noexcept
        {
            return 0.5*magnitude_squared(s);
        }
        Geometry::State2 HarmonicOscillator::derivative (const Geometry::State2& s) const noexcept
        {
            return s;
        }

    }
}