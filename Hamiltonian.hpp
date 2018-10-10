//
// Created by Panagiotis Zestanakis on 03/10/18.
//

#ifndef HAMILTONIANS_HAMILTONIAN_HPP
#define HAMILTONIANS_HAMILTONIAN_HPP

#include "State.hpp"
namespace Integrators
{
    namespace Hamiltonian
    {

        class HarmonicOscillator
        {
         public:
          double value(const Geometry::State2& s) const noexcept;
          Geometry::State2 derivative (const Geometry::State2& s) const noexcept;
        };

    }
}

#endif //HAMILTONIANS_HAMILTONIAN_HPP
