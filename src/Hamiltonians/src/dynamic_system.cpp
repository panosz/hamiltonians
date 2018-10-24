//
// Created by Panagiotis Zestanakis on 22/10/18.
//

#include "dynamic_system.hpp"

namespace Integrators
{
    namespace Dynamics
    {
        template class DynamicSystem<Hamiltonian::HarmonicOscillator>;
        template class DynamicSystem<Hamiltonian::DuffingHamiltonian>;
    }
}