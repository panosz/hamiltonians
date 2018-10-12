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

        class HarmonicOscillator {
         public:
          double value (const Geometry::State2& s) const noexcept;
          Geometry::State2 derivative (const Geometry::State2& s) const noexcept;
        };

        class DuffingHamiltonian {
         private:
          double omega_ = 1.5;
          double omega0_ = 1;
          double e_alpha_ = 0.05;
          double e_gamma_ = 2.5;
          constexpr double e_Omega() const noexcept
          {return omega0_*omega0_ - omega_*omega_;};
         public:
          DuffingHamiltonian()= default;
          DuffingHamiltonian (double omega, double omega0, double e_alpha, double e_gamma);
          double get_omega () const noexcept ;
          void set_omega (double omega) noexcept ;
          double get_omega0 () const noexcept ;
          void set_omega0 (double omega0) noexcept ;
          double get_e_alpha () const noexcept ;
          void set_e_alpha (double e_alpha) noexcept ;
          double get_e_gamma () const noexcept ;
          void set_e_gamma (double e_gamma) noexcept ;

          double value(const Geometry::State2& s) const noexcept;

          Geometry::State2 derivative(const Geometry::State2& s) const noexcept;

        };
    }
}

#endif //HAMILTONIANS_HAMILTONIAN_HPP
