//
// Created by Panagiotis Zestanakis on 03/10/18.
//

#ifndef HAMILTONIANS_HAMILTONIAN_HPP
#define HAMILTONIANS_HAMILTONIAN_HPP

#include "State.hpp"
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/math/constants/constants.hpp>


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

        class PendulumHamiltonian
        {
          double F_; // = m * g * h
          double G_; // = 1/( m * h^2)
          inline double two_kappa_square (double energy) const
          {
            return 1 + energy/F_;
          };

          inline double K( double kappa) const
          {
            return boost::math::ellint_1( kappa);
          };

          inline double Epsilon(double kappa) const
          {
            return boost::math::ellint_2(kappa);
          };


          inline double R_times_eight_over_pi () const
          {
            constexpr double eight_over_pi = 16 * boost::math::double_constants::one_div_two_pi;

            return eight_over_pi * std::sqrt(F_/G_);
          }

         public:
          PendulumHamiltonian (double F, double G);
          double F () const noexcept ;
          double G () const noexcept ;

          double value(const Geometry::State2 & s) const noexcept ;

          Geometry::State2 derivative(const Geometry::State2& s) const noexcept;

          double analytical_action(double energy) const
          {

            const auto kappa_square = 0.5*two_kappa_square(energy);
            const auto kappa = std::sqrt(kappa_square);

            if (kappa < 1 ) // not <= , because of ill-defined K(1)
                return R_times_eight_over_pi()* (Epsilon(kappa) - (1-kappa_square) * K(kappa));
            else
              return R_times_eight_over_pi() * 0.5 * kappa * Epsilon(1/kappa);


          }

          double analytical_action(const Geometry::State2& s) const
          {
            return analytical_action(value(s));
          }
        };
    }
}

#endif //HAMILTONIANS_HAMILTONIAN_HPP
