//
// Created by Panagiotis Zestanakis on 03/10/18.
//
#include <boost/math/special_functions/pow.hpp>
#include <cmath>
#include "Hamiltonian.hpp"
namespace Integrators
{
    namespace Hamiltonian
    {

        double HarmonicOscillator::value (const Geometry::State2& s) const noexcept
        {
          return 0.5 * magnitude_squared(s);
        }
        Geometry::State2 HarmonicOscillator::derivative (const Geometry::State2& s) const noexcept
        {
          return s;
        }

        DuffingHamiltonian::DuffingHamiltonian (double omega, double omega0, double e_alpha, double e_gamma)
            : omega_(omega), omega0_(omega0), e_alpha_(e_alpha), e_gamma_(e_gamma)
        { }
        double DuffingHamiltonian::get_omega () const noexcept
        {
          return omega_;
        }
        void DuffingHamiltonian::set_omega (double omega) noexcept
        {
          omega_ = omega;
        }
        double DuffingHamiltonian::get_omega0 () const noexcept
        {
          return omega0_;
        }
        void DuffingHamiltonian::set_omega0 (double omega0)noexcept
        {
          omega0_ = omega0;
        }
        double DuffingHamiltonian::get_e_alpha () const noexcept
        {
          return e_alpha_;
        }
        void DuffingHamiltonian::set_e_alpha (double e_alpha) noexcept
        {
          e_alpha_ = e_alpha;
        }
        double DuffingHamiltonian::get_e_gamma () const noexcept
        {
          return e_gamma_;
        }
        void DuffingHamiltonian::set_e_gamma (double e_gamma) noexcept
        {
          e_gamma_ = e_gamma;
        }
        double DuffingHamiltonian::value (const Geometry::State2& s) const noexcept
        {
          using boost::math::pow;

          const double q = s.q();

          const double hypot_sq = magnitude_squared(s);

          return -(e_Omega() * hypot_sq
                   + 3 * e_alpha_ / 8 * pow<2>(hypot_sq)
                   - 2 * e_gamma_ * q) / (4 * omega_);
        }
        Geometry::State2 DuffingHamiltonian::derivative (const Geometry::State2& s) const noexcept
        {
          using boost::math::pow;

          const double q = s.q();
          const double p = s.p();

          const double hypot_sq = magnitude_squared(s);

          const double dHdq = -(e_Omega() * q + 3 * e_alpha_ / 4 * hypot_sq * q - e_gamma_) / (2 * omega_);
          const double dHdp = -(e_Omega() * p + 3 * e_alpha_ / 4 * hypot_sq * p) / (2 * omega_);

          return Geometry::State2{dHdq,dHdp};

        }

        PendulumHamiltonian::PendulumHamiltonian (double FF, double GG)
            : F_(FF), G_(GG)
        { }
        double PendulumHamiltonian::F () const noexcept
        {
          return F_;
        }
        double PendulumHamiltonian::G () const noexcept
        {
          return G_;
        }


        double PendulumHamiltonian::value (const Geometry::State2& s) const noexcept
        {
          using boost::math::pow;
          using std::cos;

          const auto q = s.q();
          const auto p = s.p();

          return 0.5*G_*pow<2>(p) - F_*cos(q);
        }
        Geometry::State2 PendulumHamiltonian::derivative (const Geometry::State2& s) const noexcept
        {
          using std::sin;

          const auto q = s.q();
          const auto p = s.p();

          return Geometry::State2{F_*sin(q),p};
        }

        double FreeParticle::value (const Geometry::State2& s) const noexcept
        {
          using boost::math::pow;
          return 0.5*pow<2>(s.p());
        }
        Geometry::State2 FreeParticle::derivative (const Geometry::State2& s) const noexcept
        {
          return Integrators::Geometry::State2{0,s.p()};
        }

    }
}