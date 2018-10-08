//
// Created by Panagiotis Zestanakis on 03/10/18.
//

#include "State2.hpp"

namespace Integrators
{
    namespace Geometry
    {
        double operator* (const State2& s1, const State2& s2)
        {
          return s1.q() * s2.q() + s1.p() * s2.p();
        }

        State2 operator/ (const State2& s1, const State2& s2)
        {
          return State2{s1.q() / s2.q(), s1.p() / s2.p()};
        }
        State2 abs (const State2& s)
        {
          return State2{std::abs(s.q()), std::abs(s.p())};
        }
        std::ostream& operator<< (std::ostream& out, const State2& s)
        {
          out << s.q() << " " << s.p();
          return out;
        }

        double magnitude_squared (const State2& s)
        {
          return s * s;
        }
        State2 operator- (const State2& other)
        {
          return State2{-other.q(), -other.p()};
        }

        State2::State2 (double q, double p) noexcept
            : v_{{q, p}}
        {
        }
        double State2::q () const noexcept
        {
          return v_[0];
        }

        double State2::p () const noexcept
        {
          return v_[1];
        }
        State2& State2::operator+= (double d) noexcept
        {
          v_[0] += d;
          v_[1] += d;

          return *this;
        }
        State2& State2::operator*= (double d) noexcept
        {
          v_[0] *= d;
          v_[1] *= d;
          return *this;
        }
        State2& State2::operator+= (const State2& other) noexcept
        {
          v_[0] += other.v_[0];
          v_[1] += other.v_[1];
          return *this;
        }
        State2& State2::operator/= (double d)
        {
          v_[0] /= d;
          v_[1] /= d;
          return *this;        }
    }
}
