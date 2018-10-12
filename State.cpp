//
// Created by Panagiotis Zestanakis on 03/10/18.
//

#include "State.hpp"
#include <boost/range/algorithm/transform.hpp>
#include <boost/math/special_functions/pow.hpp>

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
          return *this;
        }
        State2::State2 (const State3& s3) noexcept
        :v_{{s3.q(),s3.p()}}
        {

        }
        State2& State2::operator-= (const State2& other) noexcept
        {
          v_[0]-= other.v_[0];
          v_[1]-= other.v_[1];
          return *this;
        }

        State3::State3 (double q, double p, double t) noexcept
            : v_{{q, p, t}}
        {

        }
        double State3::q () const noexcept
        {
          return v_[0];
        }
        double State3::p () const noexcept
        {
          return v_[1];
        }
        double State3::t () const noexcept
        {
          return v_[2];
        }
        State3& State3::operator+= (double d) noexcept
        {
          boost::range::transform(v_,
                                  std::begin(v_),
                                  [d] (auto x)
                                  { return x + d; });
          return *this;
        }
        State3& State3::operator*= (double d) noexcept
        {
          boost::range::transform(v_,
                                  std::begin(v_),
                                  [d] (auto x)
                                  { return x * d; });
          return *this;
        }
        State3& State3::operator/= (double d)
        {
          boost::range::transform(v_,
                                  std::begin(v_),
                                  [d] (auto x)
                                  { return x / d; });
          return *this;
        }
        State3& State3::operator+= (const State3& other) noexcept
        {
          boost::range::transform(v_, other.v_,
                                  std::begin(v_),
                                  [] (auto x1, auto x2)
                                  { return x1 + x2; });
          return *this;
        }
        State3::State3 (const State2& s2, double t) noexcept
        :v_{{s2.q(),s2.p(),t}}
        {

        }
        State3& State3::operator-= (const State3& other) noexcept
        {
          v_[0]-= other.v_[0];
          v_[1]-= other.v_[1];
          v_[2]-= other.v_[2];
          return *this;
        }

        State3 operator- (const State3& other)
        {
          return State3(-other.q(), -other.p(), -other.t());
        }
        double operator* (const State3& s1, const State3& s2)
        {
          return s1.q() * s2.q()
                 + s1.p() * s2.p()
                 + s1.t() * s2.t();
        }
        double magnitude_squared (const State3& s)
        {
          using boost::math::pow;
          return pow<2>(s.q()) + pow<2>(s.p()) + pow<2>(s.t());
        }
        State3 abs (const State3& s)
        {
          using std::abs;

          return State3(abs(s.q()), abs(s.p()), abs(s.t()));
        }
        std::ostream& operator<< (std::ostream& out, const State3& s)
        {
          out << s.q() << ' ' << s.p() << ' ' << s.t();
          return out;
        }
    }
}
