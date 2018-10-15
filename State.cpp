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
          v_ += d;

          return *this;
        }
        State2& State2::operator*= (double d) noexcept
        {
          v_ *= d;
          return *this;
        }
        State2& State2::operator+= (const State2& other) noexcept
        {
          v_ += other.v_;
          return *this;
        }
        State2& State2::operator/= (double d)
        {
          v_ /= d;
          return *this;
        }
//        State2::State2 (const State2_Extended& s_e) noexcept
//            : v_{s_e.q(), s_e.p()}
//        {
//
//        }

        State2::State2 (const State2_Action& s_a) noexcept
        :  v_{s_a.q(), s_a.p()}
        {

        }

        State2& State2::operator-= (const State2& other) noexcept
        {
          v_ -= other.v_;
          return *this;
        }
        double State2::operator* (const State2& other) const noexcept
        {

          return arma::as_scalar(v_.t() * other.v_);

        }
        State2 State2::operator/ (const State2& other) const noexcept
        {
          State2 ret{*this};
          ret.v_/=other.v_;
          return ret;
        }


        State2_Action abs (const State2_Action& s)
        {
          State2_Action ret{s};
          ret.v_=arma::abs(ret.v_);
          return ret;
        }
        std::ostream& operator<< (std::ostream& out, const State2_Action& s)
        {
          out<< s.v_.t();
          return out;
        }


        State2_Extended::State2_Extended (double q, double p, double J, double t) noexcept
            : v_{{q, p, t, J}}
        {

        }
        double State2_Extended::q () const noexcept
        {
          return v_[0];
        }
        double State2_Extended::p () const noexcept
        {
          return v_[1];
        }
        double State2_Extended::J() const noexcept
        {
          return v_[2];
        }

        double State2_Extended::t () const noexcept
        {
          return v_[3];
        }
        State2_Extended& State2_Extended::operator+= (double d) noexcept
        {
          v_+=d;
          return *this;
        }
        State2_Extended& State2_Extended::operator*= (double d) noexcept
        {
          v_*=d;
          return *this;
        }
        State2_Extended& State2_Extended::operator/= (double d)
        {
         v_/=d;
          return *this;
        }
        State2_Extended& State2_Extended::operator+= (const State2_Extended& other) noexcept
        {
          v_+=other.v_;
          return *this;
        }
        State2_Extended::State2_Extended (const State2_Action& s2, double t) noexcept
            : v_{{s2.q(), s2.p(), s2.J(), t}}
        {

        }
        State2_Extended& State2_Extended::operator-= (const State2_Extended& other) noexcept
        {
         v_-=other.v_;
          return *this;
        }
        double State2_Extended::operator* (const State2_Extended& other) const noexcept
        {
          return arma::as_scalar(v_.t() * other.v_);
        }


        State2_Extended State2_Extended::operator/ (const State2_Extended& other) const noexcept
        {
          State2_Extended ret{*this};
          ret.v_/=other.v_;
          return ret;
        }

        State2_Extended operator- (const State2_Extended& other)
        {
          State2_Extended ret{other};
          ret.v_=-ret.v_;
          return ret;
        }



        double magnitude_squared (const State2_Extended& s)
        {
          return s*s;
        }
        State2_Extended abs (const State2_Extended& s)
        {
          State2_Extended ret{s};
          ret.v_=arma::abs(ret.v_);
          return ret;
        }
        std::ostream& operator<< (std::ostream& out, const State2_Extended& s)
        {
          out <<s.v_.t();
          return out;
        }



        State2_Action::State2_Action (const State2_Extended& s_e) noexcept
            :v_{s_e.q(),s_e.p(),s_e.J()}
        {

        }
    }
}
