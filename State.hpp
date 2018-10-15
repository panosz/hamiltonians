//
// Created by Panagiotis Zestanakis on 03/10/18.
//

#ifndef HAMILTONIANS_STATE_HPP
#define HAMILTONIANS_STATE_HPP

#include <iostream>
#include <array>
#include <armadillo>

#include <boost/operators.hpp>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>

namespace Integrators
{
    namespace Geometry
    {

        class State2_Extended;
        class State2_Action;

        class State2 : boost::additive<State2, boost::additive<State2, double,
            boost::multiplicative<State2, double> > > {
          using vector_type = arma::colvec2;
          vector_type v_{};
         public:
          State2 () = default;
          State2 (double q, double p) noexcept;
//           explicit State2 (const State2_Extended& s_e) noexcept;
          explicit State2 (const State2_Action& s_a) noexcept;

          double q () const noexcept;
          double p () const noexcept;

          State2& operator+= (double d) noexcept;

          State2& operator*= (double d) noexcept;

          State2& operator/= (double d);

          State2& operator+= (const State2& other) noexcept;

          State2& operator-= (const State2& other) noexcept;

          double operator* (const State2& other) const noexcept;

          State2 operator/ (const State2& other) const noexcept;

          inline double inf_norm () const noexcept
          {
            return arma::norm(v_, "inf");
          }

        };

        State2 operator- (const State2& other);

        double magnitude_squared (const State2& s);

        State2 abs (const State2& s);

        std::ostream& operator<< (std::ostream& out, const State2& s);

        // State2_Action
        class State2_Action : boost::additive<State2_Action, boost::additive<State2_Action, double,
            boost::multiplicative<State2_Action, double> > > {
          using vector_type = arma::colvec3;
          vector_type v_{};
         public:
          State2_Action () = default;
          State2_Action (double q, double p, double J = 0) noexcept
              : v_{q, p, J}
          {

          };

          State2_Action (const State2& s2, double J = 0) noexcept
              : v_{s2.q(), s2.p(), J}
          { };

          explicit State2_Action(const State2_Extended & s_e) noexcept;


          inline double q () const noexcept
          {
            return v_[0];
          }

          inline double p () const noexcept
          {
            return v_[1];
          }

          inline double J () const noexcept
          {
            return v_[2];
          }

          inline State2_Action& operator+= (double d) noexcept
          {
            v_ += d;
            return *this;
          }

          inline State2_Action& operator*= (double d) noexcept
          {
            v_ *= d;
            return *this;

          }

          inline State2_Action& operator/= (double d)
          {
            v_ /= d;
            return *this;

          }

          inline State2_Action& operator+= (const State2_Action& other) noexcept
          {
            v_ += other.v_;
            return *this;

          }

          inline State2_Action& operator-= (const State2_Action& other) noexcept
          {
            v_ -= other.v_;
            return *this;

          }

          inline double operator* (const State2_Action& other) const noexcept
          {
            return arma::as_scalar(v_.t() * other.v_);
          }

          inline State2_Action operator/ (const State2_Action& other) const noexcept
          {
            State2_Action ret{*this};
            ret.v_ /= other.v_;
            return ret;
          }

          inline double inf_norm () const noexcept
          {
            return arma::norm(v_, "inf");
          }

          friend State2_Action abs (const State2_Action& s);
          friend std::ostream& operator<< (std::ostream& out, const State2_Action& s);
        };

        inline State2_Action operator- (const State2_Action& other)
        {
          return State2_Action{-other.q(), -other.q(), -other.J()};
        }
//
//        double magnitude_squared (const State2_Action& s);
//
        State2_Action abs (const State2_Action& s);
//
        std::ostream& operator<< (std::ostream& out, const State2_Action& s);



        // State2_Action end





        class State2_Extended : boost::additive<State2_Extended, boost::additive<State2_Extended, double,
            boost::multiplicative<State2_Extended, double> > > {
          using vector_type = arma::colvec4;
          vector_type v_{};
         public:
          State2_Extended () = default;
          State2_Extended (double q, double p, double J = 0, double t = 0) noexcept;
          State2_Extended (const State2_Action& s2_action, double t = 0) noexcept;

          double q () const noexcept;
          double p () const noexcept;
          double J () const noexcept;
          double t () const noexcept;

          State2_Extended& operator+= (double d) noexcept;

          State2_Extended& operator*= (double d) noexcept;

          State2_Extended& operator/= (double d);

          State2_Extended& operator+= (const State2_Extended& other) noexcept;

          State2_Extended& operator-= (const State2_Extended& other) noexcept;

          double operator* (const State2_Extended& other) const noexcept;

          State2_Extended operator/ (const State2_Extended& other) const noexcept;

          inline double inf_norm () const noexcept
          {
            return arma::norm(v_, "inf");
          }

          friend State2_Extended abs (const State2_Extended& s);
          friend State2_Extended operator- (const State2_Extended& other);
          friend double magnitude_squared (const State2_Extended& s);
          friend std::ostream& operator<< (std::ostream& out, const State2_Extended& s);

        };

        State2_Extended operator- (const State2_Extended& other);

        double magnitude_squared (const State2_Extended& s);

        State2_Extended abs (const State2_Extended& s);

        std::ostream& operator<< (std::ostream& out, const State2_Extended& s);
    }
}

namespace boost
{
    namespace numeric
    {
        namespace odeint
        {
            template<>
            struct vector_space_norm_inf<Integrators::Geometry::State2> {
                typedef double result_type;
                double operator() (const Integrators::Geometry::State2& s) const
                {
                  return s.inf_norm();
                }
            };

            template<>
            struct vector_space_norm_inf<Integrators::Geometry::State2_Action> {
                typedef double result_type;
                double operator() (const Integrators::Geometry::State2_Action& s) const
                {
                  return s.inf_norm();
                }
            };

//            template<>
//            struct vector_space_norm_inf<Integrators::Geometry::State2_Extended> {
//                typedef double result_type;
//                double operator() (const Integrators::Geometry::State2_Extended& s) const
//                {
//                  return s.inf_norm();
//                }
//            };
        }
    }
}
#endif //HAMILTONIANS_STATE_HPP
