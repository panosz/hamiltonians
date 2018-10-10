//
// Created by Panagiotis Zestanakis on 03/10/18.
//

#ifndef HAMILTONIANS_STATE_HPP
#define HAMILTONIANS_STATE_HPP

#include <iostream>
#include <array>
#include <boost/operators.hpp>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>

namespace Integrators
{
    namespace Geometry
    {

        class State3;

        class State2 : boost::additive<State2, boost::additive<State2, double,
            boost::multiplicative<State2, double> > > {
          using vector_type = std::array<double, 2>;
          vector_type v_{};
         public:
          State2 () = default;
          State2 (double q, double p) noexcept;
          explicit State2(const State3& s3) noexcept ;

          double q () const noexcept;
          double p () const noexcept;

          State2& operator+= (double d) noexcept;

          State2& operator*= (double d) noexcept;

          State2& operator/= (double d);

          State2& operator+= (const State2& other) noexcept;

        };

        State2 operator- (const State2& other);

        double operator* (const State2& s1, const State2& s2);

        double magnitude_squared (const State2& s);
        State2 operator/ (const State2& s1, const State2& s2);

        State2 abs (const State2& s);

        std::ostream& operator<< (std::ostream& out, const State2& s);

        class State3 : boost::additive<State3, boost::additive<State3, double,
            boost::multiplicative<State3, double> > > {
          using vector_type = std::array<double, 3>;
          vector_type v_{};
         public:
          State3 () = default;
          State3 (double q, double p, double t) noexcept;
          State3 (const State2& s2, double t) noexcept ;

          double q () const noexcept;
          double p () const noexcept;
          double t () const noexcept;

          State3& operator+= (double d) noexcept;

          State3& operator*= (double d) noexcept;

          State3& operator/= (double d);

          State3& operator+= (const State3& other) noexcept;

        };

        State3 operator- (const State3& other);

        double operator* (const State3& s1, const State3& s2);

        double magnitude_squared (const State3& s);
        State3 operator/ (const State3& s1, const State3& s2);

        State3 abs (const State3& s);

        std::ostream& operator<< (std::ostream& out, const State3& s);
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
                  using std::max;
                  using std::abs;
                  return max(abs(s.q()), abs(s.p()));
                }
            };

            template<>
            struct vector_space_norm_inf<Integrators::Geometry::State3> {
                typedef double result_type;
                double operator() (const Integrators::Geometry::State3& s) const
                {
                  using std::max;
                  using std::abs;
                  return max(
                      max(abs(s.q()), abs(s.p())), abs(s.t()));
                }
            };
        }
    }
}
#endif //HAMILTONIANS_STATE_HPP
