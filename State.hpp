//
// Created by Panagiotis Zestanakis on 03/10/18.
//

#ifndef HAMILTONIANS_STATE_HPP
#define HAMILTONIANS_STATE_HPP

#include <iostream>
#include <initializer_list>
#include <type_traits>
#include <armadillo>

#include <boost/operators.hpp>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>

namespace Integrators
{
    namespace Geometry
    {

        template<unsigned N, typename = typename std::enable_if<(N >= 2)>::type>
        class State : boost::additive<State<N>, boost::additive<State<N>, double,
            boost::multiplicative<State<N>, double> > > {

          template<unsigned M, typename>
          friend
          class State;

          using vector_type = arma::colvec::fixed<N>;

          template<bool B, typename T>
          using Enable_if = typename std::enable_if<B, T>::type;

          vector_type v_{arma::fill::zeros};

         public:
          State () = default;
          State (std::initializer_list<double> l) noexcept
              : v_{l}
          { };

          /// \brief constructor from State with more elements.
          /// \tparam M the dimension of the other state
          /// \param other the other state
          ///
          /// see https://stackoverflow.com/a/17842695/6060982
          template<unsigned M>
          explicit State (const State<M>& other, typename std::enable_if<(N < M)>::type * = 0) noexcept
              : v_{other.v_(arma::span(0, N - 1))}
          {
          }

          /// \brief constructor from State with fewer elements.
          /// \tparam M the dimension of the other state
          /// \param other the other state
          ///
          /// see https://stackoverflow.com/a/17842695/6060982
          template<unsigned M>
          explicit State (const State<M>& other, typename std::enable_if<(N > M)>::type * = 0) noexcept
              : v_{}
          {
            v_(arma::span(0, M - 1)) = other.v_;
          }

          template<unsigned M>
          State& operator= (const State<M>& other)
          {

            if constexpr (N < M)
              v_ = other.v_(arma::span(0, N - 1));
            else
              {
                v_(arma::span(0, M - 1)) = other.v_;
                v_(arma::span(M, N - 1)).zeros();
              }

            return *this;
          }

          double q () const noexcept
          {
            return v_[0];
          }

          double& q () noexcept
          {
            return v_[0];
          }

          double p () const noexcept
          {
            return v_[1];
          }

          double& p () noexcept
          { return v_[1]; }

          template<unsigned DIM = N>
          Enable_if<(DIM >= 3), double> J () const noexcept
          { return v_[2]; }

          template<unsigned DIM = N>
          Enable_if<(DIM >= 3), double>& J () noexcept
          { return v_[2]; }

          template<unsigned DIM = N>
          Enable_if<(DIM >= 4), double> t () const noexcept
          { return v_[3]; }

          template<unsigned DIM = N>
          Enable_if<(DIM >= 4), double>& t ()  noexcept
          { return v_[3]; }

          inline State& operator+= (double d) noexcept
          {
            v_ += d;
            return *this;
          }

          inline State& operator*= (double d) noexcept
          {
            v_ *= d;
            return *this;

          }

          inline State& operator/= (double d)
          {
            v_ /= d;
            return *this;

          }

          inline State& operator+= (const State& other) noexcept
          {
            v_ += other.v_;
            return *this;

          }

          inline State& operator-= (const State& other) noexcept
          {
            v_ -= other.v_;
            return *this;

          }

          inline double operator* (const State& other) const noexcept
          {
            return arma::as_scalar(v_.t() * other.v_);
          }

          inline State operator/ (const State& other) const noexcept
          {
            State ret{*this};
            ret.v_ /= other.v_;
            return ret;
          }

          inline double inf_norm () const noexcept
          {
            return arma::norm(v_, "inf");
          }

          State abs () const noexcept
          {
            State ret{};
            ret.v_ = arma::abs(v_);
            return ret;
          }

          template<unsigned DIM>
          friend std::ostream& operator<< (std::ostream& out, const State<DIM>& s);

        };

        template<unsigned N>
        auto abs (const State<N>& s)
        {
          return s.abs();
        }

        template<unsigned N>
        double magnitude_squared (const State<N>& s) noexcept
        {
          return s * s;
        }

        template<unsigned N>
        std::ostream& operator<< (std::ostream& out, const State<N>& s)
        {
          for (const auto& coord : s.v_)
            out << coord << ' ';
          return out;
        }

        using State2 = State<2>;
        using State2_Action = State<3>;
        using State2_Extended = State<4>;
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

            template<>
            struct vector_space_norm_inf<Integrators::Geometry::State2_Extended> {
                typedef double result_type;
                double operator() (const Integrators::Geometry::State2_Extended& s) const
                {
                  return s.inf_norm();
                }
            };
        }
    }
}
#endif //HAMILTONIANS_STATE_HPP
