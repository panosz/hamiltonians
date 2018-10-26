//
// Created by Panagiotis Zestanakis on 26/10/18.
//

#ifndef HAMILTONIANS_INTEGRATIONTIMEINTERVAL_HPP
#define HAMILTONIANS_INTEGRATIONTIMEINTERVAL_HPP
#include <optional>

namespace Integrators
{
    /// \brief To be used as input for the integration functions
    class TimeInterval {
     private:
      double t_begin_ = 0;
      double t_end_ = 0;
      std::optional<double> dt_max_{};
     public:
      TimeInterval (double t_begin, double t_end) noexcept;
      TimeInterval (double t_begin, double t_end, double dt_max)noexcept;
      double t_begin () const noexcept;
      void set_t_begin (double t_begin)noexcept;
      double t_end () const noexcept;
      void set_t_end (double t_end) noexcept;
      const std::optional<double>& dt_max () const noexcept;
      void set_dt_max (double dt_max) noexcept;
      void unset_dt_max() noexcept;

    };
}
#endif //HAMILTONIANS_INTEGRATIONTIMEINTERVAL_HPP
