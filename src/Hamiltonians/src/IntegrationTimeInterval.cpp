//
// Created by Panagiotis Zestanakis on 26/10/18.
//

#include "IntegrationTimeInterval.hpp"

namespace Integrators
{

    TimeInterval::TimeInterval (double t_begin, double t_end) noexcept
        : t_begin_{t_begin}, t_end_{t_end}
    { }
    TimeInterval::TimeInterval (double t_begin, double t_end, double dt_max) noexcept
        : t_begin_{t_begin}, t_end_{t_end}, dt_max_{std::optional{dt_max}}
    { }
    double TimeInterval::t_begin () const noexcept
    {
      return t_begin_;
    }
    void TimeInterval::set_t_begin (double t_begin) noexcept
    {
      t_begin_ = t_begin;
    }
    double TimeInterval::t_end () const noexcept
    {
      return t_end_;
    }
    void TimeInterval::set_t_end (double t_end) noexcept
    {
      t_end_ = t_end;
    }
    const std::optional<double>& TimeInterval::dt_max () const noexcept
    {
      return dt_max_;
    }
    void TimeInterval::set_dt_max (double dt_max) noexcept
    {
      dt_max_ = dt_max;
    }
    void TimeInterval::unset_dt_max () noexcept
    {
      dt_max_.reset();

    }
}