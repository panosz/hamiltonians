//
// Created by Panagiotis Zestanakis on 03/10/18.
//

#include "line.hpp"

#include <stdexcept>
#include <cmath>

namespace Integrators
{
    namespace Geometry
    {
        Line::Line (const State2& position,State2 perpendicular_vector)
            : s_{perpendicular_vector}, c_{-(position*perpendicular_vector)}
        {
          if (magnitude_squared(s_)==0)
            throw std::invalid_argument("Line: at least one of a and b must be nonzero");
        }

        double Line::operator() (const State2& position) const noexcept
        {
          return s_*position + c_;
        }

        std::ostream& operator<< (std::ostream& os, const Line& line)
        {
          const auto a = line.s_.q();
          const auto b = line.s_.p();
          const auto c = line.c_;

          os << a << "*x"
             << ((b >= 0) ? " + " : " - ")
             << std::abs(b) << "*y"
             << ((c>= 0) ? " + " : " - ")
             << std::abs(c) << " = 0";

          return os;
        }
        State2 Line::perpendicular_vector () const noexcept
        {
          return s_;
        }

        LineCrossObserver::LineCrossObserver (Line line) noexcept
            : line_(std::move(line))
        { }


        inline bool crossZeroPositiveDirectionPredicate (double current_value, double previous_value)
        {
          return (current_value >= 0 && previous_value< 0);
        }


        bool LineCrossObserver::operator() (const State2& next_point) const noexcept
        {
          const double next_distance = line_(next_point);
          const bool crossed_line = crossZeroPositiveDirectionPredicate(next_distance,distance_);
          distance_ = next_distance;
          return crossed_line;
        }
        double LineCrossObserver::distance () const noexcept
        {
          return distance_;
        }

    }
}