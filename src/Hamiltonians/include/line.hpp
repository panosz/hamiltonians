//
// Created by Panagiotis Zestanakis on 03/10/18.
//

#ifndef HAMILTONIANS_LINE_HPP
#define HAMILTONIANS_LINE_HPP
#include <iostream>
#include "State.hpp"

namespace Integrators
{
    namespace Geometry
    {






        class Line {
          friend std::ostream& operator<< (std::ostream&, const Line&);
          State2 s_{};
          double c_{0};
         public:
          double operator() (const State2& position) const noexcept;
          Line (const State2& position,State2 perpendicular_vector);
          State2 perpendicular_vector() const noexcept;
        };

        std::ostream& operator<< (std::ostream&, const Line&);

        class LineCrossObserver
        {
          Line line_;
          mutable double distance_=0;
         public:
          explicit LineCrossObserver (Line line) noexcept ;
          bool operator()(const State2& next_point) const noexcept ;
          double distance() const noexcept ;
        };
    }

}




#endif //HAMILTONIANS_LINE_HPP
