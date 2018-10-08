//
// Created by Panagiotis Zestanakis on 03/10/18.
//

#ifndef HAMILTONIANS_LINE_HPP
#define HAMILTONIANS_LINE_HPP
#include <iostream>
#include "State2.hpp"

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

        };

        std::ostream& operator<< (std::ostream&, const Line&);
    }

}




#endif //HAMILTONIANS_LINE_HPP
