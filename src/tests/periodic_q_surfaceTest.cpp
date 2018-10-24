//
// Created by Panagiotis Zestanakis on 24/10/18.
//
#include <iostream>
#include <vector>
#include <boost/range/adaptor/transformed.hpp>

#include "periodic_q_surface.hpp"
#include "myUtilities.hpp"

int main ()
{

  const auto periodic2 = Integrators::PeriodicQSurface{2};

  for (const auto p:PanosUtilities::linspace(0, 10, 2))
    {
      std::cout << p << ' ' << periodic2.distance(p) << '\n';
    }

  std::vector<double> zero_crossings{};

  const auto value_range = PanosUtilities::linspace(0, 10, 30) |
                           boost::adaptors::transformed([per = periodic2] (auto d)
                                                        { return per.distance(d); });
  std::cout << "zero crossings" << '\n';

  std::ostream_iterator<double> out_it(std::cout, "\n");
  Integrators::zero_cross(value_range.begin(), value_range.end(), out_it, 2.0);

  return 0;
}