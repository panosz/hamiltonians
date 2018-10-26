//
// Created by Panagiotis Zestanakis on 24/10/18.
//
#include <iostream>
#include <vector>
#include <algorithm>

#include <boost/range/adaptor/transformed.hpp>
#include <periodic_q_surface.hpp>

#include "details/periodic_q_distance.hpp"
#include "myUtilities/myUtilities.hpp"

int main ()
{

  const auto periodic2 = Integrators::Internals::PeriodicQDistance{0};

  const auto x_range =PanosUtilities::linspace(0, 10, 30);

  for (const auto & x: x_range)
    {
      std::cout << x << ' ' << periodic2.distance(x) << '\n';
    }

  const auto per_surf = Integrators::Geometry::PeriodicQSurfaceCrossObserver{Integrators::Geometry::State2{0,1}};

  const auto first_crossing  =  std::find_if(x_range.begin(),x_range.end(),
                                             [per_surf](auto x){return per_surf(Integrators::Geometry::State2{x,0});});
  std::cout << "zero crossing " << *first_crossing<< '\n';


  return 0;
}