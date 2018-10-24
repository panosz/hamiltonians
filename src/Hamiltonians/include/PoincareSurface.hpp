//
// Created by Panagiotis Zestanakis on 12/10/18.
//

#ifndef HAMILTONIANS_POINCARESURFACE_HPP
#define HAMILTONIANS_POINCARESURFACE_HPP

#include "State.hpp"
#include "line.hpp"
#include "dynamic_system.hpp"

namespace Integrators
{

    template<typename Ham, typename Stepper>
    class PoincareSurface {
     private:
      Ham ham_;
      Geometry::Line cross_line_;
      Geometry::State2 direction_;
      mutable Stepper stepper_;

     public:

      std::pair<Geometry::State2_Action, double> step_back (Geometry::State2_Action s, double t, double current_distance) const
      {

        auto state_extended = Geometry::State2_Extended{s};

        state_extended.t()=t;

        auto df = [this] (const Geometry::State2_Extended& s_extended, Geometry::State2_Extended& dsdt_extended, double /*time*/)
        {

            dsdt_extended = Integrators::Dynamics::dynamic_system_along_direction_impl(ham_,
                                                                                       direction_,
                                                                                       Geometry::State2_Action{
                                                                                           s_extended});
        };

        stepper_.do_step(
            df,
            state_extended,
            t,
            -current_distance);

        s = Geometry::State2_Action{state_extended};
        t = state_extended.t();

        return std::make_pair(s, t);
      };

      Geometry::Line cross_line () const
      { return cross_line_; };

      PoincareSurface (Ham&& ham, Geometry::Line cross_line, Stepper&& stepper)
          : ham_{std::forward<Ham>(ham)},
            cross_line_{cross_line},
            direction_{cross_line.perpendicular_vector()},
            stepper_{std::forward<Stepper>(stepper)}
      { };
    };

    template<typename Ham, typename Stepper>
    inline auto make_PoincareSurface (Ham&& hamiltonian, Geometry::State2 s_start, Stepper&& stepper)
    {
      const auto start_direction = Integrators::Dynamics::dynamic_system_impl(hamiltonian, s_start);

      const auto start_line = Integrators::Geometry::Line(s_start, start_direction);
      return PoincareSurface<Ham, Stepper>(std::forward<Ham>(hamiltonian), start_line, std::forward<Stepper>(stepper));

    }
}

#endif //HAMILTONIANS_POINCARESURFACE_HPP
