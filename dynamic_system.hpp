//
// Created by Panagiotis Zestanakis on 03/10/18.
//

#ifndef HAMILTONIANS_DYNAMIC_SYSTEM_HPP
#define HAMILTONIANS_DYNAMIC_SYSTEM_HPP

#include "Hamiltonian.hpp"

namespace Integrators
{
    namespace Dynamics
    {

        /// \brief system implements the dynamic system generated by the Hamiltonian Ham
        /// \tparam Ham the Hamiltonian generating the dynamic system
        /// \param direction the direction along which the normalization takes place
        /// \param s the phase space position
        /// \param dsdt the time derivatives are returned in this parameter
        template<typename Ham>
        inline Geometry::State2 dynamic_system (const Ham& ham, const Geometry::State2& s)
        {
          const auto dHds = ham.derivative(s);
          return Geometry::State2(dHds.p(), -dHds.q());
        }


        /// \brief system_along_direction  normalizes the dynamic system generated by the Hamiltonian Ham so that the  normalized
        /// derivative of the quantity s' = s*direction is equal to 1.
        /// \tparam Ham the Hamiltonian type
        /// \param ham the hamiltonian
        /// \param direction the direction along which the normalization takes place
        /// \param s the phase space position
        /// \return the normalized time derivatives of the extended phase space
        ///
        /// It is assumed that direction is not perpendicular to the dynamics flow. No check is carried out for this.
        /// Otherwise, division by zero takes place. Typically direction is chosen as the vector perpendicular to a
        /// Poincare cut in phase space.
        template<typename Ham>
        inline Geometry::State3 dynamic_system_along_direction (const Ham& ham, const Geometry::State2& direction, const Geometry::State2& s)
        {
          const auto dsdt = dynamic_system(ham, s);
          const auto deriv_along_direction = dsdt*direction;
          return Geometry::State3{dsdt,1}/deriv_along_direction;
        }

        /// \brief SystemAlongDirection
        /// \tparam Ham the Hamiltonian generating the dynamic system
        ///
        /// It is assumed that direction is not perpendicular to the dynamics flow. No check is carried out for this.
        /// Otherwise, division by zero takes place. Typically direction is chosen as the vector perpendicular to a
        /// Poincare cut in phase space.
        template<typename Ham>
        class
        SystemAlongDirection {
         private :
          Ham ham_;
          Geometry::State2 direction_;

         public:
          ///
          /// \param direction
          explicit SystemAlongDirection (Ham ham,const Geometry::State2& direction) noexcept
              : ham_{std::move(ham)}, direction_{direction}
          { };

          /// calculate the normalized time derivatives
          /// \param s the phase space position
          /// \param dsdt the time derivatives are returned in this parameter
          /// \param t time
          inline Geometry::State3 operator() (const Geometry::State3& s_extended, Geometry::State3& dsdt_extended, double t) const
          {

            Geometry::State2 s{s_extended};
            dsdt_extended = dynamic_system_along_direction(ham_, direction_, s);

          };
        };
    }
}

#endif //HAMILTONIANS_DYNAMIC_SYSTEM_HPP
