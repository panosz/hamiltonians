//
// Created by Panagiotis Zestanakis on 08/10/18.
//

#include "Integration.hpp"
namespace Integrators
{

    template std::vector<Geometry::State2_Extended>
    calculate_crossings<Hamiltonian::FreeParticle> (const Hamiltonian::FreeParticle& hamiltonian,
                                                    const Geometry::State2& s_start,
                                                    const Geometry::PeriodicQSurfaceCrossObserver& periodicQSurfaceCrossObserver,
                                                    const TimeInterval& integrationTime,
                                                    const IntegrationOptions& options);

    template std::vector<Geometry::State2_Extended>
    calculate_crossings<Hamiltonian::HarmonicOscillator> (const Hamiltonian::HarmonicOscillator& hamiltonian,
                                                          const Geometry::State2& s_start,
                                                          const Geometry::PeriodicQSurfaceCrossObserver& periodicQSurfaceCrossObserver,
                                                          const TimeInterval& integrationTime,
                                                          const IntegrationOptions& options);

    template std::vector<Geometry::State2_Extended>
    calculate_crossings<Hamiltonian::PendulumHamiltonian> (const Hamiltonian::PendulumHamiltonian& hamiltonian,
                                                           const Geometry::State2& s_start,
                                                           const Geometry::PeriodicQSurfaceCrossObserver& periodicQSurfaceCrossObserver,
                                                           const TimeInterval& integrationTime,
                                                           const IntegrationOptions& options);
    template std::vector<Geometry::State2_Extended>
    calculate_crossings<Hamiltonian::DuffingHamiltonian> (const Hamiltonian::DuffingHamiltonian& hamiltonian,
                                                          const Geometry::State2& s_start,
                                                          const Geometry::PeriodicQSurfaceCrossObserver& periodicQSurfaceCrossObserver,
                                                          const TimeInterval& integrationTime,
                                                          const IntegrationOptions& options);

    template std::vector<Geometry::State2_Extended>
    calculate_crossings<Hamiltonian::FreeParticle> (const Hamiltonian::FreeParticle& hamiltonian,
                                                    const Geometry::State2& s_start,
                                                    const Geometry::Line& cross_line,
                                                    const TimeInterval& integrationTime,
                                                    const IntegrationOptions& options);

    template std::vector<Geometry::State2_Extended>
    calculate_crossings<Hamiltonian::HarmonicOscillator> (const Hamiltonian::HarmonicOscillator& hamiltonian,
                                                          const Geometry::State2& s_start,
                                                          const Geometry::Line& cross_line,
                                                          const TimeInterval& integrationTime,
                                                          const IntegrationOptions& options);

    template std::vector<Geometry::State2_Extended>
    calculate_crossings<Hamiltonian::PendulumHamiltonian> (const Hamiltonian::PendulumHamiltonian& hamiltonian,
                                                           const Geometry::State2& s_start,
                                                           const Geometry::Line& cross_line,
                                                           const TimeInterval& integrationTime,
                                                           const IntegrationOptions& options);

    template std::vector<Geometry::State2_Extended>
    calculate_crossings<Hamiltonian::DuffingHamiltonian> (const Hamiltonian::DuffingHamiltonian& hamiltonian,
                                                          const Geometry::State2& s_start,
                                                          const Geometry::Line& cross_line,
                                                          const TimeInterval& integrationTime,
                                                          const IntegrationOptions& options);

    template
    Geometry::State2_Extended
    calculate_first_crossing<Hamiltonian::FreeParticle> (const Hamiltonian::FreeParticle& hamiltonian,
                                                         const Geometry::State2& s_start,
                                                         const Geometry::Line& cross_line,
                                                         const TimeInterval& integrationTime,
                                                         const IntegrationOptions& options);
    template
    Geometry::State2_Extended
    calculate_first_crossing<Hamiltonian::HarmonicOscillator> (const Hamiltonian::HarmonicOscillator& hamiltonian,
                                                               const Geometry::State2& s_start,
                                                               const Geometry::Line& cross_line,
                                                               const TimeInterval& integrationTime,
                                                               const IntegrationOptions& options);
    template
    Geometry::State2_Extended
    calculate_first_crossing<Hamiltonian::PendulumHamiltonian> (const Hamiltonian::PendulumHamiltonian& hamiltonian,
                                                                const Geometry::State2& s_start,
                                                                const Geometry::Line& cross_line,
                                                                const TimeInterval& integrationTime,
                                                                const IntegrationOptions& options);
    template
    Geometry::State2_Extended
    calculate_first_crossing<Hamiltonian::DuffingHamiltonian> (const Hamiltonian::DuffingHamiltonian& hamiltonian,
                                                               const Geometry::State2& s_start,
                                                               const Geometry::Line& cross_line,
                                                               const TimeInterval& integrationTime,
                                                               const IntegrationOptions& options);

    template
    Geometry::State2_Extended
    come_back_home<Hamiltonian::FreeParticle> (const Hamiltonian::FreeParticle& hamiltonian,
                                               const Geometry::State2& s_start,
                                               const TimeInterval& integrationTime,
                                               const IntegrationOptions& options);

    template
    Geometry::State2_Extended
    come_back_home<Hamiltonian::HarmonicOscillator> (const Hamiltonian::HarmonicOscillator& hamiltonian,
                                                     const Geometry::State2& s_start,
                                                     const TimeInterval& integrationTime,
                                                     const IntegrationOptions& options);
    template
    Geometry::State2_Extended
    come_back_home<Hamiltonian::PendulumHamiltonian> (const Hamiltonian::PendulumHamiltonian& hamiltonian,
                                                      const Geometry::State2& s_start,
                                                      const TimeInterval& integrationTime,
                                                      const IntegrationOptions& options);
    template
    Geometry::State2_Extended
    come_back_home<Hamiltonian::DuffingHamiltonian> (const Hamiltonian::DuffingHamiltonian& hamiltonian,
                                                     const Geometry::State2& s_start,
                                                     const TimeInterval& integrationTime,
                                                     const IntegrationOptions& options);



}