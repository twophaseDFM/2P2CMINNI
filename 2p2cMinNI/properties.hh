// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup TwoPTwoCTests
 * \brief The properties of the problem where air is injected under a low permeable layer in a depth
 *        of 2700m.
 */

#ifndef DUMUX_INJECTION_PROPERTIES_HH
#define DUMUX_INJECTION_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>
//#include <dune/alugrid/grid.hh>
#include <dune/grid/uggrid.hh>
#include <dumux/discretization/elementsolution.hh>

#include <dumux/discretization/cellcentered/mpfa/omethod/staticinteractionvolume.hh>
#include <dumux/discretization/ccmpfa.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>

#include <dumux/porousmediumflow/2p2c/model.hh>
#include <dumux/porousmediumflow/2pncmin/model.hh>

//#include <dumux/material/fluidsystems/brineco2.hh>
#include <dumux/material/components/calcite.hh>
#include <dumux/material/components/granite.hh>
#include <dumux/material/solidsystems/compositionalsolidphase_dejian.hh>
#include <dumux/material/fluidsystems/brineco2_min.hh>
#include "co2tables.hh"

#ifndef ENABLECACHING
#define ENABLECACHING 0
#endif

#include "spatialparams.hh"
#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct Injection { using InheritsFrom = std::tuple<TwoPNCMinNI>; };
struct InjectionBox { using InheritsFrom = std::tuple<Injection, BoxModel>; };
struct InjectionCCTpfa { using InheritsFrom = std::tuple<Injection, CCTpfaModel>; };
struct InjectionCCMpfa { using InheritsFrom = std::tuple<Injection, CCMpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::Injection>
{
	using Scalar = GetPropType<TypeTag, Properties::Scalar>;
	using type = Dune::YaspGrid<2, Dune::TensorProductCoordinates<Scalar, 2>>;
};
//struct Grid<TypeTag, TTag::Injection> { using type = Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>; };
//struct Grid<TypeTag, TTag::Injection> { using type = Dune::UGGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::Injection> { using type = InjectionProblem<TypeTag>; };

// Set fluid configuration
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Injection>
{
    using type = FluidSystems::BrineCO2Min<GetPropType<TypeTag, Properties::Scalar>,
											 HeterogeneousCO2Tables::CO2Tables,
											 Components::TabulatedComponent<Components::H2O<GetPropType<TypeTag, Properties::Scalar>>>,
											 FluidSystems::BrineCO2MinDefaultPolicy</*constantSalinity=*/1, /*simpleButFast=*/1>>;
};

// Set solid configuration
template<class TypeTag>
struct SolidSystem<TypeTag, TTag::Injection>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ComponentOne = Components::Calcite<Scalar>;
    using ComponentThree = Components::Granite<Scalar>;
    using ComponentTwo = Components::NaCl<Scalar>;
    static constexpr int numInertComponents = 1;
    using type = SolidSystems::CompositionalSolidPhase_dejian<Scalar, ComponentOne, ComponentTwo, ComponentThree, numInertComponents>;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Injection>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = InjectionSpatialParams<GridGeometry, Scalar>;
};

// Define whether mole(true) or mass (false) fractions are used
template<class TypeTag>
struct UseMoles<TypeTag, TTag::Injection> { static constexpr bool value = true; };

// Enable caching or not (reference solutions created without caching)
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::Injection> { static constexpr bool value = ENABLECACHING; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::Injection> { static constexpr bool value = ENABLECACHING; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::Injection> { static constexpr bool value = ENABLECACHING; };

// use the static interaction volume around interior vertices in the mpfa test
template<class TypeTag>
struct PrimaryInteractionVolume<TypeTag, TTag::InjectionCCMpfa>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using NodalIndexSet = GetPropType<TypeTag, Properties::DualGridNodalIndexSet>;

    // structured two-d grid
    static constexpr int numIvScvs = 4;
    static constexpr int numIvScvfs = 4;

    // use the default traits
    using Traits = CCMpfaODefaultStaticInteractionVolumeTraits< NodalIndexSet, Scalar, numIvScvs, numIvScvfs >;
public:
    using type = CCMpfaOStaticInteractionVolume< Traits >;
};

} // end namespace Dumux::Properties

#endif
