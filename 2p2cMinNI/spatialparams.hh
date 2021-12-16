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
 * \brief Definition of the spatial parameters for the injection problem
 *        which uses the isothermal two-phase two-component fully implicit model.
 */

#ifndef DUMUX_INJECTION_SPATIAL_PARAMS_HH
#define DUMUX_INJECTION_SPATIAL_PARAMS_HH

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/brookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/vangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/permeabilitykozenycarman.hh>
#include <dumux/material/fluidmatrixinteractions/porosityprecipitation.hh>
#include <dumux/material/spatialparams/gstatrandomfield.hh>

namespace Dumux {

/*!
 * \ingroup TwoPTwoCTests
 * \brief Definition of the spatial parameters for the injection problem
 *        which uses the isothermal two-phase two-component fully implicit model.
 */
template<class GridGeometry, class Scalar>
class InjectionSpatialParams
: public FVSpatialParams<GridGeometry, Scalar,
                         InjectionSpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVSpatialParams<GridGeometry, Scalar, InjectionSpatialParams<GridGeometry, Scalar>>;

    static constexpr int dimWorld = GridView::dimensionworld;

//    using PcKrSwCurve = FluidMatrix::BrooksCoreyDefault<Scalar>;
    using PcKrSwCurve = FluidMatrix::VanGenuchtenDefault<Scalar>;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    //! Export the type used for the permeability
    using PermeabilityType = Scalar;

    InjectionSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , pcKrSwCurve_("SpatialParams.material")
    {
        // intrinsic permeabilities
        // porosities
        porosity_ 						= getParam<Scalar>("SpatialParams.porosity");
        initialMinVolumeFraction_ 		= getParam<Scalar>("Initial.initialMinVolumeFraction");

		initRandomField(*gridGeometry);
    }

    template<class SolidSystem, class ElementSolution>
    Scalar inertVolumeFraction(const Element& element,
                               const SubControlVolume& scv,
                               const ElementSolution& elemSol,
                               int compIdx) const
    { return 1.0 - porosity(element, scv, elemSol); }

    /*!
     * *\ porosity variation due to the chemical reaction
     * the porosity variaiton is the moleFraction of CaIonIdx (from 2pncmin)
     */
    template<class ElementSolution>
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const
    {
    	using std::abs;
        auto priVars = evalSolution(element, element.geometry(), elemSol, scv.center(), /*ignoreState=*/true);
		return porosity_ + abs(initialMinVolumeFraction_ - priVars[4] - priVars[5]);
    }

    /*!
     * permeability variation due to the change of porosity
     */
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {
//    	const auto& globalPos = scv.center();
//    	using std::min;
//    	auto priVars = evalSolution(element, element.geometry(), elemSol, scv.center(), /*ignoreState=*/true);
		using Dune::power;
		auto poro = porosity(element, scv, elemSol);
		auto factor = power((1 - porosity_)/(1 - poro), 2) * power(poro/porosity_ , 3);
    	return randomPermeability_[scv.dofIndex()] * factor;
    }

    /*!
     * \brief Returns the parameter object for the capillary-pressure/
     *        saturation material law
     *
     * \param globalPos The global position
     */
    const auto fluidMatrixInteractionAtPos(const GlobalPosition& globalPos) const
    {
        return makeFluidMatrixInteraction(pcKrSwCurve_);
    }

    /*!
     * \brief Function for defining which phase is to be considered as the wetting phase.
     *
     * \param globalPos The position of the center of the element
     * \return The wetting phase index
     */
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    { return FluidSystem::BrineIdx; }
//    { return FluidSystem::H2OIdx; }

	void initRandomField(const GridGeometry& gg)
	{
		const auto& gridView = gg.gridView();
		const auto& elementMapper = gg.elementMapper();
		const auto gStatControlFile = getParam<std::string>("Gstat.ControlFile");
		const auto gStatInputFile = getParam<std::string>("Gstat.InputFile");
		const auto outputFilePrefix = getParam<std::string>("Gstat.OutputFilePrefix");

		// create random permeability object
		using RandomField = GstatRandomField<GridView, Scalar>;
		RandomField randomApertureField(gridView, elementMapper);
		randomApertureField.create(gStatControlFile,
								   gStatInputFile,
								   outputFilePrefix + ".dat",
								   RandomField::FieldType::log10,
								   true);
		randomPermeability_.resize(gridView.size(0), 0.0);

		// copy vector from the temporary gstat object
		randomPermeability_ = randomApertureField.data();
	}

    const std::vector<Scalar>& getPermField() const
    { return randomPermeability_; }


private:
    Scalar porosity_;
    Scalar referencePorosity_ ;
    Scalar initialMinVolumeFraction_;

    const PcKrSwCurve pcKrSwCurve_;

    PermeabilityKozenyCarman<PermeabilityType> permLaw_;

    std::vector<Scalar> randomPermeability_;
};

} // end namespace Dumux

#endif
