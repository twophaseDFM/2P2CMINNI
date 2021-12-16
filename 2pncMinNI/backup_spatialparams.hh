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
    , finePcKrSwCurve_("SpatialParams.FineMaterial")
    , coarsePcKrSwCurve_("SpatialParams.CoarseMaterial")
//    , deltaPorosity_(gridGeometry->gridView().size(0),0.0)
//    , poro_(gridGeometry->gridView().size(0),0.0)
    {
        layerBottom_ 		= 0;
//        layerBottom_ = 22.5;

        // intrinsic permeabilities
        fineK_ 				= getParam<Scalar>("SpatialParams.fineK");
        coarseK_ 			= getParam<Scalar>("SpatialParams.coarseK");

        // porosities
        finePorosity_ 		= getParam<Scalar>("SpatialParams.finePorosity");;
        coarsePorosity_ 	= getParam<Scalar>("SpatialParams.coarsePorosity");

		RockXMin_ 			= getParam<Scalar>("Rock.RockXmin");
		RockXMax_ 			= getParam<Scalar>("Rock.RockXmax");

//		molarDensityH2O_ 	= getParam<Scalar>("Rock.molarDensityH2O");
//		molarDensityCaIon_  = getParam<Scalar>("Rock.molarDensityCaIon");

//		initialPorosity_ 	= getParam<Scalar>("SpatialParams.initialPorosity");

//      deltaPorosity_.resize(gridGeometry->gridView().size(0));

//		for (const auto& element : elements(this->gridGeometry().gridView()))
//		{
//        	auto fvGeometry = localView(*gridGeometry);
//        	fvGeometry.bind(element);
//
//        	for (const auto& scv : scvs(fvGeometry))
//			{
//        		const auto idx = scv.dofIndex();
//        		poro_[idx] = initialPorosity_;
////        		poro_[idx] = (deltaPorosity_[idx]) * initialPorosity_;
//			}
//		}
    }

    /*!
     * \brief Returns the intrinsic permeability tensor \f$[m^2]\f$
     *
     * \param globalPos The global position
     */
//    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
//    {
//        if (isRock_(globalPos))
//            return fineK_;
//        return coarseK_;
//    }

    /*!
     * \brief Returns the porosity \f$[-]\f$
     *
     * \param globalPos The global position
     */
//    Scalar porosityAtPos(const GlobalPosition& globalPos) const
//    {
//        if (isRock_(globalPos))
//            return finePorosity_;
//        return coarsePorosity_;
//    }

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
    	const auto& globalPos = scv.center();
        auto priVars = evalSolution(element, element.geometry(), elemSol, scv.center(), /*ignoreState=*/true);
		if (isRock_(globalPos))
			return finePorosity_ + priVars[3];
		else
			return coarsePorosity_ + priVars[3];
    }

    /*!
     * permeability variation due to the change of porosity
     */
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {
    	const auto& globalPos = scv.center();
    	using std::min;
//		auto priVars = evalSolution(element, element.geometry(), elemSol, scv.center(), /*ignoreState=*/true);
    	auto poro = porosity(element, scv, elemSol);
        if (isRock_(globalPos))
			return permLaw_.evaluatePermeability(fineK_, finePorosity_, poro);
        else
			return permLaw_.evaluatePermeability(coarseK_, coarsePorosity_, poro);
    }


    /*!
     * \brief Returns the parameter object for the capillary-pressure/
     *        saturation material law
     *
     * \param globalPos The global position
     */
    const auto fluidMatrixInteractionAtPos(const GlobalPosition& globalPos) const
    {
        if (isRock_(globalPos))
            return makeFluidMatrixInteraction(finePcKrSwCurve_);
        return makeFluidMatrixInteraction(coarsePcKrSwCurve_);
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

//    Scalar deltaPorosity(const Element& element,
//    				  const SubControlVolume& scv) const
//    { return deltaPorosity_[scv.dofIndex()]; }
//
//    void setdeltaPorosity(const std::vector<Scalar> s)
//    { deltaPorosity_ = s; }
//
//    const std::vector<Scalar>& getPermField() const
//    {return poro_; }


private:
    bool isRock_(const GlobalPosition &globalPos) const
    { return globalPos[0] > RockXMax_ || globalPos[0] < RockXMin_; }

    Scalar fineK_;
    Scalar coarseK_;
    Scalar layerBottom_;
    Scalar RockXMax_ ;
	Scalar RockXMin_ ;

    Scalar finePorosity_;
    Scalar coarsePorosity_;
    Scalar referencePorosity_ ;

    const PcKrSwCurve finePcKrSwCurve_;
    const PcKrSwCurve coarsePcKrSwCurve_;

    PermeabilityKozenyCarman<PermeabilityType> permLaw_;

//    std::vector<Scalar> deltaPorosity_ ;
//    std::vector<Scalar> poro_ ;
//    Scalar initialPorosity_ ;
//    Scalar deltaPorosity_;

//    Scalar molarDensityH2O_ , molarDensityCaIon_ ;


};

} // end namespace Dumux

#endif
