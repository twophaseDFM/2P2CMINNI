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
 * \brief Problem where air is injected under a low permeable layer in a depth of 2700m.
 */

#ifndef DUMUX_INJECTION_PROBLEM_HH
#define DUMUX_INJECTION_PROBLEM_HH

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

/*!
 * \ingroup TwoPTwoCTests
 * \brief Problem where air is injected under a low permeable layer in a depth of 2700m.
 *
 * The domain is sized 60m times 40m and consists of two layers, a moderately
 * permeable one (\f$ K=10e-12\f$) for \f$ y<22m\f$ and one with a lower
 * permeablility (\f$ K=10e-13\f$) in the rest of the domain.
 *
 * A mixture of Nitrogen and Water vapor, which is composed according to the
 * prevailing conditions (temperature, pressure) enters a water-filled aquifer.
 * This is realized with a solution-dependent Neumann boundary condition at the
 * right boundary (\f$ 5m<y<15m\f$). The aquifer is situated 2700m below sea level.
 * The injected fluid phase migrates upwards due to buoyancy.
 * It accumulates and partially enters the lower permeable aquitard.
 *
 * The model is able to use either mole or mass fractions. The property useMoles
 * can be set to either true or false in the problem file. Make sure that the
 * according units are used in the problem set-up.
 * The default setting for useMoles is true.
 *
 * This problem uses the \ref TwoPTwoCModel.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_box2p2c</tt> or
 * <tt>./test_cc2p2c</tt>
 */
template <class TypeTag>
class InjectionProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using SolidSystem = GetPropType<TypeTag, Properties::SolidSystem>;

    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;

    // primary variable indices
    enum
    {
        pressureIdx = Indices::pressureIdx,
        switchIdx = Indices::switchIdx
    };

    // phase indices
    enum
    {
        nPhaseIdx = FluidSystem::CO2Idx,
		wPhaseIdx = FluidSystem::liquidPhaseIdx,
        H2OIdx = FluidSystem::BrineIdx,
        CO2Idx = FluidSystem::CO2Idx,
		CaIonIdx = FluidSystem::CaIonIdx,
    };

    // Solid phase
    enum
	{
    	sPhaseIdx = SolidSystem::comp0Idx,
		CaCO3Idx = Indices::conti0EqIdx + FluidSystem::numComponents,
	};

    // phase presence
    enum { wPhaseOnly = Indices::firstPhaseOnly };

    // equation indices
    enum
    {
        contiH2OEqIdx = Indices::conti0EqIdx + H2OIdx,
        contiCO2EqIdx = Indices::conti0EqIdx + CO2Idx,

        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
    };

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
	using FluidState = GetPropType<TypeTag, Properties::FluidState>;
	using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
	using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;

    //! Property that defines whether mole or mass fractions are used
    static constexpr bool useMoles = ModelTraits::useMoles();

public:
    InjectionProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
    	name_              		 = getParam<std::string>("Problem.Name");
        InjectionTemperature_    = getParam<Scalar>("Problem.InjectionTemperature");
        InjectionPressure_ 		 = getParam<Scalar>("Problem.InjectionPressure");
        InjectionRate_           = getParam<Scalar>("Problem.InjectionRate");
        ProductionRate_          = getParam<Scalar>("Problem.ProductionRate");
        depthBOR_           	 = getParam<Scalar>("Problem.DepthBOR");
        DomainHeight_ 			 = getParam<Scalar>("Problem.DomainHeight");
        InjectionHeight_ 		 = getParam<Scalar>("Problem.InjectionHeight");
        InjectionX_ 			 = getParam<Scalar>("Well.InjectionX");
        InjectionY_ 			 = getParam<Scalar>("Well.InjectionY");
        ProductionX_ 			 = getParam<Scalar>("Well.ProductionX");
        ProductionY_ 			 = getParam<Scalar>("Well.ProductionY");
        Wellradical_			 = getParam<Scalar>("Well.Wellradical");
        InjectionRateSource_     = getParam<Scalar>("Problem.InjectionRateSource");
        ProductionRateSource_    = getParam<Scalar>("Problem.InjectionRateSource");
        eleNi_  				 = getParam<Scalar>("Problem.eleNi");
        Swr_  				     = getParam<Scalar>("Problem.Swr");
        useSource_ 				 = getParam<Scalar>("Problem.useSource");
        useDirichlet_			 = getParam<Scalar>("Problem.useDirichlet");
        initialCaIon_ 			 = getParam<Scalar>("Rock.initialCaIon");

        // stating in the console whether mole or mass fractions are used
        if(useMoles)
            std::cout<<"problem uses mole-fractions"<<std::endl;
        else
            std::cout<<"problem uses mass-fractions"<<std::endl;

        unsigned int codim = GetPropType<TypeTag, Properties::GridGeometry>::discMethod == DiscretizationMethod::box ? dim : 0;

        permeability_.resize(gridGeometry->gridView().size(codim));
        deltaPorosity_.resize(gridGeometry->gridView().size(codim));
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the problem name
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string& name() const
    { return name_; }

    /*!
     * \brief The current time.
     */
    void setTime( Scalar time )
    {
        time_ = time;
    }

    /*!
     * \brief The time step size.
     *
     * This is used to calculate the source term.
     */
    void setTimeStepSize( Scalar timeStepSize )
     {
        timeStepSize_ = timeStepSize;
     }

    void setX (SolutionVector result)
    {
    	result_ = result;
    }

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment
     *
     * \param globalPos The global position
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;
        if (useSource_)
        {
        	if (globalPos[0] < eps_ || globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_)
				bcTypes.setAllDirichlet();
			else
				bcTypes.setAllNeumann();
        }
        else
        {
//        	        if (globalPos[0] < eps_ && !isProduction(globalPos))
//        	            bcTypes.setAllDirichlet();
			if (globalPos[0] < eps_)
				bcTypes.setAllDirichlet();
			else
				bcTypes.setAllNeumann();
        }
        return bcTypes;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet boundary segment
     *
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);
        values = initial_(globalPos);

        if (useDirichlet_)
        {
			if (globalPos[0] < eps_ && globalPos[1] < 20 + eps_)
			{
				values.setState(Indices::bothPhases);
				values[pressureIdx] = InjectionPressure_;
				values[Indices::switchIdx] = 0;
				values[Indices::temperatureIdx] = InjectionTemperature_;
			}
        }

//		if (globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_)
//		{
//			values.setState(Indices::bothPhases);
//			values[pressureIdx] = 0.9*values[pressureIdx];
//		}

    	return values;
    }

//    template<class ElementVolumeVariables>
//    PrimaryVariables dirichlet(const Element &element,
//    						   const SubControlVolumeFace &scvf,
//							   const ElementVolumeVariables &elemVolVars) const
//    {
//    	PrimaryVariables values(0.0);
//    	const auto& globalPos = scvf.ipGlobal();
//
//    	values = initial_(globalPos);
//
//
//    	return values;
//    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann
     *        boundary segment in dependency on the current solution.
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param elemVolVars All volume variables for the element
     * \param elemFluxVarsCache Flux variables caches for all faces in stencil
     * \param scvf The sub-control volume face
     *
     * This method is used for cases, when the Neumann condition depends on the
     * solution and requires some quantities that are specific to the fully-implicit method.
     * The \a values store the mass flux of each phase normal to the boundary.
     * Negative values indicate an inflow.
     */
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);

        const auto& globalPos = scvf.ipGlobal();
        const auto& initialValues = initialAtPos(globalPos);
        Scalar InjectionPressure = initialValues[pressureIdx];
        const auto& volVars = elemVolVars[scvf.insideScvIdx()];

		FluidState fs;

		fs.setPressure(nPhaseIdx, InjectionPressure); // assume pressure equality here
		fs.setTemperature(nPhaseIdx, InjectionTemperature_);

		if (useSource_)
		{}
		else
		{
			if (isInjection(globalPos))
			{
				values[contiCO2EqIdx] = -InjectionRate_/FluidSystem::molarMass(CO2Idx); //mole/(m^2*s) -> kg/(s*m^2)
				values[Indices::temperatureIdx] = values[contiCO2EqIdx] * FluidSystem::componentEnthalpy(fs, nPhaseIdx, CO2Idx) * FluidSystem::molarMass(CO2Idx);
				values[CaIonIdx] = volVars.porosity()
								 * FluidSystem::molarDensity(fs,nPhaseIdx)
								 * volVars.moleFraction(nPhaseIdx, H2OIdx)
								 * volVars.moleFraction(wPhaseIdx,CaIonIdx)
								 / timeStepSize_ ;

			}

			if (isProduction(globalPos))
			{
				values[contiCO2EqIdx] =
									ProductionRate_ /*m³/(s*m²)*/ * volVars.saturation(nPhaseIdx)
									* FluidSystem::molarDensity(volVars.fluidState(), nPhaseIdx)
									* volVars.moleFraction(nPhaseIdx, CO2Idx);
				values[contiCO2EqIdx] +=
									ProductionRate_ * volVars.saturation(wPhaseIdx)
									* FluidSystem::molarDensity(volVars.fluidState(), wPhaseIdx)
									* volVars.moleFraction(wPhaseIdx, CO2Idx);

				values[contiH2OEqIdx] =
									ProductionRate_ /*m³/(s*m²)*/ * volVars.saturation(nPhaseIdx)
									* FluidSystem::molarDensity(volVars.fluidState(), nPhaseIdx)
									* volVars.moleFraction(nPhaseIdx, H2OIdx);
				values[contiH2OEqIdx] +=
									ProductionRate_ * volVars.saturation(wPhaseIdx)
									* FluidSystem::molarDensity(volVars.fluidState(), wPhaseIdx)
									* volVars.moleFraction(wPhaseIdx, H2OIdx);

				values[Indices::energyEqIdx]  = ProductionRate_ * volVars.saturation(nPhaseIdx) * FluidSystem::density(volVars.fluidState(), nPhaseIdx) * FluidSystem::enthalpy(volVars.fluidState(), nPhaseIdx);
				values[Indices::energyEqIdx] += ProductionRate_ * volVars.saturation(wPhaseIdx) * FluidSystem::density(volVars.fluidState(), wPhaseIdx) * FluidSystem::enthalpy(volVars.fluidState(), wPhaseIdx);
			}
		}
        return values;
    }

    NumEqVector source(const Element &element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume &scv) const
    {
        NumEqVector source(0.0);
        const auto& volVars = elemVolVars[scv];
        const auto& globalPos = scv.center();
        const auto initialValues = initialAtPos(globalPos);
		Scalar InjectionPressure = initialValues[pressureIdx];

		/* *****************  REACTION  ***********************/
        // TODO determine the kinetics of the reaction
        Scalar k0 = getParam<Scalar>("Reaction.k0");
        Scalar E = getParam<Scalar>("Reaction.AE");
        Scalar R = getParam<Scalar>("Reaction.R");
        Scalar CaMax = getParam<Scalar>("Reaction.MaxCaConcentration");
//
        Scalar K0 = k0 * exp(-E/R/volVars.temperature()) * volVars.saturation(wPhaseIdx) * volVars.moleFraction(wPhaseIdx,CO2Idx) * volVars.pressure(CO2Idx);
////        Scalar K0 = k0;
        Scalar r1 = K0 * (1 - volVars.moleFraction(wPhaseIdx, CaIonIdx)/CaMax);
//
//		Scalar r2 = volVars.moleFraction(wPhaseIdx, CaIonIdx) > CaMax ?
//				    volVars.porosity() * volVars.molarDensity(wPhaseIdx)
//									   * volVars.saturation(wPhaseIdx)
//									   * abs(volVars.moleFraction(wPhaseIdx, CaIonIdx) - CaMax)
//									   / timeStepSize_
//									   : 0;
//
//		Scalar r = r1 - r2;

        source[CO2Idx] += -r1;
//        if (isReactionZone (globalPos))
		source[CaIonIdx] += r1;
		source[CaCO3Idx] += r1;

		/* ********************** INJECTION AND PRODUCTION  *******************/
		Scalar Volume = 4/3 * Wellradical_ * Wellradical_ * Wellradical_ * 3.14;

		if (useSource_)
		{
			if (useDirichlet_)
			{}
			else
			{
				FluidState fs;
				fs.setPressure(nPhaseIdx, InjectionPressure); // assume pressure equality here
				fs.setTemperature(nPhaseIdx, InjectionTemperature_);

				if (isInjectionWell (globalPos))
				{
					source[contiCO2EqIdx] += InjectionRateSource_ * FluidSystem::CO2::gasMolarDensity(InjectionTemperature_, InjectionPressure) / Volume; // kg/(s*m^2) flow rate * density
					source[Indices::energyEqIdx] = source[contiCO2EqIdx] * FluidSystem::componentEnthalpy(fs, nPhaseIdx, CO2Idx) * FluidSystem::molarMass(CO2Idx)/*J/kg*/; // W/(m^2)

					source[CaIonIdx] += - source[contiCO2EqIdx] * volVars.moleFraction(nPhaseIdx, H2OIdx) * volVars.moleFraction(wPhaseIdx,CaIonIdx);

//					source[CaIonIdx] += -volVars.porosity() * volVars.saturation(nPhaseIdx)
//									 * FluidSystem::molarDensity(fs,nPhaseIdx)
//									 * volVars.moleFraction(nPhaseIdx, H2OIdx)
//									 * volVars.moleFraction(wPhaseIdx,CaIonIdx)
//									 / timeStepSize_ ;

//					source[CaIonIdx] += - InjectionRateSource_ * FluidSystem::molarDensity(volVars.fluidState(),wPhaseIdx)
//									 * (volVars.moleFraction(wPhaseIdx, CaIonIdx) - 1e-5) / std::max(volVars.saturation(wPhaseIdx),Swr_)
//									 / volVars.porosity() / volVars.moleFraction(nPhaseIdx, H2OIdx);

	//				source[CaIonIdx] += - (volVars.moleFraction(wPhaseIdx, CaIonIdx) - 1e-5)
	//								    * FluidSystem::molarDensity(volVars.fluidState(), wPhaseIdx)
	//							    	/ timeStepSize_;

	//								   (volVars.saturation(wPhaseIdx) - 1)
	//								 * FluidSystem::molarDensity(volVars.fluidState(),wPhaseIdx)
	//								 * (volVars.moleFraction(wPhaseIdx, CaIonIdx) - 1e-5)
	//								 / timeStepSize_;

	//								 - InjectionRateSource_ /*m³/(s*m²)*/ * volVars.saturation(nPhaseIdx)
	//								 * FluidSystem::molarDensity(volVars.fluidState(), wPhaseIdx)
	//								 * volVars.moleFraction(wPhaseIdx, CaIonIdx);
				}

				if (isProductionWell(globalPos))
				{
					source[contiCO2EqIdx] =
										- ProductionRateSource_ /*m³/(s*m²)*/ * volVars.saturation(nPhaseIdx)
										* FluidSystem::molarDensity(volVars.fluidState(), nPhaseIdx)
										* volVars.moleFraction(nPhaseIdx, CO2Idx) / Volume;
					source[contiCO2EqIdx] +=
										- ProductionRateSource_ * volVars.saturation(wPhaseIdx)
										* FluidSystem::molarDensity(volVars.fluidState(), wPhaseIdx)
										* volVars.moleFraction(wPhaseIdx, CO2Idx) / Volume;

					source[contiH2OEqIdx] =
										- ProductionRateSource_ /*m³/(s*m²)*/ * volVars.saturation(nPhaseIdx)
										* FluidSystem::molarDensity(volVars.fluidState(), nPhaseIdx)
										* volVars.moleFraction(nPhaseIdx, H2OIdx) / Volume;
					source[contiH2OEqIdx] +=
										- ProductionRateSource_ * volVars.saturation(wPhaseIdx)
										* FluidSystem::molarDensity(volVars.fluidState(), wPhaseIdx)
										* volVars.moleFraction(wPhaseIdx, H2OIdx)/ Volume;

//					source[contiCO2EqIdx] =
//										- ProductionRateSource_ * volVars.saturation(nPhaseIdx) * FluidSystem::molarDensity(volVars.fluidState(), nPhaseIdx) / Volume;
//
//					source[contiH2OEqIdx] =
//										- ProductionRateSource_ * volVars.saturation(wPhaseIdx) * FluidSystem::molarDensity(volVars.fluidState(), wPhaseIdx) / Volume;

					source[Indices::energyEqIdx]  = -ProductionRateSource_ * volVars.saturation(nPhaseIdx) * FluidSystem::density(volVars.fluidState(), nPhaseIdx) * FluidSystem::enthalpy(volVars.fluidState(), nPhaseIdx) / Volume;
					source[Indices::energyEqIdx] += -ProductionRateSource_ * volVars.saturation(wPhaseIdx) * FluidSystem::density(volVars.fluidState(), wPhaseIdx) * FluidSystem::enthalpy(volVars.fluidState(), wPhaseIdx) / Volume;

					source[CaIonIdx] =
									 - ProductionRateSource_ /*m³/(s*m²)*/ * volVars.saturation(wPhaseIdx)
									 * FluidSystem::molarDensity(volVars.fluidState(), wPhaseIdx)
									 * volVars.moleFraction(wPhaseIdx, CaIonIdx) / Volume;
				}
			}
		}

        return source;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluates the initial values for a control volume.
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    { return initial_(globalPos); }

    /*!
     *
     */
    /*
     *  / Calculate the output results
     */
    void calculateOutput(const GridVariables& gridVariables, const SolutionVector& sol, Scalar time)
    {
        NumEqVector Source(0.0);
        PrimaryVariables InitialValues;
        Scalar overPin = 0;
        Scalar overPon = 0;
        Scalar eleNi = 0;
        Scalar eleNp = 0;
        Scalar outputSn = 0;
        Scalar outputCO2 = 0;
        Scalar outputCaIon = 0;
        Scalar rho11 = 0;
        Scalar rho10 = 0;
        Scalar rho0 = 0.0;
        Scalar rho_try = 0;
        Scalar wphasei, wphaseo, molarDi, molarDo;
		Scalar Volume = 4/3 * Wellradical_ * Wellradical_ * Wellradical_ * 3.14;

        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
			auto fvGeometry = localView(this->gridGeometry());
			fvGeometry.bind(element);

			for (auto&& scv : scvs(fvGeometry))
			{
				const auto& globalPos = scv.center();
				auto elemVolVars = localView(gridVariables.curGridVolVars());
				elemVolVars.bind(element, fvGeometry, sol);
				const auto& volVars = elemVolVars[scv];
				Source = this->source(element, fvGeometry, elemVolVars, scv);
				InitialValues = initialAtPos(globalPos);

				if (isProductionWell(globalPos))
				{
					eleNp += 1;
					outputSn += volVars.saturation(nPhaseIdx);
					outputCO2 += volVars.saturation(nPhaseIdx) * volVars.moleFraction(nPhaseIdx, CO2Idx) + volVars.saturation(wPhaseIdx) * volVars.moleFraction(wPhaseIdx, CO2Idx);
					outputCaIon += volVars.saturation(wPhaseIdx) * volVars.moleFraction(wPhaseIdx, FluidSystem::CaIonIdx);
					overPon += volVars.pressure(nPhaseIdx) - InitialValues[pressureIdx];
					rho11 = Source[CO2Idx];
					rho10 = Source[H2OIdx];
					wphaseo = volVars.saturation(wPhaseIdx);
					molarDo = FluidSystem::molarDensity(volVars.fluidState(), wPhaseIdx);
				}

				if (isInjectionWell(globalPos))
				{
					eleNi += 1;
					overPin += volVars.pressure(nPhaseIdx) - InitialValues[pressureIdx];
					rho0 = InjectionRateSource_ * FluidSystem::CO2::gasMolarDensity(InjectionTemperature_, 2.3e7) / Volume;
					rho_try = - InjectionRateSource_ * volVars.saturation(wPhaseIdx) * FluidSystem::molarDensity(volVars.fluidState(), wPhaseIdx) / Volume;
					wphasei = volVars.saturation(wPhaseIdx);
					molarDi = FluidSystem::CO2::gasMolarDensity(InjectionTemperature_, 2.3e7);
				}
			}
        }

        outputSn = outputSn / eleNp;
        outputCO2 = outputCO2 / eleNp;
        outputCaIon = outputCaIon / eleNp;
        overPin = overPin / 1e6;
        overPon = overPon / 1e6;

        std::cout << "outputSn:" << outputSn << " outputCO2:" << outputCO2 << " outputCaIon:" << outputCaIon << '\n';
        std::cout << "overPin:" << overPin << " overPon:" << overPon << '\n';
        std::cout << "eleNp:" << eleNp << " eleNi:" << eleNi << '\n';
        std::cout << "contiCO2EqIdx:" << contiCO2EqIdx << " contiH2OEqIdx:" << contiH2OEqIdx << " energyEqIdx:" << Indices::energyEqIdx << " CaIonIdx:" << FluidSystem::CaIonIdx <<'\n';
        std::cout << "rho10 = " << rho10 << " rho11 = " << rho11 << " rho0 = " << rho0 << " rho_try = " << rho_try << '\n' ;
        std::cout << "wphasei = " << wphasei << " wphaseo = " << wphaseo << '\n';
        std::cout << "molarDi = " << molarDi << " molarDo = " << molarDo << '\n';
        std::cout << "Production = " << ProductionRateSource_ * wphaseo * molarDo / Volume << '\n' ;
        std::cout << "Injection = " << InjectionRateSource_ * molarDi / Volume << '\n' ;
//        fileSn_.open("1outputSn_3D.log",std::ios::app);
//        fileSn_ << time << " " << OutputSn << " " << InjectionSn << " " << OutputE << " " << OutputPressure << " " << InjectionPressure << " " << totalMassCO2 << std::endl;
//        fileSn_.close();
    }

    const std::vector<Scalar>& getPermeability()
    {
        return permeability_;
    }

    const std::vector<Scalar>& getPorosity()
    {
        return deltaPorosity_;
    }

    void updateVtkOutput(const SolutionVector& curSol)
    {
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            const auto elemSol = elementSolution(element, curSol, this->gridGeometry());

            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bindElement(element);

            for (auto&& scv : scvs(fvGeometry))
            {
                VolumeVariables volVars;
                volVars.update(elemSol, *this, element, scv);
                const auto dofIdxGlobal = scv.dofIndex();
                permeability_[dofIdxGlobal] = volVars.permeability();
                deltaPorosity_[dofIdxGlobal] = volVars.porosity();
            }
        }
    }

    // \}

private:
    /*!
     * \brief Evaluates the initial values for a control volume.
     *
     * The internal method for the initial condition
     *
     * \param globalPos The global position
     */
    PrimaryVariables initial_(const GlobalPosition &globalPos) const
    {
        PrimaryVariables priVars(0.0);
        priVars.setState(wPhaseOnly);
        Scalar initialTemperature = 413.15;

        Scalar densityW = FluidSystem::H2O::liquidDensity(initialTemperature, 1e5);

        Scalar pl = 1e5 - densityW*this->spatialParams().gravity(globalPos)[1]*(depthBOR_ - globalPos[1]);
//        Scalar moleFracLiquidN2 = pl*0.95/BinaryCoeff::H2O_N2::henry(temperature_);
        Scalar moleFracLiquidN2 = 0;
        Scalar moleFracLiquidH2O = 1.0 - moleFracLiquidN2;

        Scalar meanM =
            FluidSystem::molarMass(H2OIdx)*moleFracLiquidH2O +
            FluidSystem::molarMass(CO2Idx)*moleFracLiquidN2;
        if(useMoles)
        {
            //mole-fraction formulation
            priVars[switchIdx] = moleFracLiquidN2;
        }
        else
        {
            //mass fraction formulation
            Scalar massFracLiquidN2 = moleFracLiquidN2*FluidSystem::molarMass(CO2Idx)/meanM;
            priVars[switchIdx] = massFracLiquidN2;
        }
        priVars[pressureIdx] = pl;
        priVars[Indices::temperatureIdx] = 413.15;
        priVars[CaIonIdx] = initialCaIon_ ;
        priVars[CaCO3Idx] = 0.0;
        return priVars;
    }

    bool isInjection (const GlobalPosition globalPos) const
    {
    	return globalPos[1] < (DomainHeight_ + InjectionHeight_)/2
    			&& globalPos[1] > (DomainHeight_ - InjectionHeight_)/2
				&& globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_;
    }

    bool isProduction (const GlobalPosition globalPos) const
    {
    	return globalPos[1] < (DomainHeight_ + InjectionHeight_)/2
    			&& globalPos[1] > (DomainHeight_ - InjectionHeight_)/2
				&& globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_;
    }

    bool isInjectionWell (const GlobalPosition& globalPos) const
    {
        static const GlobalPosition WellCenter = {InjectionX_, InjectionY_};

        static const Scalar WellRadius = Wellradical_;

        const auto connecVec = (WellCenter-globalPos);

        return connecVec.two_norm() < WellRadius;
    }

    bool isProductionWell (const GlobalPosition& globalPos) const
    {
        static const GlobalPosition WellCenter = {ProductionX_, ProductionY_};

        static const Scalar WellRadius = Wellradical_;

        const auto connecVec = (WellCenter-globalPos);

        return connecVec.two_norm() < WellRadius;
    }

    bool isReactionZone (const GlobalPosition globalPos) const
    {
    	return globalPos[0] < 580;
    }



    Scalar depthBOR_, DomainHeight_, InjectionHeight_;
    static constexpr Scalar eps_ = 1e-6;
    Scalar InjectionRate_, InjectionTemperature_, InjectionPressure_ ;
    Scalar ProductionRate_ ;
    Scalar InjectionX_, InjectionY_, ProductionX_, ProductionY_, Wellradical_;
    Scalar InjectionRateSource_, ProductionRateSource_ ;
    bool useSource_, useDirichlet_;
    Scalar eleNi_;

    std::string name_;
    Scalar time_ = 0.0;
	Scalar timeStepSize_ = 0.0;
	SolutionVector result_ ;
	Scalar Swr_ ;
	Scalar initialCaIon_;

	std::vector<Scalar> permeability_, deltaPorosity_;
};

} // end namespace Dumux

#endif
