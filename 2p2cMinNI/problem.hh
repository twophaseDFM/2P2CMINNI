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

#include <dumux/porousmediumflow/velocity.hh>

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
        nPhaseIdx = FluidSystem::CO2Idx,    // 1
		wPhaseIdx = FluidSystem::liquidPhaseIdx,    // 0
        H2OIdx = FluidSystem::BrineIdx,
//		H2OIdx = FluidSystem::H2OIdx,   // 0
        CO2Idx = FluidSystem::CO2Idx,   // 1
		CaIonIdx = FluidSystem::CaIonIdx,   // 2
		NaIonIdx = FluidSystem::NaIonIdx,  // 3

    };

    // Solid phase
    enum
	{
    	sCaCO3Idx = SolidSystem::comp0Idx,   // 0
		sNaClIdx = SolidSystem::comp1Idx,   // 1
		CaCO3Idx = Indices::conti0EqIdx + FluidSystem::numComponents,  // 4
		NaClIdx = Indices::conti0EqIdx + FluidSystem::numComponents + 1,   // 5
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
    using Velocity = Dune::FieldVector<double, 2>;
    using VelocityVector = std::vector<Velocity>;

    //! Property that defines whether mole or mass fractions are used
    static constexpr bool useMoles = ModelTraits::useMoles();

public:
    InjectionProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
    	name_              		 = getParam<std::string>("Problem.Name");
    	depthBOR_           	 = getParam<Scalar>("Problem.DepthBOR");

        InjectionTemperature_    = getParam<Scalar>("Injection.InjectionTemperature");
        InjectionPressure_ 		 = getParam<Scalar>("Injection.InjectionPressure");
        InjectionRate_           = getParam<Scalar>("Injection.InjectionRate");
        InjectionYMax_ 			 = getParam<Scalar>("Injection.InjectionYMax");
        InjectionYMin_ 			 = getParam<Scalar>("Injection.InjectionYMin");

        ProductionYMax_ 		 = getParam<Scalar>("Production.ProductionYMax");
        ProductionYMin_		 	 = getParam<Scalar>("Production.ProductionYMin");
        ProductionLowLimit_	 	 = getParam<Scalar>("Production.ProductionLowLimit");
        elementX_ 				 = getParam<Scalar>("Production.elementX");
        ProductionBoundary_ 	 = getParam<Scalar>("Production.ProductionBoundary");
        ProductionXMax_			 = getParam<Scalar>("Production.ProductionXMax");

        initSalinity_            = getParam<Scalar>("Initial.InitialSalinity");
        initialCaIon_ 			 = getParam<Scalar>("Initial.initialCaIon");
        initialCaCO3_			 = getParam<Scalar>("Initial.initialCaCO3");
        initialNaCl_			 = getParam<Scalar>("Initial.initialNaCl");

        porosity_		 	     = getParam<Scalar>("SpatialParams.porosity");
        // stating in the console whether mole or mass fractions are used
        if(useMoles)
            std::cout<<"problem uses mole-fractions"<<std::endl;
        else
            std::cout<<"problem uses mass-fractions"<<std::endl;

        unsigned int codim = GetPropType<TypeTag, Properties::GridGeometry>::discMethod == DiscretizationMethod::box ? dim : 0;

        permeability_.resize(gridGeometry->gridView().size(codim));
        deltaPermeability_.resize(gridGeometry->gridView().size(codim));
        deltaPorosity_.resize(gridGeometry->gridView().size(codim));
        Ca_.resize(gridGeometry->gridView().size(codim));
        Na_.resize(gridGeometry->gridView().size(codim));
//        velocity_.resize(gridGeometry->gridView().size(codim));
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

		if (globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_ && globalPos[1] > ProductionLowLimit_ - eps_)
			bcTypes.setAllDirichlet();
		else
			bcTypes.setAllNeumann();

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

//        values[pressureIdx] += - 5e6;
//		values[switchIdx] = 0.0;
//        values[CaIonIdx] = initialCaIon_ ;
//        values[CaCO3Idx] = 0.0;

    	return values;
    }

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


		if (isInjection(globalPos))
		{
			values[contiCO2EqIdx] = -InjectionRate_/FluidSystem::molarMass(CO2Idx); //mole/(m^2*s) -> kg/(s*m^2)
			values[Indices::temperatureIdx] = values[contiCO2EqIdx] * FluidSystem::componentEnthalpy(fs, nPhaseIdx, CO2Idx) * FluidSystem::molarMass(CO2Idx);
			values[CaIonIdx] = - values[contiCO2EqIdx] * volVars.moleFraction(nPhaseIdx, H2OIdx) * volVars.moleFraction(wPhaseIdx,CaIonIdx);
			values[NaIonIdx] = - values[contiCO2EqIdx] * volVars.moleFraction(nPhaseIdx, H2OIdx) * volVars.moleFraction(wPhaseIdx,NaIonIdx);
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
//        const auto& globalPos = scv.center();

		/* *****************  REACTION CaIon = CaCO3 ***********************/
        // TODO determine the kinetics of the reaction
		Scalar R = getParam<Scalar>("Reaction.R");
        Scalar Ca_k0 = getParam<Scalar>("Reaction.Ca_k0");
        Scalar Ca_E = getParam<Scalar>("Reaction.Ca_E");
        Scalar Ca_Max = getParam<Scalar>("Reaction.Ca_Max");

        Scalar Ca_K = Ca_k0 * exp(-Ca_E/R/volVars.temperature()) * volVars.saturation(wPhaseIdx) * volVars.moleFraction(wPhaseIdx,CO2Idx) * volVars.pressure(CO2Idx);
        Scalar Ca_r = Ca_K * (1 - volVars.moleFraction(wPhaseIdx, CaIonIdx)/Ca_Max);

        source[CO2Idx] += -Ca_r;
		source[CaIonIdx] += Ca_r;
		source[CaCO3Idx] += -Ca_r;

		/* *****************  REACTION NaIon = NaCl ***********************/
        Scalar Na_k0 = getParam<Scalar>("Reaction.Na_k0");
		Scalar Na_E = getParam<Scalar>("Reaction.Na_E");
		Scalar Na_Max = getParam<Scalar>("Reaction.Na_Max");

        using std::abs;
        Scalar Na_K = Na_k0 * exp(-Na_E/R/volVars.temperature());
        Scalar Na_r = Na_K * (1 - volVars.moleFraction(wPhaseIdx, NaIonIdx)/Na_Max);  //disolution rate, positive means disolution

        if (volVars.solidVolumeFraction(sNaClIdx) < 0 && Na_r > 0)
        	Na_r = 0;

//        std::cout << " solidVolumeFraction = " << volVars.solidVolumeFraction(sNaClIdx) <<'\n' ;
//        std::cout << " Na_r = " << Na_r <<'\n' ;

        source[NaIonIdx] += Na_r;
        source[NaClIdx] += -Na_r;
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
        Scalar index_p = 0, index_pf = 0, index_if = 0;
        Scalar outputSn = 0;
        Scalar outputCO2 = 0;
        Scalar outputCaIon = 0, outputNaIon = 0;
        Scalar outputT = 0, outputOverPw = 0, outputOverPn = 0;
        Scalar flux_nw = 0, flux_w = 0, flux_t = 0, flux_inj = 0;
        Scalar Area_inj = 0;

        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
			auto fvGeometry = localView(this->gridGeometry());
			fvGeometry.bind(element);
			auto elemVolVars = localView(gridVariables.curGridVolVars());
			elemVolVars.bind(element, fvGeometry, sol);
	        auto elemFluxVarsCache = localView(gridVariables.gridFluxVarsCache());
	        elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);

			for (auto&& scv : scvs(fvGeometry))
			{
				const auto& globalPos = scv.center();
				const auto& volVars = elemVolVars[scv];
				Source = this->source(element, fvGeometry, elemVolVars, scv);
				InitialValues = initialAtPos(globalPos);

				if (isProduction(globalPos))
				{
					index_p += 1;
					outputSn += volVars.saturation(nPhaseIdx);
					outputCO2 += volVars.saturation(nPhaseIdx) * volVars.moleFraction(nPhaseIdx, CO2Idx) + volVars.saturation(wPhaseIdx) * volVars.moleFraction(wPhaseIdx, CO2Idx);
					outputCaIon += volVars.saturation(wPhaseIdx) * volVars.moleFraction(wPhaseIdx, FluidSystem::CaIonIdx);
					outputNaIon += volVars.saturation(wPhaseIdx) * volVars.moleFraction(wPhaseIdx, FluidSystem::NaIonIdx);
					outputT += volVars.temperature();
					outputOverPw += volVars.pressure(wPhaseIdx) - InitialValues[pressureIdx];
					outputOverPn += volVars.pressure(nPhaseIdx) - InitialValues[pressureIdx];
				}
			}

/* ************** Calculate production Flux  ****************** */
			for (const auto& scvf : scvfs(fvGeometry))
			{
				if (scvf.boundary())
					continue;

				const auto& globalPos = scvf.ipGlobal();
				const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
                const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
                const bool isOnProduction = insideScv.center()[0] < ProductionBoundary_ && outsideScv.center()[0] > ProductionBoundary_ && globalPos[1] > ProductionLowLimit_ - eps_;
				using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
				FluxVariables fluxVars;
				fluxVars.init(*this, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
				if (isOnProduction)
				{
					const auto gasPhaseIdx = GetPropType<TypeTag, Properties::FluidSystem>::comp1Idx;
					const auto liquidPhaseIdx = GetPropType<TypeTag, Properties::FluidSystem>::comp0Idx;
					auto upwindTerm_nw = [gasPhaseIdx] (const auto& volVars)
					{ return volVars.mobility(gasPhaseIdx)*volVars.density(gasPhaseIdx); };
					auto upwindTerm_w = [liquidPhaseIdx] (const auto& volVars)
					{ return volVars.mobility(liquidPhaseIdx)*volVars.density(liquidPhaseIdx); };

					index_pf += 1;
					flux_nw += (fluxVars.advectiveFlux(gasPhaseIdx, upwindTerm_nw) * scvf.area());   // unit: mole/s
					flux_w += (fluxVars.advectiveFlux(liquidPhaseIdx, upwindTerm_w) * scvf.area());
				}
			}

		/* *************** Calculate Injection Flux ******************/
			for (const auto& scvf: scvfs(fvGeometry))
			{
				const auto& globalPos = scvf.ipGlobal();
				if (isInjection(globalPos))
				{
					if (scvf.boundary())
					 {
						flux_inj += (- neumann(element, fvGeometry, elemVolVars, elemFluxVarsCache, scvf)[contiCO2EqIdx]*scvf.area());
						Area_inj += scvf.area();
						index_if += 1;
					 }
				}
			}
        }

        outputSn /= index_p;
        outputCO2 /= index_p;
        outputCaIon /= index_p;
        outputNaIon /= index_p;
        outputT /= index_p;
        outputOverPw /= 1e6;   // transfer to MPa
        outputOverPn /= 1e6;
        flux_t = flux_w + flux_nw;

        std::cout << "outputSn = " << outputSn << " outputCO2 = " << outputCO2 << '\n';
        std::cout << "outputCaIon = " << outputCaIon << " outputNaIon = " << outputNaIon << '\n';
        std::cout << "outputT = " << outputT <<  " outputOverPw = " << outputOverPw << " outputOverPn = " << outputOverPn << '\n';
        std::cout << "flux_nw = " << flux_nw << " flux_w = " << flux_w << " flux_t = " << flux_t << '\n';
        std::cout << "flux_injection = " << flux_inj << '\n';
        std::cout << " Area_injection = " << Area_inj << '\n';

//        fileSn_.open("1outputSn_3D.log",std::ios::app);
//        fileSn_ << time << " " << OutputSn << " " << InjectionSn << " " << OutputE << " " << OutputPressure << " " << InjectionPressure << " " << totalMassCO2 << std::endl;
//        fileSn_.close();
    }

    template<class VTKWriter>
    void addVtkFields(VTKWriter& vtk)
    {
    	vtk.addField(permeability_, "permeability");
    	vtk.addField(deltaPermeability_, "deltaPermeability");
    	vtk.addField(deltaPorosity_, "deltaPorosity");
    	vtk.addField(Ca_,"concentration_Ca2+");
    	vtk.addField(Na_,"concentration_Na2+");
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
                using Dune::power;
                permeability_[dofIdxGlobal] = volVars.permeability();
                deltaPermeability_[dofIdxGlobal] = 1 / power((1 - porosity_)/(1 - volVars.porosity()), 2) * power(volVars.porosity()/porosity_ , 3);
                deltaPorosity_[dofIdxGlobal] = (initialCaCO3_ + initialNaCl_) - (volVars.solidVolumeFraction(sNaClIdx) + volVars.solidVolumeFraction(sCaCO3Idx));
                Ca_[dofIdxGlobal] = volVars.moleFraction(wPhaseIdx, CaIonIdx) * volVars.saturation(wPhaseIdx) * volVars.moleFraction(wPhaseIdx, H2OIdx) * scv.volume();
                Na_[dofIdxGlobal] = volVars.moleFraction(wPhaseIdx, NaIonIdx) * volVars.saturation(wPhaseIdx) * volVars.moleFraction(wPhaseIdx, H2OIdx) * scv.volume();

            }
        }
    }

private:
    /*!
     * \brief Returns the molality of NaCl (mol NaCl / kg water) for a given mole fraction.
     *
     * \param XwNaCl The XwNaCl [kg NaCl / kg solution]
     */
    static Scalar massToMoleFrac_(Scalar XwNaCl)
    {
       const Scalar Mw = 18.015e-3; //FluidSystem::molarMass(H2OIdx); /* molecular weight of water [kg/mol] */
       const Scalar Ms = 58.44e-3;  //FluidSystem::molarMass(NaClIdx); /* molecular weight of NaCl  [kg/mol] */

       const Scalar X_NaCl = XwNaCl;
       /* XwNaCl: conversion from mass fraction to mol fraction */
       auto xwNaCl = -Mw * X_NaCl / ((Ms - Mw) * X_NaCl - Ms);
       return xwNaCl;
    }

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
		priVars[switchIdx] = 0.0;
        priVars[pressureIdx] = pl;
        priVars[Indices::temperatureIdx] = 413.15;
        priVars[CaIonIdx] = initialCaIon_ ;
        priVars[CaCO3Idx] = initialCaCO3_;
        priVars[NaIonIdx] = massToMoleFrac_(initSalinity_);
        priVars[NaClIdx] = initialNaCl_; // [kg/m^3]
        return priVars;
    }

    bool isProduction (const GlobalPosition globalPos) const
    {
    	return globalPos[1] < ProductionYMax_ + eps_ &&
			   globalPos[1] > ProductionYMin_ - eps_ &&
			   globalPos[0] < ProductionXMax_ + eps_ &&
			   globalPos[0] > ProductionXMax_ - elementX_ - eps_;
    }

    bool isInjection (const GlobalPosition globalPos) const
    {
    	return globalPos[1] < InjectionYMax_ &&
    		   globalPos[1] > InjectionYMin_ &&
			   globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_;
    }

    Scalar depthBOR_, InjectionYMax_, InjectionYMin_, ProductionYMax_, ProductionYMin_, ProductionLowLimit_;
    Scalar InjectionRate_, InjectionTemperature_, InjectionPressure_;
    Scalar time_ = 0.0, timeStepSize_ = 0.0;
	Scalar porosity_;
	Scalar elementX_, ProductionBoundary_ , ProductionXMax_;
	Scalar initialCaIon_, initialCaCO3_ , initialNaCl_;
	Scalar initSalinity_ ;
	std::vector<Scalar> permeability_, deltaPorosity_, deltaPermeability_, Ca_, Na_ ;

	std::string name_;
	static constexpr Scalar eps_ = 1e-6;
};

} // end namespace Dumux

#endif
