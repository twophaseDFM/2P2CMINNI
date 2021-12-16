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
 * \brief Test for the two-phase two-component CC model.
 */
#include <config.h>

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>

#include <dumux/io/grid/gridmanager.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>

#include <dumux/linear/amgbackend.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/assembly/fvassembler.hh>

#include <dumux/io/vtkoutputmodule.hh>
//#include <dumux/io/grid/gridmanager_yasp.hh>

#include <dumux/io/loadsolution.hh>

#include <dumux/porousmediumflow/velocity.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    // define the type tag for this problem
    using TypeTag = Properties::TTag::TYPETAG;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // try to create a grid (from the given grid file or the input file)
    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);
    gridGeometry->update();

    // the problem (initial and boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    // get some time loop parameters
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    // check if we are about to restart a previously interrupted simulation
    Scalar restartTime = getParam<Scalar>("Restart.Time", 0);

    // the solution vector
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());
    if (restartTime > 0)
    {
        using IOFields = GetPropType<TypeTag, Properties::IOFields>;
        using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
        using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
        using SolidSystem = GetPropType<TypeTag, Properties::SolidSystem>;
        const auto fileName = getParam<std::string>("Restart.File");
        const auto pvName = createPVNameFunction<IOFields, PrimaryVariables, ModelTraits, FluidSystem, SolidSystem>();
        loadSolution(x, fileName, pvName, *gridGeometry);
    }
    else
        problem->applyInitialSolution(x);
    auto xOld = x;

    // initialize the phase saturation
//    Scalar deltaPorosity_ = 0;
//    std::vector<Scalar> deltaPorosity_(gridGeometry->numScv(), 0.0);
//    Scalar molarDensityH2O_ = getParam<Scalar>("Rock.molarDensityH2O");
//    Scalar molarDensityCaIon_ = getParam<Scalar>("Rock.molarDensityCaIon");
//    Scalar initialCaIon_ = getParam<Scalar>("Rock.initialCaIon");

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    // intialize the vtk output module
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    using VelocityOutput = GetPropType<TypeTag, Properties::VelocityOutput>;
    vtkWriter.addVelocityOutput(std::make_shared<VelocityOutput>(*gridVariables));
    IOFields::initOutputModule(vtkWriter); // Add model specific output fields
//    vtkWriter.addField(problem->getPermeability(), "Permeability");
//    vtkWriter.addField(problem->getPorosity(), "deltaPorosity");
//    vtkWriter.addField(problem->getCa(), "concentration_Ca"); // mol/m3
    problem->addVtkFields(vtkWriter);
    problem->updateVtkOutput(x);
    vtkWriter.write(restartTime);

    // instantiate time loop
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(restartTime, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    // the assembler with time loop for instationary problem
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, xOld);

    // the linear solver
    using LinearSolver = AMGBiCGSTABBackend<LinearSolverTraits<GridGeometry>>;
    auto linearSolver = std::make_shared<LinearSolver>(leafGridView, gridGeometry->dofMapper());

    // the non-linear solver
    using NewtonSolver = NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

//    using Velocity = Dune::FieldVector<double, 7>;
//    using VelocityVector = std::vector<Velocity>;
//    VelocityVector velocity_;
//    velocity_.resize(gridGeometry->gridView().size(codim));
//
//    using LocalResidual = GetPropType<TypeTag, Properties::LocalResidual>;
//    LocalResidual localResidual(&*problem, nullptr);
//
//    const Scalar elementX_ = getParam<Scalar>("Problem.elementX");
//    const auto& bBoxMaxX_ = gridGeometry->bBoxMax()[0];
//    const Scalar ProductionLowLimit_ = getParam<Scalar>("Problem.ProductionLowLimit");
//    const auto eps_ = 1e-6;
//    VelocityVector outputFlux(gridGeometry->gridView().size(0));

//    double volume = 0.0;

    // time loop
    timeLoop->start(); do
    {
        problem->setTime( timeLoop->time() + timeLoop->timeStepSize() );
        problem->setTimeStepSize( timeLoop->timeStepSize() );

        // solve the non-linear system with time step control
        nonLinearSolver.solve(x, *timeLoop);

        // make the new solution the old solution
        xOld = x;
        gridVariables->advanceTimeStep();

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // set new dt as suggested by the newton solver
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));

//        for (const auto& element : elements(leafGridView))
//        {
//        	volume += element.geometry().volume();
//
//        	auto fvGeometry = localView(*gridGeometry);
//        	fvGeometry.bind(element);
//
//        	auto elemVolVars = localView(gridVariables->curGridVolVars());
//        	elemVolVars.bind(element, fvGeometry, x);
//
//			auto elemFluxVarsCache = localView(gridVariables->gridFluxVarsCache());
//			elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);
//
//        	for (const auto& scvf : scvfs(fvGeometry))
//			{
//        		const auto& globalPos = scvf.ipGlobal();
//        		if (globalPos[1] > ProductionLowLimit_ - eps_ &&
//					globalPos[0] > bBoxMaxX_ - 2 * elementX_ - eps_ && globalPos[0] < bBoxMaxX_ - elementX_ - eps_)
//        		{
////        			outputFlux = localResidual.computeFlux(*problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
//        			calculateVelocity(velocity_, element, fvGeometry, elemVolVars, elemFluxVarsCache, 1);
//        		}
//			}
//        }
        // calculte output results
        problem->calculateOutput(*gridVariables, x, timeLoop->time());
        problem->updateVtkOutput(x);

//        problem->spatialParams().setdeltaPorosity(deltaPorosity_);

        // write vtk output
        vtkWriter.write(timeLoop->time());

    } while (!timeLoop->finished());

    timeLoop->finalize(leafGridView.comm());

    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////

    // print dumux end message
    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

    return 0;
} // end main
