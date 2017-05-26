/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#include <cstdlib>

#include "itkVTKPolyDataReader.h"
#include "itkVTKPolyDataWriter.h"
#include "itkThinShellDemonsMetric.h"
#include "itkTranslationTransform.h"
#include "itkLevenbergMarquardtOptimizer.h"
#include "itkMeshToMeshRegistrationMethod.h"
#include "itkMeshDisplacementTransform.h"

class CommandIterationUpdate : public itk::Command
{
public:
	typedef  CommandIterationUpdate   Self;
	typedef  itk::Command             Superclass;
	typedef itk::SmartPointer<Self>   Pointer;
	itkNewMacro( Self );
protected:
	CommandIterationUpdate() {};
public:
	typedef itk::LevenbergMarquardtOptimizer     OptimizerType;
	typedef const OptimizerType *                OptimizerPointer;
	void Execute(itk::Object *caller, const itk::EventObject & event)
	{
		Execute( (const itk::Object *)caller, event);
	}
	void Execute(const itk::Object * object, const itk::EventObject & event)
	{
		OptimizerPointer optimizer = static_cast< OptimizerPointer >( object );
		if( ! itk::IterationEvent().CheckEvent( &event ) )
		{
			return;
		}
		//std::cout << "Value = " << optimizer->GetCachedValue() << std::endl;
		std::cout << "Position = "  << optimizer->GetCachedCurrentPosition();
		std::cout << std::endl << std::endl;
	}
};

int itkEmptyTest( int , char * [])
{

	const unsigned int Dimension = 3;
	typedef itk::Mesh<double, Dimension>   MeshType;
	typedef itk::VTKPolyDataReader< MeshType >  ReaderType;

	/*
	Initialize fixed mesh polydata reader
	*/
	ReaderType::Pointer  fixedPolyDataReader = ReaderType::New();

	typedef ReaderType::PointType   PointType;

	fixedPolyDataReader->SetFileName("fixedMesh.vtk");

	try
	{
		fixedPolyDataReader->Update();
	}
	catch( itk::ExceptionObject & excp )
	{
		std::cerr << "Error during Fixed Mesh Update() " << std::endl;
		std::cerr << excp << std::endl;
		return EXIT_FAILURE;
	}

	std::cout << "fixedPolyDataReader:" << std::endl;
	std::cout << fixedPolyDataReader << std::endl;

	MeshType::Pointer fixedMesh = fixedPolyDataReader->GetOutput();

	/*
	Initialize moving mesh polydata reader
	*/
	ReaderType::Pointer  movingPolyDataReader = ReaderType::New();

	typedef ReaderType::PointType   PointType;

	movingPolyDataReader->SetFileName("movingMesh.vtk");

	try
	{
		movingPolyDataReader->Update();
	}
	catch( itk::ExceptionObject & excp )
	{
		std::cerr << "Error during Moving Mesh Update() " << std::endl;
		std::cerr << excp << std::endl;
		return EXIT_FAILURE;
	}

	std::cout << "movingPolyDataReader:" << std::endl;
	std::cout << movingPolyDataReader << std::endl;

	MeshType::Pointer movingMesh = movingPolyDataReader->GetOutput();

	/*
		Initialize Thin Shell Demons metric
	*/
	typedef itk::ThinShellDemonsMetric<
		MeshType,
		MeshType>
		MetricType;
	MetricType::Pointer  metric = MetricType::New();

	/*
		Initialize Thin Shell Demons transformation
	*/
	typedef itk::TranslationTransform< double, Dimension >      TransformType;
	TransformType::Pointer transform = TransformType::New();

	typedef itk::MeshDisplacementTransform< double, Dimension >    TransformTestType;
	TransformTestType::Pointer transformTest = TransformTestType::New();
	
	transformTest->SetMeshTemplate(movingMesh);
	transformTest->Initialize();
	std::cout<<transformTest->GetNumberOfParameters()<<std::endl;
	/*
		Initialize Thin Shell Demons optimizer
	*/
	typedef itk::LevenbergMarquardtOptimizer OptimizerType;
	OptimizerType::Pointer      optimizer     = OptimizerType::New();
	optimizer->SetUseCostFunctionGradient(false);
	typedef itk::MeshToMeshRegistrationMethod<
		MeshType,
		MeshType >
		RegistrationType;
	// Scale the translation components of the Transform in the Optimizer
	OptimizerType::ScalesType scales( transformTest->GetNumberOfParameters() );
	scales.Fill( 0.01 );
	unsigned long   numberOfIterations =  100;
	double          gradientTolerance  =  1e-5;    // convergence criterion
	double          valueTolerance     =  1e-5;    // convergence criterion
	double          epsilonFunction    =  1e-6;   // convergence criterion
	optimizer->SetScales( scales );
	optimizer->SetNumberOfIterations( numberOfIterations );
	optimizer->SetValueTolerance( valueTolerance );
	optimizer->SetGradientTolerance( gradientTolerance );
	optimizer->SetEpsilonFunction( epsilonFunction );
	// Start from an Identity transform (in a normal case, the user
	// can probably provide a better guess than the identity...
	transformTest->SetIdentity();

	/*
		Initialize registration
	*/
	typedef itk::MeshToMeshRegistrationMethod<
		MeshType,
		MeshType >
		RegistrationType;
	RegistrationType::Pointer   registration  = RegistrationType::New();

	registration->SetMetric(        metric        );
	registration->SetOptimizer(     optimizer     );
	registration->SetTransform(     transformTest     );
	registration->SetInitialTransformParameters( transformTest->GetParameters() );
	registration->SetFixedMesh( fixedMesh );
	registration->SetMovingMesh( movingMesh );

	// Connect an observer
	CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
	optimizer->AddObserver( itk::IterationEvent(), observer );
	try
	{
		registration->Update();
	}
	catch( itk::ExceptionObject & e )
	{
		std::cout << e << std::endl;
		return EXIT_FAILURE;
	}
	std::cout << "Solution = " << transform->GetParameters() << std::endl;

	/*
		output mesh
	*/
	typedef itk::VTKPolyDataWriter<MeshType>   WriterType;
	WriterType::Pointer writer = WriterType::New();
	MeshType::ConstPointer registeredMesh = registration->GetMovingMesh();
	writer->SetInput( registeredMesh );
	writer->SetFileName( "registeredMesh.vtk" );
	writer->Write();
	return EXIT_SUCCESS;
}
