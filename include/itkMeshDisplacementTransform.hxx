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
#ifndef itkMeshDisplacementTransform_hxx
#define itkMeshDisplacementTransform_hxx

#include "itkMeshDisplacementTransform.h"
#include "itkMath.h"
#include "itkMath.h"

namespace itk
{

template<typename TParametersValueType, unsigned int NDimensions>
MeshDisplacementTransform<TParametersValueType, NDimensions>
	::MeshDisplacementTransform()
{
	m_MeshTemplate = ITK_NULLPTR;
	this->SpaceDimension = NDimensions;
	this->ParametersDimension = 0;
}


template<typename TParametersValueType, unsigned int NDimensions>
MeshDisplacementTransform<TParametersValueType, NDimensions>
::~MeshDisplacementTransform()
{
}


template<typename TParametersValueType, unsigned int NDimensions>
void
MeshDisplacementTransform<TParametersValueType, NDimensions>
::SetParameters(const ParametersType & parameters)
{
  // Save parameters. Needed for proper operation of TransformUpdateParameters.
	if( &parameters != &(this->m_Parameters) )
	{
		this->m_Parameters = parameters;
	}

	MeshDeformationType::PixelType  deformationVector;
	unsigned int numberOfPoints = ParametersDimension/SpaceDimension;
	bool modified = false;

	for (unsigned int i = 0; i < numberOfPoints; i++ ){
		deformationVector = this->m_MeshDeformation->GetPointData(i);
		for (unsigned int j = 0; j < SpaceDimension; j++)
			if( Math::NotExactlyEquals(deformationVector[j], parameters[i*SpaceDimension+j]) )
			{
				deformationVector[j] = parameters[i*SpaceDimension+j];
				this->m_MeshDeformation->SetPointData(i,deformationVector);
				modified = true;
			}
	}
	if( modified ){
		this->Modified();
	}
}

template<typename TParametersValueType, unsigned int NDimensions>
const typename MeshDisplacementTransform<TParametersValueType, NDimensions>::ParametersType &
MeshDisplacementTransform<TParametersValueType, NDimensions>
::GetParameters() const
{
	MeshDeformationType::PixelType  deformationVector;
	unsigned int numberOfPoints = ParametersDimension/SpaceDimension;

	for (unsigned int i = 0; i < numberOfPoints; i++ ){
		deformationVector = this->m_MeshDeformation->GetPointData(i);
		for (unsigned int j = 0; j < SpaceDimension; j++)
			this->m_Parameters[i*SpaceDimension + j] = deformationVector[j];
	}

	return this->m_Parameters;
}

template<typename TParametersValueType, unsigned int NDimensions>
void
	MeshDisplacementTransform<TParametersValueType, NDimensions>
	::SetIdentity()
{
	if (!m_Mesh)
	{
		itkExceptionMacro(<< "Mesh template is not present");
	}
	if (ParametersDimension == 0)
	{
		itkExceptionMacro(<< "Mesh template has zero vertex");
	}

	MeshDeformationType::PixelType  deformationVector;
	MeshDeformationType::PointType  point;

	PointIterator pointItr = this->GetMeshTemplate()->GetPoints()->Begin();
	PointIterator pointEnd = this->GetMeshTemplate()->GetPoints()->End();

	unsigned int pointId =  0;
	while ( pointItr != pointEnd )
	{
		deformationVector[0] =  0;
		deformationVector[1] =  0;
		deformationVector[2] =  0;  

		typename Superclass::InputPointType inputPoint;
		inputPoint.CastFrom( pointItr.Value() );

		m_MeshDeformation->SetPoint( pointId, inputPoint );
		m_MeshDeformation->SetPointData( pointId, deformationVector );
		pointId++;
	}
}

template<typename TParametersValueType, unsigned int NDimensions>
void
	MeshDisplacementTransform<TParametersValueType, NDimensions>
	::Initialize()
{
	if (!m_MeshTemplate)
	{
		itkExceptionMacro(<< "FixedMesh is not present");
	}

	m_MeshDeformation = MeshDeformationType::New();
	MeshDeformationType::PixelType  deformationVector;
	MeshDeformationType::PointType  point;

	PointIterator pointItr = this->GetMeshTemplate()->GetPoints()->Begin();
	PointIterator pointEnd = this->GetMeshTemplate()->GetPoints()->End();

	unsigned int pointId =  0;
	while ( pointItr != pointEnd )
	{
		deformationVector[0] =  0;
		deformationVector[1] =  0;
		deformationVector[2] =  0;  

		typename Superclass::InputPointType inputPoint;
		inputPoint.CastFrom( pointItr.Value() );

		m_MeshDeformation->SetPoint( pointId, inputPoint );
		m_MeshDeformation->SetPointData( pointId, deformationVector );
		pointId++;
	}

	this->m_Parameters.SetParametersObject(m_MeshDeformation->GetPointData());
	this->ParametersDimension = pointId * this->SpaceDimension;

	//Superclass()
	//itkStaticConstMacro(SpaceDimension, unsigned int, NDimensions);
	//itkStaticConstMacro(ParametersDimension, unsigned int, NDimensions);
}

template<typename TParametersValueType, unsigned int NDimensions>
void
MeshDisplacementTransform<TParametersValueType, NDimensions>
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

//  os << indent << "Offset: " << m_Offset << std::endl;
}


// template<typename TParametersValueType, unsigned int NDimensions>
// void
// MeshDisplacementTransform<TParametersValueType, NDimensions>
// ::Compose(const Self *other, bool)
// {
//   this->Translate(other->m_Offset);
// }


// template<typename TParametersValueType, unsigned int NDimensions>
// void
// MeshDisplacementTransform<TParametersValueType, NDimensions>
// ::Translate(const OutputVectorType & offset, bool)
// {
//   ParametersType newOffset(SpaceDimension);
// 
//   for( unsigned int i = 0; i < SpaceDimension; i++ )
//     {
//     newOffset[i] = m_Offset[i] + offset[i];
//     }
//   this->SetParameters(newOffset);
// }


template<typename TParametersValueType, unsigned int NDimensions>
typename MeshDisplacementTransform<TParametersValueType, NDimensions>::OutputPointType
MeshDisplacementTransform<TParametersValueType, NDimensions>
::TransformPoint(const InputPointType & point) const
{
  return point;
}


template<typename TParametersValueType, unsigned int NDimensions>
typename MeshDisplacementTransform<TParametersValueType, NDimensions>::OutputVectorType
MeshDisplacementTransform<TParametersValueType, NDimensions>
::TransformVector(const InputVectorType & vect) const
{
  return vect;
}


template<typename TParametersValueType, unsigned int NDimensions>
typename MeshDisplacementTransform<TParametersValueType, NDimensions>::OutputVnlVectorType
MeshDisplacementTransform<TParametersValueType, NDimensions>
::TransformVector(const InputVnlVectorType & vect) const
{
  return vect;
}


template<typename TParametersValueType, unsigned int NDimensions>
typename MeshDisplacementTransform<TParametersValueType, NDimensions>::OutputCovariantVectorType
MeshDisplacementTransform<TParametersValueType, NDimensions>
::TransformCovariantVector(const InputCovariantVectorType & vect) const
{
  return vect;
}


// template<typename TParametersValueType, unsigned int NDimensions>
// bool
// MeshDisplacementTransform<TParametersValueType, NDimensions>
// ::GetInverse(Self *inverse) const
// {
//   if( !inverse )
//     {
//     return false;
//     }
// 
//   inverse->SetFixedParameters(this->GetFixedParameters());
//   inverse->m_Offset   = -m_Offset;
//   return true;
// }


template<typename TParametersValueType, unsigned int NDimensions>
typename MeshDisplacementTransform<TParametersValueType, NDimensions>::InverseTransformBasePointer
MeshDisplacementTransform<TParametersValueType, NDimensions>
::GetInverseTransform() const
{
  Pointer inv = New();

  //return GetInverse(inv) ? inv.GetPointer() : ITK_NULLPTR;
  return inv;
}


template<typename TParametersValueType, unsigned int NDimensions>
void
MeshDisplacementTransform<TParametersValueType, NDimensions>::ComputeJacobianWithRespectToParameters(
  const InputPointType &,
  JacobianType & jacobian) const
{
  // the Jacobian is constant for this transform, and it has already been
  // initialized in the constructor, so we just need to return it here.
  jacobian = this->m_IdentityJacobian;
}


template<typename TParametersValueType, unsigned int NDimensions>
void
MeshDisplacementTransform<TParametersValueType, NDimensions>
::ComputeJacobianWithRespectToPosition(const InputPointType &,
                                       JacobianType & jac) const
{
  jac.SetSize( NDimensions, NDimensions );
  jac.Fill(0.0);
  for( unsigned int dim = 0; dim < NDimensions; dim++ )
    {
    jac[dim][dim] = 1.0;
    }
}


} // end namespace itk

#endif
