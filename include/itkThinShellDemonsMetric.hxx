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
#ifndef itkThinShellDemonsMetric_hxx
#define itkThinShellDemonsMetric_hxx

#include "itkThinShellDemonsMetric.h"
#include "itkImageRegionConstIteratorWithIndex.h"

namespace itk
{

template< typename TFixedMesh, typename TMovingMesh, typename TDistanceMap >
ThinShellDemonsMetric< TFixedMesh, TMovingMesh, TDistanceMap >
::ThinShellDemonsMetric() :
  m_ComputeSquaredDistance( false )
{

}

template< typename TFixedMesh, typename TMovingMesh, typename TDistanceMap >
unsigned int
ThinShellDemonsMetric< TFixedMesh, TMovingMesh, TDistanceMap >
::GetNumberOfValues() const
{
  MovingMeshConstPointer movingMesh = this->GetMovingMesh();

  if ( !movingMesh )
    {
    itkExceptionMacro(<< "Moving point set has not been assigned");
    }

  return movingMesh->GetPoints()->Size();
}

template< typename TFixedMesh, typename TMovingMesh, typename TDistanceMap >
typename ThinShellDemonsMetric< TFixedMesh, TMovingMesh, TDistanceMap >::MeasureType
ThinShellDemonsMetric< TFixedMesh, TMovingMesh, TDistanceMap >
::GetValue(const TransformParametersType & parameters) const
{
  FixedMeshConstPointer fixedMesh = this->GetFixedMesh();

  if ( !fixedMesh )
    {
    itkExceptionMacro(<< "Fixed point set has not been assigned");
    }

  MovingMeshConstPointer movingMesh = this->GetMovingMesh();

  if ( !movingMesh )
    {
    itkExceptionMacro(<< "Moving point set has not been assigned");
    }

  MovingPointIterator pointItr = movingMesh->GetPoints()->Begin();
  MovingPointIterator pointEnd = movingMesh->GetPoints()->End();

  MeasureType measure;
  measure.set_size( movingMesh->GetPoints()->Size() );

  this->SetTransformParameters(parameters);

  unsigned int identifier = 0;
  while ( pointItr != pointEnd )
    {
    typename Superclass::InputPointType inputPoint;
    inputPoint.CastFrom( pointItr.Value() );
    typename Superclass::OutputPointType transformedPoint =
      this->m_Transform->TransformPoint(inputPoint);

    double minimumDistance = NumericTraits< double >::max();
    bool   closestPoint = false;

    // Try to use the distance map to solve the closest point
    if ( m_DistanceMap )
      {
      // Check if the point is inside the distance map
      typename DistanceMapType::IndexType index;
      if ( m_DistanceMap->TransformPhysicalPointToIndex(transformedPoint, index) )
        {
        minimumDistance = m_DistanceMap->GetPixel(index);
        // In case the provided distance map was signed,
        // we correct here the distance to take its absolute value.
        if ( minimumDistance < 0.0 )
          {
          minimumDistance = -minimumDistance;
          }
        closestPoint = true;
        }
      }

    // If the closestPoint has not been found, go through the list of fixed
    // points and find the closest distance
    if ( !closestPoint )
      {
      FixedPointIterator pointItr2 = fixedMesh->GetPoints()->Begin();
      FixedPointIterator pointEnd2 = fixedMesh->GetPoints()->End();

      while ( pointItr2 != pointEnd2 )
        {
        double dist = pointItr2.Value().SquaredEuclideanDistanceTo(transformedPoint);

        if ( !m_ComputeSquaredDistance )
          {
          dist = std::sqrt(dist);
          }

        if ( dist < minimumDistance )
          {
          minimumDistance = dist;
          }
        pointItr2++;
        }
      }

    measure.put(identifier, minimumDistance);

    ++pointItr;
    identifier++;
    }

  return measure;
}

template< typename TFixedMesh, typename TMovingMesh, typename TDistanceMap >
void
ThinShellDemonsMetric< TFixedMesh, TMovingMesh, TDistanceMap >
::GetDerivative( const TransformParametersType & itkNotUsed(parameters),
                 DerivativeType & itkNotUsed(derivative) ) const
{}

template< typename TFixedMesh, typename TMovingMesh, typename TDistanceMap >
void
ThinShellDemonsMetric< TFixedMesh, TMovingMesh, TDistanceMap >
::GetValueAndDerivative(const TransformParametersType & parameters,
                        MeasureType & value, DerivativeType  & derivative) const
{
  value = this->GetValue(parameters);
  this->GetDerivative(parameters, derivative);
}

template< typename TFixedMesh, typename TMovingMesh, typename TDistanceMap >
void
ThinShellDemonsMetric< TFixedMesh, TMovingMesh, TDistanceMap >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "DistanceMap: " << m_DistanceMap << std::endl;
  if ( m_ComputeSquaredDistance )
    {
    os << indent << "m_ComputeSquaredDistance: True" << std::endl;
    }
  else
    {
    os << indent << "m_ComputeSquaredDistance: False" << std::endl;
    }
}
} // end namespace itk

#endif
