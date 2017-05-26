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
#ifndef itkThinShellDemonsMetric_h
#define itkThinShellDemonsMetric_h

#include "itkMeshToMeshMetric.h"
#include "itkCovariantVector.h"
#include "itkMesh.h"
#include "itkImage.h"

namespace itk
{
/** \class ThinShellDemonsMetric
 * \brief Computes the minimum distance between a moving point-set
 *  and a fixed point-set. A vector of minimum closest point distance is
 *  created for each point in the moving point-set.
 *  No correspondance is needed.
 *  For speed consideration, the point-set with the minimum number of points
 *  should be used as the moving point-set.
 *  If the number of points is high, the possibility of setting a distance map
 *  should improve the speed of the closest point computation.
 *
 *  Reference: "A Method for Registration of 3-D Shapes",
 *             IEEE PAMI, Vol 14, No. 2, February 1992
 *
 * \ingroup RegistrationMetrics
 * \ingroup ITKRegistrationCommon
 */
template< typename TFixedMesh, typename TMovingMesh,
          typename TDistanceMap =
            ::itk::Image< unsigned short, TMovingMesh::PointDimension > >
class ITK_TEMPLATE_EXPORT ThinShellDemonsMetric:
  public MeshToMeshMetric< TFixedMesh, TMovingMesh >
{
public:

  /** Standard class typedefs. */
  typedef ThinShellDemonsMetric                                Self;
  typedef MeshToMeshMetric< TFixedMesh, TMovingMesh > Superclass;

  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ThinShellDemonsMetric, Object);

  /** Types transferred from the base class. */
  typedef typename Superclass::TransformType           TransformType;
  typedef typename Superclass::TransformPointer        TransformPointer;
  typedef typename Superclass::TransformParametersType TransformParametersType;
  typedef typename Superclass::TransformJacobianType   TransformJacobianType;

  typedef typename Superclass::MeasureType                MeasureType;
  typedef typename Superclass::DerivativeType             DerivativeType;
  typedef typename Superclass::FixedMeshType          FixedMeshType;
  typedef typename Superclass::MovingMeshType         MovingMeshType;
  typedef typename Superclass::FixedMeshConstPointer  FixedMeshConstPointer;
  typedef typename Superclass::MovingMeshConstPointer MovingMeshConstPointer;

  typedef typename Superclass::FixedPointIterator     FixedPointIterator;
  typedef typename Superclass::FixedPointDataIterator FixedPointDataIterator;

  typedef typename Superclass::MovingPointIterator     MovingPointIterator;
  typedef typename Superclass::MovingPointDataIterator MovingPointDataIterator;

  typedef typename Superclass::InputPointType InputPointType;
  typedef itk::MapContainer<int, InputPointType> TargetMapType;

  /** Get the number of values, i.e. the number of points in the moving set. */
  unsigned int GetNumberOfValues() const ITK_OVERRIDE;

  /** Get the derivatives of the match measure. */
  void GetDerivative(const TransformParametersType & parameters,
                     DerivativeType & Derivative) const ITK_OVERRIDE;

  /**  Get the match measure, i.e. the value for single valued optimizers. */
  MeasureType GetValue(const TransformParametersType & parameters) const ITK_OVERRIDE;

  /**  Get value and derivatives for multiple valued optimizers. */
  void GetValueAndDerivative(const TransformParametersType & parameters,
                             MeasureType & Value, DerivativeType & Derivative) const;

  /** Initialize the Metric by making sure that all the components
   *  are present and plugged together correctly     */
  virtual void Initialize(void) throw ( ExceptionObject ) ITK_OVERRIDE;
protected:
  ThinShellDemonsMetric();
  virtual ~ThinShellDemonsMetric() {}

  void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(ThinShellDemonsMetric);

  bool               m_TargetPositionComputed;
  TargetMapType targetMap;
  void ComputeTargetPosition();
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkThinShellDemonsMetric.hxx"
#endif

#endif
