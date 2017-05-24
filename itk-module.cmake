set(DOCUMENTATION "An ITK-based implementation of Wavelet transforms using isotropic wavelets.

A more detailed description can be found in the Insight Journal article::

Cerdan, P.H. \"ITK Wavelet Module\".
  http://hdl.handle.net/10380/3558
  September, 2016.
")

# itk_module() defines the module dependencies in ExternalTemplate
# ExternalTemplate depends on ITKCommon
# The testing module in ExternalTemplate depends on ITKTestKernel
# and ITKMetaIO(besides ExternalTemplate and ITKCore)
# By convention those modules outside of ITK are not prefixed with
# ITK.

# define the dependencies of the include module and the tests
itk_module(ExternalTemplate
  DEPENDS
    ITKCommon
	ITKMesh
	ITKRegistrationCommon
  TEST_DEPENDS
    ITKTestKernel
	ITKVtkGlue
    ITKMetaIO
  DESCRIPTION
    "${DOCUMENTATION}"
  ${_EXCLUDE}
)
