set(DOCUMENTATION "An ITK-based implementation of Wavelet transforms using isotropic wavelets.

A more detailed description can be found in the Insight Journal article::

Cerdan, P.H. \"ITK Wavelet Module\".
  http://hdl.handle.net/10380/3558
  September, 2016.
")


# define the dependencies of the include module and the tests
itk_module(ExternalTemplate
  ENABLE_SHARED
  COMPILE_DEPENDS
    ITKCommon
	  ITKMesh
	  ITKRegistrationCommon
  	ITKVtkGlue
  TEST_DEPENDS
    ITKTestKernel
    ITKMetaIO
  DESCRIPTION
    "${DOCUMENTATION}"
  ${_EXCLUDE}
)
