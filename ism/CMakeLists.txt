SET(ISM_SRCS
	ContourCriticalPointsGSL.cpp
	ContourCriticalPointsGSL.h
	ContourCritical.cpp
	ContourCritical.h
	contractors.cpp
	contractors.h
	CriticalMesh.cpp
	CriticalMesh.h
	GaussJordan.cpp
	GaussJordan.h
	Jacobi.cpp
	Jacobi.h
	Newton.cpp
	Newton.h
	NewtonClassify.cpp
	NewtonClassify.h
	NewtonCritical.cpp
	NewtonCritical.h
	NewtonDegenerate.cpp
	NewtonDegenerate.h
	SearchCritical.cpp
	SearchCritical.h
	SearchDegenerate.cpp
	SearchDegenerate.h
	SingularityClassification.cpp
	SingularityClassification.h
)

INCLUDE_DIRECTORIES(${BOOST_INCLUDE_PATH} ${GPS_SOURCE_DIR})

ADD_LIBRARY(ism
  ${ISM_SRCS}
  )
 