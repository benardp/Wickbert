# This macro merges elements in sorted lists ALIST and BLIST and stored the result in OUTPUT
MACRO(MERGE ALIST BLIST OUTPUT)
SET(BTEMP ${BLIST})
FOREACH(A ${ALIST})
	SET(SORTED)
	SET(UNINSERTED 1)
	FOREACH(B ${BTEMP})
		IF(${UNINSERTED})
			IF(${A} STRLESS ${B})
				SET(SORTED ${SORTED} ${A})
				SET(UNINSERTED 0)
			ENDIF(${A} STRLESS ${B})
		ENDIF(${UNINSERTED})
		SET(SORTED ${SORTED} ${B})
	ENDFOREACH(B ${BLIST})
	IF(${UNINSERTED})
		SET(SORTED ${SORTED} ${A})
	ENDIF(${UNINSERTED})
	SET(BTEMP ${SORTED})
ENDFOREACH(A ${ALIST})
SET(${OUTPUT} ${BTEMP})
ENDMACRO(MERGE ALIST BLIST OUTPUT)

MACRO(FOLDERFILES DIR NAME)
	FILE(GLOB DUD_SRCS ${DIR}/*.cpp)
	FILE(GLOB DUD_H ${DIR}/*.h)
	MERGE("${DUD_H}" "${DUD_SRCS}" DUD_SRCS)
	SOURCE_GROUP(${NAME} FILES ${DUD_SRCS})
	SET(SURFACE_SRCS ${SURFACE_SRCS} ${DUD_SRCS})
ENDMACRO(FOLDERFILES DIR NAME)


FILE(GLOB SURFACE_SRCS
	Surface.cpp Surface.h
	SurfParam.cpp SurfParam.h
	Interval.h itoa.cpp itoa.h
	Box.cpp Box.h
	IMatrix.cpp Imatrix.h
	KDTree.h
	gmTNTconvert.cpp gmTNTconvert.h
	RandomStream.cpp RandomStream.h
	ThinPlateSpline.cpp ThinPlateSpline.h
)
SOURCE_GROUP(Surface FILES ${SURFACE_SRCS})

FOLDERFILES("OpenMesh" "OpenMesh")

FILE(GLOB IMPLICIT_SRCS
	Implicit/Icosahedral.cpp Implicit/Icosahedral.h
	Implicit/Implicit.cpp Implicit/Implicit.h
	Implicit/ImpFileManager.cpp Implicit/ImpFileManager.h
	Implicit/ITKImplicit.cpp Implicit/ITKImplicit.h
	Implicit/VTKImplicit.cpp Implicit/VTKImplicit.h
	Implicit/polygonizer.cpp Implicit/polygonizer.h
	Implicit/pixmaps.h
)
SOURCE_GROUP("Implicit" FILES ${IMPLICIT_SRCS})
SET(SURFACE_SRCS ${SURFACE_SRCS} ${IMPLICIT_SRCS})

FOLDERFILES("Implicit/Algebraic" "Algebraic")
FOLDERFILES("Implicit/Geometric" "Geometric")
FOLDERFILES("Implicit/Blob" "Blob")
# ADF is broke
#FOLDERFILES("Implicit/ADF" "ADF")
FOLDERFILES("Implicit/MLS" "MLS")
FOLDERFILES("Implicit/Operator" "Operator")
FOLDERFILES("Implicit/Variational" "Variational")

# FOLDERFILES("Implicit/DistanceField" "DistanceField")
FILE(GLOB DF_SRCS 
	Implicit/DistanceField/*.h	
	Implicit/DistanceField/adf.cpp
	Implicit/DistanceField/DistanceField.cpp
	Implicit/DistanceField/sparsescalarlattice.cxx
	Implicit/DistanceField/cachingsparsescalarlattice.cxx	
	Implicit/DistanceField/trianglemesh2.cxx
	Implicit/DistanceField/scalarlattice.cxx
	Implicit/DistanceField/triangle.cxx
	Implicit/DistanceField/boundingvolume.cxx
	Implicit/DistanceField/signeddistancetransform.cxx
	Implicit/DistanceField/trianglenearestpoint.cxx
	Implicit/DistanceField/opengl_utils.cxx)
SOURCE_GROUP("  DistanceField" FILES ${DF_SRCS})
SET(SURFACE_SRCS ${SURFACE_SRCS} ${DF_SRCS})



INCLUDE_DIRECTORIES(
	${BOOST_INCLUDE_PATH} 
	${GPS_SOURCE_DIR} 
	${GPS_SOURCE_DIR}/Surface
	${GPS_SOURCE_DIR}/Surface/Implicit/DistanceField
	${GPS_SOURCE_DIR}/Particles
	${GPS_SOURCE_DIR}/Particles/Attributes
	)


ADD_LIBRARY(Surface ${SURFACE_SRCS})
 
