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
	SET(PARTICLES_SRCS ${PARTICLES_SRCS} ${DUD_SRCS})
ENDMACRO(FOLDERFILES DIR NAME)

FILE(GLOB PARTICLES_SRCS
	Particles.cpp Particles.h
	ParticleStuff.cpp ParticleStuff.h
	ParticleStuffParameters.cpp ParticleStuffParameters.h
	ParticleAttribute.cpp ParticleAttribute.h
	ParticleBehavior.cpp ParticleBehavior.h
	ParticleShader.cpp ParticleShader.h
	ParticleSystem.cpp ParticleSystem.h
	ParFileManager.cpp ParFileManager.h
	nv_pbuffer.cpp nv_pbuffer.h
	pstools.cpp pstools.h
	svd.cpp svd.h
)
SOURCE_GROUP(Particles FILES ${PARTICLES_SRCS})

FOLDERFILES("Attributes" "Attributes")
FOLDERFILES("Behaviors" "Behaviors")
FOLDERFILES("Shaders" "Shaders")

INCLUDE_DIRECTORIES(${GPS_SOURCE_DIR} ${BOOST_INCLUDE_PATH} ${GPS_SOURCE_DIR}/Particles ${GPS_SOURCE_DIR}/Particles/Attributes)

ADD_LIBRARY(Particles ${PARTICLES_SRCS})
  
