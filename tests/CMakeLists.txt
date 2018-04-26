ADD_EXECUTABLE(test_interval 
  test_interval.cpp
)

INCLUDE_DIRECTORIES(
	${GPS_SOURCE_DIR} 
	${GPS_SOURCE_DIR}/Surface
	${GPS_SOURCE_DIR}/Particles
	${GPS_SOURCE_DIR}/Particles/Attributes
	)

IF(APPLE)
    TARGET_LINK_LIBRARIES(test_interval 
	
    Surface 
	gm 
	ism
    )
ENDIF(APPLE)

#for UNIX, this should be changed according to the specific system

IF(UNIX AND NOT APPLE)
   
	TARGET_LINK_LIBRARIES(test_interval
	
	Surface 
	
	gm 
	
	ism 
     	)

ENDIF(UNIX AND NOT APPLE)


IF(WIN32)
TARGET_LINK_LIBRARIES(test_interval
    gm
    ism
    Surface
)
ENDIF(WIN32)
