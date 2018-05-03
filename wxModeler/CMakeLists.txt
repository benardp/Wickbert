ADD_EXECUTABLE(wxModeler WIN32 MACOSX_BUNDLE
    wxModeler.cpp
    wxModeler.h
    DisplayWindow.cpp
    DisplayWindow.h
    ParticleTree.cpp
    ParticleTree.h
    SurfaceTree.cpp
    SurfaceTree.h
    LogFrame.cpp
    LogFrame.h
    ClipartFrame.cpp
    ClipartFrame.h
    SVGshader.cpp
    SVGshader.h
    )

IF(WB_USE_CGAL)
add_to_cached_list( CGAL_EXECUTABLE_TARGETS wxModeler)
ENDIF(WB_USE_CGAL)

INCLUDE_DIRECTORIES(
    ${GPS_SOURCE_DIR}
    ${GPS_SOURCE_DIR}/Surface
    ${GPS_SOURCE_DIR}/Particles
    ${GPS_SOURCE_DIR}/Particles/Attributes
    )

TARGET_LINK_LIBRARIES(wxModeler
    Surface
    Particles
    gm
    ism
    ${OPENMESH_LIBRARIES}
    ${wxWidgets_LIBRARIES}
    ${GSL_LIBRARIES}
    ${OPENGL_LIBRARIES}
    )

IF(APPLE)
    set_target_properties(wxModeler PROPERTIES MACOSX_BUNDLE_INFO_PLIST ${CMAKE_CURRENT_SOURCE_DIR}/Info.plist)

    SET(wxModeler_PATH ${EXECUTABLE_OUTPUT_PATH})
    IF(NOT wxModeler_PATH)
    SET(wxModeler_PATH ${CMAKE_CURRENT_BINARY_DIR})
    ENDIF(NOT wxModeler_PATH)

    FIND_PROGRAM(VTK_APPLE_RESOURCE Rez /Developer/Tools)
    IF(VTK_APPLE_RESOURCE)
        ADD_CUSTOM_COMMAND(TARGET wxModeler POST_BUILD
            COMMAND ${VTK_APPLE_RESOURCE}
            ARGS Carbon.r -o
            ${wxModeler_PATH}/wxModeler
            TARGET wxModeler)
    ENDIF(VTK_APPLE_RESOURCE)

ENDIF(APPLE)

IF(WIN32)
    IF(WB_USE_CG)
    TARGET_LINK_LIBRARIES(wxModeler ${CG_LIBRARY} ${CG_GL_LIBRARY})
    ENDIF(WB_USE_CG)
ENDIF(WIN32)

IF(WB_USE_VTK)  
TARGET_LINK_LIBRARIES(wxModeler ${VTK_LIBRARIES})
ENDIF(WB_USE_VTK)

IF(WB_USE_SFL)
TARGET_LINK_LIBRARIES(wxModeler sfl )
ENDIF(WB_USE_SFL)

IF(WB_USE_ITK)
TARGET_LINK_LIBRARIES(wxModeler ${ITK_LIBRARIES} )
ENDIF(WB_USE_ITK)

IF(WB_USE_CGAL)
TARGET_LINK_LIBRARIES(wxModeler ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
ENDIF(WB_USE_CGAL)
