# NiHuLinkerSettings.cmake
# Linker flags configuration

# Variables set by this script
#   CMAKE_EXE_LINKER_FLAGS
#   CMAKE_SHARED_LINKER_FLAGS

if(WIN32)
    if(MINGW)
        set(CMAKE_EXE_LINKER_FLAGS "-MD ${CMAKE_EXE_LINKER_FLAGS} -s")
        set(CMAKE_SHARED_LINKER_FLAGS "-MD ${CMAKE_SHARED_LINKER_FLAGS} -s")
    endif()
endif()
