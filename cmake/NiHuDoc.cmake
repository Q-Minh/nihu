# NiHuDoc.cmake
# Documentation CMake file for NiHu

message(STATUS "Looking for Doxygen ...")

find_package(Doxygen QUIET)

if(DOXYGEN_FOUND)
	# Check documentation path
	if(NOT DEFINED NIHU_HTML_DOC_DIR)
		set(NIHU_HTML_DOC_DIR "${CMAKE_BINARY_DIR}/doc/html")
		file(MAKE_DIRECTORY ${NIHU_HTML_DOC_DIR})
	endif()
	message(STATUS "\tDoxygen html documentation path: ${NIHU_HTML_DOC_DIR}")
	
	if(NOT NIHU_MATHJAX_DISABLE)
		message(STATUS "\tEnabling MathJax in Doxygen documentation")
		set(NIHU_MATHJAX_USE "YES")
		if(DEFINED NIHU_MATHJAX_PATH)
			file(RELATIVE_PATH NIHU_MATHJAX_RELPATH "${NIHU_HTML_DOC_DIR}" "${NIHU_MATHJAX_PATH}")
			message(STATUS "MathJax relative path: ${NIHU_MATHJAX_RELPATH}")
		else()
			set(NIHU_MATHJAX_RELPATH "http://cdn.mathjax.org/mathjax/latest")
			message(STATUS "\tWill use CDN for MathJax: ${NIHU_MATHJAX_RELPATH}")
		endif()
	else()
		message(STATUS "\tDisabling MathJax in Doxygen documentation")
		set(NIHU_MATHJAX_USE "NO")
		set(NIHU_MATHJAX_RELPATH "")
	endif()

	# Configuration of doxyfiles
	configure_file(${CMAKE_SOURCE_DIR}/../Doxyfile.in 
		${CMAKE_BINARY_DIR}/Doxyfile)

	# doc target - runs doxygen and creates doxyerror
	add_custom_target(doc
		${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/Doxyfile 2> ${CMAKE_BINARY_DIR}/DoxyError.full
		WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/..
		COMMENT "Generating API documentation with Doxygen")

	# filter the result of doxyerror
	if(UNIX)
		add_custom_command(
			TARGET doc
			POST_BUILD
			COMMAND grep -v "\"Detected potential recursive class relation\"" ${CMAKE_BINARY_DIR}/DoxyError.full 
			| grep -v "\"Member type (typedef)\"" 
			| grep -v "\"no uniquely matching class member\""
			| grep . >${CMAKE_BINARY_DIR}/DoxyError
			DEPENDS "${CMAKE_MODULE_PATH}/WriteMatlabTestRunner.cmake"
			COMMENT "Filtering Doxygen errors and warnings"
		)
	endif()
endif()