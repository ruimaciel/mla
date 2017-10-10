
INCLUDE_DIRECTORIES(
	${PROJECT_SOURCE_DIR}
)

find_package(Boost COMPONENTS system filesystem unit_test_framework REQUIRED)


INCLUDE_DIRECTORIES(
	${PROJECT_SOURCE_DIR}
	${TEST_SOURCE_DIR}/src
)


ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK) 


INCLUDE (CheckIncludeFiles)
SET(FULL_MATRIX_MARKET_FILES_PATH "${CMAKE_CURRENT_SOURCE_DIR}/MatrixMarket/")
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/config_paths.h.in ${CMAKE_CURRENT_SOURCE_DIR}/config_paths.h @ONLY)

# Helper function that simplifies the whole ordeal of adding multiple unit tests
function(MLA_add_unit_test)
	foreach( unit_test_name ${ARGV} )
		# displays message, if required
		if(CMAKE_VERBOSE_MAKEFILE)
			message(STATUS "Adding unit test ${unit_test_name}")
		endif(CMAKE_VERBOSE_MAKEFILE)

		# add the unit test framework
		set(unit_test_source_file "${unit_test_name}.c++")
		add_executable ( ${unit_test_name} ${unit_test_source_file} )
		target_link_libraries(	${unit_test_name}
			mla
			${Boost_FILESYSTEM_LIBRARY}
			${Boost_SYSTEM_LIBRARY}
			${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
		)
		add_test("unit_${unit_test_name}" ${unit_test_name})
	endforeach(unit_test_name)
endfunction(MLA_add_unit_test)


MLA_add_unit_test(
	test_vector
	test_VectorCursor
	test_vector_convert
	test_vector_convert_to_SparseCS
	test_matrix
	test_matrix_DenseRowMajor
	test_matrix_StaticDenseRowMajor
	test_matrix_Diagonal
	test_matrix_SparseDOK
	test_matrix_SparseCOO
	test_matrix_SparseCRS
	test_MatrixCursor
	test_MatrixCursor_DenseRowMajor
	test_MatrixCursor_Diagonal
	test_matrix_convert
	test_blas_level1_axpy
	test_blas_level1_scale
	test_blas_level1_dot
	test_blas_level1_asum
	test_blas_level2_gemv
	test_blas_level2_gemv_DenseRowMajor
	test_blas_level2_gemv_SparseCRS
	test_blas_level2_syr_DenseRowMajor
	test_blas_level2_syr_SparseCRS
	test_blas_level3_gemm
	test_blas_level3_syrk
	test_cuthill_mckee
	test_parser_MatrixMarket
	test_solvers_cg
	test_solvers_umfpack
)

