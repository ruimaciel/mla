INCLUDE_DIRECTORIES(
	${PROJECT_SOURCE_DIR}
)

set(MLA_SOURCE_FILES
	shared.c++
)

set(MLA_HEADER_FILES
	matrix/all.h++
	matrix/DenseRowMajor.h++
	matrix/SparseCRS.h++
	matrix/SparseDOK.h++
	matrix/SparseCOO.h++
	matrix/SparseCCS.h++
	matrix/static/StaticDenseRowMajor.h++
	matrix/Diagonal.h++
	matrix/traits.h++
	matrix/convert.h++
	vector/all.h++
	vector/SparseCS.h++
	vector/Dense.h++
	vector/traits.h++
	vector/convert.h++
	output.h++
	operations/level1/axpy.h++
	operations/level1/scale.h++
	operations/level1/dot.h++
	operations/level1/asum.h++
	operations/level2/gemv.h++
	operations/level2/syr.h++
	operations/level3/gemm.h++
	operations/level3/syrk.h++
	operations/cuthill_mckee.h++
	ProgressIndicatorStrategy.h++
	LAException.h++
	solvers/substitution.h++
	solvers/CG.h++
	solvers/SolverReturnCodes.h++
	solvers/umfpack.h++
	solvers/Cholesky.h++
	parsers/MatrixMarket.h++
)


add_library( 
	mla 
	SHARED 
	${MLA_SOURCE_FILES}
)

target_link_libraries(
	mla
	${UMFPACK_LIBRARY}
)

set_target_properties(mla
	PROPERTIES PUBLIC_HEADER "${MLA_HEADER_FILES}"
)

install(
	TARGETS mla 
	LIBRARY DESTINATION lib
	PUBLIC_HEADER DESTINATION include/libmla
)

