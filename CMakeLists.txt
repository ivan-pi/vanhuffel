cmake_minimum_required(VERSION 3.10)

project(vanhuffel 
    VERSION 0.1.0 
    LANGUAGES Fortran
    DESCRIPTION "Partial Total Least Squares Routines by Sabine Van Huffel")

add_library(vanhuffel
    src/bidiag.f
    src/cancel.f
    src/damin.f
    src/dtls.f
    src/estim.f
    src/housh.f
    src/init.f
    src/nsingv.f
    src/psvd.f
    src/ptls.f
    src/qlstep.f
    src/qrql.f
    src/qrstep.f
    src/restor.f
    src/tr2.f
    )

include(FindBLAS REQUIRED)
target_link_libraries(vanhuffel PRIVATE blas)
target_link_libraries(vanhuffel PRIVATE linpack)

add_subdirectory(examples EXCLUDE_FROM_ALL)
