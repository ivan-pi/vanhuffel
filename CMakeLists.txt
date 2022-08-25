cmake_minimum_required(VERSION 3.12)

project(vanhuffel 
    VERSION 0.1.0 
    LANGUAGES Fortran
    DESCRIPTION "Partial Total Least Squares Routines by Sabine Van Huffel")

add_library(support OBJECT     
    src/cancel.f
    src/damin.f
    src/housh.f
    src/init.f
    src/qlstep.f
    src/qrstep.f
    src/restor.f
    src/tr2.f)

add_library(vanhuffel
    src/bidiag.f
    src/dtls.f
    src/estim.f
    src/nsingv.f
    src/psvd.f
    src/ptls.f
    src/qrql.f)

add_library(vanhuffel::ptls ALIAS vanhuffel)

target_include_directories(vanhuffel INTERFACE include)
target_link_libraries(vanhuffel PRIVATE support)
target_sources(vanhuffel PRIVATE src/linpack/dsvdc.f src/linpack/dqrdc.f)

include(FindBLAS REQUIRED)
target_link_libraries(vanhuffel PRIVATE blas)
#target_link_libraries(vanhuffel PRIVATE linpack)

add_subdirectory(examples EXCLUDE_FROM_ALL)

install(TARGETS vanhuffel ARCHIVE)