add_executable(dtls-example dtls-example.f)
target_link_libraries(dtls-example PRIVATE vanhuffel)

add_executable(ptls-example ptls-example.f)
target_link_libraries(ptls-example PRIVATE vanhuffel)

add_executable(psvd-example psvd-example.f)
target_link_libraries(psvd-example PRIVATE vanhuffel)
