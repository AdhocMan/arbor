add_executable(network_random EXCLUDE_FROM_ALL network_random.cpp)
add_dependencies(examples network_random)

target_link_libraries(network_random PRIVATE arbor arborio arborenv arbor-sup ${json_library_name})
