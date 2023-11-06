configure_file(
    "${startle_SOURCE_DIR}/include/startle.hpp.in"
    "${startle_SOURCE_DIR}/include/startle.hpp"
)

add_executable(startle startle.cxx treeio.cxx starhomoplasy.cxx)

# add libraries
target_link_libraries(startle PRIVATE nlohmann_json::nlohmann_json)
target_link_libraries(startle PRIVATE pprint)
target_link_libraries(startle PRIVATE spdlog)
target_link_libraries(startle PRIVATE csv)
target_link_libraries(startle PRIVATE argparse)
