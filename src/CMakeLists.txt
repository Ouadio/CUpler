add_library(sphereLibrary sphereUtils.c utils.c)

target_include_directories(sphereLibrary PUBLIC "${PROJECT_BINARY_DIR}" ${PROJECT_SOURCE_DIR}/include)

