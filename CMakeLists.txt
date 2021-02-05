cmake_minimum_required(VERSION 3.10)

project(SphereMovement LANGUAGES C CXX)

add_subdirectory(src sphereLibrary)

add_executable(app app.c)

target_link_libraries(app sphereLibrary m)


