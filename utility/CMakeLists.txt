add_library(utility
    jansson_wrapper.cpp
)

target_link_libraries(utility
    ${JANSSON_LIBRARIES}
)
