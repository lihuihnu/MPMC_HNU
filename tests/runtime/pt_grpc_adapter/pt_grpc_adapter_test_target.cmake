include_guard(GLOBAL)

function(mpmc_add_pt_grpc_adapter_contract_test target_name)
    if(NOT TARGET mpmc::runtime_grpc)
        message(FATAL_ERROR
            "mpmc_add_pt_grpc_adapter_contract_test requires mpmc::runtime_grpc")
    endif()

    add_executable(${target_name}
        "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/pt_grpc_adapter_test.cpp"
        "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/pt_grpc_adapter_headers.cpp")
    target_link_libraries(${target_name} PRIVATE mpmc::runtime_grpc)
endfunction()
