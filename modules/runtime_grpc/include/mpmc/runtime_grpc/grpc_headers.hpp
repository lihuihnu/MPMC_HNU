#ifndef MPMC_RUNTIME_GRPC_GRPC_HEADERS_HPP
#define MPMC_RUNTIME_GRPC_GRPC_HEADERS_HPP

// gRPC's public compatibility API still declares a deprecated credential type
// in an inline std::vector constructor. MSVC diagnoses that declaration while
// parsing the third-party header even when application code never uses it.
// Scope C4996 suppression to gRPC headers only; project code remains /W4 /WX.
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable : 4996)
#endif

#include <grpc/grpc_security_constants.h>
#include <grpcpp/create_channel.h>
#include <grpcpp/grpcpp.h>
#include <grpcpp/health_check_service_interface.h>
#include <grpcpp/resource_quota.h>
#include <grpcpp/security/credentials.h>
#include <grpcpp/security/server_credentials.h>
#include <grpcpp/server.h>
#include <grpcpp/server_builder.h>
#include <grpcpp/support/channel_arguments.h>
#include <grpcpp/support/status.h>

#if defined(_MSC_VER)
#pragma warning(pop)
#endif

#endif // MPMC_RUNTIME_GRPC_GRPC_HEADERS_HPP
