#include <mpmc/model_configuration/pr76_model_registry.hpp>

#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#include <bcrypt.h>
#elif defined(__linux__)
#include <sys/random.h>
#elif defined(__APPLE__)
#include <TargetConditionals.h>
#if TARGET_OS_OSX
#include <Security/SecBase.h>
#include <Security/SecRandom.h>
#endif
#endif

namespace mpmc::model_configuration {
ModelHandleEntropy system_model_handle_entropy() {
#if defined(_WIN32)
    ModelHandleEntropy bytes{};
    // OS system-preferred CSPRNG; no provider handle or process-global PRNG.
    if (BCryptGenRandom(nullptr, bytes.data(), static_cast<ULONG>(bytes.size()),
                        BCRYPT_USE_SYSTEM_PREFERRED_RNG) == 0) {
        return bytes;
    }
#elif defined(__linux__)
    ModelHandleEntropy bytes{};
    // A small, nonblocking request fails closed if the kernel RNG is not ready,
    // interrupted or incomplete. Never downgrade to a predictable fallback.
    const auto received = ::getrandom(bytes.data(), bytes.size(), GRND_NONBLOCK);
    if (received >= 0 && static_cast<std::size_t>(received) == bytes.size()) {
        return bytes;
    }
#elif defined(__APPLE__) && TARGET_OS_OSX
    ModelHandleEntropy bytes{};
    // macOS system CSPRNG; return only a successful complete request.
    if (SecRandomCopyBytes(kSecRandomDefault, bytes.size(), bytes.data()) == errSecSuccess) {
        return bytes;
    }
#endif
    throw ModelRegistryError(ModelRegistryErrorCode::entropy_unavailable,
                             "model registry: secure handle entropy unavailable");
}
} // namespace mpmc::model_configuration
