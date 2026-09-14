#ifndef MPMC_PT_ANDROID_PORTABILITY_H
#define MPMC_PT_ANDROID_PORTABILITY_H

#include <stdint.h>

#if defined(__GNUC__) || defined(__clang__)
#define MPMC_PT_ANDROID_EXPORT __attribute__((visibility("default")))
#else
#define MPMC_PT_ANDROID_EXPORT
#endif

#ifdef __cplusplus
extern "C" {
#endif

// Version 1 is intentionally tiny: it proves that the model-neutral PT core can
// be built as an Android shared library without introducing JNI, gRPC, process,
// UI, or parameter-database policy into this boundary.
MPMC_PT_ANDROID_EXPORT uint32_t mpmc_pt_android_core_abi_version(void);

MPMC_PT_ANDROID_EXPORT const char* mpmc_pt_android_core_convention(void);

#ifdef __cplusplus
}
#endif

#endif // MPMC_PT_ANDROID_PORTABILITY_H
