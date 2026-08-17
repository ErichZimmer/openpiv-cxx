#pragma once

#if defined(EXPORT_DLL_SYMBOLS)

    #if defined(_WIN32) || defined(__CYGWIN__)
        #define DLL_EXPORT __declspec(dllexport)
    #elif defined(__GNUC__) || defined(__clang__)
        #define DLL_EXPORT __attribute__((visibility("default")))
    #else
        #define DLL_EXPORT
    #endif

#else

    #define DLL_EXPORT

#endif