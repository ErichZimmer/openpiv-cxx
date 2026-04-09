
#pragma once

#include "loaders/image_loader.h"
#include "core/dll_export.h"

namespace openpiv::core {

    /// loader of TIFF images with support for bit depths
    /// over 8-bits per channel; will sniff the input data for
    /// TIFF of 0x49 49 2a 00 or 0x4d 4d 00 2a
    class tiff_image_loader : public image_loader
    {
    public:
        tiff_image_loader();
        ~tiff_image_loader();

        DLL_EXPORT std::string name() const override;
        DLL_EXPORT int priority() const override;
        DLL_EXPORT image_loader_ptr_t clone() const override;
        DLL_EXPORT bool can_load( std::istream& ) const override;
        DLL_EXPORT bool can_save() const override;
        DLL_EXPORT size_t num_images() const override;

        DLL_EXPORT bool open( std::istream& is ) override;
        DLL_EXPORT bool extract( size_t index, image_g16& ) override;
        DLL_EXPORT bool extract( size_t index, image_gf32& ) override;
        DLL_EXPORT bool extract( size_t index, image_gf64& ) override;
        DLL_EXPORT bool extract( size_t index, image_rgba16& ) override;

        DLL_EXPORT void save( std::ostream&, const image_g16& ) const override;
        DLL_EXPORT void save( std::ostream&, const image_gf32& ) const override;
        DLL_EXPORT void save( std::ostream&, const image_gf64& ) const override;
        DLL_EXPORT void save( std::ostream&, const image_rgba16& ) const override;
        DLL_EXPORT void save( std::ostream&, const image_view_g16& ) const override;
        DLL_EXPORT void save( std::ostream&, const image_view_gf32& ) const override;
        DLL_EXPORT void save( std::ostream&, const image_view_gf64& ) const override;
        DLL_EXPORT void save( std::ostream&, const image_view_rgba16& ) const override;

    private:
        struct impl;
        std::unique_ptr<impl> impl_;
    };

}