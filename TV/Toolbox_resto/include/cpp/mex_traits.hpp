#ifndef MEX_TRAITS_HPP
#define MEX_TRAITS_HPP

#include <mex.h>
#include <matrix.h>


template <typename T>
struct mxClassIDTrait {
	static const mxClassID type = mxUNKNOWN_CLASS;
};
template <>
struct mxClassIDTrait<mxLogical> {
	static const mxClassID type = mxLOGICAL_CLASS;
};
template <>
struct mxClassIDTrait<uint8_t> {
	static const mxClassID type = mxUINT8_CLASS;
};
template <>
struct mxClassIDTrait<int8_t> {
	static const mxClassID type = mxINT8_CLASS;
};
template <>
struct mxClassIDTrait<uint16_t> {
	static const mxClassID type = mxUINT16_CLASS;
};
template <>
struct mxClassIDTrait<int16_t> {
	static const mxClassID type = mxINT16_CLASS;
};
template <>
struct mxClassIDTrait<uint32_t> {
	static const mxClassID type = mxUINT32_CLASS;
};
template <>
struct mxClassIDTrait<int32_t> {
	static const mxClassID type = mxINT32_CLASS;
};
template <>
struct mxClassIDTrait<uint64_t> {
	static const mxClassID type = mxUINT64_CLASS;
};
template <>
struct mxClassIDTrait<int64_t> {
	static const mxClassID type = mxINT64_CLASS;
};
template <>
struct mxClassIDTrait<float> {
	static const mxClassID type = mxSINGLE_CLASS;
};
template <>
struct mxClassIDTrait<double> {
	static const mxClassID type = mxDOUBLE_CLASS;
};

#endif