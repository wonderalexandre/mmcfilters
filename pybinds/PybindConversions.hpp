#pragma once

#include "../mmcfilters/trees/MorphologicalTree.hpp"
#include "../mmcfilters/utils/RegularGridAdjacency2D.hpp"
#include "../mmcfilters/utils/Image.hpp"
#include "../mmcfilters/utils/Common.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <array>
#include <cmath>
#include <concepts>
#include <limits>
#include <stdexcept>
#include <string>
#include <sstream>
#include <string_view>
#include <vector>

namespace mmcfilters::pybind_utils {

namespace py = pybind11;

/**
 * @brief Functions for converting between C++ buffers/images and NumPy arrays.
 */
enum class FloatingDType {
    Float32,
    Float64,
};

/**
 * @brief Parses floating dtype.
 *
 * @param dtype Requested NumPy floating-point type.
 * @param argumentName Argument name included in validation error messages.
 * @return Parsed floating dtype.
 */
inline FloatingDType parseFloatingDType(py::object dtype, std::string_view argumentName = "dtype") {
    if (dtype.is_none()) {
        return FloatingDType::Float32;
    }

    py::object numpy = py::module_::import("numpy");
    py::object normalized = numpy.attr("dtype")(dtype);
    const std::string name = py::str(normalized.attr("name")).cast<std::string>();
    if (name == "float32") {
        return FloatingDType::Float32;
    }
    if (name == "float64") {
        return FloatingDType::Float64;
    }
    throw std::invalid_argument(std::string(argumentName) + " must be np.float32 or np.float64");
}

/**
 * @brief Parses floating array dtype.
 *
 * @param array NumPy array validated by the operation.
 * @param argumentName Argument name included in validation error messages.
 * @return Parsed floating array dtype.
 */
inline FloatingDType parseFloatingArrayDType(const py::array& array, std::string_view argumentName) {
    py::object numpy = py::module_::import("numpy");
    py::object normalized = numpy.attr("dtype")(array.dtype());
    const std::string name = py::str(normalized.attr("name")).cast<std::string>();
    if (name == "float32") {
        return FloatingDType::Float32;
    }
    if (name == "float64") {
        return FloatingDType::Float64;
    }
    throw std::invalid_argument(std::string(argumentName) + " must be a 1D np.float32 or np.float64 array.");
}

/**
 * @brief Validates adjacency radius.
 *
 * @param radius Neighbourhood radius used by the operation.
 * @param context Operation name used in diagnostics.
 * @return Validated finite adjacency radius.
 */
inline double requireAdjacencyRadius(double radius, std::string_view context) {
    if (!std::isfinite(radius)) {
        throw std::invalid_argument(std::string(context) + " radius must be finite.");
    }
    if (radius < 1.0) {
        throw std::invalid_argument(std::string(context) + " radius must be at least 1.0 so the adjacency stencil contains a non-central neighbour.");
    }
    const double maxSafeRadius = (std::sqrt(static_cast<double>(std::numeric_limits<int>::max())) - 1.0) / 2.0;
    if (radius > maxSafeRadius) {
        throw std::invalid_argument(std::string(context) + " radius exceeds the supported integer stencil range.");
    }
    return radius;
}

/**
 * @brief Creates a two-dimensional regular-grid adjacency.
 *
 * @param rows Number of rows in the domain.
 * @param cols Number of columns in the domain.
 * @param radius Neighbourhood radius used by the operation.
 * @param context Operation name used in diagnostics.
 * @return Validated two-dimensional regular-grid adjacency.
 */
inline RegularGridAdjacency2D makeRegularGridAdjacency2D(int rows, int cols, double radius, std::string_view context) {
    const double validatedRadius = requireAdjacencyRadius(radius, context);
    try {
        return RegularGridAdjacency2D(rows, cols, validatedRadius);
    } catch (const std::invalid_argument& error) {
        throw std::invalid_argument(std::string(context) + ": " + error.what());
    }
}

/**
 * @brief Validates the shape and contiguity of a one-dimensional array.
 *
 * @param bufferInfo Value buffer represented by `bufferInfo`.
 * @param expectedSize Expected number of values.
 * @param argumentName Argument name included in validation error messages.
 */
inline void require1DArray(const py::buffer_info& bufferInfo, py::ssize_t expectedSize, std::string_view argumentName) {
    if (bufferInfo.ndim != 1) {
        std::ostringstream oss;
        oss << argumentName << " must be a 1D array";
        throw std::invalid_argument(oss.str());
    }
    if (bufferInfo.shape[0] != expectedSize) {
        std::ostringstream oss;
        oss << argumentName << " must have length " << expectedSize << ", got " << bufferInfo.shape[0];
        throw std::invalid_argument(oss.str());
    }
}

/**
 * @brief Validates node attribute array.
 *
 * @param attr Attribute requested by the operation.
 * @param tree Tree topology used by the operation.
 * @param argumentName Argument name included in validation error messages.
 * @return Image or array produced by the operation.
 */
template <class Real>
inline py::array_t<Real, py::array::c_style> requireNodeAttributeArray(py::array attr, const MorphologicalTree& tree, std::string_view argumentName = "attr") {
    const py::buffer_info info = attr.request();
    require1DArray(info, tree.getNumInternalNodeSlots(), argumentName);
    if (info.strides[0] != static_cast<py::ssize_t>(sizeof(Real))) {
        throw std::invalid_argument(std::string(argumentName) + " must be C-contiguous.");
    }
    return py::reinterpret_borrow<py::array_t<Real, py::array::c_style>>(attr);
}

/**
 * @brief Validates vector size.
 *
 * @param values Values read or written by the operation.
 * @param expectedSize Expected number of values.
 * @param argumentName Argument name included in validation error messages.
 */
template <class T> inline void requireVectorSize(const std::vector<T>& values, std::size_t expectedSize, std::string_view argumentName) {
    if (values.size() != expectedSize) {
        std::ostringstream oss;
        oss << argumentName << " must have length " << expectedSize << ", got " << values.size();
        throw std::invalid_argument(oss.str());
    }
}

/**
 * @brief Converts to numpy.
 *
 * @param image Image used by the operation.
 * @return Converted to numpy.
 */
template <typename PixelType> inline py::array_t<PixelType> toNumpy(ImagePtr<PixelType> image) {
    int numCols = image->getNumCols();
    int numRows = image->getNumRows();

    std::shared_ptr<PixelType[]> buffer = image->rawDataPtr();
    std::shared_ptr<PixelType[]> bufferCopy = buffer;

    py::capsule free_when_done(new std::shared_ptr<PixelType[]>(bufferCopy), [](void* ptr) { delete reinterpret_cast<std::shared_ptr<PixelType[]>*>(ptr); });

    // 2D shape: (numRows, numCols), row-major strides
    const py::ssize_t itemsize = sizeof(PixelType);
    const std::array<py::ssize_t, 2> shape = {static_cast<py::ssize_t>(numRows), static_cast<py::ssize_t>(numCols)};
    const std::array<py::ssize_t, 2> strides = {static_cast<py::ssize_t>(numCols) * itemsize, itemsize};

    py::array_t<PixelType> numpy(shape, strides, buffer.get(), free_when_done);

    return numpy;
}

/**
 * @brief Transfers a one-dimensional owned buffer to a NumPy array.
 *
 * @param buffer Buffer read or written by the operation.
 * @param n Requested number of elements.
 * @return NumPy array that owns the transferred one-dimensional buffer.
 */
template <typename Real> inline py::array_t<Real> toNumpyOwned(std::vector<Real>&& buffer, int n) {
    auto* owned = new std::vector<Real>(std::move(buffer));
    py::capsule free_when_done(owned, [](void* ptr) { delete reinterpret_cast<std::vector<Real>*>(ptr); });

    return py::array_t<Real>({n}, {static_cast<py::ssize_t>(sizeof(Real))}, owned->data(), free_when_done);
}

/**
 * @brief Transfers an owned buffer to a two-dimensional NumPy array.
 *
 * @param buffer Buffer read or written by the operation.
 * @param rows Number of rows in the domain.
 * @param cols Number of columns in the domain.
 * @return NumPy array that owns the transferred two-dimensional buffer.
 */
template <typename Real> inline py::array_t<Real> toNumpyOwned2D(std::vector<Real>&& buffer, int rows, int cols) {
    auto* owned = new std::vector<Real>(std::move(buffer));
    py::capsule free_when_done(owned, [](void* ptr) { delete reinterpret_cast<std::vector<Real>*>(ptr); });

    return py::array_t<Real>({rows, cols}, {static_cast<py::ssize_t>(sizeof(Real) * cols), static_cast<py::ssize_t>(sizeof(Real))}, owned->data(),
                             free_when_done);
}

/**
 * @brief Copies a NumPy array into shared C++ storage.
 *
 * @param arr NumPy array to validate and convert.
 * @return Shared storage containing a copy of the array values.
 */
template <typename Scalar, typename Array> inline std::shared_ptr<Scalar[]> toSharedPtr(const Array& arr) {
    // Capture the Python object in the deleter so the buffer stays alive.
    return std::shared_ptr<Scalar[]>(static_cast<Scalar*>(arr.request().ptr), [obj = py::object(arr)](Scalar*) mutable {
        // Keep the py::object alive until the shared_ptr is destroyed.
    });
}

} // namespace mmcfilters::pybind_utils
