#pragma once

/**
 * @brief Preferred public include for the incremental attribute API.
 *
 * @details
 * This header gathers the stable attribute-facing contracts that ordinary C++
 * consumers should use:
 * - `Attribute` / `AttributeGroup` and their flat-buffer layouts;
 * - owning result types;
 * - the public incremental computation facade.
 *
 * Implementation helpers live under `mmcfilters/attributes/detail/`. They are
 * installed because several public headers are template/header-only, but they
 * are not part of the stable public API and should not be included directly by
 * ordinary consumers.
 */

#include "AttributeNames.hpp"
#include "AttributeResultTypes.hpp"
#include "AttributeComputation.hpp"
