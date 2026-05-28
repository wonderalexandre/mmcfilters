#pragma once

/**
 * @brief Public include for local-event primitives.
 *
 * @details
 * The local-event layer computes finite-window descriptors over the
 * image-domain proper parts owned by a `MorphologicalTree`. The public surface
 * is the generic event engine. Concrete bitquad and contour computations are
 * owned by their attribute-kernel implementations under
 * `mmcfilters/attributes/computers/`.
 *
 * `EventEngine` is the public policy-based extension point for new finite
 * window computations.
 */

#include "EventEngine.hpp"
