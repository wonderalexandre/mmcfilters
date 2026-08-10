#include "support/TestSupport.hpp"

#include "mmcfilters/attributes/computers/detail/maxdist/Geometry.hpp"
#include "mmcfilters/attributes/computers/detail/maxdist/PQueue.hpp"

using namespace mmcfilters;
using namespace mmcfilters::attributes::computers::detail::maxdist;
using namespace mmcfilters::unit_tests;

namespace {

void requireShiftedBoxRoundTrip() {
    Box2D domain(10, 20, 12, 21);
    requireEqual(domain.width(), 3, "shifted Box2D width");
    requireEqual(domain.height(), 2, "shifted Box2D height");

    for (int idx = 0; idx < domain.width() * domain.height(); ++idx) {
        const Point2D point = domain.point(idx);
        requireEqual(domain.index(point), idx, "shifted Box2D index(point(idx)) round-trip");
        requireEqual(domain.index(point.x(), point.y()), idx, "shifted Box2D index(x, y) round-trip");
    }

    requireEqual(domain.index(9, 20), -1, "shifted Box2D rejects x before top-left");
    requireEqual(domain.index(13, 21), -1, "shifted Box2D rejects x after bottom-right");
    requireEqual(domain.index(10, 19), -1, "shifted Box2D rejects y before top-left");
    requireEqual(domain.index(12, 22), -1, "shifted Box2D rejects y after bottom-right");
}

void requirePQueueLifoClearsSingletonBucket() {
    PQueue queue(16, 3);
    queue.setCost(0, 5);
    queue.insert(0);

    requireEqual(queue.popMinLIFO(), 0, "PQueue singleton LIFO pop");
    require(queue.isEmpty(), "PQueue singleton LIFO pop empties queue");

    queue.setCost(1, 5);
    queue.insert(1);
    requireEqual(queue.minElemFIFO(), 1, "PQueue insert after singleton LIFO must reset bucket head");
    requireEqual(queue.popMinFIFO(), 1, "PQueue FIFO after singleton LIFO reinsertion");
    require(queue.isEmpty(), "PQueue FIFO pop after singleton LIFO reinsertion empties queue");
}

void requirePQueueLifoMaintainsRemainingTail() {
    PQueue queue(16, 4);
    queue.setCost(0, 7);
    queue.setCost(1, 7);
    queue.insert(0);
    queue.insert(1);

    requireEqual(queue.popMinLIFO(), 1, "PQueue multi-element LIFO returns newest bucket element");
    require(!queue.isEmpty(), "PQueue multi-element LIFO leaves older bucket element queued");
    requireEqual(queue.minElemFIFO(), 0, "PQueue multi-element LIFO preserves bucket head");
    requireEqual(queue.popMinFIFO(), 0, "PQueue FIFO removes remaining element after LIFO");
    require(queue.isEmpty(), "PQueue LIFO then FIFO empties queue");

    queue.setCost(2, 7);
    queue.insert(2);
    requireEqual(queue.minElemFIFO(), 2, "PQueue insert after LIFO/FIFO sequence must reset bucket head");
    requireEqual(queue.popMinLIFO(), 2, "PQueue LIFO after LIFO/FIFO reinsertion");
    require(queue.isEmpty(), "PQueue LIFO after LIFO/FIFO reinsertion empties queue");
}

void requirePQueueMixedFifoLifoMaintainsLinks() {
    PQueue queue(16, 4);
    queue.setCost(0, 9);
    queue.setCost(1, 9);
    queue.insert(0);
    queue.insert(1);

    requireEqual(queue.popMinFIFO(), 0, "PQueue FIFO removes oldest bucket element");
    requireEqual(queue.popMinLIFO(), 1, "PQueue LIFO removes remaining element after FIFO");
    require(queue.isEmpty(), "PQueue FIFO then LIFO empties queue");

    queue.setCost(2, 9);
    queue.insert(2);
    requireEqual(queue.minElemFIFO(), 2, "PQueue insert after FIFO/LIFO sequence must reset bucket head");
}

void requirePQueueAcceptsMaxBucketCost() {
    PQueue queue(5, 2);
    queue.setCost(1, 5);
    queue.insert(1);

    requireEqual(queue.minValue(), 5, "PQueue max bucket cost min value");
    requireEqual(queue.maxValue(), 5, "PQueue max bucket cost max value");
    requireEqual(queue.minElemFIFO(), 1, "PQueue max bucket cost head");
    requireEqual(queue.popMinFIFO(), 1, "PQueue max bucket cost pop");
    require(queue.isEmpty(), "PQueue max bucket cost pop empties queue");
}

} // namespace

int main() {
    requireShiftedBoxRoundTrip();
    requirePQueueLifoClearsSingletonBucket();
    requirePQueueLifoMaintainsRemainingTail();
    requirePQueueMixedFifoLifoMaintainsLinks();
    requirePQueueAcceptsMaxBucketCost();
    return 0;
}
