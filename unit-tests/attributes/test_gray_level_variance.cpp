#include "support/TestSupport.hpp"

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

namespace {

template <class Real, AltitudeValue T> void checkVariance(const ValuedMorphologicalTree<T>& valued) {
    const auto [names, values] = AttributeComputation::computeSingleAttribute<Real>(valued, GrayLevelVariance);
    for (NodeId node : valued.topology().aliveNodeIds()) {
        std::vector<long double> samples;
        for (PixelId pixel : valued.topology().nodeSupport(node))
            samples.push_back(valued.nodeAltitude(valued.topology().smallestNode(pixel)));
        long double reference = 0;
        if (!samples.empty()) {
            const long double origin = samples.front();
            long double centeredMean = 0;
            for (long double value : samples) centeredMean += value-origin;
            centeredMean /= samples.size();
            for (long double value : samples) {
                const long double deviation = (value-origin)-centeredMean;
                reference += deviation*deviation;
            }
            reference /= samples.size();
        }
        const Real actual = values[names.linearIndex(node, GrayLevelVariance)];
        require(std::isfinite(actual), "finite variance on the reference fixtures");
        const long double tolerance = 128*std::numeric_limits<Real>::epsilon()*std::max(reference,1.0L);
        require(std::abs(static_cast<long double>(actual)-reference)<=tolerance,
                "variance agrees with independent centered two-pass support enumeration");
    }
}

template <AltitudeValue T> void checkImage(std::vector<T> pixels) {
    auto image = Image<T>::fromExternal(pixels.data(),1,static_cast<int>(pixels.size()));
    for (bool isMax : {true,false}) {
        auto valued = isMax ? MorphologicalTreeFactory::createMaxTree(image,1.0)
                            : MorphologicalTreeFactory::createMinTree(image,1.0);
        checkVariance<double>(valued);
        checkVariance<float>(valued);
    }
}

} // namespace

int main() {
    checkImage<std::int32_t>({1<<30,(1<<30)+1});
    checkImage<std::int32_t>({-(1<<30),-(1<<30)+1});
    checkImage<std::int32_t>({std::numeric_limits<std::int32_t>::min(),std::numeric_limits<std::int32_t>::max()});
    checkImage<double>({0.125,0.375,0.375,-0.125,1.125,-0.5});
    checkImage<double>({1e12,1e12+0.5});
    checkImage<double>({-7,-7,3,3,9,-2,1,9});
    checkImage<float>({-1000,-1000,-1000});
    checkImage<std::uint16_t>({0,65535,0,32768,12});

    const std::vector<NodeId> parent{0,0,0,1,1};
    const std::vector<NodeId> owner{3,4,2};
    const std::vector<double> altitude{100,-100,1000,1e9,1e9+1};
    auto valued = makeValuedHierarchy<double>(parent,owner,altitude,0,MorphologicalTreeSemantics{});
    checkVariance<double>(valued);
    checkVariance<float>(valued);
    valued.mergeNodeIntoParent(3);
    checkVariance<double>(valued);
    checkVariance<float>(valued);
}
