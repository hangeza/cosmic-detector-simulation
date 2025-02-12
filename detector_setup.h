#pragma once

#include <vector>

#include "geometry_types.h"

/** @brief DetectorSetup - class for managing geometric objects of type ExtrudedObject
 * This class stores a list of ExtrudedObject detector objects for convenience.
 * An arbitrary number of detector objects with individual alignments can be supplied to the
 * constructor. Additionaly, a global rotation of the whole setup can be set.
 * @note A reference volume can be specified through the constructor. If ommited, a ref volume
 * larger than the largest possible boundary box (considering all possible detector rotations)
 * is automatically constructed.
 * @note A trigger function can be explicitely set which handles the way how the whole
 * set of detectors responds with a boolean output (coincidence) when hit by a track.
 * If not explicitely set, the trigger function is by default constructed with a logical "AND"
 * btw. all detectors, i.e. a coincidence signal is only asserted when all detectors
 * are intersected by a track. Alternatively, an AND trigger pattern can be set through the
 * set_trigger_multiplicity method, which takes as an argument the number of detectors
 * which must be hit by a track. The second argument specifies, whether the coincidence is
 * exclusive, i.e. that only the given number of detectors must respond or inclusive, i.e.
 * that at least the given number of detectors must be hit.
*/
class DetectorSetup {
public:
    struct TriggerClass {
        unsigned number { 1u };
        bool invert { false };
        std::function<bool(bool, bool)> logic_function { std::logical_or<bool> {} };
    };
    DetectorSetup() = default;
    DetectorSetup(const DetectorSetup& other);
    DetectorSetup(DetectorSetup&& other);
    DetectorSetup(const std::vector<ExtrudedObject>& detectorlist, const ExtrudedObject& ref_volume = ExtrudedObject::invalid_volume());

    auto detectors() -> std::vector<ExtrudedObject>& { return m_detectors; }
    auto detectors() const -> const std::vector<ExtrudedObject>& { return m_detectors; }
    void add_detector(const ExtrudedObject& det);
    void set_ref_volume(const ExtrudedObject& ref_volume) { m_ref_volume = ref_volume; }
    void autogenerate_ref_volume();
    auto ref_volume() const -> const ExtrudedObject& { return m_ref_volume; }
    auto bounding_box() const -> std::pair<Point, Point>;
    void rotate(const Vector& rot_axis, double rot_angle);
    void reset_rotation();
    auto get_total_volume() const -> double;

    auto intersection(const Line& path) const -> std::vector<LineSegment>;

private:
    auto get_largest_bounding_box() const -> std::pair<Point, Point>;
    std::vector<ExtrudedObject> m_detectors {};
    ExtrudedObject m_ref_volume {};
    std::string m_name {};
};
