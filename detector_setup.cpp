#include "detector_setup.h"
#include "algebra_types.h"
#include "algebra_utils.h"
#include "geometry_types.h"

#include <iostream>

DetectorSetup::DetectorSetup(const std::vector<ExtrudedObject>& detectorlist)
    : m_detectors(detectorlist)
{
//    std::cout<<"DetectorSetup::DetectorSetup(const std::vector<ExtrudedObject>&) 1"<<std::endl;
    if (!m_detectors.empty())
        m_ref_detector = m_detectors.begin();
//    std::cout<<"DetectorSetup::DetectorSetup(const std::vector<ExtrudedObject>&) 2"<<std::endl;
}

DetectorSetup::DetectorSetup(const DetectorSetup& other)
    : m_detectors(other.m_detectors)
    //, m_ref_detector(m_detectors.end())
    , m_name(other.m_name)
{
//    std::cout<<"DetectorSetup::DetectorSetup(const std::vector<DetectorSetup>&) 1"<<std::endl;
    //!TODO: also copy the iterator pointing to the reference detector
    if (!m_detectors.empty())
        m_ref_detector = m_detectors.begin();
//    std::cout<<"DetectorSetup::DetectorSetup(const std::vector<DetectorSetup>&) 2"<<std::endl;
}

void DetectorSetup::add_detector(const ExtrudedObject& det)
{
//    std::cout<<"DetectorSetup::add_detector(ExtrudedObject) 1"<<std::endl;
    m_detectors.push_back( ExtrudedObject{ det } );
    //m_detectors.emplace_back( ExtrudedObject { det } );
    if (m_detectors.size() == 1) m_ref_detector = m_detectors.begin();
//    std::cout<<"DetectorSetup::add_detector(ExtrudedObject) 2"<<std::endl;
}

void DetectorSetup::rotate(const Vector& rot_axis, double rot_angle)
{
    if (inEpsilon(rot_angle))
        return;
    for (auto& detector : m_detectors) {
        //Point pos { detector.position() };
        //pos = ::rotate(pos, rot_axis, rot_angle);
        //detector.set_position(pos);
        detector.add_rotation(rot_axis, rot_angle);
    }
}
