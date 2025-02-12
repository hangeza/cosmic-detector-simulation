#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <functional>
#include <future>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <thread>
#include <valarray>
#include <vector>

#include "algebra_types.h"
#include "algebra_utils.h"
#include "detector_setup.h"
#include "geometry_types.h"
#include "histogram.h"
#include "sampled_distribution.h"
#include "simulation_algorithms.h"
#include "utilities.h"

constexpr int g_verbosity { 3 };
constexpr double muon_flux_literature_value { 70. };

auto main() -> int
{
    std::cout << "detector MC simulator\n";

    // set up random number generator
    std::random_device rd; // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()

    // the following parameters should be adopted to the individual simulation
    constexpr std::size_t nr_events { 100'000 }; //<! the total number of tracks to be simulated
    constexpr double theta_max { toRad(90.) }; //<! the maximum theta angle taken into account
    constexpr double theta_step { toRad(1.) }; //<! the desired granularity of the simulated angular distributions
    // define rotation axis, here x-axis
    const Vector detector_rotation_axis { R3::Base::X };
    // define the rotation angle
    // set to 0, if setup should not be rotated
    constexpr double detector_rotation_angle { toRad(0.) };

    constexpr std::size_t nr_bins { static_cast<int>(theta_max / theta_step) + 1 };
    std::cout << "nr of bins: " << nr_bins << "\n";

    // define the coincidence level, i.e. the number of detectors in a setup which have to provide a signal for one event
    // -1 for auto, i.e. coinc level is set to the number of detectors
    constexpr int min_coincidence_count { -1 };

    // vector of 2d polygon vertices defining the shape (in x-y-plane) of the detector.
    // note, that the points have to be in geometrical sequential order in
    // counter-clockwise orientation (i.e. a sequence of points defining the detector outline
    // going in ccw direction). The first point is considered to close the polygon together
    // with the last element in the vector

    // definition of the large double-paddle detector in IIPI-JLU lab
    const std::vector<Point> large_paddle_points_upper {
        { -150., -87.5 },
        { 150., -87.5 },
        { 150., 87.5 },
        { -150., 87.5 }
    };
    const std::vector<Point> large_paddle_points_lower {
        { -150., -100. },
        { 150., -100. },
        { 150., 100. },
        { -150., 100. }
    };

    // definition of the MuonPi standard-size (octagon) detector
    const std::vector<Point> octagon_points {
        { -126.5, -20. },
        { -91.5, -62.5 },
        { 91.5, -62.5 },
        { 126.5, -20. },
        { 126.5, 20. },
        { 91.5, 62.5 },
        { -91.5, 62.5 },
        { -126.5, 20. }
    };

    // definition of the MuonPi half-size detector
    const std::vector<Point> half_size_detector_points {
        { -130., -20. },
        { -60., -62.5 },
        { 0., -62.5 },
        { 0., 62.5 },
        { -60., 62.5 },
        { -130., 20. }
    };

    // definition of the MuonPi hexagon (small-size) detector
    constexpr double hex_length_a { 34.64 };
    constexpr double hex_length_b { 30.0 };
    const std::vector<Point> hexagon_detector_points {
        { -hex_length_a, 0. },
        { -hex_length_a / 2, -hex_length_b },
        { hex_length_a / 2., -hex_length_b },
        { hex_length_a, 0. },
        { hex_length_a / 2, hex_length_b },
        { -hex_length_a / 2, hex_length_b }
    };

    // definition of the large detector bars for JLU cosmic detector array
    const std::vector<Point> large_bar_points {
        { -500., -50. },
        { 500., -50. },
        { 500., 50. },
        { -500., 50. }
    };

    // definition of the JLU cubical plastic neutron detector
    const std::vector<Point> neutron_cube_points {
        { -50., -50. },
        { 50., -50. },
        { 50., 50. },
        { -50., 50. }
    };

    // create 3d objects of type ExtrudedObject defined by the 2d outline,
    // a global position offset and a thickness
    ExtrudedObject detector1 { large_paddle_points_lower, { 0., 0., -100. }, 7. };
    ExtrudedObject detector2 { large_paddle_points_upper, { 0., 0., 100. }, 7. };

    ExtrudedObject neutron_detector { neutron_cube_points, { 0., 0., -50. }, 100. };

    // create 3d objects of type ExtrudedObject but using the constructor for generation of a
    // circular shape specified by a global position offset, radius, thickness and an optional
    // number of vertex points to generate the circle
    ExtrudedObject round_detector1 { { 0., 0., 0. }, 50., 10. };
    ExtrudedObject round_detector2 { { 0., 0., 100. }, 50., 10. };

    ExtrudedObject fiber { { 0., 0., -250. }, 0.5, 500., 16 };
    fiber.add_rotation(R3::Base::Y, toRad(90.));

    // Now create the detector setup consisting of one or more ExtrudedObject instances as detectors
    // and a reference volume, which countains all detectors.
    // In case the ref volume is not specified, the DetectorSetup instance will automatically infer one
    // from the bounding box containing all detectors (+100% margin in each dimension)
    DetectorSetup setup { { detector1, detector2 } /*, reference_box */ };

    // alternative example:
    // construct a detector setup with several detectors which are individually aligned
    // the setup describes a set of scintillating fibers which are rotated and shifted to
    // align along x and y axes
    /*
    constexpr struct fibertracker_param_struct {
        std::size_t nr_fibers_per_layer { 16 };
        double fiber_diameter { 1.05 };
        double fiber_radius { fiber_diameter/2 };
        double fiber_length { 100. };
        double sublayer_distance { 2. };
        double layer_distance { 50. };
        
    } fiberparams;
    DetectorSetup setup { {} };
    //setup.add_detector(trigger_detector);
    for (std::size_t i = 0; i < fiberparams.nr_fibers_per_layer; ++i) {
        ExtrudedObject det_l1x { { 0., -8. + static_cast<double>(i) * 1. + 0.1, -fiberparams.fiber_length/2 }, fiberparams.fiber_radius, fiberparams.fiber_length, 16 };
        ExtrudedObject det_l1y { { -fiberparams.sublayer_distance, -8.5 + static_cast<double>(i) * 1. + 0.1, -fiberparams.fiber_length/2 }, fiberparams.fiber_radius, fiberparams.fiber_length, 16 };
        ExtrudedObject det_l2x { { -fiberparams.layer_distance, -8. + static_cast<double>(i) * 1. + 0.1, -fiberparams.fiber_length/2 }, fiberparams.fiber_radius, fiberparams.fiber_length, 16 };
        ExtrudedObject det_l2y { { -(fiberparams.layer_distance + fiberparams.sublayer_distance), -8.5 + static_cast<double>(i) * 1. + 0.1, -fiberparams.fiber_length/2 }, fiberparams.fiber_radius, fiberparams.fiber_length, 16 };
        det_l1x.add_rotation(R3::Base::Y, toRad(90.));
        det_l1y.add_rotation(R3::Base::Y, toRad(90.));
        det_l1y.add_rotation(R3::Base::Z, toRad(90.));
        det_l2x.add_rotation(R3::Base::Y, toRad(90.));
        det_l2y.add_rotation(R3::Base::Y, toRad(90.));
        det_l2y.add_rotation(R3::Base::Z, toRad(90.));
        setup.add_detector(det_l1x);
        setup.add_detector(det_l1y);
        setup.add_detector(det_l2x);
        setup.add_detector(det_l2y);
    }
    setup.autogenerate_ref_volume();
    std::function<bool(const std::valarray<bool>&)> trigger_fn = [&fiberparams](const std::valarray<bool>& hitvector) {
        bool first_xlayer_or { or_trigger(hitvector[std::slice(0,fiberparams.nr_fibers_per_layer,4)]) };
        bool first_ylayer_or { or_trigger(hitvector[std::slice(1,fiberparams.nr_fibers_per_layer,4)]) };
        bool second_xlayer_or { or_trigger(hitvector[std::slice(2,fiberparams.nr_fibers_per_layer,4)]) };
        bool second_ylayer_or { or_trigger(hitvector[std::slice(3,fiberparams.nr_fibers_per_layer,4)]) };
        //return true;
        bool all_layer_coinc { first_xlayer_or && first_ylayer_or && second_xlayer_or && second_ylayer_or }; 
        //if (all_layer_coinc) std::cout<<"trigger lambda: t1="<<first_xlayer_or<<std::endl;
        return all_layer_coinc;
    };
    setup.set_trigger_function(trigger_fn);
    
*/

    // add a rotation to the entire system
    // the pivot point is the origin in the detector coordinate system
    // all detectors are rotated at the pivot about the given axis and a given angle
    setup.rotate(detector_rotation_axis, detector_rotation_angle);

    for (const auto& detector : setup.detectors()) {
        auto bounds { detector.bounding_box() };
        Point dimensions { bounds.second - bounds.first };
        std::cout << "** Detector **" << std::endl;
        std::cout << "detector bounds: min=" << bounds.first << " max=" << bounds.second << "\n";
        std::cout << "detector dimensions=" << dimensions << "\n";
        std::cout << "detector base area=" << detector.getBaseArea() << "\n";
        std::cout << "detector volume=" << detector.getVolume() << "\n";
    }

    // simulate the effective area (geometric aperture) at theta=0 of the detector system
    // this quantity may be used later to infer the expected detector count rate
    [[maybe_unused]] const double effective_area_sqm { simulate_geometric_aperture(setup, gen, nr_events * 10, 0, min_coincidence_count) };

    // uncomment the following block to calculate the double differential acceptance
    // as function of phi and theta
    /*
    [[maybe_unused]] const auto acceptance_phi_theta = theta_phi_scan<361, 46>(setup, gen, nr_events, 0., theta_max, -pi(), pi());
*/

    // initialize the histogram vector
    std::vector<Histogram> histos {};

    // run a scan over theta angle (uniformly distributed)
    // to record the detector acceptance, if required
    /*
    theta_scan(setup, gen, nr_events, 0., theta_max, nr_bins, &histos);
*/

    if (min_coincidence_count < 0)
        setup.set_trigger_function(and_trigger);
    else
        setup.set_trigger_multiplicity(min_coincidence_count);

    // run a full simulation and append the resulting histograms
    // to the already existing histogram vector
    // this will simulate <nr-events> MC-generated tracks and return the total acceptance
    // i.e. the ratio of detected tracks to total generated tracks
    DataItem<double> detector_acceptance { cosmic_simulation(setup, gen, nr_events, &histos, nr_bins, theta_max, min_coincidence_count) };
    if (detector_acceptance.value != 0.) {
        DataItem<double> countrate_item { detector_acceptance };
        // calculate the count rate conversion based on given reference flux value
        // scaled by the effective detector area
        double countrate_conversion {
            2. * pi() * effective_area_sqm / 3. * muon_flux_literature_value
        };
        countrate_conversion /= 1e6 * effective_area_sqm / setup.ref_volume().getBaseArea();
        countrate_item.value *= countrate_conversion;
        countrate_item.error *= countrate_conversion;
        std::cout << "** Detector acceptance and expected count rate **" << std::endl;
        std::cout << " acceptance = " << detector_acceptance.value << " +- " << detector_acceptance.error << "\n";
        std::cout << " count rate = " << countrate_item.value << " +- " << countrate_item.error << "\n";
    }

    // export each histogram into a separate file (human readable ASCII format)
    for (auto histo : histos) {
        histo.export_file(histo.getName() + ".hist");
    }

    // quit here in case an angular detector sweep is not required
    exit(0);

    // alternatively: run a sweep over angular range of detector orientation
    // return a list of acceptance vs angle including statistical errors
    auto acceptance_dataseries { cosmic_simulation_detector_sweep(setup, gen, nr_events, detector_rotation_axis, toRad(-90.), toRad(90.), 181, min_coincidence_count) };

    // define a second list which will hold count rate values calculated from acceptance
    MeasurementVector<double, double> countrate_vs_angle_dataseries {};

    // calculate count rate for every angular acceptance
    for (const auto& [angle_value, acceptance_value] : acceptance_dataseries) {
        DataItem<double> countrate_item { acceptance_value };
        // calculate the count rate conversion based on given reference flux value
        // scaled by the effective detector area
        double countrate_conversion {
            2. * pi() * effective_area_sqm / 3. * muon_flux_literature_value
        };
        countrate_conversion /= 1e6 * effective_area_sqm / setup.ref_volume().getBaseArea();
        countrate_item.value *= countrate_conversion;
        countrate_item.error *= countrate_conversion;
        countrate_vs_angle_dataseries.emplace_back(angle_value, countrate_item);
    }

    // export data series for acceptance and count rate vs angle
    export_file(acceptance_dataseries, "detector_sweep_acceptances.dat");
    export_file(countrate_vs_angle_dataseries, "detector_sweep_countrate.dat");

    exit(0);
}
