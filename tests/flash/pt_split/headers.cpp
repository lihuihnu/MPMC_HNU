#include <mpmc/flash/pr76_phase_set.hpp>
#include <mpmc/flash/pr76_phase_set.hpp>
#include <mpmc/flash/pr76_split.hpp>
#include <mpmc/flash/pt_phase_set.hpp>
#include <mpmc/flash/pt_split.hpp>
#include <mpmc/flash/pt_vle_phase_set.hpp>
#include <mpmc/flash/rachford_rice.hpp>

bool pt_split_headers() {
    const auto r=mpmc::flash::solve_rachford_rice(
        std::vector<double>{.6,.4},std::vector<double>{std::log(2.0),std::log(.25)});
    return r.status==mpmc::flash::RachfordRiceStatus::interior;
}
