#include "rand.h"
#include <chrono>

rand::Config rand::config_{};
thread_local rand::Engine rand::engine_{};

void rand::set_master_seed(uint32_t seed) {
    std::lock_guard<std::mutex> lock(config_.mutex);
    rand::config_.use_manual_seed = true;
    rand::config_.master_seed = seed;
    rand::config_.thread_counter = 0;
}

void rand::use_hardware_entropy() {
    std::lock_guard<std::mutex> lock(config_.mutex);
    rand::config_.use_manual_seed = false;
    std::random_device rd;
    rand::config_.master_seed = rd();
    rand::config_.thread_counter = 0;
}

void rand::init_thread() {
   // The extra check !engine_.initialized would avoid erroneous re-initialization,
   // which would mess up the thread_counter and maybe compromise reproducibility.
    if (!engine_.initialized) {
        std::lock_guard<std::mutex> lock(rand::config_.mutex);
        // Set the actual engine seed using master seed and thread counter
        std::seed_seq seq{rand::config_.master_seed, rand::config_.thread_counter++};
        rand::engine_.engine.seed(seq);
        rand::engine_.initialized = true;
    }
}

uint32_t rand::get_master_seed() {
   return rand::config_.master_seed;
}