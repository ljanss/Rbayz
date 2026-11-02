//
//  BayzR -- rand.h
//  Fast native C++ random number generation, prepared for multi-threading.
//  Note: When re-using this class elsewhere it can be wise to activate again the
//  initialization checks, or else to study well the needed initialization in a main
//  function.
//
//  Created by Luc Janss on 30/10/2026.

#ifndef RAND_H
#define RAND_H

#include <random>
#include <mutex>

class rand {

   // Rand class has a bunch of member variables: seed (generated or user-given),
   // some flags, the random engine, and objects for various distributions.
   // In principle I think these could all be direct member variables, but
   // they are here separated into two struct's for some cleaner organisation
   // (and maybe to allow locking config_ for multi-threading). Internal state is
   // manged with two variables config_ (of struct Config), and engine_ (of struct
   // Engine), config_ is static across all threads, and engine_
   // is thread-local and specific to each thread. The objects that are used for
   // sampling such as 'normal' of class std::normal_distribution are all hidden
   // in the class private section, and actual access to these is done through
   // public methods such as normal() etc.
   // Note: the original design had an init-check on every call to functions like
   // normal(), but I removed that for performance concerns. It means
   // the following is required in any application's main() using the Rand class:
   //  - call either rand::set_master_seed() or rand::use_hardware_entropy(),
   //    first one sets a user-given master seed, the second generates a master
   //    seed from the machine entropy.
   //  - call rand::init_thread() to properly start the random engine.
   // In case of multi-threading, the first is done once, the init_thread() is
   // done within each thread.
   // The random number stream is reproducible by saving the master seed used,
   // and feed it back in using the manual seed-setting with set_master_seed().
   // The stored master seed is not directly accessible and should be retieved
   // using get_master_seed(). For multi-threading, each thread will use a
   // different seed, but deterministically derived, and thus reproducible from,
   // the master seed.

   private:

   struct Engine {
        std::mt19937 engine;
        bool initialized = false;
        std::normal_distribution<double> normal;
        std::normal_distribution<double> stdnormal{0.0, 1.0};
        std::uniform_real_distribution<double> uniform{0.0, 1.0};
        std::uniform_int_distribution<int> uniform_int;
        std::gamma_distribution<double> gamma;
        std::exponential_distribution<double> exponential;
        std::chi_squared_distribution<double> chi_square;
        std::negative_binomial_distribution<int> negbinom;
   };

   struct Config {
        std::mutex mutex;
        bool use_manual_seed = false;
        uint32_t master_seed = 0;
        uint32_t thread_counter = 0;
   };

   static Config config_;
   static thread_local Engine engine_;
   // In a design where every call to a random number generator checks
   // initialization, init_thread() was private, but now it is in the public
   // section so that a main() function can call it.
   // When the initialization checks would be activated again, there is no
   // need to run init_thread() from a main() application, it will then
   // automatically run on the fist call (within thread) to any of the
   // random number generator functions.
   // The config_ and engine_ remain private, and also the master seed;
   // this means a getter function is needed to retrieve the master seed.

public:
   rand() = delete;
   rand(const rand&) = delete;
   rand& operator=(const rand&) = delete;
   static void init_thread();
   static uint32_t get_master_seed();
   static void set_master_seed(uint32_t seed);
   static void use_hardware_entropy();

    // High-performance samplers (Inlined for maximum execution speed).
    // I have commented out the initialization check (which would be called
    // on every call to draw a random number!) to see if it improves
    // performance, in a self-contained software like bayz_cpp it may be OK...

//  There is normal() and stdnormal(), the first takes parameters for mean
//  and stddev, the second has no parameters and makes a standard N(0,1) draw.
    static inline double normal(double mean, double stddev) {
//      if (!engine_.initialized) init_thread();
        using param_t = std::normal_distribution<double>::param_type;
        return engine_.normal(engine_.engine, param_t(mean, stddev));
    }

   static inline double stdnormal() {
//      if (!engine_.initialized) init_thread();
        return engine_.stdnormal(engine_.engine);
    }

//  uniform() is the common real-valued uniform distribution in [0, 1).
//  I don't think I currently have a need for real uniforms with other
//  bounds, and they're easy to make by multiplying, if needed.
    static inline double uniform() {
//        if (!engine_.initialized) init_thread();
        return engine_.uniform(engine_.engine);
    }

    static inline double gamma(double alpha, double beta) {
//        if (!engine_.initialized) init_thread();
        using param_t = std::gamma_distribution<double>::param_type;
        return engine_.gamma(engine_.engine, param_t(alpha, beta));
    }

    static inline double beta(double alpha, double beta_param) {
        double x = gamma(alpha, 1.0);
        double y = gamma(beta_param, 1.0);
        if (x + y == 0.0) return 0.0;
        return x / (x + y);
    }

    static inline double negbinom(int r, double p) {
//        if (!engine_.initialized) init_thread();
        using param_t = std::negative_binomial_distribution<int>::param_type;
        return engine_.negbinom(engine_.engine, param_t(r, p));
    }
    
};

#endif // RAND_H