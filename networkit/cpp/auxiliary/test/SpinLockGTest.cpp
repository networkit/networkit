/*
 * SpinLockGTest.cpp
 *
 *  Created on: 18.07.2025
 *      Author: Andreas Scharf (andreas.b.scharf@gmail.com)
 */

#include <thread>
#include <gtest/gtest.h>
#include <networkit/auxiliary/SpinLock.hpp>

class SpinlockGTest : public testing::Test {};

namespace NetworKit {

TEST_F(SpinlockGTest, testSpinLockSingleThread) {
    Aux::Spinlock lock;
    EXPECT_NO_THROW({
        lock.lock();
        lock.unlock();
    });
}

TEST_F(SpinlockGTest, testSpinLockUnderContention) {
    Aux::Spinlock lock;

    int counter = 0;
    constexpr int incrementsPerThread = 10000;
    const unsigned int numberOfThreads = std::thread::hardware_concurrency();

    auto worker = [&]() {
        for (int i = 0; i < incrementsPerThread; ++i) {
            lock.lock();
            ++counter;
            lock.unlock();
        }
    };

    std::vector<std::thread> threads;
    for (unsigned int i = 0; i < numberOfThreads; ++i) {
        threads.emplace_back(worker);
    }
    for (std::thread &thread : threads) {
        thread.join();
    }

    EXPECT_EQ(counter, numberOfThreads * incrementsPerThread);
}

} // namespace NetworKit
