/*
 * SortedListGTest.cpp
 *
 *  Created on: 22.07.2025
 *      Author: Andreas Scharf (andreas.b.scharf@gmail.com)
 */

#include <gtest/gtest.h>
#include <networkit/auxiliary/Timer.hpp>

class LoggingTimerTest : public ::testing::Test {
protected:
    std::streambuf *originalCerrBuf;
    std::stringstream capturedOutput;

    void SetUp() override {
        originalCerrBuf = std::cerr.rdbuf();
        std::cerr.rdbuf(capturedOutput.rdbuf());
    }

    void TearDown() override { std::cerr.rdbuf(originalCerrBuf); }

    std::string getCapturedOutput() { return capturedOutput.str(); }
};

TEST_F(LoggingTimerTest, testLogsWhenLabelIsUsed) {
    Aux::Log::setLogLevel("TRACE");
    {
        Aux::LoggingTimer timer("test-timer", Aux::Log::LogLevel::TRACE);
        std::this_thread::sleep_for(std::chrono::milliseconds(5));
    }
    const std::string output = getCapturedOutput();
    EXPECT_NE(output.find("Timer \"test-timer\" ran for"), std::string::npos);
}

TEST_F(LoggingTimerTest, testLogsWhenLabelIsNotUsed) {
    Aux::Log::setLogLevel("DEBUG");
    {
        Aux::LoggingTimer timer("", Aux::Log::LogLevel::DEBUG);
        std::this_thread::sleep_for(std::chrono::milliseconds(3));
    }
    const std::string output = getCapturedOutput();
    EXPECT_NE(output.find("Timer ran for"), std::string::npos);
}

TEST_F(LoggingTimerTest, testStartedTimerElapsedValuesBeforeStop) {
    Aux::StartedTimer timer;
    std::this_thread::sleep_for(std::chrono::milliseconds(10));

    const auto ms = timer.elapsedMilliseconds();
    const auto us = timer.elapsedMicroseconds();
    const auto ns = timer.elapsedNanoseconds();

    EXPECT_GE(ms, 5);
    EXPECT_GE(us, 5000);
    EXPECT_GE(ns, 5000000);
}

TEST_F(LoggingTimerTest, testStartedTimerElapsedValuesAfterStop) {
    Aux::StartedTimer timer;
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
    timer.stop();
    std::this_thread::sleep_for(std::chrono::milliseconds(50));

    const auto ms2 = timer.elapsedMilliseconds();
    const auto us2 = timer.elapsedMicroseconds();
    const auto ns2 = timer.elapsedNanoseconds();

    EXPECT_GE(ms2, 5);
    EXPECT_GE(us2, 5000);
    EXPECT_GE(ns2, 5000000);

    EXPECT_LE(ms2, us2 / 1000 + 1);
    EXPECT_LE(us2, ns2 / 1000 + 1);
}

TEST_F(LoggingTimerTest, testStartedTimerStartAndStopTimeValid) {
    Aux::StartedTimer timer;
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
    timer.stop();

    const auto start = timer.startTime();
    const auto stop = timer.stopTime();

    EXPECT_LE(start, stop);
}
