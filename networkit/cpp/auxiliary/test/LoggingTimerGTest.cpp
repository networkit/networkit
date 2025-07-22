/*
 * SortedListGTest.cpp
 *
 *  Created on: 22.07.2025
 *      Author: Andreas Scharf (andreas.b.scharf@gmail.com)
 */

#include <thread>
#include <gmock/gmock-matchers.h>
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
    EXPECT_THAT(getCapturedOutput(), testing::HasSubstr("Timer \"test-timer\" ran for"));
}

TEST_F(LoggingTimerTest, testLogsWhenLabelIsNotUsed) {
    Aux::Log::setLogLevel("DEBUG");
    {
        Aux::LoggingTimer timer("", Aux::Log::LogLevel::DEBUG);
        std::this_thread::sleep_for(std::chrono::milliseconds(3));
    }
    EXPECT_THAT(getCapturedOutput(), testing::HasSubstr("Timer ran for"));
}

TEST_F(LoggingTimerTest, testStartedTimerElapsedValuesBeforeStop) {
    Aux::StartedTimer timer;
    std::this_thread::sleep_for(std::chrono::milliseconds(10));

    EXPECT_GE(timer.elapsedMilliseconds(), 5);
    EXPECT_GE(timer.elapsedMicroseconds(), 5000);
    EXPECT_GE(timer.elapsedNanoseconds(), 5000000);
}

TEST_F(LoggingTimerTest, testStartedTimerElapsedValuesAfterStop) {
    Aux::StartedTimer timer;
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
    timer.stop();
    std::this_thread::sleep_for(std::chrono::milliseconds(50));

    EXPECT_LE(timer.elapsedMilliseconds(), timer.elapsedMicroseconds() / 1000 + 1);
    EXPECT_LE(timer.elapsedMicroseconds(), timer.elapsedNanoseconds() / 1000 + 1);
}

TEST_F(LoggingTimerTest, testStartedTimerStartAndStopTimeValid) {
    Aux::StartedTimer timer;
    std::this_thread::sleep_for(std::chrono::milliseconds(50));
    timer.stop();
    EXPECT_LE(timer.startTime(), timer.stopTime());
}
