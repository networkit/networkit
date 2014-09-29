"""
The NetworKit benchmark

This module implements a comprehensive benchmark of NetworKit's analytics algorithms



"""

from networkit import *
import pandas
import unittest
import contextlib
import sys
from unittest.signals import registerResult
import warnings
import time
import math


class BenchmarkRunner(unittest.TextTestRunner):

    def run(self, test):
        "Run the given test case or test suite."
        result = self._makeResult()
        registerResult(result)
        result.failfast = self.failfast
        result.buffer = self.buffer
        with warnings.catch_warnings():
            if self.warnings:
                # if self.warnings is set, use it to filter all the warnings
                warnings.simplefilter(self.warnings)
                # if the filter is 'default' or 'always', special-case the
                # warnings from the deprecated unittest methods to show them
                # no more than once per module, because they can be fairly
                # noisy.  The -Wd and -Wa flags can be used to bypass this
                # only when self.warnings is None.
                if self.warnings in ['default', 'always']:
                    warnings.filterwarnings('module',
                            category=DeprecationWarning,
                            message='Please use assert\w+ instead.')
            startTime = time.time()
            startTestRun = getattr(result, 'startTestRun', None)
            if startTestRun is not None:
                startTestRun()
            try:
                oldResult = result.globalDataFrame
                newResult = test(result)
                #print(newResult)
                #print(type(newResult))
                #print(newResult.globalDataFrame)
            finally:
                stopTestRun = getattr(result, 'stopTestRun', None)
                if stopTestRun is not None:
                    stopTestRun()
            stopTime = time.time()
        timeTaken = stopTime - startTime
        result.printErrors()
        if hasattr(result, 'separator2'):
            self.stream.writeln(result.separator2)
        run = result.testsRun
        self.stream.writeln("Ran %d test%s in %.3fs" %
                            (run, run != 1 and "s" or "", timeTaken))
        self.stream.writeln()

        expectedFails = unexpectedSuccesses = skipped = 0
        try:
            results = map(len, (result.expectedFailures,
                                result.unexpectedSuccesses,
                                result.skipped))
        except AttributeError:
            pass
        else:
            expectedFails, unexpectedSuccesses, skipped = results

        infos = []
        if not result.wasSuccessful():
            self.stream.write("FAILED")
            failed, errored = len(result.failures), len(result.errors)
            if failed:
                infos.append("failures=%d" % failed)
            if errored:
                infos.append("errors=%d" % errored)
        else:
            self.stream.write("OK")
        if skipped:
            infos.append("skipped=%d" % skipped)
        if expectedFails:
            infos.append("expected failures=%d" % expectedFails)
        if unexpectedSuccesses:
            infos.append("unexpected successes=%d" % unexpectedSuccesses)
        if infos:
            self.stream.writeln(" (%s)" % (", ".join(infos),))
        else:
            self.stream.write("\n")
        return result



class BenchmarkSuite(unittest.TestSuite):
    def run(self, result, debug=False):
        topLevel = False
        if getattr(result, '_testRunEntered', False) is False:
            result._testRunEntered = topLevel = True
        #print(enumerate(self))
        for index, test in enumerate(self):
            #print(str(index)+': '+str(test))
            if result.shouldStop:
                break

            if _isnotsuite(test):
                self._tearDownPreviousClass(test, result)
                self._handleModuleFixture(test, result)
                self._handleClassSetUp(test, result)
                result._previousTestClass = test.__class__

                if (getattr(test.__class__, '_classSetupFailed', False) or
                    getattr(result, '_moduleSetUpFailed', False)):
                    continue

            if not debug:
                oldData = result.globalDataFrame
                data = test(result)
                #print(type(data))
                newData = data.globalDataFrame
                newData.append(oldData)
                #print(newData)
                result.globalDataFrame = newData
                #oldData.append(data.globalDataFrame)
                #print("return type of result is: {0}".format(type(data)))
                #print(data.globalDataFrame)
                #data.append(oldData)
                #print(newResult)
#                result.appendDataFrame(data)
                #result.appendDataFrame = data.globalDataFrame
                #if type
                #    resu
            else:
                test.debug()

            if self._cleanup:
                self._removeTestAtIndex(index)

        if topLevel:
            self._tearDownPreviousClass(None, result)
            self._handleModuleTearDown(result)
            result._testRunEntered = False
        #result.globalDataFrame = newResult
        #print(oldData)
        return result

    def debug(self):
        """Run the tests without collecting errors in a TestResult"""
        debug = _DebugResult()
        self.run(debug, True)

def _isnotsuite(test):
    "A crude way to tell apart testcases and suites with duck-typing"
    try:
        iter(test)
    except TypeError:
        return True
    return False

class BenchmarkResult(unittest.TextTestResult):
    def __init__(self, stream=sys.stdout, descriptions=True, verbosity=1):
        super(BenchmarkResult, self).__init__(stream, descriptions, verbosity)
        self.globalDataFrame = pandas.DataFrame()

    def getDataFrame(self):
        return self.globalDataFrame

    def appendDataFrame(self, dataFrame):
        self.globalDataFrame.append(dataFrame)

class SkipTest(Exception):
    """
    Raise this exception in a test to skip it.

    Usually you can use TestCase.skipTest() or one of the skipping decorators
    instead of raising this directly.
    """

class _ShouldStop(Exception):
    """
    The test should stop.
    """

class _UnexpectedSuccess(Exception):
    """
    The test was supposed to fail, but it didn't!
    """

class _Outcome(object):
    def __init__(self, result=None):
        self.expecting_failure = False
        self.result = result
        self.result_supports_subtests = hasattr(result, "addSubTest")
        self.success = True
        self.skipped = []
        self.expectedFailure = None
        self.errors = []

    @contextlib.contextmanager
    def testPartExecutor(self, test_case, isTest=False):
        old_success = self.success
        self.success = True
        try:
            yield
        except KeyboardInterrupt:
            raise
        except SkipTest as e:
            self.success = False
            self.skipped.append((test_case, str(e)))
        except _ShouldStop:
            pass
        except:
            exc_info = sys.exc_info()
            if self.expecting_failure:
                self.expectedFailure = exc_info
            else:
                self.success = False
                self.errors.append((test_case, exc_info))
            # explicitly break a reference cycle:
            # exc_info -> frame -> exc_info
            exc_info = None
        else:
            if self.result_supports_subtests and self.success:
                self.errors.append((test_case, None))
        finally:
            self.success = self.success and old_success


class BenchmarkTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def run(self, result=None):
        orig_result = result
        if result is None:
            print("got no result object")
            result = self.defaultTestResult()
            startTestRun = getattr(result, 'startTestRun', None)
            if startTestRun is not None:
                startTestRun()

        result.startTest(self)

        testMethod = getattr(self, self._testMethodName)
        if (getattr(self.__class__, "__unittest_skip__", False) or
            getattr(testMethod, "__unittest_skip__", False)):
            # If the class or method was skipped.
            try:
                skip_why = (getattr(self.__class__, '__unittest_skip_why__', '')
                            or getattr(testMethod, '__unittest_skip_why__', ''))
                self._addSkip(result, self, skip_why)
            finally:
                result.stopTest(self)
            return
        expecting_failure = getattr(testMethod,
                                    "__unittest_expecting_failure__", False)
        outcome = _Outcome(result)
        try:
            self._outcome = outcome

            with outcome.testPartExecutor(self):
                self.setUp()
            if outcome.success:
                outcome.expecting_failure = expecting_failure
                with outcome.testPartExecutor(self, isTest=True):
                    # benchmark test cases return a dataframe containing data.
                    #print("result content: {0}".format(str(result.getDataFrame())))
                    #print("returning the pandas value works?")
                    data = testMethod()
                    #print(type(data))
                    oldData = result.globalDataFrame
                    data.append(oldData)
                    result.globalDataFrame = data
                    print(data)
                    #result.appendDataFrame(data)
                outcome.expecting_failure = False
                with outcome.testPartExecutor(self):
                    self.tearDown()

            self.doCleanups()
            for test, reason in outcome.skipped:
                self._addSkip(result, test, reason)
            self._feedErrorsToResult(result, outcome.errors)
            if outcome.success:
                if expecting_failure:
                    if outcome.expectedFailure:
                        self._addExpectedFailure(result, outcome.expectedFailure)
                    else:
                        self._addUnexpectedSuccess(result)
                else:
                    result.addSuccess(self)
            # return the test result as well as the result
            #print("is the pandas appended to the result object?")
            #oldFrame = result.globalDataFrame
            #data.append(oldFrame)
            #result.globalDataFrame = data
            #print(result)
            #print(result.getDataFrame())
            #print("---------------")
            return result
        finally:
            result.stopTest(self)
            if orig_result is None:
                stopTestRun = getattr(result, 'stopTestRun', None)
                if stopTestRun is not None:
                    stopTestRun()

            # explicitly break reference cycles:
            # outcome.errors -> frame -> outcome -> outcome.errors
            # outcome.expectedFailure -> frame -> outcome -> outcome.expectedFailure
            outcome.errors.clear()
            outcome.expectedFailure = None

            # clear the outcome, no more needed
            self._outcome = None


class Benchmark_ConnectedComponents(BenchmarkTestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_init(self):
        g = loadGraph('PGPgiantcompo')
        t = stopwatch.Timer()
        cc = properties.ConnectedComponents(g)
        elapsed = t.stop()
        data = pandas.DataFrame(data={'name': ['ConnectedComponents.init'], 'time':[elapsed]})
        return data

    def test_run(self):
        g = loadGraph('PGPgiantcompo')
        cc = properties.ConnectedComponents(g)
        t = stopwatch.Timer()
        cc.run()
        elapsed = t.stop()
        data = pandas.DataFrame(data={'name': ['ConnectedComponents.run'], 'time':[elapsed]})
        return data

class Benchmark_PLM(BenchmarkTestCase):
    def setUp(self):
        pass
    def tearDown(self):
        pass
    def test_init(self):
        t = stopwatch.Timer()
        cc = community.PLM()
        elapsed = t.stop()
        data = pandas.DataFrame(data={'name': ['PLM.init'], 'time':[elapsed]})
        return data

    def test_run(self):
        g = loadGraph('PGPgiantcompo')
        plm = community.PLM()
        t = stopwatch.Timer()
        plm.run(g)
        elapsed = t.stop()
        data = pandas.DataFrame(data={'name': ['PLM.run'], 'time':[elapsed]})
        return data


# helper function

def loadGraph(key, basePath):
    (fileName, formatName) = networks[key]
    G = readGraph(os.path.join(basePath, fileName), formatName)
    return G


class Timer:
    """ Use the Python with-statement to time your code
    with this timer. """

    def __enter__(self):
        self.start = time.clock()
        return self

    def __exit__(self, *args):
        self.end = time.clock()
        self.elapsed = round(self.end - self.start, 6)



# what is a  test

# get graph
# init
# multiple times
#       run
#       take time
#


# settings

nRuns = 5   # how many runs for representative results

# collection of networks

networks = {
            "PGPgiantcompo" : ("PGPgiantcompo.metis.graph", Format.METIS),
            "power" : ("power.metis.graph", Format.METIS),
            "caidaRouterLevel" : ("caidaRouterLevel.graph", Format.METIS),
            "as-22july06" : ("as-22july06.metis.graph", Format.METIS),
            "coAuthorsDBLP" : ("coAuthorsDBLP.graph", Format.METIS),
            "uk-2007-05" : ("uk2007-05.metis.graph", Format.METIS),
            "uk-2002" : ("uk-2002.metis.graph", Format.METIS),
            "fb-Texas84" : ("Texas84.edgelist", Format.EdgeListTabZero),
            "fb-Caltech36" : ("Caltech36.edgelist", Format.EdgeListTabZero),
            "fb-MIT8" : ("MIT8.edgelist", Format.EdgeListTabZero),
            "fb-Smith60" : ("Smith60.edgelist", Format.EdgeListTabZero),
            "con-fiber_big" : ("con-fiber_big.metis.graph", Format.METIS),
            }



# - connected components (properties.ConnectedComponents, properties.ParallelConnectedComponents)

selected = ["PGPgiantcompo", "power"]
collectionDir = os.path.expanduser("~/workspace/Data/NwkBenchmark")


class Algo:
    """ runner for an algorithm"""
    def run(self, G):
        raise Exception("Not implementedd")

class ConnectedComponents_(Algo):
    name = "ConnectedComponents"

    def run(self, G):
        cc = properties.ConnectedComponents(G)
        cc.run()
        # TODO: validate result?

class ParallelConnectedComponents_(Algo):
    name = "ParallelConnectedComponents"

    def run(self, G):
        cc = properties.ParallelConnectedComponents(G)
        cc.run()


class CoreDecomposition_(Algo):
    name = "CoreDecomposition"

    def run(self, G):
        cd = properties.CoreDecomposition(G)
        cd.run()


class Fail_(Algo):
    name = "Fail"

    def run(self, G):
        raise Exception("FAIL!")



def main():
    logging.info("start benchmark")

    data = benchmark(ConnectedComponents_(), ["PGPgiantcompo", "power"])
    data = benchmark(CoreDecomposition_(), ["PGPgiantcompo", "power"])
    #data = benchmark(_Fail(), ["PGPgiantcompo", "power    "])


def benchmark(algo, graphs):
    table = []  # list of dictionaries, to be converted to a DataFrame

    for graphName in graphs:
        try:
            print("")
            logging.info("loading {0}".format(graphName))
            G = loadGraph(graphName, basePath=collectionDir)
            try:
                for i in range(nRuns):
                    row = {}    # benchmark data row
                    with Timer() as t:
                        logging.info("running {algo.name}".format(**locals()))
                        algo.run(G)
                    logging.info("took {0} s".format(t.elapsed))
                    # store data
                    row["algo"] = algo.name
                    row["graph"] = graphName
                    row["time"] = t.elapsed
                    table.append(row)
            except Exception as ex:
                logging.error("algorithm {algo.name} failed with exception: {ex}".format(**locals()))
        except Exception as ex:
            logging.error("loading graph {graphName} failed with exception: {ex}".format(**locals()))

    return pandas.DataFrame(table)

# - degree distribution power-law estimation (properties.powerLawExponent)
# - k-core decomposition (properties.CoreDecomposition)
# - BFS & Dijkstra (graph.BFS, graph.Dijkstra)
# - clustering coefficients (average local), exact (properties.ClusteringCoefficient.avgLocal) and approximated (properties.ClusteringCoefficient.approxAvgLocal)
# - degree assortativity (properties.degreeAssortativity)
# - diameter, exact (properties.Diameter.exactDiameter) and estimate (properties.Diameter.estimatedDiameterRange)
# - community detection (community.PLM, community.PLP)
# - centrality
# 	- PageRank (centrality.PageRank, centrality.SciPyPageRank)
# 	- Eigenvector centrality (centrality.EigenvectorCentrality, centrality.SciPyEVZ)
# 	- betweenness,  exact (centrality.Betweenness) and approximated (centrality.ApproxBetweenness, centrality.ApproxBetweenness2)
# - generators
# 	- Erdös-Renyi (generators.ErdosRenyiGenerator)
# 	- Barabasi-Albert (generators.BarabasiAlbertGenerator)
# 	- Chung-Lu (generators.ChungLuGenerator)
# 	- RMAT (generators.RmatGenerator)

if __name__ == "__main__":
    main()
