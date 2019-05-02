
#include <vector>
#include <string>
#include <cmath>
#define M_PI       3.14159265358979323846   // pi

#include <girgs/Node.h>
#include <girgs/girgs_api.h>

namespace girgs {

GIRGS_API double calculateRadius(int n, double alpha, double T, int deg);
GIRGS_API double hyperbolicDistance(double r1, double phi1, double r2, double phi2);

GIRGS_API std::vector<double> sampleRadii(int n, double alpha, double R, int weightSeed);
GIRGS_API std::vector<double> sampleAngles(int n, int positionSeed);

static double radiusToGirgWeight(double r, double R) { return std::exp((R - r) / 2); }
static double girgWeightToRadius(double w, double R, double scaling = 1.0) { return R - 2 * std::log(w / scaling); }

static double angleToGirgPosition(double angle) { return angle / 2 / M_PI; }
static double girgPositionToAngle(double position) { return position * 2 * M_PI; }

} // namespace girgs
