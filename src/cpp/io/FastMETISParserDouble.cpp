/*
 * FastMETISParser.cpp
 *
 *  Created on: 04.10.2013
 *      Author: cls
 */

#include "FastMETISParserDouble.h"

#include <fstream>

namespace NetworKit {

FastMETISParserDouble::FastMETISParserDouble() {
	// TODO Auto-generated constructor stub

}

FastMETISParserDouble::~FastMETISParserDouble() {
	// TODO Auto-generated destructor stub
}

static inline int isupper(char c) {
    return (c >= 'A' && c <= 'Z');
}

static inline int isalpha(char c) {
    return ((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z'));
}

static inline int isspace(char c) {
    return (c == ' ' || c == '\t' || c == '\n' || c == '\12');
}

static inline int isdigit(char c) {
    return (c >= '0' && c <= '9');
}

//static std::pair<count,index> myStrtoul(std::string &str, index i) {	
static count myStrtoul(std::string &str, index &i) {
	register uint64_t acc;
	register int c;
	register uint64_t cutoff;
	register int neg = 0, any, cutlim;

	/*
	 * See strtol for comments as to the logic used.
	 */
	// ignore white spaces
	do {
		c = str[i++];
	} while (isspace(c));
	// parse sign bit
	if (c == '-') {
		neg = 1;
		c = str[i++];
	} else if (c == '+')
		c = str[i++];

	// depending on base, determine allowed characters
	// assume base is decimal, therefore ignore the code
	register int base = 10;
/*	if ((base == 0 || base == 16) &&
		c == '0' && (str[i] == 'x' || str[i] == 'X')) {
		c = s[1];
		s += 2;
		base = 16;
	} else if ((base == 0 || base == 2) &&
		c == '0' && (*s == 'b' || *s == 'B')) {
		c = s[1];
		s += 2;
		base = 2;
	}
	if (base == 0)
		base = c == '0' ? 8 : 10;*/
	cutoff = (uint64_t)std::numeric_limits<uint64_t>::max() / (uint64_t)base;
	cutlim = (uint64_t)std::numeric_limits<uint64_t>::max() % (uint64_t)base;
	for (acc = 0, any = 0; ; c = str[i++]) { // i < end necessary?
		if (isdigit(c))
			c -= '0';
// base is 10, so only accept digits
//		else if (isalpha(c))
//			c -= isupper(c) ? 'A' - 10 : 'a' - 10;
		else
			break;
		if (c >= base)
			break;
		//explanation: if any of the "break" cases happens, set any=-1 so max_val will be returned
		if (any < 0 || acc > cutoff || acc == cutoff && c > cutlim)
			any = -1;
		else {
			any = 1;
			acc *= base;
			acc += c;
		}
	}
	if (any < 0) {
		acc = std::numeric_limits<uint64_t>::max();
//		errno = ERANGE;
	} else if (neg)
		acc = -acc;
//	return std::make_pair(acc,i);
	return acc;
}

//static std::pair<double, index> myStrtod(std::string &str, index i) {
static double myStrtod(std::string &str, index &i) {
	double maxExponent = 511;
	static double powersOf10[] = {
		10.,
		100.,
		1.0e4,
		1.0e8,
		1.0e16,
		1.0e32,
		1.0e64,
		1.0e128,
		1.0e256
	};
	index oldi = i;
	int sign, expSign = 0;
	double fraction, dblExp, *d;
	//register CONST char *p;
//	char p;
	register int c;
	int exp = 0;		/* Exponent read from "EX" field. */
	int fracExp = 0;		/* Exponent that derives from the fractional
						* part.  Under normal circumstatnces, it is
						* the negative of the number of digits in F.
						* However, if I is very long, the last digits
						* of I get dropped (otherwise a long I with a
						* large negative exponent could cause an
						* unnecessary overflow on I alone).  In this
						* case, fracExp is incremented one for each
						* dropped digit. */
	int mantSize;		/* Number of digits in mantissa. */
	int decPt;			/* Number of mantissa digits BEFORE decimal
						 * point. */
	uint64_t pExp;		/* Temporarily holds location of exponent
						* in string. */
	/*
	* Strip off leading blanks and check for a sign.
	*/
	//p = string;
	while (isspace(str[i])) {
		++i;
	}
	if (str[i] == '-') {
		sign = 1;
		++i;;
	} else {
		if (str[i] == '+') {
			++i;
		}
		sign = 0;
	}
	/*
	 * Count the number of digits in the mantissa (including the decimal
	 * point), and also locate the decimal point.
	 */
	decPt = -1;
	for (mantSize = 0; ; mantSize += 1) {
		c = str[i];
		if (!isdigit(c)) {
			if ((c != '.') || (decPt >= 0)) {
				break;
			}
		decPt = mantSize;
		}
		++i;
	}
	/*
	 * Now suck up the digits in the mantissa.  Use two integers to
	 * collect 9 digits each (this is faster than using floating-point).
	 * If the mantissa has more than 18 digits, ignore the extras, since
	 * they can't affect the value anyway.
	 */
	pExp  = i;
	i -= mantSize;
	if (decPt < 0) {
		decPt = mantSize;
	} else {
		mantSize -= 1;			/* One of the digits was the point. */
	}
	if (mantSize > 18) {
		fracExp = decPt - 18;
		mantSize = 18;
	} else {
		fracExp = decPt - mantSize;
	}
	if (mantSize == 0) {
		fraction = 0.0;
		i = 0; // this is not necessarily correct
		goto done;
	} else {
		int frac1, frac2;
		frac1 = 0;
		for ( ; mantSize > 9; mantSize -= 1) {
			c = str[i];
			++i;
			if (c == '.') {
				c = str[i];
				++i;
			}
			frac1 = 10*frac1 + (c - '0');
		}
		frac2 = 0;
		for (; mantSize > 0; mantSize -= 1) {
			c = str[i];
			++i;
			if (c == '.') {
				c = str[i];
				++i;
			}
			frac2 = 10*frac2 + (c - '0');
		}
		fraction = (1.0e9 * frac1) + frac2;
	}
	/*
	 * Skim off the exponent.
	*/
	i = pExp;
	if ((str[i] == 'E') || (str[i] == 'e')) {
		++i;
		if (str[i] == '-') {
			expSign = 1;
			++i;
		} else {
			if (str[i] == '+') {
				++i;
		}
			expSign = 0;
		}
		if (!isdigit(str[i])) {
			i = pExp;
			goto done;
		}
		while (isdigit(str[i])) {
			exp = exp * 10 + (str[i] - '0');
			++i;
		}
	}
	if (expSign) {
		exp = fracExp - exp;
	} else {
		exp = fracExp + exp;
	}
	/*
	 * Generate a floating-point number that represents the exponent.
	 * Do this by processing the exponent one bit at a time to combine
	 * many powers of 2 of 10. Then combine the exponent with the
	 * fraction.
	 */
	if (exp < 0) {
		expSign = 1;
		exp = -exp;
	} else {
		expSign = 0;
	}
	if (exp > maxExponent) {
		exp = maxExponent;
//	errno = ERANGE;
	}
	dblExp = 1.0;
	for (d = powersOf10; exp != 0; exp >>= 1, d += 1) {
		if (exp & 01) {
			dblExp *= *d;
		}
	}
	if (expSign) {
		fraction /= dblExp;
	} else {
		fraction *= dblExp;
	}
	done:
//    if (endPtr != NULL) {
//	*endPtr = (char *) p;
//    }
	if (sign) {
		//return std::make_pair(-fraction,i);
		return -fraction;
	}
	//return std::make_pair(fraction,i);
	return fraction;
}
#if 1
static inline std::tuple<count, count, int> parseHeader(const std::string& header) {
	count n;
	count m;
	int flag;

	std::vector<std::string> parts = Aux::StringTools::split(header);
	n = std::stoi(parts[0]);
	m = std::stoi(parts[1]);
	if (parts.size() > 2) {
		flag = std::stoi(parts[2]);
	} else {
		flag = 0;
	}


	return std::make_tuple(n, m, flag);
}


NetworKit::Graph FastMETISParserDouble::parse(const std::string& path) {
	std::ifstream stream(path);
	std::string line;
	std::getline(stream, line); // get header
	count n;
	count m;
	int flag; // weight flag
	std::tie(n, m, flag) = parseHeader(line);
	if (flag > 1) return Graph(0); // return empty graph in case of weighted nodes
	bool weighted = flag % 10;
	Graph G(n, weighted);

	std::string graphName = Aux::StringTools::split(Aux::StringTools::split(path, '/').back(), '.').front();

	G.setName(graphName);

	node u = 0;

	// XXX: if the first character of a line is the % char, comment lines will be ignored

	std::string current;
	node v = 0;
	if (!weighted) {
		DEBUG("reading unweighted graph");
		// unweighted edges
		while (std::getline(stream, line)) {
			if (line.empty() || (!line.empty() && line[0] == '%') ) {
				continue;
			}
			// determine index of last character of last node
			count end = line.length()-1;
			while (line[end] == ' ' && end > 0) {
				--end;
			}

			index i = 0;
			while(i < end) {
				//std::tie(v,i) = myStrtoul(line,i);
				v = myStrtoul(line,i);
				if ( u < v ) G.addEdge(u,v-1);
			}
			++u;
		}
	} else {
		DEBUG("reading weighted graph");
		double weight = 0.0;
		// weighted edges - WARNING: supports only non-negative integer weights
		while (std::getline(stream, line)) {
			if (line.empty() || (!line.empty() && line[0] == '%') ) {
				continue;
			}
			// determine index of last character of last node
			count end = line.length()-1;
			while (line[end] == ' ' && end > 0) {
				--end;
			}

			index i = 0;
			while(i < end) {
				//std::tie(v,i) = myStrtoul(line,i);
				//DEBUG("checkpoint: parsing oul was successfull");
				//std::tie(weight,i) = myStrtod(line,i);
				//DEBUG("checkpoint: parsing double was successfull with i=", i, " while end=", end);
				v = myStrtoul(line,i);
				weight = myStrtod(line,i);
				if ( u < v ) G.addEdge(u,v-1,weight);
			}
			++u;
		}
	}


	return G;
}

#else 
static inline uint64_t fast_string_to_integer(std::string::iterator it, const std::string::iterator& end) {
	uint64_t val = *it - '0';	// works on ASCII code points
	++it;
	while (it != end) {
		val *= 10;
		val += *it  - '0';
		++it;
	}
	return val;
}

static inline double fast_string_to_double(std::string::iterator it, const std::string::iterator& end) {
	double val = 0.0;
	bool negative = false;
	if (*it == '-') {
		negative = true;
		++it;
	}
	while (*it >= '0' && *it <= '9' && it != end) { //it != end
		val = (val*10.0) + (*it - '0');
		++it;
	}
	if (*it == '.' && it != end) {
		double afterComma = 0.0;
		int exp = 0;
		++it;
		while (*it >= '0' && *it <= '9' && it != end) { //it != end
			afterComma = (afterComma*10.0) + (*it - '0');
			++exp;
			++it;
		}
		val += afterComma / std::pow(10.0, exp);
	}
	if (negative) {
		val = -val;
	}
	return val;
}

static inline std::tuple<count, count, int> parseHeader(const std::string& header) {
	count n;
	count m;
	int flag;

	std::vector<std::string> parts = Aux::StringTools::split(header);
	n = std::stoi(parts[0]);
	m = std::stoi(parts[1]);
	if (parts.size() > 2) {
		flag = std::stoi(parts[2]);
	} else {
		flag = 0;
	}


	return std::make_tuple(n, m, flag);
}


NetworKit::Graph FastMETISParserDouble::parse(const std::string& path) {
	std::ifstream stream(path);
	std::string line;
	std::getline(stream, line); // get header
	count n;
	count m;
	int flag; // weight flag
	std::tie(n, m, flag) = parseHeader(line);
	bool weighted = flag % 10;
	Graph G(n, weighted);
	node u = 0;

	// TODO: handle comment lines

	if (flag == 0) {
		DEBUG("reading unweighted graph");
		// unweighted edges
		while (std::getline(stream, line)) {
			if (line.empty()) {
				continue;
			}
			auto it1 = line.begin();
			auto end = line.end();
			auto it2 = std::find(it1, end, ' ');
			if (line.back() == ' ') { // if line ends with one white space, do this trick...
				--end;
			}
			while (true) {
				node v = (fast_string_to_integer(it1, it2) - 1);
				if (u < v) {
					G.addEdge(u, v);
				}
				if (it2 == end) {
					break;
				}
				++it2;
				it1 = it2;
				it2 = std::find(it1, end, ' ');
			}
			++u; // next node
		}
	} else if (flag == 1) {
		DEBUG("reading weighted graph");
		// weighted edges - WARNING: supports only non-negative integer weights
		while (std::getline(stream, line)) {
			if (line.empty()) {
				continue;
			}
			std::stringstream linestream(line);
			node v;
			double weight;
			while (linestream >> v >> weight) {
				if (u < v) {
					G.addEdge(u,v-1, (edgeweight) weight);
				}
			}
/*			auto it1 = line.begin();
			auto end = line.end();
			auto it2 = std::find(it1, end, ' ');
			if (line.back() == ' ') { // if line ends with one white space, do this trick...
				--end;
			}
			while (true) {
				TRACE("try to read new node");
				node v = (fast_string_to_integer(it1, it2) - 1);
				TRACE("read node " , v);
				// advance
				++it2;
				it1 = it2;
				it2 = std::find(it1, end, ' ');
				TRACE("char at it1: ",*it1," and at it2: ",*it2," and their difference: ",&*it2-&*it1);
				// get weight
				double weight = (fast_string_to_double(it1, it2));
				//double weight = std::stod(line.sub);
				TRACE("read weight " , weight);
				++it2;
				it1 = it2;
				TRACE("find new it2");
				it2 = std::find(it1, end, ' ');
				TRACE("found new it2");
				TRACE("char at it1: ",*it1," and at it2: ",*it2," and their difference: ",&*it2-&*it1);
			
				if (u < v) {
					G.addEdge(u, v, (edgeweight) weight);
				}
				TRACE("new node added");
				if (it2 == end) {
					break;
				}
			}*/
			++u; // next node
		}
		
	} else if (flag == 11) {
		ERROR("weighted nodes are not supported");
		throw std::runtime_error("weighted nodes are not supported");
	} else {
		ERROR("invalid weight flag in header of METIS file: " , flag);
		throw std::runtime_error("invalid weight flag in header of METIS file");
	}


	return G;
}
#endif

} /* namespace NetworKit */
