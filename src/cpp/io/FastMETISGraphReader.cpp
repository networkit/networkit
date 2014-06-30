/*
 * FastMETISGraphReader.cpp
 *
 *  Created on: 24.04.2014
 *      Author: Maximilian Vogel
 */

#include "FastMETISGraphReader.h"
#include "../auxiliary/Log.h"

#include <fstream>

namespace NetworKit {

// a couple of checks on a character
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

/**
 * this is a custom version of an implementation of C's function 
 * 'strtoul' (string to unsigned long).
 * However, it has been adapted to work on the reference and also 
 * writes the first index after the number back to i 
 */
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
		if (any < 0 || acc > cutoff || (acc == cutoff && c > cutlim))
			any = -1;
		else {
			any = 1;
			acc *= base;
			acc += c;
		}
	}
	if (any < 0) {
		acc = std::numeric_limits<uint64_t>::max();
	} else if (neg)
		acc = -acc;
	return acc;
}

/**
 * this is a custom version of an implementation of C's function 
 * 'strtoul' (string to unsigned long).
 * However, it has been adapted to work on the reference and also 
 * writes the first index after the number back to i 
 */
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
	int sign, expSign = 0;
	double fraction, dblExp, *d;
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
		c = str[i++];
		while(!isspace(c)) c = str[i++]; // this is not necessarily correct
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
			c = str[i++];
			while(!isspace(c)) c = str[i++]; // this is not necessarily correct
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
	if (sign) {
		return -fraction;
	}
	return fraction;
}

static inline std::tuple<count, count, int> parseHeader(const std::string& header) {
	count n;
	count m;
	int flag;
	// NOTE: it's possible to use the custom strtoul function here aswell
	// however, it's not used, since parsing two unsigned longs is not performance relevant
	std::vector<std::string> parts = Aux::StringTools::split(header);
	n = std::stoul(parts[0]);
	m = std::stoul(parts[1]);

	if (parts.size() > 2) {
		flag = std::stoi(parts[2]);
	} else {
		flag = 0;
	}


	return std::make_tuple(n, m, flag);
}


Graph FastMETISGraphReader::read(const std::string& path) {
	std::ifstream stream(path);
	if (!stream.is_open()) {
		ERROR("invalid graph file: " , path);
		throw std::runtime_error("invalid graph file");
	}
	std::string line;
	count n;
	count m;
	int flag; // weight flag
	do {
		std::getline(stream, line); // get header
	} while (line[0] == '%');
	std::tie(n, m, flag) = parseHeader(line);
	if (flag > 1) return Graph(0); // return empty graph in case of weighted nodes
	bool weighted = flag % 10;
	Graph G(n, weighted);

	// set the name of the graph
	std::string graphName = Aux::StringTools::split(Aux::StringTools::split(path, '/').back(), '.').front();
	G.setName(graphName);

	// XXX: if the first character of a line is the '%' char, ignore the line as it is a comment
	node u = 0;
	node v = 0;
	if (!weighted) {
		DEBUG("reading unweighted graph");
		// unweighted edges
		while (std::getline(stream, line)) {
			if (!line.empty() && line[0] != '%') {
				count end = line.length();
				if (end > 1) --end; // only decrement if we have more than one edge
				// determine index of last character of last node
				while (line[end] == ' ' && end > 0) {
					--end;
				}
				index i = 0;
				while(i < end) {
					v = myStrtoul(line,i);
					if ( u < v && v > 0) G.addEdge(u,v-1);
				}
				DEBUG("node ", u, " completely parsed");
			}
			++u;
		}
	} else {
		DEBUG("reading weighted graph");
		double weight = 0.0;
		// weighted edges - reads edge weights as double
		while (std::getline(stream, line)) {
			// if line is note empty and not a comment line
			if (!line.empty() && line[0] != '%') {
				count end = line.length();
				if (end > 1) --end; // should never play a role.
				// determine index of last character of last node
				while (line[end] == ' ' && end > 0) {
					--end;
				}
				index i = 0;
				while(i < end) {
					v = myStrtoul(line,i);
					weight = myStrtod(line,i);
					if ( u < v && v > 0) G.addEdge(u,v-1,weight);
				}
				DEBUG("node ", u, " completely parsed");
			}
			++u;
		}
	}

	G.shrinkToFit();
	return G;
}

} /* namespace NetworKit */
