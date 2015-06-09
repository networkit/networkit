/*
 * MatrixMarketReader.cpp
 *
 *  Created on: 25.07.2014
 *      Author: dhoske
 */

#include "MatrixMarketReader.h"

#include <cctype>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <iterator>
#include <unordered_set>

namespace NetworKit {

namespace {
	static const char COMMENT_CHAR = '%';
	static const std::string MAGIC = "%%matrixmarket";

	std::string tolower(const std::string& str) {
		std::string out;
		std::transform(str.begin(), str.end(), std::back_inserter(out), ::tolower);
		return out;
	}
}

CSRMatrix MatrixMarketReader::read(const std::string &path) {
	std::ifstream in(path);
	if (!in) {
		throw std::runtime_error("could not open: " + path);
	}
	return read(in);
}

CSRMatrix MatrixMarketReader::read(std::istream& in) {
	enum { FIRST_LINE, HEADER, ENTRIES } state = FIRST_LINE;

	count nrows, ncols, nzeroes;
	std::vector<std::pair<index,index>> positions;
	std::vector<double> values;
	std::string line;
	bool weighted = true;
	while (getline(in, line)) {
		std::istringstream in_line(line);

		if (state == FIRST_LINE) {
			// First line: magic <keyword>*

			// Read id of file format at start of file
			std::string magic, object, format, data, qualifier;
			in_line >> magic >> object >> format >> data >> qualifier;
			magic = tolower(magic);
			object = tolower(object);
			format = tolower(format);
			data = tolower(data);
			qualifier = tolower(qualifier);
			if (magic != MAGIC) {
				throw std::runtime_error(MAGIC + " not found");
			}
			if (format != "coordinate") {
				throw std::runtime_error("Unsupported format: " + format);
			}
			if (data == "pattern") {
				weighted = false;
			} else if (data != "real") {
				throw std::runtime_error("Unsupported data type: " + data);
			}

			state = HEADER;
		} else if (line.empty() || line[0] == COMMENT_CHAR) {
			// Skip line if empty or comment
			continue;
		} else if (state == HEADER) {
			// Header: nrows, ncols, nzeroes
			in_line >> nrows >> ncols >> nzeroes;
			if (in_line.fail()) {
				throw std::runtime_error("expected three non-negative integers in header line:" + line);
			}

			positions = std::vector<std::pair<index,index>>(nzeroes);
			values = std::vector<double>(nzeroes);
			state = ENTRIES;
		} else {
			// Entry: row col value
			count i, j;
			double val;
			in_line >> i >> j;
			i--;
			j--;
			if (weighted) {
				in_line >> val;
			} else {
				val = 1.0;
			}

			if (in_line.fail()) {
				throw std::runtime_error("expected three entries in entry line: " + line);
			}
			if (i >= nrows) {
				throw std::runtime_error("invalid index: " + std::to_string(i));
			}
			if (j >= ncols) {
				throw std::runtime_error("invalid index: " + std::to_string(j));
			}

			positions[nzeroes] = std::make_pair(i,j);
			values[nzeroes] = val;
			nzeroes--;
		}
	}

	if (nzeroes != 0) {
		WARN("wrong number of non-zeroes given");
	}

	return CSRMatrix(nrows, ncols, positions, values);
}

}
