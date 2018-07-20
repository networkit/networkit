/*
* MemoryMappedFileGTest.h
*
*  Created on: 28.07.2018
*      Author: Manuel Penschuck (networkit@manuel.jetzt)
*/


#include <algorithm>
#include <fstream>
#include <limits>
#include <memory>
#include <random>
#include <string>
#include <stdexcept>

#include "MemoryMappedFileGTest.h"
#include "../MemoryMappedFile.h"

namespace NetworKit {
	// This object generates a temporary file and fills it with a request number
	// of random bytes. On destruction of the elemen, the file is deleted automatically.
	// Make sure that no application holds a handle to the file at this point, since
	// on Windows removal of the file fails silently in this case.
	class TemporaryFile {
	public:
		using value_type = char;
		using const_iterator = const value_type*;

		TemporaryFile(size_t bytes, std::default_random_engine& prng) : 
			dataBegin(new value_type[bytes]),
			dataEnd(dataBegin.get() + bytes)
		{
			// obtain temporary file name (C++17 offers std::filesystem, which we currently cannot rely on)
			while (true) {
				std::uniform_int_distribution<unsigned> distr;
				path = "tempfile_" + std::to_string(distr(prng));

				std::ifstream f(path);
				if (!f.good()) break;
			}

			// fill array with random data
			{
				// we use int (rather than char), as char is not an integer type for MSVC
				std::uniform_int_distribution<int> distr{
					std::numeric_limits<value_type>::min(),
					std::numeric_limits<value_type>::max() };

				std::generate(dataBegin.get(), dataEnd, [&distr, &prng] {
					return static_cast<value_type>(distr(prng)); });
			}

			// write array to file
			{
				std::ofstream file(path, std::ios::binary);
				if (bytes)
					file.write(reinterpret_cast<char*>(dataBegin.get()), bytes);
			}
		}

		~TemporaryFile() {
			if (!path.empty()) {
				std::remove(path.c_str());
			}
		}

		const std::string& filename() const { return path; }
		const_iterator cbegin() const { return dataBegin.get(); }
		const_iterator cend() const { return dataEnd; }
		size_t size() const { return dataEnd - dataBegin.get(); }

	private:
		std::unique_ptr<value_type> dataBegin;
		value_type* dataEnd;
		std::string path;
	};

	TEST_P(MemoryMappedFileGTest, ReadFile) {
		const size_t bytes = GetParam();
		std::default_random_engine prng(static_cast<unsigned>(bytes));
		TemporaryFile tmpFile(bytes, prng);
		ASSERT_EQ(bytes, tmpFile.size());

		MemoryMappedFile mmf(tmpFile.filename());
		ASSERT_EQ(bytes, mmf.size());

		for (size_t i = 0; i < bytes; i++) {
			ASSERT_EQ(*(mmf.cbegin() + i), *(tmpFile.cbegin() + i)) << "i=" << i;
		}
	}

	INSTANTIATE_TEST_CASE_P(MemoryMappedFileGTest, MemoryMappedFileGTest,
		::testing::Values(0, 1,
			(1 <<  2) - 1, (1 <<  2), (1 <<  2) + 1,
			(1 << 10) - 1, (1 << 10), (1 << 10) + 1,
			(1 << 16) - 1, (1 << 16), (1 << 16) + 1,
			(1 << 22) - 1, (1 << 22), (1 << 22) + 1
		));
}
