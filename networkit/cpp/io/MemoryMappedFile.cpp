#include "MemoryMappedFile.h"

namespace NetworKit {
	MemoryMappedFile::MemoryMappedFile() {}

	MemoryMappedFile::MemoryMappedFile(const std::string &path) { open(path); }

	MemoryMappedFile::~MemoryMappedFile() { close(); }
}

#ifdef NETWORKIT_WINDOWS
#include <Windows.h>
#include <fileapi.h>
#include <cassert>

namespace NetworKit {
	struct MemoryMappedFileState {
		HANDLE hFile{ nullptr };
		HANDLE hMap{ nullptr };
	};

	void MemoryMappedFile::open(const std::string& path) {
		if (!state) {
			state.reset(new MemoryMappedFileState{});
		}
		else {
			close();
		}

		state->hFile = CreateFile(path.c_str(), GENERIC_READ, 0, NULL, OPEN_EXISTING, 0, 0);
		if (state->hFile == INVALID_HANDLE_VALUE)
			throw std::runtime_error("Unable to open file");

		size_t size;
		{
			DWORD low, high;
			low = GetFileSize(state->hFile, &high);
			size = (static_cast<size_t>(high) << 32) | low;
		}

		// If the file is empty Mapping is not necessary
		if (!size) {
			assert(beginIt == nullptr && endIt == nullptr);
			return;
		}

		state->hMap = CreateFileMapping(state->hFile, NULL, PAGE_READONLY, 0, 0, NULL);
		if (state->hMap == NULL)
			throw std::runtime_error("Could not map file");

		auto window = MapViewOfFile(state->hMap, FILE_MAP_READ, 0, 0, 0);

		beginIt = reinterpret_cast<char *>(window);
		endIt = beginIt + size;
	}

	void MemoryMappedFile::close() {
		if (beginIt) {
			UnmapViewOfFile(static_cast<LPCVOID>(beginIt));
			beginIt = nullptr;
			endIt = nullptr;
		}

		if (state) {
			if (state->hMap) {
				CloseHandle(state->hMap);
				state->hMap = nullptr;
			}

			if (state->hFile) {
				CloseHandle(state->hFile);
				state->hFile = nullptr;
			}
		}
	}
}

#else

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>

namespace NetworKit {

	struct MemoryMappedFileState {};

	void MemoryMappedFile::open(const std::string& path) {
		if (beginIt) close();

		auto fd = ::open(path.c_str(), O_RDONLY);
		if (fd < 0)
			throw std::runtime_error("Unable to open file");

		struct stat st;
		if (fstat(fd, &st))
			throw std::runtime_error("Could not obtain file stats");

		// It does not really matter if we use a private or shared mapping.
		auto window = mmap(nullptr, st.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
		if (window == reinterpret_cast<void *>(-1))
			throw std::runtime_error("Could not map file");

		if (::close(fd))
			throw std::runtime_error("Error during close()");

		beginIt = reinterpret_cast<char *>(window);
		endIt = beginIt + st.st_size;
	}

	void MemoryMappedFile::close() {
		if (!beginIt)
			return;

		if (munmap(reinterpret_cast<void*>(const_cast<value_type*>(beginIt)),
			std::distance(beginIt, endIt))) {
			throw std::runtime_error("Could not unmap file");
		}
	}
}

#endif
