/*
 * MemoryMappedFile.h
 *
 *  Created on: 16.07.2018
 *      Author: Manuel Penschuck (networkit@manuel.jetzt)
 */

#ifndef MEMORY_MAPPED_FILE_H_
#define MEMORY_MAPPED_FILE_H_

#include <iterator>
#include <memory>
#include <string>

namespace NetworKit {
struct MemoryMappedFileState;

/**
 * This class is a wrapper to os-dependend file mapping implementations.
 * It supports to map a read-only file into virtual address space.
 * After successful opening, pointers to the begin/end of the virtual
 * address space are accessible via cbegin() / cend(). All file and mapping
 * handlers are freed automatically.
 */
class MemoryMappedFile {
public:
    using value_type = char;
    using const_iterator = const value_type*;

	MemoryMappedFile();
	explicit MemoryMappedFile(const std::string& path);
	~MemoryMappedFile();

	MemoryMappedFile(const MemoryMappedFile&) = delete;
	MemoryMappedFile& operator=(const MemoryMappedFile&) = delete;

    void open(const std::string&);
    void close();

    const_iterator cbegin() const {return beginIt;}
    const_iterator cend() const {return endIt;}
    size_t size() const { return std::distance(beginIt, endIt); }

private:
    const_iterator beginIt {nullptr};
    const_iterator endIt {nullptr};

    std::unique_ptr<MemoryMappedFileState> state; // use by windows implementation to keep handlers around.
};

}


#endif // MEMORY_MAPPED_FILE_H_