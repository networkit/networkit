/**
 * BinaryPartitionWriter.h
 *
 * @author Michael Hamann <michael.hamann@kit.edu>
 */

#ifndef BINARYPARTITIONWRITER_H_
#define BINARYPARTITIONWRITER_H_

#include "../structures/Partition.h"
#include <string>

namespace NetworKit {

/**
* @ingroup io
* Writes a partition to a file to contains a binary list of partition ids.
* Partition ids are unsigned integers.
*/
class BinaryPartitionWriter {
public:
	/**
	 * Constructs the BinaryPartitionWriter class using unsigned integers
	 * of width @a width for the written partition ids.
	 *
	 * @param[in] width The width of the written integers (supported values: 4, 8, default: 4).
	 */
	BinaryPartitionWriter(uint8_t width = 4);

	/**
	 * Write the given partition @a zeta to the given @a path.
	 *
	 * @param[in] zeta The partition to write.
         * @param[in] path The path to write to.
	 */
	virtual void write(const Partition& zeta, const std::string& path) const;
private:
	uint8_t width;
};
}

#endif
