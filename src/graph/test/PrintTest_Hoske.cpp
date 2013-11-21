#ifndef NOGTEST

#include "PrintTest_Hoske.h"

namespace NetworKit {

PrintTest_Hoske::PrintTest_Hoske() {
}

PrintTest_Hoske::~PrintTest_Hoske() {
}


/* Convert a graph to tgf format. */
static void test_convert_hoske(std::string name) {
	METISGraphReader reader;
    Graph G = reader.read("input/" + name + ".graph");
    std::ofstream out("input/" + name + "_hoske.tgf");
    print_hoske(out, G);
}

TEST_F(PrintTest_Hoske, testConvert_Hoske) {
	test_convert_hoske("add20");
	test_convert_hoske("data");
	test_convert_hoske("3elt");
	test_convert_hoske("uk");
	test_convert_hoske("add32");
}

}

#endif /*NOGTEST */
