#include "bandit/bandit/bandit.h"
//BANDIT NEEDS TO BE INCLUDED FIRST!!!

#include "../../src/cdoDebugOutput.h"
#include "../../src/modules.h"
#include "../../src/operator_help.h"
#include "../../src/process.h"
#include <iostream>
void *Oper1(void *test) { return test; }
void *Oper2(void *test) { return test; }
void *Oper3(void *test) { return test; }
void *Oper4(void *test) { return test; }

go_bandit([]() {
  bandit::describe("Process creation", []() {
    // oper1: child oper 2
    // oper2: child oper3 and oper4:
    // oper3: infile: 1
    // oper4: infile 2
    // outfile: filled by oper1
    std::vector<const char *> test_argv{"-oper1",      "-oper2", "-oper3",
                                        "input_file1", "-oper4", "input_file2",
                                        "output_file"};

    std::vector<int> expectedInputs{1,2,1,1};
    std::vector<int> expectedOutputs{1,1,1,1};


    /*clang-format off*/
    //          Name     Func  Help   oper    mod    in
    //                                           num    out
    add_module("Oper1", {Oper1, {}, {"oper1"}, 1, 0, 1, 1});
    add_module("Oper2", {Oper2, {}, {"oper2"}, 1, 0, 2, 1});
    add_module("Oper3", {Oper3, {}, {"oper3"}, 1, 0, 1, 1});
    add_module("Oper4", {Oper4, {}, {"oper4"}, 1, 0, 1, 1});
    /*clang-format on*/

    createProcesses(test_argv.size(), &test_argv[0]);

    for (unsigned int i = 0; i < Process.size(); i++) {
        auto process = Process.at(i);
      bandit::it("created inputs for:"+ std::string(Process.at(i).operatorName), [&]() {
        AssertThat(
            process.childProcesses.size() + process.inputStreams.size(),
            snowhouse::Equals(expectedInputs[i]));
      });
      bandit::it("created outputs for: " + std::string(Process.at(i).operatorName) , [&]() {
        AssertThat(
            process.parentProcesses.size() + process.outputStreams.size(),
            snowhouse::Equals(expectedOutputs[i]));
      });

    }
    CdoDebug::CdoEndMessage();
  });
});
int main(int argc, char **argv) {
  CdoDebug::outfile = "chained_debug.txt";
  CdoDebug::print_to_seperate_file = true;
  CdoDebug::CdoStartMessage();
  return bandit::run(argc, argv);
}
