#include "bandit/bandit/bandit.h"
// BANDIT NEEDS TO BE INCLUDED FIRST!!!

#include "../../src/cdoDebugOutput.h"
#include "../../src/modules.h"
#include "../../src/operator_help.h"
#include "../../src/process_int.h"
#include <iostream>
void *in1_out1(void *test) { return test; }
void *in2_out1(void *test) { return test; }
void *inVariable_out1(void *test) { return test; }

go_bandit([]() {
  add_module("in1_out1", {in1_out1, {}, {"in1_out1"}, 1, 0, 1, 1});
  add_module("in2_out1", {in2_out1, {}, {"in2_out1"}, 1, 0, 2, 1});
  add_module("inVariable_out1",
             {inVariable_out1, {}, {"inVariable_out1"}, 1, 0, -1, 0});

  // this test checks if operators with non variable input numbers can be
  // chained
  bandit::describe("Process creation for non variable operators", []() {
    // in1_out1: child oper 2
    // in2_out1: child in1_out1 and in1_out1:
    // in1_out1: infile: 1
    // in1_out1: infile 2
    // outfile: filled by in1_out1
    std::vector<const char *> test_argv{
        "-in1_out1", "-in2_out1",   "-in1_out1",  "input_file1",
        "-in1_out1", "input_file2", "output_file"};

    std::vector<unsigned int> expectedInputs{1, 2, 1, 1};
    std::vector<unsigned int> expectedOutputs{1, 1, 1, 1};

    /*clang-format off*/
    //          Name     Func  Help   oper    mod    in
    //                                           num    out
    /*clang-format on*/

    createProcessesFromInput(test_argv.size(), &test_argv[0]);

    int i;
    for (i = 0; i < processNums(); i++) {
      auto process = getProcess(i);
      bandit::it("created inputs for:" +
                     std::string(getProcess(i)->operatorName),
                 [&]() {
                   AssertThat(process->inputStreams.size(),
                              snowhouse::Equals(expectedInputs[i]));
                 });
      bandit::it("created outputs for: " +
                     std::string(getProcess(i)->operatorName),
                 [&]() {
                   AssertThat(process->outputStreams.size(),
                              snowhouse::Equals(expectedOutputs[i]));
                 });
    }
    bandit::it("created right amount of processes",
               [&]() { AssertThat(i, snowhouse::Equals(processNums())); });
  });

  clearProcesses();

  // this test checks if multiple operators can be chained if the first operator
  // is part of a module with variable number of input streams
  bandit::describe(
      "Process creation containing operators with variable number of input "
      "streams",
      []() {
        // in1_out1: child oper 2
        // in2_out1: child in1_out1 and in1_out1:
        // in1_out1: infile: 1
        // in1_out1: infile 2
        // outfile: filled by in1_out1
        std::vector<const char *> test_argv{
            "-inVariable_out1", "-in2_out1",   "-in1_out1", "input_file1",
            "-in1_out1",        "input_file2", "-in1_out1", "input_file3",
            "-in1_out1",        "input_file4"};

        std::vector<unsigned int> expectedInputs{3, 2, 1, 1, 1, 1};
        std::vector<unsigned int> expectedOutputs{0, 1, 1, 1, 1, 1};

        createProcessesFromInput(test_argv.size(), &test_argv[0]);

        int i;
        for (i = 0; i < processNums(); i++) {
          auto process = getProcess(i);
          std::string runInfo = std::string(getProcess(i)->operatorName) +
                                " in run: " + std::to_string(i + 1);
          bandit::it("created inputs for:" + runInfo, [&]() {
            AssertThat(process->inputStreams.size(),
                       snowhouse::Equals(expectedInputs[i]));
          });
          bandit::it("created outputs for: " + runInfo, [&]() {
            AssertThat(process->outputStreams.size(),
                       snowhouse::Equals(expectedOutputs[i]));
          });
        }
        bandit::it("created right amount of processes",
                   [&]() { AssertThat(i, snowhouse::Equals(processNums())); });
      });

});
int main(int argc, char **argv) {
  CdoDebug::outfile = "chainedOperators.debug";
  CdoDebug::print_to_seperate_file = true;

  CdoDebug::CdoStartMessage();
  CdoDebug::PROCESS = 1;
  CdoDebug::PSTREAM = 1;
  int result = bandit::run(argc, argv);
  CdoDebug::CdoEndMessage();

  return result;
}
