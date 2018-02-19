#include "bandit/bandit/bandit.h"
//BANDIT NEEDS TO BE INCLUDED FIRST!!!

#include "../../src/cdoDebugOutput.h"
#include "../../src/modules.h"
#include "../../src/operator_help.h"
#include "../../src/process_int.h"
#include <iostream>

void *Oper1(void *test) { return test; }
void *Oper2(void *test) {return test; }

go_bandit([]() {
  bandit::describe("Process creation", []() {
    // oper1: child oper 2
    // oper2: child oper3 and oper4:
    // oper3: infile: 1
    // oper4: infile 2
    // outfile: filled by oper1
    std::vector<std::vector<const char *>> test_argvs = {{"-oper1","in1","in2"},{"-oper1","in1"}, {"-oper2","in1","in2", "ofile1"},{"-oper2","in1","ofile2"}};

    const int numberOfRuns = 4;
    int expectedInputForRun[numberOfRuns] = {2,1,2,1};
    int expectedOutputs[numberOfRuns] = {0,0,1,1};

    /*clang-format off*/
    //          Name     Func  Help   oper    mod    in
    //                                           num    out
    add_module("Oper1", {Oper1, {}, {"oper1"}, 1, 0, -1, 0});
    add_module("Oper2", {Oper2, {}, {"oper2"}, 1, 0, -1, 1});
    /*clang-format on*/


    unsigned int i;
    for (i = 0; i < numberOfRuns; i++) {
        createProcesses(test_argvs[i].size(),&test_argvs[i][0]);
        auto process = getProcess(0);
        std::string runInfo = std::string(getProcess(0)->operatorName) + " in run " + std::to_string(i + 1);
      bandit::it(
              "created inputs for:"+ runInfo, [&]() {
        AssertThat(
            process->childProcesses.size() + process->inputStreams.size(),
            snowhouse::Equals(expectedInputForRun[i]));
      });
      bandit::it("created outputs for: " + runInfo, [&]() {
        AssertThat(
            process->parentProcesses.size() + process->outputStreams.size(),
            snowhouse::Equals(expectedOutputs[i]));
      });
     bandit::it("created right amount of processes",[&]() {
        AssertThat(1, snowhouse::Equals(processNums()));
        });
      clearProcesses();
    }
  });
});
int main(int argc, char **argv) {
  CdoDebug::outfile = "variable_inputs_operators.debug";
  CdoDebug::print_to_seperate_file = true;

  CdoDebug::CdoStartMessage();
  CdoDebug::PROCESS = 1;
  CdoDebug::PSTREAM = 1;
  int result = bandit::run(argc, argv);
  CdoDebug::CdoEndMessage();

  return result;
}
