#include "bandit/bandit/bandit.h"

#include "../../src/modules.h"
#include "../../src/operator_help.h"
#include "../../src/process.h"
#include <iostream>


void *Test(void *ptr){}
std::vector<std::string> TestHelp = {"TEST", "HELP"};
std::vector<char *> test_argv{"-test", "in_file","out_file" };

go_bandit([]() {
  bandit::describe("Process: 1 input, 1 output", []() {

    add_module("Test", {Test, TestHelp, {"test"}, 1,0,1, 2});

    createProcesses(test_argv.size(), &test_argv[0]);

    process_t test_process = Process.find(0)->second;

    bandit::it("should have appropriate number of in streams", [&]() {
      AssertThat(test_process.getInStreamCnt(), snowhouse::Equals(1));
    });

    bandit::it("should have appropriate number of out streams", [&]() {
      AssertThat(test_process.getOutStreamCnt(), snowhouse::Equals(2));
    });

  });
});

int main(int argc, char **argv) { return bandit::run(argc, argv); }
