#include "bandit/bandit/bandit.h"

#include "../../src/modules.h"
#include "../../src/operator_help.h"
#include "../../src/process.h"
#include <iostream>


void *Info(void *test){return test;}

go_bandit([]() {
  bandit::describe("Process creation", []() {

    std::vector<const char *> test_argv{"-info", "some_test_bs"};
    add_module("Info", {Info, InfoHelp, {"info"}, 1, 0, -1, 0});
    process_t *test_process = processCreate(test_argv[0]);

    bandit::it("should have the name of the operator", [&]() {
      AssertThat(test_process->operatorName, snowhouse::Equals("info"));
    });

    bandit::it("should update the count for active processes", [&]() {
      AssertThat(processNumsActive(), snowhouse::Equals(1));
    });

    bandit::it("should update the cound for existing processes",
               [&]() { AssertThat(processNums(), snowhouse::Equals(1)); });

    bandit::it("new size of Process should be updated",
               [&]() { AssertThat(Process.size(), snowhouse::Equals(1ul)); });

    bandit::it("ID should be set right",
               [&]() { AssertThat(test_process->m_ID, snowhouse::Equals(0)); });

  });
});

int main(int argc, char **argv) { return bandit::run(argc, argv); }
