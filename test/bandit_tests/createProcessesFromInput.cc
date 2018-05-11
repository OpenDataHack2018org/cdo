#include "bandit/bandit/bandit.h"
// BANDIT NEEDS TO BE INCLUDED FIRST!!!

#include "../../src/cdoDebugOutput.h"
#include "../../src/modules.h"
#include "../../src/operator_help.h"
#include "../../src/process_int.h"
#include <iostream>
/* clang-format off */
void * in0_out1(void *test) { return test; }
void * in1_out1(void *test) { return test; }
void * in2_out1(void *test) { return test; }
void * inVariable_out1(void *test) { return test; }

using namespace snowhouse;

go_bandit([]() {
//==============================================================================
  add_module("in0_out1", { in0_out1, {}, { "in0_out1" }, 1, 0, 0, 1 });
  add_module("in1_out1", { in1_out1, {}, { "in1_out1" }, 1, 0, 1, 1 });
  add_module("in2_out1", { in2_out1, {}, { "in2_out1" }, 1, 0, 2, 1 });
  add_module("inVariable_out1", { inVariable_out1, {}, { "inVariable_out1" }, 1, 0, -1, 0 });

  /* clang-format on */

  ParseStatus parseStatus;
  ProcessStatus processStatus;

  int result_parse;
  int expected_parse;
  int result_process;
  int expected_process;

  bandit::before_each([&]() { clearProcesses(); });

  //-----------------------------Test_01------------------------------------------
  //------------------------------------------------------------------------------
  bandit::describe("Negative test for unprocessed inputs", [&]() {
    /* clang-format off */
    std::vector<const char *> argv_unprocessedInput{
        "-in2_out1", "-in0_out1", "-in0_out1", "-in0_out1", "out" };
    /* clang-format on */
    parseStatus = createProcessesFromInput(argv_unprocessedInput.size(),
                                           &argv_unprocessedInput[0]);

    result_parse = static_cast<int>(parseStatus);
    expected_parse = static_cast<int>(ParseStatus::UnprocessedInput);

    bandit::it("has unprocessed Input", [&]() {
      AssertThat(result_parse, snowhouse::Equals(expected_parse));
    });
  });

  //-----------------------------Test_02------------------------------------------
  //------------------------------------------------------------------------------
  bandit::describe("Negative test for miss placement of brackets", [&]() {
    /* clang-format off */
    std::vector<const char *> argv_missingCloseBracket{
        "-in2_out1", "[", "-in0_out1", "-in0_out1", "out"
        };
    /* clang-format on */
    parseStatus = createProcessesFromInput(argv_missingCloseBracket.size(),
                                           &argv_missingCloseBracket[0]);

    result_parse = static_cast<int>(parseStatus);
    expected_parse = static_cast<int>(ParseStatus::ClosingBracketMissing);
    bandit::it("it misses a ']'", [&]() {
      AssertThat(result_parse, snowhouse::Equals(expected_parse));
    });
  });

  //-----------------------------Test_03------------------------------------------
  //------------------------------------------------------------------------------
  bandit::describe("Negative test for miss placement of brackets", [&]() {
    /* clang-format off */
    std::vector<const char *> argv_missingOpenBracket{
        "-in2_out1", "-in0_out1", "-in0_out1", "]", "out"
        };
    /* clang-format on */
    parseStatus = createProcessesFromInput(argv_missingOpenBracket.size(),
                                           &argv_missingOpenBracket[0]);

    result_parse = static_cast<int>(parseStatus);
    expected_parse = static_cast<int>(ParseStatus::OpenBracketMissing);
    bandit::it("it misses a '['", [&]() {
      AssertThat(result_parse, snowhouse::Equals(expected_parse));
    });
  });

  //-----------------------------Test_04------------------------------------------
  //------------------------------------------------------------------------------
  bandit::describe("Negative test for miss placement of brackets", [&]() {
    /* clang-format off */
    std::vector<const char *> argv_wrongBracketTooMany{
        "-in2_out1", "[", "-in0_out1", "-in0_out1", "-in0_out1", "]", "out"
    };
    /* clang-format on */
    parseStatus = createProcessesFromInput(argv_wrongBracketTooMany.size(),
                                           &argv_wrongBracketTooMany[0]);
    processStatus = getProcess(0)->checkStreamCnt();

    result_parse = static_cast<int>(parseStatus);
    expected_parse = static_cast<int>(ParseStatus::Ok);
    result_process = static_cast<int>(processStatus);
    expected_process = static_cast<int>(ProcessStatus::TooManyStreams);

    bandit::it("this one has too many inputs", [&]() {
      AssertThat(result_parse, Is().EqualTo(expected_parse));
      AssertThat(result_process, Is().EqualTo(expected_process));
    });
  });

  //-----------------------------Test_05------------------------------------------
  //------------------------------------------------------------------------------
  bandit::describe("Negative test for miss placement of brackets", [&]() {
    /* clang-format off */
    std::vector<const char *> argv_wrongBracketTooFew{ 
        "-in2_out1", "[", "-in0_out1", "]", "out"
    };
    /* clang-format on */
    parseStatus = createProcessesFromInput(argv_wrongBracketTooFew.size(),
                                           &argv_wrongBracketTooFew[0]);
    processStatus = getProcess(0)->checkStreamCnt();

    bandit::it("this one has too few inputs", [&]() {
      AssertThat(static_cast<int>(parseStatus),
                 Is().EqualTo(static_cast<int>(ParseStatus::Ok)));
      AssertThat(static_cast<int>(processStatus),
                 Is().EqualTo(static_cast<int>(ProcessStatus::TooFewStreams)));
    });
  });
});

//==============================================================================
int
main(int argc, char **argv)
{
  CdoDebug::outfile = "createProcessFromInput.debug";
  CdoDebug::print_to_seperate_file = true;

  CdoDebug::CdoStartMessage();
  CdoDebug::PROCESS = 1;
  CdoDebug::PSTREAM = 1;
  int result = bandit::run(argc, argv);
  CdoDebug::CdoEndMessage();

  return result;
}
