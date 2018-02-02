#include "bandit/bandit/bandit.h"

#include "../../src/cdoDebugOutput.h"
#include "../../src/util_wildcards.h"
#include <fstream>
#ifdef HAVE_CONFIG_H
#include "../../src/config.h"
#endif

//==========================================================================
go_bandit([]() {

  // File name rules
  const std::string prefix = "testFile";
  const std::string suffix = ".banditTestFile";
  std::vector<std::string> wildcards = {"*", "?"};
  std::vector<std::string> testArgv;
  std::string toBeExpanded;

  const int fileCount = 5;

  // +1 because there has to be another string at the binning of testArgv since
  // the first entry will be ignored by expand wildcards
  unsigned int expectedFileCnt = fileCount * wildcards.size() + 1;

  // storage for temp files;
  std::ofstream files[fileCount];

  //--------------------------------------------------------------------------
  bandit::before_each([&]() {
    std::string fileName;
    // Creating test files
    for (int fileID = 0; fileID < fileCount; fileID++) {
      fileName = prefix + std::to_string(fileID) + suffix;
      files[fileID].open(fileName);
      files[fileID].close();
    }
    // see comment at definition of expectedFileCnt
    testArgv.push_back("wildcards");
    // setting the actual filenames
    for (int wildCardIDX = 0; wildCardIDX < wildcards.size(); wildCardIDX++) {
      testArgv.push_back(prefix + wildcards[wildCardIDX] + suffix);
    }
  });

  //--------------------------------------------------------------------------
  bandit::describe("Expanding wildcard", [&]() {
    bandit::it("whatever", [&]() {
#ifdef HAVE_WORDEXP_H
      std::vector<std::string> expandedWildCards = expandWildCards(testArgv);
      for (auto i : expandedWildCards) {
        std::cout << i << " ";
      }
      AssertThat(expandedWildCards.size(), snowhouse::Equals(expectedFileCnt));
#else
    AssertThat(0, snowhouse::Equals(0));
#endif
    //--------------------------------------------------------------------------
    });
  });
});

int main(int argc, char **argv) {
  CdoDebug::outfile = "wildcards.debug";
  CdoDebug::print_to_seperate_file = true;

  CdoDebug::CdoStartMessage();
  CdoDebug::PROCESS = 1;
  CdoDebug::PSTREAM = 1;
  int result = bandit::run(argc, argv);
  CdoDebug::CdoEndMessage();

  return result;
}
