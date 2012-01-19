#include "Reporter.h"
namespace hemelb
{
  namespace reporting
  {
    Reporter::Reporter(const std::string &apath, const std::string &inputFile) :
        path(apath), snapshotCount(0), imageCount(0), stability(true), dictionary("Reporting dictionary")
    {
      dictionary.SetValue("CONFIG", inputFile);
    }

    void Reporter::Image()
    {
      imageCount++;
    }

    void Reporter::Snapshot()
    {
      snapshotCount++;
    }

    void Reporter::AddReportable(Reportable* reportable)
    {
      reportableObjects.push_back(reportable);
    }

    void Reporter::Write(const std::string &ctemplate, const std::string &as)
    {
      std::string output;
      ctemplate::ExpandTemplate(ctemplate, ctemplate::STRIP_BLANK_LINES, &dictionary, &output);
      std::string to = path + "/" + as;
      std::fstream file(to.c_str(), std::ios_base::out);
      file << output << std::flush;
      file.close();
    }
    ;

    void Reporter::FillDictionary()
    {
      dictionary.SetIntValue("IMAGES", imageCount);
      dictionary.SetIntValue("SNAPSHOTS", snapshotCount);

      if (!stability)
      {
        dictionary.AddSectionDictionary("UNSTABLE");
      }

      for (std::vector<Reportable*>::iterator reporters = reportableObjects.begin();
          reporters != reportableObjects.end(); reporters++)
      {
        (*reporters)->Report(dictionary);
      }
    }
  }
}
