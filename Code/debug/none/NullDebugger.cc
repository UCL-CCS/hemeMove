#include <iostream>
#include "debug/none/NullDebugger.h"

namespace hemelb
{
  namespace debug
  {
    // Do nothing!
    void NullDebugger::Attach(void) {}
    void NullDebugger::BreakHere(void){}
    
    NullDebugger::NullDebugger(char* executable) : Debugger(executable) {}
    
    Debugger* PlatformDebuggerFactory(char *executable) {
      return new NullDebugger(executable);
    }
  }
}
