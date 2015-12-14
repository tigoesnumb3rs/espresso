#include "TclCommand.hpp"
#include <iostream>
#include <sstream>

#ifdef HAVE_CXX11
#include <type_traits>
#endif

using namespace std;

static int TclCommand_wrapper (ClientData data, Tcl_Interp *interp, int argc, char *argv[]);
  
void TclCommand::create_command(const string &command) {
  cout << "TclCommand::create_command(" << command << ")" << endl;
  cout << this << endl;
  Tcl_CreateCommand(interp, command.c_str(), reinterpret_cast<Tcl_CmdProc *>(TclCommand_wrapper), reinterpret_cast<ClientData>(this), NULL);    
}

static int TclCommand_wrapper (ClientData data, Tcl_Interp *interp, int argc, char *argv[]) {
  std::cout << "TclCommand_wrapper()" << std::endl;
  TclCommand *p = reinterpret_cast<TclCommand *>(data);
  cout << "p = " << p << endl;
  
  if(argc > 1) {
    list<string> args(argv + 1, argv + argc);

    try {
      p->parse_from_string(args);
      if(!args.empty()) {
        throw std::string("Unknown argument '").append(args.front()).append("'");
      }
    } catch(std::string &err) {
      Tcl_AppendResult(interp, err.c_str(), 0);
      return TCL_ERROR;
    }
  } else {
    Tcl_AppendResult(interp, p->print_to_string().c_str(), 0);
  }
  return TCL_OK;
}