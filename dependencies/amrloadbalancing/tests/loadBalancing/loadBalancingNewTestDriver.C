#define CATCH_CONFIG_RUNNER
#include "catch.H"
#include "error.H"

#include <mpi.h>
#include "Pstream.H"

int argcG = 0;
char **argvG = nullptr;

bool inParallel = false;

using namespace Catch::clara;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    // OpenFOAM settings
    Foam::Warning.level = 0;
    Foam::FatalError.throwExceptions();

    // Grab a catch session
    Catch::Session session;

    // Hook to arguments parser to allow parallel with one dash
    bool parallel = false;
    auto cli = session.cli()
        | Opt(parallel)
          ["-p"]["--parallel"]
          ("Ignored by catch, used by Pstream");

    // Update CLI args
    session.cli(cli);

    // Parse the CLI args
    char *newargs[argc+1];
    if (parallel)
    {
        inParallel = true;
        char* flag = "-parallel";
        argcG = argc + 1;
        for (size_t i = 0; i < argc; i++) {
            newargs[i+1] = argv[i];
        }
        newargs[0] = flag;
        argvG = newargs;
    } else {
        argcG = argc;
        argvG = argv;
    }
    auto result = session.applyCommandLine(argc, argv);
    if (result != 0) return result;

    // Init MPI comms if requested
    if (parallel)
    {
        Foam::Pstream::init(argcG,argvG);
    }

    // Run tests and return error code
    auto code = session.run();

    // Finialize MPI comms;
    if (parallel)
    {
        Foam::Pstream::waitRequests();
        MPI_Finalize();
    }

    return code;
}


// ************************************************************************* //
