# Testing framework for LoadBalanced-AMR

## Design principles

- Tests are created and run with [Catch2 v2.x](https://github.com/catchorg/Catch2/tree/v2.x)
  ("Header-only" library)
    - Catch2 v3 is a better choice but requires C++14 and no longer "Header-only"
    - [Catch2 docs](https://github.com/catchorg/Catch2/tree/v2.x/docs)
      and [Phil's talk at CppCon 2018](https://github.com/catchorg/Catch2/tree/v2.x/docs)
      are the best places to learn more
- [FakeIt](https://github.com/eranpeer/FakeIt) (also Header-only) is used for mocking.
- Each library has its own test driver
    - Convention: The name of the test binary should end with `TestDriver`
    - Each class's tests suite should be written in a separate `.C` file
- The test driver acts like an OpenFOAM solver and must run on an OpenFOAM case
  directory. For a fully-in-memory testing setup, it's enough to move the test case
  to `/dev/shm` on most systems
- Tests for the same class in serial and parallel runs belong in separate files.
- Supported reporters:
    - `stdout` for output at the console (Best for local and CI testing)
    - `compact` for status overview in a single line
    - `xml` and `junit` for Data collectors

## Why write (Unit) tests?

- To make sure new functionality works as expected
- To make sure further development does not introduce bugs, or change older
  behaviour.
- To better document the expected usage of code.
- To gauge performance if needed.

## Best practices and guidelines for tests

### General guidelines

- BDD-style tests (GIVEN, WHEN, THEN, AND_THEN ... etc) are preferred although
    you can still use simple test-cases with sections.
- Avoid relying on dictionaries stored on disk; It's recommended to write the
    dictionaries programmatically immediately before using them and never
    writing time directories to disk.
- Each test should be manually tagged for the following
  (Task left to the use for now):
    - `[Integration]` or `[Unit]` for integration and unit tests respectively
    - `[Parallel]` for tests which are meant to test parallel functionality,
    or `[Serial]` otherwise.
- OpenFOAM results should be in ASCII uncompressed format.
- Avoid using `runTimeSelectionTable` features if you plan to use Catch2
    `GENERATORS`
- The main focus is to test the interface of classes
    - But if there is a need to test a private method, use macros defined
      in `include/memberStealer.H`
- Use `CAPTURE()` to record useful case parameters.
- Standard output is consumed by Catch2, print to standard error to see
  the logs at the console.
- Some convenient macros can be found in `include/testMacros.H`

### Strongly recommended coding guidelines

> Break these guidelines and you'll have a hard time debugging the tests!

1. Create only **ONE** Time object at the start of each test,
   at global scope; enclosed in a unique namespace. For example:
   ```cpp
   namespace multiCriteriaSerial
   {
       // Requirement 0-0: Time for serial runs
       #include "createTestTime.H"
   
       // Requirement 0-1: Pointer for mesh
       Foam::autoPtr<Foam::dynamicPolyMultiRefinementFvMesh> meshPtr;
   
   } // End of namespace multiCriteriaSerial
   
   using namespace Foam;
   using namespace multiCriteriaSerial;

   SCENARIO(" ... ") { ... }
   TEST_CASE(" ... ") { ... }
   ```

2. As the previous example outlined, it's recommended to have a single
   pointer to a (**concrete**) mesh class, which we then reset over and
   over:
   ```cpp
   namespace multiCriteriaSerial
   {
       // Requirement 0-0: Time for serial runs
       #include "createTestTime.H"
   
       // Requirement 0-1: Pointer for mesh
       Foam::autoPtr<Foam::dynamicPolyMultiRefinementFvMesh> meshPtr;
   
   } // End of namespace multiCriteriaSerial
   
   using namespace Foam;
   using namespace multiCriteriaSerial;

   TEST_CASE(" ... ")
   {
        // Generate a variable, test case runs 10 times
        auto number = GENERATE( range(1,10) );

        // Reset the mesh (macro from testMacros.H)
        resetMeshPointer
        (
            runTime, // The time object
            meshPtr, // Pointer name
            dynamicPolyMultiRefinementFvMesh, // Mesh class
            dynamicFvMesh::defaultRegion // Mesh region
        );
        auto& mesh = meshPtr();

        // Do what's needed, and
        meshPtr->clear();
   }
   ```
> Global pointers are hated for a resounding reason;
> Never try to access these mesh pointers from other files please!

3. Parallel tests are tricky when you want to use advanced Catch2 features,
   especially `GENERATORS`; which want the time object to be in global scope,
   but MPI comms don't work for global variables (So, `Pstream::myProcNo()`
   will always evaluate to 0 in global scope).

The `memberStealer` template becomes very useful in overcoming this
problem. By overwriting `processorCase_` and `case_` in the time object,
we can point it to the correct processor directory inside the test case:

```cpp
namespace multiCriteriaParallel
{
    // Requirement 0-0: Time for parallel runs
    #include "createTestTime.H"

} // End of namespace multiCriteriaParallel

using namespace Foam;
using namespace multiCriteriaParallel;

// A macro to do the necessary specializations
prepareTimePaths();

TEST_CASE(" ... ", "[Parallel]")
{
     // Alter processorCase_ and case_ for parallel cases
     word newCase = "processor"+Foam::name(Pstream::myProcNo());
     modifyTimePaths(runTime, true, newCase);

     // Generate a variable
     auto number = GENERATE( range(1,10) );
     CAPTURE(runTime.caseName(), number);

     // ...

     // Wait for requests
     Pstream::waitRequests();
}
```
Note the `waitRequests` call at each generated case end.

## Catch2 Output

On failure, Catch2 displays a message containing `TEST_CASE` name,
or the BDD logic (SCENARIO, GIVEN, WHEN, THEN) with the captured variables.

> With parallel runs, this information is printed for each processor which
> has some tests running

If you pass `-s` to the test driver, the same information is printed
for successful tests too:
```
-------------------------------------------------------------------------------
Scenario: Multi-criteria adaptive poly refinement
      Given: A time object, and a dynamicMeshDict
       When: Two refinement criterions are applied,  with maxRefinementLevels
             (1,3) each based on its own field
       Then: Interface refines, then unrefine previous interface region
        And: Second refinement doesn't unrefine previouly refined cells
-------------------------------------------------------------------------------
dynamicPolyMultiRefinementFvMesh/dynamicPolyMultiRefinementFvMeshTest.C:297
...............................................................................

dynamicPolyMultiRefinementFvMesh/dynamicPolyMultiRefinementFvMeshTest.C:299: PASSED:
  REQUIRE( boxDidntUnrefine )
with expansion:
  true
with messages:
  maxRefLevel1 := 1
  maxRefLevel2 := 3
```

## Integration with CI and reporting

When running serial and parallel tests, `$CATCH_SERIAL_OPTIONS` and
`$CATCH-PARALLEL_OPTIONS` are passed respectively to the test driver.
To collect reports on serial tests, one could use:
```bash
export CATCH_SERIAL_OPTIONS="-r junit -o serialTests.xml"
```

Note that, in parallel tests, you are not supposed to include a `-o file`
option because `--outputfile parallelTests.xml` is already passed
to `mpirun`, but you can still choose the reporter.

> This produces parallelTests.xml.1.* files on the caller machine.
> If you want to run tests on multiple machines, you have to get the files
> to a center location for reporting.

### Jenkins configuration

0. Make sure Jenkins is running as a user which has access to Foam-Extend
1. Install Junit plugin for Jenkins
2. For the build step, execute something like this:
   ```bash
    #!/bin/bash
    source /path/to/foam/etc/bashrc
    cd /path/to/repo/amrloadbalancing # Or checkout ...
    ./Allwmake
    export CATCH_SERIAL_OPTIONS="-s -d yes -r junit -o serialTests.xml"
    export CATCH_PARALLEL_OPTIONS="-s -d yes -r junit"
    ./Alltest
   ```
3. To publish the Junit report, select `**/*Tests.xml, **/parallelTests.xml*` files
   as inputs.
