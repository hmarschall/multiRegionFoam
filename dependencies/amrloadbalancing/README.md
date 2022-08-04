A library for load-balanced adaptive mesh refinement based on Foam-Extend.

## Installation

With a sourced working foam-extend 4.0 installation, clone the repository and build the library and solvers:

```bash
./Allwmake
```

To use the library, don't forget to add the following to your
`system/controlDict`:

```cpp
libs( "libbalancedRefinementFvMesh.so" );
```

## Features

The `run` directory has tutorials to showcase the library's functionality.
Please take a look there for examples of usage.

### Adaptive polyhedral mesh refinement

The `dynamicPolyMultiRefinementFvMesh` offers multi-criteria refinement algorithm
selection, collecting them from a `dynamicMeshDict.dynamicRefineFvMeshCoeffs.refinements` list.

In addition to the default refinement selection algorithms in Foam-Extend, the library provides:

- `constrainedFieldBoundsRefinement`
- `discretisationErrorRefinement`
- `fieldBoundsNoBoundariesRefinement`
- `fieldCurlRefinement`
- `fieldGradientRefinement`
- `interfaceRefinement`

Each refinement selection criterion has its own custom parameters,
and shares default ones with the others:

> PS: the refinement engine carries out these refinement operations
> in the order they appear in the dictionary.

```cpp
dynamicRefineFvMeshCoeffs
{
    // Global refinement configuration
    // Maximum refinement level
    maxRefinementLevel   3;
    ...

    refinements
    (
        unrefineBasedOnAlpha1
        {
            // Refinement-specific configuration, has priority
            // Maximum refinement level
            maxRefinementLevel   2;
	        ...
            // Refinement selection criteria
            refinementSelection
            {
                // Refines around field interface
                type        interfaceRefinement;
		        ...
            }
        }
   	    refineBasedOnSomethingElse
        {
            // maxRefinementLevel  is considered 3 if not overwritten here
	        ...
            // Refinement selection criteria
            refinementSelection
            {
                // Refines around field interface
                type        interfaceRefinement;
		        ...
            }
        }
    );
}
```

#### User-supplied code for refinement selection

If the available refinement selection algorithms are not suitable
for your needs; there is also a `codedFieldBoundsRefinement` which
enables the user to manipulate a (scalar) field and a set of bounds
in order to select refinement/unrefinement candidates (refining if
field value is between lower and upper bound).

To use the `codedFieldBoundsRefinement`:

1. Add the following to your `controlDict`:
   ```cpp
    InfoSwitches
    {
            allowSystemOperations 1;
    }
   ```
2. Point `FOAM_CODE_TEMPLATE` environment variable to `etc/codeTemplates`
   from this repository.

3. Add a refinement selection with a `code` entry:
   ```cpp
   basedOnAlpha
   {
       refineInterval   1;
       unrefineInterval 1;
       // ...

       refinementSelection
       {
           type        codedFieldBoundsRefinement;
           fieldName   target; // Name of the scalar field to create
           // Initial bounding values
           lowerBound  0;
           upperBound  1;

           cellPointCellSmoothing on;

           // Include your own headers
           codeInclude
           #{
               #include "additionalHeader.H"
           #};

           // Link to your own libraries
           codeLibs
           #{
               -lsomeLibrary
           #};

           // Compiler arguments for the dynamic code
           codeOptions
           #{
               -I$(LIB_SRC)/someLibrary/lnInclude
           #};

           // Refine if alpha1 > 80% of MaxAlpha and if x > 0.5m, unrefine otherwise
           // The following will be in its own scope, so everything
           // you define here will be lost as soon as we leave this scope
           code
           #{
               Info << "** dynamicCode **" << endl;
               auto alpha1 = mesh().lookupObject<volScalarField>("alpha1");
               forAll(field_, ci)
               {
                   field_[ci] = mesh.C()[ci].x() > 0.5 ? alpha1[ci] : 0.0 ;
               }
               lowerBound_ = 0.8*gMax(alpha1);
               // upperBound_ of 1 is good
               // Also dict() returns a copy of the basedOnAlpha dictionary
           #};

           localCode
           #{
                // Optional code section, executed just before selecting
                // unrefinement candidates
           #};
       }
   }
   ```
   Generated code for the above snippet can be found in
   `/path/to/case/dynamicCode/codedRefinement_basedOnAlpha`
   after the first call to `updateMesh` (or the solver)

> Note that each `codedFieldBoundsRefinement` will be compiled only once
> so it produces no overhead other than an extra virtual call; and you
> can use more than one `codedFieldBoundsRefinement` in multi-criteria
> refinement setups.

### Adaptive polyhedral mesh refinement with load balancing

Use `loadBalanceDynamicPolyRefinementFvMesh` as the `dynamicFvMesh`
for load-balanced adaptive polyhedral mesh refinement.

The load balancing happens after all refinement operations and its 
configuration is supplied in a separate
`dynamicMeshDict.loadBalanceFvMeshCoeffs`.

## Contributing

The `dev` branch is the corner stone of the development, please branch all of your feature/bugFix branches off of it, And "rebase" your branches on it before issuing a pull request. When your branch gets merged, it's considered a "best-practice" to delete your feature branch and start a fresh one.
