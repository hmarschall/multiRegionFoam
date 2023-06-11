# multiRegionFoam

Developed by the Computational Multiphase Flow research group.

* [Research Group Website](https://www.mathematik.tu-darmstadt.de/cmf/)
* [Report Bug](https://bitbucket.org/hmarschall/multiregionfoam/issues?status=new&status=open)
* [Request Feature](https://bitbucket.org/hmarschall/multiregionfoam/issues?status=new&status=open)

## About The Project

Unified framework for solving multiphysics problems of the multi-region coupling 
type within OpenFOAM (FOAM-extend). This framework is intended to supersede the
existing solver with the same name. The design of the new framework is modular, 
allowing users to assemble a multiphysics problem region-by-region and coupling 
conditions interface-by-interface. The present approach allows users to choose 
between deploying either monolithic or partitioned interface coupling for each 
individual transport equation. The formulation of boundary conditions is 
generalised in the sense that their implementation is based on the mathematical 
jump/transmission conditions in the most general form for tensors of any rank.

---


## Getting Started

### Prerequisites

* Working installation of foam-extend 4.1:
For installing foam-extend 4.1 please refer to and follow the [installation instructions](https://openfoamwiki.net/index.php/Installation/Linux/foam-extend-4.1) step-by-step!
    * [Installation/Linux/foam-extend-4.1/Ubuntu](https://openfoamwiki.net/index.php/Installation/Linux/foam-extend-4.1/Ubuntu)
    * [Installation/Linux/foam-extend-4.1/CentOS](https://openfoamwiki.net/index.php/Installation/Linux/foam-extend-4.1/CentOS)

### Installation

With a sourced working foam-extend 4.1 installation, clone the repository and 
build the library and solvers:

```bash
./Allwmake
```


## Usage

The `tutorials` directory has tutorials to showcase the library's functionality. 
Please take a look there for examples of usage. There are `Allrun` scripts 
demonstrating the proper usage.


## Contributing

We invite everyone in the FOAM community to collaborate with us and jointly 
develop multiRegionFoam. For this, start by forking the repository to your own 
account. Then, clone the forked repository to your local machine and create a 
new branch for your changes. Make the necessary modifications, commit your 
changes, and push the branch to your forked repository. Finally, open a pull 
request from your branch to the original repository, and providing a clear 
description of your changes. Collaborate with reviewers, address feedback, 
and once approved, your contributions can be merged.

The `dev` branch is the corner stone of the development, please branch all of 
your feature/bugFix branches off of it, And "rebase" your branches on it before 
issuing a pull request. When your branch gets merged, it's considered a 
"best-practice" to delete your feature branch and start a fresh one.


## License

Released under the GNU Public License - see code headers for details.


## Contact

[Email @hmarschall](mailto:holger.marschall@tu-darmstadt.de)


## Acknowledgements

The devopment of this project has been partly funded by

* Hessian Ministry of Higher Education, Research, Science and the Arts,
* National High Performance Computing Center for Computational Engineering Science (NHR4CES)

