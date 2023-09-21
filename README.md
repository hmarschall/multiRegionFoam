# multiRegionFoam

[![DOI](https://zenodo.org/badge/652142737.svg)](https://zenodo.org/badge/latestdoi/652142737)

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
develop multiRegionFoam. 
For this, the follwoing steps need to be carried out

1. Fork this repository to your own account.
1. Clone the forked repository to your local machine:  
`git clone <URL_of_your_fork>`
1. Also add the original multiRegionFoam as a remote to be able to regularly 
pull updates from the original repository:  
`git remote add OriginalMultiRegionFoam git@bitbucket.org:hmarschall/multiregionfoam.git`
    * To pull changes from the `dev` branch of the original multiRegionFoam repo 
    run:  
    `git pull OriginalMultiRegionFoam dev`
1. Create a new branch (The `dev` branch is the corner stone of the development,
please branch all of your feature/bugFix branches off of it):  
`git checkout -b ＜name_of_your_new_branch＞`
1. Push your new branch to your forked remote repo:  
`git push -u origin ＜name_of_your_new_branch＞`
1. Make the necessary modifications and commit them to your branch while 
providing descriptive commit messages (see this 
[Link](https://www.atlassian.com/git/tutorials/saving-changes) 
on how to work with git).
    * Try to keep your branch up to date with the new developments and bug fixes
    in the original multiRegionFoam repo. For this, do the following
        1. Switch from your branch back to the `dev` branch:  
        `git checkout dev`
        1. Pull the updates from the original repository:  
        `git pull OriginalMultiRegionFoam dev`
        1. Switch back to your branch:  
        `git checkout ＜name_of_your_new_branch＞`
        1. Merge the updated `dev` branch into your branch:  
        `git merge dev`
1. Push your changes regularly to your forked repo on bitbuket:  
`git push origin ＜name_of_your_new_branch＞`
1. If you are satisfied with your new developments and all your changes are 
pushed to your remote repository yoou can finally create a pull request from 
your branch in your remote repository on bitbucket to the original 
multiRegionFoam reposity, while providing a clear description of your changes 
(see this [Link](https://support.atlassian.com/bitbucket-cloud/docs/create-a-pull-request/) 
on how to create a pull request on bibucket). Collaborate with reviewers, address
feedback, and once approved, your contributions can be merged.
    * Make sure your branch is up to date with the most recent changes in the 
    original multiRegionFoam repo before creating the pull request 
    * When your branch gets merged, it's considered a "best-practice" to delete 
    your feature branch and start a fresh one.


## License

Released under the GNU Public License - see code headers for details.


## Contact

[Email @hmarschall](mailto:holger.marschall@tu-darmstadt.de)


## Acknowledgements

The devopment of this project has been partly funded by

* Hessian Ministry of Higher Education, Research, Science and the Arts,
* National High Performance Computing Center for Computational Engineering Science (NHR4CES)

