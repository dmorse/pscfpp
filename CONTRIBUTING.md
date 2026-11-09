
# Contributing to PSCF 

PSCF has thus far been developed primarily at the University of Minnesota 
by David Morse and student co-authors, but we welcome contributions from
elsewhere.  It is currently maintained by David Morse (morse012@umn.edu) in 
repository github.com/dmorse/pscfpp.  Users are encouraged to contribute to 
PSCF by submitting github issues to report bugs, by submitting github pull 
requests, and by participating in discussions.

PSCF is free, open source software. It is distributed under the terms 
of the GNU General Public License (GPL) as published by the Free Software 
Foundation, either version 3 of the License or (at your option) any 
later version. PSCF is distributed without any warranty, without even 
the implied warranty of merchantability or fitness for a particular 
purpose.  See the LICENSE file in the root directory of this repository 
or the <a href=http://www.gnu.org/licenses/>gnu web page</a> for details.

## Issues, Pull Requests and Discussions:

The github interfaces for issues, pull requests and discussions provide 
complementary methods of discussing bugs, problems, and proposed 
improvements:

   - Submit github Issues to report bugs. 

   - Submit github pull requests for proposed contributions to the source 
     code and/or documentation. 

   - Use the github Discussions page to ask questions about usage and to 
     carry on conversations about limitations of the sofware and possible 
     improvements.

Users may also directly email David Morse (morse012@umn.edu) to discuss 
usage or proposed improvements. 

## Issues:

Any user may submit an Issue via github to report an error:

  - Issues should be used to report both bugs in the source code and 
    errors in any documentation.

  - If you are unsure whether a problem that you encounter when using
    PSCF is the result of a bug, consider first sending an email to 
    morse012@umn.edu or using the github Discussions page to discuss 
    the problem. 

  - The maintainer(s) of PSCF often use github Issues to document plans 
    for new features and other proposed changes. We ask, however, that 
    all others instead use email and/or github discussions to discuss 
    proposed changes that are not bug or error reports.

## Pull Requests:

To submit a pull request, a user must have a github account and must
have created a private version of pscfpp in their github account by 
forking the main PSCF repository, dmorse/pscfpp.  Pull requests may 
then be submitted via the github web interface. 

  - The first pull request from a new user should be the one in which
    the user accepts the contributor licensing agreement. See 
    discussion of this below.

  - All pull requests should be directed to the "feature-pull" branch 
    of the pscfpp repository, and not to the "devel" or "main" branch. 
    To do this, use the menu provided in the github web interface for 
    submitting a pull request that allowed you choose the branch to
    which the pull request should be applied.  If a pull request is 
    accepted, the maintainer will first merge it into the "feature-pull" 
    branch, and then do some testing in that branch before merging it 
    into the devel branch.  The "feature-pull" branch is generally 
    kept synchronized with the devel branch. 

  - Before submitting a pull request that introduces a new feature, 
    we recommend that you discuss the proposal with the maintainer,
    either by sending an email to morse012@umn.edu, or by starting a 
    discussion in the github Discussions space. More generally, we 
    recommend prior discussion of most nontrivial design decisions.

  - Before submitting a pull request, merge the "devel" branch of the
    main repository into the branch from which you will be submitting 
    the pull request and resolve any conflicts or problems that arise.

  - Contributors are expected to thoroughly test all code before 
    submitting a pull request. At a minimum, you must run all relevant 
    PSCF unit tests. If you modify any code that could be used by the 
    GPU-enabled version of pscf_rp, you must run all unit tests that 
    are compiled when both C++ and CUDA variants of that program have 
    been enabled (e.g., by invoking ./setopts -c1 -u1). Running tests
    that use CUDA requires access to an NVIDIA CPU.  It is necessary 
    to run all GPU-enabled tests if you modify any code in the 
    src/pscf or src/rp directories, since code in these directories 
    is accessible to the GPU-enabled version of pscf_rp. 

  - Pull requests that add new features should include new unit tests 
    that validate the correctness of any new features.

  - A commit that fixes an issue should include the text 
    "Fix #(number)" in the brief commit message, where (number) 
    denotes the number assigned to the relevant issue, which should
    be preceded by a bash (#) symbol. For example, the brief 
    message for a commit that fixes issue 268 should contain the 
    text "Fix #268".

  - Each pull request should be documented by a comment that 
    summarizes the purpose of the pull request and the nature 
    of associated changes.

## Contributor Licensing Agreement (CLA):

Contributors to PSCF must agree to a contributor licensing agreement (CLA) 
that is located in the file docs/cla/agreement_v1.txt within the PSCF
repository. Pull requests will only be accepted from contributors who 
have accepted this agreement. Acceptance of the CLA is not required in 
order to submit issues or to participate in discussions.

The PSCF agreement allows you to retain the copyright for any code or 
documentation that you contribute to the project, but simply gives the 
maintainers of PSCF permission to redistribute it as part of material 
that is distributed under the GNU GPL copyleft license or another open 
source license. Please see the file docs/cla/agreement_v1.txt for full 
details.  

To accept the PSCF CLA, please send a pull request to feature-pull that 
contains a single commit in which you should:

  - Add your own name to the list of contributors given in the file 
    docs/cla/contributors.txt, followed by the date. This is a list of 
    contributors who have executed the CLA, in which the date listed 
    after each contributors name is the date on which the agreement 
    was executed by that contributor. Do not change anything else in 
    this file. 

  - Give this commit a brief comment that says "Accepting CLA"

  - Add a longer comment to that commit that says "I accept the terms of 
    the contributor licensing agreement given in the file 
    docs/cla/agreement_v1.txt." 

The git commit and pull request that are used to accept the CLA 
should not contain changes to any files other than the file 
docs/cla/contributors.txt.

Each contributor only needs to go through this process once to execute
the agreement. The agreement is written so as to apply to all past and
future contributions to PSCF by an individual contributor.

