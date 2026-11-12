
# Contributing to PSCF 

PSCF has thus far been developed primarily at the University of Minnesota,
but we welcome contributions from elsewhere.  It is maintained by 
David Morse (morse012@umn.edu) in repository github.com/dmorse/pscfpp. 
Users are encouraged to contribute to PSCF by submitting github issues 
to report bugs, submitting github pull requests, or participating in 
the discussion forum.

PSCF is free, open source software. It is distributed under the terms of 
the GNU General Public License (GPL) as published by the Free Software 
Foundation, either version 3 of the License or (at your option) any 
later version. PSCF is distributed without any warranty, without even 
the implied warranty of merchantability or fitness for a particular 
purpose.  See the LICENSE file in the root directory of this repository 
or the <a href=http://www.gnu.org/licenses/>gnu web page</a> for details.

## Issues, Pull Requests and Discussions:

The github interfaces for issues, pull requests and discussions provide 
complementary methods of discussing bugs, problems, and proposed 
improvements:

   - Use the PSCF discussions forum (linked via the github repository)
     to ask questions about usage and to discuss possible improvements.

   - Submit issues to report bugs. 

   - Submit pull requests for proposed contributions to the source 
     code and/or documentation. 

Users may also email David Morse (morse012@umn.edu) to discuss usage 
problems or proposed improvements. 

## Issues:

Any user may submit an issue via the github web interface to report an 
error. The following guidelines apply:

  - Issues may be submitted to report bugs in the source code and
    errors in the documentation.

  - If you are unsure whether a problem that you encounter when using
    PSCF is the result of a bug, consider first sending an email to 
    morse012@umn.edu or using the discussions form to discuss the 
    problem.

  - The maintainer of PSCF often uses github issues to describe plans 
    for new features and other proposed changes. We ask, however, that 
    other users instead use email and/or github discussions to discuss 
    proposals for improvement that are not bug or error reports.

## Pull Requests:

To submit a pull request, a user must have a github account and must
have created a clone of pscfpp in their github account by forking the 
main PSCF repository, dmorse/pscfpp.  Pull requests may then be 
submitted via the github web interface. 

  - All pull requests should be directed to the "feature-pull" branch 
    of the pscfpp repository, and not to the "devel" or "main" branch. 
    To do this, use the menu provided in the github web interface for 
    submitting a pull request that allows you choose the branch to
    which the pull request should be applied.  If a pull request is 
    accepted, the maintainer will first appply it to the "feature-pull" 
    branch, and later merge that branch into the devel branch.  The 
    "feature-pull" branch is generally kept synchronized with the 
    devel branch. 

  - The first pull request from a new contributor should be one in which 
    they accept the contributor licensing agreement. See discussion of 
    this below.

  - Before submitting a pull request that introduces a new feature, we 
    recommend that you discuss the proposal with the maintainer, either 
    by sending an email to morse012@umn.edu, or by starting a discussion 
    in the github Discussions space. More generally, we recommend prior 
    discussion of contributions that involve nontrivial design decisions.

  - Pull requests should be created by starting from a branch that is 
    identical to a recent version of the devel branch and making the
    minimum necessary modifications. We recommend that you merge the 
    "devel" branch of the main repository back into the branch that 
    you are using to create a pull request shortly before submitting 
    the request, in order to incorporate any recent changes to this 
    branch made by others.

  - Contributors are expected to thoroughly test modified code before 
    submitting a pull request. At a minimum, you must run all relevant 
    PSCF unit tests, and confirm that they all pass. If you modify any 
    code that could be used by the GPU-enabled version of pscf_rp, you 
    should run all unit tests that are compiled when both C++ and CUDA 
    specializations of that program have been enabled, e.g., by invoking 
    ./setopts -c1 -u1 before compiling. To run tests that use CUDA, you 
    will need access to computer that has the NVIDIA developer toolkit 
    and an appropriate NVIDIA GPU.

  - Pull requests that add new features should include new unit tests 
    that validate the correctness of any new features. Users are also
    encouraged to contribute unit tests that test any corrections that
    they made to fix bugs.

  - A commit that fixes an open issue should include the text 
    "Fix #(number)" in the brief commit message, where (number) denotes 
    the number assigned to the relevant issue, preceded by a bash (#) 
    symbol. For example, the brief message for a commit that fixes 
    Issue 268 should contain the text "Fix #268".

  - A pull request may contain changes created by several related 
    commits. You may squash commits if you wish, but we do not require
    it.

  - Each pull request should be documented by a comment that summarizes 
    the purpose of the pull request and the nature of associated changes.
    This comment can be provided in the web interface for comments on 
    pull requests.

## Contributor Licensing Agreement (CLA):

Contributors to PSCF must agree to a contributor licensing agreement (CLA) 
that is located in the file docs/cla/agreement_v1.txt within the PSCF
repository. Pull requests will only be accepted from contributors who have 
accepted this agreement. Acceptance of the CLA is not required in order to 
submit issues or to participate in discussions.

The PSCF agreement allows you to retain the copyright for any code or 
documentation that you contribute to the project, and simply gives the 
maintainers of PSCF permission to redistribute it as part of material that 
is distributed under the GNU GPL copyleft license or other open source 
license. Please see the file docs/cla/agreement_v1.txt for full details.  

To execute the PSCF CLA, please send a pull request to feature-pull that 
contains a single commit in which you should:

  - Add your own name to the list of contributors given in the file 
    docs/cla/contributors.txt, followed by the date. This is a list 
    of contributors who have executed the CLA. The date listed after 
    each name is the date on which the agreement was executed by that 
    contributor. Do not change anything else in this file.

  - Give the commit a brief message that says "Accepting CLA"

  - Add a longer message to that commit that says "I agree to the terms 
    of the PSCF contributor licensing agreement given in the file 
    docs/cla/agreement_v1.txt within this repository." 

The git commit and pull request that are used to accept the CLA should 
not make any changes to the repository other than to add the new
contributors name and a date to the file docs/cla/contributors.txt.

Each contributor only needs to go through this process once to execute
the agreement. The agreement is written so as to apply to all past and
future contributions to PSCF by an individual contributor.

