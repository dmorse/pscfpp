
# Contributing to PSCF 

PSCF has thus far been developed primarily at the University of Minnesota 
by David Morse and student co-authors, but we welcome contributions from 
elsewhere.  It is currently maintained by David Morse (morse012@umn.edu)
in repository github.com/dmorse/pscfpp.  Users are encouraged to contribute 
to PSCF by submitting github issues to report bugs, submitting github pull 
requests, and participating in discussions.

PSCF is free, open source software. It is distributed under the terms 
of the GNU General Public License (GPL) as published by the Free Software 
Foundation, either version 3 of the License or (at your option) any 
later version. PSCF is distributed without any warranty, without even 
the implied warranty of merchantability or fitness for a particular 
purpose.  See the LICENSE file in this repository or the
<a href=http://www.gnu.org/licenses/>gnu web page</a> for details.

## Issues, Pull Requests and Discussions:

The github interfaces for issues, pull request and discussions provide 
complementary methods of discussing bugs, problems, and proposed changes.

   - Submit github Issues to report bugs. If you are unsure whether an 
     apparent problem is a bug, consider first sending an email to 
     morse012@umn.edu or using github Discussions page to discuss the 
     problem. 

   - Submit github pull requests for proposed contributions to the source 
     code and/or documentation. 

   - Use the github Discussions page to ask questions about usage and to 
     carry on conversations about limitations of the sofware and possible 
     improvements.

Users may also directly email David Morse (morse012@umn.edu) to discuss 
usage or proposed improvements. The maintainer(s) of PSCF may use github 
issues to document plans for new features and other proposed changes that
are not bugs, but others should submit issues only to document bugs, and 
use other mechanisms to discuss proposed improvements.

## Contributor Licensing Agreement (CLA):

Contributors to PSCF must agree to a contributor licensing agreement (CLA) 
that is located in the file docs/cla/agreement_v1.txt. Pull requests will 
only be accepted from contributors who have accepted this agreement. 
Acceptance of the CLA is not required in order to submit issues or to 
participate in discussions. 

The PSCF agreement allows you to retain the copyright for any code or 
documentation that you contribute to the project, but gives the maintainer 
of PSCF permission to redistribute it as part of material that is 
distributed under the GNU GPL copyleft license or another open source
license.  Please see the file docs/cla/agreement_v1.txt for full details.  

To accept the PSCF CLA, please include a commit in your first pull request
in which you:

  - Add your own name to the list of contributors given in the file 
    docs/cla/contributors.txt. This is a list of contributors who have
    accepted the agreement.

  - Give this commit a brief comment that says "Accepting CLA"

  - Add a longer comment to that commit that says "I accept the terms 
    of contributor licensing agreement given in file 
    docs/cla/agreement_v1.txt." 

The commit that is used to accept the CLA should not change any files 
other than the file docs/cla/contributors.txt, though it may be 
included in the same pull request with other commits that change
other files.

Each contributor only needs to go through this processs once to accept 
the agreement.

## Pull Requests:

To submit a pull request, a user must have a github account and a
github repository that is a fork of the PSCF repository. Pull requests 
are submitted via Pull requests may be submitted via the github web
interface.

  - All pull requests should be directed to the "pull-request" branch 
    of the pscfpp repository, and not to the "devel" or "main" branch. 
    If a pull request is accepted, the maintainer will first accept it 
    into the "pull-request" branch, and then do some testing in that
    branch before merging it into the devel branch.  The "pull-request" 
    branch is generally kept synchronized with the devel branch.

  - Before submitting a pull request that introduces a new feature, 
    we recommend that you discuss the proposal with the maintainer,
    either by sending an email to morse012@umn.edu, or by starting 
    a discussion in the github Discussions space. More generally, we 
    recommend prior discussion of most nontrivial design decisions.

  - Thoroughly test all code before submitting a pull request. At a
    minimum, you must run all relevant PSCF unit tests, including 
    those that require access to an an NVIDIA GPU.

  - Pull requests that add a new feature should also include new unit
    tests that test the correctness of the new feature.

  - A pull request that fixes an issue should include the text 
    "Fix #(number)" in the brief commit message, where (number) denotes
    the number assigned to the relevant issue, which must be preceded
    by a bash (#) symbol. The brief message for a commit that fixes 
    issue 268 should include the text "Fix #268".

  - Please consider including a detailed message that explains the
    purpose of the commit and the nature of any changes in any commit
    that is submitted as part of a pull request. 

