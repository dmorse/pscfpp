
.. _group-page:

*************
Space Groups
*************

The symbol for a space group may be entered as the value of "space_group" in
the BASIS section of the parameter file. The tables below list the allowed
space group symbols.

1D Space Groups
===============

The only possible nontrivial symmetry for a one-dimensional lamellar
phase is inversion symmetry. There are thus only two possible groups:
The centrosymmetric group (group -1) and the non-centrosymmetric group
(group 1). Fields for a centrosymmetric lamellar phase are expanded
using a basis of cosine waves, while fields for a non-centrosymmetric
phase are expanded using a basis that contains both cosines and
sine waves.

======== ======  =================
Number   Symbol  Comments
======== ======  =================
1        -1      Inversion symmetry
2         1      No symmetry
======== ======  =================


2D Space Groups
===============

The names of all 17 possible 2D plane groups are given below in the text format
expected by PSCF. The format used in PSCF for both 2D and 3D space group names
is based on the names used in the international tables of crystallography, but
allows space group names to be written as simple ascii text strings, which
contain spaces between elements of the space group name.

 ====== ======== ==============
 Number Symbol   Lattice System
 ====== ======== ==============
 1      p 1      oblique
 2      p 2      oblique
 3      p m      rectangular
 4      p g      rectangular
 5      c m      rectangular
 6      p 2 m m  rectangular
 7      p 2 m g  rectangular
 8      p 2 g g  rectangular
 9      c 2 m m  rectangular
 10     p 4      square
 11     p 4 m m  square
 12     p 4 g m  square
 13     p 3      hexagonal
 14     p 3 m 1  hexagonal
 15     p 3 1 m  hexagonal
 16     p 6      hexagonal
 17     p 6 m m  hexagonal
 ====== ======== ==============


3D Space Groups
===============

The names of all possible 3D space groups are given below in the text format
expected by PSCF. These names are based on the names given in Hermann-Mauguin
or "international" notation used the international tables of crystallography,
but are given in a format that allows space group names to be written as simple
ascii text strings, with no special symbols or subscripts. In this format, for
example, the space group :math:`Ia\overline{3}d` of the gyroid phase (space
group 230) is written as "I a -3 d".

**Name conventions**

The following conventions are used to convert standard Hermann-Mauguin space
group symbols into text strings:

   * A single space is introduced between different elements of the space
     group name, with a few exceptions described below.

   * Integers with overbars in the Hermann-Mauguin symbol, which indicate
     inversion (:math:`\overline{1}`) or a 3-, 4- or 6-fold rotoinversion
     axis (:math:`\overline{3}`, :math:`\overline{4}`, or :math:`\overline{6}`),
     are indicated in the PSCF text string by placing a "-" sign before
     the integer. Thus for, example, :math:`\overline{3}` is replaced by
     "-3" in the ascii identifier "I a -3 d".

   * Integers with subscripts, such as :math:`4_2`, which represent screw
     axes, are indicated in the text representation by placing the two
     integers directly adjacent, with no intervening white space. Thus,
     for example, :math:`4_2` is replaced by "42".

   * Symbols that are separated by a slash appear with no space on either
     side of the slash.

   * Different "settings" of the same space group, which correspond to
     different definitions of the origin of space in the definition of
     the symmetry elements, are indicated by a colon followed by an
     integer label at the end of the space group.

**Space Groups**

 ========  =================
  Number   Symbol
 ========  =================
    1      P_1
    2      P_-1
    3_     P_1_2_1
    4      P_1_21 1
    5      C_1_2_1
    6      P_1_m_1
    7      P_1_c_1
    8      C_1_m_1
    9      C_1_c_1
   10      P_1_2%m_1
   11      P_1_21%m_1
   12      C_1_2%m_1
   13_     P_1_2%c_1
   14      P_1_21%c_1
   15      C_1_2%c_1
   16      P_2_2_2
   17      P_2_2_21
   18      P_21_21_2
   19      P_21_21_21
   20      C_2_2_21
   21      C_2_2_2
   22      F_2_2_2
   23_     I_2_2_2
   24      I_21_21_21
   25      P_m_m_2
   26      P_m_c_21
   27      P_c_c_2
   28      P_m_a_2
   29      P_c_a_21
   30      P_n_c_2
   31_     P_m_n_21
   32      P_b_a_2
   33_     P_n_a_21
   34      P_n_n_2
   35      C_m_m_2
   36      C_m_c_21
   37      C_c_c_2
   38      A_m_m_2
   39      A b_m_2
   40      A_m_a_2
   41_     A b_a_2
   42_     F_m_m_2
   43_     F_d_d 2
   44      I_m_m_2
   45      I_b_a_2
   46      I_m_a_2
   47      P_m_m_m
   48      P_n_n_n:2
   48      P_n_n_n:1
   49      P_c_c_m
   50      P_b_a_n:2
   50      P_b_a_n:1
   51      P_m_m_a
   52      P_n_n_a
   53_     P_m_n_a
   54      P_c_c_a
   55      P_b_a_m
   56      P_c_c_n
   57      P_b_c_m
   58      P_n_n_m
   59      P_m_m_n:2
   59      P_m_m_n:1
   60      P_b_c_n
   61      P_b_c_a
   62      P_n_m_a
   63_     C_m_c_m
   64      C_m_c_a
   65      C_m_m_m
   66      C_c_c_m
   67      C_m_m_a
   68      C_c_c_a:2
   68      C_c_c_a:1
   69      F_m_m_m
   70      F_d_d_d_:2
   70      F_d_d_d_:1
   71      I_m_m_m
   72      I_b_a_m
   73_     I_b_c_a
   74      I_m_m_a
   75      P_4
   76      P_41
   77      P_42
   78      P_43
   79      I_4
   80      I_41
   81      P_-4
   82      I_-4
   83_     P_4%m
   84      P_42%m
   85      P_4%n:2
   85      P_4%n:1
   86      P_42%n:2
   86      P_42%n:1
   87      I_4%m
   88      I_41%a:2
   88      I_41%a:1
   89      P_4_2_2
   90      P_4_21_2
   91      P_41_2_2
   92      P_41_21_2
   93_     P_42_2_2
   94      P_42_21_2
   95      P_43_2_2
   96      P_43_21_2
   97      I_4_2_2
   98      I_41_2_2
   99      P_4_m_m
  100      P_4 b_m
  101      P_42_c_m
  102      P_42_n_m
  103_     P_4_c_c
  104      P_4_n_c
  105      P_42_m_c
  106      P_42_b_c
  107      I_4_m_m
  108      I_4_c_m
  109      I_41_m_d
  110      I_41_c_d
  111      P_-4_2_m
  112      P_-4_2_c
  113_     P_-4_21_m
  114      P_-4_21_c
  115      P_-4_m_2
  116      P_-4_c_2
  117      P_-4 b 2
  118      P_-4_n_2
  119      I_-4_m_2
  120      I_-4_c_2
  121      I_-4_2_m
  122      I_-4_2_d
  123_     P_4%m_m_m
  124      P_4%m_c_c
  125      P_4%n_b_m:2
  125      P_4%n_b_m:1
  126      P_4%n_n_c:2
  126      P_4%n_n_c:1
  127      P_4%m_b_m
  128      P_4%m_n_c
  129      P_4%n_m_m:2
  129      P_4%n_m_m:1
  130      P_4%n_c_c:2
  130      P_4%n_c_c:1
  131_     P_42%m_m_c
  132      P_42%m_c_m
  133_     P_42%n_b_c:2
  133_     P_42%n_b_c:1
  134      P_42%n_n_m:2
  134      P_42%n_n_m:1
  135      P_42%m_b_c
  136      P_42%m_n_m
  137      P_42%n_m_c:2
  137      P_42%n_m_c:1
  138      P_42%n_c_m:2
  138      P_42%n_c_m:1
  139      I_4%m_m_m
  140      I_4%m_c_m
  141_     I_41%a_m_d_:2
  141_     I_41%a_m_d_:1
  142_     I_41%a_c_d_:2
  142_     I_41%a_c_d_:1
  143_     P_3
  144      P_31
  145      P_32
  146      R_3:H
  146      R_3:R
  147      P_-3
  148      R_-3:H
  148      R_-3:R
  149      P_3_1_2
  150      P_3_2_1
  151      P_31_1_2
  152      P_31_2_1
  153_     P_32_1_2
  154      P_32_2_1
  155      R_3_2:H
  155      R_3_2:R
  156      P_3_m_1
  157      P_3_1_m
  158      P_3_c_1
  159      P_3_1_c
  160      R_3_m:H
  160      R_3_m:R
  161      R_3_c:H
  161      R_3_c:R
  162      P_-3_1_m
  163_     P_-3_1_c
  164      P_-3_m_1
  165      P_-3_c_1
  166      R_-3_m:H
  166      R_-3_m:R
  167      R_-3_c:H
  167      R_-3_c:R
  168      P_6
  169      P_61
  170      P_65
  171      P_62
  172      P_64
  173_     P_63
  174      P_-6
  175      P_6%m
  176      P_63%m
  177      P_6_2_2
  178      P_61_2_2
  179      P_65_2_2
  180      P_62_2_2
  181      P_64_2_2
  182      P_63_2_2
  183_     P_6_m_m
  184      P_6_c_c
  185      P_63_c_m
  186      P_63_m_c
  187      P_-6_m_2
  188      P_-6_c_2
  189      P_-6 2_m
  190      P_-6 2_c
  191      P_6%m_m_m
  192      P_6%m_c_c
  193_     P_63%m_c_m
  194      P_63%m_m_c
  195      P_2_3
  196      F_2_3
  197      I_2_3
  198      P_21 3
  199      I_21 3
  200      P_m_-3
  201      P_n_-3:2
  201      P_n_-3:1
  202      F_m_-3
  203_     F_d_-3:2
  203_     F_d_-3:1
  204      I_m_-3
  205      P_a_-3
  206      I_a_-3
  207      P_4_3_2
  208      P_42_3_2
  209      F_4_3_2
  210      F_41_3_2
  211      I_4_3_2
  212      P_43_3_2
  213_     P_41_3_2
  214      I_41_3_2
  215      P_-4_3_m
  216      F_-4_3_m
  217      I_-4_3_m
  218      P_-4_3_n
  219      F_-4_3_c
  220      I_-4_3_d
  221      P_m_-3_m
  222      P_n_-3_n:2
  222      P_n_-3_n:1
  223_     P_m_-3_n
  224      P_n_-3_m:2
  224      P_n_-3_m:1
  225      F_m_-3_m
  226      F_m_-3_c
  227      F_d_-3_m_:2
  227      F_d_-3_m_:1
  228      F_d_-3_c:2
  228      F_d_-3_c:1
  229      I_m_-3_m
  230      I a_-3_d
 ========  =================

.. _group-symmetry-sec:

Symmetry Elements
=================

A list of all of the symmetry elements of any space group can be output to file by placing a "OUTPUT_GROUP" command in the parameter file at any point after the "BASIS" section.

Every space group symmetry can be expressed mathematically as an operation

.. math::

   \textbf{r} \rightarrow \textbf{A}\textbf{r}
                    + \textbf{t}

Here, :math:`\textbf{r} = [r_{1}, \ldots, r_{D}]^{T}` is a dimensionless
D-element column vector containing the components of a position within
the unit cell in a basis of Bravais lattice vectors, :math:`\textbf{A}`
is a :math:`D \times D` matrix that represents a point group symmetry
operation (e.g., identity, inversion, rotation about an axis, or
reflection through a plane), and :math:`\textbf{t}` is a dimenionless
D-element colummn vector that (if not zero) represents a translation
by a fraction of a unit cell. Every group contains an identity element in
which :math:`\textbf{A}` is the identity matrix and :math:`\textbf{t}=0`.

The elements of the column vectors :math:`\textbf{r}` and :math:`\textbf{t}`
in the above are dimensionless components defined using a basis of Bravais
basis vectors. The position :math:`\textbf{r} = [1/2, 1/2, 1/2]^{T}` thus
always represents the center of a 3D unit cell. The Cartesian representation
of a position vector is instead given by a sum

.. math::

   \sum_{i=1}^{D} r_{i}\textbf{a}_{i}


in which :math:`\textbf{a}_{i}` is the Cartesian representation of
Bravais lattice vector number i. The elements of the dimensionless
translation vector :math:`\textbf{t}` are always either zero or
simple fractions such as 1/2, 1/4, or 1/3. For example, a symmetry
element in a 3D BCC lattice in which :math:`\textbf{A}` is the identity
matrix and :math:`\textbf{t} = [1/2, 1/2, 1/2]^{T}` represents the
purely translational symmetry that relates the two equivalent positions
per cubic unit cell in a BCC lattice. Similarly, a glide plane in
a 3D crystal is represented by a diagonal :math:`\textbf{A}` matrix
with values of :math:`\pm 1` on the diagonal that represents
inversion through a plane and a translation vector that represents
a translation by half a unit cell within that plane.

The OUTPUT_GROUP command outputs a list of symmetry elements in
which each element is displayed by showing the elements of the
matrix :math:`\textbf{A}` followed by elements of the associated
column vector :math:`\textbf{t}`.

The Bravais lattice vectors used internally by PSCF for cubic, tetragonal,
and orthorhombic 3D systems are orthogonal basis vectors for the simple
cubic, tetragonal, or orthorhombic unit cells, which are aligned along
the x, y, and z axes of a Cartesian coordinate system. Similarly, the
basis vectors used for the 2D square and rectangular space groups are
orthogonal vectors which form a basis for a cubic or rectangular
unit cell. The grid used to solve the modified diffusion equation is
based on the same choice of unit cell and, thus for example, uses a
regular grid within a cubic unit cell to represent fields in a BCC or
FCC lattice.  For body-centered and space-centered lattice systems,
it is worth nothing that this unit cell not a primitive (minimum
size) unit cell of the crystal: For example, a cubic unit cell actually
contains 2 equivalent primitive unit cells of a BCC lattice or 4
primitive cells of an FCC lattice.

One consequence of the fact that PSCF does not always use a primitive
unit cell is that, in the Fourier expansion of the omega and rho fields,
the Fourier coefficients associated with some sets of symmetry-related
wavevectors (some "stars") are required to vanish in order to satisfy
the requirement that the field be invariant under all elements of the
specified space group. The rules regarding which stars must have
vanishing Fourier coefficients are the same as the rules for systematic
cancellations of Bragg reflections in X-ray or neutron scattering from
a crystal of the specified space group. The procedure used by PSCF to
construct symmetry adapted basis functions automatially identifies and
accounts for these systematic cancellations.

