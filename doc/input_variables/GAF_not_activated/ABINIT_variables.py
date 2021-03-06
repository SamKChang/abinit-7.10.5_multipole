variables={
'accesswff': {
'definition': "ACCESS to WaveFunction Files ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a> """,
'vartype': "integer parameter ",
'default': """0. However, if mpi_io is available, <b>accesswff</b> will be set to 1 for the datasets for which <a href="varpar.html#paral_kgb">paral_kgb</a>=1, while an explicit mention of <b>accesswff</b> in the input file will override this intermediate default.""",
'text': """Governs the method of access to the
internal wavefunction files. Relevant only for the wavefunctions
files for which the corresponding "mkmem"-type variable is zero, that
is, for the wavefunctions that are not kept in core memory.
<ul>
<li>0 =&gt; Use standard Fortran IO routines</li>
<li>1 =&gt; Use MPI/IO routines</li>
<li>2 =&gt; Directly use NetCDF routines (this option is not available)</li>
<li>3 =&gt; Use ETSF_IO routines, creating NetCDF files according to the ETSF specification.</li>
</ul>
<br>
In case <b>accesswff</b>=1, note the following. MPI/IO routines might be much more efficient than usual Fortran IO
routines in the case of a large number of processors, with a pool of
disks attached globally to the processors, but not one disk attached
to each processor. For a cluster of workstations, where each processor
has his own temporaries, the use of <b>accesswff</b>=0 might be perfectly
allright. This option is useful only if one is using the band-FFT parallelism.
MPI/IO routines are available in the MPI-2 library, but usually not in the MPI-1 library. So, perhaps you cannot
use <b>accesswff</b>=1.
<br>
In case <b>accesswff</b>=3, note that not only the wavefunctions will be written using the ETSF_IO routines,
but also, the same input variable governs the writing of the density and potential, that can also be
written using ETSF_IO routines. In order to use <b>accesswff</b>=3, you need to have the plug-in library ETSF_IO
working (see the documentation of the build system).
References :
<ul>
<li>
"Specification of an extensible and portable file format for electronic structure and crystallographic data",
X. Gonze, C.-O. Almbladh, A. Cucca, D. Caliste, C. Freysoldt, M. Marques, V. Olevano, Y. Pouillon, M.J. Verstraete,
Comput. Mat. Science 43, 1056 (2008)
</li>
<li>
"Sharing electronic structure and crystallographic data with ETSF_IO",
D. Caliste, Y. Pouillon, M.J. Verstraete, V. Olevano, X. Gonze,
Comput. Physics Communications 179, 748 (2008)
</li>
<li>
see also <a href="http://www.etsf.eu/fileformats">http://www.etsf.eu/fileformats</a>.
</li>
</ul>
"""
},
'acell': {
'definition': "CELL lattice vector scaling",
'section': "varbas",
'category': """<a href="../users/abinit_help.html#evolving">EVOLVING</a>,<a href="../users/abinit_help.html#length">LENGTH</a> """,
'vartype': """real array acell(3), represented internally as acell(3,<a href="varrlx.html#nimage">nimage</a>)""",
'default': "3*1 (in Bohr).",
'text': """Gives the length scales by which
dimensionless primitive translations (in <a href="varbas.html#rprim">rprim</a>) are
to be multiplied.  By default, given in Bohr atomic units
(1 Bohr=0.5291772108 Angstroms), although Angstrom can be specified,
if preferred, since <b>acell</b> has the
'<a href="../users/abinit_help.html#dimensions">LENGTH</a>' characteristics.
See further description of <b>acell</b> related to the
<a href="varbas.html#rprim">rprim</a> input variable,
the <a href="varbas.html#scalecart">scalecart</a> input variable,
and the associated internal <a href="varbas.html#rprimd">rprimd</a> input variable.
<p>
Note that <b>acell</b> is NOT the length of the conventional orthogonal basis vectors, but the scaling factors of the primitive vectors.
Use <a href="varbas.html#scalecart">scalecart</a> to scale the cartesian coordinates. """
},
'algalch': {
'definition': "ALGorithm for generating ALCHemical pseudopotentials",
'section': "vargs",
'category': " ",
'vartype': """integer array algalch(<a href="vargs.html#ntypalch">ntypalch</a>)  """,
'default': "1 for all indices",
'text': """<p> Used for the generation of alchemical pseudopotentials, that is,
when <a href="vargs.html#ntypalch">ntypalch</a> is non-zero.
<p>Give the algorithm to be used to
generate the <a href="vargs.html#ntypalch">ntypalch</a> alchemical potentials
from the different <a href="vargs.html#npspalch">npspalch</a> pseudopotentials
dedicated to this use.
<p> Presently, <b>algalch</b> can only have the value 1, that is :
<ul>
<li>the local potentials are mixed, thanks to the <a href="vargs.html#mixalch">mixalch</a>
mixing coefficients</li>
<li>the form factors of the non-local projectors are all preserved, and all considered
to generate the alchemical potential</li>
<li>the scalar coefficients of the non-local projectors are multiplied by the proportion
of the corresponding type of atom that is present in <a href="vargs.html#mixalch">mixalch</a></li>
<li>the characteristic radius for the core charge is a
linear combination of the characteristic radii of the core charges,
build with the <a href="vargs.html#mixalch">mixalch</a> mixing
coefficients</li>
<li>the core charge function f(r/rc) is a linear combination
of the core charge functions, build with the <a href="vargs.html#mixalch">mixalch</a>
mixing coefficients</li>
</ul>
Later, other algorithms for the mixing might be included.
<p>Note that alchemical mixing cannot be used with PAW.
"""
},
'amu': {
'definition': "Atomic Mass Units ",
'section': "varrlx",
'category': "EVOLVING ",
'vartype': """real array amu(<a href="varbas.html#ntypat">ntypat</a>) """,
'default': "provided by a database of atomic masses.",
'text': """A database of atomic masses is provided, giving
default values.
Note that the default database uses mixed isotope masses (for Carbon
the natural occurence of Carbon 13 is taken into account).
The values are those recommended by the commission on Atomic Weights
and
Isotopic Abundances, Inorganic Chemistry Division, IUPAC, in
<i>Pure Appl. Chem.</i> <b>60</b>, 841 (1988).
For Tc, Pm, Po to Ac, Pa and beyond U,
none of the isotopes has a half-life greater than 3.0d10 years, and
the values provided in the database do not come from that source.
<p>
For alchemical pseudoatoms, the masses of the constituents
atoms are mixed, according to the alchemical miwing
coefficients <a href="vargs.html#mixalch">mixalch</a>
<p> In most cases, the use of <b>amu</b>
will be as a static (non-evolving) variable. However, the possibility to have
different values of <b>amu</b> for different images has been coded. A population of
cells with different atomic characteristics can thus be considered,
and can be made to evolve, e.g. with a genetic algorithm (not coded in v7.0.0 though).
"""
},
'angdeg': {
'definition': "ANGles in DEGrees  ",
'section': "varbas",
'category': """<a href="../users/abinit_help.html#not_internal">NOT INTERNAL</a> """,
'vartype': "real array angdeg(3) ",
'default': """No Default (use <a href="varbas.html#rprim">rprim</a> as Default).""",
'text': """Gives the angles between directions of
primitive vectors of the unit cell (in degrees),
as an alternative to the input array <a href="varbas.html#rprim">rprim</a> .
Will be used to set up <a href="varbas.html#rprim">rprim</a>,
that, together with the array <a href="varbas.html#acell">acell</a>, will be used to define the
primitive vectors.
<ul>
<li><b>angdeg</b>(1) is the angle between the 2nd and 3rd vectors,</li>
<li><b>angdeg</b>(2) is the angle between the 1st and 3rd vectors,</li>
<li><b>angdeg</b>(3) is the angle between the 1st and 2nd vectors,</li>
</ul>
If the three angles are equal within 1.0d-12 (except if they are exactly 90 degrees),
the three primitive
vectors are chosen so that the trigonal symmetry that exchange
them is along the z cartesian axis :
<pre>
R1=( a  ,           0,c)
R2=(-a/2, sqrt(3)/2*a,c)
R3=(-a/2,-sqrt(3)/2*a,c)
</pre>
where a<sup>2</sup>+c<sup>2</sup>=1.0d0
<br>
If the angles are not all equal (or if they are all 90 degrees), one will have the following
generic form :
<ul>
<li>R1=(1,0,0)</li>
<li>R2=(a,b,0)</li>
<li>R3=(c,d,e)</li>
</ul>
where each of the vectors is normalized,
and form the desired angles with the others."""
},
'atvshift': {
'definition': "ATomic potential (V) energy SHIFTs ",
'section': "varff",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a> """,
'vartype': """real array atvshift (<a href="varff.html#natvshift">natvshift</a>, <a href="varbas.html#nsppol">nsppol</a>, <a href="varbas.html#natom">natom</a>) """,
'default': "a set of 0.0d0.",
'text': """Defines for each atom and each spin channel (at present, can only be used
with <a href="varbas.html#nsppol">nsppol</a>=1 or 2, like the +U scheme),
a possible potential shift, for the d
(with <a href="varpaw.html#lpawu">lpawu</a>=2,
<a href="varff.html#natvshift">natvshift</a>=5),
or f states
(with <a href="varpaw.html#lpawu">lpawu</a>=3,
<a href="varff.html#natvshift">natvshift</a>=7).
In the case of d states, and 2 spin channels, a set of 10 numbers for
each atom must be defined.
The first set of 5 numbers corresponds to real spherical harmonics
m=-2 to m=+2 for the spin-up channel,
the second set of 5 numbers corresponds to real spherical harmonics
m=-2 to m=+2 for the spin-down channel.
In the case of f states, the same ordering applies, for sets of 7 numbers,
corresponding to m=-3 to m=+3.
<br>
<a href="varpaw.html#usepawu">usepawu</a> should be non-zero,
<a href="varpaw.html#lpawu">lpawu</a> should be 2 or 3.
"""
},
'awtr': {
'definition': "evaluate the Adler-Wiser expression of $\chi^{(0)}_{KS}$ assuming Time-Reversal ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>  """,
'vartype': "integer parameter  ",
'default': "1 (was 0 prior to 5.8)",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=3, that is, screening calculations.
<p>
This input variable defines whether the irreducible polarizability $\chi^{(0)}_{KS}$ is evaluated
taking advantage of time-reversal symmetry or not.

<ul>
<li>  0 =&gt; Use the "standard" Adler-Wiser expression without assuming time-reversal symmetry.
In this case, the irreducible polarizability is calculated summing over all possible electronic
transitions (both resonant and antiresonant).
<li>  1 =&gt; Take advantage of time-reversal symmetry to halve the number of transitions to be
explicitly considered. This method leads to a decrease in the CPU time by a factor two with respect
to the <b>awtr</b>=0 case.
</ul>

<p>
Note that the parallel algorithm <a href="varpar.html#gwpara">gwpara</a>=2 is not compatible with
the choiche <b>awtr</b>=0.

"""
},
'bandpp': {
'definition': "BAND Per Processor ",
'section': "varpar",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a> """,
'vartype': "integer parameter ",
'default': "1.",
'text': """Control the size of the block in the LOBPCG algorithm.
This keyword works only with <a href="varpar.html#paral_kgb">paral_kgb</a>=1 and has to be either 1 or a multiple of 2.
<br><br>-- With <a href="varpar.html#npband">npband</a>=1:
<ul>
<li>1 =&gt; band-per-band algorithm</li>
<li>n =&gt; The minimization is performed using <a href="varbas.html#nband">nband</a>/n blocks of n bands.</li>
</ul>
Note: <a href="varbas.html#nband">nband</a>/n has to be an integer.
<br><br>-- With <a href="varpar.html#npband">npband</a>/=1:
<ul>
<li>1 =&gt; The minimization is performed using <a href="varbas.html#nband">nband</a>/<a href="varpar.html#npband">npband</a> blocks of <a href="varpar.html#npband">npband</a> bands.</li>
<li>n =&gt; The minimization is performed using <a href="varbas.html#nband">nband</a>/(<a href="varpar.html#npband">npband</a>*n) blocks of <a href="varpar.html#npband">npband</a>*n bands.</li>
</ul>
Note: <a href="varbas.html#nband">nband</a>/(<a href="varpar.html#npband">npband</a>*n) has to be an integer.
<br><br>By minimizing a larger number of bands together in LOBPCG, we increase the convergency of the residual.
The better minimization procedure (as concerns the convergency, but not as concerns the speed) is generally
performed by using <b>bandpp</b>*<a href="varpar.html#npband">npband</a>=<a href="varbas.html#nband">nband</a>.
Put <b>bandpp</b>=2 when <a href="vardev.html#istwfk">istwfk</a>=2 (the time spent in FFTs is divided by two).
"""
},
'bdberry': {
'definition': "BanD limits for BERRY phase ",
'section': "varff",
'category': " ",
'vartype': "integer array bdberry(4)  ",
'default': "4*0.",
'text': """<p> Used for non-zero values of <a href="varff.html#berryopt">berryopt</a>.
<p>Give the lower band and the upper band of the set of bands
for which the Berry phase must be computed.
Irrelevant if <a href="varff.html#nberry">nberry</a> is not positive.
When <a href="varbas.html#nsppol">nsppol</a> is 1 (no spin-polarisation),
only the two first numbers, giving the lower and highest
bands, are significant. Their occupation number is assumed to be 2.
When <a href="varbas.html#nsppol">nsppol</a> is 2 (spin-polarized calculation),
the two first numbers give the lowest and highest
bands for spin up, and the third and fourth numbers
give the lowest and highest bands for spin down.
Their occupation number is assumed to be 1 .
<p> Presently, <b>bdberry</b> MUST be initialized by the user
in case of Berry phase calculation: the above-mentioned
default will cause an early exit.
"""
},
'bdeigrf': {
'definition': "BanD for second-order EIGenvalues from Response-Function ",
'section': "varrf",
'category': """<a href="../users/abinit_help.html#respfn">RESPFN</a>  """,
'vartype': "integer parameters ",
'default': "-1.",
'text': """Only relevant if <a href="varrf.html#ieig2rf">ieig2rf</a> = 1 or 2, that is, if the user is performing second-order eigenvalue calculations using response-functions.
<br><br>
The variable <b>bdeigrf</b> is the maximum number of bands for which the second-order eigenvalues must be calculated: the full number of bands is still used during the computation of these corrections.
<br><br>
If <b>bdeigrf</b> is set to -1, the code will automatically set <b>bdeigrf</b> equal to nband.
"""
},
'bdgw': {
'definition': "BanDs for GW calculation ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>  """,
'vartype': """integer array bdgw(2,<a href="vargw.html#nkptgw">nkptgw</a>, <a href="varbas.html#nsppol">nsppol</a>,) """,
'default': "all 0's ",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=4, that is, sigma calculations.
<p>
For each k-point with number ikptgw in the range (1:<a href="vargw.html#nkptgw">nkptgw</a>) and each
spin index isppol, <b>bdgw(1,ikptgw,isppol)</b> is the number of the lowest band for which the GW computation must be done,
and <b>bdgw(2,ikptgw,isppol)</b> is the number of the highest band for which the GW computation must be done.
<p>
When <a href="vargw.html#gwcalctyp">gwcalctyp</a> &#62&#61 20,
the quasiparticle wavefunctions are computed and represented as linear combination of some Kohn-Sham wavefunctions.
In this case <b>bdgw</b> designates the KS wavefunctions used as basis set.
For each k-point, indeed, the quasiparticle wavefunctions are expanded considering only the KS states between
<b>bdgw(1,ikptgw,isppol)</b> and <b>bdgw(2,ikptgw,isppol)</b>.
<p>
Note that the initial values given in the input file might be changed inside the code so that all
the degenerate states at a given k-point and spin are included.
This might happen when <a href="vargw.html#symsigma">symsigma</a>=1 is used or in the case of
self-consistent GW calculations.
<p>
When <a href="vargw.html#symsigma">symsigma</a>=1, indeed, the diagonal matrix elements of the
self-energy are obtained by averaging the unsymmetrized results in the degenerate subspace.
<p>
For self-consistent calculations, on the other hand, the basis set used to expand the GW
wavefunctions should include all the degenerate states belonging to the same irreducible
representation. Only in this case, indeed, the initial symmetries and energy degenerations are preserved.
"""
},
'berryopt': {
'definition': "BERRY phase OPTions ",
'section': "varff",
'category': " ",
'vartype': "integer berryopt ",
'default': "0",
'text': """Specifies the use of Berry phase for the computation
of either the polarization, the derivatives with respect to the
wavevector, or finite electric field calculations.
<ul>
<li>0 => no computation of expressions relying on a Berry phase (default)</li>
<li>1 => the computation of Berry phases is activated (berryphase routine)</li>
<li>2 => the computation of derivatives with respect to the wavevector,
thanks to the Berry phase finite-difference formula, is activated (uderiv routine)</li>
<li>3 => same as option 1 and 2 together </li>
<li>4 => finite unreduced electric field calculation (Ground state as well as phonon)</li>
<li>5 => <b> Highly experimental!</b> Finite magnetic field calculation.  </li>
<li>6 => <b> Only for experts!</b> Finite unreduced electric displacement field calculation (out loop)</li>
<li>7 => <b> Only for experts!</b> Finite unreduced electric displacement field calculation (in loop)</li>
<li>14 => finite reduced electric field calculation (Ground state as well as phonon)</li>
<li>16 => <b> Only for experts!</b> Finite reduced electric displacement field calculation (out loop)</li>
<li>17 => <b> Only for experts!</b> Finite reduced electric displacement field calculation (in loop)</li>
<li>-1 => alternative computation of Berry phases (berryphase_new routine)</li>
<li>-2 => alternative computation of derivatives with respect to the wavevector,
thanks to the Berry phase finite-difference formula (berryphase_new routine)</li>
<li>-3 => same as option -1 and -2 together </li>
<li>-5 => <b> Highly experimental!</b> Computation of magnetization in zero external field, based on Berrys Phase formula. </li>
</ul>
<p>
The other related input variables are :
<ul>
<li>in case of <b>berryopt</b>=1,2, or 3 : <a href="varff.html#bdberry">bdberry</a>
and <a href="varff.html#kberry">kberry</a>; also, <a href="varff.html#nberry">nberry</a>
must be larger than 0;</li>
<li>in case of <b>berryopt</b>=-1,-2, or -3 : the variable
<a href="varrf.html#rfdir">rfdir</a> must be used to specify the primitive
vector along which the projection of the polarization or the ddk will be computed.
For example if <b>berryopt</b>=1 and <a href="varrf.html#rfdir">rfdir</a>=1 0 0,
the projection of the polarization along the reciprocal lattice vector
G_1 is computed. In case <a href="varrf.html#rfdir">rfdir</a>=1 1 1,
ABINIT computes the projection of P along G_1, G_2 and G_3 and transforms the results
to cartesian coordinates;</li>
<li>in case of <b>berryopt</b> is negative, <a href="varff.html#berrystep">berrystep</a>
allow a computation of multiple-step Berry phase in order to accelerate the convergence.<li>
<li><a href="varff.html#efield">efield</a>,
<a href="varrf.html#rfdir">rfdir</a> in case of <b>berryopt</b>=4 ;</li>
<li><a href="varff.html#bfield">bfield</a>,
<a href="varrf.html#rfdir">rfdir</a> in case of <b>berryopt</b>=5 ;</li>
</ul>
<p>
The cases <b>berryopt</b>=-1,-2,-3, 4, 6, 7, 14, 16, 17 and -5 and 5 have to be used with
<a href="varbas.html#nsppol">nsppol</a>=1,
<a href="vargs.html#nspinor">nspinor</a>=1, and
<a href="varbas.html#occopt">occopt</a>=1.
<p>
The cases <b>berryopt</b>=-1 and 4, 6, 7, 14, 16, 17 are compatible with PAW, howevever, if in these cases one uses
<a href="varbas.html#kptopt">kptopt</a>/=3, one must also use only symmorphic symmetries (either because the space group is
symmorphic or the variable <a href="vargw.html#symmorphi">symmorphi</a> is set to zero).
<p>
For a phonon calculation under a finite electric field, respect the following procedure.
<ul>
<li>1. Run a scf ground-state calculation at zero electric field
to get wavefunctions to initialize the ground-state calculation in finite electric fields.</li>
<li>2. Run a scf ground-state calculation in finite electric field. The
electric field is controlled by the input variable <a href="varff.html#efield">efield</a>.
<b>berryopt</b> should be 4.
The input variable <a href="varbas.html#kptopt">kptopt</a> should be set to be 2.</li>
<li>3. Based on the wave functions obtained in step (2), perform phonon
calculation by setting <b>berryopt</b>=4, <a href="varbas.html#kptopt">kptopt</a>=3 and
The same value of <a href="varff.html#efield">efield</a> than in step 2.
<a href="varbas.html#nsym">nsym</a> should be set to 1 currently but this restriction may be
removed later . The other
parameters are the same as phonon calculation at zero electric field.</li>
<li>Note : the choice of k-point sampling N x N x N should be the same in the three runs
and N should be an even number.
</ul>
<p>
At present, note that <b>berryopt</b>=5 is not compatible with non-zero <a href="varrlx.html#optcell">optcell</a>.
<p>
In case of finite electric (displacement) field calculations (<b>berryopt</b>=4,6,7,14,16,17), see also the input variables
<a href="varff.html#berrysav">berrysav</a>, <a href="varff.html#dfield">dfield</a>,
<a href="varff.html#red_dfield">red_dfield</a>, <a href="varff.html#red_efield">red_efield</a>,
<a href="varff.html#ddamp">ddamp</a>

"""
},
'berrysav': {
'definition': "BERRY SAVe ",
'section': "varff",
'category': " ",
'vartype': "integer berrysav ",
'default': "0",
'text': """<ul>
<li>0 => for finite electric field calculation (<a href="varff.html#berryopt">berryopt</a>=4/14), the polarization branch will be chosen on each iteration from (-pi, pi);
for finite electric displacement field calculation(<a href="varff.html#berryopt">berryopt</a>=6/7/16/17), the polarization will be chosen to minimize the internal energy.
</li>
<li>1 => the polarization will be kept in the same branch on each iteration. A new file "POLSAVE" will be saved.
</li>
</ul>

"""
},
'berrystep': {
'definition': "BERRY phase : multiple STEP ",
'section': "varff",
'category': " ",
'vartype': "integer berrystep ",
'default': "1",
'text': """If <a href="varff.html#berryopt">berryopt</a> is negative,
this variable is used to compute multiple-steps berry phase. The single-step berry phase
is the standard calculation using strings of k-points based on overlap of Bloch function separated by dk, while the two-step
berry phase use strings use overlaps based on dk and 2*dk, the three-step use overlaps based on dk, 2*dk and 3*dk...
<p>
The default value of this variable is 1, meaning that only the single-step berry phase calculation is done.
If a larger value is set, ABINIT will compute all the multiple-step berry phase from the single-step to the
<b>berrystep</b>-step, and use the large-step values of berry phase to correct the single-step berry phase.
As of now (2011), experience is still to be gained with this procedure, but it is promising.
"""
},
'bfield': {
'definition': "magnetic B FIELD ",
'section': "varff",
'category': " ",
'vartype': "real array bfield(3)  ",
'default': "3*0.0 .",
'text': """<b>Highly experimental!</b> In case <a href="varff.html#berryopt">berryopt</a>=5,
a finite magnetic field calculation is performed. The value
of this magnetic field, and its direction is determined by <b>bfield</b>.
It must be given in atomic units in cartesian coordinates. As related to SI units
for magnetic field, 1 a.u. of magnetic field is hbar/(e_Cb a0**2), or 2.35e5 Tesla.
"""
},
'bmass': {
'definition': "Barostat mass",
'section': "varrlx",
'category': "",
'vartype': "real parameter ",
'default': "10.",
'text': """bmass is the mass of the barostat when <a
href="varrlx.html#ionmov">ionmov</a>=13 (constant pressure)
"""
},
'boxcenter': {
'definition': "BOX CENTER ",
'section': "vargs",
'category': "",
'vartype': "real array <b>boxcenter</b>(3) ",
'default': "0.5 0.5 0.5 .",
'text': """Defines the center of the box, in reduced coordinates.
At present, this information is only used in the case of
Time-Dependent DFT computation of the oscillator strength.
One must take boxcenter such as to be roughly the center of
the cluster or molecule. The default is sensible when
the vacuum surrounding the cluster or molecule has xred 0 or 1.
On the contrary, when the cluster or molecule is close to
the origin, it is better to take <b>boxcenter</b>=(0 0 0).
"""
},
'boxcutmin': {
'definition': "BOX CUT-off MINimum ",
'section': "vargs",
'category': "",
'vartype': "real ",
'default': "2.0 .",
'text': """The box cut-off ratio is the ratio between the wavefunction plane wave sphere
radius, and the radius of the sphere that can be inserted in the
FFT box, in reciprocal space. In order for the density to be exact
(in the case of plane wave, not PAW), this ratio should be at least two.
If one uses a smaller ratio, one will gain speed, at the expense of accuracy.
In case of pure ground state calculation (e.g. for the determination
of geometries), this is sensible. However,
the wavefunctions that are obtained CANNOT be used for starting response function
calculation.
"""
},
'brvltt': {
'definition': "BRaVais LaTTice type ",
'section': "vargeo",
'category': """<a href="../users/abinit_help.html#symmetriser">SYMMETRISER</a> """,
'vartype': "integer parameter ",
'default': "0.",
'text': """Set the type of Bravais lattice,
needed only if <a href="vargeo.html#spgroup">spgroup</a>/=0 .
In this case, the cell defined by <a href="varbas.html#acell">acell</a>
and <a href="varbas.html#rprim">rprim</a> or <a href="varbas.html#angdeg">angdeg</a>
should be the CONVENTIONAL cell.
<p>If brvltt=0, the code will assign
brvltt from the space group information
<a href="vargeo.html#spgroup">spgroup</a>,
and produce the symmetry operations for the conventional unit cell.
If the conventional cell is not primitive, the user should
set <a href="vargs.html#chkprim">chkprim</a>=0.
<p>
If brvltt=-1, the code will assign brvltt from
the space group information, then reduce the unit cell
to a primitive unit cell. The echo of <a href="varbas.html#acell">acell</a>
and <a href="varbas.html#rprim">rprim</a> might thus differ from those
derived directly from the input variables.
Also, the input variable
<a href="varbas.html#xred">xred</a> will refer to the
CONVENTIONAL unit cell, but its echo will refer to the
preprocessed PRIMITIVE unit cell.
There is of course no problem with
<a href="varbas.html#xangst">xangst</a> and
<a href="varbas.html#xcart">xcart</a>, as they are independent
of the unit cell.
<p>The echo of <b>brvltt</b> in the output file will be one
of the following Bravais lattices:
<br>
<ul>
<li>1 = Primitive with no associated translations</li>
<li>2 = Inner centered with (a/2 + b/2 + c/2)
associated translation </li>
<li>3 = Face centered with (a/2 + b/2; b/2 + c/2; c/2 + a/2)
associated translations</li>
<li>4 = C - centered with (a/2 + b/2) associated translation</li>
<li>5 = A - centered with (b/2 + c/2) associated translation</li>
<li>6 = B - centered with (c/2 + a/2) associated translation</li>
<li>7 = Rhombohedral lattice.</li>
</ul>
The user might also input directly these values, although
they might not be consistent
with <a href="vargeo.html#spgroup">spgroup</a>.
<p>
The space groups 146, 148, 155, 160, 161, 166, 167, when used
with <a href="vargeo.html#spgaxor">spgaxor</a>=1 (hexagonal axes) will have <b>brvltt</b>=7
and two associated translations: (2/3, 1/3, 1/3) and
(1/3, 2/3, 2/3).
<br>For more details see the space group
<a href="../users/spacegrouphelpfile.html">help file</a>.
"""
},
'bs_algorithm': {
'definition': "Bethe-Salpeter ALGORITHM ",
'section': "vargw",
'category': "BS  ",
'vartype': "integer ",
'default': "1 ",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=99, that is, Bethe-Salpeter calculations.
<p>
The bs_algorithm input variable defines the algorithm employed to calculate the macroscopic dielectric function.
Possible values are 1, 2 or 3:

<ul>
<li> 1 =&gt;
The macroscopic dielectric is obtained by performing a direct diagonalization of the exciton Hamiltonian.
Advantages: It gives direct access to the excitonic eigenvalues as well as to the oscillator strengths.
Drawbacks: It is a very CPU- and memory-consuming approach as the size of the Hamiltonian scales as (nk*nc*nv)**2.
where nk is the number of k-point in the FULL Brillouin zone, and nc and nv are the number of
conduction and valence states, respectively.
Pros: It can be used both for resonant-only and resonant+coupling calculations.
<li> 2 =&gt;
Haydock iterative method. The macroscopic dielectric function
is obtained by iterative applications of the Hamiltonian on a set of vectors in the electron-hole space.
Advantages: It is less memory demanding and usually faster
than the direct diagonalization provided that <a href="vargw.html#zcut">zcut</a> is larger than the typical
energy spacing of the eigenvalues. Drawbacks: It's an iterative method therefore
the convergence with respect to bs_haydock_niter should be checked.
It is not possible to have direct information on the exciton spectrum, oscillator strengths and excitonic wave functions.
For the time being <b>bs_algorithm</b>=2 cannot be used for calculations in which the coupling term is included.
<li> 3 =&gt;
Conjugate-gradient method
[STILL UNDER DEVELOPMENT]
Only available for resonant calculations.
</ul>
"""
},
'bs_calctype': {
'definition': "Bethe-Salpeter CALCulation TYPE ",
'section': "vargw",
'category': "BS  ",
'vartype': "integer ",
'default': "1 ",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=99, that is, Bethe-Salpeter calculations.
<p>
Possible values are 1,2,3
<ul>
<li> 1 =&gt; use the KS eigenvalues and wave functions stored in the KSS file
to construct the transition space

<li> 2 =&gt; The transition space is constructed with Kohn-Sham orbitals
but the energies are read from the external GW file

<li> 3 =&gt; QP amplitudes and energies will be read from the QPS file and used to construct H_ex.
Not coded yet
Besides <\psi|r|\psj>^QP should be calculated taking into account the non-locality of the
the self-energy in the commutator [H,r].
</ul>
"""
},
'bs_coulomb_term': {
'definition': "Bethe-Salpeter COULOMB TERM ",
'section': "vargw",
'category': "BS  ",
'vartype': "integer ",
'default': "11 ",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=99, that is, Bethe-Salpeter calculations.
This variable governs the choice among the different options that are available for the treatment
of Coulomb term of the BS Hamiltonian.
<b>bs_coulomb_term</b> is the concatenation of two digits, labelled (A) and (B).
<p>
The first digit (A) can assume the values 0,1,2:

<ul>
<li> 0 =&gt; The Coulomb term is not computed. This choice is equivalent to computing the RPA
spectrum but using the representation in transition space instead of the
more efficient approach based on the sum over states.
<li> 1 =&gt; The Coulomb term is computed using the screened interaction read from an external SCR file
(standard excitonic calculation).
<li> 2 =&gt; The Coulomb term is computed using a model screening function
(useful for convergence studies or for reproducing published results).
</ul>
<p>
The second digit (B) can assume the values: 0,1
<ul>
<li> 0 =&gt;
Use a diagonal approximation for W_GG' (mainly used for accelerating convergence studies).
<li> 1 =&gt;
The Coulomb term is correctly evaluated using the truly non-local W(r,r').
</ul>

"""
},
'bs_coupling': {
'definition': "Bethe-Salpeter COUPLING ",
'section': "vargw",
'category': "BS ",
'vartype': "integer ",
'default': "1 ",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=99, that is, Bethe-Salpeter calculations.
<p>

The <b>bs_coupling</b> input variable defines the treatment of the coupling block of the
BS Hamiltonian. Possible values are 0,1.

<ul>
<li> 0 =&gt;
The coupling block is neglected. The code runs faster and the Hamiltonian matrix requires
less memory (factor 4). It is a good approximation for the absorption spectrum
which only requires the knowledge of Im(\epsilon). The reliability of this
approximation should be tested in the case of EELF calculations.
<li> 1 =&gt;
The coupling term is included.
</ul>
"""
},
'bs_eh_cutoff': {
'definition': "Bethe-Salpeter Electron-Hole CUTOFF ",
'section': "vargw",
'category': "BS  ",
'vartype': "integer(3) ",
'default': "(-inf,+inf,+inf) ",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=99, that is, Bethe-Salpeter calculations.
<p>
It is used to define a cutoff in the e-h basis set. Only those transitions whose
energy is between bs_eh_window(1) and bs_eh_window(2) will be considered during
the construction of the e-h Hamiltonian. bs_eh_window(3) (named stripecut in EXC)
are used to define an additional cutoff in the exciton Hamiltonian. Let it and itp
indicate two e-h basis set elements with itp >= it.Only the matrix elements
corresponding to the pairs (it,itp) such that Etrans_{itp}-Etrans_{it} are considered.
"""
},
'bs_exchange_term': {
'definition': "Bethe-Salpeter EXCHANGE TERM ",
'section': "vargw",
'category': "BS  ",
'vartype': "integer ",
'default': "1 ",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=99, that is, Bethe-Salpeter calculations.
<p>

<ul>
<li> 0 =&gt;
The exchange term is not calculated. This is equivalent to neglecting local field effects
in the macroscopic dielectric function.
<li> 1 =&gt;
The exchange term is calculated and added to the excitonic Hamiltonian.
</ul>
"""
},
'bs_freq_mesh': {
'definition': "Bethe-Salpeter FREQuency MESH ",
'section': "vargw",
'category': """BS, <a href="../users/abinit_help.html#energy">ENERGY</a>  """,
'vartype': "real(3) ",
'default': "(0.0,0.0,0.01) eV ",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=99, that is, Bethe-Salpeter calculations.
<p>
bs_freq_mesh(1) defines the first frequency for the calculation of the macroscopic dielectric function.
bs_freq_mesh(2) gives the last frequency for the calculation of the macroscopic dielectric function.
If zero, bs_omegae is automatically calculated inside the code such that the frequency mesh covers the entire set of e-h transitions.
bs_freq_mesh(3) gives the step of the linear mesh used for evaluating the macroscopic dielectric function.
"""
},
'bs_hayd_term': {
'definition': "Bethe-Salpeter HAYdock TERMinator ",
'section': "vargw",
'category': "BS  ",
'vartype': "integer ",
'default': "1 ",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=99 and <a href="vargw.html#bs_algorithm">bs_algorithm</a>=2
that is, Bethe-Salpeter calculations with the Haydock iterative method.
<p>
Defines how to terminate the continued fraction expression for the dielectic function.
The terminator reduces the number of iterations needed to converge by smoothing the oscillation
in the high energy part of the spectrum

<ul>
<li> 0 =&gt;
No terminator. The contribution given by the terms missing in the Lanczos chain are set to zero.
<li> 1 =&gt; Use the terminator function. The particular expression depends on the type of calculation:
In the resonant-only case, the a_i and b_i coefficients for i &gt; niter,
are replaced by their values at i=niter.
Ehen the coupling block is included, the terminator function described in D. Rocca, R. Gebauer, Y. Saad, S. Baroni, J. Chem. Phys.
128, 154105 (2008) is used.
</ul>

"""
},
'bs_haydock_niter': {
'definition': "Bethe-Salpeter HAYDOCK Number of Iterations ",
'section': "vargw",
'category': "BS ",
'vartype': "integer ",
'default': "100 ",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=99 and <a href="vargw.html#bs_algorithm">bs_algorithm</a>=2
that is, Bethe-Salpeter calculations with the Haydock iterative method.
<p>
<b>bs_haydock_niter</b> defines the maximum number of iterations used to calculate the macroscopic dielectric function.
The iterative algorithm stops when the difference between two consecutive evaluations of the optical
spectra is less than <a href="vargw.html#bs_haydock_tol">bs_haydock_tol</a>.
"""
},
'bs_haydock_tol': {
'definition': "Bethe-Salpeter HAYDOCK TOLerance ",
'section': "vargw",
'category': "BS  ",
'vartype': "real(2) ",
'default': "(0.02,0) ",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=99 and <a href="vargw.html#bs_algorithm">bs_algorithm</a>=2
that is, Bethe-Salpeter calculations with the Haydock iterative method.
<p>
Defines the convergence criterion for the Haydock iterative method.
The iterative algorithm stops when the difference between two consecutive evaluations of the
macroscopic dielectric function is less than <b>bs_haydock_tol(1)</b>.
The sign of <b>bs_haydock_tol(1)</b> defines how to estimate the convergence error.
A negative value signals that the converge should be reached for each frequency (strict criterion),
while a positive value indicates that the converge error is estimated
by averaging over the entire frequency range (mild criterion).
<p>
bs_haydock_tol(2) defines the quantity that will be checked for convergence:
<ul>
<li>0 &rarr; both the real and the imaginary part must converge
<li>1 &rarr; only the real part
<li>2 &rarr; only the imaginary part
</ul>
"""
},
'bs_loband': {
'definition': "Bethe-Salpeter Lowest Occupied BAND",
'section': "vargw",
'category': "BS  ",
'vartype': "integer",
'default': "0",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=99, that is, Bethe-Salpeter calculations.
<p>
This variable defines the index of the lowest occupied band used for the construction of the electron-hole basis set
An additional cutoff energy can be applied by means of the bs_eh_window input variable.
"""
},
'bs_nstates': {
'definition': "Bethe-Salpeter Number of States ",
'section': "vargw",
'category': "BS ",
'vartype': "integer ",
'default': "0 ",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=99 and <a
href="vargw.html#bs_algorithm">bs_algorithm</a>=2,3
that is, Bethe-Salpeter calculations with the Haydock iterative approach or the Conjugate gradient method
<p>
<b>bs_nstates</b> defines the maximum number of excitonic states calculated in
the direct diagonalization of the excitonic matrix or in the conjugate-gradient method.
The number of states should be sufficiently large
for a correct description of the optical properties in the frequency range of interest.
"""
},
'bxctmindg': {
'definition': "BoX CuT-off MINimum for the Double Grid (PAW) ",
'section': "varpaw",
'category': "",
'vartype': "real parameter ",
'default': "2.0",
'text': """Relevant only when
<a href="varint.html#usepaw">usepaw</a>=1.
<br>The box cut-off ratio is the ratio between the wavefunction plane wave
sphere radius, and the radius of the sphere that can be inserted in the FFT box,
in reciprocal space.
<br>If the density was generated only from wavefunctions,
this ratio should at least two in order for the density to be exact. If one
uses a smaller ratio, one will gain speed, at the expense of accuracy. In case
of pure ground state calculation (e.g. for the determination of geometries),
this is sensible. However, the wavefunctions that are obtained CANNOT be used
for starting response function calculation.
<br>However, some augmentation charge is always added in PAW, and even with the box cut-off
ratio larger than two, the density is never exact. Sometimes, this ratio must be
much larger than two for the computation to be converged at the
required level of accuracy.
"""
},
'cd_customnimfrqs': {
'definition': "Contour Deformation Custom Imaginary Frequencies",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a> """,
'vartype': "integer ",
'default': "0 (off)",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=3 or 4, that is, screening or sigma calculations,
and <a href="vargw.html#gwcalctyp">gwcalctyp</a>=2, 12 or 22.
<p>
<b>cd_customnimfrqs</b> lets the user define the gridpoints along the imaginary axis by hand. Set this
to the number of frequencies you want. The frequencies are specified with <a href="vargw.html#cd_imfrqs">cd_imfrqs</a>.
"""
},
'cd_frqim_method': {
'definition': "Contour Deformation Imaginary Frequency integration Method",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a> """,
'vartype': "integer ",
'default': "1",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=4, that is, sigma (self-energy)
calculations, and <a href="vargw.html#gwcalctyp">gwcalctyp</a>=2, 12 or 22 (contour deformation
full-frequency treatment).
<p>
<b>cd_frqim_method</b> defines the choice of integration method along the imaginary frequency axis
for Contour Deformation calculations. The default method is very robust, fast and optimal for the
vast majority of cases. However, for very accurate (&ldquo;paranoid level&rdquo;) convergence studies,
ABINIT offers the possibility of a variety of methods and grids. Note that as one starts to change
the defaults, one needs to carefully consider the grid used. Therefore we recommend that in
addition to reading the infomation below, the user reads the description of the
input variables <a href="vargw.html#freqim_alpha">freqim_alpha</a>,
<a href="vargw.html#nfreqim">nfreqim</a>, <a href="vargw.html#ppmfrq">ppmfrq</a>,
<a href="vargw.html#gw_frqim_inzgrid">gw_frqim_inzgrid</a>.
<p>
The integration to be performed for each matrix element of the self energy along the imaginary axis is of the form:
</p>
<p align="center">
<img style="width: 461px; height: 49px;" src="./vargw_img/self_energy_cd.png">
</p>
<p>
Where <span style="font-family:Times,Serif"><i>&omega;</i></span> is the frequency point along the real axis,
<span style="font-family:Times,Serif"><i>&epsilon;<sub>s</sub></i></span> is an eigenvalue,
and <span style="font-family:Times,Serif"><i>i&omega;'</i></span> is the variable along the
imaginary axis. Thus the function to be integrated
is a Lorentzian weight function centred on the origin (whose FWHM is decided by
|<span style="font-family:Times,Serif"><i>&omega; - &epsilon;<sub>s</sub></i></span>|),
times a function. The function is related to the inverse dielectric matrix. It might have a peaked structure near the origin and is very smooth otherwise. the function decays
asymptotically as <span style="font-family:Times,Serif">1 / <i>i&omega;'</i></span>, so the whole
integral converges as this to the third power.
<ul>
<li><b>cd_frqim_method = 1 - Histogram:</b> This is the <b>default</b> method where the function
<span style="font-family:Times,Serif"><i>f(i&omega;')</i></span> is approximated by a histogram,
and the Lorentzian is integrated analytically in each sub-interval. See the section on grids
below for a description of the default grid. This method combined with the default grid is the
fastest and optimised for the use of few points along the imaginary axis.

<li><b>cd_frqim_method = 2 - Trapezoid:</b> The next step up from the histogram approximation in the previous
method. The integration region is transformed
<span style="font-family:Times,Serif"><i>[0,&#8734;[ &rarr; [0,1]</i></span> with a proper
weight depending on the width of the Lorentian. In this space
<span style="font-family:Times,Serif"><i>f(i&omega;')</i></span>
is approximated by a linear function between grid points (trapezoids), and the integrand is
integrated analytically in each sub-interval. This method tends to slightly overestimate contributions
while the default method tends to slightly underestimate them, so the results from methods
1 and 2 should bracket the converged values. The asymptotic behaviour is explicitly taken
into account by a fit using the last two grid points.

<li><b>cd_frqim_method = 3, 4, 5 - Natural Spline:</b> The function is transfomed <span style="font-family:Times,Serif"><i>[0,&#8734;[ &rarr; [0,1]</i></span>. In this space <span style="font-family:Times,Serif"><i>f(i&omega;')</i></span> is approximated by a natural spline function whose starting and ending sections are linear. This transform is chosen so that the function should approach a linear function asymptotically as the integration interval approaches 1, so that the asymptotic behaviour is automatically taken into account. For each Lorenzian width
(determined by |<span style="font-family:Times,Serif"><i>&omega; - &epsilon;<sub>s</sub></i></span>|)
the integrand is appropriately scaled in the interval
<span style="font-family:Times,Serif"><i>[0,1]</i></span>, and a nested Gauss-Kronrod (GK)
numerical integration rule is performed. The integrand is evaluated at the GK nodes by means
of a spline-fit. The order of the GK rule is controlled by the index of the method:
<ul>
<li><b>3 =&gt; Gauss  7 point, Kronrod 15 point rule </b>
<li><b>4 =&gt; Gauss 11 point, Kronrod 23 point rule </b>
<li><b>5 =&gt; Gauss 15 point, Kronrod 31 point rule </b>
</ul>
There is rarely any difference to machine precision between these rules, and the code will
issue a warning if a higher-order rule is recommended.

</ul>

<p>
<b>Grids for the integral along the imaginary axis:</b>
<p>
All the methods above should execute no matter what grid is used along the imaginary axis, so this
is very much under the user's control. The only requirement is that the grid be strictly increasing.
The point at zero frequency is assumed to lie on the real axis, so the calculation of that point
is controlled by <a href="vargw.html#nfreqim">nfreqre</a> and corresponding variables. We highly
recommend extracting various elements of the dielectric matrix from the _SCR file using
the <B>Mrgscr</B> utility and plotting them for visual inspection.
<p>
<ul>
<li><b>Default</b> - The default grid is an exponentially increasing grid given by the formula:
<p align="center">
<img style="width: 309px; height: 46px;" src="./vargw_img/cd_default_grid.png">
</p>Here <span style="font-family:Times,Serif"><i>&omega;<sub>p</sub></i></span> is the
plasma frequency (by default determined by the average density of the system, but this
can be overridden by setting <a href="vargw.html#ppmfrq">ppmfrq</a>).
<span style="font-family:Times,Serif"><i>N</i></span> is the total number of gridpoints
(set by <a href="vargw.html#nfreqim">nfreqim</a>).
<span style="font-family:Times,Serif"><i>&alpha;</i></span> is a parameter which determines
how far out the final grid point will lie. The final point will be at
<span style="font-family:Times,Serif"><i>&alpha;&middot;&omega;<sub>p<sub></i></span>
(the default is <span style="font-family:Times,Serif"><i>&alpha; = 5</i></span>,
and was hard-coded in older versions of ABINIT).
This grid is designed so that approximately half the grid points are always distributed
to values lower than the plasma frequency, in order to resolve any peaked structure. If one
seeks to increase the outermost reach by increasing <a href="vargw.html#ppmfrq">ppmfrq</a>
one must simultaneously take care to increase <a href="vargw.html#nfreqim">nfreqim</a>
in order to have the appropriate resolution for the low-frequency region.
In more recent versions of ABINIT one can also simply adjust the parameter
<span style="font-family:Times,Serif"><i>&alpha;</i></span> by using
<a href="vargw.html#freqim_alpha">freqim_alpha</a>. This grid is optimised for
speed and accurate results with few gridpoints for <b>cd_frqim_method = 1</b>.
<li><b>Inverse z transform</b> - This grid is activated by the use of the variable
<a href="vargw.html#gw_frqim_inzgrid">gw_frqim_inzgrid</a>. This is the standard
<span style="font-family:Times,Serif"><i>[0,&#8734;[ &rarr; [0,1]</i></span> transform
using the formula:
<p align="center">
<img style="width: 122px; height: 38px;" src="./vargw_img/cd_inzgrid.png">
</p>Here <span style="font-family:Times,Serif"><i>&omega;<sub>p</sub></i></span> is the
plasma frequency (default can be overridden by setting <a href="vargw.html#ppmfrq">ppmfrq</a>).
The gridpoints are then picked by an equidistant grid (number of points set by
<a href="vargw.html#nfreqim">nfreqim</a>) in the interval
<span style="font-family:Times,Serif"><i>z &sub; [0,1]</i></span>. This grid can easily
be uniquely converged by just increasing <a href="vargw.html#nfreqim">nfreqim</a>. Again
the points are distributed so that approximately half of them lie below the plasma frequency.
<li><b>User defined</b> - The user can also define their own grid using the variables
<a href="vargw.html#nfreqim">cd_customnimfrqs</a> and
<a href="vargw.html#nfreqim">cd_imfrqs</a>. <i>With great power comes great responsibility!</i>
</ul>
<p>
The <B>Mrgscr</B> utility is handy in optimising the numerical effort expended in convergence studies.
By estimating the densest grid one can afford to calculate in the SCR file, and succesively removing
frequencies from a single file (using the utility), one only needs to perform the screening
calculation <b>once</b> on the dense mesh for a given convergence study. One can also use the utility to
merge independent screening calculations over q-points and frequency sections.
"""
},
'cd_full_grid': {
'definition': "Contour Deformation Full Grid in complex plane",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a> """,
'vartype': "integer ",
'default': "0 (off)",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=3, that is, screening calculations,
and <a href="vargw.html#gwcalctyp">gwcalctyp</a>=2, 12 or 22.
<p>
<b>cd_full_grid</b> enables the calculation of the screening [both chi0 and epsilon^(-1)] on a grid in the first
quadrant of the complex plane. The grid is determined by the (tensor) product of the grid in real frequency and
the grid in imaginary frequency. In the SUS and SCR files the grid points are stored as follows:
<pre>
<b>Index:</b>  1   . . .   nfreqre   nfrqre+1 . . . nfreqre+nfreqim   nfreqre+nfreqim+1 . . . nfreqre*nfreqim
<b>Entry:</b> | purely real freq.  |     purely imaginary freq.     |      gridpoints in complex plane        |
</pre>
The grid in the complex plane is stored looping over the real dimension as the inner loop and the imaginary as
the outer loop. The contents of the generated SUS and SCR files can be extracted for visualisation and further
analysis with the <B>Mrgscr</B> utility.
"""
},
'cd_halfway_freq': {
'definition': "Contour Deformation tangent grid Halfway Frequency ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>, <a href="../users/abinit_help.html#energy">ENERGY</a> """,
'vartype': "real parameter",
'default': "100.0 eV",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=3 or 4, that is, screening or sigma calculations,
and <a href="vargw.html#gwcalctyp">gwcalctyp</a>=2, 12 or 22.
<p>
<b>cd_halfway_freq</b> determines the frequency where half of the number of points defined in
<a href="vargw.html#nfreqre">nfreqre</a> are used up. The tangent transformed grid is approximately
linear up to this point. To be used in conjunction with <a href="vargw.html#gw_frqre_tangrid">gw_frqre_tangrid</a>.
"""
},
'cd_imfrqs': {
'definition': "Contour Deformation Imaginary Frequencies ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a> """,
'vartype': "real array cd_imfrqs(cd_customnimfrqs)",
'default': "none",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=3, that is, screening calculations,
and <a href="vargw.html#gwcalctyp">gwcalctyp</a>=2, 12 or 22.
<p>
<b>cd_imfrqs</b> Specifies the grid points for the imginary axis. Only activated if
<a href="vargw.html#cd_customnimfrqs">cd_customnimfrqs</a> is not equal to 0. The number of frequencies
is set by the value of <a href="vargw.html#cd_customnimfrqs">cd_customnimfrqs</a>. Example:
<pre>
cd_customnimfrqs   5
nfreqim            5
cd_imfrqs          0.1  0.2  0.5  1.0  5.0
</pre>
If <a href="vargw.html#nfreqim">nfreqim</a> is not equal to
<a href="vargw.html#cd_customnimfrqs">cd_customnimfrqs</a> a warning will be issued.
</p>
<p>
<b>Use at own risk!</b> The use of a custom grid makes it your responsibility that the SUS and
SCR files are valid in self-energy (i.e. <a href="vargs.html#optdriver">optdriver</a>=4)
calculations, so caution is advised. Note that frequencies have to be strictly increasing, and the
point at zero frequency is <b>not</b> considered to be part of the imaginary grid, but rather
the grid along the real axis. The calculation of that point should be controlled by
<a href="vargw.html#nfreqim">nfreqre</a> and related variables.
</p><p>
"""
},
'cd_max_freq': {
'definition': "Contour Deformation grid Maximum Frequency ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>, <a href="../users/abinit_help.html#energy">ENERGY</a> """,
'vartype': "real parameter",
'default': "1000.0 eV",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=3 or 4, that is, screening or sigma calculations,
and <a href="vargw.html#gwcalctyp">gwcalctyp</a>=2, 12 or 22.
<p>
<b>cd_max_freq</b> determines the frequency where all the points defined in
<a href="vargw.html#nfreqre">nfreqre</a> are used up. To be used in conjunction with
<a href="vargw.html#gw_frqre_tangrid">gw_frqre_tangrid</a>.
"""
},
'cd_subset_freq': {
'definition': "Contour Deformation grid calculate Subset of Frequencies ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a> """,
'vartype': "integer array cd_subset_freq(2)",
'default': "full set",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=3, that is, screening calculations,
and <a href="vargw.html#gwcalctyp">gwcalctyp</a>=2, 12 or 22.
<p>
<b>cd_subset_freq</b> Specifies that only a subset of the frequencies defined by
<a href="vargw.html#nfreqre">nfreqre</a> are to be calculated. The first index is the start and
the second the end, with index number 1 always being the origin. For example a calculation with
<b><a href="vargw.html#nfreqre">nfreqre</a>=100</b> could be separated into two datasets with:
<pre>
subset_freq1   1   50
subset_freq2   51  100
</pre>
Any resulting susceptibiltiy (_SUS) and screening (_SCR) files can then be merged with the <b>mrgscr</b> utility.
</p><p>
<br>
Note that this does <b>NOT</b>
have to be used in conjunction with <a href="vargw.html#gw_frqre_tangrid">gw_frqre_tangrid</a>.
"""
},
'charge': {
'definition': "CHARGE ",
'section': "vargs",
'category': " ",
'vartype': "real number  ",
'default': "0.",
'text': """Used to establish charge balance between
the number of electrons filling the bands and the
nominal <b>charge</b> associated with the atomic cores.
<br>The code adds up the number of valence electrons
provided by the pseudopotentials of each type
(call this "zval"), then add <b>charge</b>, to get the
number of electrons per unit cell,
<a href="varint.html#nelect">nelect</a>.
<br>
Then, if <a href="varbas.html#iscf">iscf</a> is positive,
the code adds up the band occupancies (given in
array <a href="vargs.html#occ">occ</a>) for all bands at each k point,
then multiplies
by the k point weight <a href="varbas.html#wtk">wtk</a> at each k point.
Call this sum "nelect_occ" (for the number of electrons
from occupation numbers).  It is then
required that:
<br>nelect_occ = nelect
<br>
To treat a neutral
system, which is desired in nearly all cases, one must
use <b>charge</b>=0.  To treat a system missing one electron
per unit cell, set <b>charge</b>=+1.
"""
},
'chkexit': {
'definition': "CHecK whether the user want to EXIT ",
'section': "vargs",
'category': "",
'vartype': "integer parameter ",
'default': "0 (before v5.3, it was 2 for sequential version of ABINIT, 1 for parallel version of ABINIT.)",
'text': """If <b>chkexit</b> is 1 or 2, ABINIT
will check whether the user wants to interrupt the run (using the keyword
"exit" on the top of the input file or creating a file
named "abinit.exit": see the
<a href="../users/abinit_help.html#chkexit">end of section 3.2</a> of abinit_help).
<p>
If <b>chkexit</b>=0, the check is not performed at all
<p>
If <b>chkexit</b>=1, the check is not performed frequently (after each SCF step)
<p>
If <b>chkexit</b>=2, the check is performed frequently
(after a few bands, at each k point)
<p>In all cases, the check is performed at most every 2 seconds of CPU time.
"""
},
'chkgwcomp': {
'definition': "Check GW Completeness",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>""",
'vartype': "integer  ",
'default': "0",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=3, that is, screening calculations.
<p>
<b>chkgwcomp</b> performs a check on the completeness relation that is assumend to be fulfilled by the oscillator strenghts
in the limit of very large number of bands.
If <b>chkgwcomp==1</b>, a _DELI file is generated, containing the incompleteness matrix in reciprocal space for each q-point.
This matrix is the difference between the identity and the projection over all bands (the oscillator strengths),
summed over all occupied states. It is relevant in case one performs a GW + PAW calculation, especially if one is
willing to use the extrapolar approximation (<b><a href="vargw.html#gwencomp">gwencomp</a>=1</b>).
Note that this closure relation is not fulfilled in PAW, due to the incompleteness of the basis.
"""
},
'chkprim': {
'definition': "CHecK whether the cell is PRIMitive ",
'section': "vargs",
'category': """<a href="../users/abinit_help.html#symmetry_finder">SYMMETRY FINDER</a> """,
'vartype': "integer parameter ",
'default': "1.",
'text': """If the symmetry finder is used
(see <a href="varbas.html#nsym">nsym</a>), a non-zero
value of <b>chkprim</b> will make the code stop if a non-primitive
cell is used. If <b>chkprim</b>=0, a warning is issued, but the run
does not stop.
<p>If you are generating the atomic and cell geometry using
<a href="vargeo.html#spgroup">spgroup</a>, you might
generate a PRIMITIVE cell using
<a href="vargeo.html#brvltt">brvltt</a>=-1 .
"""
},
'chksymbreak': {
'definition': "CHecK SYMmetry BREAKing  ",
'section': "vargs",
'category': """<a href="../users/abinit_help.html#not_internal">NOT INTERNAL</a> """,
'vartype': "integer parameter ",
'default': "1.",
'text': """This variable governs the behaviour of the code when there are potential
source of symmetry breaking, related e.g. to the k point grid or the presence of
non-symmorphic translations which might not be coherent with the exchange-correlation grid.
<p>
When <b>chksymbreak</b>=1, the code stops if :
<ul>
<li>(1) The k point grid is non-symmetric, in case <a href="varbas.html#nsym">kptopt</a>=1, 2, or 4 ;</li>
<li>(2) The non-symmorphic translation part of the symmetry operations has components that are not zero,
or simple fractions, with 2, 3, 4, 6, 8 or 12 as denominators.</li>
</ul>
<br>When <b>chksymbreak</b> is zero, there is no such check.
<br>When <b>chksymbreak</b> is minus 1, the code stops if the condition (1) is met, but in case the condition (2) is met, there will be a trial
to shift the atomic coordinates such as to obtain symmetry operations with the adequate non-symmorphic part.

<p>Explanation :
<br>In the ground-state calculation, such breaking of the symmetry is usually harmless. However, if the user is doing a
calculation of phonons using DFPT (<a href="varrf.html#rfphon">rfphon</a>=1), the convergence with respect to the number of k points will
be much worse with a non-symmetric grid than with a symmetric one. Also, if the user is doing a GW calculation, the
presence of non-symmorphic translations that are not coherent with the FFT grid might cause
problems.
In the GW part, indeed, one needs to reconstruct the wavefunctions in the full Brillouin zone for
calculating both the polarizability and the self-energy.
The wavefunctions in the full Brillouin zone are obtained from the irreducible wedge by applying the symmetry
operations of the space group of the crystal.
In the present implementation, the symmetrization of the wavefunctions is done in real space on the FFT mesh
that, therefore, has to be coherent both with the rotational part as well as with the fractional translation
of each symmetry operation.
If the condition (2) is met, the GW code won't be able to find a symmetry-preserving FFT mesh.
<br>So, it was decided to warn the user about these possible problems already at the level of the ground state calculations,
although such warning might be irrelevant.
<br>If you encounter a problem outlined above, you have two choices : change your atomic positions (translate them) such that the origin
appears as the most symmetric point ; or ignore the problem, and set <b>chksymbreak</b>=0 .
"""
},
'cineb_start': {
'definition': "Climbing-Image Nudged Elastic Band: STARTing iteration",
'section': "varrlx",
'category': "",
'vartype': "integer parameter",
'default': "<b>cineb_start</b>=7",
'text': """Relevant only when <a href="varrlx.html#imgmov">imgmov</a>=5 (Nudged Elastic Band) and
<a href="varrlx.html#neb_algo">neb_algo</a>=2 (CI-NEB).<br>
Gives the index of the first CI-NEB iteration..<br>
The CI-NEB method constitutes a small modification to the NEB method allowing a rigorous
convergence to the saddle point. As the image with the highest energy has to be identified,
the calculation begins with several iterations of the standard NEB algorithm.
The effective CI-NEB begins at the <b>cineb_start</b> iteration.<br>
<i>See: J. Chem. Phys. 113, 9901 (2000).</i>
"""
},
'cmlfile': {
'definition': "Chemical Markup Language FILE ",
'section': "varfil",
'category': """<a href="../users/abinit_help.html#not_internal">NOT INTERNAL</a> """,
'vartype': "character string ",
'default': "no file.",
'text': """Used to import some of the data from one or more Chemical Markup Language 2 (CML2) file(s)
(one per dataset). Unlike most of the other input variables, it refers
to a character string, e.g. :
<br>
<pre>
cmlfile ../t67.in_CML.xml
</pre>
The file is preprocessed, and the relevant information is translated in order
to be used as an alternative to the usual input variables. Note that the input variables
directly defined in the usual input file have precedence over the CML data :
the latter are used only when there is no occurence of the corresponding
keyword in the input file...
<br> The ABINIT CML parser is still quite primitive. The mechanism
followed to parse the CML file is described afterwards.
<br> The ABINIT CML parser will localize in the CML file the first occurence
of a 'molecule' markup section. It will ignore all other occurences of
'molecule'. Inside this 'molecule' section, it will localize the first occurences
of the 'crystal', 'symmetry' and 'atomArray' sections. All other occurences,
and all other sections, are ignored.
<br> The following ABINIT input variables will be extracted from these
sections of the CML file (if the data is available) :
<ul>
<li> <a href="varbas.html#acell">acell</a> from the first 'scalar title="a"','scalar title="b"',
and 'scalar="c"' sections (all three must be present if one is present)
in the 'crystal' section, expecting the
data in Angstrom;</li>
<li> <a href="varbas.html#angdeg">angdeg</a> from the first 'scalar title="alpha"','scalar title="beta"',
and 'scalar title="gamma"' sections (all three must be present if one is present)
in the 'crystal' section, expecting the
data in degrees; </li>
<li> <a href="varbas.html#nsym">nsym</a>,
<a href="varbas.html#symrel">symrel</a> and <a href="varbas.html#tnons">tnons</a> from
the content of 'matrix' sections in the 'symmetry' section;</li>
<li> <a href="varbas.html#natom">natom</a> from the number of items in the first
'atom' sections in the 'atomArray' section;</li>
<li> <a href="varbas.html#ntypat">ntypat</a> from the number of types of atoms over all input CML files;</li>
<li> <a href="varbas.html#typat">typat</a> from the attribute 'elementType' in the
'atom' sections in the 'atomArray' section, with identification
of the pseudopotentials that have the correct nuclear charge, according to the atomic symbol
(the first pseudopotential with the correct nuclear charge, from the pseudopotential list, will be used);</li>
<li> Coordinates of the atoms :
<ul>
<li> <a href="varbas.html#xred">xred</a> from the attributes 'xFract', 'yFract', and 'zFract'
(all three must be present if one is present) in the 'atom' sections
in the 'atomArray' section.</li>
<ul><B>OR</B></ul>
<li> <a href="varbas.html#xcart">xcart</a> from the attributes 'x3', 'y3', and 'z3'
(all three must be present if one is present) in the 'atom' sections
in the 'atomArray' section.</li>
<li> <a href="varbas.html#znucl">znucl</a> for all of the CML files read in.</li>
</ul>
</ul>
<br> These limited parsing capabilities are enough for ABINIT to read the CML files it has created
thanks to the use of the <a href="varfil.html#prtcml">prtcml</a> input variable. Mixing CML files
with normal input or <a href="vargeo.html#xyzfile">xyzfile</a> will not work.
"""
},
'cpu': {
'definition': "CPU time limit in: Seconds, Minutes, Hours   ",
'section': "vargs",
'category': """<a href="../users/abinit_help.html#no_multi">NO MULTI</a> ; for cpum and cpuh : <a href="../users/abinit_help.html#not_internal">NOT INTERNAL</a>""",
'vartype': "real parameters  ",
'default': "0.0d0.",
'text': """One of these three real parameters can be
defined in the input file, to set up a CPU time limit.
When the job reaches that limit, it will try to end smoothly.
However, note that this might still take some time.
If the user want a firm CPU time limit, the present
parameter must be reduced sufficiently. Intuition
about the actual margin to be taken into account
should come with experience ...
<br>Note that only one of these three parameters can be defined
in a single input file.
A zero value has no action of the job.
<br>Internally, only <b>cpus</b> is used in the dtset array: adequate
conversion factors are used to generate it from <b>cpum</b> or
<b>cpuh</b>."""
},
'ddamp': {
'definition': "electric Displacement field DAMPing parameter ",
'section': "varff",
'category': " ",
'vartype': "real parameter  ",
'default': "0.1 .",
'text': """In case <a href="varff.html#berryopt">berryopt</a>=6 or 7,
the electric field is updated after each SCF iteration according to
E_{n+1}= <b>ddamp</b>*(D - 4*pi*P_{n}) + (1-<b>ddamp</b>)*E_{n}
where P_{n} and E_{n} are the polarization and electric field after nth SCF iteration.
<b>ddamp</b> is a damping parameter used to control the convergence speed.
"""
},
'delayperm': {
'definition': "DELAY between trials to PERMUTE atoms ",
'section': "varrlx",
'category': "",
'vartype': "integer ",
'default': "0.",
'text': """Delay (number of time steps) between trials to permute
two atoms, in
view of accelerated search of minima. Still in development. See the
routine moldyn.F90. See also <a href="varrlx.html#signperm">signperm</a>.
When <b>delayperm</b> is zero, there is not permutation trials.
"""
},
'densty': {
'definition': "initial DENSity for each TYpe of atom ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a> """,
'vartype': """real array densty(<a href="varbas.html#ntypat">ntypat</a>) """,
'default': "0.0d0.",
'text': """Gives a rough description
of the initial GS density, for each type of atom.
This value is only used to create
the first exchange and correlation potential,
and is not used anymore afterwards.
For the time being, it corresponds to an average
radius (a.u.) of the density, and is used to generate
a gaussian density. If set to 0.0d0, an optimized value is used.
<br>No meaning for RF calculations."""
},
'dfield': {
'definition': "Displacement FIELD ",
'section': "varff",
'category': " ",
'vartype': "real array dfield(3)  ",
'default': "3*0.0 .",
'text': """In case <a href="varff.html#berryopt">berryopt</a>=6 or 7,
a unreduced finite electric displacement field calculation is performed. The value
of this displacement field, and its direction is determined by <b>dfield</b>.
It must be given in atomic units.
"""
},
'diecut': {
'definition': "DIElectric matrix Energy CUToff ",
'section': "vargs",
'category': """<a href="../users/abinit_help.html#energy">ENERGY</a> """,
'vartype': "real parameter ",
'default': "2.2d0 Ha.",
'text': """Kinetic energy cutoff that controls the number
of planewaves used to represent the dielectric matrix:
<br>(1/2)[(2 Pi)*(Gmax)]<sup>2</sup>=<a href="varbas.html#ecut">ecut</a> for Gmax.
<br>Can be specified in Ha (the default), Ry, eV or Kelvin, since
<b>diecut</b> has the
'<a href="../users/abinit_help.html#dimensions">ENERGY</a>' characteristics.
(1 Ha=27.2113845 eV)
<br>All planewaves inside this "basis sphere" centered
at G=0 are included in the basis.
This is useful only when <a href="vargs.html#iprcel">iprcel</a>>=21, which means that
a preconditioning scheme based on the dielectric matrix
is used.
<br>NOTE : a negative <b>diecut</b> will define the same dielectric
basis sphere as the corresponding positive value,
but the FFT grid will be identical to the one used
for the wavefunctions.
The much smaller FFT grid, used when <b>diecut</b> is positive,
gives exactly the same results.
<br>No meaning for RF calculations yet."""
},
'diegap': {
'definition': "DIElectric matrix GAP ",
'section': "vargs",
'category': """<a href="../users/abinit_help.html#energy">ENERGY</a> """,
'vartype': "real parameter ",
'default': "0.1 Ha.",
'text': """Gives a rough estimation of the dielectric gap
between the highest energy level computed in the run,
and the set of bands not represented.
Used to extrapolate dielectric matrix when <a href="vargs.html#iprcel">iprcel</a> >= 21.
<br>Can be specified in Ha (the default), Ry, eV or Kelvin, since
<b>diegap</b> has the
'<a href="../users/abinit_help.html#dimensions">ENERGY</a>' characteristics.
(1 Ha=27.2113845 eV)
<br>No meaning for RF calculations yet."""
},
'dielam': {
'definition': "DIElectric matrix LAMbda ",
'section': "vargs",
'category': "",
'vartype': "real parameter between 0 and 1 ",
'default': "0.5 .",
'text': """Gives the amount of occupied states with mean energy given by the
highest level computed in the run, included
in the extrapolation of the dielectric matrix.
Used when <a href="vargs.html#iprcel">iprcel</a> >= 21.
<br>No meaning for RF calculations yet."""
},
'dielng': {
'definition': "model DIElectric screening LeNGth ",
'section': "vargs",
'category': " ",
'vartype': "real parameter ",
'default': "1.0774841d0 (Bohr), for historical reasons.",
'text': """Used for screening length (in Bohr) of the model
dielectric function, diagonal in reciprocal space.
By default, given in Bohr atomic units
(1 Bohr=0.5291772108 Angstrom), although Angstrom can be specified,
if preferred, since <b>dielng</b> has the
'<a href="../users/abinit_help.html#dimensions">LENGTH</a>' characteristics.
<br>
This model dielectric function is as follows :
<pre>
(     1        + <b>dielng</b><sup>2</sup> * K<sup>2</sup> )
diel(K)= --------------------------------------------
( 1/<a href="vargs.html#diemac">diemac</a> + <b>dielng</b><sup>2</sup> * K<sup>2</sup> ) * <a href="vargs.html#diemix">diemix</a>
</pre>
The inverse of this model dielectric function will be
applied to the residual, to give the preconditioned
change of potential. Right at K=0, diel(K) is imposed to be 1.
<p>If the preconditioning were perfect,
the change of potential would lead to an exceedingly fast solution
of the self-consistency problem (two or three steps).
The present model dielectric function is excellent for
rather homogeneous unit cells.
<br>When K->0 , it tends to the macroscopic dielectric
constant, eventually divided by the mixing factor <a href="vargs.html#diemix">diemix</a>
(or <a href="vargs.html#diemix">diemixmag</a> for magnetization).
<br>For metals, simply put <a href="vargs.html#diemac">diemac</a> to a very large value (10^6 is OK)
<br>The screening length <b>dielng</b> governs the length scale
to go from the macroscopic regime to the microscopic
regime, where it is known that the dielectric function
should tend to 1. It is on the order of 1 Bohr for
metals with medium density of states at the Fermi level,
like Molybdenum, and for Silicon. For metals with a
larger DOS at the Fermi level (like Iron),
the screening will be more effective, so that <b>dielng</b>
has to be decreased by a factor of 2-4.
<br>This works for GS and RF calculations.

TO BE IMPROVED : it is not clear what the variable K is and how to deal it, David Waroquiers, 090831.
"""
},
'diemac': {
'definition': "model DIElectric MACroscopic constant ",
'section': "vargs",
'category': " ",
'vartype': "real parameter ",
'default': "10<sup>6</sup> (metallic damping). ",
'text': """A rough knowledge of the macroscopic dielectric constant <b>diemac</b>
of the system is a useful help to speed-up the SCF procedure:
a model dielectric function,
see the keyword <a href="vargs.html#dielng">dielng</a>, is used for that
purpose.  It is especially
useful for speeding up the treatment of rather homogeneous unit cells.
<p>Some hint :
<br>The value of <b>diemac</b> should usually be bigger than 1.0d0,
on physical grounds.
<br>For metals, simply put <b>diemac</b> to a very large value (the default 10<sup>6</sup> is OK)
<br>For silicon, use 12.0 . A similar value is likely to work well for
other semiconductors
<br>For wider gap insulators, use 2.0 ... 4.0
<br>For molecules in an otherwise empty big box, try 1.5 ... 3.0
<br>Systems that combine a highly polarisable part and some vacuum are rather
badly treated by the model dielectric function. One has to use the
"extrapolar" technique, activated by the  input variable
<a href="vargs.html#iprcel">iprcel</a>.
<br>In sufficiently homogeneous systems, you might have to experiment
a bit to find the best <b>diemac</b>. If you let <b>diemac</b>
to its default value, you might even never obtain the self-consistent convergence !
<br>For response function calculations, use the same
values as for GS. The improvement in speed can be considerable
for small (but non-zero) values of the wavevector.
"""
},
'diemix': {
'definition': "model DIElectric MIXing factor ",
'section': "vargs",
'category': " ",
'vartype': "real parameter ",
'default': """1.0 (norm-conserving psps or <a href="vargs.html#iprcel">iprcel</a>/=0) or 0.7 (PAW and <a href="vargs.html#iprcel">iprcel</a>=0).""",
'text': """Gives overall factor of the preconditioned
residual density/potential to be transferred in the SCF cycle.
<br>It should be between 0.0 and 1.0 .
<br>If the model dielectric function were perfect, <b>diemix</b>
should be 1.0 . By contrast, if the model dielectric function
does nothing (when <a href="vargs.html#diemac">diemac</a>=1.0d0 or <a href="vargs.html#dielng">dielng</a>
is larger than the
size of the cell), <b>diemix</b> can be used
to damp the amplifying factor inherent to the SCF loop.
<br>For molecules, a value on the order 0.5 or 0.33 is rather usual.
<br>When mod(<a href="varbas.html#iscf">iscf</a>,10)=3, 4 ,5 or 7, <b>diemix</b>
is only important at the
few first iterations when anharmonic effects are important,
since these schemes compute their own mixing factor
for self-consistency.<br>
Also note that a different value of diemix can be used for the magnetization
(see <a href="vargs.html#diemixmag">diemixmag</a>).
"""
},
'diemixmag': {
'definition': "model DIElectric MIXing factor for the MAGgnetization",
'section': "vargs",
'category': " ",
'vartype': "real parameter ",
'default': """diemixmag= <a href="vargs.html#diemixmag">diemix</a> when <a href="vargs.html#iprcel">iprcel</a>=0 or 70<<a href="vargs.html#iprcel">iprcel</a><80 or <a href="varbas.html#iscf">iscf</a><10 (SCF mixing on potential), else diemixmag= - <a href="vargs.html#diemixmag">diemix</a>""",
'text': """Gives overall factor of the preconditioned
residual magnetization/magnetic field to be transferred in the SCF cycle (see <a href="vargs.html#diemixmag">diemix</a> for further information).<br>
For the time being, apply only when the SCF mixing is done on the density
(<a href="varbas.html#iscf">iscf</a>>=10).<br><br>
A negative value of diemixmag means that magnetization is only preconditionned by ABS(diemixmag),
without the use of any preconditionner.<br><br>
When SCF cycle has some difficulties to converge, changing the value of <b>diemixmag</b>
can have a positive effect.<br>In particular <b>diemixmag</b>=-4 is a good choice
(i.e. diemixmag=4, no other preconditionner on magnetization).
"""
},
'diismemory': {
'definition': "Direct Inversion in the Iterative Subspace MEMORY ",
'section': "varrlx",
'category': "",
'vartype': "integer ",
'default': "8.",
'text': """Gives the maximum number of "time" steps for which the
forces and stresses are stored, and taken into account in the
DIIS algorithm (<a href="varrlx.html#ionmov">ionmov</a>=20)
to find zero-force and stress configurations.
"""
},
'dilatmx': {
'definition': "DILATation : MaXimal value ",
'section': "varrlx",
'category': "",
'vartype': "real parameter ",
'default': "1.0 .",
'text': """Gives the maximal permitted scaling of
the lattice parameters when the cell shape and
dimension is varied (see variable <a href="varrlx.html#optcell">optcell</a>).
It is used to define the sphere of plane waves
and FFT box coherent with the possible modifications
of the cell (<a href="varrlx.html#ionmov">ionmov</a>==2 and <a
href="varrlx.html#optcell">optcell</a>/=0).
For these definitions, it is equivalent
to changing <a href="varbas.html#ecut">ecut</a> by multiplying it by <b>dilatmx</b><sup>2</sup>
(the result is an "effective ecut", called internally "ecut_eff",
other uses of <a href="varbas.html#ecut">ecut</a> being not modified
when <b>dilatmx</b>&gt;1.0 .
<br>
Using <b>dilatmx</b>&lt;1.0 is equivalent to changing <a
href="varbas.html#ecut">ecut</a>
in all its uses. This is allowed, although its meaning
is no longer related to a maximal expected scaling.
<br>
Setting <b>dilatmx</b> to a large value leads to waste
of CPU time and memory. Supposing you think that the
optimized <a href="varbas.html#acell">acell</a> values might be 10%
larger than your
input values, use simply <b>dilatmx</b> 1.1 . This will already
lead to an increase of the number of planewaves by a factor
(1.1)<sup>3</sup>=1.331 , and a corresponding increase in CPU time
and memory.
<br>
It is possible to use <b>dilatmx</b> when <a
href="varrlx.html#optcell">optcell</a>=0, but
a value larger than 1.0 will be a waste.
"""
},
'dmatpawu': {
'definition': "initial Density MATrix for PAW+U ",
'section': "varpaw",
'category': "",
'vartype': """real array dmatpawu(2*max(<a href="varpaw.html#lpawu">lpawu</a>)+1,2*max(<a href="varpaw.html#lpawu">lpawu</a>)+1, max(<a href="varbas.html#nsppol">nsppol</a>,<a href="vargs.html#nspinor">nspinor</a>), %<a href="varint.html#natpawu">natpawu</a>) <br>where %<a href="varint.html#natpawu">natpawu</a> is the number of atoms on which LDA/GGA+U is applied. """,
'default': "-10. (not defined)",
'text': """Relevant only when
<a href="varint.html#usepaw">usepaw</a>=1,
<a href="varpaw.html#usepawu">usepawu</a>=1,
<a href="varpaw.html#usedmatpu">usedmatpu</a>/=0 for Ground-State calculations.<br>
Gives the value of an initial density matrix used in LDA+U and kept
fixed during the first abs(<a href="varpaw.html#usedmatpu">usedmatpu</a>) SCF iterations.
<br>Only components corresponding to <a href="varpaw.html#lpawu">lpawu</a> angular momentum are requested.
<br>Restriction: In order to use dmatpawu, <a href="varpaw.html#lpawu">lpawu</a> must be identical for all atom types (or -1).
<br>The occupation matrix is in the basis of real spherical harmonics Slm (note that this differs from the choice made for
<a href="varfil.html#prtdosm">prtdosm</a>, that is in the basis of complex spherical harmonics).
Their are ordered by increasing m, and are defined e.g. in
the article
"Evaluation of the rotation matrices in the basis of real spherical harmonics",
by Miguel A. Blancoa, M. Floreza, M. Bermejo,
Journal of Molecular Structure (Theochem) 419, 19 (1997), that can be downloaded from
<a href="http://azufre.quimica.uniovi.es/articles/Theochem419-19-ov-BF97-rotation-matrices.pdf">
<!-- <a href="http://www.unioviedo.es/qcg/art/Theochem419-19-ov-BF97-rotation-matrices.pdf"> -->
the author Web site</a>.
For the case l=2 (d states), the five columns corresponds respectively to (the normalisation factor has been dropped)
<br>
<ul>
<li>m=-2, xy
<li>m=-1, yz
<li>m=0, 3z^2-r^2
<li>m=1, xz
<li>m=2, x^2-y^2
</ul>
<br>
<b>dmatpawu</b> must always be given as a "spin-up" occupation matrix (and eventually a "spin-down" matrix).
Be aware that its physical meaning depends on the magnetic properties imposed to the system
(with <a href="varbas.html#nsppol">nsppol</a>,
<a href="vargs.html#nspinor">nspinor</a>,
<a href="vargs.html#nspden">nspden</a>):<br>
<ul>
<li><b>Non-magnetic system</b>
(<a href="varbas.html#nsppol">nsppol</a>=1,
<a href="vargs.html#nspinor">nspinor</a>=1,
<a href="vargs.html#nspden">nspden</a>=1):<br>
One (2lpawu+1)x(2lpawu+1) <b>dmatpawu</b> matrix is given for each atom on which +U is applied.<br>
It contains the "spin-up" occupations.</li>
<li><b>Ferromagnetic spin-polarized (collinear) system</b>
(<a href="varbas.html#nsppol">nsppol</a>=2,
<a href="vargs.html#nspinor">nspinor</a>=1,
<a href="vargs.html#nspden">nspden</a>=2):<br>
Two (2lpawu+1)x(2lpawu+1) <b>dmatpawu</b> matrices are given for each atom on which +U is applied.<br>
They contain the "spin-up" and "spin-down" occupations.</li>
<li><b>Anti-ferromagnetic spin-polarized (collinear) system</b>
(<a href="varbas.html#nsppol">nsppol</a>=1,
<a href="vargs.html#nspinor">nspinor</a>=1,
<a href="vargs.html#nspden">nspden</a>=2):<br>
One (2lpawu+1)x(2lpawu+1) <b>dmatpawu</b> matrix is given for each atom on which +U is applied.<br>
It contains the "spin-up" occupations.</li>
<li><b>Non-collinear magnetic system</b>
(<a href="varbas.html#nsppol">nsppol</a>=1,
<a href="vargs.html#nspinor">nspinor</a>=2,
<a href="vargs.html#nspden">nspden</a>=4):<br>
Two (2lpawu+1)x(2lpawu+1) <b>dmatpawu</b> matrices are given for each atom on which +U is applied.<br>
They contains the "spin-up" and "spin-down" occupations
(defined as n_up=(n+|m|)/2 and n_dn=(n-|m|)/2), where m is the integrated magnetization vector).<br>
The direction of the magnetization (which is also the direction of n_up and n_dn) is given by
<a href="vargs.html#spinat">spinat</a>.<br>
<i>Warning: unlike collinear case, atoms having the same magnetization magnitude with different directions
must be given the same occupation matrix;<br>
the magnetization will be oriented by the value of <a href="vargs.html#spinat">spinat</a>
(this is the case for antiferro-magnetism).</i></li>
<li><b>Non-collinear magnetic system with zero magnetization</b>
(<a href="varbas.html#nsppol">nsppol</a>=1,
<a href="vargs.html#nspinor">nspinor</a>=2,
<a href="vargs.html#nspden">nspden</a>=1):<br>
Two (2lpawu+1)x(2lpawu+1) <b>dmatpawu</b> matrices are given for each atom on which +U is applied.<br>
They contain the "spin-up" and "spin-down" occupations;<br>
But, as "spin-up" and "spin-down" are constrained identical, the "spin-down" one is ignored by the code.</li>
</ul>
"""
},
'dmatpuopt': {
'definition': "Density MATrix for PAW+U OPTion",
'section': "varpaw",
'category': "",
'vartype': "integer parameter",
'default': "2",
'text': """Relevant only when
<a href="varint.html#usepaw">usepaw</a>=1 and
<a href="varpaw.html#usepawu">usepawu</a>=1<br>
This option governs the way occupations of localized atomic levels are computed:<br>
<ul>
<li> <b>dmatpuopt</b>=1: atomic occupations are projections on atomic orbitals (Eq. (6) of PRB 77, 155104 (2008)). <br>
</li>
<li> <b>dmatpuopt</b>=2: atomic occupations are integrated values in PAW spheres
of angular-momentum-decomposed charge densities (Eq. (7) of PRB 77, 155104 (2008)). <br>
</li>
<li> <b>dmatpuopt</b>=3: only for tests<br>
</li>
</li>
<li> <b>dmatpuopt</b>=4: Extrapolations of occupancies outside the PAW-sphere. This Definition gives normalized operator for occupation.<br>
</li>
</ul>
In the general case <b>dmatpuopt</b>=2 is suitable. The use of <b>dmatpuopt</b>=1
is restricted to PAW datasets in which the first
atomic wavefunction of the correlated subspace is a normalized atomic eigenfunction.
"""
},
'dmatudiag': {
'definition': "Density MATrix for paw+U, DIAGonalization ",
'section': "varpaw",
'category': "",
'vartype': "real parameter",
'default': "0.4943",
'text': """Relevant only when
<a href="varint.html#usepaw">usepaw</a>=1 and
(<a href="varpaw.html#usepawu">usepawu</a>=1  or
<a href="vardev.html#usedmft">usedmft</a>=1) and
<a href="varpaw.html#lpawu">lpawu</a>=3) .<br>
This gives the ratio of Slater Integrals F6 and F2.
It is used with <a href="varpaw.html#f4of2_sla">f4of2_sla</a>=3) .<br>
in DFT+U or DFT+DMFT for the calculation of the orbital
dependent screened coulomb interaction.
"""
},
'dmft_dc': {
'definition': "Dynamical Mean Fied Theory: Double Counting",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>  """,
'vartype': "integer  ",
'default': "1 ",
'text': """<br> Value of double counting used for DMFT. Only value 1 is activated for the moment and is the FLL double counting.
"""
},
'dmft_iter': {
'definition': "Dynamical Mean Fied Theory: number of ITERation ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>  """,
'vartype': "integer  ",
'default': "0 ",
'text': """<br> Number of iterations for the DMFT inner loop.
"""
},
'dmft_mxsf': {
'definition': "Dynamical Mean Fied Theory: MiXing parameter for the SelF energy ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>  """,
'vartype': "real  ",
'default': "0.3 ",
'text': """<br> Mixing parameter for the simple mixing of the self-energy.
"""
},
'dmft_nwli': {
'definition': "Dynamical Mean Fied Theory: Number of frequency omega (W) in the LInear mesh ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>  """,
'vartype': "integer  ",
'default': "0  ",
'text': """<br> Number of Matsubara frequencies (linear  mesh)
"""
},
'dmft_nwlo': {
'definition': "Dynamical Mean Fied Theory: Number of frequency omega (W) in the log mesh ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>  """,
'vartype': "integer  ",
'default': "0  ",
'text': """<br> Number of frequencies in the log mesh.
"""
},
'dmft_read_occnd': {
'definition': "Dynamical Mean Fied Theory: Read Occupations (Non Diagonal) ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>  """,
'vartype': "integer  ",
'default': "0  ",
'text': """<br> Flag to read/write Occupations as computed in DMFT. This flag is useful
to restart a DFT+DMFT calculation with self-consistency over electronic density.
The occupations are written each time a DMFT loop is finished. So if the calculations stops
because the time limit is reached, this option offers the possibility to restart the self-consistent loop
over density at the point where it stopped.
<ul>
<li> 0=&gt;  Occupations are written but never read.
<li> 1=&gt;  Occupations are read from I_DMFTOCCND, where I is the root for input files.
<li> 2=&gt;  Occupations are read from O_DMFTOCCND, where O is the root for output files.
</ul>

"""
},
'dmft_rslf': {
'definition': "Dynamical Mean Fied Theory: Read SeLF energy ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>  """,
'vartype': "integer  ",
'default': "0  ",
'text': """<br> Flag to read/write Self-Energy. If put to one, self-energy is written and read at each LDA iteration.
"""
},
'dmft_solv': {
'definition': "Dynamical Mean Fied Theory: choice of SOLVer ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>  """,
'vartype': "real  ",
'default': "0  ",
'text': """<br> Choice of solver for the Impurity model.
<ul>
<li> 1=&gt; LDA+U self-energy is used (for testing purpose)
<li> 2=&gt; Hubbard one solver.
</ul>
<br> WARNING: Quantum Monte Carlo (QMC) solvers are not yet interfaced with the code. The present version of LDA+DMFT implemented here
should NOT be used for correlated metals. Even of correlated (Mott) insulators, QMC is expected to be much more
precise: Hubbard one is a approximation !
"""
},
'dmft_tollc': {
'definition': "Dynamical Mean Fied Theory: Tolerance on Local Charge for convergency of the DMFT loop ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>  """,
'vartype': "real  ",
'default': "0.00001  ",
'text': """<br> Tolerance for the variation of Local Charge during iterations of the DMFT Loop.
<br> The default value is good for fast calculations. However, to obtain good convergency of the DFT Loop,
the DMFT Loop needs a better convergence criterion.
"""
},
'dmftbandf': {
'definition': "(to be described) ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>  """,
'vartype': "(to be described)  ",
'default': "0 ",
'text': """<br> dmftbandi and dmftbandf are the first and last bands taken into account in the Projected Local
Orbitals scheme of LDA+DMFT. They thus define the energy window used to define Wannier Functions.
(see  Amadon, B., Lechermann, F., Georges, A., Jollet, F., Wehling, T. O., and Lichtenstein, A. I. Phys. Rev. B 77(20), (2008).)

"""
},
'dmftbandi': {
'definition': "(to be described) ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>  """,
'vartype': "(to be described)  ",
'default': "0 ",
'text': """<br> dmftbandi and dmftbandf are the first and last bands taken into account in the Projected Local
Orbitals scheme of LDA+DMFT. They thus define the energy window used to define Wannier Functions.
(see  Amadon, B., Lechermann, F., Georges, A., Jollet, F., Wehling, T. O., and Lichtenstein, A. I. Phys. Rev. B 77(20), (2008).)

"""
},
'dmftcheck': {
'definition': "Dynamical Mean Fied Theory: CHECKs ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>  """,
'vartype': "integer  ",
'default': "0 ",
'text': """<br> (Introduced by B. Amadon, v6.1.0)
"""
},
'dosdeltae': {
'definition': "DOS Delta in Energy ",
'section': "vargs",
'category': """<a href="../users/abinit_help.html#energy">ENERGY</a> """,
'vartype': "real parameter ",
'default': "0.0 .",
'text': """Defines the linear grid resolution (energy increment) to be used for the
computation of the Density-Of-States, when <a href="varfil.html#prtdos">prtdos</a>
is non-zero.
<br>If <b>dosdeltae</b> is set to zero (the default value), the actual
increment is 0.001 Ha if <a href="varfil.html#prtdos">prtdos</a>=1, and
the much smaller value 0.00005 Ha if <a href="varfil.html#prtdos">prtdos</a>=2.
This different default value arises because the <a href="varfil.html#prtdos">prtdos</a>=1 case,
based on a smearing technique, gives a quite smooth DOS, while the DOS from the
tetrahedron method, <a href="varfil.html#prtdos">prtdos</a>=2, is rapidly varying.
"""
},
'dtion': {
'definition': "Delta Time for IONs ",
'section': "varrlx",
'category': "",
'vartype': "real parameter ",
'default': "100.",
'text': """Used for controlling ion time steps.
If <a href="varrlx.html#ionmov">ionmov</a> is set to 1, 6 or 7, then
molecular dynamics is&nbsp;
used to update atomic positions in response to
forces. The parameter <b>dtion</b> is a time step in
atomic units of time. (One atomic time unit is
2.418884e-17 seconds, which is the value of
Planck's constant in hartree*sec.)
In this case the atomic masses, in amu (given in array "<a
href="varrlx.html#amu">amu</a>"),
are used in Newton's equation and the viscosity (for <a
href="varrlx.html#ionmov">ionmov</a>=1)
and number of time steps are provided to the code using input
variables "<a href="varrlx.html#vis">vis</a>" and "<a
href="varrlx.html#ntime">ntime</a>".
The code actually converts
from masses in amu to masses in atomic units (in units
of electron masses) but the user enters masses in <a
href="varrlx.html#amu">amu</a>.
(The conversion from amu to atomic units (electron
masses) is 1822.88851 electron masses/amu.)
<br>
A typical good value for <b>dtion</b> is about 100.
The user must try several values
for <b>dtion</b> in order to establish the stable and efficient
choice for the accompanying amu, atom types and positions,
and <a href="varrlx.html#vis">vis</a> (viscosity).
<br>
For quenched dynamics (<a href="varrlx.html#ionmov">ionmov</a>=7), a
larger time step might
be taken, for example 200.
<br>
No meaning for RF calculations.
"""
},
'dynimage': {
'definition': "DYNamics of the IMAGE ",
'section': "varrlx",
'category': "",
'vartype': "integer dynimage(nimage)",
'default': """<b>dynimage(:)</b>=1; if <a href="varrlx.html#imgmov">imgmov</a>=2 or 5 (String Method, NEB), <b>dynimage(1)</b>=0 and <b>dynimage(<a href="varrlx.html#nimage">nimage</a>)</b>=0.""",
'text': """This input variable is relevant when sets of images are activated (see
<a href="varrlx.html#imgmov">imgmov</a>). Not all images might be required to evolve from one time step to the other.
Indeed, in the String Method or the Nudged Elastic Band, one might impose that the extremal configurations of the string are fixed.
In case the <b>dynimage</b>(iimage)=0, the image with index "iimage" will be consider as fixed.
Thus, there is no need to compute forces and stresses for this image at each time step. The purpose of
nevertheless defining extremal images is to make the input/output easier.<br><br>
In order to save CPU time, the computation of properties of static images (<b>dynimage</b>(iimage)=0)
can be avoided: see <a href="varrlx.html#istatimg">istatimg</a> keyword.
"""
},
'ecut': {
'definition': "Energy CUToff ",
'section': "varbas",
'category': """<a href="../users/abinit_help.html#energy">ENERGY</a> """,
'vartype': "real parameter ",
'default': "",
'text': """Used for kinetic energy cutoff
which controls number
of planewaves at given k point by:
<br>(1/2)[(2 Pi)*(k+Gmax)]<sup>2</sup>=<b>ecut</b> for Gmax.
<br>All planewaves inside this "basis sphere" centered
at k are included in the basis (except if <a href="varrlx.html#dilatmx">dilatmx</a>
is defined).
<br>Can be specified in Ha (the default), Ry, eV or Kelvin, since
<b>ecut</b> has the
'<a href="../users/abinit_help.html#dimensions">ENERGY</a>' characteristics.
(1 Ha=27.2113845 eV)
<br>This is the single parameter which can have an enormous
effect on the quality of a calculation; basically the larger
<b>ecut</b> is, the better converged the calculation is.  For fixed
geometry, the total energy MUST always decrease as <b>ecut</b> is
raised because of the variational nature of the problem.
<p>
<i>Usually one runs at least several calculations at various <b>ecut</b>
to investigate the convergence needed for reliable results.</i>
<p>
For k-points whose coordinates are build from 0 or 1/2,
the implementation of time-reversal symmetry that links
coefficients of the wavefunctions in reciprocal space
has been realized. See the input variable <a href="vardev.html#istwfk">istwfk</a>.
If activated (which corresponds to the Default mode),
this input variable <a href="vardev.html#istwfk">istwfk</a> will allow to
divide the number of plane wave (npw) treated explicitly
by a factor of two. Still, the final result should be identical with
the 'full' set of plane waves.
<p>See the input variable <a href="varrlx.html#ecutsm">ecutsm</a>, for the
smoothing of the kinetic energy, needed to optimize unit cell parameters.
"""
},
'ecuteps': {
'definition': "Energy CUT-off for EPSilon (the dielectric matrix) ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>, <a href="../users/abinit_help.html#energy">ENERGY</a> """,
'vartype': "real parameter",
'default': "0.0 ",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=3 or 4, that is, screening or sigma calculations.
<p>
<b>ecuteps</b> determines the cut-off energy of the planewave set used to represent the
independent-particle susceptibility $\chi^{(0)}_{KS}$, the dielectric matrix $\epsilon$, and its inverse.
<br>
It is not worth to take <b>ecuteps</b> bigger than four times <a href="vargw.html#ecutwfn">ecutwfn</a>,
this latter limit corresponding to the highest Fourier components of a wavefunction convoluted with itself.
Usually, even twice the value of <a href="vargw.html#ecutwfn">ecutwfn</a> might overkill. A value of <b>ecuteps</b>
between 5 and 10 Hartree often leads to converged results (at the level of 0.01 eV for the energy gap).
In any case, a convergence study is worth.
<p>
This set of planewaves can also be determined by the other input variables
<a href="vargw.html#npweps">npweps</a> and <a href="vargw.html#nsheps">nsheps</a>,
but these are much less convenient to use for general systems, than the selection
criterion based on a cut-off energy.
"""
},
'ecutsigx': {
'definition': "Energy CUT-off for SIGma eXchange ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>, <a href="../users/abinit_help.html#energy">ENERGY</a>""",
'vartype': "real parameter",
'default': "0.0 .",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=4, that is, sigma calculations.
<p>
<b>ecutsigx</b> determines the cut-off energy of the planewave set used to generate the
exchange part of the self-energy operator. For norm-conserving calculations, it is pointless
to have <b>ecutsigx</b> bigger than 4*<a href="varbas.html#ecut">ecut</a>,
while for PAW calculations, the maximal useful value is
<a href="varpaw.html#pawecutdg">pawecutdg</a>. Thus, if you do not care about CPU time, please use these
values. If you want to spare some CPU time, you might try to use a value between <a href="varbas.html#ecut">ecut</a> and these
upper limits.
<p>
This set of planewaves can also be determined by the other input variables
<a href="vargw.html#npwsigx">npwsigx</a> and <a href="vargw.html#nshsigx">nshsigx</a>,
but these are much less convenient to use for general systems, than the
selection criterion based on the cut-off energy (<b>ecutsigx</b> has to be 0.0 for using these).
"""
},
'ecutsm': {
'definition': "Energy CUToff SMearing ",
'section': "varrlx",
'category': """<a href="../users/abinit_help.html#energy">ENERGY</a> """,
'vartype': "real parameter (in Hartree) ",
'default': "0.d0",
'text': """This input variable is important when performing
relaxation of unit cell
size and shape (non-zero <a href="varrlx.html#optcell">optcell</a>).
Using a non-zero
<b>ecutsm</b>, the total energy curves as a function of
<a href="varbas.html#ecut">ecut</a>, or <a href="varbas.html#acell">acell</a>,
can be smoothed,
keeping consistency with
the stress (and automatically including the Pulay stress). The
recommended
value is 0.5 Ha. Actually, when <a href="varrlx.html#optcell">optcell</a>/=0,
ABINIT requires
<b>ecutsm</b> to be larger than zero. If you want to optimize cell
shape and size without
smoothing the total energy curve (a dangerous thing to do), use a very
small <b>ecutsm</b>,
on the order of one microHartree.

<p>
Technical information :
<br>
See Bernasconi et al, J. Phys. Chem. Solids 56, 501 (1995)
for a related method.
<br>
<b>ecutsm</b> allows to define an effective kinetic energy for plane
waves, close to, but
lower than the
maximal kinetic energy <a href="varbas.html#ecut">ecut</a>. For
kinetic
energies less than <a href="varbas.html#ecut">ecut</a>-<b>ecutsm</b>,
nothing is modified,
while between <a href="varbas.html#ecut">ecut</a>-<b>ecutsm</b> and <a
href="varbas.html#ecut">ecut</a>,
the kinetic energy is multiplied by:
<br>
1.0 / ( x<sup>2</sup> (3+x-6x<sup>2</sup>+3x<sup>3</sup>))
<br>
where x = (<a href="varbas.html#ecut">ecut</a> - kinetic_energy)/<b>ecutsm</b>
<br>
Note that x<sup>2</sup> ( 3+x-6x<sup>2</sup>+3x<sup>3</sup>) is 0 at
x=0, with vanishing derivative,
and that at x=1 , it is 1, with also vanishing derivative.
<br>
If <b>ecutsm</b> is zero, the unmodified kinetic energy is used.
<br>
<b>ecutsm</b> can be specified in Ha (the default), Ry, eV or Kelvin,
since
<b>ecutsm</b> has the
'<a href="../users/abinit_help.html#dimensions">ENERGY</a>'
characteristics.
(1 Ha=27.2113845 eV).
<br>
A few test for Silicon (diamond structure, 2 k-points) have
shown 0.5 Ha to be largely enough for <a href="varbas.html#ecut">ecut</a>
between 2Ha and 6Ha,
to get smooth curves. It is likely that this value is OK
as soon as <a href="varbas.html#ecut">ecut</a> is larger than 4Ha.
Note that prior to v4.4 (in v4.3 and older versions), the smoothing
function was
<br>
1.0 / ( x<sup>2</sup> ( 3-2*x) )
<br>
However, the second derivative was not continuous at the cut-off
energy. This
leads to problems for elastic constant calculations using response
functions.
"""
},
'ecutwfn': {
'definition': "Energy CUT-off for WaveFunctions",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>, <a href="../users/abinit_help.html#energy">ENERGY</a>""",
'vartype': "real parameter",
'default': """0.0, although, if a screening or self-energy calculation is done (<a href="vargs.html#optdriver">optdriver</a>=3 or 4), <b>ecutwfn</b> is automatically set to <a href="varbas.html#ecut">ecut</a> """,
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=3 or 4, that is, screening and sigma calculations.
<p>
<b>ecutwfn</b> determines the cut-off energy of the planewave set used to represent the wavefunctions
in the formula that generates the independent-particle susceptibility $\chi^{(0)}_{KS}$
(for <a href="vargs.html#optdriver">optdriver</a>=3), or the
self-energy (for <a href="vargs.html#optdriver">optdriver</a>=4).
<br>
Usually, <b>ecutwfn</b> is smaller than <a href="varbas.html#ecut">ecut</a>,
so that the wavefunctions are filtered, and some components are ignored.
As a side effect, the wavefunctions are no more normalized, and also, no more orthogonal.
Also, the set of plane waves can be much smaller for <a href="vargs.html#optdriver">optdriver</a>=3,
than for <a href="vargs.html#optdriver">optdriver</a>=4, although a convergence
study is needed to choose correctly both values.
<p>
This set of planewaves can also be determined by the other input variables
<a href="vargw.html#npwwfn">npwwfn</a> and <a href="vargw.html#nshwfn">nshwfn</a>,
but these are much less convenient to use for general systems, than the
selection criterion based on the cut-off energy (<b>ecutwfn</b> has to be set to 0.0 for using these).
"""
},
'effmass': {
'definition': "EFFective MASS ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>  """,
'vartype': "real number  ",
'default': "one.",
'text': """<br>This parameter allows to change the electron mass, with respect to its
experimental value.
"""
},
'efield': {
'definition': "Electric FIELD ",
'section': "varff",
'category': " ",
'vartype': "real array efield(3)  ",
'default': "3*0.0 .",
'text': """In case <a href="varff.html#berryopt">berryopt</a>=4,
a finite electric field calculation is performed. The value
of this electric field, and its direction is determined by <b>efield</b>.
It must be given in atomic units (1 a.u. of electric field= 514220624373.482 V/m, see note below),
in cartesian coordinates.
<p>
References for the calculation under electric field (based on multi k point Berry phase) :
<ul>
<li> Nunes and Vanderbilt, PRL 73, 712 (1994) : real-space version of the finite-field Hamiltonian </li>
<li> Nunes and Gonze, PRB 63, 155107 (2001) : reciprocal-space version of the finite-field Hamiltonian
(the one presently implemented), and extensive theoretical analysis </li>
<li> Souza, Iniguez and Vanderbilt, PRL 89, 117602 (2003) : implementation of the finite-field Hamiltonian
(reciprocal-space version)</li>
</ul>
See also Umari, Gonze, Pasquarello, PRL 90, 027401 (2003).
<p>
The atomic unit of electric field strength is :
e_Cb/(4 pi eps0 a0**2), where e_Cb is the electronic charge in Coulomb (1.60217653e-19),
eps0 is the electric constant (8.854187817d-12 F/m), and a0 is the Bohr radius
in meter (0.5291772108e-10).
"""
},
'elph2_imagden': {
'definition': "ELectron-PHonon interaction at 2nd order : IMAGina y shoft of the DENominator ",
'section': "varrf",
'category': """<a href="../users/abinit_help.html#respfn">RESPFN</a>,'<a href="../users/abinit_help.html#dimensions">ENERGY</a>'""",
'vartype': "real parameter ",
'default': "0.0 Ha.",
'text': """Only relevant if <a href="varrf.html#ieig2rf">ieig2rf</a> is non-zero, that is, if the user is performing performing second-order eigenvalue calculations using response-functions. <br><br>The variable <b>elph2_imagden</b> determines the imaginary shift of the denominator of the sum-over-states
in the perturbation denominator, (e_{nk}-e_{n'k'}+i <b>elph2_imagden</b>).
One should use a width comparable with the Debye frequency or the maximum phonon frequency.<br>Can be
specified in Ha (the default), Ry, eV or Kelvin, since
<b>ecut</b> has the
'<a href="../users/abinit_help.html#dimensions">ENERGY</a>' characteristics.
(1 Ha=27.2113845 eV)
"""
},
'enunit': {
'definition': "ENergy UNITs ",
'section': "vargs",
'category': " ",
'vartype': "integer parameter  ",
'default': "0 (eigenvalues in hartree and phonon frequencies in hartree and cm-1).",
'text': """Governs the units to be used for
output of eigenvalues (and eventual phonon frequencies)
<ul>
<li>0=>print eigenvalues in hartree;</li>
<li>1=>print eigenvalues in eV; </li>
<li>2=>print eigenvalues in both hartree and eV.  </li>
</ul>
If phonon frequencies are to be computed :
<ul>
<li>0=> phonon frequencies in Hartree and cm-1; </li>
<li>1=> phonon frequencies in eV and THz; </li>
<li>2=> phonon frequencies in hartree, eV, cm-1, Thz and Kelvin.</li>
</ul>"""
},
'eshift': {
'definition': "Energy SHIFT ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>, <a href="../users/abinit_help.html#energy">ENERGY</a>  """,
'vartype': "real number  ",
'default': "zero.",
'text': """Used only if <a href="vardev.html#wfoptalg">wfoptalg</a>=3 .
<b>eshift</b> gives the shift of the energy used in the
shifted Hamiltonian squared.
The algorithm will determine eigenvalues and eigenvectors centered
on <b>eshift</b>.
<br>Can be specified in Ha (the default), Ry, eV or Kelvin, since
<b>ecut</b> has the
'<a href="../users/abinit_help.html#dimensions">ENERGY</a>' characteristics.
(1 Ha=27.2113845 eV)
"""
},
'esmear': {
'definition': "Eigenvalue SMEARing ",
'section': "varrf",
'category': """<a href="../users/abinit_help.html#respfn">RESPFN</a>,'<a href="../users/abinit_help.html#dimensions">ENERGY</a>'   """,
'vartype': "real parameter ",
'default': "0.04 Ha.",
'text': """Only relevant if <a href="varrf.html#smdelta">smdelta</a> = 1-5, that is, if the user is performing simulations of the electronic lifetimes induced by the electron-phonon coupling.
<br><br>
The variable <b>esmear</b> determines the width of the functions approximating the delta function, \delta(e_{nk}-e_{n'k'}),
present in the expression of the lifetimes. One should use a width comparable with the Debye frequency or the maximum phonon frequency.<br>Can be specified in Ha (the default), Ry, eV or Kelvin, since
<b>ecut</b> has the
'<a href="../users/abinit_help.html#dimensions">ENERGY</a>' characteristics.
(1 Ha=27.2113845 eV)
"""
},
'etsfgroups': {
'definition': "ETSF I/O additional GROUPS of variables",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a> """,
'vartype': "integer parameter ",
'default': "0.",
'text': """<p>NOTE : NOT USED AT PRESENT (v5.3.0)
<p>This variable is a bit-wise combination of what will be written
into&nbsp;/&nbsp;read from a special WFK/DEN/POT file. The contents of the file
follow the <a href="http://www.etsf.eu/fileformats">Nanoquanta/ETSF file format specifications</a>.</p>
<p>Please check the "etsf_io" module of the ETSF I/O library for possible
values.</p>
"""
},
'etsfmain': {
'definition': "ETSF I/O MAIN variable",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a> """,
'vartype': "integer parameter ",
'default': "0.",
'text': """<p>NOTE : NOT USED AT PRESENT (v5.3.0)
<p>This variable tells what will be written into&nbsp;/&nbsp;read from a
special WFK/DEN/POT file. The contents of the file follow the
<a href="http://www.etsf.eu/fileformats">Nanoquanta/ETSF file format specifications</a>.</p>
<p>Please check the "etsf_io" module of the ETSF I/O library for possible
values.</p>
"""
},
'exchmix': {
'definition': "EXCHange MIXing",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a> """,
'vartype': "real number ",
'default': "0.25 ",
'text': """<b>exchmix</b> allows to tune the ratio of exact exchange when
<a href="varpaw.html#useexexch">useexexch</a> is used. The default value of 0.25 corresponds to PBE0.
"""
},
'exchn2n3d': {
'definition': "EXCHange N2 and N3 Dimensions",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a> """,
'vartype': "integer parameter ",
'default': "0. ",
'text': """If <b>exchn2n3d</b> is 1, the internal representation of the FFT arrays
in reciprocal space will be array(n1,n3,n2), where the second and
third dimensions have been switched. This is to allow to be coherent with the
<a href="vardev.html#exchn2n3d">exchn2n3d</a>=4xx FFT treatment.
"""
},
'f4of2_sla': {
'definition': "F4 Over F2 ratio of Slater integrals ",
'section': "varpaw",
'category': "",
'vartype': "real parameter",
'default': "0.4943",
'text': """Relevant only when
<a href="varint.html#usepaw">usepaw</a>=1 and
(<a href="varpaw.html#usepawu">usepawu</a>=1  or
<a href="vardev.html#usedmft">usedmft</a>=1) and
<a href="varpaw.html#lpawu">lpawu</a>=3) .<br>
This gives the ratio of Slater Integrals F6 and F2.
It is used with <a href="varpaw.html#f4of2_sla">f4of2_sla</a>=3) .<br>
in DFT+U or DFT+DMFT for the calculation of the orbital
dependent screened coulomb interaction.
"""
},
'f6of2_sla': {
'definition': "F6 Over F2 ratio of Slater integrals ",
'section': "varpaw",
'category': "",
'vartype': "real parameter",
'default': "0.4943",
'text': """Relevant only when
<a href="varint.html#usepaw">usepaw</a>=1 and
(<a href="varpaw.html#usepawu">usepawu</a>=1  or
<a href="vardev.html#usedmft">usedmft</a>=1) and
<a href="varpaw.html#lpawu">lpawu</a>=3) .<br>
This gives the ratio of Slater Integrals F6 and F2.
It is used with <a href="varpaw.html#f4of2_sla">f4of2_sla</a>=3) .<br>
in DFT+U or DFT+DMFT for the calculation of the orbital
dependent screened coulomb interaction.
"""
},
'fband': {
'definition': "Factor for the number of BANDs ",
'section': "vargs",
'category': """<a href="../users/abinit_help.html#not_internal">NOT INTERNAL</a> """,
'vartype': "real parameter, positive or zero  ",
'default': """0.125 in case <a href="varbas.html#occopt">occopt</a>==1 (insulating case), 0.500 for other values of <a href="varbas.html#occopt">occopt</a> (metallic case) and 0. in wavelet case. Not used in case <a href="varbas.html#occopt">occopt</a>==0 or 2.""",
'text': """Governs the number of bands to be used in the code in the case
the parameter <a href="varbas.html#nband">nband</a> is not defined in the input file
(which means that <a href="varbas.html#occopt">occopt</a> is not equal to 0 or 2).
<p>In case <b>fband</b> is 0.0d0, the code computes from
the pseudopotential files and the geometry data
contained in the input file, the number of electrons
present in the system. Then, it computes the minimum
number of bands that can accomodate them, and use
that value for <a href="varbas.html#nband">nband</a>.
<br>In case <b>fband</b> differs from
zero, other bands will be added, just
larger than <b>fband</b> times the number of atoms.
This parameter is not echoed in the top of the main
output file, but only the parameter <a href="varbas.html#nband">nband</a> that it allowed
to compute. It is also not present in the dtset array (no internal).
<br>The default values are chosen such as to give naturally some
conduction bands. This improves the robustness of the code,
since this allows to identify lack of convergence coming from
(near-)degeneracies at the Fermi level. In the metallic
case, the number of bands generated might be too small
if the smearing factor is large. The occupation numbers
of the higher bands should be small enough such as to
neglect higher bands. It is difficult to automate
this, so a fixed default value has been chosen."""
},
'fermie_nest': {
'definition': "FERMI Energy for printing the NESTing function ",
'section': "vardev",
'category': " ",
'vartype': "real parameter  ",
'default': "0",
'text': """This input variable is only effective when <a href="vardev.html#prtnest">prtnest</a>=1. The energy is relative to the calculated fermi energy.
"""
},
'fft_opt_lob': {
'definition': "Fast Fourier Transform parallelisation - OPTion for LOB algorithm ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a> """,
'vartype': "integer ",
'default': """1. If <a href="varpar.html#paral_kgb">paral_kgb</a>=1, default is 2.""",
'text': """Option for LOB algorithm, used in the band/FFT/k-point parallelisation,
see <a href="varpar.html#npband">npband</a>, <a href="varpar.html#npfft">npfft</a>,
<a href="varpar.html#npkpt">npkpt</a>, and <a href="varpar.html#paral_kgb">paral_kgb</a>.
<ul>
<li> =1 : old implementation </li>
<li> =2 : new implementation : the
calls to getghc are made in parallel on a set of bands
<a href="vardev.html#nbdblock">nbdblock</a> :
the aim is to reduce the number of collective communications. This is
not yet implemented in lobpcgwf. </li>
</ul>
"""
},
'fftalg': {
'definition': "Fast Fourier Transform ALGorithm ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a> """,
'vartype': "integer parameter ",
'default': """112, except for VPP Fujitsu, for which the Default is 111, and for NEC, for which the default is 200. Moreover, if the FFTW3 library has been enabled, the default becomes 312, EXCEPT if <a href="vardev.html#usedmft">usedmft</a> is non-zero for at least one dataset. <br> Finally, If <a href="varpar.html#paral_kgb">paral_kgb</a>=1, <b>fftalg</b> is automatically set to 401, with the hightest precedence.""",
'text': """This keyword is <b>irrelevant</b> when Fast Fourier Transforms are done using <b>Graphics Processing Units</b> (GPU),
i.e. when <a href="varpar.html#use_gpu_cuda">use_gpu_cuda</a>=1 (in that case, it is ignored).
<br><br>Allows to choose the algorithm
for Fast Fourier Transforms. These have to be used
when applied to wavefunctions (routine fourwf.f),
as well as when
applied to densities and potentials (routine fourdp.f).
Presently, it is the concatenation of three digits,
labelled (A), (B) and (C).
<br>
<br>The first digit (A) is to be chosen among 1, 2, 3 and 4 :
<ul>
<li> 1=&gt; use FFT routines written by S. Goedecker.</li>
<li> 2=&gt; use machine-dependent FFT algorithm, taken from the vendor
library, if it exists and if it has been implemented.
The bare <b>fftalg</b>=200 has little chance to be
faster than <b>fftalg</b>=112,
but it might be tried. Implementing library
subroutines with <b>fftalg</b>/=200 has not yet been done.
Currently implemented library subroutines (<b>fftalg</b>=200)
are:
<ul>
<li>      on HP, z3dfft from Veclib; </li>
<li>      on DEC Alpha, zfft_3d from DXML;</li>
<li>      on NEC, ZFC3FB from ASL lib;</li>
<li>      on SGI, zfft3d from complib.sgimath</li>
</ul>
</li>

<li>
3=&gt; use serial or multi-threaded FFTW fortran routines (<a href="http://www.fftw.org">http://www.fftw.org</a>).
Currently implemented with <b>fftalg</b>=300.</li>
<li> 4=&gt; use FFT routines written by S. Goedecker, 2002 version, that will
be suited for MPI and OpenMP parallelism.</li>
</ul>
The second digit (B) is related to fourdp.f :
<ul>
<li> 0=&gt; only use Complex-to-complex FFT</li>
<li> 1=&gt; real-to-complex is also allowed (only coded for A==1)</li>
</ul>
The third digit (C) is related to fourwf.f :
<ul>
<li> 0=&gt; no use of zero padding </li>
<li> 1=&gt; use of zero padding (only coded for A==1 and A==4)</li>
<li> 2=&gt; use of zero padding, and also combines actual
FFT operations (using 2 routines from S. Goedecker)
with important pre- and post-processing
operations, in order to maximize cache data reuse.
This is very efficient for cache architectures.
(coded for A==1 and A==4, but A==4 is not yet sufficiently tested)</li>
</ul>
Internal representation as <a href="vargs.html#ngfft">ngfft</a>(7)."""
},
'fftcache': {
'definition': "Fast Fourier Transform CACHE size ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a> """,
'vartype': "integer parameter ",
'default': "16. Not yet machine-dependent.",
'text': """Gives the cache size of the current
machine, in Kbytes.
<br>Internal representation as <a href="vargs.html#ngfft">ngfft</a>(8)."""
},
'fftgw': {
'definition': "FFT for GW calculation",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a> """,
'vartype': "integer parameter  ",
'default': "21.",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=3 or 4, that is, screening or sigma calculations.
<p>
The basic ingredients needed to perform both a screening and a sigma calculation are the so-called
oscillator matrix elements defined as

<br><br> $<
<B>k-q</B>, b1 | e^{-i (<B>q+G</B>).<B>r</B>} | <B>k</B> b2 >$ &nbsp;&nbsp;&nbsp;1)
<br><br>

In reciprocal space, Eq.1 is given by a convolution in which the number of reciprocal
lattice vectors employed to describe the wavefunctions is given
by <a href="vargw.html#ecutwfn">ecutwfn</a>
Tn the case of screening calculations, the number of <B>G</B> vectors in Eq.1 is defined
by <a href="vargw.html#ecuteps">ecuteps</a>,
while <a href="vargw.html#ecutsigx">ecutsigx</a> defined the number of <B>G</B> used in sigma calculations.
To improve the efficiency of the code, the oscillator matrix elements are evaluated
in real space through FFT techniques, and the <b>fftgw</b> input variable is used to select the FFT
mesh to be used.
<p>
<b>fftgw</b> is the concatenation of two digits, labelled (A) and (B) whose value is internally used
to define the value of  <a href="vargs.html#ngfft">ngfft</a>(1:3) (see the setmesh.F90 routine).
<p>
The first digit (A) defines the augmentation of the FFT grid. Possible values are 1, 2 and 3.

<ul>
<li> 0 =&gt; Use the FFT grid specified by the user through <a href="vargs.html#ngfft">ngfft</a>(1:3) </li>
<li> 1 =&gt; Use a coarse FFT grid which encloses a sphere in reciprocal space whose radius depends
on the largest value between
<a href="vargw.html#ecutwfn">ecutwfn</a> and <a href="vargw.html#ecuteps">ecuteps</a> </li>
<li> 2 =&gt;  Use a slightly augmented FFT which is sufficient for the correct treatment of the
convolution
<li> 3 =&gt; Doubled FFT grid (same mesh as that used for GS calculations).</li>
</ul>

<p>
The second digit (B) can be chosen between 0 and 1. It defines whether a FFT grid compatible with all
the symmetries of the space group must be enforced or not:

<ul>
<li> 0 =&gt; Use the smallest FFT mesh which is compatible with the FFT library (faster, save memory
but is less accurate)</li>
<li> 1 =&gt; Enforce a FFT grid which is compatible with all the symmetry operations of the space
group. This method leads to an increase both of CPU time and memory, but the matrix elements
are more accurate since the symmetry properties of the system are preserved.</li>
</ul>

The behaviour of ABINIT before v5.5 corresponds to the default value 11.
"""
},
'freqim_alpha': {
'definition': "FREQuencies along the IMaginary axis ALPHA parameter",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a> """,
'vartype': "real parameter  ",
'default': "5.0",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=4, that is, self-energy calculations.
<p>
<b>freqim_alpha</b> is used only for numerical integration of the GW self-energy
(<a href="vargw.html#gwcalctyp">gwcalctyp</a>= 2, 12, 22, 9, 19, 29).
<br>
<b>freqim_alpha</b> determines the location of the maximum frequency point along the imaginary axis
if the default grid is used in Contour Deformation (numerical integration) calculations. It is set
as <span style="font-family:Times,Serif"><i>&alpha;&middot;&omega;<sub>p<sub></i></span>, where
<span style="font-family:Times,Serif"><i>&omega;<sub>p<sub></i></span> is the plasma frequency
determined by the average density of the system (this can be set by hand by using the variable <a href="vargw.html#ppmfrq">ppmfrq</a>). See the section on grids in the desciptive text for <a href="vargw.html#cd_frqim_method">cd_frqim_method</a> for a detailed description of the formula.

"""
},
'freqremax': {
'definition': "FREQuencies along the Real axis MAXimum",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>, <a href="../users/abinit_help.html#energy">ENERGY</a>""",
'vartype': "real parameter  ",
'default': "0.0",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=3, that is, screening calculations.
<p>
<b>freqremax</b> is used only for numerical integration of the GW self-energy
(<a href="vargw.html#gwcalctyp">gwcalctyp</a>= 2, 12, 22, 9, 19, 29).
<br>
<b>freqremax</b> sets the maximum real frequency used to calculate the dielectric matrix in order
to perform the numerical integration of the GW self-energy.
<b>freqremax</b>, <a href="vargw.html#freqremin">freqremin</a> and <a href="vargw.html#nfreqre">nfreqre</a>
define the spacing of the frequency mesh along the real axis.
"""
},
'freqremin': {
'definition': "FREQuencies along the Real axis MINimum",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>, <a href="../users/abinit_help.html#energy">ENERGY</a>""",
'vartype': "real parameter  ",
'default': "0.0",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=3, that is, screening calculations.
<p>
<b>freqremin</b> is used only for numerical integration of the GW self-energy
(<a href="vargw.html#gwcalctyp">gwcalctyp</a>= 2, 12, 22, 9, 19, 29).
<br>
<b>freqremin</b> sets the minimum real frequency used to calculate the dielectric matrix in order
to perform the numerical integration of the GW self-energy.
<b>freqremin</b> can be used to split a wide frequency interval into smaller subintervals that
can be calculated independently.
The different subintervals can then be merged together with the <B>Mrgscr</B> utility thus obtaining
a single screening file that can used for self-energy calculations.

Note that <a href="vargw.html#freqremax">freqremin</a>, <a href="vargw.html#freqremin">freqremin</a>
and <a href="vargw.html#nfreqre">nfreqre</a> define the spacing of the frequency mesh along the real axis.

"""
},
'freqspmax': {
'definition': "FREQuencies for the SPectral function MAXimum",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>, <a href="../users/abinit_help.html#energy">ENERGY</a>""",
'vartype': "real parameter  ",
'default': "0.0",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=4, that is, sigma calculations.
<p>
<b>freqspmax</b> sets the maximum real frequency used to calculate the spectral function
from the GW Green's function. <a href="vargw.html#freqspmin">freqspmin</a>, <b>freqspmax</b> and
<a href="vargw.html#nfreqsp">nfreqsp</a> define the spacing of an equidistant frequency mesh along
the real axis. Alternatively, the variables <a href="vargw.html#gw_customnfreqsp">gw_customnfreqsp</a> and
<a href="vargw.html#gw_freqsp">gw_freqsp</a> can be used to make a user-defined grid.

"""
},
'freqspmin': {
'definition': "FREQuencies for the SPectral function MINimum",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>, <a href="../users/abinit_help.html#energy">ENERGY</a>""",
'vartype': "real parameter  ",
'default': """-<a href="vargw.html#freqspmax">freqspmax</a>""",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=4, that is, sigma calculations.
<p>
<b>freqspmin</b> sets the minimum real frequency used to calculate the spectral function
from the GW Green's function. <b>freqspmin</b> is set to -<a href="vargw.html#freqspmax">freqspmax</a>
if left undefined. <b>freqspmin</b>, <a href="vargw.html#freqspmax">freqspmax</a>, and
<a href="vargw.html#nfreqsp">nfreqsp</a> define the spacing of an equidistant frequency mesh along
the real axis. Alternatively, the variables <a href="vargw.html#gw_customnfreqsp">gw_customnfreqsp</a> and
<a href="vargw.html#gw_freqsp">gw_freqsp</a> can be used to make a user-defined grid.
"""
},
'freqsusin': {
'definition': "FREQuencies for the SUSceptibility matrix : the INcrement ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a> """,
'vartype': "real parameter, positive or zero ",
'default': "0.0",
'text': """Define, with
<a href="vardev.html#freqsuslo">freqsuslo</a>, the series of imaginary frequencies at which
the susceptibility matrix should be computed.
<br>This is still under development.
"""
},
'freqsuslo': {
'definition': "FREQuencies for the SUSceptibility matrix : the LOwest frequency ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a> """,
'vartype': "real parameter, positive or zero ",
'default': "0.0",
'text': """Define, with
<a href="vardev.html#freqsusin">freqsusin</a>,
the series of imaginary frequencies at which
the susceptibility matrix should be computed.
<br>This is still under development.
"""
},
'friction': {
'definition': "internal FRICTION coefficient ",
'section': "varrlx",
'category': "",
'vartype': "real parameter ",
'default': "0.001 ",
'text': """The equation of motion is :<br>
M<sub>I</sub> d<sup>2</sup>R<sub>I</sub>/dt<sup>2</sup>= F<sub>I</sub>
- <b>friction</b> M<sub>I</sub> dR<sub>I</sub>/dt - F_random<sub>I</sub>
<br>
where F_random<sub>I</sub> is a Gaussian random force with average
zero,
and variance 2 <b>friction</b> M<sub>I</sub> kT.
<br>
The atomic unit of friction is
hartrees*electronic mass*(atomic time units)/Bohr<sup>2</sup>. See J.
Chelikowsky, J. Phys. D : Appl Phys. 33(2000)R33.
"""
},
'frzfermi': {
'definition': "FReeZe FERMI energy",
'section': "varrf",
'category': "        ",
'vartype': "integer parameter",
'default': "0.",
'text': """Can be used to suppress artificially the first-order change of
Fermi energy, in case of Response Function calculation
for metals at Q=0. This change is needed, but was not computed prior to v4.4 .
Its calculation has been implemented by DHamann. The input variable <b>frzfermi</b>,
if set to 1, allows to recover the previous, incorrect behaviour.
"""
},
'fxcartfactor': {
'definition': "Forces to (X) CARTesian coordinates FACTOR ",
'section': "varrlx",
'category': "",
'vartype': "real parameter ",
'default': "(Bohr^2)/Hartree ",
'text': """The forces multiplied
by <b>fxcartfactor</b> will be treated like difference in cartesian coordinates in the
process of optimization. This is a simple preconditioner.
<br> TO BE UPDATED See (<a href="varrlx.html#ionmov">ionmov</a>=2,
non-zero
<a href="varrlx.html#optcell">optcell</a>).
For example, the stopping criterion defined by
<a href="varrlx.html#tolmxf">tolmxf</a> relates to these scaled
stresses.
"""
},
'ga_fitness': {
'definition': "Genetic Algorithm FITNESS function selection",
'section': "varrlx",
'category': "",
'vartype': "integer",
'default': "1 ",
'text': """Different methodologies to perform the roulette-wheel selection of parents. Even though, the
objective function is the crystalline enthalpy (H_i), the weight of the population elements to be chosen from
in a roulette-wheel selection can be given through different functions. We consider the following cases.
<br>
1. F = H_i / Sum H_i
<br>
2. F = exp(-(H_i-H_min)) / Sum exp(-(H_i-H_min))
<br>
3. F = (1/n_i) / Sum (1/n_i). Where n_i is the position in the ordered list of enthalpies
"""
},
'ga_n_rules': {
'definition': "Genetic Algorithm Number of RULES",
'section': "varrlx",
'category': "",
'vartype': "integer",
'default': "1 ",
'text': """Different genetic rules have been implemented and the user has the change to choose between any of them.
Right now we have 4 rules. See <a href="varbas.html#rprim">ga_rules</a>
"""
},
'ga_opt_percent': {
'definition': "Genetic Algorithm OPTIMAL PERCENT",
'section': "varrlx",
'category': "",
'vartype': "double",
'default': "0.2 ",
'text': """Percentage of the population that according to the fitness function passes
to the following iteration.
"""
},
'ga_rules': {
'definition': "Genetic Algorithm RULES",
'section': "varrlx",
'category': "",
'vartype': "integer vector",
'default': "1 ",
'text': """Different genetic rules have been implemented and the user has the change to choose between any of them.
The chosen number of rules have been defined in <a href="varbas.html#rprim">ga_n_rules</a>
<br>
<br>
Implemented rules are
<br>
1) crossover. Two parents are randomly chosen and two springs are mixed from the two by (a) chosing randomly (through
Fitness function) two parents and then randomly rotating and shifting the coordinates withing that particular cell.
(b) Slice every one of the unit cell of the parents along a random direction and creating the spring offs from the
pieces of the two parents.
<br>
2) Vector flip mutation. From the coordinates from a given parent, a piece of it is inverted.
<br>
3) random strain.  A random anisotropic deformation is given to the unit cell.
<br>
4) Coordinates mutation of 1/4 of the whole coordinates.
<br>
"""
},
'genafm': {
'definition': "GENerator of the translation for Anti-FerroMagnetic space group",
'section': "vargeo",
'category': """<a href="../users/abinit_help.html#symmetriser">SYMMETRISER</a>  """,
'vartype': "real genafm(3) ",
'default': "3*0.",
'text': """This input variable might be used to define a Shubnikov type IV magnetic space group (anti-ferromagnetic
space group). The user is advised to consult
"The mathematical theory of symmetry in solids,
Representation theory for point groups and space groups, 1972,
C.J. Bradley and A.P. Cracknell, Clarendon Press, Oxford."
<br>A Shubnikov type IV magnetic space group might be defined by its Fedorov space group
(set of spatial symmetries, that do not change the magnetization), and
one translation associated with a change of magnetization.
<b>genafm</b> is precisely this translation, in reduced coordinates (like <a href="varbas.html#xred">xred</a>)
<br>Thus, one way to specify a Shubnikov IV magnetic space group, is to define both
<a href="vargeo.html#spgroup">spgroup</a> and <b>genafm</b>. Alternatively, one might
define <a href="vargeo.html#spgroup">spgroup</a> and <a href="vargeo.html#spgroupma">spgroupma</a>,
or define by hand the set of symmetries, using <a href="varbas.html#symrel">symrel</a>,
<a href="varbas.html#tnons">tnons</a> and <a href="vargs.html#symafm">symafm</a>
"""
},
'get1den': {
'definition': "GET the first-order density from _DEN file ",
'section': "varfil",
'category': " ",
'vartype': "integer parameter ",
'default': "0.",
'text': """Relevant only when <a href="vargs.html#optdriver">optdriver</a>=5.
Indicate the files from which first-order densities must be obtained,
in multi-dataset mode (in single dataset mode, use
<a href="varfil.html#irdden">ird1den</a>).
<br>NOTE : a negative value of a "get" variable indicates the number of datasets to go backwards;
it is not the number to be subtracted from the current dataset to find the proper dataset.
As an example : <pre>ndtset 3   jdtset 1 2 4  getXXX -1</pre> refers to dataset 2 when dataset 4 is initialized.
"""
},
'get1wf': {
'definition': "GET the first-order wavefunctions from _1WF file ",
'section': "varfil",
'category': " ",
'vartype': "integer parameter ",
'default': "0.",
'text': """Eventually used when <a href="varbas.html#ndtset">ndtset</a>>0
(in the multi-dataset mode), to indicate
starting wavefunctions, as an alternative to <a href="varfil.html#irdwfk">irdwfk</a>,
<a href="varfil.html#irdwfk">irdwfq</a>,
<a href="varfil.html#irdwfk">ird1wf</a>,
or <a href="varfil.html#irdwfk">irdddk</a>. One should first read the
explanations given for these latter variables.
<br>
The <b>getwfk</b>, <b>getwfq</b>, <b>get1wf</b> and <b>getddk</b> variables are typically
used to chain the calculations in the multi-dataset mode,
since they describe from which dataset the OUTPUT
wavefunctions are to be taken, as INPUT wavefunctions
of the present dataset.
<br><br>We now focus on the <b>getwfk</b> input variable (the only
one used in ground-state calculations), but
the rules for <b>getwfq</b> and <b>get1wf</b> are similar, with _WFK
replaced by _WFQ or _1WF.
<br>If <b>getwfk</b>==0, no use of previously computed output
wavefunction file appended with _DSx_WFK is done.
<br>If <b>getwfk</b> is positive, its value gives the index of the dataset
for which the output wavefunction file appended with _WFK
must be used.
<br>If <b>getwfk</b> is -1, the output wf file with _WFK
of the previous dataset must be taken,
which is a frequently occuring case.
<br>If <b>getwfk</b> is a negative number, it indicates the number
of datasets to go backward to find the needed wavefunction file.
In this case, if one refers to a non existent data set (prior
to the first), the wavefunctions are not initialised from
a disk file, so that it is as if <b>getwfk</b>=0 for that
initialisation.
Thanks to this rule, the use of <b>getwfk</b> -1 is rather
straightforward : except for the first wavefunctions, that
are not initialized by reading a disk file, the output
wavefunction of one dataset is input of the next one.
<br>In the case of a ddk calculation in a multidataset
run, in order to compute
correctly the localisation tensor, it is mandatory to
declare give getddk the value of the current dataset
(i.e. getddk3 3 ) - this is a bit strange and
should be changed in the future.
<br>NOTE : a negative value of a "get" variable indicates the number of datasets to go backwards;
it is not the number to be subtracted from the current dataset to find the proper dataset.
As an example : <pre>ndtset 3   jdtset 1 2 4  getXXX -1</pre> refers to dataset 2 when dataset 4 is initialized.
"""
},
'getbscoup': {
'definition': "GET the Bethe-Salpeter COUPling block from ... ",
'section': "varfil",
'category': " ",
'vartype': "integer ",
'default': "0.",
'text': """Eventually used when <a href="varbas.html#ndtset">ndtset</a>>0
(multi-dataset mode) and, in the case of a Bethe-Salpeter
calculation to indicate that the starting coupling block of the excitonic Hamiltonian will be taken from the output of a previous dataset.
It is used to chain the calculations, since it describes from which dataset the OUTPUT coupling
block is to be taken, as INPUT of the present dataset.
<br>If <b>getbscoup</b>==0, no such use of previously computed coupling block file is done.
<br>If <b>getbscoup</b> is positive, its value gives the index of the dataset to be used as input.
<br>If <b>getbscoup</b> is -1, the output of the previous dataset must be taken, which is a frequently occuring case.
<br>If <b>getbscoup</b> is a negative number, it indicates the number
of datasets to go backward to find the needed file.
In this case, if one refers to a non existant data set (prior to the first), the coupling block is not initialised from
a disk file, so that it is as if <b>getbscoup</b>=0 for that initialisation.
"""
},
'getbseig': {
'definition': "GET the Bethe-Salpeter EIGenstates from ... ",
'section': "varfil",
'category': " ",
'vartype': "integer ",
'default': "0.",
'text': """Eventually used when <a href="varbas.html#ndtset">ndtset</a>>0
(multi-dataset mode) and, in the case of a Bethe-Salpeter
calculation to indicate that the starting excitonic eigenstates are to be taken from the output of a previous dataset.
It is used to chain the calculations, since it describes from which dataset the OUTPUT eigenstates
are to be taken, as INPUT eigenstates of the present dataset.
<br>If <b>getbseig</b>==0, no such use of previously computed output eigenstates file is done.
<br>If <b>getbseig</b> is positive, its value gives the index of the dataset
from which the output states is to be used as input.
<br>If <b>getbseig</b> is -1, the output eigenstates of the previous dataset
must be taken, which is a frequently occuring case.
<br>If <b>getbseig</b> is a negative number, it indicates the number
of datasets to go backward to find the needed file.
In this case, if one refers to a non existant data set (prior to the first), the eigenstates are not initialised from
a disk file, so that it is as if <b>getbseig</b>=0 for that initialisation.
"""
},
'getbsreso': {
'definition': "GET the Bethe-Salpeter RESOnant block from ... ",
'section': "varfil",
'category': " ",
'vartype': "integer ",
'default': "0.",
'text': """Eventually used when <a href="varbas.html#ndtset">ndtset</a>>0
(multi-dataset mode) and, in the case of a Bethe-Salpeter
calculation to indicate that the starting resonant block of the excitonic Hamiltonian will be taken from the output of a previous dataset.
It is used to chain the calculations, since it describes from which dataset the OUTPUT resonant
block is to be taken, as INPUT of the present dataset.
<br>If <b>getbsreso</b>==0, no such use of previously computed resonant block file is done.
<br>If <b>getbsreso</b> is positive, its value gives the index of the dataset to be used as input.
<br>If <b>getbsreso</b> is -1, the output of the previous dataset must be taken, which is a frequently occuring case.
<br>If <b>getbsreso</b> is a negative number, it indicates the number
of datasets to go backward to find the needed file.
In this case, if one refers to a non existant data set (prior to the first), the resonant block is not initialised from
a disk file, so that it is as if <b>getbsreso</b>=0 for that initialisation.
"""
},
'getcell': {
'definition': "GET CELL parameters from ... ",
'section': "varrlx",
'category': "",
'vartype': "integer parameter, an instance of a 'get' variable. ",
'default': "0.",
'text': """This variable is typically used to chain the
calculations,
in the multi-dataset mode (<a href="varbas.html#ndtset">ndtset</a>&gt;0),
since it describes from which dataset <a href="varbas.html#acell">acell</a>
and
<a href="varbas.html#rprim">rprim</a> are to be taken, as input of the
present
dataset. The cell parameters are EVOLVING variables,
for which such a chain of calculations is useful.
<br>
If ==0, no use of previously computed values must occur.
<br>
If it is positive, its value gives the index of the dataset
from which the data are to be used as input data.
It must be the index of a dataset already computed in the
SAME run.
<br>
If equal to -1, the output data of the previous dataset
must be taken, which is a frequently occuring case.
However, if the first dataset is treated, -1 is equivalent
to 0, since no dataset has yet been computed in the same run.
<br>
If another negative number, it indicates the number
of datasets to go backward to find the needed data
(once again, going back beyond the first dataset is equivalent
to using a null get variable).
"""
},
'getddk': {
'definition': "GET the ddk wavefunctions from _1WF file ",
'section': "varfil",
'category': " ",
'vartype': "integer parameter ",
'default': "0.",
'text': """Eventually used when <a href="varbas.html#ndtset">ndtset</a>>0
(in the multi-dataset mode), to indicate
starting wavefunctions, as an alternative to <a href="varfil.html#irdwfk">irdwfk</a>,
<a href="varfil.html#irdwfk">irdwfq</a>,
<a href="varfil.html#irdwfk">ird1wf</a>,
or <a href="varfil.html#irdwfk">irdddk</a>. One should first read the
explanations given for these latter variables.
<br>
The <b>getwfk</b>, <b>getwfq</b>, <b>get1wf</b> and <b>getddk</b> variables are typically
used to chain the calculations in the multi-dataset mode,
since they describe from which dataset the OUTPUT
wavefunctions are to be taken, as INPUT wavefunctions
of the present dataset.
<br><br>We now focus on the <b>getwfk</b> input variable (the only
one used in ground-state calculations), but
the rules for <b>getwfq</b> and <b>get1wf</b> are similar, with _WFK
replaced by _WFQ or _1WF.
<br>If <b>getwfk</b>==0, no use of previously computed output
wavefunction file appended with _DSx_WFK is done.
<br>If <b>getwfk</b> is positive, its value gives the index of the dataset
for which the output wavefunction file appended with _WFK
must be used.
<br>If <b>getwfk</b> is -1, the output wf file with _WFK
of the previous dataset must be taken,
which is a frequently occuring case.
<br>If <b>getwfk</b> is a negative number, it indicates the number
of datasets to go backward to find the needed wavefunction file.
In this case, if one refers to a non existent data set (prior
to the first), the wavefunctions are not initialised from
a disk file, so that it is as if <b>getwfk</b>=0 for that
initialisation.
Thanks to this rule, the use of <b>getwfk</b> -1 is rather
straightforward : except for the first wavefunctions, that
are not initialized by reading a disk file, the output
wavefunction of one dataset is input of the next one.
<br>In the case of a ddk calculation in a multidataset
run, in order to compute
correctly the localisation tensor, it is mandatory to
declare give getddk the value of the current dataset
(i.e. getddk3 3 ) - this is a bit strange and
should be changed in the future.
<br>NOTE : a negative value of a "get" variable indicates the number of datasets to go backwards;
it is not the number to be subtracted from the current dataset to find the proper dataset.
As an example : <pre>ndtset 3   jdtset 1 2 4  getXXX -1</pre> refers to dataset 2 when dataset 4 is initialized.
"""
},
'getden': {
'definition': "GET the DENsity from ... ",
'section': "varfil",
'category': " ",
'vartype': "integer parameter ",
'default': "0.",
'text': """Eventually used when <a href="varbas.html#ndtset">ndtset</a>>0
(multi-dataset mode) and, in the case of a ground-state
calculation, if <a href="varbas.html#iscf">iscf</a>&lt;0 (non-SCF calculation),
to indicate that the starting density is to be taken
from the output of a previous dataset.
It is used to chain the calculations,
since it describes from which dataset the OUTPUT density
are to be taken, as INPUT density of the present dataset.
<br>If <b>getden</b>==0, no such use of previously computed output
density file is done.
<br>If <b>getden</b> is positive, its value gives the index of the dataset
from which the output density is to be used as input.
<br>If <b>getden</b> is -1, the output density of the previous dataset
must be taken, which is a frequently occuring case.
<br>If <b>getden</b> is a negative number, it indicates the number
of datasets to go backward to find the needed file.
In this case, if one refers to a non existant data set (prior
to the first), the density is not initialised from
a disk file, so that it is as if <b>getden</b>=0 for that
initialisation.
Thanks to this rule, the use of <b>getden</b> -1 is rather
straightforward : except for the first density, that
is not initialized by reading a disk file, the output
density of one dataset is input of the next one.
<br>Be careful : the output density file of a run with
non-zero <a href="varrlx.html#ionmov">ionmov</a> does not have the proper name (it has a "TIM"
indication) for use as an input of an <a href="varbas.html#iscf">iscf</a>&lt;0 calculation.
<br>One should use the output density of a <a href="varrlx.html#ionmov">ionmov</a>==0 run.
<br>NOTE : a negative value of a "get" variable indicates the number of datasets to go backwards;
it is not the number to be subtracted from the current dataset to find the proper dataset.
As an example : <pre>ndtset 3   jdtset 1 2 4  getXXX -1</pre> refers to dataset 2 when dataset 4 is initialized.
"""
},
'getgam_eig2nkq': {
'definition': "GET the GAMma phonon data EIG2NKQ from dataset ",
'section': "vardev",
'category': "",
'vartype': "integer parameter ",
'default': "0.  ",
'text': """Only relevant if <a href="varrf.html#ieig2rf">ieig2rf</a> is non-zero, that is, if the user is performing performing second-order eigenvalue calcul
ations using response-functions. Also, relevant only for non-zero wavevectors <a href="vargs.html#qpt">qpt</a>
<br><br>From the electron-phonon matrix elements at some wavevector only, it is not possible to determine the Debye-Wallercontribution : one has to know also the q=Gamma electron-phonon matrix elements.
<br>
The variable <b>getgam_eig2nkq</b> allows to transmit the information about the second-order derivatives of the
eigenvalues for q=Gamma from the dataset where the calculation at Gamma was done, to the datasets
for other wavevectors.
"""
},
'gethaydock': {
'definition': "GET the Haydock restart file from ... ",
'section': "varfil",
'category': " ",
'vartype': "integer ",
'default': "0.",
'text': """Eventually used when <a href="varbas.html#ndtset">ndtset</a>>0
(multi-dataset mode) and, in the case of a Bethe-Salpeter
calculation to indicate that the Haydock iterative technique will be restarted from the output of a previous dataset.
<br>If <b>gethaydock</b>==0, no such use of previously computed coupling block file is done.
<br>If <b>gethaydock</b> is positive, its value gives the index of the dataset to be used as input.
<br>If <b>gethaydock</b> is -1, the output of the previous dataset must be taken, which is a frequently occuring case.
<br>If <b>gethaydock</b> is a negative number, it indicates the number
of datasets to go backward to find the needed file.
In this case, if one refers to a non existant data set (prior to the first), the coupling block is not initialised from
a disk file, so that it is as if <b>gethaydock</b>=0 for that initialisation.
"""
},
'getkss': {
'definition': "GET Kohn-Sham Structure from ... ",
'section': "varfil",
'category': "GW  ",
'vartype': "integer parameter ",
'default': "0.",
'text': """Used when <a href="varbas.html#ndtset">ndtset</a>>0
(multi-dataset mode) and <a href="vargs.html#optdriver">optdriver</a>=3 or 4
(a GW calculation),
to indicate that the KSS wavefunction file is to be taken
from the output of a previous dataset.
It is used to chain the calculations,
since it describes from which dataset the OUTPUT wavefunctions
are to be taken, as INPUT of the present dataset.
<br>If <b>getkss</b>==0, no such use of previously computed output
KSS file is done.
<br>If <b>getkss</b> is positive, its value gives the index of the dataset
from which the output KSS file is to be used as input.
<br>If <b>getkss</b> is -1, the output KSS file of the previous dataset
must be taken, which is a frequently occuring case.
<br>If <b>getkss</b> is a negative number, it indicates the number
of datasets to go backward to find the needed file.
In this case, if one refers to a non existent data set (prior
to the first), the KSS file is not initialised from
a disk file, so that it is as if <b>getkss</b>=0 for that
initialisation.
<br>NOTE : a negative value of a "get" variable indicates the number of datasets to go backwards;
it is not the number to be subtracted from the current dataset to find the proper dataset.
As an example : <pre>ndtset 3   jdtset 1 2 4  getXXX -1</pre> refers to dataset 2 when dataset 4 is initialized.
"""
},
'getocc': {
'definition': "GET OCC parameters from ... ",
'section': "varfil",
'category': " ",
'vartype': "integer parameter, an instance of a 'get' variable ",
'default': "0.",
'text': """This variable is typically used to chain the calculations,
in the multi-dataset mode (<a href="varbas.html#ndtset">ndtset</a>>0),
since it describes from which dataset the array <a href="vargs.html#occ">occ</a>
is to be taken, as input of the present
dataset. The occupation numbers are EVOLVING variables,
for which such a chain of calculations is useful.
<br>If ==0, no use of previously computed values must occur.
<br>If it is positive, its value gives the index of the dataset
from which the data are to be used as input data.
It must be the index of a dataset already computed in the
SAME run.
<br>If equal to -1, the output data of the previous dataset
must be taken, which is a frequently occuring case.
However, if the first dataset is treated, -1 is equivalent
to 0, since no dataset has yet been computed in the same run.
<br>If another negative number, it indicates the number
of datasets to go backward to find the needed data
(once again, going back beyond the first dataset is equivalent
to using a null get variable).
<br>NOTE that a non-zero <b>getocc</b> MUST be used with <a href="varbas.html#occopt">occopt</a>==2,
so that the number of bands has to be initialized for
each k point. Of course, these numbers of bands must be
identical with the numbers of bands of the dataset from which
<a href="vargs.html#occ">occ</a> will be copied. The same is true for the number of k points.
<br>NOTE : a negative value of a "get" variable indicates the number of datasets to go backwards;
it is not the number to be subtracted from the current dataset to find the proper dataset.
As an example : <pre>ndtset 3   jdtset 1 2 4  getXXX -1</pre> refers to dataset 2 when dataset 4 is initialized.
"""
},
'getqps': {
'definition': "GET QuasiParticle Structure ",
'section': "varfil",
'category': "GW ",
'vartype': "integer parameter ",
'default': "0.",
'text': """Used when <a href="varbas.html#ndtset">ndtset</a>>0
(multi-dataset mode) and <a href="vargs.html#optdriver">optdriver</a>=3, or 4
(screening or sigma step of a GW calculation),
to indicate that the eigenvalues and possibly the wavefunctions have to be
taken from a previous quasiparticle calculation (instead of the usual LDA starting
point). This is to achieve quasiparticle self-consistency.
See also <a href="varfil.html#irdqps">irdqps</a>
<br>NOTE : a negative value of a "get" variable indicates the number of datasets to go backwards;
it is not the number to be subtracted from the current dataset to find the proper dataset.
As an example : <pre>ndtset 3   jdtset 1 2 4  getXXX -1</pre> refers to dataset 2 when dataset 4 is initialized.
"""
},
'getscr': {
'definition': "GET SCReening (the inverse dielectric matrix) from ... ",
'section': "varfil",
'category': "GW ",
'vartype': "integer parameter ",
'default': "0.",
'text': """Used when <a href="varbas.html#ndtset">ndtset</a>>0
(multi-dataset mode) and <a href="vargs.html#optdriver">optdriver</a>=4
(sigma step of a GW calculation),
to indicate that the dielectric matrix (_SCR file) is to be taken
from the output of a previous dataset.
It is used to chain the calculations,
since it describes from which dataset the OUTPUT dielectric matrix
is to be taken, as INPUT of the present dataset.
<br>If <b>getscr</b>==0, no such use of previously computed output
_SCR file is done.
<br>If <b>getscr</b> is positive, its value gives the index of the dataset
from which the output _SCR file is to be used as input.
<br>If <b>getscr</b> is -1, the output _SCR file of the previous dataset
must be taken, which is a frequently occuring case.
<br>If <b>getscr</b> is a negative number, it indicates the number
of datasets to go backward to find the needed file.
In this case, if one refers to a non existent data set (prior
to the first), the _SCR file is not initialised from
a disk file, so that it is as if <b>getscr</b>=0 for that
initialisation.
<br>NOTE : a negative value of a "get" variable indicates the number of datasets to go backwards;
it is not the number to be subtracted from the current dataset to find the proper dataset.
As an example : <pre>ndtset 3   jdtset 1 2 4  getXXX -1</pre> refers to dataset 2 when dataset 4 is initialized.
"""
},
'getsuscep': {
'definition': "GET SUSCEPtibility (the irreducible polarizability) from ... ",
'section': "varfil",
'category': "GW ",
'vartype': "integer parameter ",
'default': "0.",
'text': """Used when <a href="varbas.html#ndtset">ndtset</a>>0
(multi-dataset mode) and <a href="vargs.html#optdriver">optdriver</a>=4
(sigma step of a GW calculation),
to indicate that the irreducible polarizability (_SUSC file) is to be taken
from the output of a previous dataset.
It is used to chain the calculations, since it describes from which dataset the OUTPUT susceptibility
is to be taken, as INPUT of the present dataset.
Performing a GW calculations starting from the _SUSC file instead of the _SCR file presents
the advantage that starting from the irreducible polarizability, one can calculate the screened interaction
using different expressions without having to perform a screening calculation from scratch.
For example, it is possible to apply a cutoff to the Coulomb interaction in order
to facilitate the convergence of the GW correction with respect to the size of the supercell
(see <a href="vargw.html#vcutgeo">vcutgeo</a> and <a href="vargw.html#icutcoul">icutcoul</a>)

<br>If <b>getsuscep</b>==0, no such use of previously computed output _SUSC file is done.
<br>If <b>getsuscep</b> is positive, its value gives the index of the dataset
from which the output _SUSC file is to be used as input.
<br>If <b>getsuscep</b> is -1, the output _SUSC file of the previous dataset
must be taken, which is a frequently occuring case.
<br>If <b>getsuscep</b> is a negative number, it indicates the number
of datasets to go backward to find the needed file.
In this case, if one refers to a non existent data set (prior
to the first), the _SUSC file is not initialised from
a disk file, so that it is as if <b>getsuscep</b>=0 for that
initialisation.
<br>NOTE : a negative value of a "get" variable indicates the number of datasets to go backwards;
it is not the number to be subtracted from the current dataset to find the proper dataset.
As an example : <pre>ndtset 3   jdtset 1 2 4  getXXX -1</pre> refers to dataset 2 when dataset 4 is initialized.
"""
},
'getvel': {
'definition': "GET VEL from ... ",
'section': "varrlx",
'category': "",
'vartype': "integer parameters, instances of 'get' variables ",
'default': "0.",
'text': """These variables are typically used to chain the
calculations,
in the multi-dataset mode (<a href="varbas.html#ndtset">ndtset</a>&gt;0)
since they describe from which dataset the corresponding
output variables are to be taken, as input of the present
dataset. The atomic positions and velocities are EVOLVING variables,
for which such a chain of calculation is useful.
<br>
Note that the use of <b>getxcart</b> and <b>getxred</b> differs when
<a href="varbas.html#acell">acell</a> and <a href="varbas.html#rprim">rprim</a>
are different from one dataset
to the other.
<br>
If ==0, no use of previously computed values must occur.
<br>
If it is positive, its value gives the index of the dataset
from which the data are to be used as input data.
It must be the index of a dataset already computed in the
SAME run.
<br>
If equal to -1, the output data of the previous dataset
must be taken, which is a frequently occuring case.
However, if the first dataset is treated, -1 is equivalent
to 0, since no dataset has yet been computed in the same run.
<br>
If another negative number, it indicates the number
of datasets to go backward to find the needed data
(once again, going back beyond the first dataset is equivalent
to using a null get variable).
<br>
Note : <b>getxred</b> and <b>getxcart</b> cannot be simultaneously
non-zero for the same dataset. On the other hand the use of
<b>getvel</b> with <b>getxred</b> is allowed, despite the different
coordinate system.
"""
},
'getwfk': {
'definition': "GET the wavefunctions from _WFK file ",
'section': "varfil",
'category': " ",
'vartype': "integer parameter ",
'default': "0.",
'text': """Eventually used when <a href="varbas.html#ndtset">ndtset</a>>0
(in the multi-dataset mode), to indicate
starting wavefunctions, as an alternative to <a href="varfil.html#irdwfk">irdwfk</a>,
<a href="varfil.html#irdwfk">irdwfq</a>,
<a href="varfil.html#irdwfk">ird1wf</a>,
or <a href="varfil.html#irdwfk">irdddk</a>. One should first read the
explanations given for these latter variables.
<br>
The <b>getwfk</b>, <b>getwfq</b>, <b>get1wf</b> and <b>getddk</b> variables are typically
used to chain the calculations in the multi-dataset mode,
since they describe from which dataset the OUTPUT
wavefunctions are to be taken, as INPUT wavefunctions
of the present dataset.
<br><br>We now focus on the <b>getwfk</b> input variable (the only
one used in ground-state calculations), but
the rules for <b>getwfq</b> and <b>get1wf</b> are similar, with _WFK
replaced by _WFQ or _1WF.
<br>If <b>getwfk</b>==0, no use of previously computed output
wavefunction file appended with _DSx_WFK is done.
<br>If <b>getwfk</b> is positive, its value gives the index of the dataset
for which the output wavefunction file appended with _WFK
must be used.
<br>If <b>getwfk</b> is -1, the output wf file with _WFK
of the previous dataset must be taken,
which is a frequently occuring case.
<br>If <b>getwfk</b> is a negative number, it indicates the number
of datasets to go backward to find the needed wavefunction file.
In this case, if one refers to a non existent data set (prior
to the first), the wavefunctions are not initialised from
a disk file, so that it is as if <b>getwfk</b>=0 for that
initialisation.
Thanks to this rule, the use of <b>getwfk</b> -1 is rather
straightforward : except for the first wavefunctions, that
are not initialized by reading a disk file, the output
wavefunction of one dataset is input of the next one.
<br>In the case of a ddk calculation in a multidataset
run, in order to compute
correctly the localisation tensor, it is mandatory to
declare give getddk the value of the current dataset
(i.e. getddk3 3 ) - this is a bit strange and
should be changed in the future.
<br>NOTE : a negative value of a "get" variable indicates the number of datasets to go backwards;
it is not the number to be subtracted from the current dataset to find the proper dataset.
As an example : <pre>ndtset 3   jdtset 1 2 4  getXXX -1</pre> refers to dataset 2 when dataset 4 is initialized.
"""
},
'getwfq': {
'definition': "GET the wavefunctions from _WFQ file ",
'section': "varfil",
'category': " ",
'vartype': "integer parameter ",
'default': "0.",
'text': """Eventually used when <a href="varbas.html#ndtset">ndtset</a>>0
(in the multi-dataset mode), to indicate
starting wavefunctions, as an alternative to <a href="varfil.html#irdwfk">irdwfk</a>,
<a href="varfil.html#irdwfk">irdwfq</a>,
<a href="varfil.html#irdwfk">ird1wf</a>,
or <a href="varfil.html#irdwfk">irdddk</a>. One should first read the
explanations given for these latter variables.
<br>
The <b>getwfk</b>, <b>getwfq</b>, <b>get1wf</b> and <b>getddk</b> variables are typically
used to chain the calculations in the multi-dataset mode,
since they describe from which dataset the OUTPUT
wavefunctions are to be taken, as INPUT wavefunctions
of the present dataset.
<br><br>We now focus on the <b>getwfk</b> input variable (the only
one used in ground-state calculations), but
the rules for <b>getwfq</b> and <b>get1wf</b> are similar, with _WFK
replaced by _WFQ or _1WF.
<br>If <b>getwfk</b>==0, no use of previously computed output
wavefunction file appended with _DSx_WFK is done.
<br>If <b>getwfk</b> is positive, its value gives the index of the dataset
for which the output wavefunction file appended with _WFK
must be used.
<br>If <b>getwfk</b> is -1, the output wf file with _WFK
of the previous dataset must be taken,
which is a frequently occuring case.
<br>If <b>getwfk</b> is a negative number, it indicates the number
of datasets to go backward to find the needed wavefunction file.
In this case, if one refers to a non existent data set (prior
to the first), the wavefunctions are not initialised from
a disk file, so that it is as if <b>getwfk</b>=0 for that
initialisation.
Thanks to this rule, the use of <b>getwfk</b> -1 is rather
straightforward : except for the first wavefunctions, that
are not initialized by reading a disk file, the output
wavefunction of one dataset is input of the next one.
<br>In the case of a ddk calculation in a multidataset
run, in order to compute
correctly the localisation tensor, it is mandatory to
declare give getddk the value of the current dataset
(i.e. getddk3 3 ) - this is a bit strange and
should be changed in the future.
<br>NOTE : a negative value of a "get" variable indicates the number of datasets to go backwards;
it is not the number to be subtracted from the current dataset to find the proper dataset.
As an example : <pre>ndtset 3   jdtset 1 2 4  getXXX -1</pre> refers to dataset 2 when dataset 4 is initialized.
"""
},
'getxcart': {
'definition': "GET XCART from ... ",
'section': "varrlx",
'category': "",
'vartype': "integer parameters, instances of 'get' variables ",
'default': "0.",
'text': """These variables are typically used to chain the
calculations,
in the multi-dataset mode (<a href="varbas.html#ndtset">ndtset</a>&gt;0)
since they describe from which dataset the corresponding
output variables are to be taken, as input of the present
dataset. The atomic positions and velocities are EVOLVING variables,
for which such a chain of calculation is useful.
<br>
Note that the use of <b>getxcart</b> and <b>getxred</b> differs when
<a href="varbas.html#acell">acell</a> and <a href="varbas.html#rprim">rprim</a>
are different from one dataset
to the other.
<br>
If ==0, no use of previously computed values must occur.
<br>
If it is positive, its value gives the index of the dataset
from which the data are to be used as input data.
It must be the index of a dataset already computed in the
SAME run.
<br>
If equal to -1, the output data of the previous dataset
must be taken, which is a frequently occuring case.
However, if the first dataset is treated, -1 is equivalent
to 0, since no dataset has yet been computed in the same run.
<br>
If another negative number, it indicates the number
of datasets to go backward to find the needed data
(once again, going back beyond the first dataset is equivalent
to using a null get variable).
<br>
Note : <b>getxred</b> and <b>getxcart</b> cannot be simultaneously
non-zero for the same dataset. On the other hand the use of
<b>getvel</b> with <b>getxred</b> is allowed, despite the different
coordinate system.
"""
},
'getxred': {
'definition': "GET XRED from ... ",
'section': "varrlx",
'category': "",
'vartype': "integer parameters, instances of 'get' variables ",
'default': "0.",
'text': """These variables are typically used to chain the
calculations,
in the multi-dataset mode (<a href="varbas.html#ndtset">ndtset</a>&gt;0)
since they describe from which dataset the corresponding
output variables are to be taken, as input of the present
dataset. The atomic positions and velocities are EVOLVING variables,
for which such a chain of calculation is useful.
<br>
Note that the use of <b>getxcart</b> and <b>getxred</b> differs when
<a href="varbas.html#acell">acell</a> and <a href="varbas.html#rprim">rprim</a>
are different from one dataset
to the other.
<br>
If ==0, no use of previously computed values must occur.
<br>
If it is positive, its value gives the index of the dataset
from which the data are to be used as input data.
It must be the index of a dataset already computed in the
SAME run.
<br>
If equal to -1, the output data of the previous dataset
must be taken, which is a frequently occuring case.
However, if the first dataset is treated, -1 is equivalent
to 0, since no dataset has yet been computed in the same run.
<br>
If another negative number, it indicates the number
of datasets to go backward to find the needed data
(once again, going back beyond the first dataset is equivalent
to using a null get variable).
<br>
Note : <b>getxred</b> and <b>getxcart</b> cannot be simultaneously
non-zero for the same dataset. On the other hand the use of
<b>getvel</b> with <b>getxred</b> is allowed, despite the different
coordinate system.
"""
},
'goprecon': {
'definition': "Geometry Optimization PREconditioner equations",
'section': "varrlx",
'category': "",
'vartype': "integer parameter",
'default': "0 ",
'text': """Set the kind of preconditioner to be used for Geometry Optimization
<br>
(Note : Under development now (2011.05.20))
<ul>
<li><b>goprecon</b>=0 : No preconditioner</li>
<li><b>goprecon</b>=[1-9] : Linear preconditioner</li>
<li><b>goprecon</b>=[11-19] : Non-linear preconditioner</li>
</ul>
<br>
"""
},
'goprecprm': {
'definition': "Geometry Optimization PREconditioner PaRaMeters equations",
'section': "varrlx",
'category': "",
'vartype': "real array of 3 elements",
'default': "0 ",
'text': """Set the paramenters use by the preconditioner to be
used for Geometry Optimization
<br>
(Note : Under development now (2011.06.06))
"""
},
'gw_customnfreqsp': {
'definition': "GW CUSTOM SPectral FREQuencies",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a> """,
'vartype': "integer ",
'default': "0 (off)",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=4, that is, sigma calculations,
and <a href="vargw.html#gwcalctyp">gwcalctyp</a>=2, 12 or 22.
<p>
<b>gw_customnfreqsp</b> lets the user define the gridpoints along the real frequency axis by hand
for the calculation of the self-energy along the real axis. Set this to the number of frequencies
you want. The frequencies are specified with <a href="vargw.html#gw_freqsp">gw_freqsp</a>.
"""
},
'gw_eet': {
'definition': "GW Effective Energy Technique ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>  """,
'vartype': "integer ",
'default': "-1 ",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=3 or 4, that is, screening or sigma calculations.
<p>
<b>gw_eet</b> governs the use of the effective energy technique (EET).
With the EET GW calculations can be performed without summing over the empty states, both in the calculation of the screening and the self-energy.
This is achieved by using a single effective energy that takes into account the contributions of all the empty states.
See J.A. Berger, L. Reining, and, F. Sottile, Phys. Rev. B (82), 041103(R) (2010) for a more detailed description of the method.
<ul>
<li> <b>gw_eet</b> = -1, no effective energy technique is used
<li> <b>gw_eet</b> = 0, 1, or 2 an effective energy is used. These numbers correspond to the order of the approxiation for the effective energy as given in Eqs. (10), (11), and (12), respectively, in J.A. Berger, L. Reining, and, F. Sottile, Phys. Rev. B (82), 041103(R) (2010). The higher the order, the more accurate (in principle) is the result. However, also the calculation time of the effective energy will increase.
</ul>
"""
},
'gw_eet_inclvkb': {
'definition': "GW Effective Energy Technique INCLude VKB ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>  """,
'vartype': "integer ",
'default': "0 ",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=3 or 4, that is, screening or sigma calculations
and if <a href="vargw.html#gw_eet">gw_eet</a>=1 or 2
<p>
<b>gw_eet_inclvkb</b> governs the use of the commutator of the non-local part of the pseudo-potential in determining the effective energy of the EET.
<ul>
<li> <b>gw_eet_inclvkb</b> = 0 the commutator is not included
<li> <b>gw_eet_inclvkb</b> = 1 the commutator is included
</ul>
"""
},
'gw_eet_nband': {
'definition': "GW Effective Energy Technique Number of BANDs ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>  """,
'vartype': "integer ",
'default': "-1 ",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=3 or 4, that is, screening or sigma calculations
and if <a href="vargw.html#gw_eet">gw_eet</a>=0, 1, or 2
<p>
<b>gw_eet_nband</b> governs the number of bands to be used in the effective energy technique (EET). It will only be taken into account if it is higher than the number of occupied bands. If this is not the case, summations will only include occupied states.
"""
},
'gw_eet_scale': {
'definition': "GW Effective Energy Technique SCALE ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>  """,
'vartype': "real parameter ",
'default': "1 ",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=3 or 4, that is, screening or sigma calculations
and if <a href="vargw.html#gw_eet">gw_eet</a>=0, 1, or 2
<p><b>gw_eet_scale</b> sets the prefactor with which the effective energy can be scaled to
obtain a faster convergence with respect to the number of empty states.
"""
},
'gw_freqsp': {
'definition': "GW SPectral FREQuencies ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a> """,
'vartype': "real array gw_freqsp(gw_customnfreqsp)",
'default': "none",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=4, that is sigma calculations.
<p>
<b>gw_freqsp</b> specifies the grid points for the real frequency axis when the real and imaginary
(spectral funtion) parts of sigma are calculated explicitly for post-processing or plotting.
Only activated if
<a href="vargw.html#gw_customnfreqsp">gw_customnfreqsp</a> is not equal to 0. The number of frequencies
is set by the value of <a href="vargw.html#gw_customnfreqsp">gw_customnfreqsp</a>. Example:
<pre>
gw_customnfreqsp   5
nfreqsp            5
gw_freqsp         -0.5  -0.1  0.0  1.0  10.0 eV
</pre>
If <a href="vargw.html#nfreqsp">nfreqsp</a> is not equal to
<a href="vargw.html#gw_customnfreqsp">gw_customnfreqsp</a> a warning will be issued.
</p>
"""
},
'gw_frqim_inzgrid': {
'definition': "Contour Deformation Imaginary Frequencies Inverse Z Grid",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a> """,
'vartype': "integer ",
'default': "0 (off)",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=3 or 4, that is, screening or sigma calculations,
and <a href="vargw.html#gwcalctyp">gwcalctyp</a>=2, 12 or 22.
<p>
<b>gw_frqim_inzgrid</b> creates gridpoints along the <b>imaginary</b> frequency axis by using an
equidistant grid in the variable
<span style="font-family:Times,Serif"><i>z &sub; [0,1]</i></span>
where the transform is:
<p align="center">
<img style="width: 122px; height: 38px;" src="./vargw_img/cd_inzgrid.png">
</p>
<p>Here <span style="font-family:Times,Serif"><i>&omega;<sub>p</sub></i></span> is the
plasma frequency (default can be overridden by setting <a href="vargw.html#ppmfrq">ppmfrq</a>).
The equidistant grid in z is determined uniquely by <a href="vargw.html#nfreqim">nfreqim</a>) and
the points are distributed so that half of them lie below the plasma frequency.
"""
},
'gw_frqre_inzgrid': {
'definition': "Contour Deformation Real Frequencies Inverse Z Grid",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a> """,
'vartype': "integer ",
'default': "0 (off)",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=3 or 4, that is, screening or sigma calculations,
and <a href="vargw.html#gwcalctyp">gwcalctyp</a>=2, 12 or 22.
<p>
<b>gw_frqre_inzgrid</b> creates gridpoints along the <b>real</b> frequency axis by using an
equidistant grid in the variable
<span style="font-family:Times,Serif"><i>z &sub; [0,1]</i></span>
where the transform is:
<p align="center">
<img style="width: 122px; height: 38px;" src="./vargw_img/cd_inzgrid_re.png">
</p>
<p>Here <span style="font-family:Times,Serif"><i>&omega;<sub>p</sub></i></span> is the
plasma frequency (default can be overridden by setting <a href="vargw.html#ppmfrq">ppmfrq</a>).
The equidistant grid in z is determined uniquely by <a href="vargw.html#nfreqim">nfreqre</a>) and
the points are distributed so that half of them lie below the plasma frequency. This is useful in
conjuction with <a href="vargw.html#gw_frqim_inzgrid">gw_frqim_inzgrid</a> if one needs to use a grid
which maps <span style="font-family:Times,Serif"><i>[0,&#8734;[ &rarr; [0,1]</i></span>. Note that
typically <i>many</i> more points are needed along the real axis in order to properly resolve
peak structures. In contrast, both the screening and self-energy are very smooth along the
imaginary axis. Also, please note that this is <b>not</b> an efficient grid for <b>standard</b>
Contour Deformation calculations, where typically only a smaller range of frequencies near the
origin is required. The maximum value needed along the real frequency axis is output in the
logfile during Contour Deformation sigma calculations.
"""
},
'gw_frqre_tangrid': {
'definition': "Contour Deformation Use Tangent Grid ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a> """,
'vartype': "integer ",
'default': "0 (off)",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=3 or 4, that is, screening or sigma calculations,
and <a href="vargw.html#gwcalctyp">gwcalctyp</a>=2, 12 or 22.
<p>
<b>gw_frqre_tangrid</b> defines a nonuniform grid to be used in frequency, with stepsize increasing
proportional to tan(x). This makes the grid approximately linear to start with, with a rapid increase
towards the end. Also, this is the grid which gives equal importance to each point used in the integration
of a function which decays as 1/x^2. To be used in conjuction with <a href="vargw.html#nfreqre">nfreqre</a>,
<a href="vargw.html#cd_max_freq">cd max_freq</a> and <a href="vargw.html#cd_halfway_freq">cd halfway_freq</a>
which determine the parameters of the transformed grid.
"""
},
'gw_nqlwl': {
'definition': "GW, Number of Q-points for the Long Wave-Length Limit ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>  """,
'vartype': "integer ",
'default': "1 ",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=3,4,99 that is, screening, sigma or BS
calculations, althought the actual meaning of the variable depends on the particular run-level (see
discussion below).
<p>
<b>gw_nqlwl</b> defines the number of directions in reciprocal space used to describe the non-analytical behaviour
of the heads (G = G'=0) and the wings (G=0 or G'=0) of the dielectric matrix in the optical limit
(i.e. for q tending to zero).
The number of directions is specified by the additional variable <a href="vargw.html#gw_qlwl">gw_qlwl</a>.
<p>
When <a href="vargs.html#optdriver">optdriver</a>=3, <b>gw_nqlwl</b> and <b>gw_qlwl</b> define
the set of "small" q that will be calculated and stored in the final SCR file. Therefore, the two
variables can be used to analyze how the optical spectra depend on the direction
of the incident phonon (useful especially in anisotropic systems).
<p>
When <a href="vargs.html#optdriver">optdriver</a>=4, <b>gw_nqlwl</b> and <b>gw_qlwl</b> can be used
to specify the heads and the wings to be used to perform the quadrature of the correlated
part of the self-energy in the small region around the origin.
(NB: not yet available, at present the quadrature is performed using a single direction in q-space)
<p>
When <a href="vargs.html#optdriver">optdriver</a>=99, <b>gw_nqlwl</b> and <b>gw_qlwl</b> define
the set of directions in q-space along which the macroscopic dielectric function is evaluated.
By default the BS code calculates the macroscopic dielectric function using six different
directions in q-space (the three basis vectors of the reciprocal lattice and the three Cartesian axis).
"""
},
'gw_nstep': {
'definition': "GW Number of self-consistent STEPS ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>  """,
'vartype': "integer ",
'default': "30 ",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=8, that is, self-consistent GW calculations.
<p>
Gives the maximum number of self-consistent GW cycles (or "iterations").
in which G and/or W will be updated until the quasi-particle energied are converged
within <a href="vargw.html#gw_toldfeig">gw_toldfeig</a>.
<a href="vargw.html#gwcalctyp">gwcalctyp</a> and <a href="vargw.html#gw_sctype">gw_sctype</a> are
used to define the type of self-consistency.
"""
},
'gw_qlwl': {
'definition': "GW, Q-points for the Long Wave-Length Limit ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>  """,
'vartype': """real array gw_qlwl(3,<a href="vargw.html#gw_nqlwl">gw_nqlwl</a>) """,
'default': "0.00001 0.00002 0.00003",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=3,4 that is, screening or sigma calculations.
<p>
When <a href="vargs.html#optdriver">optdriver</a>=3,
<b>gw_qlwl</b> defines the set of q-points around Gamma that are considered during the evaluation
of the non-analytical behaviour of the dielectrix matrix. Optical spectra (with and without
non-local field effects) are evaluated for each direction specified by <b>gw_qlwl</b>.
<p>
"""
},
'gw_sctype': {
'definition': "GW, Self-Consistency TYPE ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>  """,
'vartype': "integer ",
'default': "0 so that no choice is done",
'text': """<p>
This variable is used to partially define the kind of self-consistency for GW calculations.
The other piece of information is given by <a href="vargw.html#gwcalctyp">gwcalctyp</a>
that defines the particular approximation for the self-energy operator as well as whether the
wavefunctions have to replaced by quasi-particle amplitudes.
<p>
If <b>gw_sctype</b> is specified in the input file, the code will perform an iterative update
of the quantities entering the GW equations until the quasi-particle energies are converged
within <a href="vargw.html#gw_toldfeig">gw_toldfeig</a>. The maximum number of iterations is
specified by <a href="vargw.html#gw_nstep">gw_nstep</a>.
Possible values are:
<ul>
<li>  1 =&gt; standard one-shot method (one screening calculation followed by a single sigma run)
<li>  2 =&gt; self-consistency only on W (iterative update of W followed by a sigma run in which G
is approximated with the Kohn-Sham independent-particle Green's function G0)
<li>  3 =&gt; self-consistency only of G (a single screening calculation to obtain the Kohn-Sham
polarizability followed by an iterative update of the Green's functions in the self-energy)
<li>  4 =&gt; fully self-consistent algorithm (iterative update of both G and W)
</ul>

It is possible to initialize the self-consistent procedure by reading a previously calculated
SCR or SUSC file via the variables <a href="varfil.html#getscr">getscr</a> or
<a href="varfil.html#getsuscep">getsuscep</a>, respectively.
<a href="varfil.html#getqps">getqps</a> can be used to read a previous QPS file thus initializing
the Green's functions to be used in the first self-consistent iteration.

"""
},
'gw_sigxcore': {
'definition': "GW, treatment of the ... ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>  """,
'vartype': "integer ",
'default': "0 ",
'text': """<p>
Only available for PAW and relevant if <a href="vargs.html#optdriver">optdriver</a>=4 that is, sigma calculations.
<p>
Theoretical introduction: GW calculations perfomed on top of electronic calculations relying
when the frozen-core approximation is used to separate inner-core electrons from valence electrons,
only the contribution to the self-energy arising from valence electrons is explicitly accounted for.
In the standard approach based on pseudopotentials the contribution to the self-energy due to core electrons is approximated
by means of the KS exchange-correlation potential generated by the core density.
In the case of GW calculations employing the PAW method, the core contribution to the self-energy
can be more accurately estimated in terms of the Fock operator generated by the core wavefunctions.
In the simplest approach, the only ingredients required for this more refined treatment are the wave
functions of the core electrons in the reference atomic configuration that are calculated during the
generation of the PAW setup. This is a good approximation provided that the core wave functions
are strictly localized inside the PAW spheres.
<p>
<b>gw_sigxcore</b> defines the approximation used to evaluate the core contribution to sigma.
<ul>
<li> <b>gw_sigxcore</b> = 0, standard approach, the core contribution is approximated with vxc.
<li> <b>gw_sigxcore</b> = 1, the core term is approximated with the Fock operator inside the PAW spheres.
</ul>
<p>
"""
},
'gw_toldfeig': {
'definition': "GW TOLerance on the DiFference of the EIGenvalues ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>, <a href="../users/abinit_help.html#energy">ENERGY</a> """,
'vartype': "REAL ",
'default': "0.1 eV ",
'text': """<p>
Sets a tolerance for absolute differences of QP energies that will cause one self-consistent GW cycle to stop.
<br>Can be specified in Ha (the default), Ry, eV or Kelvin, since
<b>toldfe</b> has the '<a href="../users/abinit_help.html#dimensions">ENERGY</a>' characteristics (1 Ha=27.2113845 eV)
<br>Effective only when self-consistent GW calculations are done  (<a href="vargs.html#optdriver">optdriver</a>=8).
<p>
"""
},
'gwcalctyp': {
'definition': "GW CALCulation TYPe",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>  """,
'vartype': "integer parameter",
'default': "0 ",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=3 or 4, that is, screening or sigma calculations.
<p>
<b>gwcalctyp</b> governs the choice between the different capabilities of the GW code.

<ul>
<li>  0 &#60&#61 <b>gwcalctyp</b> &#60&#61  9 : standard "1 shot" quasiparticle method
<li> 10 &#60&#61 <b>gwcalctyp</b> &#60&#61 19 : self-consistent quasiparticle method on energies only
<li> 20 &#60&#61 <b>gwcalctyp</b> &#60&#61 29 : self-consistent quasiparticle method on energies and wavefunctions
</ul>

<p>
<ul>
<li> <b>gwcalctyp</b> &#61 0, 10, or 20 : standard Plasmon-Pole model GW calculation
<li> <b>gwcalctyp</b> &#61 2, 12, or 22 : GW calculation using numerical integration
(contour deformation method, see e.g. S. Lebegue <i> et al.</i> PRB <b>67</b>, 155208 (2003).)
<li> <b>gwcalctyp</b> &#61 5, 15, or 25 : Hartree-Fock calculation
<li> <b>gwcalctyp</b> &#61 6, 16, or 26 : Screened Exchange calculation
<li> <b>gwcalctyp</b> &#61 7, 17, or 27 : COHSEX calculation
<li> <b>gwcalctyp</b> &#61 8, 18, or 28 : model GW calculation
following S. Faleev <i>et al.</i> PRL <b>93</b>, 126406 (2004) using a Plasmon-Pole model
<li> <b>gwcalctyp</b> &#61 9, 19, or 29 : model GW calculation
following S. Faleev <i>et al.</i> PRL <b>93</b>, 126406 (2004) using numerical integration
(contour deformation method)
</ul>
"""
},
'gwcomp': {
'definition': "GW COMPletness ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>  """,
'vartype': "integer ",
'default': "0 ",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=3 or 4, that is, screening or sigma calculations.
<p>
<b>gwcomp</b> governs the use of an extrapolar approximation. If <b>gwcomp</b>==1, one improves the completeness
in a truncated sum over states. In practice, this permits one to reduce quite much the number of
bands required in the calculation of the screening or of the self-energy.
The energy parameter
needed in the extrapolar approximation is set by <a href="vargw.html#gwencomp">gwencomp</a>.
See F. Bruneval, X. Gonze, Phys. Rev. B 78, 085125 (2008) for a description of the methodology.
"""
},
'gwencomp': {
'definition': "GW Energy for COMPletness ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>  """,
'vartype': "energy ",
'default': "2. ",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=3 or 4, that is, screening or sigma calculations,
when <a href="vargw.html#gwcomp">gwcomp</a> is equal to 1.
<p>
<b>gwencomp</b> sets the energy parameter used in the extrapolar approximation used to improve
completeness and make the convergence against the number of bands much faster.
<p>
See F. Bruneval, X. Gonze, Phys. Rev. B 78, 085125 (2008) for a description of the methodology.
"""
},
'gwgamma': {
'definition': "GW Gamma ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>  """,
'vartype': "integer ",
'default': "0 ",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=4, that is, sigma calculations.
<p>
If <b>gwgamma</b> is 1, the vertex correction will be included leading to what is known as "GWGamma" approximation.
see R. Del Sole, L. Reining, and R. W. Godby, Phys. Rev. B  49, 8024 (1994).
Note that, in order to include the vertex correction in W, one has to start the sigma calculation
from the susceptibility file_SUSC instead of the _SCR file (see <a href="varfil.html#getsuscep">getsuscep</a>&nbsp;&nbsp;
and <a href="varfil.html#irdsuscep">irdsuscep</a>&nbsp;&nbsp;)
Not available for PAW calculations.
<p>
"""
},
'gwmem': {
'definition': "GW MEMory ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>  """,
'vartype': "integer ",
'default': "11 ",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=3,4, that is, sigma calculations.
<p>
<b>gwmem</b> governs the memory strategy during a screening and/or a sigma run.
<ul>
<li> <b>gwmem</b> = 1x , the screening matrix are read for all q-vectors and stored in the memory.
<li> <b>gwmem</b> = 0x , the screening matrix are read just a q-vector after another.
<br>
<li> <b>gwmem</b> = x1 , the real-space wavefunctions are stored in the memory.
<li> <b>gwmem</b> = x0 , the real-space wavefunctions are not stored, but rather recalculated on-fly each abinit needs them using FFTs.
</ul>
The default is <b>gwmem</b> = 11, which is the fastest, but also the most memory consuming.
When experiencing memory shortage, one should try <b>gwmem</b> = 0.
The first digit is only meaningful when performing sigma calculations.
<p>
"""
},
'gwpara': {
'definition': "GW PARAllelization level ",
'section': "varpar",
'category': "GW, PARALLEL",
'vartype': "integer ",
'default': "1 <B>TODO: default should be 2</b>.",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=3 or 4, that is, screening or sigma calculations.
<p>
<b>gwpara</b> is used to choose between the two different parallelization levels
available in the GW code. The available options are:

<ul>
<li> =1 =&gt; parallelisation on k points </li>
<li> =2 =&gt; parallelisation on bands    </li>
</ul>

<p>
Additional notes:
<br>
In the present status of the code, only the parallelization over bands (<b>gwpara</b>=2)
allows to reduce the memory allocated by each processor.
<br>
Using <b>gwpara</b>=1, indeed, requires the same amount of memory as a sequential run,
irrespectively of the number of CPU's used.
<p>
A reduction of the requireed memory can be achieved by opting for an out-of-core solution
(<a href="varfil.html#mkmem">mkmem</a>=0, only coded for <a href="vargs.html#optdriver">optdriver</a>=3)
at the price of a drastic worsening of the performance.
"""
},
'gwrpacorr': {
'definition': "GW RPA CORRelation energy ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>  """,
'vartype': "integer ",
'default': "0 ",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=3 and if <a href="vargw.html#gwcalctyp">gwcalctyp</a>=x1, i.e. a screening calculation for imaginary frequencies.
<p>
<b>gwrpacorr</b> governs the calculation of the RPA correlation energy.
<ul>
<li> <b>gwrpacorr</b> = 0, no RPA correlation energy is calculated
<li> <b>gwrpacorr</b> = 1, the RPA correlation energy is calculated using an exact integration over the coupling constant: it requires one diagonalization of the polarizability matrix
<li> <b>gwrpacorr</b> = <i>n</i> &gt 1, the RPA correlation energy is calculated using <i> n </i> values for the coupling constant: it requires <i>n</i> inversions of the polarizability matrix
</ul>
<p>
"""
},
'iatcon': {
'definition': "Indices of AToms in CONstraint equations",
'section': "varrlx",
'category': """<a href="../users/abinit_help.html#no_multi">NO MULTI</a>, <a href="../users/abinit_help.html#not_internal">NOT INTERNAL</a> """,
'vartype': """integer array iatcon(<a href="varrlx.html#natcon">natcon</a>,<a href="varrlx.html#nconeq">nconeq</a>)""",
'default': "0 ",
'text': """Gives the indices of the atoms appearing in each of the
<a href="varrlx.html#nconeq">nconeq</a>
independent equations constraining the motion of
atoms during structural optimization or molecular dynamics (see <a
href="varrlx.html#nconeq">nconeq</a>, <a href="varrlx.html#natcon">natcon</a>,
and <a href="varrlx.html#wtatcon">wtatcon</a>).
<br>
(Note : combined with wtatcon to give internal representation of the
latter - this should be described)
"""
},
'iatfix': {
'definition': "Indices of AToms that are FIXed ",
'section': "varrlx",
'category': """<b>iatfixx</b>,<b>iatfixy</b> and <b>iatfixz</b> are <a href="../users/abinit_help.html#not_internal">NOT INTERNAL</a> """,
'vartype': """integer arrays of length <a href="varrlx.html#natfix">natfix</a>, <a href="varrlx.html#natfixx">natfixx</a>, <a href="varrlx.html#natfixy">natfixy</a> or <a href="varrlx.html#natfixz">natfixz</a> """,
'default': """No Default (ignored unless <a href="varrlx.html#natfix">natfix</a>,<a href="varrlx.html#natfixx">natfixx</a>, <a href="varrlx.html#natfixy">natfixy</a> or <a href="varrlx.html#natfixz">natfixz</a> &gt; 0).""",
'text': """There is no harm in fixing one atom in the three
directions
using <b>iatfix</b>, then fixing it again in other directions
by mentioning it in <b>iatfixx</b>, <b>iatfixy</b> or <b>iatfixz</b>.
<br>
The internal representation of these input data is done
by the mean of one variable <b>iatfix</b>(3,<a href="varbas.html#natom">natom</a>),
defined
for each direction and each atom, being 0 if the atom is
not fixed along the direction, and 1 if the atom is fixed
along the direction.
When some atoms are fixed along 1 or 2 directions, the
use of symmetries is restricted to symmetry operations
whose (3x3) matrices <a href="varbas.html#symrel">symrel</a> are
diagonal.
<br>
If the geometry builder is used, <b>iatfix</b> will be related
to the preprocessed set of atoms, generated by the
geometry builder. The user must thus foresee the effect
of this geometry builder (see <a href="vargeo.html#objarf">objarf</a>).
"""
},
'iatfixx': {
'definition': "Indices of AToms that are FIXed along the X direction ",
'section': "varrlx",
'category': """<b>iatfixx</b>,<b>iatfixy</b> and <b>iatfixz</b> are <a href="../users/abinit_help.html#not_internal">NOT INTERNAL</a> """,
'vartype': """integer arrays of length <a href="varrlx.html#natfix">natfix</a>, <a href="varrlx.html#natfixx">natfixx</a>, <a href="varrlx.html#natfixy">natfixy</a> or <a href="varrlx.html#natfixz">natfixz</a> """,
'default': """No Default (ignored unless <a href="varrlx.html#natfix">natfix</a>,<a href="varrlx.html#natfixx">natfixx</a>, <a href="varrlx.html#natfixy">natfixy</a> or <a href="varrlx.html#natfixz">natfixz</a> &gt; 0).""",
'text': """There is no harm in fixing one atom in the three
directions
using <b>iatfix</b>, then fixing it again in other directions
by mentioning it in <b>iatfixx</b>, <b>iatfixy</b> or <b>iatfixz</b>.
<br>
The internal representation of these input data is done
by the mean of one variable <b>iatfix</b>(3,<a href="varbas.html#natom">natom</a>),
defined
for each direction and each atom, being 0 if the atom is
not fixed along the direction, and 1 if the atom is fixed
along the direction.
When some atoms are fixed along 1 or 2 directions, the
use of symmetries is restricted to symmetry operations
whose (3x3) matrices <a href="varbas.html#symrel">symrel</a> are
diagonal.
<br>
If the geometry builder is used, <b>iatfix</b> will be related
to the preprocessed set of atoms, generated by the
geometry builder. The user must thus foresee the effect
of this geometry builder (see <a href="vargeo.html#objarf">objarf</a>).
"""
},
'iatfixy': {
'definition': "Indices of AToms that are FIXed along the Y direction ",
'section': "varrlx",
'category': """<b>iatfixx</b>,<b>iatfixy</b> and <b>iatfixz</b> are <a href="../users/abinit_help.html#not_internal">NOT INTERNAL</a> """,
'vartype': """integer arrays of length <a href="varrlx.html#natfix">natfix</a>, <a href="varrlx.html#natfixx">natfixx</a>, <a href="varrlx.html#natfixy">natfixy</a> or <a href="varrlx.html#natfixz">natfixz</a> """,
'default': """No Default (ignored unless <a href="varrlx.html#natfix">natfix</a>,<a href="varrlx.html#natfixx">natfixx</a>, <a href="varrlx.html#natfixy">natfixy</a> or <a href="varrlx.html#natfixz">natfixz</a> &gt; 0).""",
'text': """There is no harm in fixing one atom in the three
directions
using <b>iatfix</b>, then fixing it again in other directions
by mentioning it in <b>iatfixx</b>, <b>iatfixy</b> or <b>iatfixz</b>.
<br>
The internal representation of these input data is done
by the mean of one variable <b>iatfix</b>(3,<a href="varbas.html#natom">natom</a>),
defined
for each direction and each atom, being 0 if the atom is
not fixed along the direction, and 1 if the atom is fixed
along the direction.
When some atoms are fixed along 1 or 2 directions, the
use of symmetries is restricted to symmetry operations
whose (3x3) matrices <a href="varbas.html#symrel">symrel</a> are
diagonal.
<br>
If the geometry builder is used, <b>iatfix</b> will be related
to the preprocessed set of atoms, generated by the
geometry builder. The user must thus foresee the effect
of this geometry builder (see <a href="vargeo.html#objarf">objarf</a>).
"""
},
'iatfixz': {
'definition': "Indices of AToms that are FIXed along the Z direction ",
'section': "varrlx",
'category': """<b>iatfixx</b>,<b>iatfixy</b> and <b>iatfixz</b> are <a href="../users/abinit_help.html#not_internal">NOT INTERNAL</a> """,
'vartype': """integer arrays of length <a href="varrlx.html#natfix">natfix</a>, <a href="varrlx.html#natfixx">natfixx</a>, <a href="varrlx.html#natfixy">natfixy</a> or <a href="varrlx.html#natfixz">natfixz</a> """,
'default': """No Default (ignored unless <a href="varrlx.html#natfix">natfix</a>,<a href="varrlx.html#natfixx">natfixx</a>, <a href="varrlx.html#natfixy">natfixy</a> or <a href="varrlx.html#natfixz">natfixz</a> &gt; 0).""",
'text': """There is no harm in fixing one atom in the three
directions
using <b>iatfix</b>, then fixing it again in other directions
by mentioning it in <b>iatfixx</b>, <b>iatfixy</b> or <b>iatfixz</b>.
<br>
The internal representation of these input data is done
by the mean of one variable <b>iatfix</b>(3,<a href="varbas.html#natom">natom</a>),
defined
for each direction and each atom, being 0 if the atom is
not fixed along the direction, and 1 if the atom is fixed
along the direction.
When some atoms are fixed along 1 or 2 directions, the
use of symmetries is restricted to symmetry operations
whose (3x3) matrices <a href="varbas.html#symrel">symrel</a> are
diagonal.
<br>
If the geometry builder is used, <b>iatfix</b> will be related
to the preprocessed set of atoms, generated by the
geometry builder. The user must thus foresee the effect
of this geometry builder (see <a href="vargeo.html#objarf">objarf</a>).
"""
},
'iatsph': {
'definition': "Index for the ATomic SPHeres of the atom-projected density-of-states",
'section': "vargs",
'category': " ",
'vartype': """integer array iatsph(1:<a href="vargs.html#natsph">natsph</a>)  """,
'default': """1, 2, ... <a href="vargs.html#natsph">natsph</a> """,
'text': """This input variable is active only in the
<a href="varfil.html#prtdos">prtdos</a>=3 case or if
<a href="varpaw.html#pawfatbnd">pawfatbnd</a>=1 or 2.
<br>It gives the number of the <a href="vargs.html#natsph">natsph</a> atoms around which the sphere
for atom-projected density-of-states will be build,
in the <a href="varfil.html#prtdos">prtdos</a>=3 case.
The radius of these spheres is given by <a href="vargs.html#ratsph">ratsph</a>.
<br> If <a href="varpaw.html#pawfatbnd">pawfatbnd</a>=1 or 2, it gives the number of the <a href="vargs.html#natsph">natsph</a> atoms around which atom-projected band structure will be built.
"""
},
'iboxcut': {
'definition': "",
'section': "varpaw",
'category': " ",
'vartype': "integer parameter ",
'default': "0 ",
'text': """Concern all summations in the reciprocal space and is
allowed in PAW and norm-conserving.
<ul>
<li>if set to 0 all reciprocal space summations are done in a sphere
contained in the FFT box.</li>
<li>if set to 1 all reciprocal space summations are done in the whole
FFT box (useful for tests).</li>
</ul>
"""
},
'icoulomb': {
'definition': "Coulomb TReaTMenT",
'section': "vargs",
'category': " ",
'vartype': "integer ",
'default': "0",
'text': """<p>Defines the type of computation used for Hartree potential, local part of pseudo-potential and ion-ion interaction:</p>
<ul>
<li><b>icoulomb</b>=0 : usual reciprocal space computation, using 1 / g^2 for the Hartree potential and using Ewald correction. </li>
<li><b>icoulomb</b>=1 : free boundary conditions are used when the Hartree potential is computed, real space expressions of pseudo-potentials are involved (restricted to GTH pseudo-potentials) and simple coulombian interaction gives the ion-ion energy. </li>
</ul>
"""
},
'icutcoul': {
'definition': "Integer that governs the CUT-off for COULomb interaction ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>  """,
'vartype': "integer parameter ",
'default': "3",
'text': """<p>
(This documentation is still very primitive, and should be checked)
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=3 or 4, that is, a screening or self-energy calculation.
<p>
Many-body calculations for isolated systems present a slow convergence with respect to the size of
the supercell due to the long ranged Coulomb interaction and the high degree of non-locality of the
operators involved. A similar issue occur also in fully periodic systems due to the presence
of the integrable Coulomb singularity at G=0 that hinders the convergence with respect to the number
of q-points used to sample the Brillouin zone.
The convergence can be accelerated by replacing the true bare Coulomb interaction with other
expressions
<b>icutcoul</b> defines the particular expression to be used for the Coulomb term in reciprocal
space. The choice of <b>icutcoul</b>  depends on the dimensionality of the system.
Possible values of <b>icutcoul</b> are 0 to 6. To be complemented by values of <a href="vargw.html#vcutgeo">vcutgeo</a>
and <a href="vargw.html#rcut">rcut</a>
<ul>
<li>0 => sphere (molecules but also 3D-crystals) </li>
<li>1 => cylinder (nanowires, nanotubes)</li>
<li>2 => surface </li>
<li>3 => 3D crystal (no cut-off, G=0 excluded, default value)</li>
<li>4 => ERF, long-range only Coulomb interaction </li>
<li>5 => ERFC, short-range only Coulomb interaction (e.g. as used in the HSE functional) </li>
<li>6 => auxiliary function integration for 3D systems from P. Carrier <i>et al.</i>, PRB <b>75</b>, 205126 (2007).</li>
</ul>
Note that Spencer and Alavi PRB <b>77</b>, 193110 (2008) showed that the spherical cutoff can efficiently be used also for 3D systems.
In the latter case, use a negative value for the cutoff radius of the sphere (<a href="vargw.html#rcut">rcut</a>&lt0),
which is automatically calculated so that the volume enclosed in the sphere is equal to the volume of the solid.
"""
},
'idyson': {
'definition': "Integer giving the choice of method for the DYSON equation",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a> """,
'vartype': "integer parameter ",
'default': "1.",
'text': """Choice for the method used to solve the Dyson equation in the calculation
of the interacting susceptibility matrix or/and in the calculation of the ACFD exchange-correlation energy:
<ul>
<li><b>idyson</b>=1 : Solve the Dyson equation by direct matrix inversion</li>
<li><b>idyson</b>=2 : Solve the Dyson equation as a first-order differential equation
with respect to the coupling constant lambda - only implemented for the RPA at the
present stage (see header of dyson_de.f for details)</li>
<li><b>idyson</b>=3 : Calculate only the diagonal of the interacting susceptibility matrix
by self-consistently computing the linear density change in response to a set
of perturbations. Only implemented for the RPA at the present stage,
and entirely experimental (see dyson_sc.f for details).</li>
</ul>
"""
},
'ieig2rf': {
'definition': "Integer for second-order EIGenvalues from Response-Function ",
'section': "varrf",
'category': """<a href="../users/abinit_help.html#respfn">RESPFN</a>  """,
'vartype': "integer parameters  ",
'default': "0.",
'text': """If <b>ieig2rf</b> is greater then 0, the code will produce a file, named with the trailing suffix _EIGR2D, containing the second-order electronic eigenvalues for the perturbation. These files are used in the calculation of the thermal correction to the electronic eigenvalues.
<br><br>
If <b>ieig2rf</b> is set to 1, the second-order electronic eigenvalues will be calculated from the DFPT method. <br>
If <b>ieig2rf</b> is set to 2, the second-order electronic eigenvalues will be calculated from the Allen-Cardona method.
<br>Related variables : <a href="varrf.html#bdeigrf">bdeigrf</a>,<a href="varrf.html#elph2_imagden">elph2_imagden</a>,<a href="vardev.html#getgam_eig2nkq">getgam_eig2nkq</a>,<a href="varrf.html#smdelta">smdelta</a>
"""
},
'ikhxc': {
'definition': "Integer option for KHXC = Hartree XC kernel ",
'section': "vardev",
'category': "",
'vartype': "integer parameter ",
'default': "1.",
'text': """Define the HXC kernel, in the cases for which it can be
dissociated with the choice of the HXC functional given by
<a href="varbas.html#ixc">ixc</a>, namely the TD-DFT computation of excited
states (<a href="varbas.html#iscf">iscf</a>=-1), and the computation of the
susceptibility matrix (for ACFD purposes). Options 2 to 6 are for the
ACFD only.
<ul>
<li>0 =&gt; RPA for the TDDFT but no kernel for the ACFD (testing purposes).</li>
<li>1 =&gt; RPA for the TDDFT and ACFD.</li>
<li>2 =&gt; ALDA (PW92) for the ACFD</li>
<li>3 =&gt; PGG for the ACFD [M. Petersilka, U.J. Gossmann and E.K.U. Gross, PRL 76,1212 (1996)]</li>
<li>4 =&gt; BPG for the ACFD. This amounts to half the PGG kernel plus half
the ALDA kernel for spin-compensated systems [K. Burke, M. Petersilka and E.K.U. Gross,
in "Recent Advances in Density Functional Methods", Vol. III, edited by P. Fantucci and A. Bencini
(World Scientific, Singapore, 2002)]</li>
<li>5 =&gt; Linear energy optimized kernel [J. Dobson and J. Wang, PRB 62, 10038 (2000)]</li>
<li>6 =&gt; Non-linear energy optimized kernel [J. Dobson and J. Wang, PRB 62, 10038 (2000)]</li>
</ul>
<br>
For ACFD-ALDA, BPG and energy optimized kernels are highly experimental and not tested yet !!!
For ACFD calculations, a cut-off density has been defined for the ALDA, BPG and
energy optimized kernels : let rhomin = userre*rhomax (where rhomax is the maximum density
in space) ; then the actual density used to calculate the local part of these kernels
at point r is max(rho(r),rhomin.
"""
},
'imgmov': {
'definition': "IMaGe MOVEs ",
'section': "varrlx",
'category': "",
'vartype': "integer parameter ",
'default': "0.",
'text': """<br>
No meaning for RF calculations.
"""
},
'inclvkb': {
'definition': "INCLude VKB ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>  """,
'vartype': "integer parameter ",
'default': "1 <b>TODO: should be changed to 0 or 2</b>",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=3,99 that is, a screening or BS calculations.
<p>
Possible values of <b>inclvkb</b> are 0,1,2.
If <b>inclvkb</b> is 1 or 2, the commutator of the non-local part of the pseudopotential
with the position operator is correctly included in the q =&gt; 0 contribution. This is
unfortunately time-consuming. When <b>inclvkb</b> is 0, this contribution is incorrectly
omitted, but the computation is much faster.
<p>
The importance of this contribution depends on the number of k points. Turning off
<b>inclvkb</b> should be made by experienced users only.
<p>
The use of <b>inclvkb</b>=2 is strongy recommended for the calculation of optical properties.
"""
},
'intexact': {
'definition': "INTegration using an EXACT scheme ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a> """,
'vartype': "integer parameter ",
'default': "0.",
'text': """Relates to the ACFD xc functionals only.
If <b>intexact</b> &gt; 0, the integration over the coupling constant
will be performed analytically in the RPA and in the two-electron PGG
approximation for the ACFD exchange-correlation energy.
Otherwise, the integration over the coupling constant will be performed
numerically (also see
<a href="vardev.html#ndyson">ndyson</a> and
<a href="vardev.html#idyson">idyson</a>.
Note that the program will stop in <b>intexact</b> &gt; 0 and
<a href="vardev.html#ikhxc">ikhxc</a>/=1 (RPA) or
<a href="vardev.html#ikhxc">ikhxc</a>/=3 (PGG, with two electrons)
"""
},
'intxc': {
'definition': "INTerpolation for eXchange-Correlation ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a> """,
'vartype': "integer parameter ",
'default': "0.",
'text': """<ul>
<li>0=&gt; do "usual" xc quadrature on fft grid</li>
<li>1=&gt; do higher accuracy xc quadrature using fft grid
and additional points at the centers of each cube
(doubles number of grid points)--the high accuracy version
is only valid for boxcut&gt;=2.  If boxcut &lt; 2, the code stops.</li>
</ul>
<br>For RF calculations only <b>intxc</b>=0 is allowed yet. Moreover,
the GS preparation runs (giving the density file and zero-order
wavefunctions) must be done with <b>intxc</b>=0
<p> Prior to ABINITv2.3, the choice <b>intxc</b>=1 was favoured (it was the default),
but the continuation of the development of the code lead to prefer
the default <b>intxc</b>=0 . Indeed, the benefit of <b>intxc</b>=1 is
rather small, while making it available for all cases is a
non-negligible development effort. Other targets are prioritary...
You will notice that many automatice tests use <b>intxc</b>=1. Please,
do not follow this historical choice for your production runs.
"""
},
'ionmov': {
'definition': "IONic MOVEs ",
'section': "varrlx",
'category': "",
'vartype': "integer parameter ",
'default': "0.",
'text': """<br>
No meaning for RF calculations.
"""
},
'iprcch': {
'definition': "Integer for PReConditioning of CHarge response ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a> """,
'vartype': "integer parameter ",
'default': """2.If <a href="varrlx.html#ionmov">ionmov</a>=4 and <a href="varbas.html#iscf">iscf</a>=5,<b>iprcch</b> is automatically put to 3. If <a href="varpar.html#paral_kgb">paral_kgb</a>=1, <b>iprcch</b> is automatically put to 6.""",
'text': """Used when <a href="varbas.html#iscf">iscf</a>>0, to define:<br>
- the way a change of density is derived from a change of atomic position,<br>
- the way forces are corrected when the SCF cycle is not converged.<br>
<br>Supported values :
<ul>
<li>0 =&gt; density not changed (fixed charge), forces not corrected </li>
<li>1 =&gt; density not changed, forces corrected with rigid ion hypothesis (atomic charge moved with atom)</li>
<li>2 =&gt; density changed and forces corrected with rigid ion hypothesis (atomic charge moves with atom)</li>
<li>3 =&gt; density changed and forces corrected with a different implementation of the rigid ion hypothesis</li>
<li>4 =&gt; density not changed, forces corrected with the use of Harris functional formula (*)</li>
<li>5 =&gt; density changed using D. Alfe 2nd-order algorithm (**), forces not corrected</li>
<li>6 =&gt; density changed using D. Alfe 2nd-order algorithm (**) and forces corrected with the use of Harris functional formula (*)</li>
</ul>
No meaning for RF calculations.<br>

<br>For the time being,<br>
- the choice 3 must be used with <a href="varrlx.html#ionmov">ionmov</a>=4
and <a href="varbas.html#iscf">iscf</a>=5.<br>
- the choices 4, 5 or 6 must be used when band-FFT parallelism is selected.<br>
Otherwise, use the choice 2.<br>

<br><b>(*)</b><U>Note concerning the use of <b>iprcch</b>=4 or 6 (correction of forces)</U>:<br>
The force on the atom located at R is corrected by the addition of the following term:<br>
<i>F_residual=Int[dr.V_residual.dRho_atomic/dR]</i>,  where Rho_atomic is an atomic (spherical) density.<br>
- When such an atomic density (Rho_atomic) is found in the pseudopotential or PAW file, it is used. If not, a gaussian density
(defined by <a href="vardev.html#densty">densty</a> parameter) is used.<br>
- When SCF mixing is done on the density (<a href="varbas.html#iscf">iscf</a>>=10), the potential residual (V_residual)
is obtained from the density residual with the first order formula <i>V_residual=dV/drho.Rho_residual</i>
and uses the exchange-correlation kernel <i>dVxc/drho=Kxc</i> which computation is time-consuming for GGA functionals.
By default the LDA exchange-correlation kernel is used (even for GGA, for which it seems to give a reasonable accuracy).
Using the exact GGA exchange correlation kernel is always possible by giving a negative value to <b>iprcch</b>.

<br><br><b>(**)</b><U>Note concerning the use of <b>iprcch</b>=5 or 6 (density prediction)</U>:<br>
The algorithm is described in <i>Computer Physics Communications <b>118</b> (1999) 31-33</i>.
It uses an atomic (spherical) density. When such an atomic density is found in the pseudopotential or PAW file, it is used. If not, a gaussian density
(defined by <a href="vardev.html#densty">densty</a> parameter) is used.<br>
Also note that, to be efficient, this algorithm requires a minimum convergency of the SCF cycle;
Typically, vres2 (or nres2) has to be small enough (10<sup>-4</sup>...10<sup>-5</sup>)."""
},
'iprcel': {
'definition': "Integer for PReConditioning of ELectron response ",
'section': "vargs",
'category': " ",
'vartype': "integer parameter ",
'default': "0.",
'text': """Used when <a href="varbas.html#iscf">iscf</a>>0, to define the SCF preconditioning scheme.
Potential-based preconditioning schemes for the SCF loop
(electronic part) are still a subject of active research.
The present parameter (electronic part) describes the way the
change of potential is derived from the residual.
<br>The possible values of <b>iprcel</b> correspond to :
<ul>
<li>0 => model dielectric function described by <a href="vargs.html#diemac">diemac</a>,
<a href="vargs.html#dielng">dielng</a>
and <a href="vargs.html#diemix">diemix</a>.</li>
<li>larger or equal to 21 => will compute the dielectric matrix
according to <a href="vargs.html#diecut">diecut</a>, <a href="vargs.html#dielam">dielam</a>,
<a href="vargs.html#diegap">diegap</a>. This methodology is described in
P.-M. Anglade, X. Gonze, Phys. Rev. B 78, 045126 (2008).</li>
<li>Between 21 and 29 => for the first few steps
uses the same as option 0 then compute RPA dielectric function,
and use it as such.</li>
<li>Between 31 and 39 => for the first few steps
uses the same as option 0 then compute RPA dielectric function,
and use it, with the mixing factor <a href="vargs.html#diemix">diemix</a>.</li>
<li>Between 41 and 49 => compute the RPA dielectric matrix
at the first step, and recompute it at a later step,
and take into account the mixing factor <a href="vargs.html#diemix">diemix</a>.</li>
<li>Between 51 and 59 => same as between 41 and 49, but compute
the RPA dielectric matrix by another mean</li>
<li>Between 61 and 69 => same as between 41 and 49, but compute
the electronic dielectric matrix instead of the RPA one.</li>
<li>Between 71 and 78 => STILL UNDER DEVELOPMENT -- NOT USABLE ; Use the modified Kerker preconditioner with a real-space formulation (basic formulation is shown at <a href="vargs.html#dielng">dielng</a>). The dielectric matrix is approximated thanks to  <a href="vargs.html#diemac">diemac</a> and <a href="vargs.html#dielng">dielng</a>.  Note that <a href="vargs.html#diemix">diemix</a> is also used.  </li>
<li> 79 => STILL UNDER DEVELOPMENT -- NOT USABLE ; same as previous but with an alternate algorithm. </li>
<li> 141 to 169 => same as Between 41 and 69 (but, the dielectric matrix is also recomputed every iprcel modulo 10 step). </li>
</ul>
<br>

The computation of the dielectric matrix (for 0 [100]< <b>iprcel</b> < 70 [100]) is based on the  <b>extrapolar</b> approximation. This approximation can be tuned with <a href="vargs.html#diecut">diecut</a>, <a href="vargs.html#dielam">dielam</a>,
and <a href="vargs.html#diegap">diegap</a>. Yet its accuracy mainly depends on the number of conduction bands included in the system. Having 2 to 10 empty bands in the calculation is
usually enough (use <a href="varbas.html#nband">nband</a>).
<br>
<br>
NOTES:
<ul>

<li>The step at which the dielectric matrix is computed or
recomputed is determined by modulo(<b>iprcel</b>,10). The recomputation happens
just once in the calculation for <b>iprcel</b> < 100.
<li>For non-homogeneous relatively large cells <b>iprcel</b>=45
will likely give a large improvement over <b>iprcel</b>=0.
<li>In case of PAW and <b>iprcel</b>>0, see <a href="varpaw.html#pawsushat">pawsushat</a> input variable. By default,
an approximation (which can be suppressed) is done for the computation of susceptibility matrix.
<li> For extremely large inhomogeneous cells where computation of the full dielectric matrix takes too many weeks, 70 < <b>iprcel</b> < 80 is advised.
<li>For <a href="varbas.html#nsppol">nsppol</a>=2 or
<a href="vargs.html#nspinor">nspinor</a>=2 with metallic <a href="varbas.html#occopt">occopt</a>,
only <b>mod(iprcel,100)</b><50 is allowed.
<li>No meaning for RF calculations yet.
<li> The exchange term in the full dielectric matrix diverges for vanishing densities.
Therefore the values of <b>iprcel</b> beyond 60 must not be used for cells containing vacuum,
unless ones computes this matrix for every step (<b>iprcel</b>=161).
</ul>
"""
},
'iprcfc': {
'definition': "Integer for PReConditioner of Force Constants ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a> """,
'vartype': "integer parameter ",
'default': "0.",
'text': """Used when <a href="varbas.html#iscf">iscf</a>>0, to define the SCF preconditioning scheme.
Potential-based preconditioning schemes for the SCF loop
are still under development.
<br>The present parameter (force constant part)
describes the way a change of force
is derived from a change of atomic position.
<br>Supported values :
<ul>
<li>0 =&gt; hessian is the identity matrix</li>
<li>1 =&gt; hessian is 0.5 times the identity matrix</li>
<li>2 =&gt; hessian is 0.25 times the identity matrix</li>
<li>-1=&gt; hessian is twice the identity matrix</li>
<li>... (simply corresponding power of 2 times the identity matrix)</li>
</ul>
No meaning for RF calculations."""
},
'iprctfvw': {
'definition': "Integer for PReConditioning (of electron response) based on Thomas - Fermi - von Weizs&auml;cker approximations of the kinetic energy ",
'section': "vargs",
'category': " ",
'vartype': "integer parameter ",
'default': "0.",
'text': """Used when <a href="varbas.html#iscf">iscf</a>>0, to use the TFvW preconditioner. This is still in an early DEVELOPMENT stage and is not usable as is.

<ul>
<li>0 => not TFvW preconditioning applied ; </li>
<li>1 => TFvW prc. is applied on the residuals left by other preconditioners ;</li>
<li>2 => Same as previously but with a linear response formulation of TFvW ; </li>
<li>3 => (starting at abinit 5.3.?) TFvW prc. is applied before any other preconditioner.</li>
</ul>

<br>For <a href="varbas.html#nsppol">nsppol</a> >=  2 only <b>iprcel</b>=0 is allowed.
<br>No meaning for RF calculations yet.
<br>Compatible only with potential mixing : <a href="varbas.html#iscf">iscf</a><10.

<p>
<br>References:
<ul>
<li> D. Raczkowski, A Canning and L.W. Wang PRB 64 121101 (2001)</li>
</ul>
</p>
"""
},
'iqpt': {
'definition': "Index for QPoinT generation ",
'section': "vargs",
'category': """<a href="../users/abinit_help.html#not_internal">NOT INTERNAL</a> """,
'vartype': "integer parameter ",
'default': "0 .",
'text': """Only used if <a href="vargs.html#nqpt">nqpt</a>=1, and
<a href="vargs.html#qptopt">qptopt</a>=1 to 4.
<p>Defines the index of the Q point to be selected in the list of q points generated by
<a href="vargs.html#ngqpt">ngqpt</a>,
<a href="vargs.html#qptrlatt">qptrlatt</a>,
<a href="vargs.html#nshiftq">nshiftq</a>,
and
<a href="vargs.html#shiftq">shiftq</a>.
<p>If <b>iqpt</b>=0, then the q point is Gamma (0 0 0).
<p>The usual working mode is to define a series of values for <b>iqpt</b>,
starting with <b>iqpt</b>=0 or 1 (so through the definition of <b>iqpt:</b>),
and increasing it by one for each dataset (thanks to <b>iqpt+</b>).
</p>
"""
},
'irandom': {
'definition': "Integer for the choice of the RANDOM number generator",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a> """,
'vartype': "integer parameter ",
'default': "3.",
'text': """For the time being, only used when <a href="varrlx.html#imgmov">imgmov</a>=9
(Langevin Path-Integral Molecular Dynamics).<br>
<b>irandom</b> defines the random number generator.<br>
<br>Supported values :
<ul>
<li>1 =&gt; "uniformrandom", delivered with ABINIT package (initially comes from numerical receipies).</li>
<li>2 =&gt; intrinsic Fortran 90 random number generator.</li>
<li>3 =&gt; "ZBQ" non-deterministic random number generator by R. Chandler and P. Northrop.
(This was delivered for free at http://www.ucl.ac.uk/~ucakarc/work/index.html#code", but
the latter address does not seem to work anymore. In any case, the initial copyright is not violated
by the files/documentation present in the ABINIT package).
</li>
</ul>
<b>irandom</b>=3 is strongly adviced when performing Molecular Dynamics restarts (avoids bias).
"""
},
'ird1den': {
'definition': "Integer that governs the ReaDing of 1st-order DEN file ",
'section': "varfil",
'category': " ",
'vartype': "integer parameter ",
'default': "0 when iscf &#62 0, 1 when iscf &#60 0.",
'text': """<br>
<br>Ground state calculation:
Start the ground-state calculation from the density file of a previous run.
When iscf &#60 0, the reading of a DEN file is always enforced.
<p>A non-zero value of <b>irdden</b> is treated in the same way
as other "ird" variables,
see the <a href="../users/abinit_help.html#4">section 4</a> of abinit_help.
<br>
<br>Response function calculation:
If first order density is needed in single dataset mode
(for example in nonlinear optical response),
use <b>ird1den</b> to read first-order densities from _DENx files produced in
other calculations. In multi-dataset mode use <a href="varfil.html#getden">get1den</a>.
"""
},
'ird1wf': {
'definition': "Integer that governs the ReaDing of _1WF files ",
'section': "varfil",
'category': " ",
'vartype': "integer parameter ",
'default': "0.",
'text': """Indicates eventual starting
wavefunctions. As alternative, one can use the
input variables <a href="varfil.html#getwfk">getwfk</a>,
<a href="varfil.html#getwfq">getwfq</a>, <a href="varfil.html#get1wf">get1wf</a>
or <a href="varfil.html#getddk">getddk</a>.
<br>
<br>Ground-state calculation :
<ul>
<li>only <b>irdwfk</b> and <a href="varfil.html#getwfk">getwfk</a> have a meaning</li>
<li>at most one of <b>irdwfk</b> or <a href="varfil.html#getwfk">getwfk</a> can be non-zero</li>
<li>if <b>irdwfk</b> and <a href="varfil.html#getwfk">getwfk</a> are both zero,
initialize wavefunctions with
random numbers for ground state calculation. </li>
<li>if <b>irdwfk</b> = 1 : read ground state wavefunctions
from a disk file appended with _WFK , produced in a
previous ground state calculation (see the
<a href="../users/abinit_help.html#4">section 4</a> of abinit_help).</li>
</ul>
Response-function calculation :
<ul>
<li>one and only one of <b>irdwfk</b> or <a href="varfil.html#getwfk">getwfk</a> MUST be non-zero</li>
<li>if <b>irdwfk</b> = 1 : read ground state k -wavefunctions
from a disk file appended with _WFK , produced in a
previous ground state calculation (see the
<a href="../users/abinit_help.html#4">section 4</a> of abinit_help).</li>
<li>only one of <b>irdwfq</b> or <a href="varfil.html#getwfq">getwfq</a> can be non-zero,
if both of them are non-zero,
use as k + q file the one defined by
<b>irdwfk</b> and/or <a href="varfil.html#getwfk">getwfk</a></li>
<li>if <b>irdwfq</b> = 1 : read ground state k+q -wavefunctions
from a disk file appended with _WFQ , produced in a
previous ground state calculation (see the
<a href="../users/abinit_help.html#4">section 4</a> of abinit_help).</li>
<li>at most one of <b>ird1wf</b> or <a href="varfil.html#get1wf">get1wf</a> can be non-zero</li>
<li>if both are zero, initialize first order wavefunctions
to 0's.</li>
<li>if <b>ird1wf</b> = 1 : read first-order wavefunctions
from a disk file appended with _1WFx , produced in a
previous response function calculation (see the
<a href="../users/abinit_help.html#4">section 4</a> of abinit_help).</li>
<li>at most one of <b>irdddk</b> or <a href="varfil.html#getddk">getddk</a> can be non-zero</li>
<li>one of them must be non-zero if an homogeneous
electric field calculation is done
(presently, a ddk calculation in the same dataset
is not allowed)</li>
<li>if <b>irdddk</b> = 1 : read first-order ddk wavefunctions
from a disk file appended with _1WFx , produced in a
previous response function calculation (see the
<a href="../users/abinit_help.html#4">section 4</a> of abinit_help).</li>
</ul>
"""
},
'irdbscoup': {
'definition': "Integer that governs the ReaDing of COUPling block ",
'section': "varfil",
'category': " ",
'vartype': "integer ",
'default': "0.",
'text': """Start the Bethe-Salpeter calculation from the BSC file containing the coupling block produced in a previous run.
"""
},
'irdbseig': {
'definition': "Integer that governs the ReaDing of BS_EIG file ",
'section': "varfil",
'category': " ",
'vartype': "integer ",
'default': "0.",
'text': """Start the Bethe-Salpeter calculation from the BS_EIG contining the exciton eigenvectors produced in a previous run.
"""
},
'irdbsreso': {
'definition': "Integer that governs the ReaDing of RESOnant block ",
'section': "varfil",
'category': " ",
'vartype': "integer ",
'default': "0.",
'text': """Start the Bethe-Salpeter calculation from the BSR file containing the resonat block produced in a previous run.
"""
},
'irdddk': {
'definition': "Integer that governs the ReaDing of DDK wavefunctions, in _1WF files ",
'section': "varfil",
'category': " ",
'vartype': "integer parameter ",
'default': "0.",
'text': """Indicates eventual starting
wavefunctions. As alternative, one can use the
input variables <a href="varfil.html#getwfk">getwfk</a>,
<a href="varfil.html#getwfq">getwfq</a>, <a href="varfil.html#get1wf">get1wf</a>
or <a href="varfil.html#getddk">getddk</a>.
<br>
<br>Ground-state calculation :
<ul>
<li>only <b>irdwfk</b> and <a href="varfil.html#getwfk">getwfk</a> have a meaning</li>
<li>at most one of <b>irdwfk</b> or <a href="varfil.html#getwfk">getwfk</a> can be non-zero</li>
<li>if <b>irdwfk</b> and <a href="varfil.html#getwfk">getwfk</a> are both zero,
initialize wavefunctions with
random numbers for ground state calculation. </li>
<li>if <b>irdwfk</b> = 1 : read ground state wavefunctions
from a disk file appended with _WFK , produced in a
previous ground state calculation (see the
<a href="../users/abinit_help.html#4">section 4</a> of abinit_help).</li>
</ul>
Response-function calculation :
<ul>
<li>one and only one of <b>irdwfk</b> or <a href="varfil.html#getwfk">getwfk</a> MUST be non-zero</li>
<li>if <b>irdwfk</b> = 1 : read ground state k -wavefunctions
from a disk file appended with _WFK , produced in a
previous ground state calculation (see the
<a href="../users/abinit_help.html#4">section 4</a> of abinit_help).</li>
<li>only one of <b>irdwfq</b> or <a href="varfil.html#getwfq">getwfq</a> can be non-zero,
if both of them are non-zero,
use as k + q file the one defined by
<b>irdwfk</b> and/or <a href="varfil.html#getwfk">getwfk</a></li>
<li>if <b>irdwfq</b> = 1 : read ground state k+q -wavefunctions
from a disk file appended with _WFQ , produced in a
previous ground state calculation (see the
<a href="../users/abinit_help.html#4">section 4</a> of abinit_help).</li>
<li>at most one of <b>ird1wf</b> or <a href="varfil.html#get1wf">get1wf</a> can be non-zero</li>
<li>if both are zero, initialize first order wavefunctions
to 0's.</li>
<li>if <b>ird1wf</b> = 1 : read first-order wavefunctions
from a disk file appended with _1WFx , produced in a
previous response function calculation (see the
<a href="../users/abinit_help.html#4">section 4</a> of abinit_help).</li>
<li>at most one of <b>irdddk</b> or <a href="varfil.html#getddk">getddk</a> can be non-zero</li>
<li>one of them must be non-zero if an homogeneous
electric field calculation is done
(presently, a ddk calculation in the same dataset
is not allowed)</li>
<li>if <b>irdddk</b> = 1 : read first-order ddk wavefunctions
from a disk file appended with _1WFx , produced in a
previous response function calculation (see the
<a href="../users/abinit_help.html#4">section 4</a> of abinit_help).</li>
</ul>
"""
},
'irdden': {
'definition': "Integer that governs the ReaDing of DEN file ",
'section': "varfil",
'category': " ",
'vartype': "integer parameter ",
'default': "0 when iscf &#62 0, 1 when iscf &#60 0.",
'text': """<br>
<br>Ground state calculation:
Start the ground-state calculation from the density file of a previous run.
When iscf &#60 0, the reading of a DEN file is always enforced.
<p>A non-zero value of <b>irdden</b> is treated in the same way
as other "ird" variables,
see the <a href="../users/abinit_help.html#4">section 4</a> of abinit_help.
<br>
<br>Response function calculation:
If first order density is needed in single dataset mode
(for example in nonlinear optical response),
use <b>ird1den</b> to read first-order densities from _DENx files produced in
other calculations. In multi-dataset mode use <a href="varfil.html#getden">get1den</a>.
"""
},
'irdhaydock': {
'definition': "Integer that governs the ReaDing of the HAYDOCK restart file ",
'section': "varfil",
'category': " ",
'vartype': "integer ",
'default': "0.",
'text': """Used to re-start the Haydock iterative technique from the HAYDR_SAVE file produced in a previous run.
"""
},
'irdkss': {
'definition': "Integer that governs the ReaDing of KSS file ",
'section': "varfil",
'category': "GW  ",
'vartype': "integer parameter ",
'default': "0.",
'text': """Relevant only when <a href="vargs.html#optdriver">optdriver</a>=3 or 4.
Indicate the file from which the dielectric matrix must be obtained.
As alternative, one can use the input variable
<a href="varfil.html#getkss">getkss</a>.
<br>When <a href="vargs.html#optdriver">optdriver</a>=3 or 4, at least one of
<b>irdkss</b> or <a href="varfil.html#getkss">getscr</a> must be non-zero.
<p>A non-zero value of <b>irdkss</b> is treated in the same way
as other "ird" variables,
see the <a href="../users/abinit_help.html#4">section 4</a> of abinit_help.
"""
},
'irdqps': {
'definition': "Integer that governs the ReaDing of QuasiParticle Structure ",
'section': "varfil",
'category': "GW  ",
'vartype': "integer parameter ",
'default': "0.",
'text': """Relevant only when <a href="vargs.html#optdriver">optdriver</a>=3 or 4.
Indicate the file from which the eigenvalues and possibly the wavefunctions must be obtained,
in order to achieve a self-consistent quasiparticle calculations.
See also <a href="varfil.html#getqps">getqps</a>
"""
},
'irdscr': {
'definition': "Integer that governs the ReaDing of the SCReening ",
'section': "varfil",
'category': "GW  ",
'vartype': "integer parameter ",
'default': "0.",
'text': """Relevant only when <a href="vargs.html#optdriver">optdriver</a>=4.
Indicate the file from which the dielectric matrix must be obtained.
As alternative, one can use the input variable
<a href="varfil.html#getscr">getscr</a>.
<br>When <a href="vargs.html#optdriver">optdriver</a>=4, at least one of
<b>irdscr</b> or <a href="varfil.html#getscr">getscr</a>
(alternatively, <a href="varfil.html#irdsuscep">irdsuscep</a> or <a href="varfil.html#getsuscep">getsuscep</a>)
must be non-zero.
<p>A non-zero value of <b>irdscr</b> is treated in the same way as other "ird" variables,
see the <a href="../users/abinit_help.html#4">section 4</a> of abinit_help.
"""
},
'irdsuscep': {
'definition': "Integer that governs the ReaDing of the SUSCEPtibility ",
'section': "varfil",
'category': "GW  ",
'vartype': "integer parameter ",
'default': "0.",
'text': """Relevant only when <a href="vargs.html#optdriver">optdriver</a>=4.
Indicate the file from which the irreducible polarizability must be obtained.
As alternative, one can use the input variable
<a href="varfil.html#getsuscep">getsuscep</a>.
<br>When <a href="vargs.html#optdriver">optdriver</a>=4, at least one of
<b>irdsuscep</b> or <a href="varfil.html#getsuscep">getsuscep</a>
(alternatively, <a href="varfil.html#irdscr">irdscr</a> or <a href="varfil.html#getscr">getscr</a>)
must be non-zero.
<p>A non-zero value of <b>irdsuscep</b> is treated in the same way as other "ird" variables,
see the <a href="../users/abinit_help.html#4">section 4</a> of abinit_help.
"""
},
'irdwfk': {
'definition': "Integer that governs the ReaDing of _WFK files ",
'section': "varfil",
'category': " ",
'vartype': "integer parameter ",
'default': "0.",
'text': """Indicates eventual starting
wavefunctions. As alternative, one can use the
input variables <a href="varfil.html#getwfk">getwfk</a>,
<a href="varfil.html#getwfq">getwfq</a>, <a href="varfil.html#get1wf">get1wf</a>
or <a href="varfil.html#getddk">getddk</a>.
<br>
<br>Ground-state calculation :
<ul>
<li>only <b>irdwfk</b> and <a href="varfil.html#getwfk">getwfk</a> have a meaning</li>
<li>at most one of <b>irdwfk</b> or <a href="varfil.html#getwfk">getwfk</a> can be non-zero</li>
<li>if <b>irdwfk</b> and <a href="varfil.html#getwfk">getwfk</a> are both zero,
initialize wavefunctions with
random numbers for ground state calculation. </li>
<li>if <b>irdwfk</b> = 1 : read ground state wavefunctions
from a disk file appended with _WFK , produced in a
previous ground state calculation (see the
<a href="../users/abinit_help.html#4">section 4</a> of abinit_help).</li>
</ul>
Response-function calculation :
<ul>
<li>one and only one of <b>irdwfk</b> or <a href="varfil.html#getwfk">getwfk</a> MUST be non-zero</li>
<li>if <b>irdwfk</b> = 1 : read ground state k -wavefunctions
from a disk file appended with _WFK , produced in a
previous ground state calculation (see the
<a href="../users/abinit_help.html#4">section 4</a> of abinit_help).</li>
<li>only one of <b>irdwfq</b> or <a href="varfil.html#getwfq">getwfq</a> can be non-zero,
if both of them are non-zero,
use as k + q file the one defined by
<b>irdwfk</b> and/or <a href="varfil.html#getwfk">getwfk</a></li>
<li>if <b>irdwfq</b> = 1 : read ground state k+q -wavefunctions
from a disk file appended with _WFQ , produced in a
previous ground state calculation (see the
<a href="../users/abinit_help.html#4">section 4</a> of abinit_help).</li>
<li>at most one of <b>ird1wf</b> or <a href="varfil.html#get1wf">get1wf</a> can be non-zero</li>
<li>if both are zero, initialize first order wavefunctions
to 0's.</li>
<li>if <b>ird1wf</b> = 1 : read first-order wavefunctions
from a disk file appended with _1WFx , produced in a
previous response function calculation (see the
<a href="../users/abinit_help.html#4">section 4</a> of abinit_help).</li>
<li>at most one of <b>irdddk</b> or <a href="varfil.html#getddk">getddk</a> can be non-zero</li>
<li>one of them must be non-zero if an homogeneous
electric field calculation is done
(presently, a ddk calculation in the same dataset
is not allowed)</li>
<li>if <b>irdddk</b> = 1 : read first-order ddk wavefunctions
from a disk file appended with _1WFx , produced in a
previous response function calculation (see the
<a href="../users/abinit_help.html#4">section 4</a> of abinit_help).</li>
</ul>
"""
},
'irdwfq': {
'definition': "Integer that governs the ReaDing of _WFQ files ",
'section': "varfil",
'category': " ",
'vartype': "integer parameter ",
'default': "0.",
'text': """Indicates eventual starting
wavefunctions. As alternative, one can use the
input variables <a href="varfil.html#getwfk">getwfk</a>,
<a href="varfil.html#getwfq">getwfq</a>, <a href="varfil.html#get1wf">get1wf</a>
or <a href="varfil.html#getddk">getddk</a>.
<br>
<br>Ground-state calculation :
<ul>
<li>only <b>irdwfk</b> and <a href="varfil.html#getwfk">getwfk</a> have a meaning</li>
<li>at most one of <b>irdwfk</b> or <a href="varfil.html#getwfk">getwfk</a> can be non-zero</li>
<li>if <b>irdwfk</b> and <a href="varfil.html#getwfk">getwfk</a> are both zero,
initialize wavefunctions with
random numbers for ground state calculation. </li>
<li>if <b>irdwfk</b> = 1 : read ground state wavefunctions
from a disk file appended with _WFK , produced in a
previous ground state calculation (see the
<a href="../users/abinit_help.html#4">section 4</a> of abinit_help).</li>
</ul>
Response-function calculation :
<ul>
<li>one and only one of <b>irdwfk</b> or <a href="varfil.html#getwfk">getwfk</a> MUST be non-zero</li>
<li>if <b>irdwfk</b> = 1 : read ground state k -wavefunctions
from a disk file appended with _WFK , produced in a
previous ground state calculation (see the
<a href="../users/abinit_help.html#4">section 4</a> of abinit_help).</li>
<li>only one of <b>irdwfq</b> or <a href="varfil.html#getwfq">getwfq</a> can be non-zero,
if both of them are non-zero,
use as k + q file the one defined by
<b>irdwfk</b> and/or <a href="varfil.html#getwfk">getwfk</a></li>
<li>if <b>irdwfq</b> = 1 : read ground state k+q -wavefunctions
from a disk file appended with _WFQ , produced in a
previous ground state calculation (see the
<a href="../users/abinit_help.html#4">section 4</a> of abinit_help).</li>
<li>at most one of <b>ird1wf</b> or <a href="varfil.html#get1wf">get1wf</a> can be non-zero</li>
<li>if both are zero, initialize first order wavefunctions
to 0's.</li>
<li>if <b>ird1wf</b> = 1 : read first-order wavefunctions
from a disk file appended with _1WFx , produced in a
previous response function calculation (see the
<a href="../users/abinit_help.html#4">section 4</a> of abinit_help).</li>
<li>at most one of <b>irdddk</b> or <a href="varfil.html#getddk">getddk</a> can be non-zero</li>
<li>one of them must be non-zero if an homogeneous
electric field calculation is done
(presently, a ddk calculation in the same dataset
is not allowed)</li>
<li>if <b>irdddk</b> = 1 : read first-order ddk wavefunctions
from a disk file appended with _1WFx , produced in a
previous response function calculation (see the
<a href="../users/abinit_help.html#4">section 4</a> of abinit_help).</li>
</ul>
"""
},
'iscf': {
'definition': "Integer for Self-Consistent-Field cycles ",
'section': "varbas",
'category': " ",
'vartype': "integer parameter ",
'default': "7 (norm-conserving), 17 (PAW) or 0 (WVL). (prior to v5.3 : 5 and 14) ",
'text': """Controls the self-consistency.
<br>Positive values =>
this is the usual choice for doing the usual ground state (GS)
calculations or for structural relaxation, where
the potential has to be determined self-consistently.
The choice between different algorithms for SCF is possible :
<ul>
<li>=0 => SCF cycle, direct minimization scheme on the gradient of the wavefunctions. This algorithm is faster than diagonalisation and mixing but is working only for systems with a gap. It is implemented only on the wavelet basis set, when <a href="varbas.html#usewvl">usewvl</a>=1.</li>
<li>=1 => get the largest eigenvalue of the SCF cycle
<br>(DEVELOP option, used with
<a href="varfil.html#irdwfk">irdwfk</a>=1 or
<a href="varfil.html#irdwfk">irdwfq</a>=1)</li>
<li>=2 => SCF cycle, simple mixing of the potential</li>
<li>=3 => SCF cycle, Anderson mixing of the potential</li>
<li>=4 => SCF cycle, Anderson mixing of the potential based on the two previous iterations</li>
<li>=5 => SCF cycle, CG based on the minim. of the energy with respect to the potential</li>
<li>=7 => SCF cycle, Pulay mixing of the potential based on the <a href="vardev.html#npulayit">npulayit</a> previous iterations</li>
<li>=12 => SCF cycle, simple mixing of the density</li>
<li>=13 => SCF cycle, Anderson mixing of the density</li>
<li>=14 => SCF cycle, Anderson mixing of the density based on the two previous iterations</li>
<li>=15 => SCF cycle, CG based on the minim. of the energy with respect to the density</li>
<li>=17 => SCF cycle, Pulay mixing of the density based on the <a href="vardev.html#npulayit">npulayit</a> previous iterations</li>
<li>Other positive values, including zero ones, are not allowed. </li>
</ul>
<p>
Such algorithms for treating the "SCF iteration history" should be coupled with accompanying algorithms
for the SCF "preconditioning". See the input variable <a href="vargs.html#iprcel">iprcel</a>.
The default value <a href="vargs.html#iprcel">iprcel</a>=0 is often a good choice, but
for inhomogeneous systems, you might gain a lot with <a href="vargs.html#iprcel">iprcel</a>=45.
<p>
(Warning : if <b>iscf</b>&gt;10, at present (v4.6), the energy printed at each SCF cycle is not variational -
this should not affect the other properties, and at convergence, all values are OK)
<p>
- In the norm-conserving case,
the default option is <b>iscf</b>=7, which is a compromise between speed and reliability.
The value <b>iscf</b>= 2 is safer but slower.
<br>
- In the PAW case, default option is <b>iscf</b>=17.
In PAW you have the possibility to mix density/potential on the fine or coarse FFT grid (see <a href="varpaw.html#pawmixdg">pawmixdg</a>).
<br>- Note that a Pulay mixing (<b>iscf</b>=7 or 17) with <a href="vardev.html#npulayit">npulayit</a>
=1 (resp. 2) is equivalent to an Anderson mixing with <b>iscf</b>=3 or 13 (resp. 4 or 14).
<br>- Also note that: <br>* when mixing is done on potential (iscf</b><10), total energy is computed by "direct" decomposition.
<br>* when mixing is done on density (iscf</b>>=10), total energy is computed by "double counting" decomposition.
<br>"Direct" and "double counting" decomposition of energy are equal when SCF cycle is converged. Note that,
when using GGA XC functionals, these decompositions of energy can be slighty different due
to imprecise computation of density gradients on FFT grid (difference decreases as size of FFT grid increases -
see <a href="varbas.html#ecut">ecut</a> for NC pseudopotentials, <a href="varpaw.html#pawecutdg">pawecutdg</a> for PAW).
<br>
<br>Other (negative) options:
<ul>
<li>= -2 =>
a non-self-consistent calculation is to be done;
in this case an electron density rho(r) on a real space grid
(produced in a previous calculation) will be read from a
disk file (automatically if <a href="varbas.html#ndtset">ndtset</a>=0, or
according to the value of <a href="varfil.html#getden">getden</a>
if <a href="varbas.html#ndtset">ndtset</a>/=0).
<br>The name of the density file must be given as indicated
in the <a href="../users/abinit_help.html#4">section 4</a> of abinit_help.
<b>iscf</b>=-2 would be used for
band structure calculations, to permit computation of
the eigenvalues of occupied and unoccupied states at
arbitrary k points in the fixed self consistent potential
produced by some integration grid of k points.
Due to this typical use, ABINIT insist that either
<a href="varfil.html#prtvol">prtvol</a>&gt;2 or
<a href="varfil.html#prteig">prteig</a> does not vanish
when there are more than 50 k points.
<br>To compute the eigenvalues
(and wavefunctions) of unoccupied states in a separate
(non-selfconsistent) run, the user should
save the self-consistent rho(r)
and then run <b>iscf</b>=-2 for the intended set of k-points and bands.
<br>To prepare a run with <b>iscf</b>=-2, a density file
can be produced using the
parameter <a href="varfil.html#prtden">prtden</a> (see its description).
When a self-consistent set of wavefunctions is already available,
abinit can be used with
<a href="varbas.html#nstep">nstep</a>=0 (see Test_v2/t47.in),
and the adequate value of <a href="varfil.html#prtden">prtden</a>.
</li>
<li>= -3 =>
like -2, but initialize <a href="vargs.html#occ">occ</a> and <a href="varbas.html#wtk">wtk</a>,
directly or indirectly (using <a href="varbas.html#ngkpt">ngkpt</a> or
<a href="vargs.html#kptrlatt">kptrlatt</a>)
depending on the value of <a href="varbas.html#occopt">occopt</a>.
<br>For GS, this option
might be used to generate Density-of-states
(thanks to <a href="varfil.html#prtdos">prtdos</a>),
or to produce STM charge density map (thanks to <a href="varfil.html#prtstm">prtstm</a>).
<br>For RF, this option is needed to compute the response to ddk perturbation.</li>

<li>= -1 => like -2, but the non-self-consistent calculation
is followed by the determination of excited states
within TDDFT. This is only possible for <a href="varbas.html#nkpt">nkpt</a>=1,
with <a href="varbas.html#kpt">kpt</a>=0 0 0.
Note that the oscillator strength needs to be defined with respect to
an origin of coordinate, thanks to the input variable
<a href="vargs.html#boxcenter">boxcenter</a>. The maximal
number of Kohn-Sham excitations to be used to build the
excited state TDDFT matrix can be defined by <a href="varrf.html#td_mexcit">td_mexcit</a>,
or indirectly by the maximum Kohn-Sham excitation energy
<a href="varrf.html#td_maxene">td_maxene</a>.
</li>
</ul>
"""
},
'isecur': {
'definition': "Integer for level of SECURity choice ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a> """,
'vartype': "integer ",
'default': "0. ",
'text': """In the presently used algorithms, there is a compromise
between speed and robustness, that can be tuned by
using <b>isecur</b>.
<br>If <b>isecur</b> =0, an extrapolation of out-of-line
data is allowed, and might save one non-SCF calculation every
two line minimisation when some stability conditions
are fulfilled (since there are 2 non-SCF calculations
per line minimisation, 1 out of 4 is saved)
<br>Using <b>isecur</b>=1 or higher integers will raise gradually
the threshold to make extrapolation.
<br>Using <b>isecur</b>=-2 will allow to save 2 non-SCF calculations
every three line minimisation, but this can make the
algorithm unstable. Lower values of <b>isecur</b> allows
for more (tentative) savings. In any case, there must
be one non-SCF computation per line minimisation.
<br>No meaning for RF calculations yet.  """
},
'istatimg': {
'definition': "Integer governing the computation of STATic IMaGes ",
'section': "varrlx",
'category': "",
'vartype': "integer parameter",
'default': "<b>istatimg</b>= 1",
'text': """This input variable is relevant when sets of images are activated (see
<a href="varrlx.html#imgmov">imgmov</a>).<br>
Not all images might be required to evolve from one time step to the other
(see<a href="varrlx.html#dynimage">dynimage</a>): these are static images.<br>
If <b>istatimg</b>=0, the total energy of static images is not computed (but static images are
used to make the dynamic images evolve). This can be useful to save CPU time.<br>
If <b>istatimg</b>=1, the total energy of static images is computed.
"""
},
'istatr': {
'definition': "Integer for STATus file repetition Rate ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>, <a href="../users/abinit_help.html#no_multi">NO_MULTI</a> """,
'vartype': "integer parameter ",
'default': "49, and 149 for Cray T3E (slow I/Os).Values lower than 10 may not work on some machines. Default <b>istatshft</b> is 1.",
'text': """Govern the rate of output of the status file.
This status file is written when the number of the
call to the status
subroutine is equal to '<b>istatshft</b>' modulo '<b>istatr</b>', so that
it is written once every '<b>istatr</b>' call.
There is also a writing for each of the 5 first calls,
and the 10th call."""
},
'istatshft': {
'definition': "Integer for STATus file SHiFT",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>, <a href="../users/abinit_help.html#no_multi">NO_MULTI</a> """,
'vartype': "integer parameter ",
'default': "49, and 149 for Cray T3E (slow I/Os).Values lower than 10 may not work on some machines. Default <b>istatshft</b> is 1.",
'text': """Govern the rate of output of the status file.
This status file is written when the number of the
call to the status
subroutine is equal to '<b>istatshft</b>' modulo '<b>istatr</b>', so that
it is written once every '<b>istatr</b>' call.
There is also a writing for each of the 5 first calls,
and the 10th call."""
},
'istwfk': {
'definition': "Integer for choice of STorage of WaveFunction at each k point ",
'section': "vardev",
'category': " ",
'vartype': """integer array istwfk(<a href="varbas.html#nkpt">nkpt</a>)  """,
'default': """0 for all k points for GS calculations. For RF calculations, the Default is not used : <b>istwfk</b> is forced to be 1 deep inside the code, for all k points. For spin-orbit calculations (<a href="vargs.html#nspinor">nspinor</a>=2), <b>istwfk</b> is also forced to be 1, for all k points.""",
'text': """Control the way the
wavefunction for each k-point is stored inside ABINIT,
in reciprocal space.
<br>For the GS calculations, in the "cg" array containing the
wavefunction coefficients, there is for each k-point
and each band, a segment cg(1:2,1:npw). The 'full' number
of plane wave is determined by <a href="varbas.html#ecut">ecut</a>.
However, if the k-point coordinates are build
only from zeroes and halves (see list below),
the use of time-reversal symmetry (that connects coefficients)
has been implemented, in order to use real-to-complex
FFTs (see <a href="vardev.html#fftalg">fftalg</a>), and to treat explicitly only half
of the number of plane waves (this being used as 'npw').
<br>For the RF calculations, there is not only the "cg"
array, but also the "cgq" and "cg1" arrays. For the
time-reversal symmetry to decrease the number of
plane waves of these arrays, the q vector MUST be (0 0 0).
Then, for each k point, the same rule as for the
RF can be applied.
<br>WARNING (991018) : for the time being, the time-reversal
symmetry cannot be used in the RF calculations.
<ul>
<li>1=&gt; do NOT take advantage of the time-reversal symmetry</li>
<li>2=&gt; use time-reversal symmetry for k=( 0   0   0 )</li>
<li>3=&gt; use time-reversal symmetry for k=(1/2  0   0 )</li>
<li>4=&gt; use time-reversal symmetry for k=( 0   0  1/2)</li>
<li>5=&gt; use time-reversal symmetry for k=(1/2  0  1/2)</li>
<li>6=&gt; use time-reversal symmetry for k=( 0  1/2  0 )</li>
<li>7=&gt; use time-reversal symmetry for k=(1/2 1/2  0 )</li>
<li>8=&gt; use time-reversal symmetry for k=( 0  1/2 1/2)</li>
<li>9=&gt; use time-reversal symmetry for k=(1/2 1/2 1/2)</li>
<li>0=&gt; (preprocessed) for each k point, choose automatically
the appropriate time-reversal option when it is allowed,
and chose <b>istwfk</b>=1 for all the other k points.</li>
</ul>
Note that the input variable "<a href="varfil.html#mkmem">mkmem</a>" also controls
the wavefunction storage, but at the
level of core memory versus disk space."""
},
'ixc': {
'definition': "Integer for eXchange-Correlation choice ",
'section': "varbas",
'category': " ",
'vartype': "integer parameter ",
'default': "<b>ixc</b>=1 (Teter parameterization). However, if all the pseudopotentials have the same value of pspxc, the initial value of <b>ixc</b> will be that common value.",
'text': """Controls the choice of exchange and correlation (xc). The list of XC functionals is given
below. Positive values are for ABINIT native library of XC functionals, while negative values are for calling
the much wider set of functionals from the ETSF LibXC library (by M. Marques), also available at the
<a href="http://www.etsf.eu/resources/software/libraries_and_tools"> ETSF library Web page </a>

<br>Note that the choice made here should agree with the choice
made in generating the original pseudopotential, except
for <b>ixc</b>=0 (usually only used for debugging).
A warning is issued if this is not the case.
However, the choices <b>ixc</b>=1, 2, 3 and 7 are fits to the same data, from
Ceperley-Alder, and are rather similar, at least for spin-unpolarized systems.
<br>The choice between the non-spin-polarized and spin-polarized case
is governed by the value of <a href="varbas.html#nsppol">nsppol</a> (see below).
<p>

<p><b> Native ABINIT XC functionals</b><p>

<br>NOTE : in the implementation of the spin-dependence of these
functionals, and in order to avoid divergences in their
derivatives, the interpolating function between spin-unpolarized
and fully-spin-polarized function has been slightly modified,
by including a zeta rescaled by 1.d0-1.d-6. This should affect
total energy at the level of 1.d-6Ha, and should
have an even smaller effect on differences of energies, or derivatives.
<br>The value <b>ixc</b>=10 is used internally : gives the difference between <b>ixc</b>=7 and
<b>ixc</b>=9, for use with an accurate RPA correlation energy.

<p>

<ul>
<li>0=> NO xc; </li>
</ul>
<ul>
<li>1=> LDA or LSD, Teter Pade parametrization (4/93, published in S. Goedecker, M. Teter, J. Huetter, Phys.Rev.B54, 1703 (1996)), which
reproduces Perdew-Wang (which reproduces Ceperley-Alder!).</li>
<li>2=> LDA, Perdew-Zunger-Ceperley-Alder (no spin-polarization)</li>
<li>3=> LDA, old Teter rational polynomial parametrization (4/91)
fit to Ceperley-Alder data (no spin-polarization)</li>
<li>4=> LDA, Wigner functional (no spin-polarization)</li>
<li>5=> LDA, Hedin-Lundqvist functional (no spin-polarization)</li>
<li>6=> LDA, "X-alpha" functional (no spin-polarization)</li>
<li>7=> LDA or LSD, Perdew-Wang 92 functional</li>
<li>8=> LDA or LSD, x-only part of the Perdew-Wang 92 functional</li>
<li>9=> LDA or LSD, x- and RPA correlation part of the Perdew-Wang 92 functional</li>
</ul>
<ul>
<li>11=> GGA, Perdew-Burke-Ernzerhof GGA functional</li>
<li>12=> GGA, x-only part of Perdew-Burke-Ernzerhof GGA functional</li>
<li>13=> GGA potential of van Leeuwen-Baerends, while for energy, Perdew-Wang 92 functional</li>
<li>14=> GGA, revPBE of Y. Zhang and W. Yang, Phys. Rev. Lett. 80, 890 (1998)</li>
<li>15=> GGA, RPBE of B. Hammer, L.B. Hansen and J.K. Norskov, Phys. Rev. B 59, 7413 (1999)</li>
<li>16=> GGA, HTCH93 of F.A. Hamprecht, A.J. Cohen, D.J. Tozer, N.C. Handy, J. Chem. Phys. 109, 6264 (1998)</li>
<li>17=> GGA, HTCH120 of A.D. Boese, N.L. Doltsinis, N.C. Handy, and M. Sprik, J. Chem. Phys 112, 1670 (1998) - The usual HCTH functional.</li>
<li>18=> (NOT AVAILABLE : used internally for GGA BLYP pseudopotentials from
M. Krack, see Theor. Chem. Acc. 114, 145 (2005),
available from the <a href="http://cvs.berlios.de/cgi-bin/viewcvs.cgi/cp2k/potentials/Goedecker/abinit/blyp/">CP2K repository</a>
- use the LibXC instead, with <b>ixc</b>=-106131.</li>
<li>19=> (NOT AVAILABLE : used internally for GGA BP86 pseudopotentials from
M. Krack, see Theor. Chem. Acc. 114, 145 (2005),
available from the <a href="http://cvs.berlios.de/cgi-bin/viewcvs.cgi/cp2k/potentials/Goedecker/abinit/bp/">CP2K repository</a>
- use the LibXC instead, with <b>ixc</b>=-106132.</li>
</ul>
<ul>
<li>20=> Fermi-Amaldi xc ( -1/N Hartree energy, where N is the
number of electrons per cell ; G=0 is not taken
into account however), for TDDFT tests.
No spin-pol. Does not work for RF.</li>
<li>21=> same as 20, except that the xc-kernel is the LDA (<b>ixc</b>=1) one,
for TDDFT tests.</li>
<li>22=> same as 20, except that the xc-kernel is the Burke-Petersilka-Gross
hybrid, for TDDFT tests.</li>
<li>23=> GGA of Z. Wu and R.E. Cohen, Phys. Rev. 73, 235116 (2006).</li>
<li>24=> GGA, C09x exchange of V. R. Cooper, PRB 81, 161104(R) (2010).</li>
<li>26=> GGA, HTCH147 of A.D. Boese, N.L. Doltsinis, N.C. Handy, and M. Sprik, J. Chem. Phys 112, 1670 (1998).</li>
<li>27=> GGA, HTCH407 of A.D. Boese, and N.C. Handy, J. Chem. Phys 114, 5497 (2001).</li>
<li>28=> (NOT AVAILABLE : used internally for GGA OLYP pseudopotentials from
M. Krack, see Theor. Chem. Acc. 114, 145 (2005),
available from the <a href="http://cvs.berlios.de/cgi-bin/viewcvs.cgi/cp2k/potentials/Goedecker/abinit/blyp/">CP2K repository</a>
- use the LibXC instead, with <b>ixc</b>=-110131.</li>
</ul>

<p><b> ETSF Lib XC functionals</b><p>
Note that you must compile ABINIT with the LibXC plug-in in order to be able to access these functionals.
<br>
The LibXC functionals are accessed by <b>negative values</b> of <b>ixc</b>.
The LibXC contains functional forms for either exchange-only functionals, correlation-only functionals,
or combined exchange and correlation functionals. Each of them is to be specified by a three-digit number.
In case of a combined exchange and correlation functional, only one such three-digit number has to be specified as value of <b>ixc</b>,
with a minus sign (to indicate that it comes from the LibXC).
In the case of separate exchange functional (let us represent its identifier by XXX) and
correlation functional (let us represent its identified by CCC),
a six-digit number will have to be specified for <b>ixc</b>, by concatenation, be it XXXCCC or CCCXXX.
As an example, <b>ixc</b>=-020 gives access to the Teter93 LDA, while
<b>ixc</b>=-101130 gives access to the PBE GGA.
In version 0.9 of LibXC (December 2008), there are 16 three-dimensional (S)LDA functionals (1 for X, 14 for C, 1 for combined XC),
and there are 41 three-dimensional GGA (23 for X, 8 for C, 10 for combined XC).
Note that for a meta-GGA, the kinetic energy density is needed. This means having <a href="vargs.html#usekden">usekden</a>=1 .
<p>
<p> (S)LDA functionals (do not forget to add a minus sign, as discussed above)<p>
<ul>

<li>001=> XC_LDA_X
[PAM Dirac, Proceedings of the Cambridge Philosophical Society 26, 376 (1930);
F Bloch, Zeitschrift fuer Physik 57, 545 (1929)
]</li>

<li>002=> XC_LDA_C_WIGNER Wigner parametrization
[EP Wigner, Trans. Faraday Soc. 34, 678 (1938)
]</li>

<li>003=> XC_LDA_C_RPA Random Phase Approximation
[M Gell-Mann and KA Brueckner, Phys. Rev. 106, 364 (1957)
]</li>

<li>004=> XC_LDA_C_HL Hedin & Lundqvist
[L Hedin and BI Lundqvist, J. Phys. C 4, 2064 (1971)
]</li>

<li>005=> XC_LDA_C_GL !  Gunnarson & Lundqvist
[O Gunnarsson and BI Lundqvist, PRB 13, 4274 (1976)
]</li>

<li>006=> XC_LDA_C_XALPHA  !  Slater's Xalpha
]</li>

<li>007=> XC_LDA_C_VWN !  Vosko, Wilk, & Nussair
[SH Vosko, L Wilk, and M Nusair, Can. J. Phys. 58, 1200 (1980)
]</li>

<li>008=> XC_LDA_C_VWN_RPA !  Vosko, Wilk, & Nussair (RPA)
[SH Vosko, L Wilk, and M Nusair, Can. J. Phys. 58, 1200 (1980)
]</li>

<li>009=> XC_LDA_C_PZ !  Perdew & Zunger
[Perdew and Zunger, Phys. Rev. B 23, 5048 (1981)
]</li>

<li>010=> XC_LDA_C_PZ_MOD !  Perdew & Zunger (Modified)
[Perdew and Zunger, Phys. Rev. B 23, 5048 (1981)
Modified to improve the matching between the low and high rs part
]</li>

<li>011=> XC_LDA_C_OB_PZ !  Ortiz & Ballone (PZ)
[G Ortiz and P Ballone, Phys. Rev. B 50, 1391 (1994) ;
G Ortiz and P Ballone, Phys. Rev. B 56, 9970(E) (1997) ;
Perdew and Zunger, Phys. Rev. B 23, 5048 (1981)
]</li>

<li>012=> XC_LDA_C_PW !  Perdew & Wang
[JP Perdew and Y Wang, Phys. Rev. B 45, 13244 (1992)
]</li>

<li>013=> XC_LDA_C_PW_MOD !  Perdew & Wang (Modified)
[JP Perdew and Y Wang, Phys. Rev. B 45, 13244 (1992) ;
Added extra digits to some constants as in the PBE routine
see <a href="http://www.chem.uci.edu/~kieron/dftold2/pbe.php">http://www.chem.uci.edu/~kieron/dftold2/pbe.php</a>
(at some point it was available at http://dft.uci.edu/pbe.php)
]</li>

<li>014=> XC_LDA_C_OB_PW !  Ortiz & Ballone (PW)
[G Ortiz and P Ballone, Phys. Rev. B 50, 1391 (1994) ;
G Ortiz and P Ballone, Phys. Rev. B 56, 9970(E) (1997) ;
JP Perdew and Y Wang, Phys. Rev. B 45, 13244 (1992)
]</li>

<li>017=> XC_LDA_C_vBH !  von Barth & Hedin
[U von Barth and L Hedin, J. Phys. C: Solid State Phys. 5, 1629 (1972)
]</li>

<li>020=> XC_LDA_XC_TETER93 !  Teter 93 parametrization
[S Goedecker, M Teter, J Hutter, PRB 54, 1703 (1996)
]</li>

<li>022=> XC_LDA_C_ML1 !  Modified LSD (version 1) of Proynov and Salahub
[EI Proynov and D Salahub, Phys. Rev. B 49, 7874 (1994)
]</li>

<li>023=> XC_LDA_C_ML2 !  Modified LSD (version 2) of Proynov and Salahub
[EI Proynov and D Salahub, Phys. Rev. B 49, 7874 (1994)
]</li>

<li>024=> XC_LDA_C_GOMBAS !  Gombas parametrization
[P. Gombas, Pseudopotentials (Springer-Verlag, New York, 1967)
]</li>


</ul>
<p> GGA functionals (do not forget to add a minus sign, as discussed above)<p>
<ul>

<li>101=> XC_GGA_X_PBE !  Perdew, Burke & Ernzerhof exchange
[JP Perdew, K Burke, and M Ernzerhof, Phys. Rev. Lett. 77, 3865 (1996) ;
JP Perdew, K Burke, and M Ernzerhof, Phys. Rev. Lett. 78, 1396(E) (1997)
]</li>

<li>102=> XC_GGA_X_PBE_R !  Perdew, Burke & Ernzerhof exchange (revised)
[Y Zhang and W Yang, Phys. Rev. Lett 80, 890 (1998)
]</li>

<li>103=> XC_GGA_X_B86 !  Becke 86 Xalfa,beta,gamma
[AD Becke, J. Chem. Phys 84, 4524 (1986)
]</li>

<li>104=> XC_GGA_X_HERMAN !  Herman Xalphabeta GGA
[F Herman, JP Van Dyke, and IB Ortenburger, Phys. Rev. Lett. 22, 807 (1969) ;
F Herman, IB Ortenburger, and JP Van Dyke, Int. J. Quantum Chem. Symp. 3, 827 (1970)
]</li>

<li>105=> XC_GGA_X_B86_MGC !  Becke 86 Xalfa,beta,gamma (with mod. grad. correction)
[AD Becke, J. Chem. Phys 84, 4524 (1986) ;
AD Becke, J. Chem. Phys 85, 7184 (1986)
]</li>

<li>106=> XC_GGA_X_B88 !  Becke 88
[AD Becke, Phys. Rev. A 38, 3098 (1988)
]</li>

<li>107=> XC_GGA_X_G96 !  Gill 96
[PMW Gill, Mol. Phys. 89, 433 (1996)
]</li>

<li>108=> XC_GGA_X_PW86 !  Perdew & Wang 86
[JP Perdew and Y Wang, Phys. Rev. B 33, 8800 (1986)
]</li>

<li>109=> XC_GGA_X_PW91 !  Perdew & Wang 91
[JP Perdew, in Proceedings of the 21st Annual International Symposium on the Electronic Structure of Solids, ed. by P Ziesche and H
Eschrig (Akademie Verlag, Berlin, 1991), p. 11. ;
JP Perdew, JA Chevary, SH Vosko, KA Jackson, MR Pederson, DJ Singh, and C Fiolhais, Phys. Rev. B 46, 6671 (1992) ;
JP Perdew, JA Chevary, SH Vosko, KA Jackson, MR Pederson, DJ Singh, and C Fiolhais, Phys. Rev. B 48, 4978(E) (1993)
]</li>

<li>110=> XC_GGA_X_OPTX !  Handy & Cohen OPTX 01
[NC Handy and AJ Cohen, Mol. Phys. 99, 403 (2001)
]</li>

<li>111=> XC_GGA_X_DK87_R1 !  dePristo & Kress 87 (version R1)
[AE DePristo and JD Kress, J. Chem. Phys. 86, 1425 (1987)
]</li>

<li>112=> XC_GGA_X_DK87_R2 !  dePristo & Kress 87 (version R2)
[AE DePristo and JD Kress, J. Chem. Phys. 86, 1425 (1987)
]</li>

<li>113=> XC_GGA_X_LG93 !  Lacks & Gordon 93
[DJ Lacks and RG Gordon, Phys. Rev. A 47, 4681 (1993)
]</li>

<li>114=> XC_GGA_X_FT97_A !  Filatov & Thiel 97 (version A)
[M Filatov and W Thiel, Mol. Phys 91, 847 (1997)
]</li>

<li>115=> XC_GGA_X_FT97_B !  Filatov & Thiel 97 (version B)
[M Filatov and W Thiel, Mol. Phys 91, 847 (1997)
]</li>

<li>116=> XC_GGA_X_PBE_SOL !  Perdew, Burke & Ernzerhof exchange (solids)
[JP Perdew, et al, Phys. Rev. Lett. 100, 136406 (2008)
]</li>

<li>117=> XC_GGA_X_RPBE !  Hammer, Hansen & Norskov (PBE-like)
[B Hammer, LB Hansen and JK Norskov, Phys. Rev. B 59, 7413 (1999)
]</li>

<li>118=> XC_GGA_X_WC !  Wu & Cohen
[Z Wu and RE Cohen, Phys. Rev. B 73, 235116 (2006)
]</li>

<li>119=> XC_GGA_X_mPW91 !  Modified form of PW91 by Adamo & Barone
[C Adamo and V Barone, J. Chem. Phys. 108, 664 (1998)
]</li>

<li>120=> XC_GGA_X_AM05 !  Armiento & Mattsson 05 exchange
[R Armiento and AE Mattsson, Phys. Rev. B 72, 085108 (2005) ;
AE Mattsson, R Armiento, J Paier, G Kresse, JM Wills, and TR Mattsson, J. Chem. Phys. 128, 084714 (2008)
]</li>

<li>121=> XC_GGA_X_PBEA !  Madsen (PBE-like)
[G Madsen, Phys. Rev. B 75, 195108 (2007)
]</li>

<li>122=> XC_GGA_X_MPBE !  Adamo & Barone modification to PBE
[C Adamo and V Barone, J. Chem. Phys. 116, 5933 (2002)
]</li>

<li>123=> XC_GGA_X_XPBE !  xPBE reparametrization by Xu & Goddard
[X Xu and WA Goddard III, J. Chem. Phys. 121, 4068 (2004)
]</li>

<li>125=> XC_GGA_X_BAYESIAN !  Bayesian best fit for the enhancement factor
[JJ Mortensen, K Kaasbjerg, SL Frederiksen, JK Norskov, JP Sethna, and KW Jacobsen, Phys. Rev. Lett. 95, 216401 (2005)
]</li>

<li>126=> XC_GGA_X_PBE_JSJR !  PBE JSJR reparametrization by Pedroza, Silva & Capelle
[LS Pedroza, AJR da Silva, and K. Capelle, Phys. Rev. B 79, 201106(R) (2009)
]</li>

<li>130=> XC_GGA_C_PBE !  Perdew, Burke & Ernzerhof correlation
[JP Perdew, K Burke, and M Ernzerhof, Phys. Rev. Lett. 77, 3865 (1996) ;
JP Perdew, K Burke, and M Ernzerhof, Phys. Rev. Lett. 78, 1396(E) (1997)
]</li>

<li>131=> XC_GGA_C_LYP !  Lee, Yang & Parr
[C Lee, W Yang and RG Parr, Phys. Rev. B 37, 785 (1988)
B Miehlich, A Savin, H Stoll and H Preuss, Chem. Phys. Lett. 157, 200 (1989)
]</li>

<li>132=> XC_GGA_C_P86 !  Perdew 86
[JP Perdew, Phys. Rev. B 33, 8822 (1986)
]</li>

<li>133=> XC_GGA_C_PBE_SOL !  Perdew, Burke & Ernzerhof correlation SOL
[JP Perdew, et al, Phys. Rev. Lett. 100, 136406 (2008)
]</li>

<li>134=> XC_GGA_C_PW91 !  Perdew & Wang 91
[JP Perdew, JA Chevary, SH Vosko, KA Jackson, MR Pederson, DJ Singh, and C Fiolhais, Phys. Rev. B 46, 6671 (1992)
]</li>

<li>135=> XC_GGA_C_AM05 !  Armiento & Mattsson 05 correlation
[ R Armiento and AE Mattsson, Phys. Rev. B 72, 085108 (2005) ;
AE Mattsson, R Armiento, J Paier, G Kresse, JM Wills, and TR Mattsson, J. Chem. Phys. 128, 084714 (2008)
]</li>

<li>136=> XC_GGA_C_XPBE !  xPBE reparametrization by Xu & Goddard
[X Xu and WA Goddard III, J. Chem. Phys. 121, 4068 (2004)
]</li>

<li>137=> XC_GGA_C_LM !  Langreth and Mehl correlation
[DC Langreth and MJ Mehl, Phys. Rev. Lett. 47, 446 (1981)
]</li>

<li>138=> XC_GGA_C_PBE_JRGX !  JRGX reparametrization by Pedroza, Silva & Capelle
[LS Pedroza, AJR da Silva, and K. Capelle, Phys. Rev. B 79, 201106(R) (2009)
]</li>

<li>139=> XC_GGA_X_OPTB88_VDW !  Becke 88 reoptimized to be used with vdW functional of Dion et al
[J Klimes, DR Bowler, and A Michaelides, J. Phys.: Condens. Matter 22, 022201 (2010)
]</li>

<li>140=> XC_GGA_X_PBEK1_VDW !  PBE reparametrization for vdW
[J Klimes, DR Bowler, and A Michaelides, J. Phys.: Condens. Matter 22, 022201 (2010)
]</li>

<li>141=> XC_GGA_X_OPTPBE_VDW !  PBE reparametrization for vdW
[J Klimes, DR Bowler, and A Michaelides, J. Phys.: Condens. Matter 22, 022201 (2010)
]</li>

<li>142=> XC_GGA_X_RGE2 !  Regularized PBE
[A Ruzsinszky, GI Csonka, and G Scuseria, J. Chem. Theory Comput. 5, 763 (2009)
]</li>

<li>143=> XC_GGA_C_RGE2 !  Regularized PBE
[A Ruzsinszky, GI Csonka, and G Scuseria, J. Chem. Theory Comput. 5, 763 (2009)
]</li>

<li>144=> XC_GGA_X_RPW86 !  refitted Perdew & Wang 86
[ED Murray, K Lee and DC Langreth, J. Chem. Theory Comput. 5, 2754-2762 (2009)
]</li>

<li>145=> XC_GGA_X_KT1 !  Keal and Tozer version 1
[TW Keal and DJ Tozer, J. Chem. Phys. 119, 3015 (2003)
]</li>

<li>146=> XC_GGA_XC_KT2 !  Keal and Tozer version 2
[TW Keal and DJ Tozer, J. Chem. Phys. 119, 3015 (2003)
]</li>

<li>147=> XC_GGA_C_WL !  Wilson & Levy
[LC Wilson and M Levy, Phys. Rev. B 41, 12930 (1990)
]</li>

<li>148=> XC_GGA_C_WI !  Wilson & Ivanov
[LC Wilson & S Ivanov, Int. J. Quantum Chem. 69, 523-532 (1998)
]</li>

<li>149=> XC_GGA_X_MB88 !  Modified Becke 88 for proton transfer
[V Tognetti and C Adamo, J. Phys. Chem. A 113, 14415-14419 (2009)
]</li>

<!-- The following functional is untested in ABINIT
<li>160=> XC_GGA_XC_LB !  van Leeuwen & Baerends
[R van Leeuwen and EJ Baerends, Phys. Rev. A. 49, 2421 (1994)
]</li>
-->

<li>161=> XC_GGA_XC_HCTH_93 !  HCTH functional fitted to  93 molecules
[FA Hamprecht, AJ Cohen, DJ Tozer, and NC Handy, J. Chem. Phys. 109, 6264 (1998)
]</li>

<li>162=> XC_GGA_XC_HCTH_120 !  HCTH functional fitted to 120 molecules
[AD Boese, NL Doltsinis, NC Handy, and M Sprik, J. Chem. Phys. 112, 1670 (2000)
]</li>

<li>163=> XC_GGA_XC_HCTH_147 !  HCTH functional fitted to 147 molecules
[AD Boese, NL Doltsinis, NC Handy, and M Sprik, J. Chem. Phys. 112, 1670 (2000)
]</li>

<li>164=> XC_GGA_XC_HCTH_407 !  HCTH functional fitted to 407 molecules
[AD Boese, and NC Handy, J. Chem. Phys. 114, 5497 (2001)
]</li>

<li>165=> XC_GGA_XC_EDF1 !  Empirical functionals from Adamson, Gill, and Pople
[RD Adamson, PMW Gill, and JA Pople, Chem. Phys. Lett. 284 6 (1998)
]</li>

<li>166=> XC_GGA_XC_XLYP !  XLYP functional
[X Xu and WA Goddard, III, PNAS 101, 2673 (2004)
]</li>

<li>167=> XC_GGA_XC_B97 !  Becke 97
[AD Becke, J. Chem. Phys. 107, 8554-8560 (1997)
]</li>

<li>168=> XC_GGA_XC_B97_1 !  Becke 97-1
[FA Hamprecht, AJ Cohen, DJ Tozer, and NC Handy, J. Chem. Phys. 109, 6264 (1998);
AD Becke, J. Chem. Phys. 107, 8554-8560 (1997)
]</li>

<li>169=> XC_GGA_XC_B97_2 !  Becke 97-2
[AD Becke, J. Chem. Phys. 107, 8554-8560 (1997)
]</li>

<li>170=> XC_GGA_XC_B97_D !  Grimme functional to be used with C6 vdW term
[S Grimme, J. Comput. Chem. 27, 1787 (2006)
]</li>

<li>171=> XC_GGA_XC_B97_K !  Boese-Martin for Kinetics
[AD Boese and JML Martin, J. Chem. Phys., Vol. 121, 3405 (2004)
]</li>

<li>172=> XC_GGA_XC_B97_3 !  Becke 97-3
[TW Keal and DJ Tozer, J. Chem. Phys. 123, 121103 (2005)
]</li>

<li>173=> XC_GGA_XC_PBE1W !  Functionals fitted for water
[EE Dahlke and DG Truhlar, J. Phys. Chem. B 109, 15677 (2005)
]</li>

<li>174=> XC_GGA_XC_MPWLYP1W !  Functionals fitted for water
[EE Dahlke and DG Truhlar, J. Phys. Chem. B 109, 15677 (2005)
]</li>

<li>175=> XC_GGA_XC_PBELYP1W !  Functionals fitted for water
[EE Dahlke and DG Truhlar, J. Phys. Chem. B 109, 15677 (2005)
]</li>

<li>176=> XC_GGA_XC_SB98_1a !  Schmider-Becke 98 parameterization 1a
[HL Schmider and AD Becke, J. Chem. Phys. 108, 9624 (1998)
]</li>

<li>177=> XC_GGA_XC_SB98_1b !  Schmider-Becke 98 parameterization 1b
[HL Schmider and AD Becke, J. Chem. Phys. 108, 9624 (1998)
]</li>

<li>178=> XC_GGA_XC_SB98_1c !  Schmider-Becke 98 parameterization 1c
[HL Schmider and AD Becke, J. Chem. Phys. 108, 9624 (1998)
]</li>

<li>179=> XC_GGA_XC_SB98_2a !  Schmider-Becke 98 parameterization 2a
[HL Schmider and AD Becke, J. Chem. Phys. 108, 9624 (1998)
]</li>

<li>180=> XC_GGA_XC_SB98_2b !  Schmider-Becke 98 parameterization 2b
[HL Schmider and AD Becke, J. Chem. Phys. 108, 9624 (1998)
]</li>

<li>181=> XC_GGA_XC_SB98_2c !  Schmider-Becke 98 parameterization 2c
[HL Schmider and AD Becke, J. Chem. Phys. 108, 9624 (1998)
]</li>

<li>183=> XC_GGA_X_OL2 !  Exchange form based on Ou-Yang and Levy v.2
[P Fuentealba and O Reyes, Chem. Phys. Lett. 232, 31-34 (1995) ; H Ou-Yang, M Levy, Int. J. of Quant. Chem. 40, 379-388 (1991)
]</li>

<li>184=> XC_GGA_X_APBE !  mu fixed from the semiclassical neutral atom
[LA Constantin, E Fabiano, S Laricchia, and F Della Sala, Phys. Rev. Lett. 106, 186406 (2011)
]</li>

<li>186=> XC_GGA_C_APBE !  mu fixed from the semiclassical neutral atom
[LA Constantin, E Fabiano, S Laricchia, and F Della Sala, Phys. Rev. Lett. 106, 186406 (2011)
]</li>

<li>202=> XC_MGGA_X_TPSS !  Tao, Perdew, Staroverov & Scuseria
[J Tao, JP Perdew, VN Staroverov, and G Scuseria, Phys. Rev. Lett. 91, 146401 (2003) ;
JP Perdew, J Tao, VN Staroverov, and G Scuseria, J. Chem. Phys. 120, 6898 (2004)
]</li>

<li>203=> XC_MGGA_X_M06L !  Zhao, Truhlar exchange
[Y Zhao and DG Truhlar, JCP 125, 194101 (2006);
Y Zhao and DG Truhlar, Theor. Chem. Account 120, 215 (2008)
]</li>

<li>204=> XC_MGGA_X_GVT4 !  GVT4 (X part of VSXC) from van Voorhis and Scuseria
[T Van Voorhis and GE Scuseria, JCP 109, 400 (1998)
]</li>

<li>205=> XC_MGGA_X_TAU_HCTH ! tau-HCTH from Boese and Handy
[AD Boese and NC Handy, JCP 116, 9559 (2002)
]</li>

<li>207=> XC_MGGA_X_BJ06 !  Becke & Johnson correction to Becke-Roussel 89
[AD Becke and ER Johnson, J. Chem. Phys. 124, 221101 (2006)
] WARNING : this Vxc-only mGGA can only be used with a LDA correlation, typically Perdew-Wang 92.</li>

<li>208=> XC_MGGA_X_TB09 !  Tran-blaha - correction to Becke & Johnson correction to Becke-Roussel 89
[F Tran and P Blaha, Phys. Rev. Lett. 102, 226401 (2009)
] WARNING : this Vxc-only mGGA can only be used with a LDA correlation, typically Perdew-Wang 92.</li>

<li>209=> XC_MGGA_X_RPP09 !  Rasanen, Pittalis, and Proetto correction to Becke & Johnson
[E Rasanen, S Pittalis & C Proetto, arXiv:0909.1477 (2009)
] WARNING : this Vxc-only mGGA can only be used with a LDA correlation, typically Perdew-Wang 92.</li>

<li>232=> XC_MGGA_C_VSXC !  VSxc from Van Voorhis and Scuseria (correlation part)
[T Van Voorhis and GE Scuseria, JCP 109, 400 (1998)
]</li>

</ul>


"""
},
'ixcpositron': {
'definition': "Integer for the eXchange-Correlation applied to the electron-POSITRON interaction",
'section': "vargs",
'category': " ",
'vartype': "integer parameter",
'default': "0",
'text': """Relevant only when <a href="vargs.html#positron">positron</a>/=0.<br>
Define the type of electron-positron correlation that is used in case
of a electron-positron two-component DFT calculation.<br>
Define also the analytical formula of the enhancement factor used to compute the electron-positron annhilation rate:<br><br>
Electron-positron correlation functional:<br>
<ul>
<b>ixcpositron=1</b>:
LDA zero positron density limit parametrized by Arponen & Pajanne and provided by Boronski & Nieminen [1,2]<br>
<b>ixcpositron=11</b>:
LDA zero positron density limit parametrized by Arponen & Pajanne and fitted by Sterne & Kaiser [1,3]<br>
<b>ixcpositron=2</b>:
LDA electron-positron correlation provided by Puska, Seitsonen, and Nieminen [1,4]<br>
<b>ixcpositron=3</b>:
GGA zero positron density limit parametrized by Arponen & Pajanne and provided by Boronski & Nieminen [1,2,5]<br>
<b>ixcpositron=31</b>:
GGA zero positron density limit parametrized by Arponen & Pajanne and fitted by Sterne & Kaiser [1,3,5]
</ul>
Annihilation rate enhancement factor:<br>
<ul>
<b>ixcpositron=1</b>:
Boronski and Nieminen full modelisation and RPA limit [1]<br>
<b>ixcpositron=11</b>:
Sterne and Kaiser [2]<br>
<b>ixcpositron=2</b>:
Puska, Seitsonen and Nieminen [3]<br>
<b>ixcpositron=3</b>:
Boronski and Nieminen full modelisation and RPA limit [1], with GGA corrections<br>
<b>ixcpositron=31</b>:
Sterne and Kaiser [2], with GGA corrections<br>
</ul>
References:
<ul>
<b>[1]</b> J. Arponen and E. Pajanne, Ann. Phys. (N.Y.) 121, 343 (1979).<br>
<b>[2]</b> Boronski and R.M. Nieminen, Phys. Rev. B 34, 3820 (1986).<br>
<b>[3]</b> P.A. Sterne and J.H. Kaiser, Phys. Rev. B 43, 13892 (1991).<br>
<b>[4]</b> M.J. Puska, A.P. Seitsonen and R.M. Nieminen, Phys. Rev. B 52, 10947 (1994).<br>
<b>[5]</b> B. Barbiellini, M.J. Puska, T. Torsti and R.M.Nieminen, Phys. Rev. B 51, 7341 (1994)<br>
</ul>
"""
},
'jdtset': {
'definition': "index -J- for DaTaSETs ",
'section': "varbas",
'category': """<a href="../users/abinit_help.html#no_multi">NO MULTI</a> """,
'vartype': """integer array <b>jdtset</b>(<a href="varbas.html#ndtset">ndtset</a>)""",
'default': """the series 1, 2, 3 ... <a href="varbas.html#ndtset">ndtset</a> .""",
'text': """Gives the dataset index
of each of the datasets. This index will be used :
<ul>
<li>to determine which input variables are specific to each
dataset, since the variable names for this
dataset will be made from the bare variable
name concatenated with this index, and only if
such a composite variable name does not exist,
the code will consider the bare variable name,
or even, the Default;</li>
<li>to characterize output variable names, if their
content differs from dataset to dataset;</li>
<li>to characterize output files ( root names appended with _DSx
where 'x' is the dataset index ).</li>
</ul>
The allowed index values are between 1 and 9999.
<br>An input variable name appended with 0 is not allowed.
<br>When <a href="varbas.html#ndtset">ndtset</a>==0, this array is not used, and moreover,
no input variable name appended with a digit is allowed.
This array might be initialized thanks to the use of
the input variable <a href="varbas.html#udtset">udtset</a>. In this case, <b>jdtset</b> cannot
be used."""
},
'jellslab': {
'definition': "include a JELLium SLAB in the cell",
'section': "vargs",
'category': " ",
'vartype': "integer parameter   ",
'default': "0 (no jellium slab).",
'text': """If set to 1, a slab of uniform positive background charge density,
that is, a jellium slab, is included in the calculation cell.
A portion of the unit cell is filled with such positive charge density distribution
which is equal to a bulk-mean value n<sub>bulk</sub> between two edges
and zero in the vacuum region if present.
<br>
For the sake of convenience the unit cell is supposed
to have the third crystal primitive lattice vector orthogonal
to the other ones so that the portion of the cell filled by the jellium slab can be defined through its edges along z.
<br>
The bulk-mean positive charge density is fixed by the input variable <a href="vargs.html#slabwsrad">slabwsrad</a>,
while the position of the slab edges along z is defined through
the input variables <a href="vargs.html#slabzbeg">slabzbeg</a> and <a href="vargs.html#slabzend">slabzend</a>."""
},
'jfielddir': {
'definition': "electric/displacement FIELD DIRection ",
'section': "varff",
'category': " ",
'vartype': "integer array jfielddir(3)   ",
'default': "3*0 .",
'text': """If set to 1, a slab of uniform positive background charge density,
When specifying mixed electric field boundary conditions ( <a href="varff.html#berryopt">berryopt</a>=17),
jfielddir controls whether reduced electric field (<b>jfielddir</b>=1)
or reduced electric displacement field (<b>jfielddir</b>=2) is chosen to be fixed,
in each of the three lattice directions (i.e., in the reduced, not the Cartesian, frame).
For example, <b>jfielddir</b>=(1 1 2) tells the code to use fixed ebar_1 and ebar_2
along the first two lattice directions and fixed d_3 along the third.
<br>
For the case of mixed electric field boundary conditions,
<a href="varff.html#red_efieldbar">red_efieldbar</a> and <a href="varff.html#red_dfield">red_dfield</a>
are used to control ebar and d, respectively.
For example, for electric boundary conditions corresponding to a material
in a parallel-plate capacitor, if you want to control d_3=d0,
while fixing ebar_1=ebar_2=0, then the input files should have
<a href="varff.html#berryopt">berryopt</a>=17, <b>jfielddir</b>=(1 1 2), <a href="varff.html#red_efieldbar">red_efieldbar</a>=(0.0 0.0 a),
and <a href="varff.html#red_dfield">red_dfield</a>=(b c d0).
Here a, b, and c are the starting values.
They can be chosen in this way: do a single run for fixed d calculation (<a href="varff.html#red_dfield">red_dfield</a>=0,0,d0),
from the final results you will have ebar_3, which is a good guess for a.
Then do another single run for fixed ebar calculation (<a href="varff.html#red_efieldbar">red_efieldbar</a>=(0 0 0)),
from the final results you will have d_1,d_2, these are good guesses for b, c.
"""
},
'jpawu': {
'definition': "value of J for PAW+U ",
'section': "varpaw",
'category': """<a href="../users/abinit_help.html#energy">ENERGY</a>""",
'vartype': """real array jpawu(<a href="varbas.html#ntypat">ntypat</a>) """,
'default': "0",
'text': """Relevant only when
<a href="varint.html#usepaw">usepaw</a>=1, and <a
href="varpaw.html#usepawu">usepawu</a>=1.<br>
<br>
Gives the value of the screened exchange interaction between correlated
electrons corresponding to <a href="varpaw.html#lpawu">lpawu</a>
for each species.<br>
In the case where <a href="varpaw.html#lpawu">lpawu</a> =-1,
the value is not used."""
},
'kberry': {
'definition': "K wavevectors for BERRY phase computation ",
'section': "varff",
'category': " ",
'vartype': """integer array kberry(3,<a href="varff.html#nberry">nberry</a>)  """,
'default': "an array of 0",
'text': """<p>Used for non-zero values of <a href="varff.html#berryopt">berryopt</a>.
<p>This array defines, for each Berry phase calculation
(the number of such calculations is defined by
<a href="varff.html#nberry">nberry</a>), the
difference of wavevector between k points for which
the overlap matrix must be computed.
The polarisation vector will be projected
on the direction of that wavevector,
and the result of the computation will be the magnitude of this
projection.
Doing more than one wavevector, with different independent
direction, allows to find the full polarisation vector.
However, note that converged results need oriented grids,
denser along the difference wavevector than usual Monkhorst-Pack
grids.
<p> The difference of wavevector is computed in the coordinate
system defined by the k-points grid
(see <a href="varbas.html#ngkpt">ngkpt</a>
and <a href="vargs.html#kptrlatt">kptrlatt</a>), so that
the values of <b>kberry</b> are integers.
Of course, such a k point grid must exist, and all the
corresponding wavefunctions must be available, so that the
computation is allowed only when <a href="varbas.html#kptopt">kptopt</a>
is equal to 3. In order to save computing time, it is suggested
to make a preliminary calculation of the wavefunctions on the
irreducible part of the grid, with <a href="varbas.html#kptopt">kptopt</a>
equal to 1, and then use these converged wavefunctions
in the entire Brillouin zone, by reading them to initialize
the <a href="varbas.html#kptopt">kptopt</a>=3 computation.
"""
},
'kpt': {
'definition': "K - PoinTs ",
'section': "varbas",
'category': " ",
'vartype': """real array kpt(3,<a href="varbas.html#nkpt">nkpt</a>)  """,
'default': "0. 0. 0. (for just one k-point, adequate for one molecule in a supercell)",
'text': """Contains the k points in terms
of reciprocal space primitive translations (NOT in
cartesian coordinates!).
<br> Needed ONLY
if <a href="varbas.html#kptopt">kptopt</a>=0, otherwise
deduced from other input variables.

<p>It contains dimensionless numbers in terms of which
the cartesian coordinates would be:
<br><tele>k_cartesian = k1*G1+k2*G2+k3*G3 </tele>
<br>where <tele>(k1,k2,k3)</tele> represent the dimensionless "reduced
coordinates" and <tele>G1, G2, G3</tele> are the cartesian coordinates
of the primitive translation vectors.  G1,G2,G3 are related
to the choice of direct space primitive translation vectors
made in <a href="varbas.html#rprim">rprim</a>.
Note that an overall norm for the k
points is supplied by <a href="varbas.html#kptnrm">kptnrm</a>.  This allows
one to avoid supplying many digits for the k points to
represent such points as (1,1,1)/3.
<br>Note: one of the algorithms used to set up the sphere
of G vectors for the basis needs components of k-points
in the range [-1,1], so the
remapping is easily done by adding or subtracting 1 from
each component until it is in the range [-1,1].  That is,
given the k point normalization <a href="varbas.html#kptnrm">kptnrm</a> described below,
each component must lie in [-<a href="varbas.html#kptnrm">kptnrm</a>,<a href="varbas.html#kptnrm">kptnrm</a>].
<br>Note: a global shift can be provided by <a href="varint.html#qptn">qptn</a>
<br>Not read if <a href="varbas.html#kptopt">kptopt</a>/=0 .
"""
},
'kptbounds': {
'definition': "K PoinTs BOUNDarieS ",
'section': "vargs",
'category': """<a href="../users/abinit_help.html#not_internal">NOT INTERNAL</a> """,
'vartype': """real array kptbounds(3,abs(<a href="varbas.html#kptopt">kptopt</a>)+1)  """,
'default': "No Default",
'text': """It is used to generate the circuit to be followed by the band structure,
when <a href="varbas.html#kptopt">kptopt</a> is negative (it is
not read if <a href="varbas.html#kptopt">kptopt</a> is zero or positive).

<p>There are abs(<a href="varbas.html#kptopt">kptopt</a>)
segments to be defined, each of which starting from
the end point of the preceeding one. Thus,
the number of points to be input is
abs(<a href="varbas.html#kptopt">kptopt</a>)+1.
They form a circuit starting
at <b>kptbounds</b>(1:3,1)/<a href="varbas.html#kptnrm">kptnrm</a>
and ending at
<b>kptbounds</b>(1:3,abs(<a href="varbas.html#kptopt">kptopt</a>)+1)/<a href="varbas.html#kptnrm">kptnrm</a>.
The number of divisions of each segment can be defined either using the array <a href="vargs.html#ndivk">ndivk</a>
or the variable <a href="vargs.html#ndivsm">ndivsm</a> that just defines the number of divisions for the smallest segment
<p>As for <a href="varbas.html#kpt">kpt</a>, <b>kptbounds</b> is specified
using the primitive vectors in reciprocal space. If your Bravais lattice is simple,
then it should be quite easy to find the coordinates of the end points.
On the other hand, for centered, body-centered, face-centered, hexagonal, and rhombohedral
Bravais lattice,
the conversion might be more difficult. See the description of <a href="varbas.html#kpt">kpt</a>
for an explanation of how to convert data from the "conventional" cartesian coordinates to
the primitive vectors in the reciprocal space. In order to help a bit, we list below a series
of typical values, for the FCC, BCC, hexagonal and rhombohedral Bravais lattices. Note : all the data below
are given in dimensionless units ; they have to be rescaled by the actual lengths defined by the
<a href="varbas.html#acell">acell</a> values. However, <b>kptbounds</b> values can be used as such,
if the values of <a href="varbas.html#rprim">rprim</a> given below are adopted.

<p>A. <b>FCC lattice</b>
<p>Suppose the primitive vectors in real space are given by
<br><pre>
rprim   0 1 1    1 0 1    1 1 0
</pre>
or
<pre>
rprim   0 1/2 1/2    1/2 0 1/2    1/2 1/2 0
</pre>
(these two possibilities only differ by a scaling factor, irrelevant for the definition
of the k points in the primitive vectors in reciprocal space).

Then, the reciprocal primitive vectors (in conventional cartesian coordinates) are
<pre>
(-1/2 1/2 1/2), (1/2 -1/2 1/2), (1/2 1/2 -1/2)
</pre>
or
<pre>
(-1 1 1), (1 -1 1), (1 1 -1)
</pre>
and, in both cases, the coordinates of several special points with respect to primitive vectors in reciprocal space are
<pre>
X (0   1/2 1/2)   (conventional cartesian coordinate 1/2 0 0)
X'(1/2 1/2 1  )   (conventional cartesian coordinate 1/2 1/2 0)  (an other instance of X, in another Brillouin zone)
L (1/2 1/2 1/2)   (conventional cartesian coordinate  1/4 1/4 1/4)
L'(1/2 0   0  )   (conventional cartesian coordinate -1/4 1/4 1/4) (an other instance of L, on another face of the BZ)
W (1/4 1/2 3/4)   (conventional cartesian coordinate 1/2 1/4 0)
U (1/4 5/8 5/8)   (conventional cartesian coordinate 1/2 1/8 1/8)
K (3/8 3/8 3/4)   (conventional cartesian coordinate 3/8 3/8 0)
</pre>
Note that K is actually equivalent to U, by spatial and translational symmetry.
So, if you want to specify a typical circuit, the following might do the work :
L-Gamma-X-W-K,U-L-W-X-K,U-Gamma with
<br>
<pre>kptbounds  1/2 0 0  0 0 0  0 1/2 1/2  1/4 1/2 3/4  3/8 3/8 3/4  1/2 1/2 1/2  1/4 1/2 3/4  1/2 1/2 1  3/8 3/8 3/4  0 0 0</pre>
<p>
The lengths of segments (this information is useful to draw the band structure, with the correct relative
scale between special points)
can be found using
the conventional cartesian coordinates :
l(L-Gamma)=sqrt(3)/4=0.433... ;
l(Gamma-X)=1/2=0.5 ;
l(X-W)=1/4=0.25 ;
l(W-K)=sqrt(2)/8=0.177... ;
l(K-L)=sqrt(6)/8=0.306... ;
l(L-W)=sqrt(2)/4=0.354... ;
l(W-X)=1/4=0.25 ;
l(X-K)=sqrt(2)/8=0.177... ;
l(K-Gamma)=sqrt(2).3/8=0.530...
<p>

<p>B. <b>BCC lattice</b>
<p>Suppose the primitive vectors in real space are given by
<br><pre>
rprim  -1 1 1    1 -1 1    1 1 -1
</pre>
(as for the FCC lattice, there is a scale invariance).
Then, the reciprocal primitive vectors (in conventional cartesian coordinates) are
(0 1/2 1/2), (1/2 0 1/2), and (1/2 1/2 0)
and the coordinates of several special points with respect to primitive vectors in reciprocal space are
<pre>
H (-1/2 1/2 1/2)   (conventional cartesian coordinate 1/2 0 0)
N ( 0   0   1/2)   (conventional cartesian coordinate 1/4 1/4 0)
P ( 1/4 1/4 1/4)   (conventional cartesian coordinate 1/4 1/4 1/4)
</pre>
So, if you want to specify a typical circuit, the following might do the work :
Gamma-H-N-Gamma-P-N-P-H
<br>
<pre>kptbounds  0 0 0  -1/2 1/2 1/2  0 0 1/2  0 0 0   1/4 1/4 1/4  0 0 1/2  1/4 1/4 1/4  -1/2 1/2 1/2 </pre>
<p>
The lengths of segments (this information is useful to draw the band structure, with the correct relative scale between special points)
can be found using the conventional cartesian coordinates :
l(Gamma-H)=1/2=0.5 ;
l(H-N)=sqrt(2)/4=0.354... ;
l(N-Gamma)=sqrt(2)/4=0.354... ;
l(Gamma-P)=sqrt(3)/4=0.433... ;
l(P-N)=1/4=0.25 ;
l(N-P)=1/4=0.25 ;
l(P-H)=sqrt(3)/4=0.433...
<p>

<p>C. <b>Hexagonal lattices</b>
<p>Suppose the primitive vectors in real space are given by
<br><pre>
rprim  1 0 0    -1/2 sqrt(0.75) 0    0 0 1
</pre>
The coordinates of several special points with respect to primitive vectors in reciprocal space are
<pre>
M (1/2 0 0) or (0 1/2 0) or (-1/2 1/2 0)
L (1/2 0 1/2) or (0 1/2 1/2) or (-1/2 1/2 1/2)
K (1/3 1/3 0) or (2/3 -1/3 0) or (-1/3 2/3 0)
H (1/3 1/3 1/2) or (2/3 -1/3 1/2) or (-1/3 2/3 1/2)
A (0 0 1/2)
</pre>
So, if you want to specify a typical circuit, the following might do the work :
K-Gamma-M-K-H-A-L-H-L-M-Gamma-A
<br>
<pre>kptbounds  1/3 1/3 0  0 0 0  1/2 0 0  1/3 1/3 0  1/3 1/3 1/2  0 0 1/2  1/2 0 1/2  1/3 1/3 1/2  1/2 0 1/2  1/2 0 0  0 0 0  0 0 1/2 </pre>
<p>
In order to find the lengths of segments
(this information is useful to draw the band structure,
with the correct relative scale between special points)
one needs to know the a and c lattice parameters. Also, in what follows, we omit the 2*pi factor sometimes
present in the definition of the reciprocal space vectors.
The reciprocal vectors are (1/a 1/(sqrt(3)*a) 0) , (0 2/(sqrt(3)*a) 0), (0 0 1/c). The lengths of
the above-mentioned segments can be computed as :
l(K-Gamma)=2/(3*a)=0.666.../a ;
l(Gamma-M)=1/(sqrt(3)*a)=0.577.../a ;
l(M-K)=1/(3*a)=0.333.../a ;
l(K-H)=1/(2*c)=0.5.../c ;
l(H-A)=2/(3*a)=0.666.../a ;
l(A-L)=1/(sqrt(3)*a)=0.577.../a ;
l(L-H)=1/(3*a)=0.333.../a ;
l(H-L)=1/(3*a)=0.333.../a ;
l(L-M)=1/(2*c)=0.5.../c ;
l(M-Gamma)=-1/(sqrt(3)*a)=0.577.../a ;
l(Gamma-A)=1/(2*c)=0.5.../c
<p>

<p>D. <b>Rhombohedral lattices</b>
<p>Rhombohedral lattices are characterised by two parameters, the length of the primitive
vectors, that we will denote a0, and the angle they form, alpha.
These can be directly input of ABINIT, as
<a href="varbas.html#acell">acell</a> and <a href="varbas.html#angdeg">angdeg</a>
<p>This will generate the primitive vectors in real space , with
<pre>
<a href="varbas.html#acell">acell</a> a0 a0 a0    and      <a href="varbas.html#rprim">rprim</a>  a 0 c    -a/2 a*sqrt(0.75) c    -a/2 -a*sqrt(0.75) c
</pre>
with a^2+c^2=1, a^2=(1-cos(alpha))*2/3, c^2=(1+2*cos(alpha))*1/3, (a/c)^2=2*(1-cos(alpha))/(1+2*cos(alpha))
and also cos(alpha)=(1-(a/c)^2/2)/(1+(a/c)^2).
Alternatively, these values of rprim might directly be the input of ABINIT (then, the balance
of the scaling factor might be adjusted between
<a href="varbas.html#acell">acell</a> and <a href="varbas.html#rprim">rprim</a>).
<p>
Unlike for the simple cubic, FCC, BCC, hexagonal (and some other) Bravais lattice,
the topology of the Brillouin zone will depend on the alpha (or a/c) value. We give below
information concerning the case when cos(alpha) is positive, that is, (a/c)^2 lower than 2.
<p>
The coordinates of several special points with respect to primitive vectors in reciprocal space will
not depend on the a/c ratio, but some others will depend on it. So, some care has to be exercised.
Notations for the Brillouin Zone special points are the same as in Phys. Rev. B 41, 11827 (1990).
<pre>
L (1/2 0 0) or (0 1/2 0) or (0 0 1/2) (or with negative signs)
T (1/2 1/2 1/2)
X (1/2 1/2 0) or (1/2 0 1/2) or (0 1/2 1/2) (or with separate negative signs)
W (5/6 - (a/c)^2/6 , 1/2 , 1/6 + (a/c)^2/6 ) = (1 0 -1)*(1-(a/c)^2/2)/3 + (1 1 1)/2
U ( (1+(a/c)^2)/6 , (8-(a/c)^2)/12 , (8-(a/c)^2)/12 ) = (-1 1/2 1/2)*(1-(a/c)^2/2)/3 + (1 1 1)/2
K (1 0 -1)*(1+(a/c)^2/4)/3
</pre>
So, if you want to specify a typical circuit, the following might do the work (the representative points on lines of symmetry are indicated - there are sometimes more than one way to go from one point to another) :
X-V-K-Sigma-Gamma-Lambda-T-Q-W-Y-L-sigma-Gamma-sigma-X . The suggestion is to sample this path
with the following coordinates
for the special points X, Gamma, T, L, Gamma, X :
<br>
<pre>kptbounds  1/2 0 -1/2   0 0 0    1/2 1/2 1/2  1 1/2 0   1 0 0  1 1/2 1/2
</pre>
<p>
In order to find the lengths of segments
(this information is useful to draw the band structure,
with the correct relative scale between special points)
one needs to know the a and c lattice parameters. Also, in what follows, we omit the 2*pi factor sometimes
present in the definition of the reciprocal space vectors.
The reciprocal vectors are (2/(3*a) 0 1/(3*c)) , -(1/(3*a) 1/(sqrt(3)*a) 1/(3*c), -(1/(3*a) -1/(sqrt(3)*a) 1/(3*c) ). The lengths of
the above-mentioned segments can be computed as :
l(X-Gamma)=2/(sqrt(3)*a)=1.155.../a , with l(K-Gamma)=(1+(a/c)^2/4)*4/(3*sqrt(3)*a);
l(Gamma-T)=1/(2*c) ;
l(T-L)=2/(sqrt(3)*a)=1.155.../a , with l(T-W)=(1-(a/c)^2/2)*4/(3*sqrt(3)*a);
l(L-Gamma)=sqrt(4/(a^2)+1/(c^2))/3
l(Gamma-X)=sqrt(1/(a^2)+1/(c^2))*2/3
<p>

"""
},
'kptgw': {
'definition': "K-PoinTs for GW calculations ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>  """,
'vartype': """real array kptgw(3,<a href="vargw.html#nkptgw">nkptgw</a>)""",
'default': "all 0.0's ",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=4, that is, sigma calculations.
<p>
For each k-point with number igwpt in the range (1:<a href="vargw.html#nkptgw">nkptgw</a>),
<b>kptgw(1,igwpt)</b> is the reduced coordinate of the k-point where GW corrections are required.
while <a href=vargw.html#bdgw">bdgw</a>(1:2,igwpt) specifies the range of bands to be considered.
<p>
At present, not all k-points are possible. Only those corresponding to the k-point
grid defined with the same repetition parameters (<a href="vargs.html#kptrlatt" target="kwimg">kptrlatt</a>,
or <a href="varbas.html#ngkpt" target="kwimg">ngkpt</a>) than the GS one, but WITHOUT any shift, are allowed.
"""
},
'kptnrm': {
'definition': "K - PoinTs NoRMalization ",
'section': "varbas",
'category': " ",
'vartype': "real parameter  ",
'default': "1.",
'text': """Establishes a normalizing denominator
for each k point.
Needed only
if <a href="varbas.html#kptopt">kptopt</a>&lt;=0, otherwise
deduced from other input variables.
<br>The k point coordinates as fractions
of reciprocal lattice translations are therefore
<a href="varbas.html#kpt">kpt</a>(mu,ikpt)/<b>kptnrm</b>.  <b>kptnrm</b> defaults to 1 and can
be ignored by the user.  It is introduced to avoid
the need for many digits in representing numbers such as 1/3.
<br>It cannot be smaller than 1.0d0
"""
},
'kptns': {
'definition': "K-PoinTs re-Normalized and Shifted",
'section': "varint",
'category': """<a href="../users/abinit_help.html#internal">INTERNAL</a> """,
'vartype': """real array <b>kptns</b>(3,<a href="vargs.html#qpt">nkpt</a>) """,
'default': "",
'text': """If <a href="vargs.html#nqpt">nqpt</a>=0, or if one is
doing a reponse calculation,
this internal variable is derived from
<a href="varbas.html#kpt">kpt</a> and <a href="varbas.html#kptnrm">kptnrm</a>:
<b>kptns</b>(1:3,:)=
<a href="varbas.html#kpt">kpt</a>(1:3,:)/
<a href="varbas.html#kptnrm">kptnrm</a>, so that
it is <a href="varbas.html#kpt">kpt</a> renormalized by
<a href="varbas.html#kptnrm">kptnrm</a>.
<br>
If <a href="vargs.html#nqpt">nqpt</a>=1 and one is
not doing a ground-state calculation,
this internal variable is derived from
<a href="varbas.html#kpt">kpt</a>,<a href="varbas.html#kptnrm">kptnrm</a>
and <a href="varint.html#qptn">qptn</a>
<b>kptns</b>(1:3,:)=
<a href="varbas.html#kpt">kpt</a>(1:3,:)/
<a href="varbas.html#kptnrm">kptnrm</a>+
<a href="varint.html#qptn">qptn</a>(1:3), so that
it is <a href="varbas.html#kpt">kpt</a> renormalized by
<a href="varbas.html#kptnrm">kptnrm</a>, then shifted
by <a href="varint.html#qptn">qptn</a>(1:3).
"""
},
'kptopt': {
'definition': "KPoinTs OPTion ",
'section': "varbas",
'category': " ",
'vartype': "integer parameter ",
'default': "1 (WARNING : was 0 prior to 5.8).",
'text': """Controls the set up of the k-points list.
The aim will be to initialize, by straight reading
or by a preprocessing approach based on other input variables,
the following input variables, giving the k points, their number,
and their weight:
<a href="varbas.html#kpt">kpt</a>,
<a href="varbas.html#kptnrm">kptnrm</a>,
<a href="varbas.html#nkpt">nkpt</a>,
and, for <a href="varbas.html#iscf">iscf</a>/=-2,
<a href="varbas.html#wtk">wtk</a>.

<p>
Often, the k points will form a lattice in reciprocal space. In this case,
one will also aim at initializing input variables that give
the reciprocal of this k-point lattice, as well as its shift with respect
to the origin:
<a href="varbas.html#ngkpt">ngkpt</a> or
<a href="vargs.html#kptrlatt">kptrlatt</a>,
as well as on <a href="varbas.html#nshiftk">nshiftk</a> and
<a href="varbas.html#shiftk">shiftk</a>.

<p>A global additional shift can be provided by <a href="varint.html#qptn">qptn</a>

<ul>
<li>0=> read directly <a href="varbas.html#nkpt">nkpt</a>, <a href="varbas.html#kpt">kpt</a>,
<a href="varbas.html#kptnrm">kptnrm</a> and <a href="varbas.html#wtk">wtk</a>.
<li>1=> rely on <a href="varbas.html#ngkpt">ngkpt</a> or
<a href="vargs.html#kptrlatt">kptrlatt</a>, as well as on
<a href="varbas.html#nshiftk">nshiftk</a> and
<a href="varbas.html#shiftk">shiftk</a> to set up the k points.
Take fully into account the symmetry to generate the
k points in the Irreducible Brillouin Zone only.<br>
(This is the usual mode for GS calculations)</li>
<li>2=> rely on
<a href="varbas.html#ngkpt">ngkpt</a> or
<a href="vargs.html#kptrlatt">kptrlatt</a>, as well as on
<a href="varbas.html#nshiftk">nshiftk</a> and
<a href="varbas.html#shiftk">shiftk</a> to set up the k points.
Take into account only the time-reversal symmetry :
k points will be generated in half the Brillouin zone.<br>
(This is to be used when preparing or executing a
RF calculation at q=(0 0 0) ) </li>
<li>3=> rely on <a href="varbas.html#ngkpt">ngkpt</a> or
<a href="vargs.html#kptrlatt">kptrlatt</a>, as well as on
<a href="varbas.html#nshiftk">nshiftk</a> and
<a href="varbas.html#shiftk">shiftk</a> to set up the k points.
Do not take into account any symmetry :
k points will be generated in the full Brillouin zone.<br>
(This is to be used when preparing or executing a
RF calculation at non-zero q )</li>
<li>4=> rely on <a href="varbas.html#ngkpt">ngkpt</a> or
<a href="vargs.html#kptrlatt">kptrlatt</a>, as well as on
<a href="varbas.html#nshiftk">nshiftk</a> and
<a href="varbas.html#shiftk">shiftk</a> to set up the k points.
Take into account all the symmetries EXCEPT the time-reversal symmetry
to generate the k points in the Irreducible Brillouin Zone.<br>
This has to be used when performing PAW calculations including
spin-orbit coupling (<a href="varpaw.html#pawspnorb">pawspnorb</a>/=0)</li>
<li>A negative value =>
rely on <a href="vargs.html#kptbounds">kptbounds</a>,
and <a href="vargs.html#ndivk">ndivk</a>
to set up a band structure calculation along different lines
(allowed only for <a href="varbas.html#iscf">iscf</a>==-2).
The absolute value of <b>kptopt</b> gives the number of segments
of the band structure.
</li>
</ul>
In the case of a grid of k points, the auxiliary variables
<a href="vargs.html#kptrlen">kptrlen</a>,
<a href="varbas.html#ngkpt">ngkpt</a>  and
<a href="varfil.html#prtkpt">prtkpt</a> might help
you to select the optimal grid.
"""
},
'kptrlatt': {
'definition': "K - PoinTs grid : Real space LATTice ",
'section': "vargs",
'category': " ",
'vartype': "integer array kptrlatt(3,3)",
'default': "No default.",
'text': """This input variable is used only when <a href="varbas.html#kptopt">kptopt</a>
is positive. It partially defines the k point grid.
The other piece of information is contained in
<a href="varbas.html#shiftk">shiftk</a>.
<b>kptrlatt</b> cannot be used together with <a href="varbas.html#ngkpt">ngkpt</a>.

<p>The values kptrlatt(1:3,1), kptrlatt(1:3,2), kptrlatt(1:3,3)
are the coordinates of three vectors in real space, expressed
in the <a href="varbas.html#rprimd">rprimd</a> coordinate system (reduced coordinates).
They defines a super-lattice in real space.
The k point lattice is the reciprocal of this super-lattice,
possibly shifted (see <a href="varbas.html#shiftk">shiftk</a>).

<p>If neither <a href="varbas.html#ngkpt">ngkpt</a> nor <b>kptrlatt</b>
are defined, ABINIT will automatically generate a set
of k point grids, and select the best combination
of <b>kptrlatt</b> and <a href="varbas.html#shiftk">shiftk</a>
that allows to reach a sufficient value of <a href="vargs.html#kptrlen">kptrlen</a>.
See this latter variable for a complete description of this
procedure.
"""
},
'kptrlen': {
'definition': "K - PoinTs grid : Real space LENgth ",
'section': "vargs",
'category': " ",
'vartype': "real parameter ",
'default': "30.0d0 (WARNING : was 20 prior to v5.8).",
'text': """This input variable is used only when <a href="varbas.html#kptopt">kptopt</a>
is positive and non-zero.

<p>Preliminary explanation :
<br>
The k point lattice defined by <a href="varbas.html#ngkpt">ngkpt</a>
or <a href="vargs.html#kptrlatt">kptrlatt</a> is used to perform integrations
of periodic quantities in the Brillouin Zone, like
the density or the kinetic energy. One can relate the
error made by replacing the continuous integral by a sum
over k point lattice to the Fourier transform of the
periodic quantity. Erroneous contributions will appear
only for the vectors in real space that belong to the reciprocal
of the k point lattice, except the origin.
Moreover, the expected size of these
contributions usually decreases exponentially with the distance.
So, the length of the smallest of these real space vectors
is a measure of the accuracy of the k point grid.

<p>When either <a href="varbas.html#ngkpt">ngkpt</a> or
<a href="vargs.html#kptrlatt">kptrlatt</a> is defined, <b>kptrlen</b> is not
used as an input variable, but the length of the
smallest vector will be placed in this variable, and echoed
in the output file.

<p>On the other hand, when neither <a href="varbas.html#ngkpt">ngkpt</a> nor
<a href="vargs.html#kptrlatt">kptrlatt</a> are defined, ABINIT will
automatically generate a large set of possible k point grids,
and select among this set, the grids that give
a length of smallest vector LARGER than <b>kptrlen</b>,
and among these grids, the one that, when used with
<a href="varbas.html#kptopt">kptopt</a>=1, reduces to the smallest number
of k points. Note that this procedure can be time-consuming.
It is worth doing it once for a given unit cell
and set of symmetries, but not use this procedure by default.
The best is then to set <a href="varfil.html#prtkpt">prtkpt</a>=1, in order
to get a detailed analysis of the set of grids.

<p>If some layer of vacuum is detected in the unit cell
(see the input variable <a href="vargs.html#vacuum">vacuum</a>), the
computation of <b>kptrlen</b> will ignore the
dimension related to the direction perpendicular
to the vacuum layer, and generate a bi-dimensional k point grid.
If the system is confined in a tube,
a one-dimensional k point grid will be generated.
For a cluster, this procedure will only generate the Gamma point.
"""
},
'kssform': {
'definition': "Kohn Sham Structure file FORMat ",
'section': "varfil",
'category': " ",
'vartype': "integer parameter  ",
'default': "1, i.e. the KSS format",
'text': """Governs the choice of the format for the file that
contains the Kohn-Sham electronic structure information,
for use in GW calculations, see the input variables
<a href="vargs.html#optdriver">optdriver</a> and
<a href="vargw.html#nbandkss">nbandkss</a>.
<ul>
<li><b>(obsolete and, starting v5.6, not supported anymore) kssform</b>=0,
the _STA file is generated together with a _VKB
file containing information on the pseudopotential </li>
<li><b>kssform</b>=1, a single file .kss (double precision) containing
complete information on the Kohn Sham Structure (eigenstates and the
pseudopotentials used) will be generated through full diagonalization
of the complete Hamiltonian matrix.
The file has at the beginning the standard abinit header</li>
<li><b>(obsolete and not supported anymore starting v5.6) kssform</b>=2, the same as 1,
but most of the relevant informations are in single precision. </li>
<li><b>kssform</b>=3, a single file .kss (double precision) containing
complete information on the Kohn Sham Structure (eigenstates and the
pseudopotentials used) will be generated through the usual conjugate gradient
algorithm (so, a restricted number of states)
The file has at the beginning the standard abinit header</li>
</ul>
<p>Very important : for the time being, <a href="vardev.html#istwfk">istwfk</a>
must be 1 for all the k-points.</p>
"""
},
'ldgapp': {
'definition': "Lein-Dobson-Gross approximation",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>  """,
'vartype': "integer parameter  ",
'default': "0.",
'text': """<p>Concern only the ACFD computation of the correlation energy
(<a href="vargs.html#optdriver">optdriver</a>=3).
<br>If <b>ldgapp</b> &gt; 0, the Lein, Dobson and Gross first-order
approximation to the correlation energy is also computed during the ACFD run.
[See Lein, Dobson and Gross, J. Comput. Chem. 20,12 (1999)]. This
is only implemented for the RPA, for the PGG kernel and
for the linear energy optimized kernel at the present time.
"""
},
'lexexch': {
'definition': "value of angular momentum L for EXact EXCHange ",
'section': "varpaw",
'category': "",
'vartype': """integer array lexexch(<a href="varbas.html#ntypat">ntypat</a>)""",
'default': "-1",
'text': """Activated if
<a href="varpaw.html#useexexch">useexexch</a> is equal to 1.<br>
Give for each species the value of the angular momentum (only values 2 or 3 are allowed)
on which to apply the exact exchange correction.
"""
},
'localrdwf': {
'definition': "LOCAL ReaD WaveFunctions",
'section': "varpar",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>, PARALLEL """,
'vartype': "integer ",
'default': "1.",
'text': """This input variable is used only when running abinit in parallel.
If <b>localrdwf</b>=1, the input wavefunction disk file or the KSS/SCR file in case of GW
calculations, is read locally by each processor, while
if <b>localrdwf</b>=0, only one processor reads it, and
broadcast the data to the other processors.
<p> The option <b>localrdwf</b>=0 is NOT allowed when parallel I/O are activated (MPI-IO access),
i.e. when <a href="vardev.html#accesswff">accesswff</a>==1.
<p> The option <b>localrdwf</b>=0 is NOT allowed when
<a href="varfil.html#mkmem">mkmem</a>==0 (or, for RF, when
<a href="varrf.html#mkqmem">mkqmem</a>==0, or <a href="varrf.html#mk1mem">mk1mem</a>==0), that
is, when the wavefunctions are stored on disk.
This is still to be coded ...
<p> In the case of a parallel computer with a unique file system,
both options are as convenient for the user. However, if the I/O
are slow compared to communications between processors,
(e.g. for CRAY T3E machines), <b>localrdwf</b>=0 should be much more
efficient;
if you really need temporary disk storage, switch to localrdwf=1 ).
<p> In the case of a cluster of nodes, with a different file system for
each machine, the input wavefunction file must be available on all
nodes if <b>localrdwf</b>=1, while it is needed only for the
master node if <b>localrdwf</b>=0.

"""
},
'lpawu': {
'definition': "value of angular momentum L for PAW+U",
'section': "varpaw",
'category': "",
'vartype': """integer array lpawu(<a href="varbas.html#ntypat">ntypat</a>)""",
'default': "-1 ",
'text': """Activated if <a href="varpaw.html#usepawu">usepawu</a>
is equal to 1 or 2.<br>
Give for each species the value of the angular momentum (only
values 2 or 3 are allowed)&nbsp; on which to apply the LDA+U correction.<br>
<ul>
<li>If equal to 2 (d-orbitals)&nbsp; or 3 (f-orbitals), values of
<a href="varpaw.html#upawu">upawu</a> and&nbsp;
<a href="varpaw.html#jpawu">jpawu</a> are used in the calculation.</li>
<li>If equal to -1: do not apply LDA+U correction on the species.</li>
</ul>
"""
},
'macro_uj': {
'definition': "Macro variable that activates the determination of the U and J parameter (for the PAW+U calculations) ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>  """,
'vartype': "integer parameter  ",
'default': "0.",
'text': """<p>Sets proper input values for the determination of U and J i.e.
for <a href="vardev.html#pawujat">pawujat</a> (first atom treated with PAW+U),
<a href="varfil.html#irdwfk">irdwfk</a> (=1),
<a href="varbas.html#tolvrs">tolvrs</a> (=10^(-8)),
<a href="varbas.html#nstep">nstep</a> (=255),
<a href="vargs.html#diemix">diemix</a> (=0.45),
<a href="varff.html#atvshift">atvshift</a> (<a href="vardev.html#pawujat">pawujat</a>) <a href="vardev.html#pawujv">pawujv</a>). Do not overwrite these variables manually unless you know what you do.

<ul>
<li> <b>macro_uj</b>=1 (and <a href="varbas.html#nsppol">nsppol</a>=2) Standard procedure to determine U on atom pawujat through a shift of the potential on both spin channels.</li>
<li> <b>macro_uj</b>=1 (and <a href="varbas.html#nsppol">nsppol</a>=1) Non standand procedure to determine U from potential shift on atom pawujat (experimental).</li>
<li> <b>macro_uj</b>=2 (and <a href="varbas.html#nsppol">nsppol</a>=2) Non standand procedure to determine U from potential shift on atom pawujat through a shift on spin channel 1 on this atom and the response on this channel (experimental).</li>
<li> <b>macro_uj</b>=3 (and <a href="varbas.html#nsppol">nsppol</a>=2) Standand procedure to determine J from potential shift on spin channel 1 on atom pawujat and response on spin channel 2 (experimental).</li>
</ul>

Determination of U and J can be done only if the symmetry of the atomic arrangement is reduced and the atom pawujat is not connected to any other atom by symmetry relations (either input reduced symmetries manually, define concerned atom as a separate atomic species or shift concerned atom from ideal postion).

"""
},
'maxestep': {
'definition': "MAXimum Electric field STEP",
'section': "varff",
'category': "",
'vartype': "real parameter",
'default': "0.005",
'text': """This variable controls the maximum change of electric field when updating the electric field after each SCF iteration.
When the calculation is difficult to converge, try reducing this value or reducing <a href="varff.html#ddamp">ddamp</a>.
This variable is used in finite electric displacement field calculations (<a href="varff.html#berryopt">berryopt</a>=6,16,17).
"""
},
'maxnsym': {
'definition': "MAXimum Number of SYMetries",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>  """,
'vartype': "integer parameter  ",
'default': "384.",
'text': """<p>Gives the maximum number of spatial symetries allowed in the memory.<br>
The default value is sufficient for most applications; it has to be increase in the case of the use of a supercell (unit cell identically repeated).
"""
},
'mband': {
'definition': "Maximum number of BANDs ",
'section': "varint",
'category': """<a href="../users/abinit_help.html#internal">INTERNAL</a> """,
'vartype': "integer ",
'default': "",
'text': """This internal variable derives
the maximum number of bands
over all k-points and spin-polarisation from
<a href="varbas.html#nband">nband</a>(1:nkpt*nsppol).
"""
},
'mdtemp': {
'definition': "Molecular Dynamics Temperatures ",
'section': "varrlx",
'category': "",
'vartype': "real array mdtemp(2)",
'default': "<b>mdtemp</b>= (300, 300)",
'text': """Give the initial and final temperature
of the Nose-Hoover thermostat (<a href="varrlx.html#ionmov">ionmov</a>=8)
and Langevin dynamics (<a href="varrlx.html#ionmov">ionmov</a>=9), in
Kelvin.
This temperature will change linearly from the initial temperature
<b>mdtemp(1)</b> at itime=1 to
the final temperature <b>mdtemp(2)</b> at the end of the
<a href="varrlx.html#ntime">ntime</a> timesteps.
"""
},
'mdwall': {
'definition': "Molecular Dynamics WALL location ",
'section': "varrlx",
'category': "",
'vartype': "real parameter ",
'default': "10000.0 Bohr (the walls are extremely far away).",
'text': """Gives the location (atomic units) of walls
on which the atoms will bounce back.
when <a href="varrlx.html#ionmov">ionmov</a>=6, 7, 8 or 9. For each
cartesian direction idir=1, 2 or 3, there is a pair of walls with
coordinates xcart(idir)=-wall and xcart(idir)=rprimd(idir,idir)+wall .
Supposing the particle will cross the wall, its velocity normal to the
wall is reversed, so that it bounces back.
<br>
By default, given in Bohr atomic units
(1 Bohr=0.5291772108 Angstroms), although Angstrom can be specified,
if preferred, since <b>mdwall</b> has the
'<a href="../users/abinit_help.html#dimensions">LENGTH</a>'
characteristics.
"""
},
'mep_mxstep': {
'definition': "Minimal Energy Path search: MaXimum allowed STEP size",
'section': "varrlx",
'category': """<a href="../users/abinit_help.html#length">LENGTH</a>""",
'vartype': "real parameter",
'default': """<b>mep_mxstep</b>=0.4 Bohr if <a href="varrlx.html#imgmov">imgmov</a>=5 (NEB) otherwise <b>mep_mxstep</b>=100. Bohr""",
'text': """Relevant only when <a href="varrlx.html#imgmov">imgmov</a>=1 (Steepest-Descent), 2 (String Method) or 5 (Nudged Elastic Band).<br>
The optimizer used to solve the Ordinary Differential Equation (ODE) can be constrained with a maximum allowed step size for each image. By default this feature is only activated for Nudged Elastic Band (NEB) and the value is inspired by <i>J. Chem. Phys. 128, 134106 (2008)</i>.<br>
Note that the step size is defined for each image as
<i>step = SQRT[SUM(R_i dot R_i)]</i> where the <i>R_i</i> are the positions of
the atoms in the cell.
"""
},
'mep_solver': {
'definition': "Minimal Energy Path ordinary differential equation SOLVER ",
'section': "varrlx",
'category': "",
'vartype': "integer parameter",
'default': "<b>mep_solver</b>=0",
'text': """Relevant only when <a href="varrlx.html#imgmov">imgmov</a>=2 (String Method) or 5 (Nudged Elastic Band).<br>
Gives the algorithm used to solve the Ordinary Differential Equation (ODE) when searching for
a Minimal Energy Path (MEP).<br>
Possible values can be:<br>
<ul>
<li>0=&gt; <b>Steepest-Descent algorithm</b> following the (scaled) forces,
the scaling factor being <a href="varrlx.html#fxcartfactor">fxcartfactor</a>
(forward Euler method).<br>
Compatible with all MEP search methods.
</li>
<br>
<li>1=&gt; <b>Quick-min optimizer</b> following the (scaled) forces,
the scaling factor being <a href="varrlx.html#fxcartfactor">fxcartfactor</a>.
The "quick minimizer" improves upon the steepest-descent method by
accelerating the system in the direction of the forces. The velocity (of the image) is
projected long the force and cancelled if antiparallel to it.<br>
Compatible only with Nudged Elastic Band (<a href="varrlx.html#imgmov">imgmov</a>=5).<br>
<i>See, for instance: J. Chem. Phys. 128, 134106 (2008).</i><br>
</li>
<br>
<li>2=&gt; <b>Local Broyden-Fletcher-Goldfarb-Shanno (L-BFGS) algorithm</b>; each image along the
band is minimized with a different instance of the BFGS optimizer.<br>
Compatible only with Nudged Elastic Band (<a href="varrlx.html#imgmov">imgmov</a>=5).<br>
<i>See, for instance: J. Chem. Phys. 128, 134106 (2008).</i><br>
IN DEVELOPPMENT - NOT RELIABLE
</li>
<br>
<li>3=&gt; <b>Global Broyden-Fletcher-Goldfarb-Shanno (GL-BFGS) algorithm</b>; all images along the
band are minimized with a single instance of the BFGS optimizer.<br>
Compatible only with Nudged Elastic Band (<a href="varrlx.html#imgmov">imgmov</a>=5).<br>
<i>See, for instance: J. Chem. Phys. 128, 134106 (2008).</i><br>
IN DEVELOPPMENT - NOT RELIABLE
</li>
<br>
<li>4=&gt; <b>Fourth-order Runge-Kutta method</b>; the images along the band are moved
every four steps (1&lt=istep&lt=<a href="varrlx.html#ntimimage">ntimimage</a>)
following the Runge-Kutta algorithm,
the time step being <a href="varrlx.html#fxcartfactor">fxcartfactor</a>.<br>
Compatible only with Simplified String Method (<a href="varrlx.html#imgmov">imgmov</a>=2
and <a href="varrlx.html#string_algo">string_algo</a>=1 or 2).<br>
<i>See: J. Chem. Phys. 126, 164103 (2007).</i>
</li>
</ul>
All of the optimizers can be constrained with a maximum allowed step size for each image; see <a href="varrlx.html#mep_mxstep">mep_mxstep</a>. This is by default the case of the Nudged Elastic Band (<a href="varrlx.html#imgmov">imgmov</a>=5).
"""
},
'mffmem': {
'definition': "Maximum number of FFt grids in MEMory ",
'section': "varfil",
'category': " ",
'vartype': "integer parameter  ",
'default': "1, i.e. an in-core solution.",
'text': """Governs the choice of number of
FFT arrays that will be kept permanently in core memory.
<br>The allowed values are 0, in which case maximal use
is made of disk space, saving core memory at the expense of
execution time (not much, usually),
or 1, in which case everything is kept in core memory.
<br>More detailed explanations : if <b>mffmem</b>==0, some arrays
of size <br>double precision :: xx(nfft,<a href="varbas.html#nsppol">nsppol</a>) <br>will be saved
on disk when the wavefunctions are optimized or when
the Hartree and xc potential is computed (which can require
some sizeable memory space also).
<br>The number of these arrays is
<br> 5 if <a href="varbas.html#iscf">iscf</a>==1,
<br> 3 if <a href="varbas.html#iscf">iscf</a>==2 or 12,
<br> 4 if <a href="varbas.html#iscf">iscf</a>==3 or 3,
<br> 6 if <a href="varbas.html#iscf">iscf</a>==4 or 14,
<br> 10 if <a href="varbas.html#iscf">iscf</a>==5,
<br> (2+2*<a href="vardev.html#npulayit">npulayit</a>) if <a href="varbas.html#iscf">iscf</a>==7 or 17.
<br>The saving of memory can be appreciable especially when <a href="varbas.html#iscf">iscf</a>==5, 7 ,17 and <a href="varbas.html#nsppol">
nsppol</a>=2."""
},
'mgfft': {
'definition': "Maximum of nGFFT ",
'section': "varint",
'category': """<a href="../users/abinit_help.html#internal">INTERNAL</a> """,
'vartype': "integer ",
'default': "",
'text': """This internal variable contains the maximum of
<a href="vargs.html#ngfft">ngfft</a>(1:3).
"""
},
'mgfftdg': {
'definition': "Maximum of nGFFT for the Double Grid",
'section': "varint",
'category': """<a href="../users/abinit_help.html#internal">INTERNAL</a> """,
'vartype': "integer ",
'default': "",
'text': """This internal variable contains the maximum of
<a href="varpaw.html#ngfftdg">ngfftdg</a>(1:3).
"""
},
'mixalch': {
'definition': "MIXing coefficients for ALCHemical potentials",
'section': "vargs",
'category': "EVOLVING ",
'vartype': """integer array mixalch(<a href="vargs.html#npspalch">npspalch</a>,<a href="vargs.html#ntypalch">ntypalch</a>)""",
'default': "0.d0 (will not accepted !)",
'text': """<p> Used for the generation of alchemical pseudoatoms, that is,
when <a href="vargs.html#ntypalch">ntypalch</a> is non-zero.
<p> This array gives, for each type of alchemical pseudatom (there are
<a href="vargs.html#ntypalch">ntypalch</a> such pseudoatoms), the mixing coefficients
of the basic <a href="vargs.html#npspalch">npspalch</a> pseudopotentials for
alchemical use. For each type of alchemical pseudoatom, the sum of the
mixing coefficients must equal 1.
<p> The actual use of the mixing coefficients is defined by the input
variable <a href="vargs.html#algalch">algalch</a>. Note that the masses of the atoms,
<a href="varrlx.html#amu">amu</a>
are also mixed according to the value of <b>mixalch</b>, by default.
<p> Example 1. Suppose that we want to describe Ba(0.25) Sr(0.75) Ti O3.
<br> The input variables related to the construction of the alchemical Ba(0.25) Sr(0.75)
potential will be :
<pre>
npsp   4                 ! 4 pseudopotentials should be read.
znucl  8 40 56 38        ! The nuclear charges. Note that the two
! atoms whose pseudopotentials are to be mixed
! are mentioned at the end of the series.
ntypat  3                ! There will be three types of atoms.
ntypalch   1             ! One pseudoatom will be alchemical.
! Hence, there will be ntyppure=2 pure pseudoatoms,
! with znucl 8 (O) and 40 (Ti), corresponding to
! the two first pseudopotentials. Out of the
! four pseudopotentials, npspalch=2 are left
! for alchemical purposes, with znucl 56 (Ba)
! and 38 (Sr).
mixalch    0.25  0.75    ! For that unique pseudoatom to be
! generated, here are the mixing coeeficients,
! to be used to combine the Ba and Sr pseudopotentials.
</pre>
<p> Example 2. More complicated, and illustrate some minor drawback of the
design of input variables.
Suppose that one wants to generate Al(0.25)Ga(0.75) As(0.10)Sb(0.90).
<br> The input variables will be :
<pre>
npsp  4                  ! 4 pseudopotentials should be read
znucl  13 31 33 51       ! The atomic numbers. All pseudopotentials
! will be used for some alchemical purpose
ntypat  2                ! There will be two types of atoms.
ntypalch   2             ! None of the atoms will be "pure".
! Hence, there will be npspalch=4 pseudopotentials
!  to be used for alchemical purposes.
mixalch    0.25  0.75 0.0  0.0   ! This array is a (4,2) array, arranged in the
0.0   0.0  0.1  0.9   ! usual Fortran order.
</pre>
Minor drawback : one should not forget to fill <b>mixalch</b> with the needed zero's, in this later case.
<p> In most cases, the use of <b>mixalch</b>
will be as a static (non-evolving) variable. However, the possibility to have
different values of <b>mixalch</b> for different images has been coded. A population of
cells with different atomic characteristics can thus be considered,
and can be made to evolve, e.g. with a genetic algorithm (not coded in v7.0.0 though).
There is one restriction to this possibility : the value of <a href="varint.html#ziontypat">ziontypat</a> for the atoms that are mixed should be
identical.
"""
},
'mk1mem': {
'definition': "Maximum number of K - points for 1st order wavefunctions, kept in MEMory ",
'section': "varrf",
'category': """<a href="../users/abinit_help.html#respfn">RESPFN</a>  """,
'vartype': "integer parameters  ",
'default': """<a href="varbas.html#nkpt">nkpt</a>, i.e. in-core solution.""",
'text': """Plays a role similar to <a href="varfil.html#mkmem">mkmem</a>
but for different sets of wavefunctions : the
ground state wavefunctions at k+q and the first-order
wavefunctions. Only needed for response calculations.
<br>Internal representation as mkmems(2) and mkmems(3).
<br>Note (991019) that although the effective number of k points
can be reduced thanks to symmetry for different
perturbations, <b>mkqmem</b> and <b>mk1mem</b> are presently
still compared with the input <a href="varbas.html#nkpt">nkpt</a>. This should be changed
later."""
},
'mkmem': {
'definition': "Maximum number of K - points in MEMory ",
'section': "varfil",
'category': " ",
'vartype': "integer parameter ",
'default': """<a href="varbas.html#nkpt">nkpt</a>, i.e. in-core solution.""",
'text': """Sets the maximum number of k
points for which the ground state wavefunctions
are kept in core memory at one time.
<br>This value should either be 0, in which case an out-of-core
solution will be used, or else <a href="varbas.html#nkpt">nkpt</a>,
in which case an in-core solution will be used.
<br>
This input variable can be used during a GW calculation to reduce the memory allocated
by each processor (presently only for <a href="vargs.html#optdriver">optdriver</a>=3).
See also the related input variable <a href="varpar.html#gwpara">gwpara</a>.

<br>Internal representation as mkmems(1)"""
},
'mkqmem': {
'definition': "Maximum number of K+Q - points in MEMory ",
'section': "varrf",
'category': """<a href="../users/abinit_help.html#respfn">RESPFN</a>  """,
'vartype': "integer parameters  ",
'default': """<a href="varbas.html#nkpt">nkpt</a>, i.e. in-core solution.""",
'text': """Plays a role similar to <a href="varfil.html#mkmem">mkmem</a>
but for different sets of wavefunctions : the
ground state wavefunctions at k+q and the first-order
wavefunctions. Only needed for response calculations.
<br>Internal representation as mkmems(2) and mkmems(3).
<br>Note (991019) that although the effective number of k points
can be reduced thanks to symmetry for different
perturbations, <b>mkqmem</b> and <b>mk1mem</b> are presently
still compared with the input <a href="varbas.html#nkpt">nkpt</a>. This should be changed
later."""
},
'mpw': {
'definition': "Maximum number of Plane Waves ",
'section': "varint",
'category': """<a href="../users/abinit_help.html#internal">INTERNAL</a> """,
'vartype': "integer ",
'default': "",
'text': """This internal variable gives the maximum of the number of
plane waves over all k-points. It is computed
from <a href="varbas.html#ecut">ecut</a> and the description
of the cell, provided by
<a href="varbas.html#acell">acell</a>,
<a href="varbas.html#rprim">rprim</a>, and/or
<a href="varbas.html#angdeg">angdeg</a>.
"""
},
'mqgrid': {
'definition': "Maximum number of Q-space GRID points for pseudopotentials",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>  """,
'vartype': "integer parameter  ",
'default': "3001.",
'text': """<p>Govern the size of the one-dimensional information
related to pseudopotentials, in reciprocal space :
potentials, or projector functions.
"""
},
'mqgriddg': {
'definition': "",
'section': "varpaw",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a> """,
'vartype': "integer parameter ",
'default': "3001 ",
'text': """Maximum number of wavevectors used to sample the local part of the potential, in PAW.
Actually referred to as mqgrid_vl internally. Should change name to the latter ...
See also <a href="vardev.html#mqgrid">mqgrid</a>
"""
},
'natcon': {
'definition': "Number of AToms in CONstraint equations",
'section': "varrlx",
'category': """<a href="../users/abinit_help.html#no_multi">NO MULTI</a>, <a href="../users/abinit_help.html#not_internal">NOT INTERNAL</a> """,
'vartype': """integer array of length <a href="varrlx.html#nconeq">nconeq</a>""",
'default': "0 ",
'text': """Gives the number of atoms appearing in each of the <a
href="varrlx.html#nconeq">nconeq</a>
independent equations constraining the motion of
atoms during structural optimization or molecular dynamics (see <a
href="varrlx.html#nconeq">nconeq</a>, <a href="varrlx.html#iatcon">iatcon</a>,
and <a href="varrlx.html#wtatcon">wtatcon</a>).
"""
},
'natfix': {
'definition': "Number of Atoms that are FIXed ",
'section': "varrlx",
'category': "NOT INTERNAL ",
'vartype': "integer parameter ",
'default': "0 (no atoms held fixed).",
'text': """Gives the number of atoms (not to exceed
<a href="varbas.html#natom">natom</a>) which are to be held fixed
during a structural
optimization or molecular dynamics. <br>
When <b>natfix</b> &gt; 0,
<b>natfix</b> entries should be provided in array <a
href="varrlx.html#iatfix">iatfix</a>.
<br>
When <b>natfixx</b> &gt; 0, <b>natfixx</b> entries should be provided
in array <a href="varrlx.html#iatfixx">iatfixx</a>, and so on ...
"""
},
'natfixx': {
'definition': "Number of Atoms that are FIXed along the X direction ",
'section': "varrlx",
'category': "NOT INTERNAL ",
'vartype': "integer parameter ",
'default': "0 (no atoms held fixed).",
'text': """Gives the number of atoms (not to exceed
<a href="varbas.html#natom">natom</a>) which are to be held fixed
during a structural
optimization or molecular dynamics. <br>
When <b>natfix</b> &gt; 0,
<b>natfix</b> entries should be provided in array <a
href="varrlx.html#iatfix">iatfix</a>.
<br>
When <b>natfixx</b> &gt; 0, <b>natfixx</b> entries should be provided
in array <a href="varrlx.html#iatfixx">iatfixx</a>, and so on ...
"""
},
'natfixy': {
'definition': "Number of Atoms that are FIXed along the Y direction ",
'section': "varrlx",
'category': "NOT INTERNAL ",
'vartype': "integer parameter ",
'default': "0 (no atoms held fixed).",
'text': """Gives the number of atoms (not to exceed
<a href="varbas.html#natom">natom</a>) which are to be held fixed
during a structural
optimization or molecular dynamics. <br>
When <b>natfix</b> &gt; 0,
<b>natfix</b> entries should be provided in array <a
href="varrlx.html#iatfix">iatfix</a>.
<br>
When <b>natfixx</b> &gt; 0, <b>natfixx</b> entries should be provided
in array <a href="varrlx.html#iatfixx">iatfixx</a>, and so on ...
"""
},
'natfixz': {
'definition': "Number of Atoms that are FIXed along the Z direction ",
'section': "varrlx",
'category': "NOT INTERNAL ",
'vartype': "integer parameter ",
'default': "0 (no atoms held fixed).",
'text': """Gives the number of atoms (not to exceed
<a href="varbas.html#natom">natom</a>) which are to be held fixed
during a structural
optimization or molecular dynamics. <br>
When <b>natfix</b> &gt; 0,
<b>natfix</b> entries should be provided in array <a
href="varrlx.html#iatfix">iatfix</a>.
<br>
When <b>natfixx</b> &gt; 0, <b>natfixx</b> entries should be provided
in array <a href="varrlx.html#iatfixx">iatfixx</a>, and so on ...
"""
},
'natom': {
'definition': "Number of ATOMs ",
'section': "varbas",
'category': " ",
'vartype': "integer parameter  ",
'default': "1 ",
'text': """Gives the total number of atoms in the unit cell.
Default is 1 but you will obviously want to input this
value explicitly.
<br>Note that <b>natom</b> refers to all atoms in the unit cell, not
only to the irreducible set of atoms in the unit cell (using symmetry operations,
this set allows to recover all atoms). If you want
to specify only the irreducible set of atoms, use the
symmetriser, see the input variable <a href="vargeo.html#natrd">natrd</a>."""
},
'natpawu': {
'definition': "Number of AToms on which PAW+U is applied",
'section': "varint",
'category': """<a href="../users/abinit_help.html#internal">INTERNAL</a> """,
'vartype': "integer parameter  ",
'default': "",
'text': """This internal variable gives the number of atoms on which LDA/GGA+U method
is applied.
<br>It is determined by <a href="varbas.html#natom">natom</a>,
<a href="varpaw.html#usepawu">usepawu</a> and
<a href="varpaw.html#lpawu">lpawu</a> input keywords.
"""
},
'natrd': {
'definition': "Number of AToms ReaD ",
'section': "vargeo",
'category': """<a href="../users/abinit_help.html#geometry_builder">GEOMETRY BUILDER</a>,  SYMMETRISER  """,
'vartype': "integer parameter  ",
'default': """<a href="varbas.html#natom">natom</a>.""",
'text': """Gives the number of atoms to be
read from the input file, in the case the geometry builder
or the symmetriser is used. In this case,
<b>natrd</b> is also used to dimension
the array <a href="varbas.html#typat">typat</a>,
and the arrays <a href="varbas.html#xred">xred</a>,
<a href="varbas.html#xangst">xangst</a> and <a href="varbas.html#xcart">xcart</a>.
<br>Must take into account the vacancies (see <a href="vargeo.html#vacnum">vacnum</a> and
<a href="vargeo.html#vaclst">vaclst</a>).
<br>Despite possible vacancies, cannot be bigger than <a href="varbas.html#natom">natom</a>.
"""
},
'natsph': {
'definition': "Number of ATomic SPHeres for the atom-projected density-of-states",
'section': "vargs",
'category': " ",
'vartype': "integer parameter  ",
'default': """<a href="varbas.html#natom">natom</a> """,
'text': """This input variable is active only in the
<a href="varfil.html#prtdos">prtdos</a>=3 case or if
<a href="varpaw.html#pawfatbnd">pawfatbnd</a>=1 or 2.
<br>It gives the number of atoms around which the sphere
for atom-projected density-of-states will be built,
in the <a href="varfil.html#prtdos">prtdos</a>=3 case.
The indices of these atoms are given by <a href="vargs.html#iatsph">iatsph</a>.
The radius of these spheres is given by <a href="vargs.html#ratsph">ratsph</a>.
<br> If <a href="varpaw.html#pawfatbnd">pawfatbnd</a>=1 or 2, it gives the number of atoms around which atom-projected band structure will be
built (the indices of these atoms are given by  <a href="vargs.html#iatsph">iatsph</a>).
"""
},
'natvshift': {
'definition': "Number of ATomic potential (V) energy SHIFTs (per atom) ",
'section': "varff",
'category': "",
'vartype': "integer parameter ",
'default': "0.",
'text': """Number of atomic potential energy shifts (per atom), to be used to define the
array <a href="varff.html#atvshift">atvshift</a>.
If non-zero, only two possibilities exist : 5 for d states
(with <a href="varpaw.html#lpawu">lpawu</a>=2),
and 7 for f states (with <a href="varpaw.html#lpawu">lpawu</a>=3).
If non-zero, one should define
<a href="varpaw.html#usepawu">usepawu</a>,
<a href="varpaw.html#lpawu">lpawu</a> and
<a href="varff.html#atvshift">atvshift</a>.
"""
},
'nband': {
'definition': "Number of BANDs ",
'section': "varbas",
'category': " ",
'vartype': "integer parameter  ",
'default': "1.",
'text': """Gives number of bands, occupied plus
possibly unoccupied, for which wavefunctions are being computed
along with eigenvalues.
<br>Note : if the parameter
<a href="varbas.html#occopt">occopt</a> (see below) is not set to 2,
<b>nband</b> is a scalar integer, but
if the parameter <a href="varbas.html#occopt">occopt</a> is set to 2,
then <b>nband</b> must be an array <b>nband</b>(<a href="varbas.html#nkpt">nkpt</a>*
<a href="varbas.html#nsppol">nsppol</a>) giving the
number of bands explicitly for each k point.  This
option is provided in order to allow the number of
bands treated to vary from k point to k point.
<br>For the values of <a href="varbas.html#occopt">occopt</a> not equal to 0 or 2, <b>nband</b>
can be omitted. The number of bands will be set up
thanks to the use of the variable <a href="vargs.html#fband">fband</a>. The present Default
will not be used.
<p>If <a href="vargs.html#nspinor">nspinor</a> is 2, nband must be even for
each k point.
<p>In the case of a GW calculation (<a href="vargs.html#optdriver">optdriver</a>=3 or 4),
<b>nband</b> gives the number of bands to be treated to generate the screening (susceptibility
and dielectric matrix), as well as the self-energy. However, to generate the _KSS
file (see <a href="varfil.html#kssform">kssform</a>)
the relevant number of bands is given by <a href="vargw.html#nbandkss">nbandkss</a>.

"""
},
'nbandkss': {
'definition': "Number of BaNDs in the KSS file ",
'section': "vargw",
'category': " ",
'vartype': "integer parameter  ",
'default': "0",
'text': """<p>
This input variable is used for the preparation of a GW calculation :
it is used in a GS run (where <a href="vargs.html#optdriver">optdriver</a>=0) to generate a _KSS file.
In this run, <b>nbandkss</b> should be non-zero.
The generated _KSS file can be subsequently used to calculate the irreducible polarizabilty
$\chi^{(0)}_{KS}$ using <a href="vargs.html#optdriver">optdriver</a>=3
or to calculate GW corrections setting <a href="vargs.html#optdriver">optdriver</a>=4.
<p>

<ul>
<li>If <b>nbandkss</b>=0, no _KSS file is created</li>
<li>If <b>nbandkss</b>=-1, all the available eigenstates (energies and eigenfunctions) are stored in the
abo_KSS file at the end of the ground state calculation. The number of states is forced to be
the same for all k-points : it will be the minimum of the number of plane waves over all k-points.</li>
<li>If <b>nbandkss</b> is greater than 0, abinit stores (about) <b>nbandkss</b> eigenstates in the abo_KSS file.
This number of states is forced to be the same for all k-points.</li>
</ul>

<p>
See <a href="vargw.html#npwkss">npwkss</a> for the selection of the number of the planewave components of
the eigenstates to be stored.
<br>
The input variable <a href="vardev.html#accesswff">accesswff</a> can be used
to read and write KSS files according to different fileformat
(presently only <a href="vardev.html#accesswff">accesswff</a>=0 and 3 are available in the GW part).
<br>
The precision of the KSS file can be tuned through the input variable <a href="varfil.html#kssform">kssform</a>.
<br>
For more details about the format of the abo_KSS file, see the routine outkss.F90.
<p>
Very important : for the time being, <a href="vardev.html#istwfk">istwfk</a> must be 1 for all the k-points
in order to generate a _KSS file.
"""
},
'nbandsus': {
'definition': "Number of BANDs to compute the SUSceptibility ",
'section': "vardev",
'category': "",
'vartype': "integer parameter ",
'default': """<a href="varbas.html#nband">nband</a>.""",
'text': """Number of bands to be used in the calculation of the susceptibility matrix (ACFD only).
"""
},
'nbdblock': {
'definition': "Number of BanDs in a BLOCK ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>  """,
'vartype': "integer parameter  ",
'default': "1",
'text': """<p>In case of non-standard, blocked algorithms for the
optimization of the wavefunctions (that is, if
<a href="vardev.html#wfoptalg">wfoptalg</a>=4):
<ul>
<li>if <a href="vardev.html#wfoptalg">wfoptalg</a>=4,
<b>nbdblock</b> defines the number of blocks (the number of bands in the block is
then <a href="varbas.html#nband">nband</a>/nbdblock ).</li>
</ul>
"""
},
'nbdbuf': {
'definition': "Number of BanDs for the BUFfer",
'section': "vargs",
'category': "",
'vartype': "integer parameter ",
'default': "0. However, the default is changed to 2 in some cases, see later.",
'text': """<b>nbdbuf</b> gives the number of bands, the highest in energy, that,
among the
<a href="varbas.html#nband">nband</a> bands, are to be considered
as part of a buffer. This concept is useful in three situations:
in non-self-consistent
calculations, for the determination of the convergence tolerance ;
for response functions of metals, to avoid instabilities,
and also when finite electric fields or non-linear responses (with electric field
perturbations) are considered.
For the two first, the need of a buffer is a natural requirement
of the problem, so that the default value is changed to 2 automatically,
as explained in the following.
The third case is only for implementation convenience.

<p> In non-self-consistent GS calculations (<a href="varbas.html#iscf">iscf</a>&lt;0),
the highest levels might be
difficult to converge, if they are degenerate with another level,
that does not belong to the set of bands treated. Then, it might
take extremely long to reach <a href="varbas.html#tolwfr">tolwfr</a>, although
the other bands are already extremely well-converged, and the energy
of the highest bands (whose residual are not yet good enough), is
also rather well converged.
<br> In response to this problem, for non-zero <b>nbdbuf</b>, the
largest residual (residm), to be later compared with <a href="varbas.html#tolwfr">tolwfr</a>,
will be computed only in the set of non-buffer bands (this modification
applies for non-self-consistent as well as self-consistent calculation,
for GS as well as RF calculations).
<br> For a GS calculation, with <a href="varbas.html#iscf">iscf</a>&lt;0, supposing
<b>nbdbuf</b> is not initialized in the input file,
then ABINIT will overcome the default <b>nbdbuf</b> value,
and automatically set <b>nbdbuf</b> to 2.
<p> In metallic RF calculations, in the conjugate gradient optimisation
of first-order wavefunctions, there is an instability situation
when the q wavevector of the perturbation brings the eigenenergy of the
highest treated band at some k point higher than the lowest
untreated eigenenergy at some k+q point.
If one accepts a buffer of frozen states, this instability can be made to
disappear. Frozen states receive automatically a residual value of -0.1d0.
<br> For a RF calculation, with 3&lt;=<a href="varbas.html#occopt">occopt</a>&lt;=7,
supposing
<b>nbdbuf</b> is not initialized in the input file, then
ABINIT will overcome the default <b>nbdbuf</b> value,
and automatically set <b>nbdbuf</b> to 2. This value might be too low
in some cases.

<p>Also, the number of active bands, in all cases, is imposed
to be at least 1, irrespective of the value of <b>nbdbuf</b>.
"""
},
'nberry': {
'definition': "Number of BERRY phase computations",
'section': "varff",
'category': "",
'vartype': "integer nberry",
'default': "1",
'text': """<p> Used for non-zero values of <a href="varff.html#berryopt">berryopt</a>.
<p> Gives the number of Berry phase computations of polarisation,
or finite-difference estimations of the derivative of wavefunctions
with respect to the wavevector,
each of which might be characterized by a different change of
wavevector <a href="varff.html#kberry">kberry</a>.
<p> When equal to 0, no Berry phase calculation of polarisation
is performed. The maximal value of <b>nberry</b> is 20.
<p> Note that the computation of the polarisation for a set of bands
having different occupation numbers is meaningless (although
in the case of spin-polarized calculations, the spin up bands
might have an identical occupation number, that might differ
from the identical occupation number of spin down bands).
Although meaningless, ABINIT will perform such computation,
if required by the user. The input variable
<a href="varff.html#bdberry">bdberry</a> governs the set of bands
for which a Berry phase is computed.
<p> The computation of the Berry phase is not yet implemented
for spinor wavefunctions (<a href="vargs.html#nspinor">nspinor</a>=2).
Moreover, it is not yet implemented in the parallel version of ABINIT.
"""
},
'nconeq': {
'definition': "Number of CONstraint EQuations",
'section': "varrlx",
'category': """<a href="../users/abinit_help.html#no_multi">NO MULTI</a> """,
'vartype': "integer parameter",
'default': "0 ",
'text': """Gives the number of independent equations constraining
the motion of
atoms during structural optimization or molecular dynamics (see <a
href="varrlx.html#natcon">natcon</a>, <a href="varrlx.html#iatcon">iatcon</a>,
and <a href="varrlx.html#wtatcon">wtatcon</a>).
"""
},
'nctime': {
'definition': "NetCdf TIME between output of molecular dynamics informations ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>  """,
'vartype': "integer parameter  ",
'default': "0",
'text': """<p>When <b>nctime</b> is non-zero, the molecular dynamics information
is output in NetCDF format, every <b>nctime</b> time step. Here is the content of an example file :
<pre>
netcdf md32.outH_moldyn1 {
dimensions:
time = UNLIMITED ; // (11 currently)
DimTensor = 6 ;
DimCoord = 3 ;
NbAtoms = 32 ;
DimVector = 3 ;
DimScalar = 1 ;
variables:
double E_pot(time) ;
E_pot:units = "hartree" ;
double E_kin(time) ;
E_kin:units = "hartree" ;
double Stress(time, DimTensor) ;
Stress:units = "hartree/Bohr^3" ;
double Position(time, DimCoord, NbAtoms) ;
Position:units = "Bohr" ;
double Celerity(time, DimCoord, NbAtoms) ;
Celerity:units = "Bohr/(atomic time unit)" ;
double PrimitiveVector1(DimVector) ;
double PrimitiveVector2(DimVector) ;
double PrimitiveVector3(DimVector) ;
double Cell_Volume(DimScalar) ;
Cell_Volume:units = "Bohr^3" ;
}
</pre>

"""
},
'ndivk': {
'definition': "Number of DIVisions of K lines",
'section': "vargs",
'category': """<a href="../users/abinit_help.html#not_internal">NOT INTERNAL</a> """,
'vartype': """integer array ndivk(abs(<a href="varbas.html#kptopt">kptopt</a>))  """,
'default': "No default.",
'text': """Gives the number of divisions of each of the segments
of the band structure, whose path is determined by
<a href="varbas.html#kptopt">kptopt</a>
and
<a href="vargs.html#kptbounds">kptbounds</a>.
This is only needed when <a href="varbas.html#kptopt">kptopt</a>
is negative. In this case, the absolute value of
<a href="varbas.html#kptopt">kptopt</a> is the number of such segments.
<p>For example, suppose that the number of segment is just one
(<a href="varbas.html#kptopt">kptopt</a>=-1),
a value <b>ndivk</b>=4 will lead to the computation
of points with relative coordinates 0.0, 0.25, 0.5, 0.75 and 1.0 , along
the segment in consideration.
<p>Now, suppose that there are two segments
(<a href="varbas.html#kptopt">kptopt</a>=-2), with
<b>ndivk</b>(1)=4 and <b>ndivk</b>(2)=2, the computation of the
eigenvalues will be done at 7 points, 5 belonging to the
first segment, with relative coordinates 0.0, 0.25, 0.5, 0.75 and 1.0,
the last one being also the starting point of the next segment,
for which two other points must be computed, with relative coordinates
0.5 and 1.0 .
<p>It is easy to compute disconnected circuits (non-chained segments),
by separating
the circuits with the value <b>ndivk</b>=1 for the intermediate
segment connecting the end of one circuit with the
beginning of the next one (in which case no intermediate
point is computed along this segment).
<p>Alternatively it is possible to generate automatically the array <b>ndivk</b>
by just specifying the number of divisions for the smallest segment.
See the related input variable <a href="vargs.html#ndivsm">ndivsm</a>.
"""
},
'ndivsm': {
'definition': "Number of DIVisions for the SMallest segment",
'section': "vargs",
'category': """<a href="../users/abinit_help.html#not_internal">NOT INTERNAL</a> """,
'vartype': "integer ",
'default': "No default.",
'text': """Gives the number of divisions used to sample the smallest segment
of the circuit employed during a band structure calculation
(see related input variables
<a href="varbas.html#kptopt">kptopt</a>
and
<a href="vargs.html#kptbounds">kptbounds</a>).
If <b>ndivsm</b> is given in the input file, there is no need to specify the number of divisions
to be used for the other segments.
Indeed <a href="vargs.html#ndivk">ndivk</a> is automatically calculated inside the
code in order to generate a path where the number of divisions in each segment is proportial
to the length of the segment itself.
This option is activated only when <a href="varbas.html#kptopt">kptopt</a> is negative.
In this case, the absolute value of
<a href="varbas.html#kptopt">kptopt</a> is the number of such segments.
"""
},
'ndtset': {
'definition': "Number of DaTaSETs ",
'section': "varbas",
'category': """<a href="../users/abinit_help.html#no_multi">NO MULTI</a>  """,
'vartype': "integer parameter  ",
'default': "0 (no multi-data set).",
'text': """Gives the number of data sets to be
treated.
<br>If 0, means that the multi-data set treatment is not used,
so that the root filenames will not be appended with _DSx,
where 'x' is the dataset index defined
by the input variable <a href="varbas.html#jdtset">jdtset</a>,
and also that input names with a dataset index are not allowed.
Otherwise, <b>ndtset</b>=0 is equivalent to <b>ndtset</b>=1."""
},
'ndynimage': {
'definition': "Number of DYNamical IMAGEs ",
'section': "varint",
'category': """<a href="../users/abinit_help.html#internal">INTERNAL</a> """,
'vartype': "integer ",
'default': "",
'text': """This internal variable gives the number of dynamical images,
immediately deduced from the number of non-zero values present in
<a href="varrlx.html#dynimage">dynimage</a>.
It is used to dimension many memory-consuming arrays (one copy for each image),
in case they are not stored on disk (see <a href="varfil.html#mkmem">mkmem</a>),
e.g. the wavefunction array (cg), the density array (rho), etc .
"""
},
'ndyson': {
'definition': "Number of points to be added for the solution of the DYSON equation ",
'section': "vardev",
'category': "",
'vartype': "integer parameter ",
'default': "-1.",
'text': """Number of points to be added to lambda=0 and lambda=1 (that are always calculated
for the integration ober the coupling constant lambda in the ACFD calculation of the
exchange-correlation energy.
<ul>
<li>ndyson=-1 : let the code decide how many points to use (presently, 3 points
for <a href="vardev.html#idyson">idyson</a>=1 or 3, and 9 points for
<a href="vardev.html#idyson">idyson</a>=2)</li>
<li>ndyson=0 : only compute the non-interacting and fully-interacting
susceptibility matrices.</li>
<li>ndyson&gt;0 : use <b>ndyson</b> more points in ]0,1[</li>
</ul>
"""
},
'neb_algo': {
'definition': "Nudged Elastic Band ALGOrithm ",
'section': "varrlx",
'category': "",
'vartype': "integer parameter",
'default': "<b>neb_algo</b>=1",
'text': """Relevant only when <a href="varrlx.html#imgmov">imgmov</a>=5 (Nudged Elastic Band).<br>
Gives the variant of the NEB method used.<br>
Possible values can be:<br>
<ul>
<li>0=&gt; <b>Original NEB method</b>.<br>
<i>See: Classical and Quantum Dynamics in Condensed Phase Simulations, edited by
Berne, Ciccotti, Coker (World Scientific, Singapore, 1998), pp. 385-404</i>
</li>
<br>
<li>1=&gt; <b>NEB + improved tangent</b>.<br>
The Improved Tangent Method builds on the NEB with an improved estimate of the
tangent direction and a resulting change of the component of the spring force acting
on the images.<br>
<i>See: J. Chem. Phys. 113, 9978 (2000).</i>
</li>
<br>
<li>2=&gt; <b>Climbing-Image NEB (CI-NEB)</b>.<br>
The CI-NEB method constitutes a small modification to the NEB method.
Information about the shape of the MEP is retained, but a rigorous
convergence to a saddle point is also obtained.
By defaut the spring constants are variable (see <a href="varrlx.html#neb_spring">neb_spring</a>).
As the image with the highest energy has to be identified,
the calculation begins with several iterations of the standard NEB
algorithm. The effective CI-NEB begins at the <a href="varrlx.html#cineb_start">cineb_start</a> iteration.<br>
<i>See: J. Chem. Phys. 113, 9901 (2000).</i>
</li>
</ul>
Note that, in all cases, it is possible to define the value of the spring constant connecting images with
<a href="varrlx.html#neb_spring">neb_spring</a>, keeping it constant or allowing it
to vary between 2 values (to have higher resolution close to the saddle point).
"""
},
'neb_spring': {
'definition': "Nudged Elastic Band: SPRING constant",
'section': "varrlx",
'category': "",
'vartype': "real array neb_spring(2)",
'default': """<b>neb_spring</b>=2*0.05 Ha/bohr^2 when <a href="varrlx.html#neb_algo">neb_algo</a>/=2, <b>neb_spring</b>=(0.02, 0.15) Ha/bohr^2 when <a href="varrlx.html#neb_algo">neb_algo</a>/=2 (CI-NEB)""",
'text': """Relevant only when <a href="varrlx.html#imgmov">imgmov</a>=5 (Nudged Elastic Band).<br>
Gives the minimal and maximal values of the spring constant connecting images for the NEB method.<br>
In the standard "Nudged Elastic Band" method, the spring constant is constant along the path,
but, in order to have higher resolution close to the saddle point, it can be better
to have stronger springs close to it.<br>
Units: Hartree/Bohr^2.<br>
<i>See: J. Chem. Phys. 113, 9901 (2000).</i>
"""
},
'nelect': {
'definition': "Number of ELECTrons ",
'section': "varint",
'category': """<a href="../users/abinit_help.html#internal">INTERNAL</a> """,
'vartype': "real number ",
'default': "",
'text': """This internal variable gives the number of electrons per unit
cell, as computed from the sum of the valence electrons
related to each atom (given in the pseudopotential, where it is called
"zion"), and the input variable
<a href="vargs.html#charge">charge</a>:
<br>
<b>nelect</b>=zion-<a href="vargs.html#charge">charge</a>.
"""
},
'nfft': {
'definition': "Number of FFT points ",
'section': "varint",
'category': """<a href="../users/abinit_help.html#internal">INTERNAL</a> """,
'vartype': "integer ",
'default': "",
'text': """If space parallelisation is not used,
this internal variable gives the number of Fast Fourier Transform
points in the grid generated by
<a href="vargs.html#ngfft">ngfft</a>(1:3). It is simply the
product of the three components of <a href="vargs.html#ngfft">ngfft</a>.
<br> If space parallelisation is used, then it becomes the
number of Fast Fourier Transform points attributed to the
particular processor. It is no longer the above-mentioned simple product,
but a number usually close to this product divided by the
number of processors on which the space is shared.
"""
},
'nfftdg': {
'definition': "Number of FFT points for the Double Grid",
'section': "varint",
'category': """<a href="../users/abinit_help.html#internal">INTERNAL</a> """,
'vartype': "integer ",
'default': "",
'text': """If space parallelisation is not used,
this internal variable gives the number of Fast Fourier Transform
points in the (double) grid generated by
<a href="varpaw.html#ngfftdg">ngfftdg</a>(1:3). It is simply the
product of the three components of <a href="varpaw.html#ngfftdg">ngfftdg</a>.
<br> If space parallelisation is used, then it becomes the
number of Fast Fourier Transform points attributed to the
particular processor. It is no longer the above-mentioned simple product,
but a number usually close to this product divided by the
number of processors on which the space is shared.
"""
},
'nfreqim': {
'definition': "Number of FREQuencies along the IMaginary axis ",
'section': "vargw",
'category': " ",
'vartype': "integer parameter  ",
'default': "0",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=3, that is, screening calculations.
<p>
<b>nfreqim</b> is used only for numerical integration of the GW self-energy
(<a href="vargw.html#gwcalctyp">gwcalctyp</a>= 2, 12, 22, 9, 19, 29).
<br>
<b>nfreqim</b> sets the number of pure imaginary frequencies used to calculate
the dielectric matrix in order to perform the numerical integration of the GW self-energy.
"""
},
'nfreqmidm': {
'definition': "Nth FREQuencey moment of the IMaginary part of DM ",
'section': "vargw",
'category': " ",
'vartype': "integer parameter  ",
'default': "No default ",
'text': """at the moment only relevant if <a href="vargs.html#optdriver">optdriver</a>=4,
that is, a screening calculation.
<p>
depending on the value of <b>>nfreqmidm</b> will calculate the frequency moment of the Dielectric matrix or its inverse,

if <b>>nfreqmidm</b> positive : calculate (nth=nfreqmidm) frequency moment of the Dielectric matrix
if <b>>nfreqmidm</b> negative : calculate (nth=nfreqmidm) frequency moment of the inverse Dielectric matrix
if <b>>nfreqmidm</b> = 0 : calculate first frequency moment of the full polarizability

see M. Taut, J. Phys. C: Solid State Phys. 18 (1985) 2677-2690.

"""
},
'nfreqre': {
'definition': "Number of FREQuencies along the Real axis ",
'section': "vargw",
'category': " ",
'vartype': "integer parameter  ",
'default': "0",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=3, that is, screening calculations.
<p>
<b>nfreqre</b> is used only for numerical integration of the GW self-energy
(<a href="vargw.html#gwcalctyp">gwcalctyp</a>= 2, 12, 22, 9, 19, 29).
<br>
<b>nfreqre</b> sets the number of real frequencies used to calculate
the dielectric matrix in order to perform the numerical integration of the GW self-energy.
<p>
It can be used also in case of GW calculations with plasmon-pole models,
<i>i.e</i> <a href="vargw.html#gwcalctyp">gwcalctyp</a>&lt;10,
to reduce the number of frequencies used to evaluate the dielectric matrix from the (default) two to
one frequency (omega=0) by setting <b>nfreqre</b>=1.
This might be a good idea in case one is planning to use ppmodel&gt;1.
This will force the calculation of the screening on a single frequency (omega=0) and hence
reduce memory and disk space requirement.
The only draw back is that the user will not be able to perform self energy calculation using
<a href="vargw.html#ppmodel">ppmodel</a>=1,
since in the last case the dielectric matrix calculated on two frequencies is required.
If the user is not sure which ppmodel to use, then s/he is not advised
to use this input variable.
Using the default values, one must be able to get a screening file that can be used with any ppmodel.
"""
},
'nfreqsp': {
'definition': "Number of FREQuencies for the SPectral function ",
'section': "vargw",
'category': " ",
'vartype': "integer parameter  ",
'default': "0",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=4, that is, sigma calculations.
<p>
<b>nfreqsp</b> defines the number of real frequencies used to calculate the spectral function
of the GW Green's function.
"""
},
'nfreqsus': {
'definition': "Number of FREQuencies for the SUSceptibility matrix",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>  """,
'vartype': "integer parameter  ",
'default': "0",
'text': """If 0, no computation of frequency-dependent susceptibility matrix.
If 1 or larger, will read <a href="vardev.html#freqsuslo">freqsuslo</a> and
<a href="vardev.html#freqsusin">freqsusin</a>
to define the frequencies
(1 is currently the only value allowed)
"""
},
'ngfft': {
'definition': "Number of Grid points for Fast Fourier Transform ",
'section': "vargs",
'category': " ",
'vartype': "integer array ngfft(3)",
'default': "0 0 0 (so, automatic selection of optimal values)",
'text': """gives the size of fast Fourier transform
(fft) grid in three dimensions.  Each number must be
composed of the factors 2, 3, and 5 to be consistent with
the radices available in our fft.  If no <b>ngfft</b> is provided or
if <b>ngfft</b> is set to 0 0 0, the code will automatically provide
an optimal set of <b>ngfft</b> values, based on <a href="varbas.html#acell">acell</a>,
<a href="varbas.html#rprim">rprim</a> and <a href="varbas.html#ecut">ecut</a>
(see also <a href="vargs.html#boxcutmin">boxcutmin</a> for speed/accuracy concerns).
This is the recommended procedure, of course.
<br>The total number of FFT points
is the product:<br> <tele><b>ngfft</b>(1)*<b>ngfft</b>(2)*<b>ngfft</b>(3)=nfft </tele>. <br>
When <b>ngfft</b> is made smaller
than recommended values (e.g. by setting <a href="vargs.html#boxcutmin">boxcutmin</a>
to a value smaller than 2.0 or by setting <b>ngfft</b> manually), the code runs faster and the
equations in effect are approximated by a low pass Fourier
filter.  The code reports to standard output (unit 06) a
parameter "boxcut" which is the smallest ratio of the fft
box side to the G vector basis sphere diameter.  When
boxcut is less than 2 the Fourier filter approximation is being
used.  When boxcut gets less than about 1.5 the
approximation may be too severe for realistic results
and should be tested against larger values of <b>ngfft</b>.
When boxcut is larger than 2, <b>ngfft</b> could be reduced without
loss of accuracy. In this case, the small variations
that are observed are solely due to the
xc quadrature, that may be handled with <a href="vardev.html#intxc">intxc</a>=1
to even reduce this effect.
<p>
Internally, <b>ngfft</b> is an array of size 18. The present
components are stored in <b>ngfft</b>(1:3), while
<ul>
<li> <b>ngfft</b>(4:6) contains slightly different (larger) values,
modified for efficiency of the FFT </li>
<li> <b>ngfft</b>(7) is <a href="vardev.html#fftalg">fftalg</a> </li>
<li> <b>ngfft</b>(8) is <a href="vardev.html#fftcache">fftcache</a> </li>
<li> <b>ngfft</b>(9) is set to 0 if the parallelization of the FFT is not activated,
while it is set to 1 if it is activated.</li>
<li> <b>ngfft</b>(10) is the number of processors of the FFT group </li>
<li> <b>ngfft</b>(11) is the index of the processor in the group of processors </li>
<li> <b>ngfft</b>(12) is n2proc, the number of x-z planes, in reciprocal space, treated by the processor </li>
<li> <b>ngfft</b>(13) is n3proc, the number of x-y planes, in real space, treated by the processor </li>
<li> <b>ngfft</b>(14) is mpi_comm_fft, the handle on the MPI communicator in charge of the FFT parallelisation </li>
<li> <b>ngfft</b>(15:18) are not yet used </li>
</ul>
<br>
The number of points stored by this processor in real space is n1*n2*n3proc, while in reciprocal
space, it is n1*n2proc*n3.
"""
},
'ngfftdg': {
'definition': "Number of Grid points for Fast Fourier Transform : Double Grid ",
'section': "varpaw",
'category': "",
'vartype': "integer array ngfftdg(3) ",
'default': "0 0 0 (so, automatic selection of optimal values)",
'text': """Relevant only when
<a href="varint.html#usepaw">usepaw</a>=1.
<br> This variable has the same meaning as ngfft (gives the size of fast Fourier
transform (fft) grid in three dimensions) but concerns the "double
grid" only used for PAW calculations.
"""
},
'ngkpt': {
'definition': "Number of Grid points for K PoinTs generation ",
'section': "varbas",
'category': """<a href="../users/abinit_help.html#not_internal">NOT INTERNAL</a> """,
'vartype': "integer array ngkpt(3)",
'default': "No Default",
'text': """Used when <a href="varbas.html#kptopt">kptopt</a>>=0,
if <a href="vargs.html#kptrlatt">kptrlatt</a>
has not been defined (<a href="vargs.html#kptrlatt">kptrlatt</a>
and <b>ngkpt</b> are exclusive of each other).
<br> Its three positive components
give the number of k points of Monkhorst-Pack grids
(defined with respect to primitive axis in reciprocal space)
in each of the three dimensions.
<b>ngkpt</b> will be used to generate the
corresponding <a href="vargs.html#kptrlatt">kptrlatt</a>
input variable.
The use of <a href="varbas.html#nshiftk">nshiftk</a>
and <a href="varbas.html#shiftk">shiftk</a>, allows to generate
shifted grids, or Monkhorst-Pack grids defined
with respect to conventional unit cells.
<p>When <a href="varbas.html#nshiftk">nshiftk</a>=1,
<a href="vargs.html#kptrlatt">kptrlatt</a> is initialized
as a diagonal (3x3) matrix, whose diagonal elements
are the three values <b>ngkpt</b>(1:3). When
<a href="varbas.html#nshiftk">nshiftk</a> is greater than 1,
ABINIT will try to generate <a href="vargs.html#kptrlatt">kptrlatt</a>
on the basis of the primitive vectors of the k-lattice:
the number of shifts might be reduced, in which case
<a href="vargs.html#kptrlatt">kptrlatt</a> will not be diagonal
anymore.
<p>Monkhorst-Pack grids are usually the most efficient when
their defining integer numbers are even.
For a measure of the efficiency, see the input variable
<a href="vargs.html#kptrlen">kptrlen</a>.
"""
},
'ngqpt': {
'definition': "Number of Grid pointsfor Q PoinTs generation ",
'section': "vargs",
'category': """<a href="../users/abinit_help.html#not_internal">NOT INTERNAL</a> """,
'vartype': "integer array ngqpt(3)",
'default': "(0 0 0)",
'text': """Used when <a href="vargs.html#nqpt">nqpt</a>=1, with <a href="vargs.html#qptopt">kptopt</a>>=0,
if <a href="vargs.html#qptrlatt">qptrlatt</a>
has not been defined (<a href="vargs.html#qptrlatt">qptrlatt</a>
and <b>ngqpt</b> are exclusive of each other).
<br>At variance with <a href="varbas.html#ngkpt">ngkpt</a>, note that only one q point
is selected per dataset (see <a href="vargs.html#iqpt">iqpt</a>).
<br> Its three positive components
give the number of q points of Monkhorst-Pack grids
(defined with respect to primitive axis in reciprocal space)
in each of the three dimensions.
The use of <a href="vargs.html#nshiftq">nshiftq</a>
and <a href="vargs.html#shiftq">shiftq</a>, allows to generate
shifted grids, or Monkhorst-Pack grids defined
with respect to conventional unit cells.
<p>For more information on Monkhorst-Pack grids, see <a href="varbas.html#ngkpt">ngkpt</a>.
"""
},
'nimage': {
'definition': "Number of IMAGEs",
'section': "varrlx",
'category': "",
'vartype': "integer parameter",
'default': "1",
'text': """Give the number of images (or replicas) of the system,
for which the forces and stresses might be computed independantly,
in the context of the string method, the genetic algorithm, hyperdynamics or Path-Integral Molecular Dynamics
depending on the value of <a href="varrlx.html#imgmov">imgmov</a>).
Related input variables : <a href="varrlx.html#dynimage">dynimage</a>,
<a href="varpar.html#npimage">npimage</a>, <a href="varrlx.html#ntimimage">ntimimage</a>
and <a href="varfil.html#prtvolimg">prtvolimg</a>.
<br>
Images might differ by the position of atoms in the unit cell, their
velocity, as well as by their cell geometry. The following input variables might be used to define
the images :
<ul>
<li><a href="varbas.html#acell">acell</a>
<li><a href="varbas.html#angdeg">angdeg</a>
<li><a href="varbas.html#rprim">rprim</a>
<li><a href="varrlx.html#vel">vel</a>
<li><a href="varbas.html#xangst">xangst</a>
<li><a href="varbas.html#xcart">xcart</a>
<li><a href="varbas.html#xred">xred</a>
</ul>
These input variables, non-modified, will be used to define the image with index 1.
For the image with the last index, the input file might specify the values
of such input variables, appended with "_lastimg", e.g. one of these :
<ul>
<li><a href="varbas.html#acell">acell_lastimg</a>
<li><a href="varbas.html#angdeg">angdeg_lastimg</a>
<li><a href="varbas.html#rprim">rprim_lastimg</a>
<li><a href="varrlx.html#vel">vel_lastimg</a>
<li><a href="varbas.html#xangst">xangst_lastimg</a>
<li><a href="varbas.html#xcart">xcart_lastimg</a>
<li><a href="varbas.html#xred">xred_lastimg</a>
</ul>

By default, these values will be interpolated linearly to define values for the other images, unless
there exist specific values for some images, for which the string "last" has to be replaced by the
index of the image, e.g. for the image number 4 :
<ul>
<li><a href="varbas.html#acell">acell_4img</a>
<li><a href="varbas.html#angdeg">angdeg_4img</a>
<li><a href="varbas.html#rprim">rprim_4img</a>
<li><a href="varrlx.html#vel">vel_4img</a>
<li><a href="varbas.html#xangst">xangst_4img</a>
<li><a href="varbas.html#xcart">xcart_4img</a>
<li><a href="varbas.html#xred">xred_4img</a>
</ul>

It is notably possible to specify the starting point and the end point of the path (of images), while specifying intermediate points.<br>
<br>
It usually happen that the images do not have the same symmetries and space group.
ABINIT has not been designed to use different set of symmetries for different images.
ABINIT will use the symmetry and space group of the image number 2, that
is expected to have a low number of symmetries.
This might lead to erroneous calculations, in case some image has even less
symmetry. By contrast, there is no problem if some other image has more symmetries
than those of the second image.
"""
},
'nkpt': {
'definition': "Number of K - Points ",
'section': "varbas",
'category': " ",
'vartype': "integer parameter  ",
'default': """0 if <a href="varbas.html#kptopt">kptopt</a>/=0, and 1 if <a href="varbas.html#kptopt">kptopt</a>==0.""",
'text': """If non-zero,
<b>nkpt</b>
gives the number of k points in the k point array
<a href="varbas.html#kpt">kpt</a>. These points
are used either to sample the Brillouin zone, or to
build a band structure along specified lines.
<p>
If <b>nkpt</b> is zero, the code deduces from other input variables
(see the list in the description of <a href="varbas.html#kptopt">kptopt</a>)
the number of k points, which is possible only
when <a href="varbas.html#kptopt">kptopt</a>/=0.
If <a href="varbas.html#kptopt">kptopt</a>/=0 and
the input value of <a href="varbas.html#nkpt">nkpt</a>/=0,
then ABINIT will check that the number of k points
generated from the other input variables
is exactly the same than <b>nkpt</b>.

<p>If <a href="varbas.html#kptopt">kptopt</a> is positive,
<b>nkpt</b> must be coherent with the values
of <a href="vargs.html#kptrlatt">kptrlatt</a>,
<a href="varbas.html#nshiftk">nshiftk</a>
and <a href="varbas.html#shiftk">shiftk</a>.
<br>
For ground state calculations, one should select the
k point in the irreducible Brillouin Zone (obtained
by taking into account point symmetries and the time-reversal
symmetry).
<br>For response function calculations, one should
select k points in the full Brillouin zone, if the wavevector
of the perturbation does not vanish, or in a half of
the Brillouin Zone if q=0. The code will automatically decrease
the number of k points to the minimal set needed for
each particular perturbation.
<p>If <a href="varbas.html#kptopt">kptopt</a> is negative,
<b>nkpt</b> will be the sum of the number of points on
the different lines of the band structure.
For example,
if <a href="varbas.html#kptopt">kptopt</a>=-3, one
will have three segments; supposing
<a href="vargs.html#ndivk">ndivk</a> is 10 12 17,
the total number of k points of the circuit will be
10+12+17+1(for the final point)=40.
"""
},
'nkptgw': {
'definition': "Number of K-Points for GW corrections ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>  """,
'vartype': "integer parameter",
'default': "0 ",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=4, that is, sigma calculations.
This input variable was called "ngwpt" in versions before v4.3.
<p>
<b>nkptgw</b> gives the number of k-points for which the GW calculation must be done.
It is used to dimension <a href="vargw.html#kptgw">kptgw</a>
"""
},
'nline': {
'definition': "Number of LINE minimisations ",
'section': "vargs",
'category': " ",
'vartype': "integer parameter  ",
'default': "4.",
'text': """Gives maximum number of line minimizations
allowed in preconditioned conjugate gradient minimization
for each band. The Default, 4, is fine.
<br>Special cases, with degeneracies or near-degeneracies
of levels at the Fermi energy may require a larger value of
<b>nline</b> (5 or 6 ?)
Line minimizations will be stopped anyway when improvement
gets small. With the input variable <a href="vardev.html#nnsclo">nnsclo</a>,
governs the convergence of the wavefunctions
for fixed potential.
<br>Note that <b>nline</b>=0 can be used to diagonalize the Hamiltonian
matrix in the subspace spanned by the input wavefunctions."""
},
'nloalg': {
'definition': "Non Local ALGorithm ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a> """,
'vartype': "integer variable  ",
'default': "4 (norm-conserving psps) or 14 (PAW), except for the NEC where it is 2 (or 12).",
'text': """Allows to choose the algorithm
for non-local operator application.
On super-scalar architectures, the Default <b>nloalg</b>=4/14 is the best,
but you can save memory by using <b>nloalg</b>=-4.<br>
More detailed explanations:<br><br>
<div style="margin-left: 10px;">Units figure of <b>nloalg</b>:</div>
<div style="margin-left: 40px;">
- <b>nloalg</b>=?2 : Should be efficient on vector machines. It is
indeed the fastest algorithm for the NEC, but
actual tests on Fujitsu machine did not gave better
performances than the other options.<br>
- <b>nloalg</b>=?3 : same as <b>nloalg</b>==2, but the loop order is inverted.<br>
- <b>nloalg</b>=?4 : same as <b>nloalg</b>==3, but maximal use of registers
has been coded. This should be especially efficient on
scalar and super-scalar machines. This has been
confirmed by tests.<br>
</div><br>
<div style="margin-left: 10px;">Tens figure of <b>nloalg</b>:</div>
<div style="margin-left: 40px;">
- <b>nloalg</b><10 : (k+G) vectors are not precomputed, in order to save memory space.<br>
- <b>nloalg</b>>=10 : (k+G) vectors are precomputed, once per k-point.<br>
</div><br>
<div style="margin-left: 10px;">Sign of <b>nloalg</b>:</div>
<div style="margin-left: 40px;">
Negative values of <b>nloalg</b> correspond positive ones,
where the phase precomputation has been suppressed,
in order to save memory space: an array <i>double precision :: ph3d(2,npw,<a href="varbas.html#natom">natom</a>)</i>
is saved (typically half the space needed
for the wavefunctions at 1 k point - this corresponds
to the silicon case). However, the computation of phases
inside nonlop is somehow time-consuming.<br>
</div>
<br>
Note: internally, <b>nloalg</b> is an array <i>nloalg(1:5)</i>,
that also allows to initialize several internal variables (not documented):
<div style="margin-left: 40px;">
- <i>nloalg(1)</i>=mod(<b>nloalg</b>,10)<br>
- <i>jump</i>=nloalg(2)<br>
- <i>mblkpw</i>=nloalg(3)<br>
- <i>mincat</i>=nloalg(4)<br>
- <i>nloalg(5)</i>=<b>nloalg</b>/10<br>
</div>
However, only <i>nloalg(1)+10*nloalg(5)</i> is read as an input variable."""
},
'nnos': {
'definition': "Number of nose masses",
'section': "varrlx",
'category': "",
'vartype': "integer parameter",
'default': "0",
'text': """Give the number of oscillators in the Martyna et al.
chain of oscillators thermostats (when <a href="varrlx.html#ionmov">ionmov</a>=13).
The mass of these oscillators is given by qmass (similar to <a
href="varrlx.html#noseinert">noseinert</a>)&nbsp;
"""
},
'nnsclo': {
'definition': "Number of Non-Self Consistent LOops ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a> """,
'vartype': "integer parameter ",
'default': "0.",
'text': """Gives the maximum number of
non-self-consistent loops of <a href="vargs.html#nline">nline</a> line minimisations,
in the SCF case (when <a href="varbas.html#iscf">iscf</a> &gt;0).  In the case <a href="varbas.html#iscf">iscf</a> &lt;=0 ,
the number of non-self-consistent loops is determined
by <a href="varbas.html#nstep">nstep</a>.
<br>The Default value of 0 correspond to make
the two first fixed potential determinations
of wavefunctions have 2 non-self consistent loops,
and the next ones to have only 1 non-self consistent loop. """
},
'nobj': {
'definition': "Number of OBJects ",
'section': "vargeo",
'category': """<a href="../users/abinit_help.html#geometry_builder">GEOMETRY BUILDER</a>, <a href="../users/abinit_help.html#not_internal">NOT INTERNAL</a> """,
'vartype': "integer parameter ",
'default': "0 (no use of the geometry builder).",
'text': """Gives the number of 'objects'
to be used by the geometry builder in order to find the full
set of atoms. At present, only one or two objects
can be defined, identified as objects 'a' and 'b'.
<br>Related variables for object 'a' are : <a href="vargeo.html#objan">objan</a>,
<a href="vargeo.html#objaat">objaat</a>,
<a href="vargeo.html#objarf">objarf</a>, <a href="vargeo.html#objatr">objatr</a>,
<a href="vargeo.html#objaro">objaro</a>, <a href="vargeo.html#objaax">objaax</a>. Related variables
for object 'b' are : <a href="vargeo.html#objan">objbn</a>, <a href="vargeo.html#objaat">objbat</a>,
<a href="vargeo.html#objarf">objbrf</a>, <a href="vargeo.html#objatr">objbtr</a>,
<a href="vargeo.html#objaro">objbro</a>, <a href="vargeo.html#objaax">objbax</a>.
<br>More detailed explanation : when the geometry builder
is used (i.e. when <b>nobj</b>==1 or <b>nobj</b>==2), the code
will be given a primitive set of atoms, from which it
will have to deduce the full set of atoms.
<br>An object will be specified by the number of
atoms it includes (<a href="vargeo.html#objan">objan</a> or <a href="vargeo.html#objan">objbn</a>),
and the list of these atoms (<a href="vargeo.html#objaat">objaat</a> or <a href="vargeo.html#objaat">objbat</a>).
<br>Examples of physical realisation of an object can be a molecule,
or a group of atom to be repeated, or a part of a molecule
to be rotated.
The geometry builder can indeed repeat these objects (<a href="vargeo.html#objarf">objarf</a>
or <a href="vargeo.html#objarf">objbrf</a>), rotate them (<a href="vargeo.html#objaro">objaro</a> or
<a href="vargeo.html#objaro">objbro</a>) with respect
to an axis (<a href="vargeo.html#objaax">objaax</a> or <a href="vargeo.html#objaax">objbax</a>), and translate them
(<a href="vargeo.html#objatr">objatr</a> or <a href="vargeo.html#objatr">objbtr</a>).
After having generated a geometry
thanks to rotation, translation and repetition of objects,
it is possible to remove some atoms, in order to create
vacancies (<a href="vargeo.html#vacnum">vacnum</a> and <a href="vargeo.html#vaclst">vaclst</a>).
The number of atoms in the primitive
set, those that will be read from the input file, is
specified by the variable <a href="vargeo.html#natrd">natrd</a>. It will be always smaller
than the final number of atoms, given by the variable <a href="varbas.html#natom">natom</a>.
The code checks whether the primitive number of atoms
plus those obtained by the repetition operation is
coherent with the variable <a href="varbas.html#natom">natom</a>, taking into
account possible vacancies. <br>You should look
at the other variables for more information.
Go to <a href="vargeo.html#objan">objan</a>, for example.
<br>Not present in the dtset array (no internal).
"""
},
'nomegasf': {
'definition': "Number of OMEGA to evaluate the Spectral Function ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>, <a href="../users/abinit_help.html#energy">ENERGY</a> """,
'vartype': "integer parameter ",
'default': "0 ",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=3 and
<a href="vargw.html#spmeth">spmeth</a>/=0, that is, screening calculations
based on the spectral reprentation of the irreducible polarizability.
<p>
<b>nomegasf</b> defines the number of real frequencies used to describe the spectral function
associated to the irreducible polarizability $\chi^{(0)}_{KS}$.
The frequency mesh will cover the interval between 0 and the maximum (positive) transition
energy between occupied and empty states.
The delta function entering the expression defining the spectral function is approximated using two
different methods according to the value of the <a href="vargw.html#spmeth">spmeth</a> input variable.
<p> It is important to notice that an accurate description of the imaginary part of $\chi^{(0)}_{KS}$
requires an extremely dense frequency mesh. It should be kept in mind, however, that the memory required
grows fast with the value of <b>nomegasf</b>.
"""
},
'nomegasi': {
'definition': "Number of OMEGA(S) along the Imaginary axis ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>  """,
'vartype': "integer ",
'default': "12 ",
'text': """<p>

Only relevant for sigma calculations in which the self-energy along the real axis is obtained
by performing the analytic continuation from the imaginary axis to the full complex plane
via the Pade approximant
(<a href="vargs.html#optdriver">optdriver</a>=4 and <a href="vargw.html#gwcalctyp">gwcalctyp</a>=x1).
<p>
<b>nomegasi</b> defines the number of frequency points used to sample the self-energy along the
imaginary axis. The frequency mesh is linear and covers the interval between OMEGASIMIN=0.01 Hartree and
<a href="vargw.html#omegasimax">omegasimax</a>.
"""
},
'nomegasrd': {
'definition': "Number of OMEGA to evaluate the Sigma Real axis Derivative ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>  """,
'vartype': "integer parameter ",
'default': "9 ",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=4, that is, sigma calculations.
<p>
The number of real frequencies around the KS energy where the self-energy Sigma is evaluated.
From these values, the derivative of Sigma at the KS energy is numerically estimated through linear interpolation.
"""
},
'normpawu': {
'definition': "NORMalize atomic PAW+U projector ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a> """,
'vartype': """integer normpawu(<a href="varbas.html#ntypat">ntypat</a>) """,
'default': "0 ",
'text': """Defines whether the atomic wave function (used as projectors in PAW+U) should be renormalized to 1
within PAW sphere.
<ul>
<li><b>normpawu</b>=0 : leave projector</li>
<li><b>normpawu</b>=1 : renormalize </li>
</ul>
"""
},
'noseft': {
'definition': "",
'section': "vardev",
'category': " ",
'vartype': "",
'default': "",
'text': """TO BE DOCUMENTED
"""
},
'noseinert': {
'definition': "NOSE INERTia factor ",
'section': "varrlx",
'category': "",
'vartype': "real noseinert",
'default': "1.0d5 ",
'text': """Give the inertia factor W<sub>T</sub> of the
Nose-Hoover thermostat
(when <a href="varrlx.html#ionmov">ionmov</a>=8), in atomic
units of weight*length<sup>2</sup>, that is (electron mass)*(Bohr)<sup>2</sup>.
The equations of motion are :<br>
M<sub>I</sub> d<sup>2</sup>R<sub>I</sub>/dt<sup>2</sup>= F<sub>I</sub>
- dX/dt M<sub>I</sub> dR<sub>I</sub>/dt
<br>
and
<br>
W<sub>T</sub> d<sup>2</sup>X/dt<sup>2</sup>=
Sum(I) M<sub>I</sub> (dR<sub>I</sub>/dt)<sup>2</sup>
- 3Nk<sub>B</sub>T
<br>
where I represent each nucleus, M<sub>I</sub> is the mass of each
nucleus
(see <a href="varrlx.html#amu">amu</a>),
R<sub>I</sub> is the coordinate of each nucleus (see <a
href="varbas.html#xcart">xcart</a>),
dX/dt is a dynamical friction coefficient, and T is the temperature
of the thermostat (see <a href="varrlx.html#mdtemp">mdtemp</a>.
"""
},
'np_slk': {
'definition': "Number of mpi Processors used for ScaLapacK calls",
'section': "varpar",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>  """,
'vartype': "integer parameter",
'default': "1000000",
'text': """"""
},
'npband': {
'definition': "Number of Processors at the BAND level ",
'section': "varpar",
'category': "",
'vartype': "integer ",
'default': "1.",
'text': """Relevant only for the band/FFT parallelisation
(see the <a href="varpar.html#paral_kgb">paral_kgb</a> input variable).
<br> <b>npband</b> gives the number of processors among which the work load over the band level is shared.
<b>npband</b>, <a href="varpar.html#npfft">npfft</a>,
<a href="varpar.html#npkpt">npkpt</a> and <a href="varpar.html#npspinor">npspinor</a>
are combined to give the total number
of processors (nproc) working on the band/FFT/k-point parallelisation.<br>
See <a href="varpar.html#npfft">npfft</a>, <a href="varpar.html#npkpt">npkpt</a>,
<a href="varpar.html#npspinor">npspinor</a> and
<a href="varpar.html#paral_kgb">paral_kgb</a> for the additional information on the use of
band/FFT/k-point parallelisation.
<p>
Note : at present, <b>npband</b> has to be a divisor or equal to <a href="varbas.html#nband">nband</a>
"""
},
'npfft': {
'definition': "Number of Processors at the FFT level ",
'section': "varpar",
'category': "",
'vartype': "integer ",
'default': "nproc.",
'text': """Relevant only for the band/FFT/k-point parallelisation
(see the <a href="varpar.html#paral_kgb">paral_kgb</a> input variable).
<br> <b>npfft</b> gives the number of processors among
which the work load over the FFT level is shared.
<b>npfft</b>, <a href="varpar.html#npkpt">npkpt</a>,
<a href="varpar.html#npband">npband</a> and <a href="varpar.html#npspinor">npspinor</a>
are combined to give the total number
of processors (nproc) working on the band/FFT/k-point parallelisation.<br>
See <a href="varpar.html#npband">npband</a>, <a href="varpar.html#npkpt">npkpt</a>,
<a href="varpar.html#npspinor">npspinor</a>, and
<a href="varpar.html#paral_kgb">paral_kgb</a> for the additional information on the use of
band/FFT/k-point parallelisation.
<p>

Note : <a href="vargs.html#ngfft">ngfft</a> is automatically adjusted to <b>npfft</b>.
If the number of processor is changed from a calculation to another one,
<b>npfft</b> may change, and then <a href="vargs.html#ngfft">ngfft</a> also.

"""
},
'npimage': {
'definition': "Number of Processors at the IMAGE level ",
'section': "varpar",
'category': "",
'vartype': "integer ",
'default': """min(nproc,<a href="varint.html#ndynimage">ndynimage</a>) (see below).""",
'text': """Relevant only
when sets of images are activated (see <a href="varrlx.html#imgmov">imgmov</a>
and <a href="varrlx.html#nimage">nimage</a>.
<br> <b>npimage</b> gives the number of processors among
which the work load over the image level is shared. It is compatible with all other parallelization
levels available for ground-state calculations.<br>
Note on the <b>npimage</b> default value: this default value is crude.
It is set to the number of dynamic images (<a href="varint.html#ndynimage">ndynimage</a>)
if the number of available processors allows this choice.
If <a href="varrlx.html#ntimimage">ntimimage</a>=1, <b>npimage</b> is set to
min(nproc,<a href="varrlx.html#nimage">nimage</a>).<br>

<br><br>

<i>See <a href="varpar.html#paral_kgb">paral_kgb</a>,
<a href="varpar.html#npkpt">npkpt</a>,
<a href="varpar.html#npband">npband</a>,
<a href="varpar.html#npfft">npfft</a>
and <a href="varpar.html#npspinor">npspinor</a>
for the additional information on the use of k-point/band/FFT parallelisation.</i>
<p>
"""
},
'npkpt': {
'definition': "Number of Processors at the K-Point Level ",
'section': "varpar",
'category': "",
'vartype': "integer ",
'default': "1.",
'text': """Relevant only for the band/FFT/k-point parallelisation
(see the <a href="varpar.html#paral_kgb">paral_kgb</a> input variable).
<br> <b>npkpt</b> gives the number of processors among
which the work load over the k-point/spin-component level is shared.
<b>npkpt</b>, <a href="varpar.html#npfft">npfft</a>, <a href="varpar.html#npband">npband</a> and
<a href="varpar.html#npspinor">npspinor</a> are combined to give the total number
of processors (nproc) working on the band/FFT/k-point parallelisation.<br>
See <a href="varpar.html#npband">npband</a>, <a href="varpar.html#npfft">npfft</a>,
<a href="varpar.html#npspinor">npspinor</a> and
<a href="varpar.html#paral_kgb">paral_kgb</a> for the additional information on the use of
band/FFT/k-point parallelisation.
<p>

Note : <b>npkpt</b> should be a divisor or equal to with the number of k-point/spin-components
(<a href="varbas.html#nkpt">nkpt</a>*<a href="varbas.html#nsppol">nsppol</a>)
in order to have the better load-balancing and efficiency.
<br>
"""
},
'nppert': {
'definition': "Number of Processors at the PERTurbation level  ",
'section': "varpar",
'category': "can even be specified separately for each dataset, parameter paral_rf is necessary ",
'vartype': "integer ",
'default': "1.",
'text': """This parameter is used in connection to the parallelization over perturbations(see <a href="varpar.html#paral_rf">paral_rf</a> ), for a linear response calculation.
<b>nppert</b> gives the number of processors among which the work load over the perturbation level is shared.
<br>
"""
},
'npsp': {
'definition': "Number of PSeudoPotentials ",
'section': "vargs",
'category': """<a href="../users/abinit_help.html#no_multi">NO MULTI</a> """,
'vartype': "integer parameter ",
'default': """<a href="varbas.html#ntypat">ntypat</a>""",
'text': """Usually, the number of pseudopotentials to be read is equal
to the number of type of atoms. However,
in the case an alchemical mixing of pseudopotential is to be used,
often the number of pseudopotentials to be read will not equal the number of types of atoms.
<p> Alchemical pseudopotentials will be present
when <a href="vargs.html#ntypalch">ntypalch</a> is non-zero.
See <a href="vargs.html#ntypalch">ntypalch</a>
to understand how
to use alchemical potentials in ABINIT.
The input variables
(<a href="vargs.html#ntypalch">ntypalch</a>,
<a href="vargs.html#algalch">algalch</a>,<a href="vargs.html#mixalch">mixalch</a>)
are active, and generate alchemical potentials from the available
pseudopotentials. Also, the inner variables
(<a href="vargs.html#ntyppure">ntyppure</a>,<a href="vargs.html#npspalch">npspalch</a>)
become active. See these input variables, especially
<a href="vargs.html#mixalch">mixalch</a>, to understand how
to use alchemical potentials in ABINIT.
"""
},
'npspalch': {
'definition': """Number of PSeudoPotentials that are "ALCHemical" """,
'section': "vargs",
'category': "Inner ",
'vartype': "integer parameter, non-negative ",
'default': "",
'text': """Gives the number of pseudopotentials that are used for alchemical mixing (when <a href="vargs.html#ntypalch">ntypalch</a> is non-zero) :
<p> <b>npspalch</b>=<a href="vargs.html#npsp">npsp</a>-<a href="vargs.html#ntyppure">ntyppure</a>
"""
},
'npspinor': {
'definition': "Number of Processors at the SPINOR level ",
'section': "varpar",
'category': "",
'vartype': "integer ",
'default': "1.",
'text': """Can be 1 or 2 (if <a href="vargs.html#nspinor">nspinor</a>=2).<br>
Relevant only for the band/FFT/k-point parallelisation
(see the <a href="varpar.html#paral_kgb">paral_kgb</a> input variable).
<br> <b>npspinor</b> gives the number of processors among
which the work load over the spinorial components of wave-functions is shared.
<b>npspinor</b>, <a href="varpar.html#npfft">npfft</a>,
<a href="varpar.html#npband">npband</a> and <a href="varpar.html#npkpt">npkpt</a>
are combined to give the total number of processors (nproc)
working on the band/FFT/k-point parallelisation.<br>
See <a href="varpar.html#npkpt">npkpt</a>, <a href="varpar.html#npband">npband</a>,
<a href="varpar.html#npfft">npfft</a>, and
<a href="varpar.html#paral_kgb">paral_kgb</a> for the additional information on the use of
band/FFT/k-point parallelisation.
<br>
"""
},
'npulayit': {
'definition': "Number of PULAY ITerations for SC mixing ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a> """,
'vartype': "integer parameter ",
'default': "7.",
'text': """Gives the  number of previous iterations involved in Pulay mixing (mixing
during electronic SC iterations).
"""
},
'npweps': {
'definition': "Number of PlaneWaves for EPSilon (the dielectric matrix) ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>  """,
'vartype': "integer parameter",
'default': "0 ",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=3 or 4, that is, screening or sigma calculations.
<p>
<b>npweps</b> determines the size of the planewave set used to represent the independent-particle
susceptibility $\chi^{(0)}_{KS}$, the dielectric matrix $\epsilon$ and its inverse.
<br>See <a href="vargw.html#ecuteps">ecuteps</a> (preferred over <b>npweps</b>) for more information.
"""
},
'npwkss': {
'definition': "Number of PlaneWaves in the KSS file ",
'section': "vargw",
'category': "",
'vartype': "integer parameter  ",
'default': "0",
'text': """<p>
This input variable is used for the preparation of a GW calculation:
the GS run (where <a href="vargs.html#optdriver">optdriver</a>=1 and <b>nbandkss</b>/=0) should be followed
with a run where <a href="vargs.html#optdriver">optdriver</a>=3.
Also, if <a href="vargw.html#nbandkss">nbandkss</a>=0, no use of <b>npwkss</b>.
<p>
<b>npwkss</b> defines the number of planewave components of the Kohn-Sham states
to build the Hamiltonian, in the routine outkss.F90, and so, the size of the matrix, the size of eigenvectors,
and the number of available states, to be stored in the abo_KSS file.
If it is set to 0, then, the planewave basis set defined by the usual Ground State input variable
<a href="varbas.html#ecut">ecut</a> is used to generate the superset of all planewaves used for all k-points.
Note that this (large) planewave basis is the same for all k-points.
<p>
Very important : for the time being, <a href="vardev.html#istwfk">istwfk</a> must be 1 for all the k-points.
"""
},
'npwsigx': {
'definition': "Number of PlaneWaves for SIGma eXchange ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>  """,
'vartype': "integer parameter",
'default': "0 ",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=4, that is, sigma calculations.
<p>
<b>npwsigx</b> determines the cut-off energy of the planewave set used to generate the
exchange part of the self-energy operator.
<br>See <a href="vargw.html#ecutsigx">ecutsigx</a> (preferred over <b>npwsigx</b>) for more information.
"""
},
'npwwfn': {
'definition': "Number of PlaneWaves for WaveFunctioNs ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>  """,
'vartype': "integer parameter",
'default': "0 ",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=3 or 4, that is, screeening or sigma calculations.
<p>
<b>npwwfn</b> determines the size of the planewave set used to represent the wavefunctions
in the formula that generates the independent-particle susceptibility $\chi^{(0)}_{KS}$.
<br>
See <a href="vargw.html#ecutwfn">ecutwfn</a> (preferred over <b>nshwfn</b>) for more information.
"""
},
'nqpt': {
'definition': "Number of Q - POINTs ",
'section': "vargs",
'category': " ",
'vartype': "integer parameter ",
'default': "0.",
'text': """Determines whether one q point
must be read (See the variable <a href="varint.html#qptn">qptn</a>).
<br>Can be either 0 or 1.
<br>If 1 and used in ground-state calculation,
a global shift of all the k-points is applied, to give
calculation at k+q.
In this case, the output wavefunction will be appended
by _WFQ instead of _WFK (see the <a href="../users/abinit_help.html">section 4</a> of abinit_help)
Also, if 1 and a RF calculation is done, defines the
wavevector of the perturbation."""
},
'nqptdm': {
'definition': "Number of Q-PoinTs for the Dielectric Matrix",
'section': "vargw",
'category': " ",
'vartype': "integer parameter",
'default': "0",
'text': """<p>
Used only in the screening part, that is for <a href="vargs.html#optdriver">optdriver</a>=3.
<p>
If <b>nqptdm</b> is equal to 0, the set of q-points for computing the dielectric matrix is
determined automatically considering all the possible differences between the k-points contained
in the _KSS file. When <b>nqptdm</b> is non-zero, the list of q points is read from <a href="vargw.html#qptdm">qptdm</a>.
This allows one to split the big calculation of all the dielectric matrices into smaller
calculations that can be performed independently. The _SCR files generated in different runs
can be merged thanks to the <B>Mrgscr</B> utility.
"""
},
'nscforder': {
'definition': "SCaling Function ORDER",
'section': "vardev",
'category': " ",
'vartype': "",
'default': "16",
'text': """This variable controls the order of used scaling functions when the Hartree potential is computed using the Poisson solver (see <a href="vargs.html#icoulomb">icoulomb</a> imput variable). This variable is of seldom use since the default value is large enough. Nonetheless, possible values are 8, 14, 16, 20, 24, 30, 40, 50, 60, 100. Values greater than 20 are included in ABINIT for test purposes only.
"""
},
'nsheps': {
'definition': "Number of SHells for EPSilon (the dielectric matrix) ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>  """,
'vartype': "integer parameter",
'default': "0 ",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=3, that is, screening calculations.
<p>
<b>nsheps</b> determines the size of the planewave set used to represent the independent-particle
susceptibility $\chi^{(0)}_{KS}$, the dielectric matrix $\epsilon$ and its inverse.
<br>
See <a href="vargw.html#ecuteps">ecuteps</a> (preferred over <b>nsheps</b>) for more information.
"""
},
'nshiftk': {
'definition': "Number of SHIFTs for K point grids ",
'section': "varbas",
'category': " ",
'vartype': "integer parameter ",
'default': "1.",
'text': """This parameter
gives the number of shifted grids
to be used concurrently to generate the full grid of k points.
It can be used with primitive grids defined either from
<a href="varbas.html#ngkpt">ngkpt</a>
or
<a href="vargs.html#kptrlatt">kptrlatt</a>.
The maximum allowed value of <b>nshiftk</b> is 8.
The values of the shifts are given by <a href="varbas.html#shiftk">shiftk</a>.
"""
},
'nshiftq': {
'definition': "Number of SHIFTs for Q point grids ",
'section': "vargs",
'category': """<a href="../users/abinit_help.html#not_internal">NOT INTERNAL</a> """,
'vartype': "integer parameter ",
'default': "1.",
'text': """This parameter
gives the number of shifted grids
to be used concurrently to generate the full grid of q points.
It can be used with primitive grids defined either from
<a href="vargs.html#ngqpt">ngqpt</a>
or
<a href="vargs.html#qptrlatt">qptrlatt</a>.
The maximum allowed value of <b>nshiftq</b> is 8.
The values of the shifts are given by <a href="vargs.html#shiftq">shiftq</a>.
"""
},
'nshsigx': {
'definition': "Number of SHells for SIGma eXchange ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>  """,
'vartype': "integer parameter",
'default': "0 ",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=4, that is, sigma calculations.
<p>
<b>nshsigx</b> determines the cut-off energy of the planewave set used to generate the
exchange part of the self-energy operator.
<br>
See <a href="vargw.html#ecutsigx">ecutsigx</a> (preferred over <b>nshsigx</b>) for more information.
"""
},
'nshwfn': {
'definition': "Number of SHells for WaveFunctioNs ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>  """,
'vartype': "integer parameter",
'default': "0 ",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=3 or 4, that is, screening or sigma calculations.
<p>
<b>nshwfn</b> determines the number of shells of the planewave set used to represent the wavefunctions
in the formula that generates the independent-particle susceptibility $\chi^{(0)}_{KS}$.
<br>
See <a href="vargw.html#ecutwfn">ecutwfn</a> (preferred over <b>nshwfn</b>) for more information.
"""
},
'nspden': {
'definition': "Number of SPin-DENsity components ",
'section': "vargs",
'category': "",
'vartype': "integer parameter ",
'default': """the value of <a href="varbas.html#nsppol">nsppol</a>.""",
'text': """If <b>nspden</b>=1, no spin-magnetization : the density matrix is
diagonal, with same values spin-up and spin-down
(compatible with <a href="varbas.html#nsppol">nsppol</a>=1 only,
for both <a href="vargs.html#nspinor">nspinor</a>=1 or 2)
<p>
If <b>nspden</b>=2, scalar magnetization (the axis is arbitrarily
fixed in the z direction) : the density matrix is
diagonal, with different values for spin-up and spin-down
(compatible with <a href="vargs.html#nspinor">nspinor</a>=1,
either with <a href="varbas.html#nsppol">nsppol</a>=2 -general
collinear magnetization- or
<a href="varbas.html#nsppol">nsppol</a>=1 -antiferromagnetism)
<p>
If <b>nspden</b>=4, vector magnetization : the density matrix is full,
with allowed x, y and z magnetization
(useful only with <a href="vargs.html#nspinor">nspinor</a>=2 and
<a href="varbas.html#nsppol">nsppol</a>=1, either
because there is spin-orbit without time-reversal
symmetry - and thus spontaneous magnetization, or
with spin-orbit, if one allows for spontaneous
non-collinear magnetism). Not yet available for response functions. Also note
that, with <b>nspden</b>=4, time-reversal symmetry is not taken into account
(at present ; this has to be checked) and thus <a href="varbas.html#kptopt">kptopt</a>
has to be different from 1 or 2.

<p>
The default (<b>nspden</b>=<a href="varbas.html#nsppol">nsppol</a>)
does not suit the case of vector magnetization.
"""
},
'nspinor': {
'definition': "Number of SPINORial components of the wavefunctions ",
'section': "vargs",
'category': "",
'vartype': "integer parameter ",
'default': """1 (2 if <a href="varpaw.html#pawspnorb">pawspnorb</a>=1).""",
'text': """If <b>nspinor</b>=1, usual case : scalar wavefunction
(compatible with (<a href="varbas.html#nsppol">nsppol</a>=1,
<a href="vargs.html#nspden">nspden</a>=1) as well
as (<a href="varbas.html#nsppol">nsppol</a>=2, <a href="vargs.html#nspden">nspden</a>=2) )
<p>
If <b>nspinor</b>=2, the wavefunction is a spinor
(compatible with <a href="varbas.html#nsppol">nsppol</a>=1, with
<a href="vargs.html#nspden">nspden</a>=1 or 4,
but not with <a href="varbas.html#nsppol">nsppol</a>=2)
<p>When <b>nspinor</b> is 2, the values of <a href="vardev.html#istwfk">istwfk</a>
are automatically set to 1. Also, the number of bands, for each k-point,
should be even.
"""
},
'nsppol': {
'definition': "Number of SPin POLarization ",
'section': "varbas",
'category': " ",
'vartype': "integer parameter ",
'default': "1.",
'text': """Give the number of INDEPENDENT
spin polarisations. Can take the values
1 or 2.
<p>If <b>nsppol</b>=1, one has an unpolarized calculation
(<a href="vargs.html#nspinor">nspinor</a>=1,
<a href="vargs.html#nspden">nspden</a>=1) or
an antiferromagnetic system
(<a href="vargs.html#nspinor">nspinor</a>=1,
<a href="vargs.html#nspden">nspden</a>=2), or
a calculation in which spin up and spin down cannot be disantengled
(<a href="vargs.html#nspinor">nspinor</a>=2), that is, either
non-collinear magnetism or presence of spin-orbit coupling,
for which one needs spinor wavefunctions.
<p>If <b>nsppol</b>=2, one has a spin-polarized (collinear) calculation
with separate and different wavefunctions for up and
down spin electrons for each band and k point.
Compatible only with <a href="vargs.html#nspinor">nspinor</a>=1,
<a href="vargs.html#nspden">nspden</a>=2.
<p>In the present status of development,
with <b>nsppol</b>=1,
all values of <a href="varbas.html#ixc">ixc</a> are allowed, while
with <b>nsppol</b>=2,
some values of <a href="varbas.html#ixc">ixc</a> might not be allowed (e.g. 2, 3, 4, 5, 6, 20, 21, 22 are not allowed).
<p>See also the input variable <a href="vargs.html#nspden">nspden</a>
for the components of the density matrix with respect to
the spin-polarization.
"""
},
'nstep': {
'definition': "Number of (non-)self-consistent field STEPS ",
'section': "varbas",
'category': " ",
'vartype': "integer parameter ",
'default': "30. (was 1 before v5.3)",
'text': """Gives the maximum number of cycles (or "iterations") in a SCF or non-SCF run.
<br>Full convergence from random numbers if usually achieved in
12-20 SCF iterations. Each can take from minutes to hours.
In certain difficult cases, usually related to a small or
zero bandgap or magnetism, convergence performance may be much worse.
When the convergence tolerance <a href="varbas.html#tolwfr">tolwfr</a> on the wavefunctions
is satisfied, iterations will stop, so for well converged
calculations you should set <b>nstep</b> to a value larger than
you think will be needed for full convergence, e.g.
if 20 steps usually converges the system, set <b>nstep</b> to 30.
<br>
For non-self-consistent runs (<a name="iscf">iscf</a> &lt; 0) nstep governs
the number of cycles of convergence for the wavefunctions for a fixed density
and Hamiltonian.
<p>NOTE that a choice of <b>nstep</b>=0 is permitted; this will
either read wavefunctions from disk (with <a href="varfil.html#irdwfk">irdwfk</a>=1
or  <a href="varfil.html#irdwfq">irdwfq</a>=1,
or non-zero <a href="varfil.html#getwfk">getwfk</a>
or <a href="varfil.html#getwfq">getwfq</a> in the case
of multi-dataset) and
compute the density, the total energy and stop, or else
(with all of the above vanishing) will initialize
randomly the wavefunctions and
compute the resulting density and total energy.
This is provided for testing purposes.
<br>Also NOTE that <b>nstep</b>=0
with <a href="varfil.html#irdwfk">irdwfk</a>=1 will exactly give the same result as
the previous run only if the latter is done with <a href="varbas.html#iscf">iscf</a><10
(potential mixing).
<br>One can output the density by using <a href="varfil.html#prtden">prtden</a>.
<br>The forces and stress tensor are computed with <b>nstep</b>=0."""
},
'nsym': {
'definition': "Number of SYMmetry operations ",
'section': "varbas",
'category': """<a href="../users/abinit_help.html#symmetry_finder">SYMMETRY FINDER</a> """,
'vartype': "integer parameter ",
'default': "0.",
'text': """Gives number of space group symmetries
to be applied in this problem.  Symmetries will be input in
array "<a href="varbas.html#symrel">symrel</a>" and (nonsymmorphic) translations vectors
will be input
in array "<a href="varbas.html#tnons">tnons</a>".  If there is no symmetry in the problem
then set <b>nsym</b> to 1, because the identity is still a symmetry.

<br>In case of a RF calculation, the code is able to use
the symmetries of the system to decrease the number of
perturbations to be calculated, and to decrease of the
number of special k points to be used for the sampling of
the Brillouin zone.
After the response to the perturbations have been calculated,
the symmetries are used to generate as many as
possible elements of the 2DTE from those already
computed.

<p>SYMMETRY FINDER mode (Default mode).
If <b>nsym</b> is 0, all the atomic coordinates must be
explicitely given (one cannot use the geometry builder
and the symmetrizer): the code will then find automatically
the symmetry operations that leave the lattice and each
atomic sublattice invariant. It also checks whether the
cell is primitive (see <a href="vargs.html#chkprim">chkprim</a>).
<br>Note that the tolerance on symmetric atomic positions and
lattice is rather stringent :
for a symmetry operation to be admitted,
the lattice and atomic positions must map on themselves
within 1.0e-8 .

<p>The user is allowed to set up systems with non-primitive unit cells (i.e.
conventional FCC or BCC cells, or supercells without any distortion).
In this case, pure translations will be identified as symmetries
of the system by the symmetry finder.
Then, the combined "pure translation + usual rotation and inversion" symmetry
operations can be very numerous. For example, a conventional FCC cell
has 192 symmetry operations, instead of the 48 ones of the primitive cell.
A maximum limit of 384 symmetry operations is hard-coded. This
corresponds to the maximum number of symmetry operations of a 2x2x2
undistorted supercell. Going beyond
that number will make the code stop very rapidly. If you want
nevertheless, for testing purposes, to treat a larger number of symmetries,
change the initialization of "msym" in the abinit.F90 main routine,
then recompile the code.

<p>For GW calculation, the user might want to select only the symmetry operations whose
non-symmorphic translation vector <a href="varbas.html#tnons">tnons</a>
is zero. This can be done with the help of the input variable
<a href="vargw.html#symmorphi">symmorphi</a>
"""
},
'ntime': {
'definition': "Number of TIME steps ",
'section': "varrlx",
'category': "",
'vartype': "integer parameter ",
'default': "0.",
'text': """Gives the number of molecular dynamics time steps or
Broyden
structural optimization steps to be done if <a
href="varrlx.html#ionmov">ionmov</a>is non-zero.<br>
Note that at the present
the option <a href="varrlx.html#ionmov">ionmov</a>=1 is initialized
with four
Runge-Kutta steps which costs some overhead in the startup.
By contrast, the initialisation of
other <a href="varrlx.html#ionmov">ionmov</a> values is only one
SCF call. <br>
<b>ntime</b> is ignored if <a href="varrlx.html#ionmov">ionmov</a>=0.
"""
},
'ntimimage': {
'definition': "Number of TIME steps for IMAGE propagation ",
'section': "varrlx",
'category': "",
'vartype': "integer parameter ",
'default': "1.",
'text': """Gives the maximal number of molecular dynamics time steps or
structural optimization steps to be done
for the set of images, referred to as 'image-timesteps'. At each image-timestep,
all the images are propagated simulaneously, each according to the algorithm
determined by <a href="varrlx.html#imgmov">imgmov</a> and the usual accompanying
input variables, and then the next positions and velocities for each image
are determined from the set of results obtained for all images.
"""
},
'ntypalch': {
'definition': """Number of TYPe of atoms that are "ALCHemical" """,
'section': "vargs",
'category': "",
'vartype': "integer parameter ",
'default': "0",
'text': """<p> Used for the generation of alchemical pseudopotentials :
when <b>ntypalch</b> is non-zero, alchemical mixing
will be used.
<p> Among the <a href="varbas.html#ntypat">ntypat</a> types of atoms, the
last <b>ntypalch</b> will be "alchemical" pseudoatoms, while only
the first <b>ntyppure</b> will be uniquely associated with a pseudopotential
(the <b>ntyppure</b> first of these, actually). The
<a href="vargs.html#ntypalch">ntypalch</a> types of alchemical
pseudoatoms are to be made
from the remaining <a href="vargs.html#npspalch">npspalch</a> pseudopotentials.
<p> In this case,
the input variables
<a href="vargs.html#algalch">algalch</a>,<a href="vargs.html#mixalch">mixalch</a>
are active, and generate alchemical potentials from the available
pseudopotentials.  See these input variables, especially
<a href="vargs.html#mixalch">mixalch</a>, to understand how
to use alchemical potentials in ABINIT.
"""
},
'ntypat': {
'definition': "Number of TYPEs of atoms ",
'section': "varrlx",
'category': """<a href="../users/abinit_help.html#no_multi">NO MULTI</a> """,
'vartype': "integer parameter ",
'default': "1.",
'text': """Gives the number of types of atoms. E.g. for
a homopolar system (e.g. pure Si) <b>ntypat</b> is 1.
<br>
The code tries to read the same number pseudopotential files.
<br>
The first pseudopotential is assigned type number 1, and so
on ...
"""
},
'ntyppure': {
'definition': """Number of TYPe of atoms that are "PURe" """,
'section': "vargs",
'category': "Inner ",
'vartype': "integer parameter, non-negative ",
'default': "",
'text': """Gives the number of type of atoms that are "pure" when alchemical mixing is used (<a href="vargs.html#ntypalch">ntypalch</a> /= 0) :
<p> <b>ntyppure</b>=<a href="varbas.html#ntypat">ntypat</a>-<a href="vargs.html#ntypalch">ntypalch</a>
"""
},
'nwfshist': {
'definition': "Number of WaveFunctionS HISTory ",
'section': "vargs",
'category': "",
'vartype': "integer parameter, non-negative ",
'default': "0",
'text': """<p>In the wavelet basis set, the ground state is found by direct
minimisation. The algorithm used can be either the steepest
descent or the DIIS (Direct Inversion of Iteration
Space). When <b>nwfshist</b> = 0, the steepest descent is used
(<i>i.e.</i> there is no history storage of the previous
iterations). If <b>nwfshist</b> is strictly positive, a DIIS
is used. A typical value is 6. Using a DIIS increases the
memory required by the program since N previous wavefunctions
are stored during the electronic minimisation.</p>
"""
},
'objaat': {
'definition': "OBJect A : list of AToms, OBJect B : list of AToms ",
'section': "vargeo",
'category': """<a href="../users/abinit_help.html#geometry_builder">GEOMETRY BUILDER</a>, <a href="../users/abinit_help.html#not_internal">NOT INTERNAL</a> """,
'vartype': """integer arrays objaat(<a href="vargeo.html#objan">objan</a>) and objbat(<a href="vargeo.html#objan">objbn</a>) """,
'default': "",
'text': """Gives the list of atoms in either
object a or object b. This list is specified
by giving, for each atom, its index in the list
of coordinates (<a href="varbas.html#xred">xred</a>, <a href="varbas.html#xangst">xangst</a>
or <a href="varbas.html#xcart">xcart</a>), that also
corresponds to a type of atom (given by the array type).
These objects can be thought as molecules, or groups
of atoms, or parts of molecules, to be repeated,
rotated and translated to generate the full set
of atoms. <br>Look at <a href="vargeo.html#objarf">objarf</a> and <a href="vargeo.html#objarf">objbrf</a>
for further explanations.
<b>objaat</b> MUST be provided if <a href="vargeo.html#nobj">nobj</a>==1.
<b>objaat</b> and <b>objbat</b> MUST be provided if <a href="vargeo.html#nobj">nobj</a>==2.
<br>Not present in the dtset array (no internal).
"""
},
'objaax': {
'definition': "OBJect A : AXis, OBJect B : AXis",
'section': "vargeo",
'category': """<a href="../users/abinit_help.html#geometry_builder">GEOMETRY BUILDER</a>, <a href="../users/abinit_help.html#not_internal">NOT INTERNAL</a>, <a href="../users/abinit_help.html#length">LENGTH</a> """,
'vartype': "real arrays objaax(6) and objbax(6) ",
'default': "",
'text': """Gives, for each object, the cartesian coordinates of
two points (first point :  <b>objaax</b>(1:3) or <b>objbax</b>(1:3),
second point : <b>objaax</b>(4:6) or <b>objbax</b>(4:6) ).
<br>By default, given in Bohr atomic units
(1 Bohr=0.5291772108 Angstroms), although Angstrom can be specified,
if preferred, since these variables have the
'<a href="../users/abinit_help.html#dimensions">LENGTH</a>' characteristics.
<br>
The two points define an axis of rotation
of the corresponding object. <br>Note that the rotation of
the object is done BEFORE the object is translated.
<br>The sign of the rotation angle is positive if the
object is to be rotated clockwise when looking to it
along the axis, from point 1 (coordinates 1:3)
toward point 2 (coordinates 4:6).
<br><a href="vargeo.html#objaat">objaat</a> MUST be provided if <a href="vargeo.html#nobj">
nobj</a>==1 and one component
of <a href="vargeo.html#objaro">objaro</a> does not vanish.
<br><a href="vargeo.html#objaat">objaat</a> and <a href="vargeo.html#objaat">objbat</a>
MUST be provided if <a href="vargeo.html#nobj">nobj</a>==2 and one
component of <a href="vargeo.html#objaro">objbro</a> does not vanish.
<br>Not present in the dtset array (no internal).
"""
},
'objan': {
'definition': "OBJect A : Number of atoms, OBJect B : Number of atoms ",
'section': "vargeo",
'category': """<a href="../users/abinit_help.html#geometry_builder">GEOMETRY BUILDER</a>, <a href="../users/abinit_help.html#not_internal">NOT INTERNAL</a>""",
'vartype': "integer parameters ",
'default': "",
'text': """Gives the number of atoms in either
object a or object b. The list of atoms is given by the
variables <a href="vargeo.html#objaat">objaat</a> and <a href="vargeo.html#objaat">objbat</a>.
<br><b>objan</b> MUST be provided if <a href="vargeo.html#nobj">nobj</a>==1.
<br><b>objan</b> and <b>objbn</b> MUST be provided if <a href="vargeo.html#nobj">nobj</a>==2.
<br>Not present in the dtset array (no internal).
"""
},
'objarf': {
'definition': "OBJect A : Repetition Factors, OBJect B : Repetition Factors ",
'section': "vargeo",
'category': """<a href="../users/abinit_help.html#geometry_builder">GEOMETRY BUILDER</a>, <a href="../users/abinit_help.html#not_internal">NOT INTERNAL</a>""",
'vartype': "integer arrays objarf(3) and objbrf(3) ",
'default': "1 1 1 .",
'text': """Gives three repetition factors of the objects a or b.
<br>This gives the opportunity to generate a three-dimensional
set of repeated objects, although a simple one-dimensional
repetition will be easily obtained through the specification of<br>
nrep 1 1  &lt;r&gt;
where nrep is the 1D repetition factor.
<br>The initial rotation and translation of the object,
as well as the increment of rotation or translation
from one object to the next are specified by
the variables <a href="vargeo.html#objaro">objaro</a> and <a href="vargeo.html#objatr">objatr</a>, for object a,
and <a href="vargeo.html#objaro">objbro</a> and <a href="vargeo.html#objatr">objbtr</a>, for object b.
<br>Note that the geometry builder will generate the full set
of atoms from the primitive set of atoms using the
following order : it will process each atom in the
primitive list one by one, determine whether it belongs
to either object a or object b, and then repeat it
taking into account the proper rotation and translation,
with the fastest varying repetition factor being the first,
then the second, then the third.
<br>In the final list of atoms, one will first
find the atoms generated from atom 1 in the primitive list,
then those generated from atom 2 in the primitive list, and
so on.
<br>If the geometry builder is only used to rotate
or translate an object, without repeating it,
simply use 1 1 1, which is also the Default value.
<br>Not present in the dtset array (no internal).
"""
},
'objaro': {
'definition': "OBJect A : ROtations, OBJect B : ROtations ",
'section': "vargeo",
'category': """<a href="../users/abinit_help.html#geometry_builder">GEOMETRY BUILDER</a>, <a href="../users/abinit_help.html#not_internal">NOT INTERNAL</a>""",
'vartype': "real arrays objaro(4) and objbro(4) ",
'default': "4*0.0d0 (no rotation).",
'text': """Give, for each object, the angles of rotation in degrees
to be applied to the corresponding object.
<br>The rotation is applied before the translation,
and the axis is defined by the variables <a href="vargeo.html#objaax">objaax</a> and
<a href="vargeo.html#objaax">objbax</a>.
See the latter variables for the definition of the sign
of the rotation.
<br>The first component <b>objaro</b>(1) and <b>objbro</b>(1) gives
the angle of rotation to be applied to the first
instance of the object. The second, third or fourth
component (resp.) gives the increment of rotation angle
from one instance to the next instance, defined
by the first, second or third repetition factor (resp.) .
This allows to generate 3D arrays of molecules with
different rotation angles.
<br>Not present in the dtset array (no internal).
"""
},
'objatr': {
'definition': "OBJect A : TRanslations, OBJect B : TRanslations ",
'section': "vargeo",
'category': """<a href="../users/abinit_help.html#geometry_builder">GEOMETRY BUILDER</a>, <a href="../users/abinit_help.html#not_internal">NOT INTERNAL</a>, <a href="../users/abinit_help.html#length">LENGTH</a> """,
'vartype': "real arrays objatr(3,4) and objbtr(3,4)",
'default': "12*0.0d0 (no translation).",
'text': """Give, for each object, the vectors of translations,
in cartesian coordinates,
to be applied to the corresponding object.
By default, given in Bohr atomic units
(1 Bohr=0.5291772108 Angstroms), although Angstrom can be specified,
if preferred, since these variables have the
'<a href="../users/abinit_help.html#dimensions">LENGTH</a>' characteristics.
<br>The translation is applied after the rotation.
<br>The first vector <b>objatr</b>(3,1) and <a href="vargeo.html#objaro">objbro</a>(3,1) gives
the translation to be applied to the first
instance of the object. The second, third or fourth
component (resp.) gives the increment of translation
from one instance to the next instance, defined
by the first, second or third repetition factor (resp.) .
This allows to generate 3D arrays of molecules.
<br>In general, when the objects are repeated, a translation
vector must be given, since otherwise, the repeated objects
pack in the same region of space. As an exception, one can
have a set of molecules regularly spaced on a circle, in
which case, only rotations are needed.
<br>Not present in the dtset array (no internal).
"""
},
'occ': {
'definition': "OCCupation numbers ",
'section': "vargs",
'category': """<a href="../users/abinit_help.html#evolving">EVOLVING</a> """,
'vartype': """real array occ(<a href="varbas.html#nband">nband</a>)  """,
'default': "0's.",
'text': """Gives occupation numbers for all
bands in the problem. Needed if <a href="varbas.html#occopt">occopt</a>==0 or
<a href="varbas.html#occopt">occopt</a>==2.
Ignored otherwise. Also ignored when <a href="varbas.html#iscf">iscf</a>=-2.
<br>Typical band occupancy is either
2 or 0, but can be 1 for half-occupied band or other
choices in special circumstances.  <br>If <a href="varbas.html#occopt">occopt</a> is not 2,
then the occupancies must be the same for each k point.
<br>If <a href="varbas.html#occopt">occopt</a>=2, then the band occupancies must be
provided explicitly for each band, EACH k POINT,
and EACH SPIN-POLARIZATION, in an
array which runs over all bands, k points,
and spin-polarizations.
<br>The order of entries in the array would correspond to
all bands at the first k point (spin up), then all bands at the
second k point (spin up), etc, then all k-points spin down.
<br>The total number of array elements
which must be provided is <br>
( <a href="varbas.html#nband">nband</a>(1)+<a href="varbas.html#nband">nband</a>(2)+...+
<a href="varbas.html#nband">nband</a>(<a href="varbas.html#nkpt">nkpt</a>) ) *
<a href="varbas.html#nsppol">nsppol</a> .
<br>The occupation numbers evolve only for metallic occupations,
that is, <a href="varbas.html#occopt">occopt</a> &ge; 3 ."""
},
'occopt': {
'definition': "OCCupation OPTion ",
'section': "varbas",
'category': " ",
'vartype': "integer option parameter  ",
'default': "<b>occopt</b>=1.",
'text': """Controls how input
parameters <a href="varbas.html#nband">nband</a>, <a href="vargs.html#occ">occ</a>,
and <a href="varbas.html#wtk">wtk</a> are handled.
<ul>
<li>              <b>occopt</b>=0:
<br>All k points have the same number of bands
and the same occupancies of bands.  <a href="varbas.html#nband">nband</a> is given as a
single number, and <a href="vargs.html#occ">occ</a>(<a href="varbas.html#nband">nband</a>)
is an array of <a href="varbas.html#nband">nband</a>
elements, read in by the code.  <br>
The k point weights in array <a href="varbas.html#wtk">wtk</a>(<a href="varbas.html#nkpt">nkpt</a>) are
automatically normalized by the code to add to 1.</li>
<li>              <b>occopt</b>=1:
<br>Same as <b>occopt</b>=0, except that the array <a href="vargs.html#occ">occ</a> is
automatically generated by the code, to give a semiconductor.
<br>An error occurs when filling cannot be done with
occupation numbers equal to 2 or 0 in each k-point (non-spin-polarized case),
or with occupation numbers equal to 1 or 0 in each k-point (spin-polarized case).</li>
<li>             <b>occopt</b>=2:
<br>k points may optionally have different numbers of
bands and different occupancies.  <a href="varbas.html#nband">nband</a>(<a href="varbas.html#nkpt">
nkpt</a>*<a href="varbas.html#nsppol">nsppol</a>) is given
explicitly as an array of <a href="varbas.html#nkpt">nkpt</a>*<a href="varbas.html#nsppol">nsppol</a> elements.
<a href="vargs.html#occ">occ</a>() is given explicitly for all bands at each k point,
and eventually for each spin --
the total number of elements is the sum of <a href="varbas.html#nband">nband</a>(ikpt)
over all k points and spins. The k point weights <a href="varbas.html#wtk">wtk</a>
(<a href="varbas.html#nkpt">nkpt</a>) are
NOT automatically normalized under this option.</li>
<li>             <b>occopt</b>=3, 4, 5, 6 and 7
<br>Metallic occupation of levels, using different occupation
schemes (see below). The corresponding thermal
broadening, or cold smearing, is defined by
the input variable <a href="vargs.html#tsmear">tsmear</a> (see below : the variable
xx is the energy in Ha, divided by <a href="vargs.html#tsmear">tsmear</a>)
<br>Like for <b>occopt</b>=1, the variable <a href="vargs.html#occ">occ</a> is not read
<br>All k points have the same number of bands,
<a href="varbas.html#nband">nband</a> is given as a single number, read by the code.
<br>The k point weights in array <a href="varbas.html#wtk">wtk</a>(<a href="varbas.html#nkpt">nkpt</a>) are
automatically normalized by the code to add to 1.</li>
<ul>
<li> <b>occopt</b>=3:
<br>Fermi-Dirac smearing (finite-temperature metal)
Smeared delta function : 0.25d0/(cosh(xx/2.0d0)**2)</li>
<li> <b>occopt</b>=4:
<br>"Cold smearing" of N. Marzari (see his thesis work),
with a=-.5634 (minimization of the bump)
<br>Smeared delta function :
<br>  exp(-xx<sup>2</sup>)/sqrt(pi) * (1.5d0+xx*(-a*1.5d0+xx*(-1.0d0+a*xx)))</li>
<li> <b>occopt</b>=5:
<br>"Cold smearing" of N. Marzari (see his thesis work),
with a=-.8165 (monotonic function in the tail)
<br>Same smeared delta function as <b>occopt</b>=4, with different a.</li>
<li> <b>occopt</b>=6:
<br>Smearing of Methfessel and Paxton (PRB40,3616(1989))
with Hermite polynomial of degree 2, corresponding
to "Cold smearing" of N. Marzari with a=0
(so, same smeared delta function as <b>occopt</b>=4, with different a).</li>
<li> <b>occopt</b>=7:
<br>Gaussian smearing, corresponding to the 0 order
Hermite polynomial of Methfessel and Paxton.
<br>Smeared delta function : 1.0d0*exp(-xx**2)/sqrt(pi)</li>
<li> <b>occopt</b>=8:
<br>Uniform smearing (the delta function is replaced by a constant function of value one
over ]-1/2,1/2[ (with one-half value at the boundaries). Used for testing purposes only.
</li>
</ul>
</ul>

WARNING : one can use metallic occupation of levels in the
case of a molecule, in order to avoid any problem with
degenerate levels. However, it is advised NOT to use
<b>occopt</b>=6 (and to a lesser extent <b>occopt</b>=4 and 5),
since the associated number of electron
versus the Fermi energy is NOT guaranteed to be
a monotonic function. For true metals, AND a sufficiently
dense sampling of the Brillouin zone, this should not happen,
but be cautious ! As an indication of this problem,
a small variation of input parameters might lead to
a jump of total energy, because there might be two or even
three possible values of the Fermi energy, and the
bisection algorithm find one or the other. """
},
'omegasimax': {
'definition': "OMEGA to evaluate Sigma along the Imaginary axis D: MAXimal value ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>, <a href="../users/abinit_help.html#energy">ENERGY</a>""",
'vartype': "real ",
'default': "50 eV ",
'text': """<p>
Only relevant for sigma calculations in which the self-energy along the real axis is obtained
by performing the analytic continuation from the imaginary axis to the full complex plane via the Pade approximant
(<a href="vargs.html#optdriver">optdriver</a>=4 and <a href="vargw.html#gwcalctyp">gwcalctyp</a>=x1).
<p>
<b>omegasimax</b> defines the maximum frequency along the imaginary the axis.
In conjunction with <a href="vargw.html#nomegasi">nomegasi</a>,
<b>omegasimax</b> uniquely defines the linear mesh employed to sample the self-energy along the imaginary axis.

"""
},
'omegasrdmax': {
'definition': "OMEGA to evaluate the Sigma Real axis Derivative : MAXimal value ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>, <a href="../users/abinit_help.html#energy">ENERGY</a>""",
'vartype': "real ",
'default': "1.0 eV ",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=4, that is, sigma calculations.
<p>
The maximum distance from the KS energy where to evaluate Sigma. Sigma is evaluated at
[KS_energy - maxomegasrd, KS_energy + maxomegasrd] sampled <a href="vargw.html#nomegasrd">nomegasrd</a> times.
"""
},
'optcell': {
'definition': "OPTimize the CELL shape and dimensions",
'section': "varrlx",
'category': "",
'vartype': "integer parameter ",
'default': "<b>optcell</b>=0",
'text': """Presently (v3.1), one cannot restart (<a
href="varrlx.html#restartxf">restartxf</a>)
a calculation with a non-zero <b>optcell</b> value from the
(x,f) history of another run with a different non-zero <b>optcell</b>
value. There
are still a few problems at that level.
"""
},
'optdriver': {
'definition': "OPTions for the DRIVER ",
'section': "vargs",
'category': " ",
'vartype': "integer parameter  ",
'default': "<b>optdriver</b>=0 ",
'text': """For each dataset, choose
the task to be done, at the level of the "driver" routine.
<p>The choice is among :
<br><b>optdriver</b>=0 : ground-state calculation (GS), routine "gstate"
<br><b>optdriver</b>=1 : response-function calculation (RF), routine "respfn"
<br><b>optdriver</b>=2 : susceptibility calculation (SUS), routine "suscep"
<br><b>optdriver</b>=3 : susceptibility and dielectric matrix calculation (SCR), routine "screening"
<br>
(see the input variables <a href="vargw.html#ecutwfn">ecutwfn</a>,
<a href="vargw.html#ecuteps">ecuteps</a>,
<a href="vargw.html#ppmfrq">ppmfrq</a>,
<a href="varfil.html#getkss">getkss</a>,
as well as <a href="vargw.html#nbandkss">nbandkss</a> and <a href="varbas.html#nband">nband</a>)
<br><b>optdriver</b>=4 : self-energy calculation (SIG), routine "sigma"
<br><b>optdriver</b>=5 : non-linear response functions, using the 2n+1 theorem, routine "nonlinear"
<br><b>optdriver</b>=99 : Bethe-Salpeter calculation, routine "bethe_salpeter"
<p>
If one of <a href="varrf.html#rfphon">rfphon</a>, <a href="varrf.html#rfddk">rfddk</a>,
<a href="varrf.html#rfelfd">rfelfd</a>,
or <a href="varrf.html#rfstrs">rfstrs</a> is non-zero, while <b>optdriver</b>
is not defined in the input file, ABINIT will set <b>optdriver</b> to 1
automatically. These input variables (<a href="varrf.html#rfphon">rfphon</a>,
<a href="varrf.html#rfddk">rfddk</a>, <a href="varrf.html#rfelfd">rfelfd</a>,
and <a href="varrf.html#rfstrs">rfstrs</a>) must be
zero if <b>optdriver</b> is not set to 1.
"""
},
'optforces': {
'definition': "OPTions for the calculation of FORCES ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>  """,
'vartype': "integer parameter  ",
'default': """2. However, <b>optforces</b> is automatically set to 1 when <a href="varbas.html#toldff">toldff</a> or <a href="varbas.html#tolrff">tolrff</a> are non-zero.""",
'text': """Allows to choose options for the calculation of forces.
<ul>
<li><b>optforces</b>=0 : the forces are set to zero, and many steps of the
computation of forces are skipped </li>
<li><b>optforces</b>=1 : calculation of forces at each SCF iteration, allowing
to use forces as criterion to stop the SCF cycles
</li>
<li><b>optforces</b>=2 : calculation of forces at the end of the SCF iterations
(like the stresses)
</li>
</ul>
"""
},
'optfreqsus': {
'definition': "OPTion for the generation of FREQuency grids for the SUSceptibility",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>""",
'vartype': "integer parameter",
'default': "2",
'text': """<p>
Selects the type of frequency grid that will be used to compute ACFD energies,
as follows:
</p>
<ul>
<li>0: use preassigned mesh (see defs_suscep module)
<ul>
<li>nfreqsus= 2: pick-up 2 highest frequencies of H_2 mesh</li>
<li>nfreqsus= 8: pick-up 8 frequencies inside Be_2 mesh, depending on freq1</li>
<li>nfreqsus= 9: pick-up 9 frequencies inside H_2 mesh, depending on freq1</li>
<li>nfreqsus=11: pick-up 11 highest frequencies of Be_2 mesh</li>
<li>nfreqsus=16: use full He mesh</li>
<li>nfreqsus=18: use full H_2 mesh</li>
<li>nfreqsus=20: use full He mesh good up to 8 Ha</li>
<li>nfreqsus=24: use full Be_2 mesh</li>
</ul>
</li>
<li>1: create linear mesh and weights for quadrature by Taylor rule
<ul>
<li>freqsusin=starting frequency</li>
<li>freqsuslo=frequency increment</li>
</ul>
</li>
<li>2: create mesh and weights using Gauss-Legendre quadrature
<p>A first Gauss-Legendre mesh is built for interval [0,freqsuslo], then
a second one is obtained by transforming the first for the
[freqsuslo,+\infty[ interval. freqsusin may be use to compress or expand
the mesh on the second interval (a value of 1.0 is adequate for
most cases). For practical reasons, nfreqsus must be even.
</p>
</li>
</ul>

<p>
See also:
<a href="vardev.html#nfreqsus">nfreqsus</a>,
<a href="vardev.html#freqsuslo">freqsuslo</a>,
<a href="vardev.html#freqsusin">freqsusin</a>.
</p>

"""
},
'optnlxccc': {
'definition': "OPTion for the calculation of Non-Linear eXchange-Correlation Core Correction",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>  """,
'vartype': "integer parameter  ",
'default': "1.",
'text': """Allows to choose options for the calculation of non-linear XC correction.
At present, only relevant for the FHI type of pseudopotentials, with pspcod=6 .
<ul>
<li><b>optnlxccc</b>=1 : uses the old psp6cc.f routine, with inconsistent treatment of real-space derivatives of the core function (computed in this routine, while splined in the other parts of the code) </li>
<li><b>optnlxccc</b>=2 : consistent calculation derivatives, in the psp6cc_dhr.f routine from DHamann.</li>
</ul>
"""
},
'optstress': {
'definition': "OPTion for the computation of STRess",
'section': "vargs",
'category': " ",
'vartype': "integer ",
'default': "1",
'text': """If set to 1, the computation of stresses is done,
in the SCF case
(under the conditions <a href="varbas.html#iscf">iscf</a> > 0 , <a href="varfil.html#prtstm">prtstm</a>==0 ,
<a href="vargs.html#positron">positron</a>==0,
and either  <a href="varbas.html#nstep">nstep</a> >0 , or
<a href="varint.html#usepaw">usepaw</a>==0 or <a href="varfil.html#irdwfk">irdwfk</a>==1).
<br>
Otherwise, to save CPU time, if no optimization of the cell is required,
one can skip the computation of stresses. The CPU time saving might be interesting
for some PAW calculations.
"""
},
'ortalg': {
'definition': "ORThogonalisation ALGorithm ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>  """,
'vartype': "integer parameter  ",
'default': """2 when <a href="vardev.html#wfoptalg">wfoptalg</a> &lt; 10,<br> -2 when <a href="vardev.html#wfoptalg">wfoptalg</a> >=10.""",
'text': """Allows to choose the algorithm
for orthogonalisation.
<br>Positive or zero values make two projections per
line minimisation, one before the preconditioning, one
after. This is the clean application of the band-by-band
CG gradient for finding eigenfunctions.
<br>Negative values make only one projection per line mininisation.
<br>The orthogonalisation step is twice faster, but the
convergence is less good. This actually calls to
a better understanding of this effect.
<br><b>ortalg</b>=0, 1 or -1 is the conventional coding, actually
identical to the one in versions prior to 1.7
<br><b>ortalg</b>=2 or -2 try to make better use of existing registers
on the particular machine one is running.
<br>More demanding use of registers
is provided by <b>ortalg</b>=3 or -3, and so on.
<br>The maximal value is presently 4 and -4.
<br>Tests have shown that <b>ortalg</b>=2 or -2 is suitable for
use on the available platforms."""
},
'papiopt': {
'definition': "PAPI OPTion",
'section': "vardev",
'category': " ",
'vartype': "integer",
'default': "0 ",
'text': """<a href="http://icl.cs.utk.edu/papi/index.html">PAPI</a> aims to
provide the tool designer and application engineer with a
consistent interface and methodology for use of the
performance counter hardware found in most major
microprocessors. PAPI enables software engineers to see, in
near real time, the relation between software performance and
processor events.<br />

This option can be used only when ABINIT has been compiled with the
<code>--enable-papi</code> configure option.<br />

If <b>papiopt</b>=1, then PAPI counters are used instead of
the usual time() routine. All the timing output of ABINIT is
then done with PAPI values. The measurements are more accurate and
give also access to the flops of the calculation.
"""
},
'paral_atom': {
'definition': "activate PARALelization over (paw) ATOMic sites",
'section': "varpar",
'category': "",
'vartype': "integer ",
'default': "0.",
'text': """Relevant only for PAW calculations.<br>
This keyword controls the parallel distribution of memory over atomic sites. Calculations are
also distributed using the "kpt-band" communicator.<br>
Warning: use of <b>paral_atom</b> is highly experimental.<br>
Only compatible (for the moment) with ground-state calculations.
<br>
"""
},
'paral_kgb': {
'definition': "activate PARALelization over K-point, G-vectors and Bands",
'section': "varpar",
'category': "",
'vartype': "integer ",
'default': "0.",
'text': """<p>
<b>If paral_kgb is not explicitely put in the input file</b>,
ABINIT automatically detects if the job has been sent in sequential or in parallel. In this last case, it detects the number of processors on which the job has been sent and calculates values of
<a href="varpar.html#npkpt">npkpt</a>, <a href="varpar.html#npfft">npfft</a>,
<a href="varpar.html#npband">npband</a>, <a href="varpar.html#bandpp">bandpp</a> ,<a href="varpar.html#npimage">npimage</a> and <a href="varpar.html#npspinor">npspinor</a> that are compatible with the number of processors. It then set
paral_kgb to 0 or 1 (see hereunder)and launches the job.
</p>
<b>If paral_kgb=0</b>,
the parallelization over k-points only is activated. In this case,
<a href="varpar.html#npkpt">npkpt</a>, <a href="varpar.html#npfft">npfft</a> and
<a href="varpar.html#npband">npband</a> are ignored. Require compilation option --enable-mpi="yes".
</p>




<p>
<b>If paral_kgb=1</b>,
the parallelization over bands, FFTs, and k-point/spin-components is activated
(see <a href="varpar.html#npkpt">npkpt</a>, <a href="varpar.html#npfft">npfft</a> and
<a href="varpar.html#npband">npband</a>). With this parallelization, the work load is split over
three levels of parallelization. The different communications almost occur
along one dimension only. Require compilation option --enable-mpi="yes".
</p>



HOWTO fix the number of processors along one level of parallelisation:
<br>
At first, try to parallelise over the k point and spin
(see <a href="varpar.html#npkpt">npkpt</a>).
Otherwise, for unpolarized calculation at the gamma point, parallelise over the
two other levels: the band and FFT ones. For nproc<=50,
the best speed-up is achieved for
<a href="varpar.html#npband">npband</a>=nproc and
<a href="varpar.html#npfft">npfft</a>=1 (which is not yet the default).
For nproc>=50, the best speed-up is achieved for
<a href="varpar.html#npband">npband</a> >=4*<a href="varpar.html#npfft">npfft</a>.

<p>
For additional information,
download F. Bottin's presentation at the <a href="http://www.abinit.org/community/events/program3rd">ABINIT workshop 2007</a>
<p>
Suggested acknowledgments :
<br>
F. Bottin, S. Leroux, A. Knyazev and G. Zerah,
<i>Large scale ab initio calculations based on three levels of parallelization</i>,
Comput. Mat. Science <b>42</b>, 329 (2008),
available on arXiv, http://arxiv.org/abs/0707.3405 .
<p>
If the total number of processors used is compatible with the three levels of parallelization, the values for <a href="varpar.html#npband">npband</a>, <a href="varpar.html#npfft">npfft</a>, <a href="varpar.html#npband">npband</a> and <a href="varpar.html#bandpp">bandpp</a> will be filled automatically, although the repartition may not be optimal. To optimize the repartition use:
<p>

<b>If paral_kgb=-n </b>, ABINIT will test automatically if all the processor numbers between 2 and n are convenient for a parallel calculation and print the possible values in the log file. A weight is attributed to each possible processors repartition. It is adviced to select a processor repartition the weight of which is close to 1. The code will then stop after the printing. This test can be done as well with a sequential as with a parallel version of the code. The user can then choose the adequate number of processor on which he can run his job. He must put again paral_kgb=1 in the input file and put the corresponding values for <a href="varpar.html#npband">npband</a>, <a href="varpar.html#npfft">npfft</a>, <a href="varpar.html#npband">npband</a> and <a href="varpar.html#bandpp">bandpp</a> in the input file.





"""
},
'paral_rf': {
'definition': "activate PARALlelization over Response Function perturbations ",
'section': "varpar",
'category': "can even be specified separately for each dataset ",
'vartype': "integer ",
'default': "0.",
'text': """<p>
This parameter activates the parallelization over perturbations which can be used during
RF-Calculation. It is possible to use this type of parallelization in combination to the
parallelization over k-points.
</p>
<p>
Currently total energies calculated by groups, where the master process is not in, are saved
in .status_LOGxxxx files.
</p>


"""
},
'pawcpxocc': {
'definition': "PAW - use ComPleX rhoij OCCupancies ",
'section': "varpaw",
'category': "",
'vartype': "integer parameter ",
'default': """1, except for Ground-state static calculations (<a href="vargs.html#optdriver">optdriver</a>=0, <a href="varrlx.html#ionmov">ionmov</a><6) with spin-orbit coupling (<a href="varpaw.html#pawspnorb">pawspnorb</a>=1), density mixing (<a href="varbas.html#iscf">iscf</a>>=10) and time-reversal symmetry desactivated for k-points (<a href="varbas.html#kptopt">kptopt</a>/=1 or 2)""",
'text': """Relevant only when
<a href="varint.html#usepaw">usepaw</a>=1.
<br>
When <b>pawcpxocc=2</b>, PAW augmentation occupancies are treated as COMPLEX; else they are considered are REAL.<br>
This is needed when time-reversal symmetry is broken (typically when spin-orbit coupling is activated).<br><br>

Note for ground-state calculations (<a href="vargs.html#optdriver">optdriver</a>=0):<br>
The imaginary part of PAW augmentation occupancies is only used for the computation of the total energy by "direct scheme";
this only necessary when SCF mixing on potential is chosen (<a href="varbas.html#iscf">iscf</a><10).<br>
When SCF mixing on density is chosen (<a href="varbas.html#iscf">iscf</a>>=10), the "direct" decomposition
of energy is only printed out without being used. It is thus possible to use <b>pawcpxocc</b>=1 is the latter case.<br>
In order to save CPU time, when molecular dynamics is selected (<a href="varrlx.html#ionmov">ionmov</a>>=6) and
SCF mixing done on density (<a href="varbas.html#iscf">iscf</a>>=10), <b>pawcpxocc</b>=2 is (by default) set to <b>1</b>.<br>
When <b>pawcpxocc=1</b>, "direct" decomposition of total energy cannot be printed out.


"""
},
'pawcross': {
'definition': "PAW - add CROSS term in oscillator strengths ",
'section': "varpaw",
'category': "",
'vartype': "integer parameter ",
'default': "0 ",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=3 or 4, that is, screening or sigma calculations, and
<a href="varint.html#usepaw">usepaw</a>=1.
<p>
When <b>pawcross=1</b>, the overlap between the plane-wave part of one band and the on-site part of an other
is taken into account in the computation of the oscillator strengths. Hence, the completeness of the on-site basis is no longer assumed.


"""
},
'pawecutdg': {
'definition': "PAW - Energy CUToff for the Double Grid ",
'section': "varpaw",
'category': """<a href="../users/abinit_help.html#energy">ENERGY</a> """,
'vartype': "real parameter ",
'default': "-1, so pawecutdg MUST be specified for PAW calculations.",
'text': """Needed only when
<a href="varint.html#usepaw">usepaw</a>=1.
<br>Define the energy cut-off for the fine FFT grid
(the "double grid", that allows to transfer data from the normal, coarse, FFT grid to the
spherical grid around each atom).
<br><b>pawecutdg</b> must be larger or equal to
<a href="varbas.html#ecut">ecut</a>. If it is equal to it, then no fine
grid is used. The results are not very accurate, but the computations
proceed quite fast.
<br>For typical PAW computations, where <a href="varbas.html#ecut">ecut</a>
is on the order of 15 Ha, <b>pawecutdg</b> should be on the order of 40 Ha.
Choosing a larger value should not increase the accuracy, but does not slow down the
computation either, only the memory. The choice made for this variable DOES have a bearing
on the numerical accuracy of the results, and, as such, should be the object
of a convergence study. The convergence test might be made on the total energy
or derived quantities, like forces, but also on the two values of the
"Compensation charge inside spheres", a quantity written in the log file.
"""
},
'pawfatbnd': {
'definition': "PAW: print band structure in the FAT-BaND representation ",
'section': "varpaw",
'category': "",
'vartype': "integer parameter",
'default': "0",
'text': """Relevant only when
<a href="varint.html#usepaw">usepaw</a>=1
for Ground-State calculations and non self-consistent calculations.<br>
This option can be used to plot band structure. For each atom (specified by
<a href="vargs.html#natsph">natsph</a> and <a href="vargs.html#iatsph">iatsph</a>), each angular momentum, and each
spin polarisation, the band structure is written in files (such as e.g. FATBANDS_at0001_Ni_is2_l2_m-1). Each file
contains the eigenvalue, and the contribution of angular momentum L, and projection of angular momentum m,
(for the corresponding wavefunction) to the PAW density inside the PAW sphere
as a function of the index of the k-point.
The output can be readily plotted with the software <a href="http://plasma-gate.weizmann.ac.il/Grace/">xmgrace</a> (e.g xmgrace FATBANDS_at0001_Ni_is2_l2_m-1).
Relevant values are:<br>
<ul>
<li>0: desactivated.</li>
<li>1: The fatbands are only resolved in L.</li>
<li>2: The fatbands are resolved in L and M.</li>
</ul>
"""
},
'pawlcutd': {
'definition': "PAW - L angular momentum used to CUT the development in moments of the Densitites ",
'section': "varpaw",
'category': "",
'vartype': "integer parameter ",
'default': "10 ",
'text': """Relevant only when
<a href="varint.html#usepaw">usepaw</a>=1.
<br>The expansion of the densities in angular momenta is
performed up to l=<b>pawlcutd</b>.
<br>Note that, for a given system, the maximum value of <b>pawlcutd</b> is <b>2*l_max</b>,
where l_max is the maximum l of the PAW partial waves basis.
<br><br>The choice made for this variable DOES have a bearing
on the numerical accuracy of the results, and, as such, should be the object
of a convergence study. The convergence test might be made on the total energy
or derived quantities, like forces, but also on the two values of the
"Compensation charge inside spheres", a quantity written in the log file.
"""
},
'pawlmix': {
'definition': "PAW - maximum L used in the spherical part MIXing ",
'section': "varpaw",
'category': "",
'vartype': "integer parameter ",
'default': "10 ",
'text': """Relevant only when
<a href="varint.html#usepaw">usepaw</a>=1.
<br>
The choice made for this variable determine how the spherical part of the
density is mixed during electronic iterations.
<br>
<br>
Only parts of (rhoij quantities) associated with
l angular momenta up to l=pawlmix are mixed.
Other parts of augmentation occupancies are not included in the mixing process.
<br>
This option is useful to save CPU time but DOES have a bearing on the numerical accuracy of the results.
"""
},
'pawmixdg': {
'definition': "PAW - MIXing is done (or not) on the (fine) Double Grid ",
'section': "varpaw",
'category': "",
'vartype': "integer parameter ",
'default': """0 if <a href="varpar.html#npfft">npfft</a>=1 (else 1)""",
'text': """Relevant only when
<a href="varint.html#usepaw">usepaw</a>=1.
<br>
The choice made for this variable determines the grid on which the density (or potential) is mixed during the SCF cycle.
<br>
<br> - If <b>pawmixdg=1</b> the density/potential is mixed in REAL space using the fine FFT grid
(defined by <a href="varpaw.html#">pawecutdg</a> or <a href="varpaw.html#ngfftdg">ngfftdg</a>).
<br> - If <b>pawmixdg=0</b> the density/potential is mixed in RECIPROCAL space using the coarse FFT grid
(defined by <a href="varbas.html#ecut">ecut</a> or <a href="vargs.html#ngfft">ngfft</a>).
Only components of the coarse grid are mixed using the scheme defined by <a href="varbas.html#iscf">iscf</a>;
other components are only precondionned by <a href="vargs.html#diemix">diemix</a> and simply mixed.
<br>
This option is useful to save memory and does not affect numerical accuracy of converged results.
If <b>pawmixdg=1</b>, density and corresponding residual are stored for previous iterations and are
REAL arrays of size <a href="varint.html#nfftdg">nfftdg</a>.
If <b>pawmixdg=0</b>, density and corresponding residual are stored
for previous iterations and are COMPLEX arrays of size <a href="varint.html#nfft">nfft</a>.
The memory saving is particulary efficient when using the Pulay mixing
(<a href="varbas.html#iscf">iscf</a>=7 or 17).
"""
},
'pawnhatxc': {
'definition': "PAW - Flag for exact computation of gradients of NHAT density in eXchange-Correlation. ",
'section': "varpaw",
'category': "",
'vartype': "integer parameter ",
'default': "1 ",
'text': """Needed only when
<a href="varint.html#usepaw">usepaw</a>=1.
<br>Relevant only when a GGA exchange-correlation functional is used.
<br>When this flag is activated, the gradients of compensation charge density (n_hat) are exactly
computed (i.e. analytically); when it is desactivated, they are computed with a numerical scheme in
reciprocal space (which can produce inaccurate results if the compensation charge density
is highly localized).
<br>As analytical treatment of compensation charge density gradients is CPU time demanding,
it is possible to bypass it with <b>pawnhatxc</b>=0; but the numerical accuracy can be affected
by this choice. It is recommended to test the validity of this approximation before use.
"""
},
'pawnphi': {
'definition': "PAW - Number of PHI angles used to discretize the sphere around each atom. ",
'section': "varpaw",
'category': "",
'vartype': "integer parameter ",
'default': "13 ",
'text': """Needed only when
<a href="varint.html#usepaw">usepaw</a>=1.
<br>Number of phi angles (longitude) used to discretize the
data on the atomic spheres. This discretization is completely
defined by <b>pawnphi</b>
and <a href="varpaw.html#pawntheta">pawntheta</a>.
"""
},
'pawntheta': {
'definition': "PAW - Number of THETA angles used to discretize the sphere around each atom. ",
'section': "varpaw",
'category': "",
'vartype': "integer parameter ",
'default': "12 ",
'text': """Relevant only when
<a href="varint.html#usepaw">usepaw</a>=1.
<br>Number of theta angles (latitude) used to discretize the
data on the atomic spheres. This discretization is completely
defined by <b>pawntheta</b>
and <a href="varpaw.html#pawnphi">pawnphi</a>.
"""
},
'pawnzlm': {
'definition': "PAW - only compute Non-Zero LM-moments of the contributions to the density from the spheres",
'section': "varpaw",
'category': "",
'vartype': "integer parameter ",
'default': "1 ",
'text': """Relevant only when
<a href="varint.html#usepaw">usepaw</a>=1.
<br> Concerns the computation of the contributions to the density from the spheres (named rho_1 -
rho_tild_1).
<br>If set to 0, all lm-moments of the sphere contributions to the density are computed at each
electronic iteration.
<br> If set to 1, only non-zero lm-moments of the sphere contributions to the density are computed
at each electronic iteration (they are all computed at the first iteration then only
those found to be non-zero will be computed ; thus the first iteration is more cpu intensive)
"""
},
'pawoptmix': {
'definition': "PAW - OPTion for the MIXing of the spherical part",
'section': "varpaw",
'category': " ",
'vartype': "integer parameter ",
'default': "0",
'text': """When PAW is activated, the self-consistent requires the mixing of both the total potential (or density)
and the "spherical part" (in fact the augmentation occupancies rho_ij).
<br>The same mixing scheme is applied to the potential (density) and the spherical part.
It is optimized in order to minimize a residual.
<br>If <b>pawoptmix</b>=0 the residual is the potential (or density) residual.
<br>If <b>pawoptmix</b>=1 the residual is a sum of the potential (or density) residual
and the "spherical part" residual.
"""
},
'pawovlp': {
'definition': "PAW - spheres OVerLap allowed (in percentage)",
'section': "varpaw",
'category': " ",
'vartype': "real parameter ",
'default': "5.0",
'text': """When PAW is activated, a localized atomic basis is added to describe wave functions.
Spheres around atoms are defined and they are IN PRINCIPLE not allowed to overlap. However,
a small overlap can be allowed without compromising the accurary of results. Be aware that
too high overlaps can lead to unphysical results.<br>
With the <b>pawovlp</b> variable, the user can control the (voluminal) overlap percentage
allowed without stopping the execution.<br>
<b>pawovlp</b> is the value (in percentage: 0...100%) obtained by dividing
the volume of the overlap of two spheres by the volume of the smallest sphere.<br>
The following values are permitted for <b>pawovlp</b>:<br>
<div style="margin-left: 40px;">
- <b>pawovlp</b><0. : overlap is always allowed<br>
- <b>pawovlp</b>=0. : no overlap is allowed<br>
- <b>pawovlp</b>>0. and <100. : overlap is allowed only if it is less than <b>pawovlp</b> %
</div>"""
},
'pawprt_b': {
'definition': "PAW print band",
'section': "vardev",
'category': " DEVELOP ",
'vartype': "integer",
'default': "0 ",
'text': """Forces the output of the all-electron wavefunction for
only a single band. To be used in conjuction with:<b><br />
<a href="varpaw.html#pawprtwf">pawprtwf</a>=1</b> and
<a href="vardev.html#pawprt_k">pawprt_k</a>.
The indexing of the bands start with one for the lowest occupied band
and goes up from there.
"""
},
'pawprt_k': {
'definition': "PAW print k-point",
'section': "vardev",
'category': "DEVELOP",
'vartype': "integer",
'default': "0 ",
'text': """Forces the output of the all-electron wavefunction for
only a single k-point. To be used in conjuction with:<b><br />
<a href="varpaw.html#pawprtwf">pawprtwf</a>=1</b> and
<a href="vardev.html#pawprt_b">pawprt_b</a>.
The indexing follows the order in ouptput of the internal
variable <b>kpt</b> in the beginning of the run.
"""
},
'pawprtden': {
'definition': "PAW: PRinT total physical electron DENsity",
'section': "varpaw",
'category': " ",
'vartype': "integer parameter ",
'default': "0",
'text': """<b>Deprecated :</b> See the <a href="varfil.html#prtden">prtden</a>.
"""
},
'pawprtdos': {
'definition': "PAW: PRinT partial DOS contributions",
'section': "varpaw",
'category': " ",
'vartype': "integer parameter ",
'default': "0",
'text': """This input variable controls the computation and/or printing of contributions to the PAW partial DOS in _DOS file(s):<br>
<div style="margin-left: 40px;">
+ Plane-waves contribution<br>
+ "on-site" all-electron contribution (phi)<br>
- "on-site" pseudo contribution (phi_tild).<br>
</div>
If <b>pawprtdos=0:</b><br>
- The 3 contributions are computed; only the total partial DOS is output in _DOS file.<br>
If <b>pawprtdos=1:</b><br>
- The 3 contributions are computed and output in _DOS file.<br>
- In that case, integrated DOS is not output.<br>
If <b>pawprtdos=2:</b><br>
- Only "on-site" all-electron contribution is computed and output in _DOS file.<br>
- This a (very) good approximation of total DOS, provided that (1) the PAW local basis is complete, (2) the electronic charge is mostly contained in PAW spheres.<br>
- In that case, the <a href="vargs.html#ratsph">ratsph</a> variable is automatically set to the PAW radius.
"""
},
'pawprtvol': {
'definition': "PAW: PRinT VOLume",
'section': "varpaw",
'category': " ",
'vartype': "integer parameter ",
'default': "0",
'text': """Control print volume and debugging output for PAW in log file or standard output.
If set to 0, the print volume is at its minimum.<br>
<b>pawprtvol</b> can have values from -3 to 3:<br>
- <b>pawprtvol</b>=-1 or 1: matrices rho_ij (atomic occupancies) and D_ij (psp strength)
are printed at each SCF cycle with details about their contributions. <br>
- <b>pawprtvol</b>=-2 or 2: like -1 or 1 plus additional printing: moments of "on-site" densities,
details about local exact exchange.<br>
- <b>pawprtvol</b>=-3 or 3: like -2 or 2 plus additional printing: details about
PAW+U, rotation matrices of sphercal harmonics.<br>
When <b>pawprtvol</b>>=0, up to 12 components of rho_ij and D_ij matrices for the 1st and last atom are printed.<br>
When <b>pawprtvol</b><0, all components of rho_ij and D_ij matrices for all atoms are printed.<br>
"""
},
'pawprtwf': {
'definition': "PAW: PRinT WaveFunctions",
'section': "varpaw",
'category': " ",
'vartype': "integer parameter ",
'default': "0",
'text': """This input variable controls the output of the full PAW wave functions including the on-site
contribution inside each PAW sphere needed to reconstruct the correct nodal shape in the augmentation region.
<b>pawprtwf=1</b> causes the generation of a file  _AE_WFK that contains the full
wavefunctions in real space on the fine FFT grid defined by
<a href="varpaw.html#pawecutdg">pawecutdg</a> or <a href="varpaw.html#ngfftdg">ngfftdg</a>).
Limitations: At present (v6.0), <b>pawprtwf=1</b> is not compatible neither with the k-point parallelism
nor with the parallelism over G-vectors. Therefore the output of the _AE_WFK has to be done in sequential.
Moreover, in order to use this feature, one has to be enable the support for ETSF-IO at configure-time
as the _AW_WFK file is written using the NETCDF file format following the ETSF-IO specification for
wavefunctions in real space. If the code is run in entirely in serial, additional output is made of various
contributions to the all-electron avefunction. By default the full available set of bands and k-points are
ouput, but a single band and k-point index can be requested by using the variables
<a href="vardev.html#pawprt_b">pawprt_b</a> and <a href="vardev.html#pawprt_k">pawprt_k</a>.
"""
},
'pawspnorb': {
'definition': "PAW - option for SPiN-ORBit coupling",
'section': "varpaw",
'category': " ",
'vartype': "integer parameter ",
'default': "0",
'text': """When PAW is activated, the <b>spin-orbit coupling</b> can be added without the use
of specific PAW datasets (pseudopotentials).<br>
If <b>pawspnorb</b>=1, spin-orbit will be added.<br><br>
Note that only the all-electron "on-site" contribution to the Hamiltonian is taken into account;
this a very good approximation but requires the
following conditions to be fullfilled:
<br> <div style="margin-left: 20px;">1- the<span
style="position: relative; top: -8pt; left: 8pt;"></span><span
style="position: relative; top: -5pt; left: 6pt;">~</span>&#966;<sub>i</sub>
basis is complete enough<br> 2- the electronic density is mainly
contained in the PAW sphere<br>
</div><br>
Also note that, when spin-orbit coupling is activated, the time-reversal symmetry might be broken.<br>
The use of <a href="varbas.html#kptopt">kptopt</a>=1 or <a href="varbas.html#kptopt">kptopt</a>=2 is thus
forbidden. It is adviced to use <a href="varbas.html#kptopt">kptopt</a>=3 (no symmetry used to generate k-points)
or <a href="varbas.html#kptopt">kptopt</a>=4 (only spatial symmetries used to generate k-points).<br>
Be careful if you choose to use <a href="varbas.html#kptopt">kptopt</a>=0 (k-points given by hand); Time-reversal symmetry
has to be avoided.
<br>An artificial scaling of the spin-orbit can be introduced thanks to the <a href="varpaw.html#spnorbscl">spnorbscl</a> input variable.
"""
},
'pawstgylm': {
'definition': "PAW - option for the STorage of G_l(r).YLM(r)",
'section': "varpaw",
'category': " ",
'vartype': "integer parameter ",
'default': "1",
'text': """When PAW is activated, the computation of compensation charge density (so called
"hat" density) requires the computation of g_l(r).Y_lm(r) factors (and cartesian derivatives) at each point
of real space contained in PAW spheres. The number of atoms, of (l,m) quantum numbers
and the sharpness of the real FFT grid can lead to a very big {g_l.Y_lm} datastructure.
One can save memory by putting <b>pawstgylm</b>=0; but, in that case, g_l(r).Y_lm(r) factors
a re-computed each time they are needed and CPU time increases.
<br><br><div style="margin-left: 20px;">
Possible choices:<br>
- <b>pawstgylm</b>=0 : g_l(r).Y_lm(r) are not stored in memory and recomputed.<br>
- <b>pawstgylm</b>=1 : g_l(r).Y_lm(r) are stored in memory.<br>
</div>
<br><div style="margin-left: 40px;">
Note:<br>
g_l(r) are shape functions (analytically known)<br>
Y_lm(r) are real spherical harmonics
</div>
"""
},
'pawsushat': {
'definition': "PAW - SUSceptibility, inclusion of HAT (compensation charge) contribution",
'section': "varpaw",
'category': " ",
'vartype': "integer parameter ",
'default': "0",
'text': """When a sophisticated preconditioning scheme is selected for the SCF cycle of a Ground-State calculation
(<a href="vargs.html#iprcel">iprcel</a>>0), the computation of the susceptibility matrix
is required several times during the cycle. This computation is computer time consuming, especially
-- within PAW -- because of the inclusion of additional terms due to the compensation charge density.
As only a crude valuation of the susceptibilty matrix is needed (to evaluate a preconditioning matrix),
the compensation charge contribution can be neglected to save CPU time (select <b>pawsushat</b>=0).
This approximation could be unfavourable in some cases; in the latter, we advice to put <b>pawsushat</b>=1.
<br><br><div style="margin-left: 20px;">
Possible choices:<br>
- <b>pawsushat</b>=0 : only plane-wave contribution to suscep. matrix is computed.<br>
- <b>pawsushat</b>=1 : the whole suscep. matrix (PW + PAW on-site) is computed.<br>
</div>
"""
},
'pawujat': {
'definition': "PAW+macro_UJ, ATom number ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>""",
'vartype': "integer",
'default': "1, i.e. the first atom treated with PAW+U. ",
'text': """Determines the atom for which U (or J) should be determined. See also <a href="vardev.html#macro_uj">macro_uj</a>.
"""
},
'pawujrad': {
'definition': "PAW+macro_UJ, sphere RADius",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>""",
'vartype': "real ",
'default': "20 a.u. ",
'text': """The sphere radius serves to extrapolate the U value calculated at r_paw to a larger sphere radius.
See also <a href="vardev.html#macro_uj">macro_uj</a>.
As most projector functions are localized within r_paw to &asymp;80%,
20 a.u. contains &asymp;100% of the wavefunction and corresponds to r_paw &rarr; &infin;.
"""
},
'pawujv': {
'definition': "PAW+macro_UJ, potential shift (V)",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>""",
'vartype': "real,  ",
'default': "0.1 eV. ",
'text': """Amplitude of the potential shift for the determination of U (or J). See also <a href="vardev.html#macro_uj">macro_uj</a>.
"""
},
'pawusecp': {
'definition': "PAW - option for the USE of CPrj in memory (cprj=WF projected with NL projector)",
'section': "varpaw",
'category': " ",
'vartype': "integer parameter ",
'default': "1",
'text': """When PAW is activated, the computation of cprj arrays is memory and time consuming.<br>
When <b>pawusecp</b>=0, then the cprj are never kept in memory, they are recomputed when needed (this is CPU-time consuming).
When <b>pawusecp</b>=1, then the cprj are computed once and then kept in memory.<br>
Change the value of the keyword only if you are an experienced user (developper).<br>
Remember: cprj = (WF_n .dot. p_i) (WF_n=wave function, p_i=non-local projector).<br>
<br>For the time being, only activated for RF calculations.<br>
"""
},
'pawxcdev': {
'definition': "PAW - choice for eXchange-Correlation DEVelopment (spherical part) ",
'section': "varpaw",
'category': "",
'vartype': "integer parameter ",
'default': "1 ",
'text': """Relevant only when
<a href="varint.html#usepaw">usepaw</a>=1.
<ul>
<li> If set to 0, the exchange-correlation term in the spherical part of energy is
totally computed on the angular mesh </li>
<li> If set to 1, the exchange-correlation term in the spherical part of energy is
developed onto lm-moments at order 1 </li>
<li> If set to 2, the exchange-correlation term in the spherical part of energy is
developed onto lm-moments at order 2 (can be memory/CPU consuming) </li>
</ul>
<br>
Be careful: GGA requires <b>pawxcdev</b> &gt; 0
"""
},
'pimass': {
'definition': "Path Integral inertial MASSes",
'section': "varrlx",
'category': "",
'vartype': """real array pimass(<a href="varbas.html#ntypat">ntypat</a>) """,
'default': """<a href="varrlx.html#amu">amu</a>(<a href="varbas.html#ntypat">ntypat</a>).""",
'text': """Only relevant if <a href="varrlx.html#imgmov">imgmov</a>=9 or 13 (Path-Integral Molecular Dynamics).<br>
Gives the inertial masses (called "fictitious masses" in
<i>D. Marx and M. Parrinello, J. Chem. Phys. 104, 4077 (1996)</i>)
in atomic mass units for each kind of atom in cell. These masses are used
in performing molecular dynamical atomic motion.<br>
See <a href="varrlx.html#amu">amu</a> for further details.
"""
},
'pitransform': {
'definition': "Path Integral TRANSFORMation",
'section': "varrlx",
'category': "",
'vartype': "integer parameter ",
'default': "0.",
'text': """Only relevant if <a href="varrlx.html#imgmov">imgmov</a>=13 (Path-Integral Molecular Dynamics).
Allows to change the dynamics of the Path Integral chains, as described by M. Tuckerman et al, J. Chem. Phys. 104, 5579 (1996).
<p>
If equal to 1, staging transformation.
<br>
If equal to 2, normal modes transformation.
"""
},
'polcen': {
'definition': "POLarization for Centrosymmetric geometry",
'section': "varff",
'category': "",
'vartype': "real array polcen(3) ",
'default': "3*0",
'text': """When doing a finite electric displacement field calculation,
if the structure is centrosymmetric but the polarization is non-zero (such as for AlAs),
this non-zero polarization should be specified as <b>polcen</b> (in atomic units) in the input file.
"""
},
'positron': {
'definition': "POSITRON calculation ",
'section': "vargs",
'category': " ",
'vartype': "integer ",
'default': "0",
'text': """This input parameter can be positive or negative.<br>
Negative <b>positron</b> are only relevant for PAW calculation.<br>
Electron-positron correlation functional is defined by <a href="vargs.html#ixcpositron">ixcpositron</a>.<br>

<br><U>Positive values for <b>positron</b>:</U><br>
<i>For <b>positron=1 or 2</b>, will perform the calculation of positron
lifetime (and annihilation rate).</i><br>
<ul>
<li><b>positron=1</b>:<br>
Starting from a previous electronic GS density (with <b>positron=0</b>),
a positronic ground-state calculation is performed, considering that the electrons are not perturbed
by the presence of the positron.<br>
This is almost correct for a positron in a perfect bulk material.
But this approximation fails when defects exist in
the material (for instance: the positron might be trapped by a vacancy).<br>
The electronic density will be automatically read from a _DEN file (with or without
<a href="varfil.html#getden">getden</a>/<a href="varfil.html#irdden">irdden</a>
keyword).<br>
At the end of the SCF cycle, the positron lifetime and annihilation rate
are printed out.<br><br>
<i>Additional information for the use of pseudopotentials:<br><ul>
<li>PAW datasets: nothing to do; simply use usual electronic PAW datasets
<li>Norm-conserving pseudopotentials: One has to use specifics
pseudopotentials for the positron calculation.
They must be of the FHI type (pspcod=6), and must contain at their
end, the all-electrons core density generated with FHI98PP. They must have
lmax=lloc=0 (check that this works for the electronic GS !! No ghost, etc ...). Otherwise, their are similar to a usual FHI pseudopotential.<br>
</i></ul><br>

<li><b>positron=2</b>:<br>
Starting from a previous positronic GS density (with <b>positron=1</b>),
an electronic ground-state calculation is performed, keeping the positronic
density constant.<br>
The positronic density will be automatically read from a _DEN file (with or without
<a href="varfil.html#getden">getden</a>/<a href="varfil.html#irdden">irdden</a>
keyword).<br>
At the end of the SCF cycle, the positron lifetime and annihilation rate
are printed out.<br><br>
<i>Additional information for the use of pseudopotentials:<br><ul>
<li>PAW datasets: nothing to do; simply use usual electronic PAW datasets
<li>Norm-conserving pseudopotentials: One has to use specifics
pseudopotentials for the positron calculation.
They must be of the FHI type (pspcod=6), and must contain at their
end, the all-electrons core density generated with FHI98PP.<br>
</i></ul><br>

<li><b>Typical use</b>:<br>
The calculation is done in several steps:<br>
The first one is a normal GS calculation for the electrons, with <b>positron</b>=0.
The only specific thing to do is to set <a href="varfil.html#prtden">prtden</a>=1.
This will create the associated _DEN file which will used as input file for the positronic
GS calculation.<br>
The second step is the GS calculation of the positron and subsequently its
lifetime, with <b>positron</b>=1.
One has to define also <a href="vargs.html#ixcpositron">ixcpositron</a>.<br>
Then, it is possible to perform an additional step, computing
the GS electronic density in presence of the positron, with <b>positron</b>=2.<br>
and so on...<br>
This procedure can be automated (for PAW only) by the use of a negative value for <b>positron</b>.<br>
</ul>

<br><U>Negative values for <b>positron</b>:</U><br>
<i>For <b>positron&lt;0</b>, will perform an automatic calculation of electrons and positron
densities in the two-component DFT context; then will compute positron lifetime (and annihilation rate).</i><br>

<ul><li><b>positron=-1</b>:<br>
Starting from scratch, will first perform a usual electronic ground-state
calculation until convergency (controlled by the use of one of the <i>tolerance</i> keywords).<br>
Then will perform a positronic ground state calculation in presence of the electrons
and ions; then an electronic ground state calculation in presence of the positron and the ions...<br>
and so on... until the total energy is converged.<br>
The convergence of the total energy of the ions+electrons+positron system
is controlled by the use of the <a href="vargs.html#postoldfe">postoldfe</a> input keyword.<br>
With <b>positron=-1</b>, at the beginning of each new electronic/positronic step,
the wave functions are unknown.</ul>

<ul><li><b>positron=-10</b>:<br>
Same as <b>positron=-1</b> except that the electronic/positronic wave functions
are stored in memory.<br>
Consequently, the total number of iterations to reach the convergency (Etotal<<a href="vargs.html#postoldfe">postoldfe</a>) is smaller.<br>
But, this can increase the total amount of memory needed by the code.</ul><br>

References:<br>
<ul>
<b>[1]</b> J. Arponen and E. Pajanne, Ann. Phys. (N.Y.) 121, 343 (1979).<br>
<b>[2]</b> Boronski and R.M. Nieminen, Phys. Rev. B 34, 3820 (1986).<br>
<b>[3]</b> P.A. Sterne and J.H. Kaiser, Phys. Rev. B 43, 13892 (1991).<br>
<b>[4]</b> M.J. Puska, A.P. Seitsonen and R.M. Nieminen, Phys. Rev. B 52, 10947 (1994).<br>
<b>[5]</b> B. Barbiellini, M.J. Puska, T. Torsti and R.M.Nieminen, Phys. Rev. B 51, 7341 (1994)<br>
</ul>
"""
},
'posnstep': {
'definition': "POSitron calculation: max. Number of STEPs for the two-component DFT",
'section': "vargs",
'category': "",
'vartype': "integer parameter",
'default': "50",
'text': """Relevant only when <a href="vargs.html#positron">positron</a>&lt;0.<br>
Sets the maximum number of electronic/positronic iterations that, when reached, will cause the two-component DFT SCF cycle to stop.<br>
The code will first compute the electronic ground-state, then the positronic ground state in the electronic density, then the electronic ground-state in the positronic density, ...<br>
...until diff_Etotal<<a href="vargs.html#postoldfe">postoldfe</a> or diff_Forces<<a href="vargs.html#postoldff">postoldff</a>
or the number of electronic/positronic steps is <b>posnstep</b>.<br>
"""
},
'posocc': {
'definition': "POSitron calculation: OCCupation number for the positron",
'section': "vargs",
'category': "",
'vartype': "real parameter ",
'default': "1.",
'text': """Relevant only when <a href="vargs.html#positron">positron</a>/=0.<br>
Sets the occupation number for the positron. Has to be <=1.<br>
Changing <b>posocc</b> is only usefull for bulk calculation when one wants to perform lifetime computations using a small simulation cell (can avoid the use of a supercell). It simulates the dispersion of the positron in the whole crystal.<br>
"""
},
'postoldfe': {
'definition': "POSITRON calculation: TOLerance on the DiFference of total Energy",
'section': "vargs",
'category': """<a href="../users/abinit_help.html#energy">ENERGY</a> """,
'vartype': "real parameter ",
'default': "1.e-6",
'text': """Relevant only when <a href="vargs.html#positron">positron</a>&lt;0.<br>
Sets a tolerance for absolute difference of total energy (of <i>ions+electrons+positron</i> system)
that, when reached, will cause the SCF cycle to stop before the number of
steps is <a href="varbas.html#nstep">nstep</a>.<br>
<br>Can be specified in Ha (the default), Ry, eV or Kelvin, since
<b>toldfe</b> has the '<a href="../users/abinit_help.html#dimensions">ENERGY</a>' characteristics.
"""
},
'postoldff': {
'definition': "POSitron calculation: TOLerance on the DiFference of Forces",
'section': "vargs",
'category': "",
'vartype': "real parameter ",
'default': "0",
'text': """Relevant only when <a href="vargs.html#positron">positron</a>&lt;0.<br>
Sets a tolerance for absolute difference of maximum force acting on ions (due to <i>ions+electrons+positron</i> system)
that, when reached, will cause the SCF cycle to stop before the number of SCF
steps is <a href="varbas.html#nstep">nstep</a> or the number of electronic/positronic steps is <a href="vargs.html#posnstep">posnstep</a>.<br>
Only one and only one of <a href="vargs.html#postoldfe">postoldfe</a> or <a href="vargs.html#postoldff">postoldff</a> can be set.<br>
"""
},
'ppmfrq': {
'definition': "Plasmon Pole Model FReQuency ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>, <a href="../users/abinit_help.html#energy">ENERGY</a> """,
'vartype': "real parameter",
'default': "0.0 Ha ",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=3 or 4, that is, screening
calculations or sigma calculations.
Usually only effective if GW corrections are evaluated using the plasmon-pole model of Godby-Needs
(<i>i.e</i> <a href="vargw.html#ppmodel">ppmodel</a>=1).
<br><br>
<B> In plasmon-pole calculations</B>
<p>
In the present status of the GW code, the convolution in frequency space
defining the self-energy operator can be evaluated using two different approaches:
numerical integration and plasmon-pole models.
<br>
Methods based on the numerical integration (contour deformation, analytic continuation) require
the knowledge of the screened interaction for several frequencies. These approaches give
the most accurate results but at the price of an increase in the CPU time required.
<br>
Alternatively, it is possible to approximate the dynamical behaviour of the screened interaction
through simple analytical expressions, the so-called plasmon-pole models.
In the plasmon-pole model proposed by Godby-Needs (<a href="vargw.html#ppmodel">ppmodel</a>=1),
the screening must be available at zero frequency, as well as at another imaginary frequency,
of the order of the plasmon frequency (the peak in the EELS spectrum).
This information is used to model the behaviour of the dielectric matrix
for all frequencies.
During the calculation of the screening, <b>ppmfrq</b> defines the imaginary frequency where the
dielectric matrix is evaluated, in addition to the zero frequency.
During the self-energy run, <b>ppmfrq</b> can be used to define the second frequency to be used
to calculate the plasmon-pole parameters. This is particularly useful when the
SCR file contains several frequencies along the imaginary axis.
In this case the frequency whose value is the closest one to <b>ppmfrq</b> will be selected.
Note that, if the plasmon-pole approximation is good, then, the
choice of <b>ppmfrq</b> should have no influence on the final result.
One should check whether this is the case. In general, the plasmon frequencies of bulk solids
are of the order of 0.5 Hartree.
<br><br>
<B> In Contour Deformation calculations</B>
<p>
<b>ppmfrq</b> is here used to <b>override</b> the default value calculated from the average electronic
density per unit cell. This can affect the distribution of gridpoints along the imaginary and
real frequency axes. See <a href="vargw.html#cd_frqim_method">gw_frqim_method</a>, <a href="vargw.html#gw_frqim_inzgrid">gw_frqim_inzgrid</a> and <a href="vargw.html#gw_frqre_inzgrid">gw_frqre_inzgrid</a> for more details.

"""
},
'ppmodel': {
'definition': "Plasmon Pole MODEL ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>  """,
'vartype': "integer parameter ",
'default': "1 ",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=3 or 4, that is, screening calculations and self-energy calculations.

<ul>
<li> <b>ppmodel</b>=1 : PP model of Godby and Needs,
See <B> FIXME>REFERENCE MISSING!</B></li>
<li> <b>ppmodel</b>=2 : PP model of Hybertsen and Louie,
See Phys Rev B 34, 5390 (1986) </li>
<li> <b>ppmodel</b>=3 : PP model of  W. von der Linden
and P. Horsh see Phys Rev B 37, 8351 (1988) </li>
<li> <b>ppmodel</b>=4 : PP model of Farid and Engel.
See Phys Rev B47,15931 (1993)
<li> <b>ppmodel</b>=0 : no PP model, numerical integration
(contour deformation method, see e.g. S. Lebegue  et al. PRB  67, 155208 (2003).)  </li>
</ul>

Please note the difference between <b>ppmodel</b> 1 and <b>ppmodel</b> 2,3,4.
In the first case (<b>ppmodel</b>=1), the plasmon-pole parameters are determined in order to reproduce
the behaviour of the dielectric matrix at two calculated frequencies: the static limit (omega=0)
and the imaginary frequency defined by <a href="vargw.html#ppmfrq">ppmfrq</a>.
In the last three cases, instead, the plasmon-pole parameters are found
by using the dielectric matrix calculated only at omega=0 and enforcing the so-called f-sum rule.
See also <a href="vargw.html#nfreqre">nfreqre</a>.
<p>
Please note also that in the case of <b>ppmodel</b> 4, the plasmon energies are not simple mathematical parameters,
but rather have a physical meaning (at least the lowest ones).
Thus the calculated plasmon band structure (plasmon energy vs q vector) is reported in the output file for the
lowest 10 bands.
"""
},
'prepanl': {
'definition': "PREPAre Non-Linear response calculation ",
'section': "varrf",
'category': """<a href="../users/abinit_help.html#respfn">RESPFN</a>  """,
'vartype': "integer parameter  ",
'default': "0.",
'text': """The computation of third-order derivatives from the 2n+1 theorem
requires the first-order wavefunctions and densities obtained from
a linear response calculation. The standard approach in a linear
response calculation is (i) to compute only the
irreducible perturbations, and (ii) to use symmetries to
reduce the number of k-points for the k-point integration.
<br>This approach cannot be applied, presently (v4.1),
if the first-order wavefunctions are to be used to compute third-order derivatives.
First, for electric fields, the code needs the derivatives
along the three directions. Still, in case of phonons, only the
irreducible perturbations are required.
Second, for both electric fields and phonons, the wavefunctions
must be available in half the BZ (kptopt=2), or the full BZ (kptopt=3).
<br>During the linear response calculation, in order to prepare a non-linear
calculation, one should put <b>prepanl</b> to 1 in order
to force ABINIT (i) to compute the electric field perturbation
along the three directions explicitly, and (ii) to keep the full number of k-points.
"""
},
'prepgkk': {
'definition': "PREPAre GKK calculation",
'section': "varrf",
'category': """<a href="../users/abinit_help.html#respfn">RESPFN</a>  """,
'vartype': "integer parameter  ",
'default': "0.",
'text': """The calculation of electron-phonon coupling quantities requires the presence
of all the perturbations (all atoms in all directions) for the chosen set
of (irreducible) q-points. To impose this and prevent ABINIT from using
symmetry to reduce the number of perturbations, set <b>prepgkk</b> to 1.
Use in conjunction with <a href="varfil.html#prtgkk"><b>prtgkk</b></a>.
"""
},
'prepscphon': {
'definition': "PREPare Self-Consistent PHONon calculation",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>""",
'vartype': "integer",
'default': "0 ",
'text': """Print PCINFO, PHFREQ, and PHVEC files, for use with self-consistent phonon runs, after a perturbation
calculation. Only prints out files for the present q-point, and there is presently no tool to symmetrize
or merge these files, so use anaddb instead (with prtscphon input variable). The abinit input
variable is destined to someday bypass the use of anaddb for scphon calculations.
"""
},
'prt1dm': {
'definition': "PRinT 1-DiMensional potential and density ",
'section': "varfil",
'category': " ",
'vartype': "integer parameter ",
'default': "0.",
'text': """If set >= 1, provide one-dimensional projection of
potential and density, for each of the three axis.
This corresponds to averaging the potential
or the density on bi-dimensional slices of the FFT grid."""
},
'prtatlist': {
'definition': "PRinT by ATom LIST of ATom",
'section': "varrlx",
'category': """<a href="../users/abinit_help.html#no_multi">NO MULTI</a> """,
'vartype': "integer array ",
'default': "0.",
'text': """This is an array of the numbers associated to the index atoms that
the user want to print in the output or log files, this is useful when
you have a large number of atoms and you are only interested to
follow especific atoms, the numbers associated should be consistent
with the list in <a href="varbas.html#xcart">xcart</a> or
<a href="varbas.html#xcart">xred</a>

This input varible does not affect the contents of the "OUT.nc" or
"HIST.nc", those are NetCDF files containing the information about
all the atoms.
"""
},
'prtbbb': {
'definition': "PRinT Band-By-Band decomposition ",
'section': "varrf",
'category': """<a href="../users/abinit_help.html#respfn">RESPFN</a>  """,
'vartype': "integer parameter  ",
'default': "0.",
'text': """If <b>prtbbb</b> is 1, print the band-by-band decomposition of
Born effective charges and localization tensor, in case they are computed.
See Ph. Ghosez and X. Gonze, J. Phys.: Condens. Matter 12, 9179 (2000).
"""
},
'prtbltztrp': {
'definition': "PRinT output for BoLTZTRaP code",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>""",
'vartype': "integer",
'default': "0 ",
'text': """Print out geometry (_BLZTRP_GEOM) and eigenenergy (_BLZTRP_EIGEN) files for the
<a href="http://www.icams.de/content/departments/ams/madsen/boltztrap.html">BoltzTraP code</a> by Georg Madsen.
"""
},
'prtcif': {
'definition': "PRinT Crystallographic Information File",
'section': "vardev",
'category': " DEVELOP",
'vartype': "integer flag",
'default': "0 ",
'text': """If set to 1, a CIF file is output with the crystallographic data for the present run (cell size shape and atomic positions).
"""
},
'prtcml': {
'definition': "PRinT CML file ",
'section': "varfil",
'category': " ",
'vartype': "integer parameter ",
'default': "0.",
'text': """If set to 1 or 2, provide output of geometrical
parameters using CML (the Chemical Markup Language,
see papers by P. Murray-Rust and H. Rzepa, especially
J. Chem. Inf. Comput. Sci. 39, 928-942 (1998) and the Web site
<a href="http://cml.sourceforge.net">http://cml.sourceforge.net</a>).
<br>If prtcml==1, atomic positions are written in fractional coordinates (attributes
'xFract', 'yFract', and 'zFract').
<br>If prtcml==2, atomic positions are written in cartesian coordinates (angstrom units,
attributes 'x3', 'y3', and 'z3').
<br>Such file can be treated automatically by tools developed
to handle XML formatted files.
<br>
<br>Such a CML file contains :
<ul>
<li>The crystallographic information (space group number and the needed
unit cell parameters and angles)</li>
<li>The list of symmetry elements</li>
<li>The list of atoms in the cell (symbols and fractional or cartesian coordinates in angstroms)</li>
</ul>
<br>If <a href="varrlx.html#ionmov">ionmov</a>==0, the name of the CML file will be
the root output name, followed by _CML.xml .
<br>If <a href="varrlx.html#ionmov">ionmov</a>==1 or 2, the CML file will be output
at each time step, with the name being made of
<ul>
<li> the root output name,</li>
<li> followed by _TIMx , where
x is related to the timestep (see later)</li>
<li> then followed by _CML.xml </li>
</ul>
<br>No output is provided by <b>prtcml</b> different from 1 or 2
"""
},
'prtcs': {
'definition': "PRinT Chemical Shielding ",
'section': "varpaw",
'category': "",
'vartype': "integer parameter ",
'default': "0 ",
'text': """Relevant only when
<a href="varint.html#usepaw">usepaw</a>=1.
<ul>
<li> If set to 1, calculate the chemical shielding at each atomic site in the unit cell. THIS CODE
IS UNDER DEVELOPMENT AND IS NOT READY FOR USE.
</li>
</ul>

"""
},
'prtden': {
'definition': "PRinT the DENsity ",
'section': "varfil",
'category': " ",
'vartype': "integer parameter ",
'default': """1. (was 0 before v5.3), except when <a href="varrlx.html#nimage">nimage</a>>1""",
'text': """If set to 1  or a larger value , provide output of electron density
in real space rho(r), in units of electrons/Bohr^3.
<br>If <a href="varrlx.html#ionmov">ionmov</a>==0, the name of the density file will be
the root output name, followed by _DEN .
<br>If <a href="varrlx.html#ionmov">ionmov</a>==1 or 2, density files will be output
at each time step, with the name being made of
<ul>
<li> the root output name,</li>
<li> followed by _TIMx , where
x is related to the timestep (see later)</li>
<li> then followed by _DEN</li>
</ul>
The file structure of the unformatted output file is
described below, see section 6).
<br>If  <b>prtden</b> is lower than 0, two files will be printed for restart every <b>prtden</b>  step,
with the names being made of
<ul>
<li> the root temporary name,</li>
<li> followed by _DEN_x , where
x is 0000 or 0001 alternatively.</li>
<li> The most recent of the two files should be used for restart, and copied to root input name_DS2_DEN</li>
<li>To perform a restart, in a multidataset mode, use ndtset 2 and jdtset 2 3 (that is 2 datasets, numbered 2 and 3) </li>
<li>In the dataset 2, get the density you just copied (getden2 -1), perform a non selfconsistent calculation and print the wave function (prtwf2 1)</li>
<li>In the dataset 3, get  the previous wf(getwfk3 -1), and continue the calculation </li>
<li>This complicated procedure is due to the fact that reading the density is only allowed for a non sc calculation,
and also for a dataset different of 0 or the previous one, the option we choose here.</li>

</ul>

Please note that in the case of PAW (<a href="varint.html#usepaw">usepaw</a>=1) calculations, the _DEN density output
is not the full physical electron density. If what is wanted is the full physical electron density, say
for post-processing with <A href="../users/aim_help.html">AIM</A> or visualization, prtden > 1 will produce physical electron
density or other interesting quantites (see below). Nevertheless, even in the PAW case, when chaining together
calculations where the density from one calculation is to be used in a subsequent calculation, it is necessary
to use the _DEN files and <b>not</b> one of the other files produced with prtden > 1, i.e. _PAWDEN, ATMDEN_xxx or else.
Note that the usual _DEN file is allways generated as soon as prtden >= 1.

Options 2 to 6 for prtden are relevant only for <a href="varint.html#usepaw">usepaw</a>=1 and control the output of the
full electron density in the PAW case :<br><br>
<div style="margin-left: 10px; ">
<b> prtden=2</b> causes generation of a file _PAWDEN that contains the bulk <b>valence</b> charge density together with the PAW on-site contributions, and has the same format as the other density files.<br>
<b> prtden=3</b> causes generation of a file _PAWDEN that contains the bulk <b>full</b> charge density (valence+core)<br>
<b> prtden=4</b> causes generation of three files _ATMDEN_CORE, _ATMDEN_VAL and _ATMDEN_FULL which respectively contain the core, valence and full atomic protodensity (the density of the individual component atoms in vacuum superposed at the bulk atomic positions). This can be used to generate various visualizations of the bonding density.<br>
<b> prtden=5</b> options 2 and 4 taken together.<br>
<b> prtden=6</b> options 3 and 4 taken together.<br>
<b> prtden=7</b> causes the generation of all the individual contributions to the bulk <b>valence</b> charge density : n_tilde-n_hat (_N_TILDE), n_onsite (_N_ONE) and n_tilde_onsite (_NT_ONE). This is for diagnosis purposes only. <br>
</div><br>
Options 3 to 6 currently require the user to supply the atomic core and valence density in external files in the working directory. The files must be named properly; for example, the files for an atom of type 1 should be named: "core_density_atom_type1.dat" and "valence_density_atom_type1.dat". The file should be a text file, where the first line is assumed to be a comment, and the subsequent lines contain two values each, where the first one is a radial coordinate and the second the value of the density n(r). Please note that it is n(r) which should be supplied, <b>not</b> n(r)/r^2. The first coordinate point must be the origin, i.e. <b><i>r = 0</i></b>. The atomic densities are spherically averaged, so assumed to be completely spherically symmetric, even for open shells.<br><br>
NOTE: in the PAW case, <b>DO NOT</b> use _PAWDEN or _ATMDEN_xxx files produced by prtden > 1 to chain the density output from one
calculation as the input to another, use the _DEN file for that.

"""
},
'prtdensph': {
'definition': "PRinT integral of DENsity inside atomic SPHeres",
'section': "vargs",
'category': "",
'vartype': "integer parameter",
'default': """0 except when antiferro-magnetism is activated (<a href="varbas.html#nsppol">nsppol</a>=1, <a href="vargs.html#nspden">nspden</a>=2).""",
'text': """When this flag is activated, values of integral(s) of total density inside sphere(s) around each atom are printed in output file (for each spin component).
Spheres around atoms are defined by a radius given by <a href="vargs.html#ratsph">ratsph</a> keyword.<br>
Note: integral of density inside a sphere around an atom
can be used to determine a rough approximation of the local magnetic moment;
this is particulary useful for antiferromagnetic systems.
<br>
The algorithm to compute this integral is particularly primitive : the points on the FFT grids, belonging
to the interior of the sphere are determined, and the value of the functions on these points are summed,
taking into account a fixed volume attributed to each point.
In particular, the integral as a function of the radius will be a constant, except when
a new point enters the sphere, in which case a sudden jump occurs.
However, since the purpose of this output is to get a rough idea of the repartition of the density,
this is not a real problem. If you are interested in a more accurate estimation
of the density within a sphere, you should use the cut3d postprocessor.
"""
},
'prtdipole': {
'definition': "PRinT DIPOLE",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>""",
'vartype': "integer",
'default': "0 ",
'text': """Print out dipole of unit cell, calculated in real space for the primitive cell only. Under development.
"""
},
'prtdos': {
'definition': "PRinT the Density Of States ",
'section': "varfil",
'category': " ",
'vartype': "integer parameter ",
'default': "0.",
'text': """Provide output of Density of States if set to 1, 2 or 3.
Can either use a smearing technique (<b>prtdos</b>=1),
or the tetrahedron method (<b>prtdos</b>=2).
If <b>prtdos</b>=3, provide output of Local Density of States inside a sphere centered on an atom,
as well as the angular-momentum projected DOS, in the same sphere. The resolution of the linear grid of energies
for which the DOS is computed can be tuned thanks to <a href="vargs.html#dosdeltae">dosdeltae</a>.
<p>
If <b>prtdos</b>=1, the smeared density of states is obtained
from the eigenvalues, properly weighted at each k point
using <a href="varbas.html#wtk">wtk</a>, and smeared according to <a href="varbas.html#occopt">occopt</a>
and <a href="vargs.html#tsmear">tsmear</a>. All levels that are present in the calculation
are taken into account (occupied and unoccupied).
Note that <a href="varbas.html#occopt">occopt</a> must be between 3 and 7 .
<br> In order to compute the DOS of an insulator with <b>prtdos</b>=1, compute its density thanks to
a self-consistent calculation (with a non-metallic <a href="varbas.html#occopt">occopt</a>
value, 0, 1 or 2), then use <b>prtdos</b>=1, together
with <a href="varbas.html#iscf">iscf</a>=-3, and a metallic <a href="varbas.html#occopt">occopt</a>,
between 3 and 7, providing the needed smearing.
If <b>prtdos</b>=1, the name of the DOS file is the root name for the output
files, followed by "_DOS" .
<p>
If <b>prtdos</b>=2, the DOS is computed using the tetrahedron method.
As in the case of <b>prtdos</b>=1, all levels that are present in the calculation
are taken into account (occupied and unoccupied). In this case, the
k-points must have been defined using the input variable <a href="varbas.html#ngkpt">ngkpt</a>
or the input variable <a href="vargs.html#kptrlatt">kptrlatt</a>. There must be at least
two non-equivalent points in the Irreducible Brillouin Zone to use  <b>prtdos</b>=2.
It is strongly advised to use a non-shifted k-point grid
(<a href="varbas.html#shiftk">shiftk</a> 0 0 0): such grids contain naturally
more extremal points (band minima and maxima at Gamma or at the zone-boundaries) than
shifted grids, and lead to more non-equivalent points than shifted grids, for the same grid spacing.

There is no need to take care of
the <a href="varbas.html#occopt">occopt</a> or <a href="vargs.html#tsmear">tsmear</a> input variables,
and there is no subtlety to be taken into account for insulators. The computation
can be done in the self-consistent case as well as in the non-self-consistent case,
using <a href="varbas.html#iscf">iscf</a>=-3. This allows to refine the DOS at fixed
starting density.
<br>In that case, if <a href="varrlx.html#ionmov">ionmov</a>==0, the name of the potential file will be
the root output name, followed by _DOS (like in the <b>prtdos</b>=1 case).
<br>However, if <a href="varrlx.html#ionmov">ionmov</a>==1 or 2, potential files will be output
at each time step, with the name being made of
<ul>
<li> the root output name, </li>
<li> followed by _TIMx , where
x is related to the timestep (see later)</li>
<li> then followed by _DOS.</li>
</ul>
<p>
If <b>prtdos</b>=3, the same tetrahedron method as for <b>prtdos</b>=2 is used, but
the DOS inside a sphere centered on some atom is delivered, as well as the angular-momentum
projected (l=0,1,2,3,4) DOS in the same sphere. The preparation of this
case, the parameters under which the computation is to be done, and the file
denomination is similar
to the <b>prtdos</b>=2 case. However, three additional input variables might be provided,
describing the atoms that are the center of the sphere (input variables
<a href="vargs.html#natsph">natsph</a> and  <a href="vargs.html#iatsph">iatsph</a>), as well as the radius of this
sphere (input variable <a href="vargs.html#ratsph">ratsph</a>).
<br>
In case of PAW, <a href="vargs.html#ratsph">ratsph</a> radius has to be greater or equal to
largest PAW radius of the atom types considered (which is read from the PAW atomic data file; see rc_sph or r_paw).
Additional printing and/or approximations in PAW mode can be controlled with <a href="varpaw.html#pawprtdos">pawprtdos</a> keyword
(in particular,<a href="varpaw.html#pawprtdos">pawprtdos</a>=2 can be used to compute quickly a very good approximation of the DOS).
<br><br>Note 1: when <b>prtdos</b>=3, it is possible to ouptut m-decomposed LDOS in _DOS file; simply use <a href="varfil.html#prtdosm">prtdosm</a> keyword.
<br>Note 2: the integrated total DOS in spheres around atoms can be obtained when
<a href="vargs.html#prtdensph">prtdensph</a> flag is activated.
It can be compared to the integrated DOS provided in _DOS file when <b>prtdos</b>=3.
"""
},
'prtdosm': {
'definition': "PRinT the Density Of States with M decomposition ",
'section': "varfil",
'category': " ",
'vartype': "integer parameter ",
'default': "0.",
'text': """Relevant only when <a href="varfil.html#prtdos">prtdos</a>=3.
<br>If set to 1, the m-decomposed LDOS is delivered in DOS file.
<br>Note that <b>prtdosm</b> computes the M-resolved partial dos for complex spherical harmonics,giving e.g.
DOS(L,M) == DOS(L,-M) (without spin-orbit). In the contrary, the LDA+U occupation matrix,
see <a href="varpaw.html#dmatpawu">dmatpawu</a> is in the real spherical harmonics basis.
"""
},
'prtefg': {
'definition': "PRint Electric Field Gradient ",
'section': "varpaw",
'category': "",
'vartype': "integer parameter ",
'default': "0 ",
'text': """Relevant only when
<a href="varint.html#usepaw">usepaw</a>=1.
<ul>
<li> If set to 1, calculate the electric field gradient at each atomic site in the unit cell. Using
this option requires <a href="varpaw.html#quadmom">quadmom</a> to be set as well.
Values written to main output file (search for Electric Field Gradient). If prtefg=1,
only the quadrupole coupling in MHz and asymmetry are reported. If prtefg=2, the full electric field
gradient tensors in atomic units are also given, showing separate contributions from the valence
electrons, the ion cores, and the PAW reconstruction. If prtefg=3, then in addition to the prtefg=2 output,
the EFG's are computed using an ionic point charge model. This is useful for comparing the accurate PAW-based
results to those of simple ion-only models. Use of prtefg=3 requires that the variable
<a href="varpaw.html#ptcharge">ptcharge</a> be set as well.
<br> The option prtefg is compatible with spin polarized calculations
(see <a href="vargs.html#nspden">nspden</a>) and also LDA+U (see <a href="varpaw.html#usepawu">usepawu</a>).

</li>
</ul>

"""
},
'prteig': {
'definition': "PRinT EIGenenergies ",
'section': "varfil",
'category': " ",
'vartype': "integer parameter ",
'default': """1 (was 0 before v5.3), except when <a href="varrlx.html#nimage">nimage</a>>1""",
'text': """If set to 1, a file *_EIG, containing the
k-points and one-electron eigenvalues is printed.
"""
},
'prtelf': {
'definition': "PRinT Electron Localization Function (ELF) ",
'section': "varfil",
'category': " ",
'vartype': "integer parameter ",
'default': "0 ",
'text': """If set to 1 or a larger value, provide output of ELF
in real space elf(r). This is a dimensionless quantity bounded between 0 and 1.
<br>The name of the ELF file will be the root output name, followed by _ELF.
<br>Like a _DEN file, it can be analyzed by cut3d. However unlike densities, in case of spin polarized
calculations, the spin down component can not be obtained by substracting the spin up component to
the total ELF. Hence when spin polarized calculations are performed the code produces also output files with
_ELF_UP and _ELF_DOWN extensions. (For technical reasons these files contain also two components but the second is zero.
So to perform analysis of _ELF_UP and _ELF_DOWN files with cut3d you have to answer "ispden= 0 ==> Total density"
when cut3d ask you which ispden to choose. Also remember that spin down component can not be obtained by using cut3d on the _ELF file.
Sorry for the inconvenience, this will be fixed in the next release.)
<br>
ELF is not yet implemented in non collinear spin case.
<br>
If prtelf is set to 2, in the case of spin polarized calculation, the total ELF is computed from an alternative approach which should better take into account the existence of spin dependent densities (see the documentation in /doc/theory/ELF of your ABINIT repository)
<br>
<br>Please note that ELF is <b>not</b> yet implemented in the case of PAW (<a href="varint.html#usepaw">usepaw</a>=1) calculations.
"""
},
'prtfc': {
'definition': "PRinT Fermi Contact term ",
'section': "varpaw",
'category': "",
'vartype': "integer parameter ",
'default': "0 ",
'text': """Relevant only when
<a href="varint.html#usepaw">usepaw</a>=1.
<ul>
<li> If set to 1,  print the Fermi contact interaction at each nuclear site, that is,
the electron density at each site. The result appears in the main output file (search for FC).
Note that this calculation is different than what is done by cut3d, because it also computes the
PAW on-site corrections in addition to the contribution from the valence pseudo-wavefunctions.
</li>
</ul>

"""
},
'prtfsurf': {
'definition': "PRinT Fermi SURFace file ",
'section': "varfil",
'category': " ",
'vartype': "integer parameter ",
'default': "0.",
'text': """If set to 1, provide Fermi surface file in the BXSF format (Xcrysden)
If <b>prtfsurf</b>=1, a _BXSF file readable by <a href="http://www.xcrysden.org">XCrySDen</a> will
be produced at the end of the calculation. The file contains information on the band structure
of the system and can be used to visualize the Fermi surface or any other energy isosurface.
<b>prtfsurf</b>=1 is compatible only with SCF calculations
(<a href="varbas.html#iscf">iscf</a> > 1) or NSCF runs in which the
occupation factors and Fermi level are recalculated once convergence is achieved (<a href="varbas.html#iscf">iscf</a> = -3).
The two methods should produce the same Fermi surface provided that the k-meshes are sufficiently dense.
The k-mesh used for the sampling of the Fermi surface can be specified using the standard
variables <a href="varbas.html#ngkpt">ngkpt</a>, (<a href="varbas.html#shiftk">shiftk</a>, and <a href="varbas.html#nshiftk">nshiftk</a>.
Note, however, that the mesh must be homogeneous and centered on gamma (multiple shifts are not supported by Xcrysden)
"""
},
'prtgden': {
'definition': "PRinT the Gradient of electron DENsity ",
'section': "varfil",
'category': " ",
'vartype': "integer parameter ",
'default': "0. ",
'text': """If set to 1  or a larger value, provide output of gradient of electron density
in real space grho(r), in units of Bohr^-(5/2).
<br>The names of the gradient of electron density files will be the root output name,
followed by _GDEN1, _GDEN2, GDEN3 for each principal direction (indeed it is a vector).
<br>Like a _DEN file, it can be analyzed by cut3d. The file structure of the unformatted output file is
described below, see section 6).
"""
},
'prtgeo': {
'definition': "PRinT the GEOmetry analysis  ",
'section': "varfil",
'category': " ",
'vartype': "integer parameter ",
'default': "0.",
'text': """If set to 1  or a larger value, provide output of geometrical analysis
(bond lengths and bond angles). The value
of <b>prtgeo</b> is taken by the code to be the
maximum coordination number of atoms in the system.
<br>It will deduce a maximum number of "nearest" and "next-nearest"
neighbors accordingly , and compute corresponding bond lengths.
<br>It will compute bond angles for the "nearest" neighbours only.
<br>If <a href="varrlx.html#ionmov">ionmov</a>==0, the name of the file will be
the root output name, followed by _GEO .
<br>If <a href="varrlx.html#ionmov">ionmov</a>==1 or 2, one file will be output
at each time step, with the name being made of
<ul>
<li> the root output name,</li>
<li> followed by _TIMx , where
x is related to the timestep (see later)</li>
<li> then followed by _GEO</li>
</ul>
The content of the file should be rather self-explanatory.
<br>No output is provided by <b>prtgeo</b> is lower than or equal to 0.
<br>If <b>prtgeo</b>>0, the maximum number of atoms (<a href="varbas.html#natom">natom</a>) is 9999."""
},
'prtgkk': {
'definition': "PRinT the GKK matrix elements file  ",
'section': "varfil",
'category': " ",
'vartype': "integer parameter ",
'default': "0.",
'text': """If set to 1  or a larger value, provide output of electron-phonon
"gkk" matrix elements, for further treatment by mrggkk utility or anaddb utility.
Additional information on electron-phonon treatment in ABINIT is given in ~abinit/doc/users/elphon_manual.ps .
"""
},
'prtkden': {
'definition': "PRinT the Kinetic energy DENsity ",
'section': "varfil",
'category': " ",
'vartype': "integer parameter ",
'default': "0. ",
'text': """If set to 1  or a larger value , provide output of kinetic energy density
in real space tau(r), in units of Bohr^-5.
<br>The name of the kinetic energy density file will be the root output name, followed by _KDEN.
<br>Like a _DEN file, it can be analyzed by cut3d. The file structure of the unformatted output file is
described below, see section 6).
<br>Note that the computation of the kinetic energy density must be activate, thanks to the input variable <a href="vargs.html#usekden">usekden</a>.
<br>Please note that kinetic energy density is <b>not</b> yet implemented in the case of PAW (<a href="varint.html#usepaw">usepaw</a>=1) calculations.
"""
},
'prtkpt': {
'definition': "PRinT the K-PoinTs sets  ",
'section': "varfil",
'category': " ",
'vartype': "integer parameter ",
'default': "0.",
'text': """If set /= 0 , proceeds to a detailed analysis
of different k point grids. Works only if
<a href="varbas.html#kptopt">kptopt</a> is positive, and neither
<a href="vargs.html#kptrlatt">kptrlatt</a>
nor
<a href="varbas.html#ngkpt">ngkpt</a> are defined.
ABINIT will stop after this analysis.

<p>Different sets of k point grids are defined,
with common values of <a href="varbas.html#shiftk">shiftk</a>.
In each set, ABINIT increases the length of vectors of the
supercell (see <a href="vargs.html#kptrlatt">kptrlatt</a>) by integer
steps. The different sets are labelled by "iset".
For each k point grid,  <a href="vargs.html#kptrlen">kptrlen</a>
and  <a href="varbas.html#nkpt">nkpt</a> are computed (the latter always
invoking  <a href="varbas.html#kptopt">kptopt</a>=1, that is, full use of
symmetries). A series is finished when the computed
<a href="vargs.html#kptrlen">kptrlen</a> is twice larger than the
input variable <a href="vargs.html#kptrlen">kptrlen</a>.
After the examination of the different sets,
ABINIT summarizes, for each <a href="varbas.html#nkpt">nkpt</a>, the
best possible grid, that is, the one with the
largest computed <a href="vargs.html#kptrlen">kptrlen</a>.

<p>Note that this analysis is also performed when
<b>prtkpt</b>=0, as soon as neither <a href="vargs.html#kptrlatt">kptrlatt</a>
nor
<a href="varbas.html#ngkpt">ngkpt</a> are defined. But, in this case,
no analysis report is given, and the code selects the grid
with the smaller <a href="varbas.html#ngkpt">ngkpt</a>
for the desired <a href="vargs.html#kptrlen">kptrlen</a>. However,
this analysis takes some times (well sometimes,
it is only a few seconds - it depends on the value of
the input <a href="vargs.html#kptrlen">kptrlen</a>), and it is better
to examine the full analysis for a given cell and set of symmetries,
<a href="varbas.html#shiftk">shiftk</a> for all the production runs.
"""
},
'prtlden': {
'definition': "PRinT the Laplacian of electron DENsity ",
'section': "varfil",
'category': " ",
'vartype': "integer parameter ",
'default': "0. ",
'text': """If set to 1  or a larger value, provide output of Laplacian of electron density
in real space grho(r), in units of Bohr^-(7/2).
<br>The name of the Laplacian of electron density file will be the root output name,
followed by _LDEN.
<br>Like a _DEN file, it can be analyzed by cut3d. The file structure of the unformatted output file is
described below, see section 6).
"""
},
'prtnabla': {
'definition': "PRint NABLA ",
'section': "varpaw",
'category': "",
'vartype': "integer parameter ",
'default': "0 ",
'text': """Relevant only when
<a href="varint.html#usepaw">usepaw</a>=1.
<ul>
<li> If set to 1, calculate the matrix elements &lt;Psi_n|-inabla|Psi_m&gt; and write it in file _OPT to be read by the code conducti. </li>
</ul>

"""
},
'prtnest': {
'definition': "PRinT NESTing function",
'section': "vardev",
'category': " DEVELOP",
'vartype': "integer flag",
'default': "0 ",
'text': """If set to 1, the nesting function for the k-point grid is printed. For the moment the path in q space for the nesting function is fixed, but will become an input as well.
"""
},
'prtposcar': {
'definition': "PRinT POSCAR file",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>""",
'vartype': "integer",
'default': "0 ",
'text': """Print out VASP-style POSCAR and FORCES files, for use with PHON or frophon codes for frozen phonon calculations.
See the associated script in ~abinit/extras/post_processing/phondisp2abi.py for further details on interfacing
with PHON, PHONOPY, etc...
"""
},
'prtpot': {
'definition': "PRinT the iotal (kohn-sham)POTential ",
'section': "varfil",
'category': " ",
'vartype': "integer parameter ",
'default': "0.",
'text': """If set >=1 , provide output of different
potentials.
<br>For <b>prtpot</b>, output the total (Kohn-Sham) potential,
sum of local pseudo-potential, Hartree potential, and xc potential.
<br>For <b>prtvha</b>, output the Hartree potential.
<br>For <b>prtvhxc</b>, output the
sum of Hartree potential and xc potential.
<br>For <b>prtvxc</b>, output the exchange-correlation potential.

<p>If <a href="varrlx.html#ionmov">ionmov</a>==0, the name of the potential file will be
the root output name, followed by _POT, _VHA, _VHXC, or _VXC .
<br>If <a href="varrlx.html#ionmov">ionmov</a>==1 or 2, potential files will be output
at each time step, with the name being made of
<ul>
<li> the root output name, </li>
<li> followed by _TIMx , where
x is related to the timestep (see later)</li>
<li> then followed by _POT, _VHA, _VHXC, or _VXC.</li>
</ul>
The file structure of this unformatted output file is
described in <a href="../users/abinit_help.html#localpotfile">section 6.6</a> of abinit_help.
No output is provided by a negative value of these variables.
"""
},
'prtspcur': {
'definition': "PRinT the SPin CURrent density ",
'section': "varfil",
'category': " ",
'vartype': "integer parameter ",
'default': "0.",
'text': """If set to 1  or a larger value, provide output of the current density of
different direction spins (x,y,z) in the whole unit cell. Should require spinorial wave functions
<a href="vargs.html#nspinor"><b>nspinor</b></a> = 2. Experimental: this does not work yet.
"""
},
'prtstm': {
'definition': "PRinT the STM density ",
'section': "varfil",
'category': " ",
'vartype': "integer parameter ",
'default': "0.",
'text': """If set to 1  or a larger value, provide output of the electron density
in real space rho(r), made only from the electrons
close to the Fermi energy, in a range of energy (positive or negative), determined
by the (positive or negative, but non-zero) value of the STM bias <a href="vargs.html#stmbias">stmbias</a>.
<br>This is a very approximate way to obtain STM profiles : one can choose an equidensity surface,
and consider that the STM tip will follow this surface. Such equidensity surface might be determined
with the help of Cut3D, and further post-processing of it (to be implemented). The big approximations
of this technique are : neglect of the finite size of the tip, and
position-independent transfer matrix elements
between the tip and the surface.
<br>The charge density is provided in units of electrons/Bohr^3.
The name of the STM density file will be the root output name, followed by _STM .
Like a _DEN file, it can be analyzed by cut3d.
The file structure of this unformatted output file is
described in <a href="../users/abinit_help.html#densoutputfile">section 6.5</a> of abinit_help.
<br>For the STM charge density to be generated, one must give, as an input file, the
converged wavefunctions obtained from a previous run, at exactly the same k-points and cut-off energy,
self-consistently determined, using the occupation numbers from <a href="varbas.html#occopt">occopt</a>=7.
<br>In the run with positive <b>prtstm</b>, one has to use :
<ul>
<li>positive <a href="varbas.html#iscf">iscf</a> </li>
<li> <a href="varbas.html#occopt">occopt</a>=7, with specification of <a href="vargs.html#tsmear">tsmear</a> </li>
<li> <a href="varbas.html#nstep">nstep</a>=1 </li>
<li> the <a href="varbas.html#tolwfr">tolwfr</a> convergence criterion</li>
<li> <a href="varrlx.html#ionmov">ionmov</a>=0 (this is the default value)</li>
<li> <a href="vargs.html#optdriver">optdriver</a>=0 (this is the default value)</li>
</ul>
<br>Note that you might have to adjust the value of <a href="varbas.html#nband">nband</a>
as well, for the treatment of unoccupied states, because the automatic determination
of <a href="varbas.html#nband">nband</a> will often not include enough unoccupied
states.
<br>When <b>prtstm</b> is non-zero, the stress tensor is set to zero.
<br>No output of _STM file is provided by <b>prtstm</b> lower or equal to 0.
<br>No other printing variables for density or potentials should be activated (e.g. <a href="varfil.html#prtden">prtden</a> has to be set to zero).
"""
},
'prtsuscep': {
'definition': "PRinT the SUSCEPtibility file (the irreducible polarizability)",
'section': "varfil",
'category': " ",
'vartype': "integer parameter ",
'default': "1.",
'text': """If set to 0, no _SUSC file will be produced after the screening calculation,
only the _SCR file will be output.
"""
},
'prtvha': {
'definition': "PRinT V_HArtree ",
'section': "varfil",
'category': " ",
'vartype': "integer parameter ",
'default': "0.",
'text': """If set >=1 , provide output of different
potentials.
<br>For <b>prtpot</b>, output the total (Kohn-Sham) potential,
sum of local pseudo-potential, Hartree potential, and xc potential.
<br>For <b>prtvha</b>, output the Hartree potential.
<br>For <b>prtvhxc</b>, output the
sum of Hartree potential and xc potential.
<br>For <b>prtvxc</b>, output the exchange-correlation potential.

<p>If <a href="varrlx.html#ionmov">ionmov</a>==0, the name of the potential file will be
the root output name, followed by _POT, _VHA, _VHXC, or _VXC .
<br>If <a href="varrlx.html#ionmov">ionmov</a>==1 or 2, potential files will be output
at each time step, with the name being made of
<ul>
<li> the root output name, </li>
<li> followed by _TIMx , where
x is related to the timestep (see later)</li>
<li> then followed by _POT, _VHA, _VHXC, or _VXC.</li>
</ul>
The file structure of this unformatted output file is
described in <a href="../users/abinit_help.html#localpotfile">section 6.6</a> of abinit_help.
No output is provided by a negative value of these variables.
"""
},
'prtvhxc': {
'definition': "PRinT V_(Hartree+XC) ",
'section': "varfil",
'category': " ",
'vartype': "integer parameter ",
'default': "0.",
'text': """If set >=1 , provide output of different
potentials.
<br>For <b>prtpot</b>, output the total (Kohn-Sham) potential,
sum of local pseudo-potential, Hartree potential, and xc potential.
<br>For <b>prtvha</b>, output the Hartree potential.
<br>For <b>prtvhxc</b>, output the
sum of Hartree potential and xc potential.
<br>For <b>prtvxc</b>, output the exchange-correlation potential.

<p>If <a href="varrlx.html#ionmov">ionmov</a>==0, the name of the potential file will be
the root output name, followed by _POT, _VHA, _VHXC, or _VXC .
<br>If <a href="varrlx.html#ionmov">ionmov</a>==1 or 2, potential files will be output
at each time step, with the name being made of
<ul>
<li> the root output name, </li>
<li> followed by _TIMx , where
x is related to the timestep (see later)</li>
<li> then followed by _POT, _VHA, _VHXC, or _VXC.</li>
</ul>
The file structure of this unformatted output file is
described in <a href="../users/abinit_help.html#localpotfile">section 6.6</a> of abinit_help.
No output is provided by a negative value of these variables.
"""
},
'prtvol': {
'definition': "PRinT VOLume ",
'section': "varfil",
'category': " ",
'vartype': "integer parameter ",
'default': "0.",
'text': """Control the volume of printed output.
<br>Standard choice is 0. Positive values print more in the
output and log files, while negative values are
for debugging (or preprocessing only), and cause
the code to stop at some point.
<ul>
<li>0 => there is a limit on the number of k-points
for which related information will be written. This
limit is presently 50. Due to some subtlety, if for <b>some</b> dataset
<b>prtvol</b> is non-zero, the limit for input and output echoes
cannot be enforced, so it is like if <b>prtvol</b>=1 for
<b>all</b> the datasets for which <b>prtvol</b> was set to 0.
Also, if <a href="varbas.html#iscf">iscf</a>=-2 and
<a href="varbas.html#kptopt">kptopt</a>&lt;=0,
the eigenvalues for all the k points are printed anyway.</li>
<li>1 => there is no such limit for the input
and output echoes, in the main output file. Also, if <a href="varbas.html#iscf">iscf</a>=-2 and
<a href="varbas.html#kptopt">kptopt</a>&lt;=0, the eigenvalues for all the k points are printed anyway.</li>
<li>2 => there is no such limit in the whole main output file.</li>
<li>3 => there is no such limit in both output and log files.</li>
<li>10 => no limit on the number of k points, and moreover,
the eigenvalues are printed for every SCF iteration,
as well as other additions (to be specified in the future...)</li>
</ul>
Debugging options :
<ul>
<li>= -1 => stop in abinit (main program), before call gstate.
Useful to see the effect of the preprocessing of
input variables (memory needed, effect of symmetries,
k points ...) without going further. Run very fast,
on the order of the second.
<li>=-2 => same as -1, except that print only the first dataset.
All the non default input variables associated to all datasets
are printed in the output file, but only for the first dataset.
Also all the input variables are written in the NetCDF file
"OUT.nc", even if the value is the default.
</li>
<li>= -3 => stop in gstate, before call scfcv, move or brdmin.
Useful to debug pseudopotentials</li>
<li>= -4 => stop in move, after completion of all loops</li>
<li>= -5 => stop in brdmin, after completion of all loops</li>
<li>= -6 => stop in scfcv, after completion of all loops</li>
<li>= -7 => stop in vtorho, after the first rho is obtained</li>
<li>= -8 => stop in vtowfk, after the first k point is treated</li>
<li>= -9 => stop in cgwf, after the first wf is optimized </li>
<li>= -10 => stop in getghc, after the Hamiltonian is applied once</li>
</ul>
This debugging feature is not yet activated in the RF routines.
Note that <a href="vardev.html#fftalg">fftalg</a> offers
another option for debugging."""
},
'prtvolimg': {
'definition': "PRinT VOLume for IMaGes ",
'section': "varfil",
'category': " ",
'vartype': "integer parameter ",
'default': "0.",
'text': """Control the volume of printed output
when an algorithm using images of the cell is used (<a href="varrlx.html#nimage">nimage</a>>1).<br>
When such an algorithm is activated, the printing volume (in output file)
can be large and difficult to read.<br>
Using <b>prtvolimg=1</b>, the printing volume, for each image,
is reduced to unit cell, atomic positions, total energy, forces, stresses, velocities
and convergence residuals.<br>
Using <b>prtvolimg=2</b>, the printing volume, for each image,
is reduced to total energy and convergence residuals only.
"""
},
'prtvxc': {
'definition': "PRinT V_XC ",
'section': "varfil",
'category': " ",
'vartype': "integer parameter ",
'default': "0.",
'text': """If set >=1 , provide output of different
potentials.
<br>For <b>prtpot</b>, output the total (Kohn-Sham) potential,
sum of local pseudo-potential, Hartree potential, and xc potential.
<br>For <b>prtvha</b>, output the Hartree potential.
<br>For <b>prtvhxc</b>, output the
sum of Hartree potential and xc potential.
<br>For <b>prtvxc</b>, output the exchange-correlation potential.

<p>If <a href="varrlx.html#ionmov">ionmov</a>==0, the name of the potential file will be
the root output name, followed by _POT, _VHA, _VHXC, or _VXC .
<br>If <a href="varrlx.html#ionmov">ionmov</a>==1 or 2, potential files will be output
at each time step, with the name being made of
<ul>
<li> the root output name, </li>
<li> followed by _TIMx , where
x is related to the timestep (see later)</li>
<li> then followed by _POT, _VHA, _VHXC, or _VXC.</li>
</ul>
The file structure of this unformatted output file is
described in <a href="../users/abinit_help.html#localpotfile">section 6.6</a> of abinit_help.
No output is provided by a negative value of these variables.
"""
},
'prtwant': {
'definition': "PRinT WANT file ",
'section': "varfil",
'category': " ",
'vartype': "integer parameter ",
'default': "0.",
'text': """Flag used to indicate that either the Wannier90 or the WanT interfaces
will be used.
<ul>
<li> =1 => Use the <b>ABINIT- WanT </b> interface.<p>
Provide an output file that can be used by the WanT postprocessing
program (see http://www.wannier-transport.org).
The value of the prtwant indicates the version of the WanT code that
can read it. Currently only the value <b>prtwant</b>=1 is implemented,
corresponding to WanT version 1.0.1, available since Oct. 22, 2004.
<p>
Notes : Several requirements must be
fulfilled by the wavefunction. Among them, two are mandatory:
<br>
<ul>
<li>1. An uniform grid of k-points, including the GAMMA point must be used.</li>
<li>2. The use of time reversal symmetry is not allowed (istwfk=1) </li>
<li>3. The list of k-points must be ordered, such that the coordinates,
namely three-components vectors has the third index varying the most
rapidly, then the second index, then the first index </li>
</ul>
If these requirement are not fulfilled, the program will stop and an error message is returned.
<p>
As an example of k-point grid in case of systems that have some 3D character (1D systems are easy) :
<pre>
nkpt 8
kpt  0   0   0
0   0   1/2
0   1/2 0
0   1/2 1/2
1/2 0   0
1/2 0   1/2
1/2 1/2 0
1/2 1/2 1/2
istwfk 8*1
</pre>
<p>
Also, in order to use WanT as a postprocessing program for ABINIT you might have to
recompile it with the appropriate flags (see ABINIT makefile). Up to now only
the -convert big-endian was found to be mandatory, for machines with little-endian default choice.
<p>
If set = 2, provide an output file that can be used by the wannier90 postprocessing
program in next versions
<li>=2 => Use the <b>ABINIT- Wannier90 </b> interface.<p>
ABINIT will produce the input files required by Wannier90  and it will
run Wannier90 to produce the Maximally-locallized Wannier functions
(see <a href="http://www.wannier.org">
http://www.wannier.org </a> ).

<p>Notes:

<ul>
<li>The files that are created can also be used by Wannier90 in
stand-alone mode.
<li>
In order to use Wannier90 as a postprocessing program for ABINIT you
might have to recompile it with the appropriate flags (see ABINIT
makefile). You might use ./configure --enable-wannier90
<li>
There are some other variables related to the interface of Wannier90
and ABINIT. See, <a href="varw90.html">VARW90</a>.
</ul>
<p>

<li>=3 => Use the <b>ABINIT- Wannier90 </b> interface after converting the
input wavefunctions to <b>GW quasiparticle</b> wavefunctions.<p>
ABINIT will produce the input files required by Wannier90  and it will
run Wannier90 to produce the Maximally-locallized Wannier functions
(see <a href="http://www.wannier.org">
http://www.wannier.org </a> ).

<p>Additional Notes:

<ul>
<li>An input file of LDA wave functions is required which is completely
consistent with the _KSS file used in the self-consistent GW calculation.
This means that <a href="varfil.html#kssform">kssform</a> 3 must be used
to create the _KSS file and the output _WFK file from the same run must
be used as input here.
<li>Wannier90 requires <a href="varbas.html#nshiftk">nshiftk</a>=1, and
<a href="varbas.html#shiftk">shiftk</a>= 0 0 0 is recommended.  The k-point
set used for the GW calculation, typically the irreducible BZ set created
using <a href="varbas.html#kptopt">kptopt</a>=1, and that for the Abinit-
Wannier90 interface must be consistent.
<li>Full-BZ wavefunctions should be generated in the run calling the
interface by setting <a href="varbas.html#kptopt">kptopt</a>=3,
<a href="varbas.html#iscf">iscf</a>=-2, and
<a href="varbas.html#nstep">nstep</a>=3.  This will simply use symmetry to
transform the input IBZ wavefunctions to the full BZ set, still consistent
with the GW _KSS input.
<li>The final _QPS file created by the self-consistent GW run is required
as input.
<li>Any value of <a href="vargw.html#gwcalctyp">gwcalctyp</a> between
between 20 and 29 should be suitable, so, for example, Hartree-Fock
maximally-localized Wannier functions could be generated setting
<a href="vargw.html#gwcalctyp">gwcalctyp</a>=25.

</ul>


</ul>

"""
},
'prtwf': {
'definition': "PRinT the WaveFunction ",
'section': "varfil",
'category': " ",
'vartype': "integer parameter ",
'default': """1, except when <a href="varrlx.html#nimage">nimage</a>>1""",
'text': """If <b>prtwf</b>=1 , provide output of wavefunction
and eigenvalue file, as described in
<a href="../users/abinit_help.html#wavefctfile">section 6.7</a> of the main abinit help file.
<br>For a standard ground-state calculation, the name of the wavefunction file will be
the root output name, followed by _WFK. If <a href="vargs.html#nqpt">nqpt</a>=1,
the root name will be followed by _WFQ. For response-function calculations,
the root name will be followed by _1WFx, where x is the number of the perturbation.
The dataset information will be added as well, if relevant.
<br>No wavefunction output is provided by <b>prtwf</b>=0.
<p>
<br>If <b>prtwf</b>=2, a file pwfn.data is produced, to be used as input for the
CASINO QMC code.
<br>To produce a wave function suitable for use as a CASINO trial wave
function, certain ABINIT parameters must be set correctly. Primarily,
CASINO (and QMC methods generally) can only take advantage of
time-reversal symmetry, and not the full set of symmetries of the crystal
structure. Therefore, ABINIT must be instructed to generate k-points not
just in the Irreducible Brillouin Zone, but in a full half of the
Brillouin Zone (using time-reversal symmetry to generate the other half).
Additionally, unless instructed otherwise, Abinit avoids the need for
internal storage of many of the coefficients of its wave functions for
k-points that have the property 2k=G_latt, where G_latt is a reciprocal
lattice vector, by making use of the property that
c_k(G)=c^*_k(-G-G_latt). Abinit must be instructed not to do this in order
to output the full set of coefficients for use in CASINO. See the ABINIT
theoretical background documents ABINIT/Infos/Theory/geometry.pdf and
ABINIT/Infos/Theory/1WF.pdf for more information.
<br>
The first of these requirements is met by setting the ABINIT input
variable kptopt to 2 (see ABINIT/Infos/varbas.html#kptopt) and the second
by setting istwfk to 1 for all the k points (see
ABINIT/Infos/vardev.html#istwfk). Since CASINO is typically run with
relatively small numbers of k-points, this is easily done by defining an
array of 1's in the input file.
<br>
For example, for the 8 k-points generated with ngkpt 2 2 2, we add the
following lines to the input file:
<pre>
# Turn off special storage mode for time-reversal k-points
istwfk 1 1 1 1 1 1 1 1
# Use only time reversal symmetry, not full set of symmetries.
kptopt 2
</pre>
Other useful input variables of relevance to the plane waves ABINIT will
produce include ecut, nshiftk, shiftk, nband, occopt, occ, spinat and
nsppol (see relevant input variable documents in ABINIT/Infos/). If ABINIT
is run in multiple dataset mode, the different wave functions for the
various datasets are exported as pwfn1.data, pwfn2.data, ..., pwfnn.data
where the numbers are the contents of the contents of the input array
jdtset (defaults to 1,2,...,ndtset).
<br>
Once the routine is incorporated into the ABINIT package it is anticipated
that there will be an input variable to control whether or not a CASINO
pwfn.data file is written.
<p>
Other issues related to <b>prtwf</b>=2.
<br>
The exporter does not currently work when ABINIT is used in parallel mode
on multiple processors if k-point parallelism is
chosen. ABINIT does not store the full wave function on each processor but
rather splits the k-points between the processors, so no one processor
could write out the whole file. Clearly this could be fixed but we haven't
done it yet. The sort of plane wave DFT calculations usually required to
generate QMC trial wave functions execute very rapidly anyway and will
generally not require a parallel machines. The outqmc routine currently
bails out with an error if this combination of modes is selected - this
will hopefully be fixed later.
<br>
There has not been very extensive testing of less common situations such
as different numbers of bands for different k-points, and more complicated
spin polarized systems, so care should be taken when using the output in
these circumstances.
<br>
If there is any doubt about the output of this routine, the first place to
look is the log file produced by ABINIT: if there are any warnings about
incorrectly normalized orbitals or non-integer occupation numbers there is
probably something set wrong in the input file.
"""
},
'prtxml': {
'definition': "PRinT an XML output",
'section': "varfil",
'category': " ",
'vartype': "0 or 1  ",
'default': "0",
'text': """Create an XML output with common values. The corresponding
DTD is distributed in sources as extras/post_processing/abinitRun.dtd. All the DTD is not
yet implemented and this one is currently restricted to ground computations
(and derivative such as geometry optimisation)."""
},
'ptcharge': {
'definition': "PoinT CHARGEs ",
'section': "varpaw",
'category': "",
'vartype': """real array ptcharge(<a href="varbas.html#ntypat">ntypat</a>) """,
'default': "0",
'text': """Relevant only when
<a href="varint.html#usepaw">usepaw</a>=1, and <a
href="varpaw.html#prtefg">prtefg</a> = 3 or greater.<br>
<ul>
<li> Array of point charges, in atomic units, of the nuclei. In the normal computation of electric field gradients
(see <a href="varpaw.html#prtefg">prtefg</a>) the ionic contribution is calculated from the core charges of the atomic
sites. Thus for example in a PAW data set for oxygen where the core is 1s2, the core charge is +6 (total nuclear charge
minus core electron charge). In point charge models, which are much less accurate than PAW calculations, all atomic sites
are treated as ions with charges determined by their valence states. In such a case oxygen almost always would have a
point charge of -2. The present variable taken together with <a href="varpaw.html#prtefg">prtefg=3</a> performs a full
PAW computation of the electric field gradient and also a simple point charge computation. The user inputs whatever
point charges he/she wishes for each atom type.
</li>
</ul>

"""
},
'ptgroupma': {
'definition': "PoinT GROUP number for the MAgnetic space group",
'section': "vargeo",
'category': """<a href="../users/abinit_help.html#symmetriser">SYMMETRISER</a>, <a href="../users/abinit_help.html#internal">INTERNAL</a> """,
'vartype': "integer parameter ",
'default': "0.",
'text': """This internal variable characterizes a Shubnikov type III magnetic space group (anti-ferromagnetic
space group). The user is advised to consult
"The mathematical theory of symmetry in solids,
Representation theory for point groups and space groups, 1972,
C.J. Bradley and A.P. Cracknell, Clarendon Press, Oxford."
<br>A Shubnikov type III magnetic space group might be defined by its Fedorov space group
(set of all spatial symmetries, irrespective of their magnetic action), and
the halving space group (only the symmetries that do not change the magnetization).
<br>The specification of the halving space group might be done by specifying, for each
point symmetry, the magnetic action. See Table 7.1 of the above-mentioned reference.
Magnetic point groups are numbered from 1 to 58.
<p> Related input variables :
<a href="vargeo.html#spgroup">spgroup</a>,
<a href="vargeo.html#spgroupma">spgroupma</a>,
<a href="vargeo.html#genafm">genafm</a>
"""
},
'qmass': {
'definition': "Q thermostat mass ",
'section': "varrlx",
'category': "",
'vartype': """real array qmass(<a href="varrlx.html#nnos">nnos</a>) """,
'default': "nnos*10.0 ",
'text': """<br>
This temperature control can be used with&nbsp; <a
href="varrlx.html#optcell">optcell</a>==0, 1 (homogeneous cell
deformation) or 2 (full cell deformation)
"""
},
'qprtrb': {
'definition': "Q-wavevector of the PERTurbation ",
'section': "varff",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>  """,
'vartype': "integer array of three values",
'default': "0 0 0.",
'text': """Gives the wavevector,
in units of reciprocal lattice primitive translations,
of a perturbing potential of strength vprtrb.  See vprtrb
for more explanation."""
},
'qpt': {
'definition': "Q PoinT ",
'section': "vargs",
'category': """<a href="../users/abinit_help.html#not_internal">NOT INTERNAL</a> """,
'vartype': "real array of 3 elements  ",
'default': "0 0 0.",
'text': """Only used if <a href="vargs.html#nqpt">nqpt</a>=1.
<p>
Combined with <a href="vargs.html#qptnrm">qptnrm</a>,
define the q vector <a href="varint.html#qptn">qptn</a>(1:3)
in the case <a href="vargs.html#qptopt">qptopt</a>=0.
<p>This input variable is not internal (<a href="varint.html#qptn">qptn</a>(1:3) is used
instead), but is used to echo the value of <a href="varint.html#qptn">qptn</a>(1:3),
with renormalisation factor one.
"""
},
'qptdm': {
'definition': "Q-PoinTs for the Dielectric Matrix",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a> """,
'vartype': """real array qptdm(3,<a href="vargw.html#nqptdm">nqptdm</a>)  """,
'default': "0. 0. 0. (for just one q-point)",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=3, that is, screening calculations,
and only if <a href="vargw.html#nqptdm">nqptdm</a>/=0.
<p>
<b>qptdm</b> contains the set of q-points used in the screening part of ABINIT,
instead of the automatic generation of the q points when <a href="vargw.html#nqptdm">nqptdm</a>=0.
These q points are given in terms of reciprocal space primitive translations (NOT in cartesian coordinates!).
For further explanation, see the input variable <a href="vargw.html#nqptdm">nqptdm</a>.
"""
},
'qptn': {
'definition': "Q-PoinT re-Normalized ",
'section': "varint",
'category': """<a href="../users/abinit_help.html#internal">INTERNAL</a> """,
'vartype': "real array <b>qptn</b>(3) ",
'default': "",
'text': """Only used if <a href="vargs.html#nqpt">nqpt</a>=1.
<br>In ground-state calculation,
the vector <b>qptn</b>(1:3) is added to
each renormalized k point (whatever the value of
<a href="varbas.html#kptopt">kptopt</a> that was used)
to generate the normalized, shifted, set of k-points
<a href="varint.html#kptns">kptns</a>(1:3,1:<b>nkpt</b>).
<br>In response-function calculations,
<b>qptn</b>(1:3)
is the wavevector of the phonon-type calculation.
<br><b>qptn</b>(1:3) can be produced on the basis of
the different methods described in <a href="vargs.html#qptopt">qptopt</a>,
like using <a href="vargs.html#qpt">qpt</a>(1:3)
with renormalisation provided by <a href="vargs.html#qptnrm">qptnrm</a>,
or using the other possibilities defined by
<a href="vargs.html#iqpt">iqpt</a>,
<a href="vargs.html#ngqpt">ngqpt</a>,
<a href="vargs.html#nshiftq">nshiftq</a>,
<a href="vargs.html#qptrlatt">qptrlatt</a>,
<a href="vargs.html#shiftq">shiftq</a>,
<br>For insulators, there is no restriction on the
q-points to be used for the perturbations. By contrast,
for metals, for the time being, it is advised to take
q points for which the k and k+q grids are the same
(when the periodicity in reciprocal space is taken
into account). <br>Tests remain to be done to see whether
other q points might be allowed (perhaps with some
modification of the code)."""
},
'qptnrm': {
'definition': "Q PoinTs NoRMalization ",
'section': "vargs",
'category': """<a href="../users/abinit_help.html#not_internal">NOT INTERNAL</a> """,
'vartype': "real parameter  ",
'default': "1.0  ",
'text': """Only used if <a href="vargs.html#nqpt">nqpt</a>=1 and <a href="vargs.html#qptopt">qptopt</a>=0
<p>
Provides re-normalization
of <a href="vargs.html#qpt">qpt</a>.
Must be positive, non-zero.
The actual q vector (renormalized) is
<a href="varint.html#qptn">qptn</a>(1:3)=
<a href="vargs.html#qpt">qpt</a>(1:3)/<b>qptnrm</b>.
"""
},
'qptopt': {
'definition': "QPoinTs OPTion ",
'section': "vargs",
'category': """<a href="../users/abinit_help.html#not_internal">NOT INTERNAL</a> """,
'vartype': "integer parameter ",
'default': "0 .",
'text': """Only used if <a href="vargs.html#nqpt">nqpt</a>=1.
<p>
Controls the set up to generate the Q point
<a href="varint.html#qptn">qptn</a>(1:3)
to be used for the specific dataset,
either as a shift of k-point grid in ground-state calculations,
or as a stand-alone phonon wavevector.
<p>There are two basic techniques to generate the Q point : either
by specifying it directly, possibly with a renormalisation factor (<b>qptopt</b>=0),
or extracting it from a grid a Q points (<b>qptopt</b>=1 to 4), using the
index <a href="vargs.html#iqpt">iqpt</a>.
At variance with the similar generation of k points, only ONE q point
can be used per dataset.
<p>With <b>qptopt</b>=1 to 4, rely on <a href="vargs.html#ngqpt">ngqpt</a> or
<a href="vargs.html#qptrlatt">qptrlatt</a>, as well as on
<a href="vargs.html#nshiftq">nshiftq</a> and
<a href="vargs.html#shiftq">shiftq</a> to set up a q point
grid, from which the q point with number <a href="vargs.html#iqpt">iqpt</a> will be selected.
The values <b>qptopt</b>=1 to 4 differ by the treatment of symmetries. Note
that the symmetries are
recomputed starting from the values of <a href="varbas.html#rprimd">rprimd</a>
<a href="varbas.html#xred">xred</a> and
<a href="vargs.html#spinat">spinat</a>. So, the explicit value of
<a href="varbas.html#symrel">symrel</a> are not used.
This is to allow doing calculations with <a href="varbas.html#nsym">nsym</a>=1,
sometimes needed for T-dependent electronic structure, still decreasing
the number of q points in the case <b>qptopt</b>=1 or <b>qptopt</b>=3.

<ul>
<li>0=> read directly <a href="vargs.html#qpt">qpt</a>,
and its (eventual) renormalisation factor <a href="vargs.html#qptnrm">qptnrm</a>.
<li>1=>
Take fully into account the symmetry to generate the
grid of q points in the Irreducible Brillouin Zone only.<br>
(This is the usual mode for RF calculations)</li>
<li>2=>
Take into account only the time-reversal symmetry :
q points will be generated in half the Brillouin zone.<br>
</li>
<li>3=>
Do not take into account any symmetry :
q points will be generated in the full Brillouin zone.<br>
</li>
<li>4=>
Take into account all the symmetries EXCEPT the time-reversal symmetry
to generate the k points in the Irreducible Brillouin Zone.<br>
</li>
</ul>
In the case of a grid of q points, the auxiliary variables
<a href="vargs.html#kptrlen">kptrlen</a>,
<a href="varbas.html#ngkpt">ngkpt</a>  and
<a href="varfil.html#prtkpt">prtkpt</a> might help
you to select the optimal grid, similarly to the case of the K point grid.
"""
},
'qptrlatt': {
'definition': "Q - PoinTs grid : Real space LATTice ",
'section': "vargs",
'category': """<a href="../users/abinit_help.html#not_internal">NOT INTERNAL</a> """,
'vartype': "integer array qptrlatt(3,3)  ",
'default': "No default.",
'text': """This input variable is used only when <a href="vargs.html#qptopt">qptopt</a>
is positive. It partially defines the q point grid.
The other piece of information is contained in
<a href="vargs.html#shiftq">shiftq</a>.
<b>qptrlatt</b> cannot be used together with <a href="vargs.html#ngqpt">ngqpt</a>.

<p>The values <b>qptrlatt</b>(1:3,1), <b>qptrlatt</b>(1:3,2),  <b>qptrlatt</b>(1:3,3)
are the coordinates of three vectors in real space, expressed
in the <a href="varbas.html#rprimd">rprimd</a> coordinate system (reduced coordinates).
They defines a super-lattice in real space.
The k point lattice is the reciprocal of this super-lattice,
possibly shifted (see <a href="vargs.html#shiftq">shiftq</a>).

<p>If neither <a href="vargs.html#ngqpt">ngqpt</a> nor <b>qptrlatt</b>
are defined, ABINIT will automatically generate a set
of k point grids, and select the best combination
of <b>qptrlatt</b> and <a href="vargs.html#shiftq">shiftq</a>
that allows to reach a sufficient value of <a href="vargs.html#kptrlen">kptrlen</a>.
See this latter variable for a complete description of this
procedure.
"""
},
'quadmom': {
'definition': "QUADrupole MOMents ",
'section': "varpaw",
'category': "",
'vartype': """real array quadmom(<a href="varbas.html#ntypat">ntypat</a>) """,
'default': "0",
'text': """Relevant only when
<a href="varint.html#usepaw">usepaw</a>=1, and <a
href="varpaw.html#prtefg">prtefg</a> = 1 or greater.<br>
<ul>
<li> Array of quadrupole moments, in barns, of the nuclei. These values are used in conjunction with the
electric field gradients computed with <a href="varpaw.html#prtefg">prtefg</a> to calculate the quadrupole
couplings in MHz, as well as the asymmetries. Note that the electric field gradient at a nuclear site is independent of
the nuclear quadrupole moment, thus the quadrupole moment of a nucleus can be input as 0, and the
option <a href="varpaw.html#prtefg">prtefg = 2</a> used to determine the electric field gradient at the site.
</li>
</ul>

"""
},
'random_atpos': {
'definition': "RANDOM ATomic POSitions ",
'section': "varrlx",
'category': "",
'vartype': "integer parameter ",
'default': "0.",
'text': """Control the inner coordinates, which can be generated randomly by using 4 different methods depending
ont its value
<br>
(0) if zero, no random generation and xred are taken as they have been introduced by the user
<br>
(1) if one, particles are generated completly random within the unit cell.
<br>
(2) if two, particles are generated randomly but the inner particle distace is always larger than a factor of the
sum of the covalent bonds between the atoms (note : this is incompatible with the definition of alchemical mixing, in which
<a href="varbas.html#ntypat">ntypat</a> differs from <a href="vargs.html#npsp">npsp</a>)
"""
},
'ratsph': {
'definition': "Radii of the ATomic SPHere(s) ",
'section': "vargs",
'category': " ",
'vartype': """real array ratsph(<a href="varbas.html#ntypat">ntypat</a>)  """,
'default': "2.0 Bohr for norm-conserving psps, PAW radius (found in PAW dataset file) for PAW",
'text': """Relevant only when
<a href="varfil.html#prtdos">prtdos</a>=3 or <a href="vargs.html#prtdensph">prtdensph</a>=1.<br>
<br>When <a href="varfil.html#prtdos">prtdos</a>=3:<br>
Provides the radius of the spheres around the <a href="vargs.html#natsph">natsph</a> atoms
of indices <a href="vargs.html#iatsph">iatsph</a>, in which the local
DOS and its angular-momentum projections will be analysed.
The choice of this radius is quite arbitrary. In a plane-wave basis set,
there is no natural definition of an atomic sphere. However, it might be wise
to use the following well-defined and physically motivated procedure
(in version 4.2, this procedure is NOT implemented, unfortunately) :
from the Bader analysis, one can define the radius of the sphere
that contains the same charge as the Bader volume. This
"Equivalent Bader charge atomic radius" might then be used to perform
the present analysis.
See the <A href="../users/aim_help.html"><B>AIM (Bader)</B></a> help file for more explanations.
Another physically motivated choice would be to rely on another
charge partitioning, like the Hirshfeld one (see the cut3d utility).
The advantage of using charge partitioning schemes comes from the fact that the
sum of atomic DOS, for all angular momenta and atoms, integrated on the
energy range of the occupied states,
gives back the total charge.
If this is not an issue, one could rely on the half of the nearest-neighbour distances, or
any scheme that allows to define an atomic radius. Note that the choice of this
radius is however critical for the balance between the s, p and d components. Indeed,
the integrated charge within a given radius, behave as a different power of the
radius, for the different channels s, p, d. At the limit of very small radii, the s component
dominates the charge contained in the sphere ...
<br><br>When <a href="vargs.html#prtdensph">prtdensph</a>=1:<br>
Provides the radius of the spheres around (all) atoms in which the total charge density will be integrated.
<br><br>In case of PAW, <a href="vargs.html#ratsph">ratsph</a> radius has to be greater or equal to
PAW radius of considered atom type (which is read from the PAW dataset file; see rc_sph or r_paw).
"""
},
'rcut': {
'definition': "Radius of the CUT-off for coulomb interaction ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>  """,
'vartype': "real parameter ",
'default': "0.0d0",
'text': """<p>
Truncation of the Coulomb interaction in real space. The meaning of <b>rcut</b> is governed by the cutoff shape option <a href="vargw.html#icutcoul">icutcoul</a>.
<p>
If <b>rcut</b> is negative, the cutoff is automatically calculated so to enclose the same volume inside the cutoff as the volume of the solid.
"""
},
'recefermi': {
'definition': "RECursion - initial guess  of the FERMI Energy",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>  """,
'vartype': "real",
'default': "0",
'text': """Used in Recursion method (<a href="vardev.html#tfkinfunc">tfkinfunc</a>=2).
In the first  SCF calculation it fixes the initial guess for the Fermi energy.
"""
},
'recgratio': {
'definition': "RECursion - Grid Ratio ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>  """,
'vartype': "integer",
'default': "1",
'text': """Used in Recursion method (<a href="vardev.html#tfkinfunc">tfkinfunc</a>=2).
It represents the ratio of the two grid step: <b>recgratio</b>=fine_step/coarse_step and
it is bigger or equal than 1.  It introduces a double-grid system which permits
to compute the electronic density on a coarse grid, using a fine grid
(defined by <a href="vargs.html#ngfft">ngfft</a>) in the
discretisation of the green kernel (see <ahref="vardev.html#recptrott">recptrott</a>).
Successively the density and the recursion coefficients are interpolated on the fine grid by
FFT interpolation.  Note that ngfft/recgratio=number of points of the
coarse grid has to be compatible with the parallelization parameters.
"""
},
'recnpath': {
'definition': "RECursion - Number of point for PATH integral calculations ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>  """,
'vartype': "integer",
'default': "500",
'text': """Used in Recursion method (<a href="vardev.html#tfkinfunc">tfkinfunc</a>=2).
Determine the number of discretisation points to compute some path
integral in the recursion method ; those path integrals are used to
compute the entropy and the eigenvalues energy. during the latest SFC
cycles. """
},
'recnrec': {
'definition': "RECursion - Number of RECursions ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>  """,
'vartype': "integer",
'default': "10",
'text': """Used in Recursion method (<a href="vardev.html#tfkinfunc">tfkinfunc</a>=2).
Determine the maximum order of recursion, that is the dimension of the
krylov space we use to compute density. If the precision setten by
<a href="vardev.html#rectolden">rectolden</a> is reached before that order, the recursion method
automatically stops.  """
},
'recptrott': {
'definition': "RECursion - TROTTer P parameter ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>  """,
'vartype': "integer",
'default': "0",
'text': """Used in Recursion method (<a href="vardev.html#tfkinfunc">tfkinfunc</a>=2).
Determine the trotter parameter used to compute the exponential of the hamiltonian in the
recursion method: exp(-beta*(-Delta + V)) ~ (exp(-beta/(4*recptrott)
V) exp(-beta/(4*recptrott) Delta) exp(-beta/(4*recptrott)
V))^(2*recptrott).
If set to 0, we use recptrott = 1/2 in the above formula.
Increasing <b>recptrott</b> improve the accuracy of the trotter formula, but
increase the dicretisation error: it may be necessary to increase
<a href="vargs.html#ngfft">ngfft</a>. The discretisation error is essentially the discretisation
error of the green kernel exp((recptrott/beta*|r|^2)) on the ngfft
grid.
"""
},
'recrcut': {
'definition': "RECursion - CUTing Radius  ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>  """,
'vartype': "integer",
'default': "0",
'text': """Used in Recursion method (<a href="vardev.html#tfkinfunc">tfkinfunc</a>=2).
Used to improve the
computational time in the case of the recursion method in a large
cell: the density at a point will be computed with taking account only of
a sphere of radius <b>recrcut</b>.
"""
},
'rectesteg': {
'definition': "RECursion - TEST on Electron Gas  ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>  """,
'vartype': "integer",
'default': "0",
'text': """Used in Recursion method (<a href="vardev.html#tfkinfunc">tfkinfunc</a>=2).
It is used to test an electron gas by putting the ion potential
equal to zero.
"""
},
'rectolden': {
'definition': "RECursion - TOLerance on the difference of electronic DENsity  ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>  """,
'vartype': "real",
'default': "0.0E00 (to change)",
'text': """Used in Recursion method (<a href="vardev.html#tfkinfunc">tfkinfunc</a>=2).
Sets a tolerance for differences of electronic density that, reached TWICE
successively, will cause one SCF cycle to stop. That electronic
density difference is computed in the infinity norm (that is, it is
computed point-by-point, and then the maximum difference is computed).
"""
},
'red_dfield': {
'definition': "REDuced Displacement FIELD ",
'section': "varff",
'category': " ",
'vartype': "real array red_dfield(3)  ",
'default': "3*0.0 .",
'text': """In case <a href="varff.html#berryopt">berryopt</a>=16 or 17,
a reduced finite electric displacement field calculation is performed. The value
of this displacement field, and its direction is determined by <b>red_dfield</b>.
It must be given in atomic units.
<p>
<b>red_dfield</b> is defined as Eq.(26) in the Suppl. of Nat. Phys. (M. Stengel, N.A. Spaldin and D. Vanderbilt, Nat. Phys. 5,304 (2009))
"""
},
'red_efield': {
'definition': "REDuced Electric FIELD ",
'section': "varff",
'category': " ",
'vartype': "real array red_efield(3)  ",
'default': "3*0.0 .",
'text': """In case <a href="varff.html#berryopt">berryopt</a>=16 or 17,
a reduced finite electric displacement field calculation is performed. The initial value
of this electric field, and its direction is determined by <b>red_efield</b>.
It must be given in atomic units.
<p>
<b>red_efield</b> is defined as Eq.(25) in the Suppl. of Nat. Phys. (M. Stengel, N.A. Spaldin and D. Vanderbilt, Nat. Phys. 5,304 (2009))
"""
},
'red_efieldbar': {
'definition': "REDuced Electric FIELD BAR ",
'section': "varff",
'category': " ",
'vartype': "real array red_efieldbar(3)  ",
'default': "3*0.0 .",
'text': """In case <a href="varff.html#berryopt">berryopt</a>=14,
a reduced finite electric field calculation is performed. The value
of this electric field, and its direction is determined by <b>red_efieldbar</b>.
It must be given in atomic units.
<p>
<b>red_efieldbar</b> is defined as Eq.(28) in the Suppl. of Nat. Phys. (M. Stengel, N.A. Spaldin and D. Vanderbilt, Nat. Phys. 5,304 (2009))
"""
},
'restartxf': {
'definition': "RESTART from (X,F) history ",
'section': "varrlx",
'category': "",
'vartype': "integer parameter ",
'default': "0.",
'text': """Control the restart of a molecular dynamics or structural
optimization job.
<br><br>
<b>restartxf>0 (Deprecated)</b>:The code reads from the input wf file,
the previous history of atomic coordinates and corresponding forces, in order
to continue the work done by the job that produced this wf file.
If <a href="varrlx.html#optcell">optcell</a>/=0, the history of
<a href="varbas.html#acell">acell</a> and
<a href="varbas.html#rprim">rprim</a> variables is also taken into account.

The code will take into consideration the whole history (if <b>restartxf</b>=1),
or discard the few first (x,f) pairs, and begin only at the
pair whose number corresponds to <b>restartxf</b>.
<br>
Works only for <a href="varrlx.html#ionmov">ionmov</a>=2 (Broyden) and
when an input wavefunction file is specified, thanks to the
appropriate values of <a href="varfil.html#irdwfk">irdwfk</a> or <a
href="varfil.html#getwfk">getwfk</a>.
<br><br>
NOTES :<br>
* The input wf file must have been produced by a run that exited cleanly.
It cannot be one of the temporary wf files that exist when a job crashed.
<br>
* One cannot restart a calculation with a non-zero <a href="varrlx.html#optcell">optcell</a>
value from the (x,f) history of another run with a different non-zero <a
href="varrlx.html#optcell">optcell</a> value. Starting a non-zero <a
href="varrlx.html#optcell">optcell</a> run from a zero
<a href="varrlx.html#optcell">optcell</a> run should work.<br>
* Deprecated, the use of the new options (-1 and -2) is prefered.
<br><br>
<b>restartxf=0 (Default)</b>: No restart procedure is enable
and will start a Molecular dynamics or structural optimization
from scratch.
<br><br>
<b>restartxf=-1 (New)</b>: Use the HIST file to reconstruct the
a partial calculation. It will reconstruct the different configurations
using the forces and stress store in the HIST file, instead of calling
the SCF procedure.
<br>
Enable <b>restartxf=-1</b> from the begining is harmless.
The only condition is to keep the input file the same in such a way
that the same predictor is used and it will predict the same structure
recorded in the HIST file.
<br>
This option will always compute extra <a href="varrlx.html#ntime">ntime</a>
iterations independent of the number of iterations recovered previously.
<br><br>
<b>restartxf=-2 (New)</b>:Read the HIST file and select the atomic
positions and cell parameters with the lowest energy. Forget all the
history and start the calculation using those values. The original
atomic coordinates and cell parameters are irrelevant in that case.
<br><br>
NOTES:
<br>
* You can use <b>restartxf=-1 or -2</b> for all predictiors that
make no use of random numbers.<br>
* You can use <b>restartxf=-1 or -2</b> to restart a calculation
that was not completed. The HIST file is written on each iteration. So
you always have something to recover from.<br>
* You can take advantage of the appropriate values
of <a href="varfil.html#irdwfk">irdwfk</a> or <a
href="varfil.html#getwfk">getwfk</a> to get a good wave function
to continue your job.
"""
},
'rf1atpol': {
'definition': "non-linear Response Function, 1st mixed perturbation : limits of ATomic POLarisations ",
'section': "varrf",
'category': "NON-LINEAR  ",
'vartype': "integer array of 2 elements  ",
'default': "1 1",
'text': """Control the range
of atoms for which displacements will be considered
in phonon calculations (atomic polarisations), or in non-linear
computations, using the 2n+1 theorem.
<br>These values are only relevant to phonon response function
calculations, or non-linear computations.
<br>May take values from 1 to <a href="varbas.html#natom">natom</a>, with <b>rfatpol</b>(1)&lt;=<b>rfatpol</b>(2).
<br>The atoms to be moved will be defined by the
<br>do-loop variable iatpol :
<br>do iatpol=<b>rfatpol</b>(1),<b>rfatpol</b>(2)
<br>For the calculation of a full dynamical matrix, use
<b>rfatpol</b>(1)=1 and <b>rfatpol</b>(2)=<a href="varbas.html#natom">natom</a>, together with
<a href="varrf.html#rfdir">rfdir</a> 1 1 1 . For selected elements of the
dynamical matrix, use different values of <b>rfatpol</b> and/or
<a href="varrf.html#rfdir">rfdir</a>. The name 'iatpol' is used for the part of the
internal variable ipert when it runs from 1 to <a href="varbas.html#natom">natom</a>. The
internal variable ipert can also assume values larger
than <a href="varbas.html#natom">natom</a>, of electric field or stress type (see respfn.help)."""
},
'rf1dir': {
'definition': "non-linear Response Function, 1st mixed perturbation : DIRections ",
'section': "varrf",
'category': "NON-LINEAR  ",
'vartype': "integer array of 3 elements  ",
'default': "0 0 0.",
'text': """Gives the directions
to be considered for response function calculations, or non-linear computations
(also for the Berry phase computation of the polarization, see
the <a href="varff.html#berryopt">berryopt</a> input variable).
<br>The three elements corresponds to the three primitive
vectors, either in real space (phonon calculations),
or in reciprocal space (d/dk, homogeneous electric field, homogeneous magnetic field
calculations). So, they generate a basis
for the generation of the dynamical matrix or
the macroscopic dielectric tensor or magnetic susceptibility and magnetic
shielding, or the effective
charge tensors.
<br>If equal to 1, response functions, as defined
by <a href="varrf.html#rfddk">rfddk</a>,
<a href="varrf.html#rfelfd">rfelfd</a>, <a href="varrf.html#rfphon">rfphon</a>, <b>rfdir</b>
and <a href="varrf.html#rfatpol">rfatpol</a>, are to be computed
for the corresponding direction. If 0, this direction
should not be considered (for non-linear computations, the corresponding input
variables should be used)."""
},
'rf1elfd': {
'definition': "non-linear Response Function, 1st mixed perturbation : ELectric FielD ",
'section': "varrf",
'category': "NON-LINEAR  ",
'vartype': "integer parameter  ",
'default': "0.",
'text': """Turns on electric field response
function calculations (or non-linear computation, including the electric field
perturbation). Actually, such calculations
requires first the non-self-consistent calculation
of derivatives with respect to k, independently of the
electric field perturbation itself.
<br>
Only <a href="varrf.html#rfelfd">rfelfd</a>
is compatible with both norm-conserving pseudopotentials as well as PAW. Higher mixed
perturbations can be used only with norm-conserving pseudopotentials.
<ul>
<li>0=>no electric field perturbation</li>
<li>1=>full calculation, with first the
derivative of ground-state wavefunction with
respect to k (d/dk calculation), by a
non-self-consistent calculation, then the generation of
the first-order response to an homogeneous
electric field</li>
<li>2=>only the derivative of ground-state wavefunctions with
respect to k</li>
<li>3=>only the generation of the first-order response
to the electric field,
assuming that the data on derivative of ground-state
wavefunction with respect to k is available on disk.</li>
</ul>
(Note : because the tolerances to be used for derivatives or
homogeneous electric field are different, one often does the
calculation of derivatives in a separate dataset, followed by
calculation of electric field response as well as phonon.
<br>The options 2 and 3 proves useful in that context ;
also, in case a scissor shift is to be used,
it is usually not applied for the d/dk response). """
},
'rf1phon': {
'definition': "non-linear Response Function, 1st mixed perturbation : PHONons ",
'section': "varrf",
'category': "NON-LINEAR  ",
'vartype': "integer parameter ",
'default': "0.",
'text': """It must be equal to 1
to run phonon response function calculations, or to include some phonon perturbation
in non-linear computations.   """
},
'rf2atpol': {
'definition': "non-linear Response Function, 2nd mixed perturbation : limits of ATomic POLarisations ",
'section': "varrf",
'category': "NON-LINEAR  ",
'vartype': "integer array of 2 elements  ",
'default': "1 1",
'text': """Control the range
of atoms for which displacements will be considered
in phonon calculations (atomic polarisations), or in non-linear
computations, using the 2n+1 theorem.
<br>These values are only relevant to phonon response function
calculations, or non-linear computations.
<br>May take values from 1 to <a href="varbas.html#natom">natom</a>, with <b>rfatpol</b>(1)&lt;=<b>rfatpol</b>(2).
<br>The atoms to be moved will be defined by the
<br>do-loop variable iatpol :
<br>do iatpol=<b>rfatpol</b>(1),<b>rfatpol</b>(2)
<br>For the calculation of a full dynamical matrix, use
<b>rfatpol</b>(1)=1 and <b>rfatpol</b>(2)=<a href="varbas.html#natom">natom</a>, together with
<a href="varrf.html#rfdir">rfdir</a> 1 1 1 . For selected elements of the
dynamical matrix, use different values of <b>rfatpol</b> and/or
<a href="varrf.html#rfdir">rfdir</a>. The name 'iatpol' is used for the part of the
internal variable ipert when it runs from 1 to <a href="varbas.html#natom">natom</a>. The
internal variable ipert can also assume values larger
than <a href="varbas.html#natom">natom</a>, of electric field or stress type (see respfn.help)."""
},
'rf2dir': {
'definition': "non-linear Response Function, 2nd mixed perturbation : DIRections ",
'section': "varrf",
'category': "NON-LINEAR  ",
'vartype': "integer array of 3 elements  ",
'default': "0 0 0.",
'text': """Gives the directions
to be considered for response function calculations, or non-linear computations
(also for the Berry phase computation of the polarization, see
the <a href="varff.html#berryopt">berryopt</a> input variable).
<br>The three elements corresponds to the three primitive
vectors, either in real space (phonon calculations),
or in reciprocal space (d/dk, homogeneous electric field, homogeneous magnetic field
calculations). So, they generate a basis
for the generation of the dynamical matrix or
the macroscopic dielectric tensor or magnetic susceptibility and magnetic
shielding, or the effective
charge tensors.
<br>If equal to 1, response functions, as defined
by <a href="varrf.html#rfddk">rfddk</a>,
<a href="varrf.html#rfelfd">rfelfd</a>, <a href="varrf.html#rfphon">rfphon</a>, <b>rfdir</b>
and <a href="varrf.html#rfatpol">rfatpol</a>, are to be computed
for the corresponding direction. If 0, this direction
should not be considered (for non-linear computations, the corresponding input
variables should be used)."""
},
'rf2elfd': {
'definition': "non-linear Response Function, 2nd mixed perturbation : ELectric FielD ",
'section': "varrf",
'category': "NON-LINEAR  ",
'vartype': "integer parameter  ",
'default': "0.",
'text': """Turns on electric field response
function calculations (or non-linear computation, including the electric field
perturbation). Actually, such calculations
requires first the non-self-consistent calculation
of derivatives with respect to k, independently of the
electric field perturbation itself.
<br>
Only <a href="varrf.html#rfelfd">rfelfd</a>
is compatible with both norm-conserving pseudopotentials as well as PAW. Higher mixed
perturbations can be used only with norm-conserving pseudopotentials.
<ul>
<li>0=>no electric field perturbation</li>
<li>1=>full calculation, with first the
derivative of ground-state wavefunction with
respect to k (d/dk calculation), by a
non-self-consistent calculation, then the generation of
the first-order response to an homogeneous
electric field</li>
<li>2=>only the derivative of ground-state wavefunctions with
respect to k</li>
<li>3=>only the generation of the first-order response
to the electric field,
assuming that the data on derivative of ground-state
wavefunction with respect to k is available on disk.</li>
</ul>
(Note : because the tolerances to be used for derivatives or
homogeneous electric field are different, one often does the
calculation of derivatives in a separate dataset, followed by
calculation of electric field response as well as phonon.
<br>The options 2 and 3 proves useful in that context ;
also, in case a scissor shift is to be used,
it is usually not applied for the d/dk response). """
},
'rf2phon': {
'definition': "non-linear Response Function, 2nd mixed perturbation : PHONons ",
'section': "varrf",
'category': "NON-LINEAR  ",
'vartype': "integer parameter ",
'default': "0.",
'text': """It must be equal to 1
to run phonon response function calculations, or to include some phonon perturbation
in non-linear computations.   """
},
'rf3atpol': {
'definition': "non-linear Response Function, 3rd mixed perturbation : limits of ATomic POLarisations ",
'section': "varrf",
'category': "NON-LINEAR  ",
'vartype': "integer array of 2 elements  ",
'default': "1 1",
'text': """Control the range
of atoms for which displacements will be considered
in phonon calculations (atomic polarisations), or in non-linear
computations, using the 2n+1 theorem.
<br>These values are only relevant to phonon response function
calculations, or non-linear computations.
<br>May take values from 1 to <a href="varbas.html#natom">natom</a>, with <b>rfatpol</b>(1)&lt;=<b>rfatpol</b>(2).
<br>The atoms to be moved will be defined by the
<br>do-loop variable iatpol :
<br>do iatpol=<b>rfatpol</b>(1),<b>rfatpol</b>(2)
<br>For the calculation of a full dynamical matrix, use
<b>rfatpol</b>(1)=1 and <b>rfatpol</b>(2)=<a href="varbas.html#natom">natom</a>, together with
<a href="varrf.html#rfdir">rfdir</a> 1 1 1 . For selected elements of the
dynamical matrix, use different values of <b>rfatpol</b> and/or
<a href="varrf.html#rfdir">rfdir</a>. The name 'iatpol' is used for the part of the
internal variable ipert when it runs from 1 to <a href="varbas.html#natom">natom</a>. The
internal variable ipert can also assume values larger
than <a href="varbas.html#natom">natom</a>, of electric field or stress type (see respfn.help)."""
},
'rf3dir': {
'definition': "non-linear Response Function, 3rd mixed perturbation : DIRections ",
'section': "varrf",
'category': "NON-LINEAR  ",
'vartype': "integer array of 3 elements  ",
'default': "0 0 0.",
'text': """Gives the directions
to be considered for response function calculations, or non-linear computations
(also for the Berry phase computation of the polarization, see
the <a href="varff.html#berryopt">berryopt</a> input variable).
<br>The three elements corresponds to the three primitive
vectors, either in real space (phonon calculations),
or in reciprocal space (d/dk, homogeneous electric field, homogeneous magnetic field
calculations). So, they generate a basis
for the generation of the dynamical matrix or
the macroscopic dielectric tensor or magnetic susceptibility and magnetic
shielding, or the effective
charge tensors.
<br>If equal to 1, response functions, as defined
by <a href="varrf.html#rfddk">rfddk</a>,
<a href="varrf.html#rfelfd">rfelfd</a>, <a href="varrf.html#rfphon">rfphon</a>, <b>rfdir</b>
and <a href="varrf.html#rfatpol">rfatpol</a>, are to be computed
for the corresponding direction. If 0, this direction
should not be considered (for non-linear computations, the corresponding input
variables should be used)."""
},
'rf3elfd': {
'definition': "non-linear Response Function, 3rd mixed perturbation : ELectric FielD ",
'section': "varrf",
'category': "NON-LINEAR  ",
'vartype': "integer parameter  ",
'default': "0.",
'text': """Turns on electric field response
function calculations (or non-linear computation, including the electric field
perturbation). Actually, such calculations
requires first the non-self-consistent calculation
of derivatives with respect to k, independently of the
electric field perturbation itself.
<br>
Only <a href="varrf.html#rfelfd">rfelfd</a>
is compatible with both norm-conserving pseudopotentials as well as PAW. Higher mixed
perturbations can be used only with norm-conserving pseudopotentials.
<ul>
<li>0=>no electric field perturbation</li>
<li>1=>full calculation, with first the
derivative of ground-state wavefunction with
respect to k (d/dk calculation), by a
non-self-consistent calculation, then the generation of
the first-order response to an homogeneous
electric field</li>
<li>2=>only the derivative of ground-state wavefunctions with
respect to k</li>
<li>3=>only the generation of the first-order response
to the electric field,
assuming that the data on derivative of ground-state
wavefunction with respect to k is available on disk.</li>
</ul>
(Note : because the tolerances to be used for derivatives or
homogeneous electric field are different, one often does the
calculation of derivatives in a separate dataset, followed by
calculation of electric field response as well as phonon.
<br>The options 2 and 3 proves useful in that context ;
also, in case a scissor shift is to be used,
it is usually not applied for the d/dk response). """
},
'rf3phon': {
'definition': "non-linear Response Function, 3rd mixed perturbation : PHONons ",
'section': "varrf",
'category': "NON-LINEAR  ",
'vartype': "integer parameter ",
'default': "0.",
'text': """It must be equal to 1
to run phonon response function calculations, or to include some phonon perturbation
in non-linear computations.   """
},
'rfasr': {
'definition': "Response Function : Acoustic Sum Rule ",
'section': "varrf",
'category': """<a href="../users/abinit_help.html#respfn">RESPFN</a>  """,
'vartype': "integer parameter  ",
'default': "0.",
'text': """Control the evaluation of the
acoustic sum rule in effective charges and dynamical matrix at Gamma
within a response function calculation (not active at the level of producing the DDB, but
at the level of the phonon eigenfrequencies output).
<ul>
<li>0 => no acoustic sum rule imposed</li>
<li>1 => acoustic sum rule imposed for dynamical matrix at Gamma, and charge neutrality imposed with
extra charge evenly distributed among atoms</li>
<li>2 => acoustic sum rule imposed for dynamical matrix at Gamma, and charge neutrality imposed with
extra charge given proportionally to those atoms with
the largest effective charge.</li>
</ul>
The treatment of the acoustic sum rule and charge neutrality sum rule is finer at the level of the ANADDB utility,
with the two independent input variables <a href="../users/anaddb_help.html#asr">asr</a>
and <a href="../users/anaddb_help.html#chneut">chneut</a>.
"""
},
'rfatpol': {
'definition': "Response Function : limits of ATomic POLarisations ",
'section': "varrf",
'category': "NON-LINEAR  ",
'vartype': "integer array of 2 elements  ",
'default': "1 1",
'text': """Control the range
of atoms for which displacements will be considered
in phonon calculations (atomic polarisations), or in non-linear
computations, using the 2n+1 theorem.
<br>These values are only relevant to phonon response function
calculations, or non-linear computations.
<br>May take values from 1 to <a href="varbas.html#natom">natom</a>, with <b>rfatpol</b>(1)&lt;=<b>rfatpol</b>(2).
<br>The atoms to be moved will be defined by the
<br>do-loop variable iatpol :
<br>do iatpol=<b>rfatpol</b>(1),<b>rfatpol</b>(2)
<br>For the calculation of a full dynamical matrix, use
<b>rfatpol</b>(1)=1 and <b>rfatpol</b>(2)=<a href="varbas.html#natom">natom</a>, together with
<a href="varrf.html#rfdir">rfdir</a> 1 1 1 . For selected elements of the
dynamical matrix, use different values of <b>rfatpol</b> and/or
<a href="varrf.html#rfdir">rfdir</a>. The name 'iatpol' is used for the part of the
internal variable ipert when it runs from 1 to <a href="varbas.html#natom">natom</a>. The
internal variable ipert can also assume values larger
than <a href="varbas.html#natom">natom</a>, of electric field or stress type (see respfn.help)."""
},
'rfddk': {
'definition': "Response Function with respect to Derivative with respect to K ",
'section': "varrf",
'category': """<a href="../users/abinit_help.html#respfn">RESPFN</a>  """,
'vartype': "integer parameter  ",
'default': "0.",
'text': """Activates computation of derivatives of ground state
wavefunctions with respect to wavevectors. This is not strictly a response
function but is a needed auxiliary quantity for example in the electric field and
magnetic field
calculations (see <a href="varrf.html#rfelfd">rfelfd</a>) The directions for the
derivatives are determined by <a href="varrf.html#rfdir">rfdir</a>.
<ul>
<li>0=>no derivative calculation</li>
<li>1=>calculation of first derivatives of wavefunctions with respect to k points
(d/dk calculation). The exact same functionality is provided by
<a href="varrf.html#rfelfd">rfelfd</a> = 2.
</li>
<li>2=>will activate second derivative calculation but this is NOT YET implemented.</li>
</ul>
"""
},
'rfdir': {
'definition': "Response Function : DIRections ",
'section': "varrf",
'category': "NON-LINEAR  ",
'vartype': "integer array of 3 elements  ",
'default': "0 0 0.",
'text': """Gives the directions
to be considered for response function calculations, or non-linear computations
(also for the Berry phase computation of the polarization, see
the <a href="varff.html#berryopt">berryopt</a> input variable).
<br>The three elements corresponds to the three primitive
vectors, either in real space (phonon calculations),
or in reciprocal space (d/dk, homogeneous electric field, homogeneous magnetic field
calculations). So, they generate a basis
for the generation of the dynamical matrix or
the macroscopic dielectric tensor or magnetic susceptibility and magnetic
shielding, or the effective
charge tensors.
<br>If equal to 1, response functions, as defined
by <a href="varrf.html#rfddk">rfddk</a>,
<a href="varrf.html#rfelfd">rfelfd</a>, <a href="varrf.html#rfphon">rfphon</a>, <b>rfdir</b>
and <a href="varrf.html#rfatpol">rfatpol</a>, are to be computed
for the corresponding direction. If 0, this direction
should not be considered (for non-linear computations, the corresponding input
variables should be used)."""
},
'rfelfd': {
'definition': "Response Function with respect to the ELectric FielD ",
'section': "varrf",
'category': "NON-LINEAR  ",
'vartype': "integer parameter  ",
'default': "0.",
'text': """Turns on electric field response
function calculations (or non-linear computation, including the electric field
perturbation). Actually, such calculations
requires first the non-self-consistent calculation
of derivatives with respect to k, independently of the
electric field perturbation itself.
<br>
Only <a href="varrf.html#rfelfd">rfelfd</a>
is compatible with both norm-conserving pseudopotentials as well as PAW. Higher mixed
perturbations can be used only with norm-conserving pseudopotentials.
<ul>
<li>0=>no electric field perturbation</li>
<li>1=>full calculation, with first the
derivative of ground-state wavefunction with
respect to k (d/dk calculation), by a
non-self-consistent calculation, then the generation of
the first-order response to an homogeneous
electric field</li>
<li>2=>only the derivative of ground-state wavefunctions with
respect to k</li>
<li>3=>only the generation of the first-order response
to the electric field,
assuming that the data on derivative of ground-state
wavefunction with respect to k is available on disk.</li>
</ul>
(Note : because the tolerances to be used for derivatives or
homogeneous electric field are different, one often does the
calculation of derivatives in a separate dataset, followed by
calculation of electric field response as well as phonon.
<br>The options 2 and 3 proves useful in that context ;
also, in case a scissor shift is to be used,
it is usually not applied for the d/dk response). """
},
'rfmeth': {
'definition': "Response Function METHod ",
'section': "varrf",
'category': """<a href="../users/abinit_help.html#respfn">RESPFN</a>  """,
'vartype': "integer parameter  ",
'default': "1.",
'text': """Selects method used in
response function calculations. Presently, only 1
is allowed. """
},
'rfphon': {
'definition': "Response Function with respect to PHONons ",
'section': "varrf",
'category': "NON-LINEAR  ",
'vartype': "integer parameter ",
'default': "0.",
'text': """It must be equal to 1
to run phonon response function calculations, or to include some phonon perturbation
in non-linear computations.   """
},
'rfstrs': {
'definition': "Response Function with respect to STRainS ",
'section': "varrf",
'category': """<a href="../users/abinit_help.html#respfn">RESPFN</a>  """,
'vartype': "integer parameter  ",
'default': "0.",
'text': """Used to run strain response-function
calculations (e.g. needed to get elastic constants). Define, with
<a href="varrf.html#rfdir">rfdir</a>, the set of perturbations.
<ul>
<li>0=>no strain perturbation</li>
<li>1=>only uniaxial strain(s) (ipert=natom+3 is activated)</li>
<li>2=>only shear strain(s) (ipert=natom+4 is activated)</li>
<li>3=>both uniaxial and shear strain(s) (both ipert=natom+3 and ipert=natom+4 are activated)</li>
</ul>
See the possible restrictions on the use of strain perturbations, in the
<a href="../users/respfn_help.html#1">respfn_help file</a>.
"""
},
'rfuser': {
'definition': "Response Function, USER-defined  ",
'section': "varrf",
'category': """<a href="../users/abinit_help.html#respfn">RESPFN</a>  """,
'vartype': "integer parameter  ",
'default': "0. ",
'text': """Available to the developpers, to activate
the use of ipert=natom+5 and ipert=natom+6, two sets of perturbations
that the developpers can define.
<br>
<ul>
<li>0=>no computations for ipert=natom+5 or ipert=natom+6</li>
<li>1=>response with respect to perturbation natom+5 will be computed </li>
<li>2=>response with respect to perturbation natom+6 will be computed </li>
<li>3=>responses with respect to perturbations natom+5 and natom+6 will be computed </li>
</ul>
<p>In order to define and use correctly the new perturbations,
the developper might have to include code lines or additional routines
at the level of the following routines :
cgwf3.f, chkph3.f, dyout3.f, d2sym3.f, eneou3.f, eneres3.f, gath3.f, insy3.f,
loper3.f, mkcor3.f, nstdy3.f, nstwf3.f, respfn.f,
scfcv3.f, syper3.f, vloca3.f, vtorho3.f, vtowfk3.f, wings3.f, .
In these routines, the developper should pay a particular
attention to the rfpert array, defined in the routine respfn.f ,
as well as to the ipert local variable.
"""
},
'rhoqpmix': {
'definition': "RHO QuasiParticle MIXing",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a> """,
'vartype': "real parameter  ",
'default': "1.0",
'text': """<p>
For self-consistent GW runs, <b>rhoqpmix</b> sets the mixing coefficient between the new and the previous
electronic densities.
This mixing damps the spurious oscillations in the Hartree potential when achieving self-consistency.
<b>rhoqpmix</b> is meaningful only when doing self-consistency on the wavefunctions with
<a href="vargw.html#gwcalctyp">gwcalctyp</a> &#62&#61 20.
"""
},
'rprim': {
'definition': "Real space PRIMitive translations ",
'section': "varbas",
'category': """<a href="../users/abinit_help.html#evolving">EVOLVING</a> (if <a href="varrlx.html#ionmov">ionmov</a>==2 and <a href="varrlx.html#optcell">optcell</a>/=0) """,
'vartype': """real array rprim(3,3), represented internally as rprim(3,3,<a href="varrlx.html#nimage">nimage</a>)""",
'default': "3x3 unity matrix.",
'text': """Give, in columnwise entry,
the three dimensionless primitive translations in real space, to be rescaled by
<a href="varbas.html#acell">acell</a> and <a href="varbas.html#scalecart">scalecart</a>.
<br>If the Default is used, that is, <b>rprim</b> is the unity matrix,
the three dimensionless primitive vectors are three
unit vectors in cartesian coordinates. Each will be (possibly) multiplied
by the corresponding <a href="varbas.html#acell">acell</a> value, then (possibly)
stretched along the cartesian coordinates by the corresponding <a href="varbas.html#scalecart">scalecart</a> value,
to give the dimensional primitive vectors, called <a href="varbas.html#rprimd">rprimd</a>.
<br>In the general case, the dimensional cartesian
coordinates of the crystal primitive translations R1p, R2p and R3p, see
<a href="varbas.html#rprimd">rprimd</a>, are
<ul>
<li>  R1p(i)=<a href="varbas.html#scalecart">scalecart</a>(i)<b>rprim</b>(i,1)*<a href="varbas.html#acell">acell</a>(1) </li>
<li>  R2p(i)=<a href="varbas.html#scalecart">scalecart</a>(i)<b>rprim</b>(i,2)*<a href="varbas.html#acell">acell</a>(2) </li>
<li>  R3p(i)=<a href="varbas.html#scalecart">scalecart</a>(i)<b>rprim</b>(i,3)*<a href="varbas.html#acell">acell</a>(3) </li>
</ul>
where i=1,2,3 is the component of the primitive translation (i.e. x, y, and z).
<br>
<br>
The <b>rprim</b> variable, scaled by <a href="varbas.html#scalecart">scalecart</a>, is thus used to define directions
of the primitive vectors, that will be multiplied (so keeping the direction unchanged) by
the appropriate length scale <a href="varbas.html#acell">acell</a>(1), <a href="varbas.html#acell">acell</a>(2),
or <a href="varbas.html#acell">acell</a>(3),
respectively to give the dimensional primitive translations
in real space in cartesian coordinates.
<br>
Presently, it is requested that the mixed product
(R1xR2).R3 is positive. If this is not the case,
simply exchange a pair of vectors.
<br>
To be more specific, keeping the default value of <a href="varbas.html#scalecart">scalecart</a>=1 to simplify the matter,
<b>rprim</b> 1 2 3 4 5 6 7 8 9 corresponds to input of the
three primitive translations R1=(1,2,3) (to be multiplied by <a
href="varbas.html#acell">acell</a>(1)), R2=(4,5,6) (to be multiplied by <a
href="varbas.html#acell">acell</a>(2)), and R3=(7,8,9) (to be multiplied by <a
href="varbas.html#acell">acell</a>(3)).
<br>
Note carefully that the first
three numbers input are the first column of <b>rprim</b>, the next
three are the second, and the final three are the third.
This corresponds with the usual Fortran order for arrays.
The matrix whose columns are the reciprocal space primitive
translations is the inverse transpose of the matrix whose
columns are the direct space primitive translations.
<p>Alternatively to <b>rprim</b>, directions of dimensionless primitive
vectors can be specified by using the input variable <a href="varbas.html#angdeg">angdeg</a>.
This is especially useful for hexagonal lattices (with 120 or 60 degrees angles).
Indeed, in order for symmetries to be recognized, rprim must be symmetric up to
<a href="vargeo.html#tolsym">tolsym</a> (10 digits by default),
inducing a specification such as
<pre>
rprim  0.86602540378  0.5  0.0
-0.86602540378  0.5  0.0
0.0            0.0  1.0
</pre>
that can be avoided thanks to <a href="varbas.html#angdeg">angdeg</a>:
<pre>
angdeg 90 90 120
</pre>
<br>Although the use of <a href="varbas.html#scalecart">scalecart</a> or <a href="varbas.html#acell">acell</a> is
rather equivalent when the primitive vectors are aligned with the cartesian directions, it is not the case for
non-orthogonal primitive vectors. In particular, beginners often make the error of trying to use <a href="varbas.html#acell">acell</a>
to define primitive vectors in face-centered tetragonal lattice, or body-centered tetragonal lattice, or similarly
in face or body-centered orthorhombic lattices. Let's take the example of a body-centered tetragonal lattice, that
might be defined using the following ("a" and "c" have to be replaced by the appropriate conventional cell vector length):
<pre>
rprim  "a"      0      0
0      "a"     0
"a/2"   "a/2"  "c/2"
acell 3*1     scalecart 3*1    !  ( These are default values)
</pre>
The following is a valid, alternative way to define  the same primitive vectors :
<pre>
rprim   1       0      0
0       1      0
1/2     1/2    1/2
scalecart  "a"  "a"  "c"
acell 3*1    !  ( These are default values)
</pre>
Indeed, the cell has been stretched along the cartesian coordinates, by "a", "a" and "c" factors.
<p>
At variance, the following is WRONG :
<pre>
rprim   1       0      0
0       1      0
1/2     1/2    1/2
acell  "a"  "a"  "c"    !   THIS IS WRONG
scalecart 3*1    !  ( These are default values)
</pre>
Indeed, the latter would correspond to :
<pre>
rprim  "a"      0      0
0      "a"     0
"c/2"   "c/2"  "c/2"
acell 3*1     scalecart 3*1    !  ( These are default values)
</pre>
Namely, the third vector has been rescaled by "c". It is not at all in the center of the tetragonal cell whose basis vectors
are defined by the scaling factor "a".
<br> As another difference between <a href="varbas.html#scalecart">scalecart</a> or <a href="varbas.html#acell">acell</a>,
note that <a href="varbas.html#scalecart">scalecart</a> is <a href="../users/abinit_help.html#not_internal">NOT INTERNAL</a> :
its content will be immediately applied to rprim, at parsing time,
and then scalecart will be set to the default values (3*1). So, in case <a href="varbas.html#scalecart">scalecart</a> is used,
the echo of <b>rprim</b> in the output file is not the value contained in the input file,
but the value rescaled by <a href="varbas.html#scalecart">scalecart</a>.
"""
},
'rprimd': {
'definition': "Real space PRIMitive translations, Dimensional ",
'section': "varbas",
'category': """<a href="../users/abinit_help.html#internal">INTERNAL</a>, <a href="../users/abinit_help.html#evolving">EVOVING</a>(if <a href="varrlx.html#ionmov">ionmov</a>==2 and <a href="varrlx.html#optcell">optcell</a>/=0) """,
'vartype': """real array rprimd(3,3), represented internally as rprimd(3,3,<a href="varrlx.html#nimage">nimage</a>)""",
'default': "",
'text': """This internal variable gives the dimensional real space primitive
vectors, computed from <a href="varbas.html#acell">acell</a>,
<a href="varbas.html#scalecart">scalecart</a>,
and <a href="varbas.html#rprim">rprim</a>.
<ul>
<li>  R1p(i)=<b>rprimd</b>(i,1)=<a href="varbas.html#scalecart">scalecart</a>(i)*<a href="varbas.html#rprim">rprim</a>(i,1)*<a href="varbas.html#acell">acell</a>(1) for i=1,2,3 (x,y,and z)</li>
<li>  R2p(i)=<b>rprimd</b>(i,2)=<a href="varbas.html#scalecart">scalecart</a>(i)*<a href="varbas.html#rprim">rprim</a>(i,2)*<a href="varbas.html#acell">acell</a>(2) for i=1,2,3</li>
<li>  R3p(i)=<b>rprimd</b>(i,3)=<a href="varbas.html#scalecart">scalecart</a>(i)*<a href="varbas.html#rprim">rprim</a>(i,3)*<a href="varbas.html#acell">acell</a>(3) for i=1,2,3</li>
</ul>
"""
},
'scalecart': {
'definition': "SCALE CARTesian coordinates",
'section': "varbas",
'category': """<a href="../users/abinit_help.html#not_internal">NOT INTERNAL</a>  """,
'vartype': "real array scalecart(3)   ",
'default': "3*1  ",
'text': """Gives the scaling factors of cartesian coordinates by which
dimensionless primitive translations (in "<a href="varbas.html#rprim">rprim</a>") are
to be multiplied.
<a href="varbas.html#rprim">rprim</a> input variable,
the <a href="varbas.html#acell">acell</a> input variable,
and the associated internal <a href="varbas.html#rprimd">rprimd</a> internal variable.
<br> Especially useful for body-centered and face-centered tetragonal lattices, as well as
body-centered and face-centered orthorhombic lattices, see <a href="varbas.html#rprimd">rprimd</a>.
<br> Note that this input variable is NOT INTERNAL : its content will be immediately applied to rprim, at parsing time,
and then scalecart will be set to the default values. So, it will not be echoed.
"""
},
'sciss': {
'definition': "SCISSor operator ",
'section': "varrf",
'category': """<a href="../users/abinit_help.html#respfn">RESPFN</a>, <a href="../users/abinit_help.html#energy">ENERGY</a> """,
'vartype': "real parameter ",
'default': "0.",
'text': """It is the value of the "scissors operator", the
shift of conduction band eigenvalues,
used in response function calculations.
<br>Can be specified in Ha (the default), Ry, eV or Kelvin, since
<b>ecut</b> has the
'<a href="../users/abinit_help.html#dimensions">ENERGY</a>' characteristics.
(1 Ha=27.2113845 eV)
<br>Typical use is for response to electric field (<a href="varrf.html#rfelfd">rfelfd</a>=3),
but NOT for d/dk (<a href="varrf.html#rfelfd">rfelfd</a>=2) and phonon responses."""
},
'scphon_supercell': {
'definition': "Self Consistent PHONon SUPERCELL",
'section': "vargs",
'category': "",
'vartype': "integer array of 3 elements",
'default': "1 1 1 ",
'text': """Give extent, in number of primitive unit cells, of the supercell being used for
a self-consistent phonon calculation. Presumes the phonon frequencies and eigenvectors
have been calculated in the original primitive unit cell, on a grid of q-points which
corresponds to the supercell in the present calculation.

TO BE IMPROVED : should contain a tutorial on how to do self-consistent phonon calculations, David Waroquiers 090831
"""
},
'scphon_temp': {
'definition': "Self Consistent PHONon TEMPerature",
'section': "vargs",
'category': """<a href="../users/abinit_help.html#energy">ENERGY</a> """,
'vartype': "real parameter",
'default': "0.0d0 ",
'text': """Temperature which is imposed on phonon distribution, in the self-consistent scheme of
Souvatzis et al. PRL <b>100</b>, 095901. Determines the extent of the finite displacements
used, and consequent anharmonic effects. Experimental.
"""
},
'shiftk': {
'definition': "SHIFT for K points ",
'section': "varbas",
'category': " ",
'vartype': """real array shift(3,<a href="varbas.html#nshiftk">nshiftk</a>)  """,
'default': "0.5 0.5 0.5 ... 0.5  ",
'text': """It is used only when <a href="varbas.html#kptopt">kptopt</a>>=0,
and must be defined if <a href="varbas.html#nshiftk">nshiftk</a> is larger than 1.
<br><b>shiftk</b>(1:3,1:<a href="varbas.html#nshiftk">nshiftk</a>) defines
<a href="varbas.html#nshiftk">nshiftk</a> shifts
of the homogeneous grid of k points
based on <a href="varbas.html#ngkpt">ngkpt</a> or
<a href="vargs.html#kptrlatt">kptrlatt</a>.
<br>The shifts induced by <b>shiftk</b> corresponds
to the reduced coordinates in the coordinate system
defining the k-point lattice. For example,
if the k point lattice is defined using <a href="varbas.html#ngkpt">ngkpt</a>,
the point whose reciprocal space reduced coordinates are
( <b>shiftk</b>(1,ii)/<a href="varbas.html#ngkpt">ngkpt</a>(1)
<b>shiftk</b>(2,ii)/<a href="varbas.html#ngkpt">ngkpt</a>(2)
<b>shiftk</b>(3,ii)/<a href="varbas.html#ngkpt">ngkpt</a>(3) )
belongs to the shifted grid number ii.

<p>The user might rely on ABINIT to suggest suitable and
efficient combinations of <a href="vargs.html#kptrlatt">kptrlatt</a>
and <a href="varbas.html#shiftk">shiftk</a>.
The procedure to be followed is described with the
input variables <a href="vargs.html#kptrlen">kptrlen</a>.
In what follows, we suggest some interesting values of the shifts,
to be used with even values of <a href="varbas.html#ngkpt">ngkpt</a>.
This list is much less exhaustive than the above-mentioned automatic
procedure.
<p>
1) When the primitive vectors of the lattice do NOT form
a FCC or a BCC lattice, the usual (shifted) Monkhorst-Pack
grids are formed by using
<a href="varbas.html#nshiftk">nshiftk</a>=1 and <b>shiftk</b> 0.5 0.5 0.5 .
This is often the preferred k point sampling.
For a non-shifted Monkhorst-Pack grid, use
<a href="varbas.html#nshiftk">nshiftk</a>=1 and <b>shiftk</b> 0.0 0.0 0.0 ,
but there is little reason to do that.
<p>
2) When the primitive vectors of the lattice form a FCC lattice,
with <a href="varbas.html#rprim">rprim</a>
<pre>
0.0 0.5 0.5
0.5 0.0 0.5
0.5 0.5 0.0
</pre>
the (very efficient) usual Monkhorst-Pack sampling will be generated by using
<a href="varbas.html#nshiftk">nshiftk</a>= 4  and <b>shiftk</b>
<pre>
0.5 0.5 0.5
0.5 0.0 0.0
0.0 0.5 0.0
0.0 0.0 0.5
</pre>
<p>
3) When the primitive vectors of the lattice form a BCC lattice,
with <a href="varbas.html#rprim">rprim</a>
<pre>
-0.5  0.5  0.5
0.5 -0.5  0.5
0.5  0.5 -0.5
</pre>
the usual Monkhorst-Pack sampling will be generated by using
<a href="varbas.html#nshiftk">nshiftk</a>= 2  and <b>shiftk</b>
<pre>
0.25  0.25  0.25
-0.25 -0.25 -0.25
</pre>
However, the simple sampling
<a href="varbas.html#nshiftk">nshiftk</a>=1 and <b>shiftk</b> 0.5 0.5 0.5
is excellent.
<p>
4) For hexagonal lattices with hexagonal axes, e.g. <a href="varbas.html#rprim">rprim</a>
<pre>
1.0  0.0       0.0
-0.5  sqrt(3)/2 0.0
0.0  0.0       1.0
</pre>
one can use
<a href="varbas.html#nshiftk">nshiftk</a>= 1  and <b>shiftk</b>  0.0 0.0 0.5
<p>
In rhombohedral axes, e.g. using <a href="varbas.html#angdeg">angdeg</a> 3*60.,
this corresponds to <b>shiftk</b>  0.5 0.5 0.5, to keep the shift along the
symmetry axis.
"""
},
'shiftq': {
'definition': "SHIFT for Q points ",
'section': "vargs",
'category': """<a href="../users/abinit_help.html#not_internal">NOT INTERNAL</a> """,
'vartype': """real array shift(3,<a href="vargs.html#nshiftq">nshiftq</a>)  """,
'default': "0.5 0.5 0.5 ... 0.5  ",
'text': """It is used only when <a href="vargs.html#qptopt">qptopt</a>>=0,
and must be defined if <a href="vargs.html#nshiftq">nshiftq</a> is larger than 1.
<br><b>shiftq</b>(1:3,1:<a href="vargs.html#nshiftq">nshiftq</a>) defines
<a href="vargs.html#nshiftq">nshiftq</a> shifts
of the homogeneous grid of q points
based on <a href="vargs.html#ngqpt">ngqpt</a> or
<a href="vargs.html#qptrlatt">qptrlatt</a>.
<p>See <a href="varbas.html#shiftk">shiftk</a> for more information on the definition,
use, and suitable values for these shifts.
"""
},
'signperm': {
'definition': "SIGN of PERMutation potential ",
'section': "varrlx",
'category': "",
'vartype': "integer ",
'default': "1.",
'text': """-1 favors segregation
"""
},
'slabwsrad': {
'definition': "jellium SLAB Wigner-Seitz RADius",
'section': "vargs",
'category': """<a href="../users/abinit_help.html#length">LENGTH</a>""",
'vartype': "real parameter ",
'default': "0.0d0. ",
'text': """Fix the bulk-mean positive charge density n<sub>bulk</sub>
of a jellium slab (if the latter is employed, e.g. <a href="vargs.html#jellslab">jellslab</a> &ne 0).
Often called "r<sub>s</sub>" [see for example N. D. Lang and W. Kohn PRB 1, 4555 (1970)],
<b>slabwsrad</b> is the radius of a sphere which has the same volume
as the average volume per particle in a homogeneous electron gas with density n<sub>bulk</sub>, so:
<pre>
1/n<sub>bulk</sub> = 4/3 Pi * <a href="vargs.html#slabwsrad">slabwsrad</a><sup>3</sup>
</pre>
For example, the bulk aluminum fcc lattice constant
is a=4.0495 Angstroms (webelements.com), each cubic centered cell
includes 4 Al atoms and each atom has 3 valence electrons,
so the average volume per electron is a<sup>3</sup>/12=37.34 Bohr<sup>3</sup>
which has to be equal to 4/3 Pi*r<sub>s</sub><sup>3</sup>.
Consequently Al has approximately r<sub>s</sub>=2.07 Bohr, while for example magnesium has r<sub>s</sub>=2.65 Bohr, sodium 3.99 Bohr.
<br>
By default, given in Bohr atomic units (1 Bohr=0.5291772108 Angstroms)."""
},
'slabzbeg': {
'definition': "jellium SLAB BEGinning edge along the Z direction",
'section': "vargs",
'category': "",
'vartype': "real parameter ",
'default': "0.0d0, 0.0d0.  ",
'text': """Define the edges of the jellium slab (if used, so if <a href="vargs.html#jellslab">jellslab</a> &ne 0)
along z, namely the slab starts at a point along z which is expressed in Bohr by <b>slabzbeg</b>
and it ends at a point expressed in Bohr by <b>slabzend</b>.
The z direction is parallel to the third crystal primitive lattice vector
which has to be orthogonal to the other ones,
so the length of the cell along z is <a href="varbas.html#rprimd">rprimd</a>(3,3).
In addition <b>slabzbeg</b> and <b>slabzend</b> have to be such that:
<pre>
0 &le <b>slabzbeg</b> &lt <b>slabzend</b> &le <a href="varbas.html#rprimd">rprimd</a>(3,3)
</pre>
Together with <a href="vargs.html#slabwsrad">slabwsrad</a>
they define the jellium positive charge density distribution n<sub>+</sub>(x,y,z) in this way:
<pre>
n<sub>+</sub>(x,y,z) = n<sub>bulk</sub>   if <b>slabzbeg</b> &le z &le <b>slabzend</b>
= 0      otherwise,
</pre>
so the positive charge density is invariant along the xy plane as well as the electrostatic potential generated by it."""
},
'slabzend': {
'definition': "jellium SLAB ENDing edge along the Z direction",
'section': "vargs",
'category': "",
'vartype': "real parameter ",
'default': "0.0d0, 0.0d0.  ",
'text': """Define the edges of the jellium slab (if used, so if <a href="vargs.html#jellslab">jellslab</a> &ne 0)
along z, namely the slab starts at a point along z which is expressed in Bohr by <b>slabzbeg</b>
and it ends at a point expressed in Bohr by <b>slabzend</b>.
The z direction is parallel to the third crystal primitive lattice vector
which has to be orthogonal to the other ones,
so the length of the cell along z is <a href="varbas.html#rprimd">rprimd</a>(3,3).
In addition <b>slabzbeg</b> and <b>slabzend</b> have to be such that:
<pre>
0 &le <b>slabzbeg</b> &lt <b>slabzend</b> &le <a href="varbas.html#rprimd">rprimd</a>(3,3)
</pre>
Together with <a href="vargs.html#slabwsrad">slabwsrad</a>
they define the jellium positive charge density distribution n<sub>+</sub>(x,y,z) in this way:
<pre>
n<sub>+</sub>(x,y,z) = n<sub>bulk</sub>   if <b>slabzbeg</b> &le z &le <b>slabzend</b>
= 0      otherwise,
</pre>
so the positive charge density is invariant along the xy plane as well as the electrostatic potential generated by it."""
},
'smdelta': {
'definition': "SMeared DELTA function ",
'section': "varrf",
'category': """<a href="../users/abinit_help.html#respfn">RESPFN</a>  """,
'vartype': "integer parameter ",
'default': "0",
'text': """When <b>smdelta</b> in non-zero, it will trigger the calculation of the imaginary part of the second-order electronic eigenvalues, which can be related to the electronic lifetimes. The delta function is evaluated using:
<br>
<ul>
<li> when <b>smdelta</b> == 1, Fermi-Dirac smearing : 0.25_dp/(cosh(xx/2.0_dp)**2 </li>
<li> when <b>smdelta</b> == 2, Cold smearing by Marzari using the parameter a=-.5634 (minimization of the bump): exp(-xx2)/sqrt(pi) * (1.5d0+xx*(-a*1.5d0+xx*(-1.0d0+a*xx))) </li>
<li> when <b>smdelta</b> == 3, Cold smearing by Marzari using the parameter a=-.8165 (monotonic function in the tail): as 2 but different a </li>
<li> when <b>smdelta</b> == 4, Smearing of Methfessel and Paxton (PRB40,3616(1989)) with Hermite polynomial of degree 2, corresponding to "Cold smearing" of N. Marzari with a=0 (so, same smeared delta function as smdelta=2, with different a). </li>
<li> when <b>smdelta</b> == 5, Gaussian smearing :  1.0d0*exp(-xx**2)/sqrt(pi) </li>
</ul>
"""
},
'so_psp': {
'definition': "Spin-Orbit treatment for each PSeudoPotential ",
'section': "vargs",
'category': " ",
'vartype': """integer array <b>so_psp</b>(<a href="vargs.html#npsp">npsp</a>) """,
'default': """<a href="vargs.html#npsp">npsp</a>*1 """,
'text': """For each type of atom (each pseudopotential), specify
the treatment of spin-orbit interaction (if <a href="vargs.html#nspinor">nspinor</a>==2).
<br> If 0 : no spin-orbit interaction, even if <a href="vargs.html#nspinor">nspinor</a>=2
<br> If 2 : treat spin-orbit in the HGH form
(not allowed for all pseudopotentials)
<br> If 3 : treat spin-orbit in the HFN form
(not allowed for all pseudopotentials)
<p>Also, <b>so_psp</b>=1 lead automatically to treatments 0, 2, or 3 according
to the data contained in the pseudopotential file
(0= there is no spin-orbit information in the psp file;
2= the spin-orbit information is of the HGH form;
3= the spin-orbit information is of the HFN form )
So, for typical usage, one will need only the values 0 or 1
<p>Note that if <a href="vargs.html#nspinor">nspinor</a>==1, the spin-orbit cannot be treated
anyhow, so the value of <b>so_psp</b> is irrelevant.
<p>Prior to v5.4, the input variable <b>so_typat</b>
was used, in place of <b>so_psp</b>. Because the values 0 and 1 have been switched
between <b>so_psp</b> and so_typat, it was dangerous to continue to allow the use of so_typat.
"""
},
'soenergy': {
'definition': "Scissor Operator ENERGY ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>, <a href="../users/abinit_help.html#energy">ENERGY</a> """,
'vartype': "real ",
'default': "0.0 ",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=3,4,99 that is, screening, sigma or Bethe-Salpeter calculations.
<p>
The Scissors operator energy added to the conductions states.
In some cases, it mimics a second iteration self-consistent GW calculation.
"""
},
'spbroad': {
'definition': "SPectral BROADening ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>, <a href="../users/abinit_help.html#energy">ENERGY</a> """,
'vartype': "real parameter ",
'default': "0.0 ",
'text': """Only relevant if <a href="vargs.html#optdriver">optdriver</a>=3 and
<a href="vargw.html#spmeth">spmeth</a>=2, that is, screening calculations
based on the spectral reprentation of the irreducible polarizability in which
the delta function is replaced by the gaussian approximant whose standard deviation
is given by <b>spbroad</b>.

"""
},
'spgaxor': {
'definition': "SPace Group : AXes ORientation ",
'section': "vargeo",
'category': """<a href="../users/abinit_help.html#symmetriser">SYMMETRISER</a>  """,
'vartype': " integer parameter",
'default': "1.",
'text': """It is taken into account only when <a href="vargeo.html#spgroup">spgroup</a>/=0; it
allows one to define the axes orientation
for the specific space groups for which this is needed.
Trigonal groups (number 146,148,155,160,161,166,167):
<ul>
<li> 1 represents the hexagonal axes</li>
<li> 2 represents the rhombohedral axes</li>
</ul>
Orthorhombic space groups :
there are six possibilities corresponding to the possible
axes permutations
<ul>
<li> 1 abc -> abc</li>
<li> 2 abc -> cab</li>
<li> 3 abc -> bca</li>
<li> 4 abc -> acb</li>
<li> 5 abc -> bac</li>
<li> 6 abc -> cba</li>
</ul>
Monoclinic : there are 3 or 9 possibilities depending on
the space group. See the space group
<a href="../users/spacegrouphelpfile.html">help file</a> for details.
In the log/output file the notation used to describe the
monoclinic groups is for example:<br>
15:c1, A2/a_c = C2/c<br>
where,
<ul>
<li>15 represents the space group number,</li>
<li>c1 the orientation as it appears on the web page, </li>
<li>A is the real Bravais type lattice, </li>
<li>2/a the existent symmetry elements, </li>
<li>_c marks the orientation of the two-fold axis or of the mirror plane, </li>
<li>C2/c represents the parent space group.</li>
</ul>
"""
},
'spgorig': {
'definition': "SPace Group : ORIGin  ",
'section': "vargeo",
'category': """<a href="../users/abinit_help.html#symmetriser">SYMMETRISER</a>  """,
'vartype': "integer parameter ",
'default': "1.",
'text': """Gives the choice of origin for the axes system,
taken into account only when <a href="vargeo.html#spgroup">spgroup</a>/=0,
<br>It is defined according to the origin choice in the
International Tables of Crystallography.
<br>It applies only to the space groups 48, 50, 59, 70, 85, 86, 88, 125,
126, 129, 130, 133, 134, 137, 141, 142, 201, 203, 222, 224, 227, 228.
<br>For details see the space group
<a href="../users/spacegrouphelpfile.html">help file</a>."""
},
'spgroup': {
'definition': "SPace GROUP number ",
'section': "vargeo",
'category': """<a href="../users/abinit_help.html#symmetriser">SYMMETRISER</a>  """,
'vartype': "integer parameter ",
'default': "0.",
'text': """Gives the number of the space group.
<br>If <b>spgroup</b> is 0, the code assumes that all the symmetries
are input through the <a href="varbas.html#symrel">symrel</a> matrices and the <a href="varbas.html#tnons">tnons</a>
vectors, or obtained from the symmetry finder (the default when
<a href="varbas.html#nsym">nsym</a>==0).
<br>It should be between 1 and 230.  This option can be
used to obtain all the atoms in the unit cell, starting
from the assymetric unit cell.
<br>The references for computing the symmetry corresponding to
the space groups are :
<ul>
<li>International Tables for Crystallography, 1983, Ed. Theo Hahn,
D. Reidel Publishing Company</li>
<li>The mathematical theory of symmetry in solids,
Representation theory for point groups and space groups, 1972,
C.J. Bradley and A.P. Cracknell, Clarendon Press, Oxford.</li></ul>
<br>For details see the space group
<a href="../users/spacegrouphelpfile.html">help file</a>.
"""
},
'spgroupma': {
'definition': "SPace GROUP number defining a MAgnetic space group",
'section': "vargeo",
'category': """<a href="../users/abinit_help.html#symmetriser">SYMMETRISER</a>, <a href="../users/abinit_help.html#not_internal">NOT INTERNAL</a> """,
'vartype': "integer parameter ",
'default': "0.",
'text': """This input variable might be used to define a Shubnikov
magnetic space group (anti-ferromagnetic space group). The user is advised to consult
"The mathematical theory of symmetry in solids,
Representation theory for point groups and space groups, 1972,
C.J. Bradley and A.P. Cracknell, Clarendon Press, Oxford."
<br>A Shubnikov type IV magnetic space group might be defined by its Fedorov space group
(set of spatial symmetries that do not change the magnetization), and an additional
magnetic space group number <b>spgroupma</b>.
<br>A Shubnikov type III magnetic space group might be defined by its Fedorov space group
(set of all spatial symmetries, irrespective of their magnetic action), and an additional
magnetic space group number <b>spgroupma</b>.
<br>For the additional number <b>spgroupma</b>, we follow the definition of Table 7.4 of the
above-mentioned Bradley and Cracknell textbook.
<br>Thus, one way to specify a Shubnikov IV magnetic space group, is to define both
<a href="vargeo.html#spgroup">spgroup</a> and <b>spgroupma</b>.
<br>For example, the group P2_1/c' has <a href="vargeo.html#spgroup">spgroup</a>=14
and <b>spgroupma</b>=78.
<br>Alternatively,
for Shubnikov IV magnetic groups, one might
define <a href="vargeo.html#spgroup">spgroup</a> and <a href="vargeo.html#genafm">genafm</a>.
For both the type III and IV, one might
define by hand the set of symmetries, using <a href="varbas.html#symrel">symrel</a>,
<a href="varbas.html#tnons">tnons</a> and <a href="vargs.html#symafm">symafm</a>
"""
},
'spinat': {
'definition': "SPIN for AToms ",
'section': "vargs",
'category': " ",
'vartype': """real array spinat(3,<a href="varbas.html#natom">natom</a>) or spinat(3,<a href="vargeo.html#natrd">natrd</a>) if the symmetriser is used  """,
'default': "0.0d0.",
'text': """Gives the initial electronic spin-magnetization
for each atom, in unit of h-bar/2.
<p>Note that if <a href="vargs.html#nspden">nspden</a>=2,
the z-component must be given
for each atom, in triplets (0 0 z-component).
<br>For example, the electron of an hydrogen atom
can be spin up (0 0 1.0) or spin down (0 0 -1.0).
<p>This value is only used to create
the first exchange and correlation potential,
and is not used anymore afterwards.
<br>It is not checked against the initial occupation numbers
<a href="vargs.html#occ">occ</a> for each spin channel.
<br>It is meant to give an easy way to break
the spin symmetry, and to allow
to find stable local spin fluctuations, for example :
antiferromagnetism, or the spontaneous spatial
spin separation of elongated H2 molecule.
<br><br><li>If the geometry builder is used, <b>spinat</b> will be related
to the preprocessed set of atoms, generated by the
geometry builder. The user must thus foresee the effect
of this geometry builder (see <a href="vargeo.html#objarf">objarf</a>).</li>
<br><li>If the geometry builder is not used, and the symmetries
are not specified by the user (<a href="varbas.html#nsym">nsym</a>=0),
spinat will be used, if present, to determine the anti-ferromagnetic
characteristics of the symmetry operations, see <a href="vargs.html#symafm">symafm</a>.<br>
In case of collinear antiferromagnetism
(<a href="varbas.html#nsppol">nsppol</a>=1,
<a href="vargs.html#nspinor">nspinor</a>=1,
<a href="vargs.html#nspden">nspden</a>=2),
these symmetries are used to symmetrize the density.<br>
In case of non-collinear magnetism
(<a href="varbas.html#nsppol">nsppol</a>=1,
<a href="vargs.html#nspinor">nspinor</a>=1,
<a href="vargs.html#nspden">nspden</a>=4),
they are also used to symmetrize the density.
In the latter case, this strongly constrains the magnetization (imposing its direction).
If the user want to let all degrees of freedom of the magnetization evolve, it is
then recommended to put <a href="varbas.html#nsym">nsym</a>=1.<br>
</li>
<br><li>If the symmetries are specified, and the irreducible set of atoms
is specified, the anti-ferromagnetic characteristics of the symmetry
operations <a href="vargs.html#symafm">symafm</a> will be used to generate
<b>spinat</b> for all the non-irreducible atoms.</li>
<br><li>In the case of PAW+U calculations using the <a href="varpaw.html#dmatpawu">dmatpawu</a>
initial occupation matrix, and if <a href="vargs.html#nspden">nspden</a>=4, <b>spinat</b> is
also used to determine the direction of the integrated magnetization matrix.</li>
"""
},
'spinmagntarget': {
'definition': "SPIN-MAGNetization TARGET ",
'section': "varff",
'category': " ",
'vartype': "real parameter  ",
'default': "-99.99d0",
'text': """This input variable is active only in the
<a href="varbas.html#nsppol">nsppol</a>=2 case.
If <b>spinmagntarget</b> is not the "magic" value of -99.99d0, the
spin-magnetization of the system will be fixed (or optimized, if it is not possible to impose it)
to the value of <b>spinmagntarget</b>.
<br>If <a href="varbas.html#occopt">occopt</a> is a metallic one, the
Fermi energies for spin up and spin down are adjusted to give the target
spin-polarisation (this is equivalent to an exchange splitting).
If <a href="varbas.html#occopt">occopt</a>=1 and <a href="varbas.html#nsppol">nsppol</a>=2,
the occupation numbers for spin up and spin down will be adjusted to give the required
spin-magnetization (occupation numbers are identical for all k-points, with <a href="varbas.html#occopt">occopt</a>=1)
<br>
If <b>spinmagntarget</b> is the default one, the spin-magnetization will not be constrained,
and will be determined self-consistently, by having the same spin up and spin down
Fermi energy in the metallic case, while for the other cases, there will be no spin-magnetization,
except for an odd number of electrons if <a href="varbas.html#occopt">occopt</a>=1 and <a href="varbas.html#nsppol">nsppol</a>=2.
<p>
Note : for the time being, only the spin down Fermi energy
is written out in the main output file. In the fixed
magnetic moment case, it differs from the
spin up Fermi energy.
"""
},
'spmeth': {
'definition': "SPectral METHod ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>  """,
'vartype': "integer ",
'default': "0  </b>",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=4, that is, sigma calculations.
<p>
The <b>spmeth</b> input variable defines the method used to calculate the irreducible
polarizability $\chi^{(0)}_{KS}$.
<p>
By default $\chi^{(0)}_{KS}$ is calculated employing the Adler-Wiser expression (<b>spmeth</b>=0)
with a CPU effort that scales linearly with the number of frequencies.
This approach is convenient when few frequencies are required, and is usually
used in conjunction with plasmon-pole models in which only one or two frequencies are calculated,
according to the value of <a href=vargw.html#ppmodel">ppmodel</a>.
<br>
Unfortunately a calculation based on the Adler-Wiser expression might be quite CPU demanding
if the matrix elements of the self-energy operator are evaluated by performing numerically
the convolution defining the self-energy.
The integrand function, indeed, has poles above and below the real axis, and
the screened interaction has to be evaluated on a dense frequency mesh in order to obtain accurate
results.
<p>
In the spectral method (<b>spmeth</b>=1 or 2) the irreducible polarizability is expressed as
the Hilbert transform of the imaginary part.
The advantage in using this approach consists in the fact that, once the spectral function is known,
the irreducible polarizability for an arbitrary frequency can be easily obtained through inexpensive
integrations. On the other hand an accurate evaluation of the imaginary part requires a dense
frequency mesh due to the presence of delta functions.
Two different approaches can be used to approximate these delta functions thus allowing the use
of affordable frequency grids.
<p>
Summarizing:

<ul>
<li> 0 =&gt; use Adler-Wiser expression to calculate $\chi^{(0)}_{KS}$  </li>
<li> 1 =&gt; use the spectral method approximating the delta function with a triangular approximant
as proposed in <B>REF TO BE ADDED</B> </li>
<li> 2 =&gt; use spectral method but approximating the delta function with a Taylor expansion of the exponential
as proposed in <B>REF TO BE ADDED</B> </li>
</ul>
"""
},
'spnorbscl': {
'definition': "SPin-ORBit SCaLing",
'section': "varpaw",
'category': "",
'vartype': "real ",
'default': "1.d0",
'text': """Relevant only when
<a href="varint.html#usepaw">usepaw</a>=1, and <a
href="varpaw.html#pawspnorb">pawspnorb</a> = 1 (or greater).<br>
Scaling of the spin-orbit interaction. The default values gives the first-principles value, while
other values are used for the analysis of the effect of the spin-orbit interaction,
but are not expected to correspond to any physical situation.

"""
},
'stmbias': {
'definition': "Scanning Tunneling Microscopy BIAS voltage ",
'section': "vargs",
'category': """<a href="../users/abinit_help.html#energy">ENERGY</a> """,
'vartype': "real parameter ",
'default': "0.00",
'text': """Gives, in Hartree, the
bias of the STM tip, with respect to the sample, in order to generate
the STM density map.
<br>Used with positive <a href="varbas.html#iscf">iscf</a>,
<a href="varbas.html#occopt">occopt</a>=7 (metallic, gaussian),
<a href="varbas.html#nstep">nstep</a>=1 ,
and positive <a href="varfil.html#prtstm">prtstm</a>, this
value is used to generate a charge density map from electrons
close to the Fermi energy, in a (positive or negative) energy range.
Positive <b>stmbias</b> will lead to the inclusion of occupied (valence) states only, while
negative <b>stmbias</b> will lead to the inclusion of unoccupied (conduction) states only.
<br>Can be specified in Ha (the default), Ry, eV or Kelvin, since
<b>stmbias</b> has the
'<a href="../users/abinit_help.html#dimensions">ENERGY</a>' characteristics
(0.001 Ha = 27.2113845 meV = 315.773 Kelvin).
With <a href="varbas.html#occopt">occopt</a>=7,
one has also to specify an independent broadening <a href="vargs.html#tsmear">tsmear</a>.
"""
},
'strfact': {
'definition': "STRess FACTor ",
'section': "varrlx",
'category': "",
'vartype': "real parameter ",
'default': "100.0 (Bohr^2)",
'text': """The stresses multiplied
by <b>strfact</b> will be treated like forces in the
process of optimization (<a href="varrlx.html#ionmov">ionmov</a>=2,
non-zero
<a href="varrlx.html#optcell">optcell</a>).
<br>
For example, the stopping criterion defined by
<a href="varrlx.html#tolmxf">tolmxf</a> relates to these scaled
stresses.
"""
},
'string_algo': {
'definition': "STRING method ALGOrithm ",
'section': "varrlx",
'category': "",
'vartype': "integer parameter",
'default': "<b>string_algo</b>=1",
'text': """Relevant only when <a href="varrlx.html#imgmov">imgmov</a>=2 (String Method).<br>
Gives the variant of the String Method method used.<br>
Possible values can be:<br>
<ul>
<li>0=&gt; <b>Original String Method</b>.<br>
NOT YET IMPLEMENTED<br>
<i>See: Phys. Rev. B 66, 052301 (2002)</i>
</li>
<br>
<li>1=&gt; <b>Simplified String Method</b> with parametrization by <b>equal arc length</b>.<br>
Instead of using the normal force (wr the band), the full force is used; the
reparametrization is enforced by keeping the points of the string equally spaced.<br>
<i>See: J. Chem. Phys. 126, 164103 (2007)</i>
</li>
<br>
<li>2=&gt; <b>Simplified String Method</b> with parametrization by <b>energy-weighted arc length</b>.<br>
A variant of the Simplified String Method (like 2-); the reparametrization is done by using
energy-weight arc-lengths, giving a finer distribution near the saddle point..<br>
<i>See: J. Chem. Phys. 126, 164103 (2007) and J. Chem. Phys. 130, 244108 (2009)</i>
</li>
</ul>
"""
},
'strprecon': {
'definition': "STRess PRECONditioner ",
'section': "varrlx",
'category': "",
'vartype': "real parameter ",
'default': "1.0 ",
'text': """This is a scaling factor to initialize the part of
the Hessian related to the treatment of the stresses (optimisation
of the unit cell). In case there is an instability, decrease the
default value, e.g. set it to 0.1 .
"""
},
'strtarget': {
'definition': "STRess TARGET ",
'section': "varrlx",
'category': "",
'vartype': "real array strtarget(6) ",
'default': "6*0.0 (Ha/Bohr**3)",
'text': """The components of the stress tensor must be stored
according to :
(1,1)-&gt;1 ; (2,2)-&gt;2 ; (3,3)-&gt;3 ; (2,3)-&gt;4 ; (3,1)-&gt;5 ;
(1,2)-&gt;6.
The conversion factor
between Ha/Bohr**3 and GPa is : 1 Ha/Bohr**3 = 29421.033d0 GPa.
<br>
Not used if <a href="varrlx.html#optcell">optcell</a>==0.
"""
},
'suskxcrs': {
'definition': "SUSceptibility times KXC treated in real space ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>  """,
'vartype': "integer",
'default': "0",
'text': """Only relevant for the ACFD calculation of total energies.
If <b>suskxcrs</b>=1, the XC kernel is not treated in
reciprocal space, but combined with the susceptibility (chi_0), to avoid
Kxc divergences where the density goes to zero (G. Onida &amp; M. Gatti !)
<p>
Not applicable for RPA (as there should be a Kxc present). Initially tested for ikhxc==2 (ALDA).
"""
},
'symafm': {
'definition': "SYMmetries, Anti-FerroMagnetic characteristics ",
'section': "vargs",
'category': " ",
'vartype': """integer array symafm(<a href="varbas.html#nsym">nsym</a>) """,
'default': """<a href="varbas.html#nsym">nsym</a>*1.""",
'text': """In case the material is magnetic (well, this is only interesting in the
case of antiferromagnetism, collinear or not), additional symmetries might appear, that
change the sign of the magnetization.
They have been introduced by Shubnikov (1951). They can be used by ABINIT
to decrease the CPU time, by using them to decrease the number of k-points.
<br> <b>symafm</b> should be set to +1 for all the usual symmetry operations,
that do not change the sign of the magnetization, while it should be
set to -1 for the magnetization-changing symmetries.
<br> If the symmetry operations are not specified by the user
in the input file, that is, if <a href="varbas.html#nsym">nsym</a>=0,
then ABINIT will use the values of <a href="vargs.html#spinat">spinat</a>
to determine the content of <b>symafm</b>.<br>
The symmetries found as "antiferro magnetic" (<b>symafm</b>=-1) are used to symmetrize density and magnetization in the following cases:<br>
- antiferromagnetism (<a href="varbas.html#nsppol">nsppol</a>=1,
<a href="vargs.html#nspinor">nspinor</a>=1,
<a href="vargs.html#nspden">nspden</a>=2)<br>
- non-collinear magnetism (<a href="varbas.html#nsppol">nsppol</a>=1,
<a href="vargs.html#nspinor">nspinor</a>=1,
<a href="vargs.html#nspden">nspden</a>=4)<br>
In other cases they are not used.
"""
},
'symchi': {
'definition': "SYMmetryze \chi_o ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>, GW """,
'vartype': "integer ",
'default': "1 ",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=3, that is, screening calculations.
<p>
The evaluation of the irreducible polarizability for a given q-point requires an
integration over the Brillouin zone (BZ) which is approximated by a discrete sum over k-points.
In principle the integrand function should be evaluated for each k-point in the BZ, however it is possible
to reduce the number of points to be explicitly considered by taking advantage of symmetry properties.
The development input variable <b>symchi</b> is used to choose between these two equivalent methods:

<ul>
<li> 0=&gt; the summation over k-points is performed considering ALL the points in the BZ (useful for testing and debugging).
<li> 1=&gt; the summation is restricted to the k-points belonging to the irreducible wedge
defined by the little group associated to the external vector q.
</ul>

"""
},
'symmorphi': {
'definition': "SYMMORPHIc symmetry operations ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>, GW """,
'vartype': "integer parameter",
'default': "1",
'text': """With <b>symmorphi</b>=1, symmetry operations with a non-symmorphic vector are allowed.
With <b>symmorphi</b>=0, they are not allowed.
In the latter case, if the symmetry operations are specified in the input file, the code
will stop and print an error message if a non-symmorphic vector is encountered.
By contrast, if the symmetry operations are to be determined automatically
(if <a href="varbas.html#nsym">nsym</a>=0), then the set of symmetries will
not include the non-symmorphic operations.
<p>
Note : this feature exist because in a previous status of the GW calculations, non-symmorphic
symmetry operations could not be exploited. Thus, the k points were restricted
to the IBZ. In order to prepare GW calculations, and to perform GW calculations,
<b>symmorphi</b>=0 was to be used, together with <a href="varbas.html#nsym">nsym</a>=0.
"""
},
'symrel': {
'definition': "SYMmetry in REaL space  ",
'section': "varbas",
'category': " ",
'vartype': """integer array symrel(3,3,<a href="varbas.html#nsym">nsym</a>)  """,
'default': "the identity matrix for one symmetry.",
'text': """Gives "<a href="varbas.html#nsym">nsym</a>" 3x3 matrices
expressing space group symmetries in terms of their action
on the direct (or real) space primitive translations.
<br>It turns out that these can always be expressed as integers.
<br>Always give the identity matrix even if no other symmetries
hold, e.g.
<b>symrel</b> 1 0 0 0 1 0 0 0 1
<br>Also note that for this array as for all others the array
elements are filled in a columnwise order as is usual for
Fortran.
<br>The relation between the above symmetry matrices <b>symrel</b>,
expressed in the basis of primitive translations, and
the same symmetry matrices expressed in cartesian coordinates,
is as follows.  Denote the matrix whose columns are the
primitive translations as R, and denote the cartesian
symmetry matrix as S.  Then <br>
<b>symrel</b> = R(inverse) * S * R <br>
where matrix multiplication is implied.
<br>When the symmetry finder is used (see <a href="varbas.html#nsym">nsym</a>), <b>symrel</b>
will be computed automatically."""
},
'symsigma': {
'definition': "SYMmetrization of SIGMA matrix elements ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a> """,
'vartype': "integer ",
'default': "0",
'text': """Only relevant if <a href="vargs.html#optdriver">optdriver</a>=4 that is sigma calculations.
<p>
This option is used to switch on the symmetrization of the self-energy matrix elements (<a href="vargw.html#symsigma">symsigma</a>=1).
In this case the BZ integration defining the self-energy matrix elements is
reduced to an appropriate irreducible wedge defined by the point group of the wave-vector k
specified in the <a href="vargw.html#kptgw">kptgw</a> list.
<p>
The symmetrized expression leads to a considerable speedup of the run but, unfortunately,
this option is not yet compatible with self-consistent GW calculations
(see <a href="vargw.html#gwcalctyp">gwcalctyp</a>).
<p>
The algorithm implemented in <a href="vargw.html#symsigma">symsigma</a>=1
constructs a symmetric invariant for the diagonal matrix elements of the self-energy
by simply averaging the GW results within the degenerate subspace.

Therefore particular care has to be taken in the presence of accidental degeneracies.
since GW calculations performed with <a href="vargw.html#symsigma">symsigma</a>=1 won't be able to remove
the initial accidental degeneracy.

"""
},
'td_maxene': {
'definition': "Time-Dependent dft : MAXimal kohn-sham ENErgy difference ",
'section': "varrf",
'category': "TDDFT  ",
'vartype': "real parameter ",
'default': "0.0",
'text': """The Matrix to be diagonalized in the Casida framework
(see "Time-Dependent Density Functional Response Theory of Molecular
systems: Theory, Computational Methods, and Functionals", by M.E. Casida,
in Recent Developments and Applications of Modern Density Functional
Theory, edited by J.M. Seminario (Elsevier, Amsterdam, 1996).)
is a NxN matrix, where, by default, N is the product of
the number of occupied states by the number of unoccupied states.
<br>
The input variable <b>td_maxene</b> allows to diminish N : it selects
only the pairs of occupied and unoccupied states for which the
Kohn-Sham energy difference is less than <b>td_maxene</b>.
The default value 0.0 means that all pairs are taken into account.
<br>See <a href="varrf.html#td_mexcit">td_mexcit</a> for an alternative
way to decrease N.
"""
},
'td_mexcit': {
'definition': "Time-Dependent dft : Maximal number of EXCITations ",
'section': "varrf",
'category': "TDDFT  ",
'vartype': "real parameter ",
'default': "0.",
'text': """The Matrix to be diagonalized in the Casida framework
(see "Time-Dependent Density Functional Response Theory of Molecular
systems: Theory, Computational Methods, and Functionals", by M.E. Casida,
in Recent Developments and Applications of Modern Density Functional
Theory, edited by J.M. Seminario (Elsevier, Amsterdam, 1996).)
is a NxN matrix, where, by default, N is the product of
the number of occupied states by the number of unoccupied states.
<br>
The input variable <b>td_mexcit</b> allows to diminish N : it selects
the first <b>td_mexcit</b> pairs of occupied and unoccupied states, ordered
with respect to increasing Kohn-Sham energy difference.
However, when <b>td_mexcit</b> is zero, all pairs are allowed.
<br>See <a href="varrf.html#td_mexcit">td_maxene</a> for an alternative
way to decrease N.
"""
},
'tfkinfunc': {
'definition': "Thomas-Fermi KINetic energy FUNCtional ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>  """,
'vartype': "integer",
'default': "0",
'text': """<ul>
<li>
<b>tfkinfunc</b>=1 : Thomas-Fermi
kinetic functional (explicit functional of the density)
is used instead of Kohn-Sham kinetic energy functional (implicit functional of the density
through Kohn-Sham wavefunctions).
<br></li>
<li> <b>tfkinfunc</b>=2 : the Recursion Method is used in order to compute
electronic density, entropy, Fermi energy and eigenvalues energy. This
method computes the density without computing any orbital, is
efficient at high temperature, with a efficient parallelization
(almost perfect scalability). When that option is in use, the
<a href="varbas.html#ecut">ecut</a> input variable is no longer a
convergence parameter ;
<a href="vargs.html#ngfft">ngfft</a> becomes
the main convergence parameter: you should adapt ecut for the ngfft
grid you need (it is not yet automatically computed). Other
convergence parameter are for the energetic values:
<a href="vardev.html#recnrec">recnrec</a>, <a href="vardev.html#recptrott">recptrott</a>,
<a href="vardev.html#recnpath">recnpath</a>.
Since the convergence of the self-consistent cycle
is determined directly by the convergence of the density:
<a href="varbas.html#toldfe">toldfe</a>,
<a href="varbas.html#toldff">toldff</a>
<a href="varbas.html#tolrff">tolrff</a>,
<a href="varbas.html#tolvrs">tolvrs</a>,
<a href="varbas.html#tolwfr">tolwfr</a> are not used, and are replaced by
<a href="vardev.html#rectolden">rectolden</a>; the energetic values, except for the fermi energy, are only
computed during the latest SFC cycle : the output file will show a
jump of the total energy at the end, but it is not because of a bad
convergence behavior. Computational speed can be improved by the use
of  <a href="vardev.html#recrcut">recrcut</a> and  <a href="vardev.html#recgratio">recgratio</a>.
The recursion method has not be tested in the case of non cubic cell
or with the use of symmetries.
In the recursion method the following variables are set to:
<br> <a href="vardev.html#useylm">useylm</a>=1,
<br> <a href="varint.html#userec">userec</a>=1.
</ul>
</li>
"""
},
'timopt': {
'definition': "TIMing OPTion ",
'section': "vargs",
'category': """<a href="../users/abinit_help.html#no_multi">NO MULTI</a> """,
'vartype': "integer parameter ",
'default': "1 for sequential code, 0 for parallel code.",
'text': """This input variable allows to modulate the use of the timing routines.
<p>
If 0 => as soon as possible, suppresses all calls to timing routines
<br>
If 1 => usual timing behaviour, with short analysis, appropriate
for sequential execution
<br>
If 2 => close to <b>timopt</b>=1, except that the analysis routine
does not time the timer, appropriate for parallel execution.
<br>
If 3 => close to <b>timopt</b>=1, except that the different parts of the lobpcg routine are timed in detail.
<br>
If 4 => close to <b>timopt</b>=1, except that the different parts of the lobpcg routine are timed in detail.
A different splitting of lobpcg than for <b>timopt</b>=-3 is provided.
<br>
If -1 => a full analysis of timings is delivered
<br>
If -2 => a full analysis of timings is delivered,
except timing the timer
<br>
If -3 => a full analysis of timings is delivered, including the detailed timing of the different parts of the lobpcg routine.
(this takes time, and is discouraged for too small runs - the timing would take more time than the run !). The timer is timed.

<br>
If -4 => a full analysis of timings is delivered, including the detailed timing of the different parts of the lobpcg routine.
A different splitting of lobpcg than for <b>timopt</b>=-3 is provided
(this takes time, and is discouraged for too small runs - the timing would take more time than the run !). The timer is timed.
The sum of the independent parts is closer to 100% than for <b>timopt</b>=-3.
"""
},
'tl_nprccg': {
'definition': "TaiL maximum Number of PReConditionner Conjugate Gradient iterations ",
'section': "vargs",
'category': " ",
'vartype': "integer parameter ",
'default': "30.",
'text': """This variable is the same than <a href="#wvl_nprccg">wvl_nprccg</a>
but for the preconditionner iterations during the tail
corrections (see <a href="#tl_radius">tl_radius</a>).

TO BE IMPROVED : all tl_* and wvl_* variables should contain a link to a tutorial, David Waroquiers 090831.
"""
},
'tl_radius': {
'definition': "TaiL expansion RADIUS ",
'section': "vargs",
'category': """<a href="../users/abinit_help.html#length">LENGTH</a> """,
'vartype': "real parameter ",
'default': "0.0d0.",
'text': """In the wavelet computation case, the linkage between the grid and the
free boundary conditions can be smoothed using an exponential
decay. This means a correction on the energy at the end on each
wavefunction optimisation run. If this parameter is set to zero,
no tail computation is done. On the contrary, put it to a
positive value makes the tail correction available. The value
correspond to a length in atomic units being the spacial expansion
with the exponential decay around the grid.
"""
},
'tnons': {
'definition': "Translation NON-Symmorphic vectors ",
'section': "varbas",
'category': " ",
'vartype': """real array tnons(3,<a href="varbas.html#nsym">nsym</a>)  """,
'default': "",
'text': """Gives the (nonsymmorphic) translation
vectors associated with the symmetries expressed
in "<a href="varbas.html#symrel">symrel</a>".
<br>These may all be 0, or may be fractional (nonprimitive)
translations expressed relative to the real space
primitive translations (so, using the "reduced" system
of coordinates, see "<a href="varbas.html#xred">xred</a>").
If all elements of the space
group leave 0 0 0 invariant, then these are all 0.
<br>When the symmetry finder is used (see <a href="varbas.html#nsym">nsym</a>), <b>tnons</b>
is computed automatically."""
},
'toldfe': {
'definition': "TOLerance on the DiFference of total Energy ",
'section': "varbas",
'category': """<a href="../users/abinit_help.html#energy">ENERGY</a>  """,
'vartype': "real parameter  ",
'default': "0.0 (stopping condition ignored)",
'text': """Sets a tolerance for absolute differences
of total energy that, reached TWICE successively,
will cause one SCF cycle to stop (and ions to be moved).
<br>Can be specified in Ha (the default), Ry, eV or Kelvin, since
<b>toldfe</b> has the
'<a href="../users/abinit_help.html#dimensions">ENERGY</a>' characteristics.
(1 Ha=27.2113845 eV)
<br>If set to zero, this stopping condition is ignored.
<br>Effective only when SCF cycles are done (<a href="varbas.html#iscf">iscf</a>>0).
<br>Because of machine precision, it is not worth to try
to obtain differences in energy that are smaller
than about 1.0d-12 of the total energy.
To get accurate stresses may be quite demanding.
<br>When the geometry is optimized (relaxation of atomic positions or primitive vectors), the use of
<b>toldfe</b> is to be avoided. The use of <a href="varbas.html#toldff">toldff</a> or
<a href="varbas.html#tolrff">tolrff</a> is by far preferable, in order to have a handle on the
geometry characteristics. When all forces vanish by symmetry (e.g. optimization of the lattice parameters
of a high-symmetry crystal), then place <b>toldfe</b> to 1.0d-12, or use (better) <a href="varbas.html#tolvrs">tolvrs</a>.
<br>Since <b>toldfe</b>, <a href="varbas.html#toldff">toldff</a>,
<a href="varbas.html#tolrff">tolrff</a>,
<a href="varbas.html#tolvrs">tolvrs</a> and <a href="varbas.html#tolwfr">tolwfr</a>
are aimed at the same goal (causing the SCF cycle to stop),
they are seen as a unique input variable at reading. Hence, it is forbidden that two of these input variables
have non-zero values for the same dataset, or generically (for all datasets).
However, a non-zero value for one such variable for one dataset will have precedence on the non-zero value for another
input variable defined generically.

"""
},
'toldff': {
'definition': "TOLerance on the DiFference of Forces ",
'section': "varbas",
'category': " ",
'vartype': "real parameter  ",
'default': "0.0 (stopping condition ignored)",
'text': """Sets a tolerance for differences of forces
(in hartree/Bohr) that, reached TWICE successively,
will cause one SCF cycle to stop (and ions to be moved).
<br>If set to zero, this stopping condition is ignored.
<br>Effective only when SCF cycles are done (<a href="varbas.html#iscf">iscf</a>>0).
This tolerance applies to any particular cartesian
component of any atom, INCLUDING fixed ones.
This is to be used when trying to equilibrate a
structure to its lowest energy configuration (<a href="varrlx.html#ionmov">ionmov</a>=2),
or in case of molecular dynamics (<a href="varrlx.html#ionmov">ionmov</a>=1)
<br>A value ten times smaller
than <a href="varrlx.html#tolmxf">tolmxf</a> is suggested (for example 5.0d-6 hartree/Bohr).
<br>This stopping criterion is not allowed for RF calculations.
<br>Since <b>toldfe</b>, <a href="varbas.html#toldff">toldff</a>,
<a href="varbas.html#tolrff">tolrff</a>,
<a href="varbas.html#tolvrs">tolvrs</a> and <a href="varbas.html#tolwfr">tolwfr</a>
are aimed at the same goal (causing the SCF cycle to stop),
they are seen as a unique input variable at reading. Hence, it is forbidden that two of these input variables
have non-zero values for the same dataset, or generically (for all datasets).
However, a non-zero value for one such variable for one dataset will have precedence on the non-zero value for another
input variable defined generically.
"""
},
'tolimg': {
'definition': "TOLerance on the mean total energy for IMaGes ",
'section': "varrlx",
'category': """<a href="../users/abinit_help.html#energy">ENERGY</a> """,
'vartype': "real parameter ",
'default': "5.0d-5 hartree",
'text': """Sets a maximal absolute energy tolerance (in hartree, averaged over dynamic images)
below which iterations on images (the one governed by the <a
href="varrlx.html#ntimimage">ntimimage</a> input variable) will stop.
<br>
This is to be used when trying to optimize a
population of structures to their lowest energy configuration,
taking into account the particular algorithm defined by <a href="varrlx.html#imgmov">imgmov</a>
<br>
A value of about 5.0d-5 hartree or smaller
is suggested (this corresponds to about 3.7d-7 eV).
<br>
No meaning for RF calculations.
"""
},
'tolmxf': {
'definition': "TOLerance on the MaXimal Force ",
'section': "varrlx",
'category': "",
'vartype': "real parameter ",
'default': "5.0d-5 hartree/Bohr.",
'text': """Sets a maximal absolute force tolerance
(in hartree/Bohr) below which
BFGS structural relaxation iterations will stop.
<br>
Can also control tolerance on stresses, when <a
href="varrlx.html#optcell">optcell</a>/=0,
using the conversion factor <a href="varrlx.html#strfact">strfact</a>.
This tolerance applies to any particular cartesian
component of any atom, excluding fixed ones.
See the parameter <a href="varrlx.html#ionmov">ionmov</a>.
<br>
This is to be used when trying to equilibrate a
structure to its lowest energy configuration (<a
href="varrlx.html#ionmov">ionmov</a>=2).
<br>
A value of about 5.0d-5 hartree/Bohr or smaller
is suggested (this corresponds to about 2.5d-3 eV/Angstrom).
<br>
No meaning for RF calculations.
"""
},
'tolrde': {
'definition': "TOLerance on the Relative Difference of Eigenenergies ",
'section': "vardev",
'category': " ",
'vartype': "real parameter  ",
'default': "0.005",
'text': """Sets a tolerance for the ratio of differences of eigenenergies
in the line minisation conjugate-gradient algorithm. It compares the
decrease of the eigenenergy due to the last line minimisation, with the
one observed for the first line minimisation.
When the ratio is lower than <b>tolrde</b>,
the next line minimisations are skipped.
<br>The number of line minimisations is limited by
<a href="vargs.html#nline">nline</a> anyhow.
<br>This stopping criterion is present for both GS and RF calculations.
In RF calculations, <b>tolrde</b> is actually doubled before comparing with the above-mentioned
ratio, for historical reasons."""
},
'tolrff': {
'definition': "TOLerance on the Relative diFference of Forces ",
'section': "varbas",
'category': " ",
'vartype': "real parameter  ",
'default': "0.0 (stopping condition ignored)",
'text': """Sets a tolerance for the ratio of differences of forces
(in hartree/Bohr) to maximum force, that, reached TWICE successively,
will cause one SCF cycle to stop (and ions to be moved) : diffor < tolrff * maxfor.
<br>If set to zero, this stopping condition is ignored.
<br>Effective only when SCF cycles are done (<a href="varbas.html#iscf">iscf</a>>0).
This tolerance applies to any particular cartesian
component of any atom, INCLUDING fixed ones.
This is to be used when trying to equilibrate a
structure to its lowest energy configuration (<a href="varrlx.html#ionmov">ionmov</a>=2),
or in case of molecular dynamics (<a href="varrlx.html#ionmov">ionmov</a>=1)
<br>A value of 0.02 is suggested.
<br>This stopping criterion is not allowed for RF calculations.
<br>Since <b>toldfe</b>, <a href="varbas.html#toldff">toldff</a>,
<a href="varbas.html#tolrff">tolrff</a>,
<a href="varbas.html#tolvrs">tolvrs</a> and <a href="varbas.html#tolwfr">tolwfr</a>
are aimed at the same goal (causing the SCF cycle to stop),
they are seen as a unique input variable at reading. Hence, it is forbidden that two of these input variables
have non-zero values for the same dataset, or generically (for all datasets).
However, a non-zero value for one such variable for one dataset will have precedence on the non-zero value for another
input variable defined generically.
"""
},
'tolsym': {
'definition': "TOLERANCE for SYMmetries ",
'section': "vargeo",
'category': "",
'vartype': "real ",
'default': "1.e-8 ",
'text': """Gives the tolerance on the atomic positions (reduced coordinates), primitive
vectors, or magnetization, to be considered equivalent, thanks to symmetry operations.
This is used in the recognition of the set of symmetries of the system,
or the application of the symmetry operations to generate from a reduced set of atoms,
the full set of atoms. Note that a value larger than 0.01 is considered to be unacceptable.
<br>
Note : ABINIT needs the atomic positions to be symmmetric to each others within 1.e-8 .
If <b>tolsym</b> is set to a larger value than 1.e-8, then the input atomic coordinates
will be automatically symmetrized by the symmetry operations that will have been found.
"""
},
'tolvrs': {
'definition': "TOLerance on the potential V(r) ReSidual ",
'section': "varbas",
'category': " ",
'vartype': "real parameter  ",
'default': "0.0 (stopping condition ignored)",
'text': """Sets a tolerance for potential
residual that, when reached, will cause one SCF cycle
to stop (and ions to be moved).
<br>If set to zero, this stopping condition is ignored.
<br>Effective only when SCF cycles are done (<a href="varbas.html#iscf">iscf</a>>0).
<br>To get accurate stresses may be quite demanding.
<p>Additional explanation : the residual of the potential is the difference between the
input potential and the output potential, when the latter is obtained from the density
determined from the eigenfunctions of the input potential. When the self-consistency
loop is achieved, both input and output potentials must be equal, and the residual
of the potential must be zero. The tolerance on the
potential residual is imposed by first subtracting the mean of the residual of the potential
(or the trace of the potential matrix, if the system is spin-polarized),
then summing the square of this function over all FFT grid points. The result should be
lower than <b>tolvrs</b>.
<br>Since <b>toldfe</b>, <a href="varbas.html#toldff">toldff</a>,
<a href="varbas.html#tolrff">tolrff</a>,
<a href="varbas.html#tolvrs">tolvrs</a> and <a href="varbas.html#tolwfr">tolwfr</a>
are aimed at the same goal (causing the SCF cycle to stop),
they are seen as a unique input variable at reading. Hence, it is forbidden that two of these input variables
have non-zero values for the same dataset, or generically (for all datasets).
However, a non-zero value for one such variable for one dataset will have precedence on the non-zero value for another
input variable defined generically.
"""
},
'tolwfr': {
'definition': "TOLerance on WaveFunction squared Residual ",
'section': "varbas",
'category': " ",
'vartype': "real parameter  ",
'default': "0.0d0 (stopping criterion ignored)",
'text': """The signification of this tolerance depends on
the basis set. In plane waves, it gives a convergence tolerance for the
largest squared "residual" (defined below) for any
given band.  The squared residual is:<br>
<pre>
&lt; nk|(H-E)<sup>2</sup>|nk &gt;,    E = &lt; nk|H|nk &gt;
</pre>
<br>which clearly is nonnegative and goes to 0 as
the iterations converge to an eigenstate.
With the squared residual expressed in
Hartrees<sup>2</sup> (Hartrees squared), the largest squared
residual (called residm) encountered over all bands
and k points must be less than <b>tolwfr</b> for iterations
to halt due to successful convergence.
<br>Note that if <a href="varbas.html#iscf">iscf</a>>0, this criterion should be replaced
by those based on <a href="varbas.html#toldfe">toldfe</a> (preferred for <a href="varrlx.html#ionmov">ionmov</a>==0),
<a href="varbas.html#toldff">toldff</a>
<a href="varbas.html#tolrff">tolrff</a> (preferred for <a href="varrlx.html#ionmov">ionmov</a>/=0), or
<a href="varbas.html#tolvrs">tolvrs</a> (preferred for theoretical reasons!).
<br>When <b>tolwfr</b> is 0.0, this criterion is ignored,
and a finite value of <a href="varbas.html#toldfe">toldfe</a>, <a href="varbas.html#toldff">toldff</a>
or <a href="varbas.html#tolvrs">tolvrs</a> must be specified.
This also imposes a restriction
on taking an ion step; ion steps are not permitted
unless the largest squared residual is less than
<b>tolwfr</b>, ensuring accurate forces.
<br>To get accurate stresses may be quite demanding.
<br>
Note that the preparatory GS calculations
before a RF calculations must be highly converged.
<br>Typical values for these preparatory runs are <b>tolwfr</b>
between 1.0d-16 and 1.0d-22.
<p> Note that <b>tolwfr</b> is often used in the test cases, but this is
<em>tolwfr</em> purely for historical reasons :
except when <a href="varbas.html#iscf">iscf</a>&lt;0, other critera
should be used.</p>
<p>In the wavelet case (see <a href="varbas.html#usewvl">usewvl</a> =
1), this criterion is the favoured one. It is based on the
norm 2 of the gradient of the wavefunctions. Typical values
range from 5*10<sup>-4</sup> to 5*10<sup>-5</sup>.</p>
<br>Since <b>toldfe</b>, <a href="varbas.html#toldff">toldff</a>,
<a href="varbas.html#tolrff">tolrff</a>,
<a href="varbas.html#tolvrs">tolvrs</a> and <a href="varbas.html#tolwfr">tolwfr</a>
are aimed at the same goal (causing the SCF cycle to stop),
they are seen as a unique input variable at reading. Hence, it is forbidden that two of these input variables
have non-zero values for the same dataset, or generically (for all datasets).
However, a non-zero value for one such variable for one dataset will have precedence on the non-zero value for another
input variable defined generically.
"""
},
'tphysel': {
'definition': "Temperature (PHYSical) of the ELectrons ",
'section': "vargs",
'category': """<a href="../users/abinit_help.html#energy">ENERGY</a> """,
'vartype': "real parameter ",
'default': "0.00",
'text': """Gives, in Hartree, the physical temperature of the
system, in case <a href="varbas.html#occopt">occopt</a>=4, 5, 6, or 7.
<br>Can be specified in Ha (the default), Ry, eV or Kelvin, since
<b>ecut</b> has the
'<a href="../users/abinit_help.html#dimensions">ENERGY</a>' characteristics
(0.001 Ha = 27.2113845 meV = 315.773 Kelvin).
One has to specify an independent broadening <a href="vargs.html#tsmear">tsmear</a>.
The combination of the two parameters
<b>tphysel</b> and <a href="vargs.html#tsmear">tsmear</a> is described
in a paper by M. Verstraete and X. Gonze, Phys. Rev. B 65, 035111 (2002).
Note that the signification of the entropy is modified with respect
to the usual entropy. The choice has been made to use
<a href="vargs.html#tsmear">tsmear</a> as a prefactor of the entropy,
to define the entropy contribution to the free energy.
"""
},
'tsmear': {
'definition': "Temperature of SMEARing ",
'section': "vargs",
'category': """<a href="../users/abinit_help.html#energy">ENERGY</a> """,
'vartype': "real parameter ",
'default': "0.04",
'text': """Gives the broadening of occupation
numbers <a href="vargs.html#occ">occ</a>, in the metallic cases
(<a href="varbas.html#occopt">occopt</a>=3, 4, 5, 6 and 7).
Can be specified in Ha (the default), eV, Ry, or Kelvin, since
<b>tsmear</b> has the
'<a href="../users/abinit_help.html#dimensions">ENERGY</a>' characteristics
(0.001 Ha = 27.2113845 meV = 315.773 Kelvin).
<br>Default is 0.04 Ha. This should be OK for a free-electron
metal like Al. For d-band metals, use 0.01 Ha.
<br>Always check the convergence of the calculation
with respect to this parameter, and simultaneously,
with respect to the sampling of k-points (see <a href="varbas.html#nkpt">nkpt</a>)
<br>If <a href="varbas.html#occopt">occopt</a>=3, <b>tsmear</b> is the
physical temperature, as the broadening is based on Fermi-Dirac statistics.
However,
if <a href="varbas.html#occopt">occopt</a>=4, 5, 6, or 7,
the broadening is not based on Fermi-Dirac statistics, and
<b>tsmear</b> is only a convergence parameter. It is still possible
to define a physical temperature, thanks to the input variable
<a href="vargs.html#tphysel">tphysel</a>. See the paper
by M. Verstraete and X. Gonze, Phys. Rev. B (2002).
"""
},
'typat': {
'definition': "TYPE of atoms ",
'section': "varbas",
'category': " ",
'vartype': """integer array typat(<a href="varbas.html#natom">natom</a>) (or : typat(<a href="vargeo.html#natrd">natrd</a>), if the geometry builder is used) """,
'default': """1 (for <a href="varbas.html#natom">natom</a>=1)""",
'text': """Array giving an integer label to every atom in the unit
cell to denote its type. <br>The different types of atoms
are constructed from the pseudopotential files.
There are at most <a href="varbas.html#ntypat">ntypat</a> types of atoms.
<br> As an example, for BaTiO3, where the pseudopotential for Ba is number 1,
the one of Ti is number 2, and the one of O is number 3, the actual
value of the <b>typat</b> array might be :
<pre>
typat 1 2 3 3 3
</pre>
<br>The array <b>typat</b> has to agree with the actual locations
of atoms given in <a href="varbas.html#xred">xred</a> , <a href="varbas.html#xcart">xcart</a> or
<a href="varbas.html#xangst">xangst</a>, and the input
of pseudopotentials has to be ordered to agree with the
atoms identified in <b>typat</b>.
<br>The nuclear charge of the
elements, given by the array <a href="varbas.html#znucl">znucl</a>, also must agree with
the type of atoms designated in "<b>typat</b>".
<br>The array <b>typat</b> is
not constrained to be increasing. An
internal representation of the list of atoms,
deep in the code (array atindx), groups the atoms of same type
together. This should be transparent to the
user, while keeping efficiency.
"""
},
'udtset': {
'definition': "Upper limit on DaTa SETs ",
'section': "varbas",
'category': " ",
'vartype': "integer array udtset(2)",
'default': "No Default (since it is not used when it is not defined).",
'text': """Used to define the set of indices in the multi-data set
mode, when a double loop is needed (see later).
<br>The values of <b>udtset</b>(1) must be between 1 and 999,
the values of <b>udtset</b>(2) must be between 1 and 9, and their
product must be equal to <a href="varbas.html#ndtset">ndtset</a>.
<br>
The values of <a href="varbas.html#jdtset">jdtset</a> are obtained by
looping on the two indices defined by <b>udtset</b>(1) and  <b>udtset</b>(2) as follows :
<pre>
do i1=1,intarr(1)
do i2=1,intarr(2)
idtset=idtset+1
dtsets(idtset)%jdtset=i1*10+i2
end do
end do
</pre>
So, <b>udtset</b>(2) sets the largest value for the unity digit, that varies between 1 and <b>udtset</b>(2).
<br>If <b>udtset</b> is used, the input variable <a href="varbas.html#jdtset">jdtset</a> cannot be used."""
},
'upawu': {
'definition': "value of U for PAW+U ",
'section': "varpaw",
'category': """<a href="../users/abinit_help.html#energy">ENERGY</a>""",
'vartype': """real array upawu(<a href="varbas.html#ntypat">ntypat</a>)""",
'default': "0",
'text': """Relevant only when
<a href="varint.html#usepaw">usepaw</a>=1, and <a
href="varpaw.html#usepawu">usepawu</a>=1.<br>
Gives the value of the
screened coulomb interaction between correlated electrons corresponding
to <a href="varpaw.html#lpawu">lpawu</a> for each species.<br>
In the case where
<a href="varpaw.html#lpawu">lpawu</a>
=-1, the value is not used.
<br>
In the case of a GW calculation, the U interaction defined by <b>upawu</b> will be REMOVED from the self energy.
In particular, for G0 W0 calculations (perturbative calculations), the energy eigenvalues obtained
after an underlying DFT+U calculation will be
<br>
E_GW = E_DFT+U + < phi | Self-energy - U | phi>
<br>
Actually, in order to perform a GW @ DFT+U calculation, one should define the same value of U in the self-energy calculation,
than the one defined in the DFT calculation.
The easiest is actually to define the value of U for the whole set of calculations (for the different datasets),
including the screening, even if the U value does not play explicitly a role in the computation of the latter (well, the input
wavefunctions will be different anyhow).
<br>
It is possible to perform calculations of the type GW+U' @ DFT+U , so keeping a U interaction  (usually smaller than the initial U)
in the GW calculation, by defining a smaller U than the one used in the DFT calculation. This value will be subtracted in the GW correction calculation,
as outlined above.
<br> Explicitly, in order to do a calculation of a material with a DFT U value of 7.5 eV, followed by a GW calculation where there is a residual
U value of 2 eV, one has to define :
<pre>
uldau1   7.5 eV   ! This is for the DFT calculation
...
optdriver4  4
uldau4   5.5 eV   ! This is for the screening calculation
</pre>
"""
},
'use_gpu_cuda': {
'definition': "activate USE of GPU accelerators with CUDA (nvidia)",
'section': "varpar",
'category': "",
'vartype': "integer ",
'default': """1 for ground-state calculations (<a href="vargs.html#optdriver">optdriver</a>=0) when ABINIT has been compiled using cuda, 0 otherwise""",
'text': """<p>
Only available if ABINIT executable has been compiled with cuda nvcc compiler.<br>
This parameter activates the use of NVidia graphic accelerators (GPU) if present.<br>
If <b>use_gp_cuda</b>=1, some parts of the computation are transmitted to the GPUs.<br>
If <b>use_gp_cuda</b>=1, no computation is done on GPUs, even if present.<br><br>
Note that, while running ABINIT on GPUs, it is recommended to use MAGMA external library
(i.e. Lapack on GPUs). The latter is activated during compilation stage (see "configure"
step of ABINIT compilation process). If MAGMA is not used, ABINIT performances on GPUs
can be poor.

"""
},
'use_slk': {
'definition': "USE ScaLapacK",
'section': "varpar",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>  """,
'vartype': "integer parameter  ",
'default': "0",
'text': """"""
},
'usedmatpu': {
'definition': "USE of an initial Density MATrix in Paw+U",
'section': "varpaw",
'category': "",
'vartype': "integer parameter ",
'default': "0",
'text': """Relevant only when
<a href="varint.html#usepaw">usepaw</a>=1, and <a
href="varpaw.html#usepawu">usepawu</a>=1.<br>
When <b>usedmatpu</b>/=0, an initial density matrix (given by <a href="varpaw.html#dmatpawu">dmatpawu</a>
keyword) is used and kept fixed during the first ABS(<b>usedmatpu</b>) SCF steps.<br>
This starting value of the density matrix can be useful to find the correct ground state.
Within LDA+U formalism, finding the minimal energy of the system is tricky; thus it is advised to test several values
of the initial density matrix.<br>
Note also that the density matrix has to respect some symmetry rules determined by the space group.
If the symmetry is not respected in the input, the matrix is however automatically symmetrised.<br><br>
The sign of <b>usedmatpu</b> has influence only when <a href="varrlx.html#ionmov">ionmov</a>/=0 (dynamics or relaxation):<br>
- When <b>usedmatpu</b>>0, the density matrix is kept constant only at first ionic step<br>
- When <b>usedmatpu</b><0, the density matrix is kept constant at each ionic step<br>
"""
},
'usedmft': {
'definition': "USE Dynamical Mean Field Theory ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>  """,
'vartype': "integer parameter  ",
'default': "0 ",
'text': """"""
},
'useexexch': {
'definition': "USE of EXact EXCHange ",
'section': "varpaw",
'category': "",
'vartype': "integer parameter ",
'default': "0",
'text': """When <b>useexexch</b>=1, the hybrid functional PBE0 is used
in PAW, inside PAW spheres only, and only for correlated orbitals given by
<a href="varpaw.html#lexexch">lexexch</a>. To change the ratio of exact exchange, see also <a href="vardev.html#exchmix">exchmix</a>.
"""
},
'usekden': {
'definition': "USE Kinetic energy DENsity",
'section': "vargs",
'category': " ",
'vartype': "integer ",
'default': "0",
'text': """<p>If <b>usekden</b>=1 the kinetic energy density will be computed during the self-consistency loop,
in a way similar to the computation of the density.
This is needed if a meta-GGA is to be used as XC functional. Otherwise (<b>usekden</b>=0), the kinetic energy
density is not computed during the self-consistency loop.</p>
"""
},
'usepaw': {
'definition': "USE Projector Augmented Waves method ",
'section': "varint",
'category': """<a href="../users/abinit_help.html#internal">INTERNAL</a> """,
'vartype': "integer parameter  ",
'default': "set by the pseudopotential files : either PAW (1) or norm-conserving (0).",
'text': """This variable is determined by the pseudopotentials files.
PAW calculations (see <a href="varpaw.html">PAW variables</a>) can only
be performed with PAW atomic data input files, while
pseudopotential calculations are performed in ABINIT with norm-conserving
pseudopotential input files. Most functionalities in ABINIT are available
with either type of calculation.
"""
},
'usepawu': {
'definition': "USE PAW+U (spherical part) ",
'section': "varpaw",
'category': "",
'vartype': "integer parameter ",
'default': "0",
'text': """Must be non-zero if a DFT+U calculation is done, or if a GW calculation following a DFT+U calculation is done (important !).
<ul>
<li> If set to 0, the LDA+U method is not used.<br>
</li>
<li> If set to 1 or 2, the LDA+U method (cf [1])  is used. The full rotationally invariant
formulation is used (see Eq. (3) of Ref [2]) for the interaction term of the energy.
Two choices are allowed concerning the double counting term: <br>
</li>
<ul>
<li> If <b>usepawu</b>=1, the Full Localized
Limit (FLL) (or Atomic limit) double counting is used (cf Eq. (4) of Ref.[2] or
Eq. (8) of Ref[3]).
<br>
</li>
<li> If <b>usepawu</b>=2, the Around Mean Field (AMF) double counting is used
(cf Eq. (7) of Ref [3]). Not valid if nspinor=2.<br>
</li>
</ul>
</ul>
If LDA+U is activated (<b>usepawu</b>=1 or 2), the <a href="varpaw.html#lpawu">lpawu</a>,
<a href="varpaw.html#upawu">upawu</a> and&nbsp; <a
href="varpaw.html#jpawu">jpawu</a>&nbsp; input variables are read. <br>
The implementation is done inside PAW augmentation regions only (cf Ref [4]). The initial density matrix can be
given in the input file (see  <a href="varpaw.html#usedmatpu">usedmatpu</a>).
The expression of the density matrix is chosen thanks to
<a href="varpaw.html#dmatpuopt">dmatpuopt</a>. See also <a href="../users/How_to_use_LDA_plus_U.txt">How_to_use_LDA_plus_U.txt</a>.
for some informations.
<br>
In the case of a GW calculation on top of a DFT+U, the absence of definition of a U value in the self-energy
will LEAVE the underlying U from the DFT calculation.
Thus, the code will actually do a GW+U @ DFT+U calculation.
Note that the screening calculation will not be affected by the presence/absence of a U value.
<br>
Actually, in order to perform a GW @ DFT+U calculation, one should define the same value of U in the self-energy calculation,
than the one defined in the DFT calculation. The code will know that the interaction corresponding to that value has to be
SUBTRACTED inside the self-energy. The easiest is actually to define the presence U for the whole set of calculations (for the different datasets),
including the screening, even if the U value does not play explicitly a role in the computation of the latter (well, the input
wavefunctions will be different anyhow).
<br>
It is possible to perform calculations of the type GW+U' @ DFT+U , so keeping a smaller U interaction in the GW calculation, by subtracting
a smaller U than the one used in the DFT calculation. See the description of the <a href="varpaw.html#upawu">upawu</a> input variable.
<br>

References: <br>
[1] V. I. Anisimov, J. Zaanen, and O. K. Andersen PRB 44, 943 (1991) <br>
[2] A.I. Lichtenstein, V.I. Anisimov and J. Zaanen  PRB 52, 5467 (1995) <br>
[3] M. T. Czyzyk and G.  A. Sawatzky PRB 49, 14211 (1994) <br>
[4] O. Bengone, M. Alouani, P. Blochl, and J. Hugel PRB 62, 16392 (2000) <br>

<br>
<br>
Suggested acknowledgment:<br>
</span>- B. Amadon, F. Jollet and M. Torrent, Phys. Rev. B 77, 155104 (2008).<br>
"""
},
'userec': {
'definition': "USE RECursion ",
'section': "varint",
'category': """<a href="../users/abinit_help.html#internal">INTERNAL</a> """,
'vartype': "integer parameter  ",
'default': "Value is 0",
'text': """This internal variable is set to 1 when the recursion method is
activated (see <a href="vardev.html#tfkinfunc">tfkinfunc</a>).
"""
},
'useri': {
'definition': "USER Integer variables A, B, C, D and E ",
'section': "vardev",
'category': " ",
'vartype': "integers ",
'default': "0 .",
'text': """These are user-definable integers which the user may
input and then utilize in subroutines of his/her own
design.  They are not used in the official versions
of the ABINIT code, and should ease independent
developments (hopefully integrated in the official
version afterwards).
<br>Internally, they are available in the dtset structured datatype,
e.g. dtset%useria .
"""
},
'userr': {
'definition': "USER Real variables A, B, C, D, and E ",
'section': "vardev",
'category': " ",
'vartype': "real numbers",
'default': "0.0 .",
'text': """These are user-definable with the same purpose as <a href="vardev.html#useri">useri</a>. """
},
'usewvl': {
'definition': "Use WaVeLet basis set ",
'section': "varbas",
'category': " ",
'vartype': "integer (0 or 1)",
'default': "0 (use plane-wave basis set).",
'text': """Used to define if the calculation is done on a
wavelet basis set or not.
<br>The values of <b>usewvl</b> must be 0 or 1. Putting <b>usewvl</b>
to 1, makes <a href="vargs.html#icoulomb">icoulomb</a>
mandatory to 1. The number of band (<a
href="varbas.html#nband">nband</a>) must be set manually to
the strict number need for an isolator system (<i>i.e.</i>
number of electron over two). The cut-off is not relevant in the
wavelet case, use <a href="varbas.html#wvl_hgrid">wvl_hgrid</a>
instead.
<br>In wavelet case, the system must be isolated systems (molecules or
clusters). All geometry optimization are available (see <a
href="varrlx.html#ionmov">ionmov</a>, especially the geometry
optimisation and the molecular dynamics.
<br>The spin computation is not currently possible with wavelets and
metalic systems may be slow to converge.
"""
},
'usexcnhat': {
'definition': "USE eXchange-Correlation with NHAT (compensation charge density) ",
'section': "varpaw",
'category': "",
'vartype': "integer parameter ",
'default': "-1",
'text': """Relevant only when <a href="varint.html#usepaw">usepaw</a>=1.<br>
This flag determines how the exchange-correlation terms are computed for the pseudo-density.<br>
When <b>usexcnhat</b>=0, exchange-correlation potential does not include the compensation charge density, i.e.
Vxc=Vxc(tild_Ncore + tild_Nvalence).<br>
When <b>usexcnhat</b>=1, exchange-correlation potential includes the compensation charge density, i.e.
Vxc=Vxc(tild_Ncore + tild_Nvalence + hat_N).<br>
When <b>usexcnhat</b>=-1,the value of <b>usexcnhat</b> is determined from the reading of the PAW dataset file
(pseudopotential file). When PAW datasets with different treatment of Vxc are used in the same
run, the code stops.

"""
},
'useylm': {
'definition': "USE YLM (the spherical harmonics) ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>  """,
'vartype': "integer parameter  ",
'default': """0 for norm-conserving pseudopotential(s), 1 for Projector Augmented-Wave (PAW), 1 when the recursion method is used (<a href="vardev.html#tfkinfunc">tfkinfunc</a>=1).""",
'text': """When this flag is activated, the non-local operator is applied using an algorithm based on spherical harmonics. Non-local projectors are used with their usual form:<br>
<ul>P<sub>lmn</sub>(r)=Y<sub>lm</sub>(r)*p<sub>ln</sub>(r)</ul><br>
<br>When <b>useylm</b>=0, the sum over Y_lm can be reduced to a Legendre polynomial form.
"""
},
'vaclst': {
'definition': "VACancies LiST ",
'section': "vargeo",
'category': """<a href="../users/abinit_help.html#geometry_builder">GEOMETRY BUILDER</a>, <a href="../users/abinit_help.html#not_internal">NOT INTERNAL</a>  """,
'vartype': """integer array vaclst(<a href="vargeo.html#vacnum">vacnum</a>) """,
'default': "No Default.",
'text': """Gives the identification number(s) of atoms
to be subtracted from the set of atoms that are obtained
after having rotated, translated and repeated the objects.
<br>Useful to created vacancies."""
},
'vacnum': {
'definition': "VACancies NUMber ",
'section': "vargeo",
'category': """<a href="../users/abinit_help.html#geometry_builder">GEOMETRY BUILDER</a>  """,
'vartype': "integer parameter  ",
'default': "0.",
'text': """Gives the number of atoms to be subtracted
from the list of atoms after the rotations, translations
and repetitions have been done. The list of these
atoms is contained in <a href="vargeo.html#vaclst">vaclst</a>."""
},
'vacuum': {
'definition': "VACUUM identification ",
'section': "vargs",
'category': """<a href="../users/abinit_help.html#not_internal">NOT INTERNAL</a> """,
'vartype': "integer array <b>vacuum</b>(3)  ",
'default': "No Default",
'text': """Establishes the presence (if 1) or absence (if 0) of a vacuum
layer, along the three possible directions normal to the
primitive axes.
<p>
This information might be used to generate k-point grids,
if <a href="varbas.html#kptopt">kptopt</a>=0 and neither
<a href="varbas.html#ngkpt">ngkpt</a> nor <a href="vargs.html#kptrlatt">kptrlatt</a>
are defined (see explanations with the input variable
<a href="varfil.html#prtkpt">prtkpt</a>).
<br> It will allow to select
a zero-, one-, two- or three-dimensional
grid of k points. The coordinate of the k points
along vacuum directions is automatically set to zero.
<p>
If <b>vacuum</b> is not defined, the input variable
<a href="vargs.html#vacwidth">vacwidth</a>
will be used to determine automatically whether the
distance between atoms is sufficient to have the
presence or absence of vacuum.
"""
},
'vacwidth': {
'definition': "VACuum WIDTH ",
'section': "vargs",
'category': """<a href="../users/abinit_help.html#not_internal">NOT INTERNAL</a>, <a href="../users/abinit_help.html#length">LENGTH</a>  """,
'vartype': "real parameter  ",
'default': "10.0",
'text': """Give a minimum "projected" distance between
atoms to be found in order to declare that there
is some <a href="vargs.html#vacuum">vacuum</a> present for each of the three
directions.
By default, given in Bohr atomic units
(1 Bohr=0.5291772108 Angstroms), although Angstrom can be specified,
if preferred, since <b>vacwidth</b> has the
'<a href="../users/abinit_help.html#dimensions">LENGTH</a>' characteristics.
<br>The precise requirement is that a slab
of width <b>vacwidth</b>, delimited by two
planes of constant reduced coordinates in the
investigated direction, must be empty of atoms.
"""
},
'vcutgeo': {
'definition': "V (potential) CUT-off GEOmetry ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>  """,
'vartype': "real array vcutgeo(3)",
'default': "3*0.0d0",
'text': """<p>
(No explicit documentation, please see
<a href="https://listes-2.sipr.ucl.ac.be/abinit.org/arc/forum/2008-11/msg00087.html">
https://listes-2.sipr.ucl.ac.be/abinit.org/arc/forum/2008-11/msg00087.html </a>
for the time being)
<p>
<b>vcutgeo</b> is used in conjunction with <a href="vargw.html#icutcoul">icutcoul</a>
to specify the geometry used to truncate the Coulomb interaction as well as the
particular approach to be used.
It has a meaning only for the cylindrical symmetry
(<a href="vargw.html#icutcoul">icutcoul</a>=1) or in the case of surfaces
<a href="vargw.html#icutcoul">icutcoul</a>=2.
For each geometry, two different definition of the cutoff region are available
(see Phys. Rev. B 73, 233103 and Phys. Rev. B 73, 205119 for a complete description of the methods)
<p>
In Beigi's method (Phys. Rev. B 73, 233103), the cutoff region is given by the Wigner-Seitz cell
centered on the axis of the cylinder.
The cutoff region is thus automatically defined by the unit cell and there's no need to specify
When <a href="vargw.html#rcut">rcut</a>.
<p>
To define a cylinder along the z-axis use the following two lines.

icutcoul 1
rprim    1 0 0
0 1 0
0 0 1
vcutgeo  0 0 1

<p>
Please note that Beigi's method is implemented only in the case if an orthorombic
Bravais lattic. For hexagonal lattices, one has to use the method of Rozzi (Phys. Rev. B 73, 205119)
In this case, the interaction is truncated in a finite cylinder.
Contrarily to the first approach, here one has to specify both the radius of the cylinder
with <a href="vargw.html#rcut">rcut</a>
as well as the length of the cylinder along the periodic dimension that should always be smaller
than the extension of the Born von Karman box.
The length of the cylinder is given in terms of the primitive vector along the periodic direction.
<p>
For example, in order to define a finite cylinder along z of radius 2.5 Bohr and length 3*R3

icutcoul 1
rprim    1 0 0
0 1 0
0 0 1
vcutgeo  0 0 -3.0 # note the minus sign
rcut     2.5

<p>
For surface calculations (<a href="vargw.html#icutcoul">icutcoul</a>=2),
<b>vcutgeo</b> is used to define the two periodic directions defining the surface.
Also in this case two different techniques are available.
In the method of Beigi, the (positive) non-zero components of vcutgeo define the periodic
directions of the infinite surface. The interaction is truncated within a slab
of width L where L is the length of the primitive vector of the lattice along the non-periodic dimension.
For example:

icutcoul 2
rprim    1 0 0
0 1 0
0 0 1
vcutgeo 1 1 0

It is also possible to define a finite surface by employing negative values
For example:

icutcoul 2
vcutgeo -3 -2 0

defines ....


"""
},
'vdw_nfrag': {
'definition': "van der Waals Number of interacting FRAGments",
'section': "vardev",
'category': "DEVELOP",
'vartype': "integer vdw_nfrag ",
'default': "1",
'text': """The absolute value of vdw_nfrag is the number of vdW interacting
fragments in the unit cell. As wannierization takes place in reciprocal space, the MLWF
center positions could be translated by some lattice vector from the cell where atoms
are placed. If vdw_nfrag >= 1 then MLWFs are translated to the original
unit cell, otherwise the program will keep the positions obtained by Wannier90. The
later is usually correct if some atoms are located at the corners or at limiting
faces of the unit cell.

Used only if <b>vdw_xc=</b>10,11.
"""
},
'vdw_supercell': {
'definition': "Van Der Waals correction from Wannier functions in SUPERCELL",
'section': "vardev",
'category': "DEVELOP",
'vartype': "integer array vdw_supercell(3) ",
'default': "0 0 0",
'text': """Set of dimensionless positive numbers which define the maximum multiples
of the primitive translations (rprimd) in the supercell construction. Each component of vdw_supercell
indicates the maximum number of cells along both positive or negative directions of the corresponding
primitive vector i.e. the components of <a href="varbas.html#rprimd">rprimd</a>. In the case of layered
systems for which vdW interactions occur between layers made of tightly bound atoms, the evaluation
of vdW corrections comming from MLWFs in the same layer (fragment) must be avoided. Both a negative or
null value for one component of <b>vdw_supercell</b>  will indicate that the  corresponding direction
is normal to the layers.


Used only if <b>vdw_xc</b>=10,11.
"""
},
'vdw_tol': {
'definition': "van der Waals TOLerance",
'section': "vardev",
'category': "DEVELOP",
'vartype': "real number ",
'default': "10^-10",
'text': """Used only when Van der Waals DFT-D2 correction is activated (<a href="vardev.html#vdw_xc">vdw_xc</a>=5).<br>
The DFT-D2 (S. Grimme approach) dispersion potential is implemented as a pair potential.
The number of pairs of atoms contributing to the potential is necessarily limited. To be included
in the potential a pair of atom must have contribution to the energy larger than <b>vdw_tol</b>.
"""
},
'vdw_typfrag': {
'definition': "van der Waals TYPe of FRAGment",
'section': "vardev",
'category': "DEVELOP",
'vartype': "integer array vdw_typfrag(natom) ",
'default': "1*natom",
'text': """This array defines the interacting fragments by assigning to each atom an
integer index from 1 to <b>vdw_nfrag</b>. The ordering of <b>vdw_typfrag</b> is the same as
<a href="varbas.html#typat">typat</a> or <a href="varbas.html#xcart">xcart</a>. Internally each MLWF is
assigned to a given fragment by computing the distance to the atoms. MLWFs belong to
the same fragment as their nearest atom. The resulting set of MLWFs in each interacting fragment
can be found in the output file in xyz format for easy visualization.

Used only if <b>vdw_xc</b>=10,11.
"""
},
'vdw_xc': {
'definition': "van der Waals eXchange-Correlation functional",
'section': "vardev",
'category': "DEVELOP",
'vartype': "integer  ",
'default': "0",
'text': """Selects a van-der-Waals density functional to
apply the corresponding correction to the exchange-correlation energy.
If set to zero, no correction will be applied.
<br>
Possible values are:
<ul>
<li>0: no correction.</li>
<li>1: apply vdW-DF1 (DRSLL) from Dion <i>et al.</i><br>
<i>doi:10.1103/PhysRevLett.92.246401</i></li>
<li>2: apply vdw-DF2 (LMKLL) from Lee <i>et al.</i><br>
<i>arXiv:1003.5255v1</i></li>
<li>5: apply vdw-DFT-D2 as proposed by S. Grimme (adding a semi-empirical dispersion potential)<br>
Available only for ground-state calculations ; see <a href="vardev.html#vdw_tol">vdw_tol</a> variable
to control convergency<br>
<i>J. Comp. Chem. 27, 1787 (2006)</i></li>
<li>10: evaluate the vdW correlation energy from maximally localized Wannier functions, as proposed by
P. L. Silvestrelli, also known as vdW-WF1 method. <i>doi:10.1103/PhysRevLett.100.053002,
doi:10.1016/j.cpc.2011.11.003</i></li>
<li>11: evaluate the vdW correlation energy from maximally localized Wannier functions, as proposed by
A. Ambrosetti and P. L. Silvestrelli, also known as vdW-WF2 method. <i>doi:10.1103/PhysRevB.85.073101</i></li>
</ul>
For vdw_xc=1 and vdw_xc=2, the implementation follows the strategy devised
in the article of Rom&aacute;n-P&eacute;rez and Soler
(doi:10.1103/PhysRevLett.103.096102).
"""
},
'vel': {
'definition': "VELocity ",
'section': "varrlx",
'category': """<a href="../users/abinit_help.html#evolving">EVOLVING</a> """,
'vartype': """real array vel(3,<a href="varbas.html#natom">natom</a>) represented internally as vel(3,<a href="varbas.html#natom">natom</a>,<a href="varrlx.html#nimage">nimage</a>)""",
'default': """3*<a href="varbas.html#natom">natom</a> 0's.""",
'text': """Gives the starting velocities
of atoms, in cartesian coordinates, in Bohr/atomic time
units (atomic time units given where <a href="varrlx.html#dtion">dtion</a>
is described).
<br>
Irrelevant unless <a href="varrlx.html#ionmov">ionmov</a> &gt; 0.
<br>
For <a href="varrlx.html#ionmov">ionmov</a>=8 (Nose thermostat),
if <b>vel</b> is not initialized, a random initial
velocity giving the right kinetic energy will be generated.
<br>
If the geometry builder is used, <b>vel</b> will be related
to the preprocessed set of atoms, generated by the
geometry builder. The user must thus foresee the effect
of this geometry builder (see <a href="vargeo.html#objarf">objarf</a>).
<br>
Velocities evolve is <a href="varrlx.html#ionmov">ionmov</a>==1.
"""
},
'vel_cell': {
'definition': "VELocity of the CELL parameters ",
'section': "varrlx",
'category': """<a href="../users/abinit_help.html#evolving">EVOLVING</a> """,
'vartype': """real array vel_cell(3,3) represented internally as vel_cell(3,3,<a href="varrlx.html#nimage">nimage</a>)""",
'default': "3*3 0's.",
'text': """Irrelevant unless <a href="varrlx.html#imgmov">imgmov</a>=9 or 13
and <a href="varrlx.html#optcell">optcell</a>>0 (Path-Integral Molecular Dynamics
with NPT algorithm).<br>
Gives the starting velocities of the dimensional cell parameters in Bohr/atomic time
units (atomic time units given where <a href="varrlx.html#dtion">dtion</a>
is described).
"""
},
'vis': {
'definition': "VIScosity ",
'section': "varrlx",
'category': "",
'vartype': "real parameter ",
'default': "100.",
'text': """The equation of motion is :<br>
M<sub>I</sub> d<sup>2</sup>R<sub>I</sub>/dt<sup>2</sup>= F<sub>I</sub>
- <b>vis</b> dR<sub>I</sub>/dt
<br>
<br>
The atomic unit of viscosity is hartrees*(atomic time units)/Bohr<sup>2</sup>.
Units are not
critical as this is a fictitious damping used to relax
structures. A typical value for silicon is 400 with
<a href="varrlx.html#dtion">dtion</a> of 350 and atomic mass 28 <a
href="varrlx.html#amu">amu</a>. Critical
damping is most desirable and is found only by
optimizing <b>vis</b> for a given situation.
"""
},
'vprtrb': {
'definition': "potential -V- for the PeRTuRBation ",
'section': "varff",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>, <a href="../users/abinit_help.html#energy">ENERGY</a>  """,
'vartype': "real array of 2 elements  ",
'default': "0.d0 0.d0.",
'text': """Gives the real and imaginary
parts of a scalar potential perturbation.
Can be specified in Ha (the default), Ry, eV or Kelvin, since
<b>ecut</b> has the
'<a href="../users/abinit_help.html#dimensions">ENERGY</a>' characteristics.
<br>This is made
available for testing responses to such perturbations.
The form of the perturbation, which is added to the local
potential, is:
<ul>
<li> (<b>vprtrb</b>(1)+I*<b>vprtrb</b>(2))/2 at G=<a href="varff.html#qprtrb">qprtrb</a>  and</li>
<li> (<b>vprtrb</b>(1)-I*<b>vprtrb</b>(2))/2 at G=-<a href="varff.html#qprtrb">qprtrb</a>
(see <a href="varff.html#qprtrb">qprtrb</a> also).</li>
</ul>"""
},
'w90iniprj': {
'definition': "Wannier90- INItial PROJections ",
'section': "varw90",
'category': " ",
'vartype': "integer ",
'default': "1 (random projections).",
'text': """In order to find the Maximally Localized Wannier Functions, the user
has to provide an initial guess. A set of localized trial orbitals
is chosen
corresponding to some rough initial guess at the
Wannier Functions, and these are projected onto the  Bloch
eigenstates. See Ivo
Souza, Nicola Marzari, and David Vanderbilt. Phys. Rev. B, 65, 035109 (2001).
<br>These initial projections are stored in a file .amn and the variable
<b> w90iniprj</b> is used to construct them:
<ul>
<li><b> w90iniprj</b>=1:
Random projections.<br>
<br>
</li>

<li><b> w90iniprj</b>=2:
The initial projections will be a linear combination of hydrogenic
atomic orbitals.
<br>The user has to define the projections in the secondary input file
wannier90.win
<br>Information about how to define them can be found in the manual of
Wannier90. See <span><a href="http://www.wannier.org">www.wannier.org</a>
</ul>
Not read if <a href="varfil.html#prtwant">prtwant</a> /= 2 or 3.
"""
},
'w90prtunk': {
'definition': "Wannier90- PRINT UNKp.s file",
'section': "varw90",
'category': " ",
'vartype': "integer ",
'default': "set to zero.",
'text': """Defines whether or not the UNKp.s file will be printed.
<ul>
<li> <b>w90prtunk</b>=0: Do not print the UNKp.s files<br>
<br>
</li>
<li> <b>w90prtunk</b>=1: Print the UNKp.s files on a fine grid<br>
<br>
</li>
<li> <b>w90prtunk</b>>1: Print the UNKp.s files on a coarse grid <br>
Instead of printing every record we will print every w90prtunk records. This is useful to reduce the size of the UNKp.s files, but, the quality is also reduced.</li>

</ul>
Not read if <a href="varfil.html#prtwant">prtwant</a> /= 2 or 3.<br>
<br>
The default is set to zero because UNKp.s files occupy a lot of
memory.
These files contain the periodic part of the bloch states represented
on a regular real space grid. They are indexed by k-point <b>p</b> (from 1 to
nkpt) and spin <b>s</b> ('1' for 'up','2' for 'down').
<br> <br> The name of the wavefunction file is assumed to have the form:
<br><br>
write(wfnname,200) <b>p</b>, <b>spin</b>
<br>200 format ('UNK',i5.5,'.',i1)
<br><br>
These file are unformatted.
The first line of each file contains 5 integers: the number of
grid points in each direction (<b>n1</b>, <b>n2</b> and <b>n3</b>), the k-point number <b>ikpt</b>
and the total number of bands mband in the file. The following rows contain the wavefunctions in real space.
<p>
These files are written in the following way for the coarse grid:</p>
<pre>
write(iun_plot) n1/w90prtunk,n2/w90prtunk,n3/w90prtunk,ikpt,nband
write(iun_plot) (((fofr(1,jj1,jj2,jj3),fofr(2,jj1,jj2,jj3),&
&      jj1=1,n1,w90prtunk),jj2=1,n2,w90prtunk),jj3=1,n3,w90prtunk)
</pre>
Where <b>fofr</b> is a double precision variable which contains the wavefunctions in real space.
Note that in order to reduce the size of the UNK files we are just
including records in the wavefunctions for 1/(w90prtunk^3) of the grid points.
That is why we divide n1, n2 and n3 by prtunk. The output .xsf files for plotting
with XCrysDen will also be on the coarse grid.  When this dosen't produce an
acceptable plot, prtunk can be set to 1 to output every grid point.
(You should try spline interpolation in XCrysDen first.)



</ul>
"""
},
'wfoptalg': {
'definition': "WaveFunction OPTimisation ALGorithm ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>  """,
'vartype': "integer parameter  ",
'default': """0 when <a href="varint.html#usepaw">usepaw</a>=0 (norm-conserving pseudopotentials),<br>10 when <a href="varint.html#usepaw">usepaw</a>=1 (PAW).<br>Default is 14 if <a href="varpar.html#paral_kgb">paral_kgb</a>=1.""",
'text': """Allows to choose the algorithm
for the optimisation of the wavefunctions.
<br>The different possibilities are :
<ul>
<li><b>wfoptalg</b>=0 : standard state-by-state conjugate gradient algorithm,
with no possibility to parallelize over the states;</li>
<li><b>wfoptalg</b>=2 : minimisation of the residual with respect
to different shifts, in order to cover the whole set of occupied
bands, with possibility to parallelize over blocks of states (or bands).
The number of states in a block
is defined in <a href="vardev.html#nbdblock">nbdblock</a>.
THIS IS STILL IN DEVELOPMENT.</li>
<li><b>wfoptalg</b>=3 : minimisation of the residual with respect
to a shift. Available only in the non-self-consistent case
<a href="varbas.html#iscf">iscf</a>=-2,
in order to find eigenvalues and wavefunctions close to a
prescribed value.</li>
<li><b>wfoptalg</b>=4 : (see also <b>wfoptalg</b>=14), a parallel code based on the Locally Optimal
Block Preconditioned Conjugate Gradient (LOBPCG) method of Knyazev.
<a href="http://dx.doi.org/10.1137/S1064827500366124">
Reference : A.V. Knyazev, "Toward the Optimal Preconditioned Eigensolver
: Locally Optimal Block Preconditioned Conjugate Gradient Method". SIAM
Journal on Scientific Computing 23, pp517-541 (2001)</a>.
The implementation rests
on the <a href="http://www.mathworks.com/matlabcentral/fileexchange/48-lobpcg-m">
matlab program by Knyazev</a>.
<a href="http://dx.doi.org/10.1137/060661624">
Reference A. V. Knyazev, I. Lashuk, M. E. Argentati, and E. Ovchinnikov,
Block Locally Optimal Preconditioned Eigenvalue Xolvers (BLOPEX) in
hypre and PETSc (2007). SIAM Journal on Scientific Computing (SISC).
25(5): 2224-2239</a>.
For more
information see
<a href="http://dx.doi.org/10.1016/j.commatsci.2007.07.019">
F. Bottin, S. Leroux, A. Knyazev, G. Zerah, Large scale
ab initio calculations based on three levels of parallelization. (2008).
Computational Material Science, 42(2), 329-336.
</a>
<li><b>wfoptalg</b>=10 : (for PAW) standard state-by-state conjugate gradient algorithm,
with no possibility to parallelize over the states, but modified
scheme described in Kresse, Furthmuller, PRB 54, 11169 (1996)
(modified kinetic energy, modified preconditionning, minimal
orthogonalization, ...) ;</li>
<li><b>wfoptalg</b>=14 :
the recommended for massively parallel code, the same as <b>wfoptalg</b>=4 except that the preconditioning of
the block vectors does not depend on the kinetic energy of each band,
and the orthogonalization after the LOBPCG algorithm is no longer
performed. The first modification increases the convergence and the
second one the efficiency.
</li>
</ul>
"""
},
'wtatcon': {
'definition': "WeighTs for AToms in CONstraint equations",
'section': "varrlx",
'category': """<a href="../users/abinit_help.html#no_multi">NO MULTI</a> """,
'vartype': """real array wtatcon(3,<a href="varrlx.html#natcon">natcon</a>,<a href="varrlx.html#nconeq">nconeq</a>)""",
'default': "0.",
'text': """Gives the weights determining how the motion of atoms
is constrained
during structural optimization or molecular dynamics (see <a
href="varrlx.html#nconeq">nconeq</a>, <a href="varrlx.html#natcon">natcon</a>,
and <a href="varrlx.html#iatcon">iatcon</a>). For each of the <a
href="varrlx.html#nconeq">nconeq</a> independent constraint equations,
wtatcon is a 3*<a href="varrlx.html#natcon">natcon</a> array giving
weights, W<sub>I</sub>,
for the x, y, and z components of each of the atoms (labeled by I) in
the list of indices <a href="varrlx.html#iatcon">iatcon</a>.
Prior to taking an atomic step, the calculated forces, F<sub>I</sub>,
are
replaced by projected forces, F'<sub>I</sub>, which satisfy the set of
constraint equations
<br>
<br>
Sum<sub>mu=x,y,z; I=1,natcon</sub>: W<sub>mu,I</sub> * F'<sub>mu,I</sub>
= 0 for each of the <a href="varrlx.html#nconeq">nconeq</a> arrays W<sub>I</sub>.
<br>
<br>
Different types of motion constraints can be implemented this way. For
example,
<br>
<br>
nconeq 1 natcon 2 iatcon 1 2 wtatcon 0 0 +1 0 0 -1
<br>
<br>
could be used to constrain the relative height difference of two
adsorbate atoms on a surface (assuming their
masses are equal), since F'<sub>z,1</sub> - F'<sub>z,2</sub> = 0
implies z<sub>1</sub> - z<sub>2</sub> = constant.
"""
},
'wtk': {
'definition': "WeighTs for K points ",
'section': "varbas",
'category': " ",
'vartype': """real array wtk(<a href="varbas.html#nkpt">nkpt</a>)  """,
'default': """<a href="varbas.html#nkpt">nkpt</a>*1.0d0 .""",
'text': """Gives the k point weights.  <br>The
k point weights will have their sum normalized to 1
(unless <a href="varbas.html#occopt">occopt</a>=2; see description of <a href="varbas.html#occopt">occopt</a>)
within the program and therefore may be input with any
arbitrary normalization.  This feature helps avoid the
need for many digits in representing fractional weights
such as 1/3.
<br><b>wtk</b> is ignored if <a href="varbas.html#iscf">iscf</a> is not positive,
except if <a href="varbas.html#iscf">iscf</a>=-3."""
},
'wvl_crmult': {
'definition': "WaVeLet Coarse grid Radius MULTiplier ",
'section': "vargs",
'category': "",
'vartype': "real parameter  ",
'default': "6.0",
'text': """This factor is used to defined the expansion of the coarse resolution
grid in the case of wavelets (see <a
href="varbas.html#usewvl">usewvl</a>). The grid is made of
points inside spheres centered on atoms. The radius of these
spheres are the product between this factor and the covalent
radius of element (read from the pseudo-potential file).<br />
This factor is responsible for the amount of used memory (see also <a href="varbas.html#wvl_hgrid">wvl_hgrid</a>).
"""
},
'wvl_frmult': {
'definition': "WaVeLet Fine grid Radius MULTiplier ",
'section': "vargs",
'category': "",
'vartype': "real parameter  ",
'default': "10.0",
'text': """This factor is used to defined the expansion of the fine resolution
grid in the case of wavelets (see <a
href="varbas.html#usewvl">usewvl</a>). This fine resolution
grid has the same grid step than the coarse one (see <a
href="vargs.html#wvl_crmult">wvl_crmult</a>), but on each
point, 8 coefficients are stored instead of one, increasing the
precision of the calculation in this area. The grid is made of
points inside spheres centered on atoms. The radius of these
spheres are the product between this factor and a value read from the pseudo-potential file.<br />
This factor is responsible for the amount of used memory (see also <a href="varbas.html#wvl_hgrid">wvl_hgrid</a>).
"""
},
'wvl_hgrid': {
'definition': "WaVeLet H step GRID ",
'section': "varbas",
'category': """<a href="../users/abinit_help.html#length">LENGTH</a>  """,
'vartype': "real parameter ",
'default': "0.5d0 .",
'text': """It gives the step size in real space for the
grid resolution in the wavelet basis set. This value is highly
responsible for the memory occupation in the wavelet
computation. The value is a length in atomic units.
"""
},
'wvl_nprccg': {
'definition': "WaVeLet maximum Number of PReConditionner Conjugate Gradient iterations",
'section': "vargs",
'category': "",
'vartype': "integer parameter  ",
'default': "5",
'text': """In the wavelet computation case, the wavefunctions are directly
minimised using a real-space preconditionner. This preconditionner
has internally some conjugate gradient iterations. This value
defines a boundary for the number of conjugate gradient
iterations on each wavefunction convergence step.
"""
},
'xangst': {
'definition': "vectors (X) of atom positions in cartesian coordinates -length in ANGSTrom- ",
'section': "varbas",
'category': """<a href="../users/abinit_help.html#not_internal">NOT INTERNAL</a> """,
'vartype': """real array xangst(3,<a href="varbas.html#natom">natom</a>) (or xangst(3,<a href="vargeo.html#natrd">natrd</a>) if the geometry builder is used)""",
'default': "",
'text': """Gives the cartesian coordinates
of atoms within unit cell, in angstrom.  This information is
redundant with that supplied by array <a href="varbas.html#xred">xred</a> or <a href="varbas.html#xcart">xcart</a>.
<br>If <a href="varbas.html#xred">xred</a> and <b>xangst</b> are ABSENT from the input file and
<a href="varbas.html#xcart">xcart</a> is
provided, then the values of <a href="varbas.html#xred">xred</a> will be computed from
the provided <a href="varbas.html#xcart">xcart</a> (i.e. the user may use xangst instead
of <a href="varbas.html#xred">xred</a> or <a href="varbas.html#xcart">xcart</a> to provide starting coordinates).
<br>One and only one of <a href="varbas.html#xred">xred</a>, <a href="varbas.html#xcart">xcart</a>
and <b>xangst</b> must be provided.
<br>The conversion factor between Bohr and Angstrom
is 1 Bohr=0.5291772108 Angstrom, see the <a href="http://physics.nist.gov/cuu/Constants/index.html">NIST site</a>.
<br>Atomic positions evolve if <a href="varrlx.html#ionmov">ionmov</a>/=0 .
In constrast with <a href="varbas.html#xred">xred</a> and
<a href="varbas.html#xcart">xcart</a>, <b>xangst</b> is not internal.
"""
},
'xc_denpos': {
'definition': "eXchange-Correlation - DENsity POSitivity value  ",
'section': "vardev",
'category': """<a href="../users/abinit_help.html#develop">DEVELOP</a>  """,
'vartype': "real",
'default': "1.0e-14 ",
'text': """For the evaluation of the exchange-correlation functionals, the density
cannot be negative, or even too small (e.g. the LDA exchange kernel
behaves like the density at power -(2/3), and the density is used at the denominator
of different factors in GGAs and metaGGAs.
<b>xc_denpos</b> is the smallest value that the density can assume at the time of the
evaluation of a XC functional, in ABINIT. When then computed density drops below <b>xc_denpos</b>
before attacking the evaluation of the XC functional, then it will be (only for that purpose)
replaced by <b>xc_denpos</b>. Note that the evaluation of the gradients or other quantities
that are density-dependent is performed before this replacement.
<p>
It has been observed that the SCF cycle of the Tran-Blaha mGGA can be quite hard to make converge, for systems
for which there is some vacuum. In this case, setting <b>xc_denpos</b> to 1.0e-7 ... 1.0e-6 has been seen
to allow good convergence. Of course, this will affect the numerical results somehow, and one should play
a bit with this value to avoid incorrect calculations.
"""
},
'xc_tb09_c': {
'definition': "Value of the c parameter in the eXchange-Correlation TB09 functional ",
'section': "vardev",
'category': "</a> ",
'vartype': "real ",
'default': "all 99.99d0",
'text': """The modified Becke-Johnson exchange-correlation functional by Tran and Blaha (Phys. Rev. Lett. 102, 226401 (2009)) reads :
<p>V_x(r) = c * V_x^{BR}(r) + (3*c - 2) * 1/pi * sqrt(5/12) * sqrt(2*kden(r)/den(r))
<p>in which V_x^{BR}(r) is the Becke-Roussel potential.
<p>In this equation the parameter c can be evaluated at each SCF step according to the following equation :
<p>c = alpha + beta * sqrt(1/V_{cell} * \int_{V_{cell}} |grad(den(r))|/den(r) d3r)
<p>The c parameter is evaluated thanks to the previous equation when xc_tb09_c is equal to the "magic" default value 99.99.
The c parameter can also be fixed to some (property-optimized or material-optimized) value by using this variable.
"""
},
'xcart': {
'definition': "vectors (X) of atom positions in CARTesian coordinates  ",
'section': "varbas",
'category': """<a href="../users/abinit_help.html#evolving">EVOLVING</a>, <a href="../users/abinit_help.html#length">LENGTH</a> """,
'vartype': """real array xcart(3,<a href="varbas.html#natom">natom</a>) (or xcart(3,<a href="vargeo.html#natrd">natrd</a>) if the geometry builder is used) """,
'default': "",
'text': """Gives the cartesian coordinates
of atoms within unit cell.  This information is
redundant with that supplied by array <a href="varbas.html#xred">xred</a> or <a href="varbas.html#xangst">xangst</a>.
By default, <b>xcart</b> is given in Bohr atomic units
(1 Bohr=0.5291772108 Angstroms), although Angstrom can be specified,
if preferred, since <b>xcart</b> has the
'<a href="../users/abinit_help.html#dimensions">LENGTH</a>' characteristics.
<br>If <a href="varbas.html#xred">xred</a> and <a href="varbas.html#xangst">xangst</a> are
ABSENT from the input file and <b>xcart</b> is
provided, then the values of <a href="varbas.html#xred">xred</a> will be computed from
the provided <b>xcart</b> (i.e. the user may use <b>xcart</b> instead
of <a href="varbas.html#xred">xred</a> or <a href="varbas.html#xangst">xangst</a> to provide starting coordinates).
<br>One and only one of <a href="varbas.html#xred">xred</a>, <a href="varbas.html#xcart">xcart</a>
and <b>xangst</b> must be provided.
<br>Atomic positions evolve if <a href="varrlx.html#ionmov">ionmov</a>/=0 ."""
},
'xclevel': {
'definition': "eXchange Correlation functional level",
'section': "vargs",
'category': """<a href="../users/abinit_help.html#internal">INTERNAL</a>""",
'vartype': "integer parameter  ",
'default': "0",
'text': """Automatically determined from the value of <a href="varbas.html#ixc">ixc</a>.
<ul>
<li>0 => No XC contribution.</li>
<li>1 => LDA functional.</li>
<li>2 => GGA functional.</li>
<li>3 => Functional for TDDFT.</li>
</ul>
"""
},
'xred': {
'definition': "vectors (X) of atom positions in REDuced coordinates  ",
'section': "varbas",
'category': """<a href="../users/abinit_help.html#evolving">EVOLVING</a> """,
'vartype': """real array xred(3,<a href="varbas.html#natom">natom</a>) (or xred(3,<a href="vargeo.html#natrd">natrd</a>) if the geometry builder is used), represented internally as xred(3,<a href="varbas.html#natom">natom</a>,<a href="varrlx.html#nimage">nimage</a>) """,
'default': "all 0.0d0",
'text': """Gives the atomic locations within
unit cell in coordinates relative to real space primitive
translations (NOT in cartesian coordinates).  Thus these
are fractional numbers typically between 0 and 1 and
are dimensionless.  The cartesian coordinates of atoms (in Bohr)
are given by:<br>
<tele>  R_cartesian = xred1*rprimd1+xred2*rprimd2+xred3*rprimd3</tele><br>
where (t1,t2,t3) are the "reduced coordinates" given in
columns of "<b>xred</b>", (rprimd1,rprimd2,rprimd3) are the columns of
primitive vectors array "<a href="varbas.html#rprimd">rprimd</a>" in Bohr.
<br>If you prefer to work only with cartesian coordinates, you
may work entirely with "<a href="varbas.html#xcart">xcart</a>" or "<a href="varbas.html#xangst">xangst</a>" and ignore <b>xred</b>, in
which case <b>xred</b> must be absent from the input file.
<br>One and only one of <a href="varbas.html#xred">xred</a>, <a href="varbas.html#xcart">xcart</a>
and <b>xangst</b> must be provided.
<br>Atomic positions evolve if <a href="varrlx.html#ionmov">ionmov</a>/=0 ."""
},
'xyzfile': {
'definition': "XYZ FILE input for geometry ",
'section': "vargeo",
'category': """<a href="../users/abinit_help.html#not_internal">NOT INTERNAL</a>""",
'vartype': "character file name  ",
'default': "No default",
'text': """Gives the name of a xyz format file, to take
<a href="varbas.html#natom">natom</a>, <a href="varbas.html#ntypat">ntypat</a>, <a href="varbas.html#typat">typat</a>, <a href="varbas.html#znucl">znucl</a>,
and <a href="varbas.html#xangst">xangst</a> from. This input can not be mixed with normal atom specifications or <a href="varfil.html#cmlfile">cmlfile</a> for other datasets."""
},
'zcut': {
'definition': "Z-CUT ",
'section': "vargw",
'category': """<a href="../users/abinit_help.html#gw">GW</a>, <a href="../users/abinit_help.html#energy">ENERGY</a> """,
'vartype': "real parameter ",
'default': "0.1 eV = 3.67493260d-03Ha ",
'text': """<p>
Only relevant if <a href="vargs.html#optdriver">optdriver</a>=3,4,99 that is, screening, sigma or BS calculations.
<p>
It is meant to avoid some divergencies that might occur during the evaluation of the Adler-Wiser expression of
the irreducible polarizability (<a href="vargs.html#optdriver">optdriver</a>=3) or during the numerical treatment
of the integrals defining the contribution to the self-energy matrix elements
(<a href="vargs.html#optdriver">optdriver</a>=4).
If the denominator becomes smaller than <b>zcut</b>, a small imaginary part (depending on <b>zcut</b>) is added,
in order to avoid the divergence.
<p>
When <a href="vargs.html#optdriver">optdriver</a>=99, <b>zcut</b> defines the small complex shift
used to avoid divergences in the expression for the macroscopic dieletric function.
It simulates the experimental uncertainty and the finite lifetime of the quasiparticles
(although the true lifetime should be k- and band-dependent).
The value of <b>zcut</b> affects the number of iteration needed to achieve convergence
in the Haydock iterative method. In this case, <b>zcut</b> should be
larger than the typical distance between the eigenvalues of the exciton Hamiltonian.
<br>
Ideally, one should make a convergence study decreasing the value of <b>zcut</b> for increasing number of k-points.
"""
},
'zeemanfield': {
'definition': "ZEEMAN FIELD",
'section': "varff",
'category': "MAGNETIC FIELD",
'vartype': "real array zeemanfield(3)  ",
'default': "0",
'text': """Give the value of the Zeeman field, H, acting on the spinorial wavefunctions.
Note that Tesla are admitted. This sets the magnitude of mu_0*H, in Tesla,
with H in Amperes/metre.
"""
},
'ziontypat': {
'definition': "Z (charge) of the IONs for the different TYPes of AToms",
'section': "varint",
'category': " ",
'vartype': "real array ziontypat(ntypat)",
'default': "value is set by the pseudopotential files.",
'text': """Charge of the pseudo-ion (=number of valence electrons
that are needed to screen exactly the pseudopotential).
"""
},
'znucl': {
'definition': "charge -Z- of the NUCLeus ",
'section': "varbas",
'category': """<a href="../users/abinit_help.html#no_multi">NO MULTI</a>  """,
'vartype': """real array znucl(<a href="vargs.html#npsp">npsp</a>)  """,
'default': "",
'text': """Gives nuclear charge for each
type of pseudopotential, in order.
<br>If <b>znucl</b> does not agree with nuclear charge,
as given in pseudopotential files, the program writes
an error message and stops.
<p>N.B. : In the pseudopotential files, <b>znucl</b> is called "zatom".
<p>For a "dummy" atom, with znucl=0 , as used in the case of calculations
with only a jellium surface, ABINIT sets arbitrarily the covalent radius to one.
"""
}
}
